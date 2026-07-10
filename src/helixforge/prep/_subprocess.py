"""Shared subprocess driver for the prep wrappers."""

from __future__ import annotations

import subprocess
import time
from collections.abc import Callable, Iterable, Sequence
from logging import Logger
from pathlib import Path
from typing import IO

from helixforge.constants import STDERR_TAIL as _STDERR_TAIL
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)

# _STDERR_TAIL: trailing stderr characters kept in a failed-tool exception message.


class ToolError(RuntimeError):
    """An external tool exited nonzero. Carries ``returncode`` + ``argv``.

    Subclasses ``RuntimeError`` so existing ``pytest.raises(RuntimeError)`` and
    callers keep working unchanged. The structured ``returncode`` is what the
    bounded-retry wrapper gates on — only an **allowlisted** exit code from an
    **idempotent** step is retried.
    """

    def __init__(
        self,
        message: str,
        *,
        returncode: int | None = None,
        argv: Sequence[str] | None = None,
    ) -> None:
        super().__init__(message)
        self.returncode = returncode
        self.argv: list[str] | None = list(argv) if argv is not None else None


def with_retries(
    call: Callable[[], object],
    *,
    retries: int = 0,
    retry_exit_codes: Iterable[int] | None = None,
    backoff: float = 0.5,
    _sleep: Callable[[float], object] = time.sleep,
    log: Logger | None = None,
) -> object:
    """Run ``call`` with bounded retry-with-backoff on a :class:`ToolError`.

    Wraps **only idempotent** external steps
    (DIAMOND, alignment, StringTie — they fully regenerate their output) so a
    transient cluster failure (FS hiccup, scratch contention, momentary OOM) does
    not kill a multi-hour run. A retry happens only while attempts remain **and**
    the exit code is allowlisted (``retry_exit_codes=None`` ⇒ any nonzero;
    a tuple ⇒ exactly those). A missing binary (``FileNotFoundError``) is **not**
    a ``ToolError`` and so is never retried. Backoff is linear (``backoff * n``).
    **Never** wrap a non-idempotent step (e.g. Mikado ``serialise``) — checkpoint
    it instead.
    """
    _log_use = log or _log
    attempt = 0
    while True:
        try:
            return call()
        except ToolError as exc:
            retry_set = set(retry_exit_codes) if retry_exit_codes is not None else None
            allowlisted = retry_set is None or exc.returncode in retry_set
            if attempt < retries and allowlisted:
                attempt += 1
                delay = backoff * attempt
                argv0 = exc.argv[0] if exc.argv else "?"
                _log_use.warning(
                    "transient failure (%s exit %s) — retry %d/%d in %.1fs",
                    argv0,
                    exc.returncode,
                    attempt,
                    retries,
                    delay,
                )
                _sleep(delay)
                continue
            raise


def output_is_fresh(
    output: str | Path,
    inputs: Iterable[str | Path] = (),
) -> bool:
    """True if ``output`` already exists, is non-empty, and is up to date.

    Skip-if-exists idempotency for the prep wrappers (Phase 22 §3.1): a re-run of
    a long maize-scale prep should not re-align every sample. ``output`` is fresh
    when it exists, is non-empty, and is **newer than** every input in ``inputs``
    that exists (an input modified after the output means the output is stale and
    must be regenerated). Missing/empty output ⇒ not fresh ⇒ regenerate.
    """
    output = Path(output)
    if not output.exists() or output.stat().st_size == 0:
        return False
    out_mtime = output.stat().st_mtime
    for inp in inputs:
        ip = Path(inp)
        if ip.exists() and ip.stat().st_mtime > out_mtime:
            return False
    return True


def _write_stderr_log(stderr_log: str | Path | None, stderr: str | None) -> None:
    """Stream a tool's *full* stderr to ``stderr_log``.

    The exception message keeps only the tail (``_STDERR_TAIL``); the head of a
    long error — often the real cause — is preserved here. Best-effort: a failure
    to write the diagnostic log never masks the underlying tool error.
    """
    if not stderr_log or not stderr:
        return
    try:
        Path(stderr_log).parent.mkdir(parents=True, exist_ok=True)
        Path(stderr_log).write_text(stderr)
    except OSError:
        pass


def _run_tool_once(
    argv: list[str],
    *,
    cwd: str | Path | None,
    stdout_path: str | Path | None,
    stderr_log: str | Path | None,
    log: Logger,
) -> subprocess.CompletedProcess[str]:
    """One attempt of :func:`run_tool`; raise :class:`ToolError` on a nonzero exit."""
    log.info("run: %s", " ".join(argv))
    t0 = time.perf_counter()
    cwd_s = str(cwd) if cwd is not None else None
    stdout_fh: IO[str] | None = (
        open(stdout_path, "w") if stdout_path is not None else None
    )
    try:
        if stdout_fh is not None:
            proc = subprocess.run(
                argv,
                check=True,
                stdout=stdout_fh,
                stderr=subprocess.PIPE,
                text=True,
                cwd=cwd_s,
            )
        else:
            proc = subprocess.run(
                argv,
                check=True,
                capture_output=True,
                text=True,
                cwd=cwd_s,
            )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"tool not found: {argv[0]!r}. Install it and put it on PATH, "
            "or pass an absolute path via the wrapper's *_bin argument."
        ) from exc
    except subprocess.CalledProcessError as exc:
        _write_stderr_log(stderr_log, exc.stderr)
        stderr_tail = (exc.stderr or "")[-_STDERR_TAIL:]
        raise ToolError(
            f"external tool {argv[0]!r} failed (exit {exc.returncode}).\n"
            f"argv: {argv}\n"
            f"stderr tail:\n{stderr_tail}",
            returncode=exc.returncode,
            argv=argv,
        ) from exc
    finally:
        if stdout_fh is not None:
            stdout_fh.close()
    _write_stderr_log(stderr_log, getattr(proc, "stderr", None))
    log.info("done in %.2fs: %s", time.perf_counter() - t0, argv[0])
    return proc


def run_tool(
    argv: Sequence[str | Path | int],
    *,
    cwd: str | Path | None = None,
    stdout_path: str | Path | None = None,
    stderr_log: str | Path | None = None,
    log: Logger | None = None,
    retries: int = 0,
    retry_exit_codes: Iterable[int] | None = None,
    backoff: float = 0.5,
    _sleep: Callable[[float], object] = time.sleep,
) -> subprocess.CompletedProcess[str]:
    """Run ``argv`` (a list) with logging + timing; raise on failure.

    ``check=True``, output captured. ``cwd`` sets the working directory (some
    tools write to CWD). ``stdout_path`` redirects stdout to that file instead
    of capturing it (stderr is still captured for error reporting). ``stderr_log``
    streams the tool's *full* stderr to that file (head + tail), while the
    exception message keeps only the ``_STDERR_TAIL`` tail. On a nonzero exit,
    re-raise as :class:`ToolError` (a ``RuntimeError``) with the argv + stderr
    tail; on a missing binary, raise a clear ``tool not found: {argv[0]}`` error.

    ``retries``: when > 0 and the step is **idempotent**, an
    allowlisted nonzero exit is retried with linear backoff (see
    :func:`with_retries`). Default ``retries=0`` reproduces the original
    single-attempt behavior byte-for-byte — never retry a non-idempotent step.
    """
    _log_use = log or _log
    argv_strs: list[str] = [str(a) for a in argv]

    def _attempt() -> subprocess.CompletedProcess[str]:
        return _run_tool_once(
            argv_strs,
            cwd=cwd,
            stdout_path=stdout_path,
            stderr_log=stderr_log,
            log=_log_use,
        )

    result = with_retries(
        _attempt,
        retries=retries,
        retry_exit_codes=retry_exit_codes,
        backoff=backoff,
        _sleep=_sleep,
        log=_log_use,
    )
    return result  # type: ignore[return-value]  # with_retries returns object; _attempt always returns CompletedProcess


def run_pipe(
    commands: Sequence[Sequence[str | Path | int]],
    *,
    cwd: str | Path | None = None,
    stdout_path: str | Path | None = None,
    log: Logger | None = None,
) -> list[subprocess.Popen[str]]:
    """Run ``commands`` (a list of argv lists) connected by OS pipes.

    ``commands[i]``'s stdout feeds ``commands[i+1]``'s stdin (``cmd0 | cmd1 |
    ...``). The last command's stdout goes to ``stdout_path`` if given, else to
    the parent's stdout. Each command's stderr is captured; any nonzero exit
    raises ``RuntimeError`` with the offending argv + stderr tail. A missing
    binary raises the same ``tool not found`` error as :func:`run_tool`.
    Returns the list of completed ``Popen`` objects.
    """
    _log_use = log or _log
    cmds: list[list[str]] = [[str(a) for a in argv] for argv in commands]
    _log_use.info("pipe: %s", " | ".join(" ".join(c) for c in cmds))
    cwd_s = str(cwd) if cwd is not None else None
    final_fh: IO[str] | None = (
        open(stdout_path, "w") if stdout_path is not None else None
    )
    procs: list[subprocess.Popen[str]] = []
    try:
        prev: subprocess.Popen[str] | None = None
        for i, argv in enumerate(cmds):
            is_last = i == len(cmds) - 1
            stdout_arg: int | IO[str] | None
            if is_last:
                stdout_arg = final_fh if final_fh is not None else None
            else:
                stdout_arg = subprocess.PIPE
            try:
                proc: subprocess.Popen[str] = subprocess.Popen(
                    argv,
                    stdin=(prev.stdout if prev is not None else None),
                    stdout=stdout_arg,
                    stderr=subprocess.PIPE,
                    text=True,
                    cwd=cwd_s,
                )
            except FileNotFoundError as exc:
                raise FileNotFoundError(
                    f"tool not found: {argv[0]!r}. Install it and put it on "
                    "PATH, or pass an absolute path via the wrapper's *_bin."
                ) from exc
            if prev is not None:
                # Let the upstream process receive SIGPIPE if we exit early.
                prev.stdout.close()  # type: ignore[union-attr]  # prev.stdout is PIPE so not None
            procs.append(proc)
            prev = proc

        errs = [proc.communicate()[1] or "" for proc in procs]
        for argv, proc, err in zip(cmds, procs, errs):
            if proc.returncode != 0:
                raise ToolError(
                    f"external tool {argv[0]!r} failed (exit {proc.returncode}).\n"
                    f"argv: {argv}\n"
                    f"stderr tail:\n{err[-_STDERR_TAIL:]}",
                    returncode=proc.returncode,
                    argv=argv,
                )
    finally:
        if final_fh is not None:
            final_fh.close()
    return procs


def tool_version(bin_name: str | Path, version_arg: str = "--version") -> str | None:
    """Return a best-effort version string for ``bin_name`` (None if unresolved).

    Runs ``bin_name <version_arg>`` and returns the first non-empty output line
    (some tools print to stdout, some to stderr — both are inspected). Never
    raises: a missing binary or nonzero exit yields ``None``. Used by
    ``helixforge doctor`` (Phase 13).
    """
    try:
        proc = subprocess.run(
            [str(bin_name), str(version_arg)],
            check=False,
            capture_output=True,
            text=True,
            errors="replace",
        )
    except (OSError, ValueError, subprocess.SubprocessError):
        # Missing binary, bad argv, or any spawn failure → unresolved, not a crash.
        # ``errors="replace"`` keeps a tool that prints non-UTF-8 version banners
        # (e.g. a latin-1 © byte) from raising UnicodeDecodeError mid-decode.
        return None
    for stream in (proc.stdout, proc.stderr):
        for line in (stream or "").splitlines():
            line = line.strip()
            if line:
                return line
    return None
