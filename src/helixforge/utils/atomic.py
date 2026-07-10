"""Atomic file writes."""

from __future__ import annotations

import os
import stat
import tempfile
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any


def _target_mode(path: str | Path) -> int:
    """Permission bits the final file should carry.

    Reuse an existing target's mode on overwrite (so ``atomic_write`` is
    transparent); otherwise the umask-respecting default ``0o666 & ~umask`` —
    ``mkstemp``'s private ``0o600`` must not leak onto a freshly created file.
    """
    p = Path(path)
    if p.exists():
        return stat.S_IMODE(p.stat().st_mode)
    cur = os.umask(0)
    os.umask(cur)
    return 0o666 & ~cur


@contextmanager
def atomic_write(
    path: str | Path,
    mode: str = "w",
    *,
    encoding: str = "utf-8",
    newline: str | None = None,
) -> Generator[IO[Any], None, None]:
    """Context manager that writes ``path`` atomically.

    Yields a writable file handle backed by a temp file in ``path``'s directory.
    On clean exit the handle is flushed + ``fsync``'d and ``os.replace``'d onto
    ``path``. On **any** exception the temp file is removed and ``path`` is left
    exactly as it was (an existing target is untouched).

    Text mode defaults to UTF-8; binary mode (``"wb"``) ignores ``encoding``.
    Permissions of an existing target are preserved.
    """
    path = Path(path)
    directory = path.parent
    directory.mkdir(parents=True, exist_ok=True)
    binary = "b" in mode
    fd, tmp_name = tempfile.mkstemp(
        dir=str(directory), prefix=f".{path.name}.", suffix=".tmp"
    )
    tmp = Path(tmp_name)
    try:
        open_kw = {} if binary else {"encoding": encoding, "newline": newline}
        with os.fdopen(fd, mode, **open_kw) as fh:  # type: ignore[call-overload]  # dynamic mode + kwargs
            yield fh
            fh.flush()
            os.fsync(fh.fileno())
        os.chmod(tmp, _target_mode(path))
        os.replace(tmp, path)
    except BaseException:
        # Leave the destination intact; drop the partial temp file.
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass
        raise
