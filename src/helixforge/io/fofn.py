"""Multi-file evidence-input expansion (repeatable + comma + FOFN)."""

from __future__ import annotations

import os
import warnings
from typing import Callable, Sequence


def expand_file_args(
    values: Sequence[str],
    list_file: str | os.PathLike[str] | None,
    *,
    label: str,
) -> list[str]:
    """Merge the three v1 multi-file forms into one path list.

    ``values`` are the raw direct-flag values (a click ``multiple=True`` tuple);
    each is split on ``,`` so the repeatable and comma forms collapse in one
    pass. ``list_file`` is the optional ``-list`` FOFN. Every resolved path is
    validated to exist — a missing path raises :class:`FileNotFoundError` naming
    ``label`` and, for FOFN entries, the list file + line number; nothing is
    silently dropped. The merge is order-preserving and de-duplicated by resolved
    absolute path (direct flags first, then FOFN entries).
    """
    merged: list[str] = []
    seen: set[str] = set()

    def _add(path: str) -> None:
        key = os.path.abspath(path)
        if key in seen:
            return
        seen.add(key)
        merged.append(path)

    # Forms 1 + 2: repeatable + comma-separated direct flags, in one pass.
    for value in values:
        for part in str(value).split(","):
            part = part.strip()
            if not part:
                continue
            if not os.path.exists(part):
                raise FileNotFoundError(f"{label} not found: {part}")
            _add(part)

    # Form 3: the -list FOFN companion (relative entries resolve against the
    # list file's own directory, not the CWD).
    if list_file:
        list_path = os.fspath(list_file)
        list_dir = os.path.dirname(os.path.abspath(list_path))
        with open(list_path) as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                resolved = line if os.path.isabs(line) else os.path.join(list_dir, line)
                if not os.path.exists(resolved):
                    raise FileNotFoundError(
                        f"{label} not found: {resolved} "
                        f"(from {list_path} line {lineno})"
                    )
                _add(resolved)

    return merged


_GTF_SUFFIXES = (".gtf", ".gff", ".gff3")


def _looks_like_gtf(path: str) -> bool:
    """True when ``path`` has a GTF/GFF extension (``.gz`` stripped first)."""
    lower = path.lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    return lower.endswith(_GTF_SUFFIXES)


def _first_real_line(path: str) -> str | None:
    """First non-blank, non-``#`` line of ``path`` (stripped), or ``None``."""
    try:
        with open(path) as fh:
            for raw in fh:
                line = raw.strip()
                if line and not line.startswith("#"):
                    return line
    except OSError:
        return None
    return None


def expand_stringtie_args(
    values: Sequence[str],
    list_file: str | os.PathLike[str] | None,
    *,
    label: str = "StringTie GTF",
    warn: Callable[[str], None] | None = None,
) -> list[str]:
    """:func:`expand_file_args` for ``--stringtie`` + a back-compat shim.

    ``--stringtie`` used to be a single sample-list file; it is now repeatable
    individual GTFs with the list-file moved to ``--stringtie-list`` (D3). The
    one-release safety net: a lone ``--stringtie`` value that is *itself* a FOFN
    — an existing non-GTF file whose first real line resolves to an existing GTF
    — is treated as ``--stringtie-list`` with a deprecation warning. ``warn`` is
    the message sink (defaults to :func:`warnings.warn` as a ``DeprecationWarning``).
    """
    values = list(values)
    if list_file is None and len(values) == 1 and "," not in values[0]:
        candidate = values[0]
        if os.path.isfile(candidate) and not _looks_like_gtf(candidate):
            first = _first_real_line(candidate)
            if first is not None:
                base = os.path.dirname(os.path.abspath(candidate))
                resolved = first if os.path.isabs(first) else os.path.join(base, first)
                if os.path.isfile(resolved) and _looks_like_gtf(resolved):
                    message = (
                        f"--stringtie was given a list-file ({candidate}); passing a "
                        f"sample-list to --stringtie is deprecated — use "
                        f"--stringtie-list instead."
                    )
                    if warn is not None:
                        warn(message)
                    else:
                        warnings.warn(message, DeprecationWarning, stacklevel=2)
                    return expand_file_args((), candidate, label=label)

    return expand_file_args(tuple(values), list_file, label=label)


__all__ = ["expand_file_args", "expand_stringtie_args"]
