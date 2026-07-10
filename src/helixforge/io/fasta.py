"""Genome FASTA access."""

from __future__ import annotations

import os
from typing import Any, cast

from pyfaidx import Fasta

from helixforge.utils.sequences import reverse_complement


class GenomeAccessor:
    """Random-access reader over a genome FASTA, 0-based half-open."""

    def __init__(self, fasta_path: str | os.PathLike[str]) -> None:
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA not found: {fasta_path}")
        self.fasta_path = str(fasta_path)
        # as_raw=True -> plain strings; sequence_always_upper handled in get_sequence.
        self._fasta: Fasta | None = Fasta(
            self.fasta_path, as_raw=True, sequence_always_upper=True
        )

    # --- context manager ---
    def __enter__(self) -> GenomeAccessor:
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.close()

    def close(self) -> None:
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None

    # --- membership / metadata ---
    def __contains__(self, seqid: object) -> bool:
        return seqid in self._fasta  # type: ignore[operator]  # pyfaidx lacks stubs

    def get_seqids(self) -> list[str]:
        """Return all sequence ids, sorted."""
        return sorted(self._fasta.keys())  # type: ignore[union-attr]  # pyfaidx lacks stubs

    def get_length(self, seqid: str) -> int:
        if seqid not in self._fasta:  # type: ignore[operator]  # pyfaidx lacks stubs
            raise KeyError(f"unknown seqid: {seqid!r}")
        return len(self._fasta[seqid])  # type: ignore[index]  # pyfaidx lacks stubs

    def get_scaffold_lengths(self) -> dict[str, int]:
        """Return ``{seqid: length}`` for every scaffold."""
        return {sid: len(self._fasta[sid]) for sid in self._fasta.keys()}  # type: ignore[union-attr,index]  # pyfaidx lacks stubs

    # --- sequence extraction ---
    def get_sequence(self, seqid: str, start: int, end: int, strand: str = "+") -> str:
        """Return uppercase sequence for ``[start, end)``; RC if ``strand == '-'``.

        Coordinates are 0-based half-open and always low→high. Strand only
        controls reverse-complement, it never changes the coordinates.
        """
        if seqid not in self._fasta:  # type: ignore[operator]  # pyfaidx lacks stubs
            raise KeyError(f"unknown seqid: {seqid!r}")
        if strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {strand!r}")
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError(
                f"coordinates must be ints, got start={start!r} end={end!r}"
            )
        if start < 0:
            raise ValueError(f"start must be >= 0, got {start} ({seqid})")
        if end <= start:
            raise ValueError(
                f"require start < end, got start={start} end={end} ({seqid})"
            )
        length = len(self._fasta[seqid])  # type: ignore[index]  # pyfaidx lacks stubs
        if end > length:
            raise IndexError(f"end {end} exceeds length {length} of {seqid!r}")
        seq: str = cast(str, self._fasta[seqid][start:end])  # type: ignore[index]  # pyfaidx lacks stubs
        if strand == "-":
            return reverse_complement(seq)
        return seq
