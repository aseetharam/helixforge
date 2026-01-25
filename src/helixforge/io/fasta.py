"""FASTA file handling for genome sequences.

This module provides efficient access to genome sequences stored in
FASTA format, using pyfaidx for indexed random access.

Features:
    - Random access to sequences by region
    - Memory-efficient streaming for large genomes
    - Sequence extraction with strand awareness
    - Coordinate validation
    - Integration with HDF5 coordinate mapping

Example:
    >>> from helixforge.io.fasta import GenomeAccessor
    >>> genome = GenomeAccessor("genome.fa")
    >>> seq = genome.get_sequence("chr1", 1000, 2000)
    >>> print(seq[:50])
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import pyfaidx

if TYPE_CHECKING:
    pass

# =============================================================================
# Type Aliases
# =============================================================================

Strand = Literal["+", "-"]

logger = logging.getLogger(__name__)

# =============================================================================
# Complement Table
# =============================================================================

COMPLEMENT_TABLE = str.maketrans(
    "ACGTacgtNnRYSWKMBDHVryswkmbdhv",
    "TGCAtgcaNnYRSWMKVHDByrswmkvhdb",
)


# =============================================================================
# Utility Functions
# =============================================================================


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement of a DNA sequence.

    Handles IUPAC ambiguity codes and preserves case.

    Args:
        sequence: DNA sequence string.

    Returns:
        Reverse complement sequence.
    """
    return sequence.translate(COMPLEMENT_TABLE)[::-1]


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        GC content as a fraction (0.0 to 1.0).
    """
    if not sequence:
        return 0.0

    seq_upper = sequence.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    # Exclude N's from the calculation
    valid_bases = len(seq_upper) - seq_upper.count("N")

    if valid_bases == 0:
        return 0.0

    return gc_count / valid_bases


# =============================================================================
# Main Accessor Class
# =============================================================================


class GenomeAccessor:
    """Indexed FASTA access using pyfaidx.

    Provides efficient random access to genome sequences with support
    for strand-aware extraction and coordinate validation.

    Attributes:
        path: Path to the FASTA file.
        fai_path: Path to the .fai index file.

    Example:
        >>> genome = GenomeAccessor("genome.fa")
        >>> print(f"Scaffolds: {list(genome.scaffold_lengths.keys())[:5]}")
        >>> seq = genome.get_sequence("chr1", 1000, 2000)
        >>> rc_seq = genome.get_sequence("chr1", 1000, 2000, strand="-")
    """

    def __init__(self, fasta_path: Path | str) -> None:
        """Initialize the genome accessor.

        Args:
            fasta_path: Path to FASTA file. Will create .fai index if needed.

        Raises:
            FileNotFoundError: If FASTA file doesn't exist.
        """
        self.path = Path(fasta_path)
        if not self.path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.path}")

        self._fasta: pyfaidx.Fasta | None = None
        self._scaffold_lengths: dict[str, int] | None = None
        self._scaffold_order: list[str] | None = None

        self._open()

    def _open(self) -> None:
        """Open the FASTA file with pyfaidx."""
        # pyfaidx will create index if it doesn't exist
        self._fasta = pyfaidx.Fasta(
            str(self.path),
            sequence_always_upper=False,  # Preserve case for soft-masking
            read_ahead=10000,  # Read-ahead buffer for efficiency
            rebuild=False,  # Don't rebuild index unless necessary
        )

        # Cache scaffold information
        self._scaffold_order = list(self._fasta.keys())
        self._scaffold_lengths = {
            seqid: len(self._fasta[seqid]) for seqid in self._scaffold_order
        }

        logger.info(
            f"Opened FASTA: {self.path.name}, "
            f"{len(self._scaffold_order)} scaffolds, "
            f"{self.total_length:,} bp total"
        )

    @property
    def scaffold_lengths(self) -> dict[str, int]:
        """Return {seqid: length} mapping."""
        if self._scaffold_lengths is None:
            raise RuntimeError("FASTA file not opened")
        return self._scaffold_lengths.copy()

    @property
    def scaffold_order(self) -> list[str]:
        """Return list of scaffold names in file order."""
        if self._scaffold_order is None:
            raise RuntimeError("FASTA file not opened")
        return self._scaffold_order.copy()

    @property
    def total_length(self) -> int:
        """Total genome size in bases."""
        if self._scaffold_lengths is None:
            raise RuntimeError("FASTA file not opened")
        return sum(self._scaffold_lengths.values())

    @property
    def n_scaffolds(self) -> int:
        """Number of scaffolds."""
        if self._scaffold_order is None:
            raise RuntimeError("FASTA file not opened")
        return len(self._scaffold_order)

    def __enter__(self) -> GenomeAccessor:
        """Context manager entry."""
        return self

    def __exit__(self, *args: Any) -> None:
        """Context manager exit."""
        self.close()

    def close(self) -> None:
        """Close the FASTA file."""
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None

    def get_sequence(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: Strand = "+",
    ) -> str:
        """Get sequence for region (0-based, half-open coordinates).

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).
            strand: Strand (+ or -). Returns reverse complement if "-".

        Returns:
            Sequence string.

        Raises:
            KeyError: If seqid not in FASTA.
            ValueError: If coordinates are invalid.
        """
        if self._fasta is None:
            raise RuntimeError("FASTA file not opened")

        # Validate scaffold
        if seqid not in self._scaffold_lengths:
            raise KeyError(f"Unknown scaffold: {seqid}")

        # Validate coordinates
        scaffold_length = self._scaffold_lengths[seqid]
        if start < 0:
            raise ValueError(f"Start position cannot be negative: {start}")
        if end > scaffold_length:
            raise ValueError(
                f"End position {end} exceeds scaffold length {scaffold_length}"
            )
        if start >= end:
            raise ValueError(f"Start ({start}) must be less than end ({end})")

        # pyfaidx uses 0-based coordinates with slicing
        # Note: pyfaidx intervals are 0-based, inclusive on both ends for fetching
        # but we want 0-based half-open, so we adjust
        sequence = str(self._fasta[seqid][start:end])

        if strand == "-":
            sequence = reverse_complement(sequence)

        return sequence

    def get_length(self, seqid: str) -> int:
        """Get the length of a scaffold.

        Args:
            seqid: Scaffold name.

        Returns:
            Scaffold length in base pairs.

        Raises:
            KeyError: If seqid not in FASTA.
        """
        if self._scaffold_lengths is None:
            raise RuntimeError("FASTA file not opened")

        if seqid not in self._scaffold_lengths:
            raise KeyError(f"Unknown scaffold: {seqid}")

        return self._scaffold_lengths[seqid]

    def get_fai_path(self) -> Path:
        """Return path to .fai index file."""
        return self.path.with_suffix(self.path.suffix + ".fai")

    def validate_region(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> bool:
        """Check if region is valid.

        Args:
            seqid: Scaffold name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            True if region is valid, False otherwise.
        """
        if self._scaffold_lengths is None:
            return False

        if seqid not in self._scaffold_lengths:
            return False

        scaffold_length = self._scaffold_lengths[seqid]

        return 0 <= start < end <= scaffold_length

    def build_offset_index(self) -> dict[str, int]:
        """Build cumulative offset index for HDF5 coordinate mapping.

        Returns:
            {seqid: cumulative_start_position}

        Used to convert HDF5 linear coordinates to genomic coordinates.
        """
        if self._scaffold_order is None or self._scaffold_lengths is None:
            raise RuntimeError("FASTA file not opened")

        offsets = {}
        cumulative = 0

        for seqid in self._scaffold_order:
            offsets[seqid] = cumulative
            cumulative += self._scaffold_lengths[seqid]

        return offsets

    def get_gc_content(
        self,
        seqid: str,
        start: int | None = None,
        end: int | None = None,
    ) -> float:
        """Calculate GC content for a region or scaffold.

        Args:
            seqid: Scaffold name.
            start: Start position (0-based). If None, use 0.
            end: End position (0-based, exclusive). If None, use scaffold end.

        Returns:
            GC content as a fraction (0.0 to 1.0).
        """
        if self._scaffold_lengths is None:
            raise RuntimeError("FASTA file not opened")

        scaffold_length = self.get_length(seqid)

        if start is None:
            start = 0
        if end is None:
            end = scaffold_length

        sequence = self.get_sequence(seqid, start, end)
        return calculate_gc_content(sequence)

    def get_dinucleotide(
        self,
        seqid: str,
        position: int,
        strand: Strand = "+",
    ) -> str:
        """Get dinucleotide at a position (for splice site analysis).

        Args:
            seqid: Scaffold name.
            position: Position of first base (0-based).
            strand: Strand (+ or -).

        Returns:
            Two-nucleotide string.

        Raises:
            ValueError: If position is out of bounds.
        """
        scaffold_length = self.get_length(seqid)

        if position < 0 or position + 2 > scaffold_length:
            raise ValueError(
                f"Position {position} out of bounds for dinucleotide extraction"
            )

        return self.get_sequence(seqid, position, position + 2, strand)

    def iter_scaffolds(self) -> "Iterator[tuple[str, str]]":
        """Iterate over all scaffolds.

        Yields:
            Tuples of (seqid, sequence).

        Warning:
            May use significant memory for large scaffolds.
        """
        if self._scaffold_order is None:
            raise RuntimeError("FASTA file not opened")

        for seqid in self._scaffold_order:
            yield seqid, self.get_sequence(seqid, 0, self.get_length(seqid))

    def __contains__(self, seqid: str) -> bool:
        """Check if scaffold exists in FASTA."""
        if self._scaffold_lengths is None:
            return False
        return seqid in self._scaffold_lengths

    def __len__(self) -> int:
        """Return number of scaffolds."""
        return self.n_scaffolds


# Import Iterator for type hint
from typing import Iterator

# =============================================================================
# Legacy Compatibility
# =============================================================================

# Alias for backward compatibility with stub
FastaReader = GenomeAccessor


# =============================================================================
# Convenience Functions
# =============================================================================


def load_genome(fasta_path: Path | str) -> GenomeAccessor:
    """Convenience function to load a genome.

    Args:
        fasta_path: Path to FASTA file.

    Returns:
        GenomeAccessor instance.
    """
    return GenomeAccessor(fasta_path)


def extract_sequences_for_regions(
    genome: GenomeAccessor,
    regions: list[tuple[str, int, int, str]],
) -> list[str]:
    """Extract sequences for multiple regions.

    Args:
        genome: GenomeAccessor instance.
        regions: List of (seqid, start, end, strand) tuples.

    Returns:
        List of sequences in same order as regions.
    """
    return [
        genome.get_sequence(seqid, start, end, strand)
        for seqid, start, end, strand in regions
    ]
