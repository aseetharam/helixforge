"""Genomic region parsing and validation utilities.

This module provides utilities for parsing and validating genomic region
strings used in CLI arguments and task generation for parallel execution.

Coordinate conventions:
    - CLI input: 1-based inclusive (standard genomic convention)
    - Internal storage: 0-based half-open (Python convention)
    - GFF3 files: 1-based inclusive
    - BED files: 0-based half-open
    - HDF5 arrays: 0-based

Example:
    >>> from helixforge.utils.regions import parse_region, GenomicRegion
    >>> region = parse_region("chr1:1000-2000")
    >>> print(region.seqid)  # 'chr1'
    >>> print(region.start)  # 999 (0-based)
    >>> print(region.end)    # 2000 (half-open)
    >>> print(region.length) # 1001
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from helixforge.io.fasta import GenomeAccessor


class GenomicRegion(NamedTuple):
    """Parsed genomic region with 0-based half-open coordinates.

    Attributes:
        seqid: Scaffold/chromosome name.
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
    """

    seqid: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive

    def __str__(self) -> str:
        """Return string representation in 1-based inclusive format."""
        return f"{self.seqid}:{self.start + 1}-{self.end}"

    @property
    def length(self) -> int:
        """Get region length in base pairs."""
        return self.end - self.start

    def contains(self, seqid: str, position: int) -> bool:
        """Check if a position is contained within this region.

        Args:
            seqid: Scaffold name.
            position: Position to check (0-based).

        Returns:
            True if position is within region.
        """
        return self.seqid == seqid and self.start <= position < self.end

    def overlaps(self, other: GenomicRegion) -> bool:
        """Check if this region overlaps another.

        Args:
            other: Another GenomicRegion.

        Returns:
            True if regions overlap.
        """
        if self.seqid != other.seqid:
            return False
        return self.start < other.end and other.start < self.end

    def contains_region(self, other: GenomicRegion) -> bool:
        """Check if this region fully contains another.

        Args:
            other: Another GenomicRegion.

        Returns:
            True if this region fully contains the other.
        """
        if self.seqid != other.seqid:
            return False
        return self.start <= other.start and other.end <= self.end


# Regex pattern for region parsing
# Handles: chr1:1000-2000, chr1:1000..2000, scaffold_123:100-200
_REGION_PATTERN = re.compile(r"^(.+):(\d+)[-.]\.?(\d+)$")


def parse_region(region_str: str) -> GenomicRegion:
    """Parse region string into GenomicRegion.

    Supported formats:
        chr1:1000-2000      (1-based, inclusive - standard)
        chr1:1000..2000     (1-based, inclusive - GFF style)
        scaffold_1:0-1000   (handles underscores in names)

    The input coordinates are assumed to be 1-based inclusive
    (standard genomic convention). They are converted to 0-based
    half-open for internal storage.

    Args:
        region_str: Region string in format seqid:start-end.

    Returns:
        GenomicRegion with 0-based, half-open coordinates.

    Raises:
        ValueError: If format is invalid or coordinates are invalid.

    Example:
        >>> region = parse_region("chr1:1000-2000")
        >>> region.seqid
        'chr1'
        >>> region.start  # 0-based
        999
        >>> region.end    # half-open
        2000
    """
    match = _REGION_PATTERN.match(region_str)

    if not match:
        raise ValueError(
            f"Invalid region format: '{region_str}'. "
            "Expected format: seqid:start-end (e.g., chr1:1000-2000)"
        )

    seqid = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))

    # Validate 1-based input coordinates
    if start < 1:
        raise ValueError(f"Start position must be >= 1, got {start}")
    if end < start:
        raise ValueError(f"End must be >= start: {start}-{end}")

    # Convert to 0-based half-open
    # Input is 1-based inclusive (standard genomic convention)
    start_0based = start - 1
    end_0based = end  # half-open, so end stays same

    return GenomicRegion(seqid, start_0based, end_0based)


def validate_region(
    region: GenomicRegion,
    genome: GenomeAccessor,
) -> None:
    """Validate region against genome.

    Checks that:
    - The scaffold exists in the genome
    - The region coordinates are within scaffold bounds

    Args:
        region: GenomicRegion to validate.
        genome: GenomeAccessor for the reference genome.

    Raises:
        ValueError: If region is invalid for genome.
    """
    scaffold_lengths = genome.scaffold_lengths

    if region.seqid not in scaffold_lengths:
        available = list(scaffold_lengths.keys())[:5]
        suffix = "..." if len(scaffold_lengths) > 5 else ""
        raise ValueError(
            f"Scaffold '{region.seqid}' not found in genome. "
            f"Available: {available}{suffix}"
        )

    scaffold_len = scaffold_lengths[region.seqid]
    if region.end > scaffold_len:
        raise ValueError(
            f"Region end ({region.end}) exceeds scaffold length ({scaffold_len})"
        )


def region_to_str(region: GenomicRegion, one_based: bool = True) -> str:
    """Convert region to string representation.

    Args:
        region: GenomicRegion (0-based internally).
        one_based: If True, output 1-based inclusive coordinates.
            If False, output 0-based half-open coordinates.

    Returns:
        Region string in requested format.

    Example:
        >>> region = GenomicRegion("chr1", 999, 2000)
        >>> region_to_str(region, one_based=True)
        'chr1:1000-2000'
        >>> region_to_str(region, one_based=False)
        'chr1:999-2000'
    """
    if one_based:
        return f"{region.seqid}:{region.start + 1}-{region.end}"
    return f"{region.seqid}:{region.start}-{region.end}"


def region_from_scaffold(seqid: str, length: int) -> GenomicRegion:
    """Create a region covering an entire scaffold.

    Args:
        seqid: Scaffold name.
        length: Scaffold length in base pairs.

    Returns:
        GenomicRegion covering the entire scaffold.

    Example:
        >>> region = region_from_scaffold("chr1", 100000)
        >>> region.start
        0
        >>> region.end
        100000
    """
    return GenomicRegion(seqid, 0, length)


def format_region_for_cli(
    seqid: str,
    start: int,
    end: int,
    from_zero_based: bool = True,
) -> str:
    """Format coordinates as CLI region string.

    Args:
        seqid: Scaffold name.
        start: Start coordinate.
        end: End coordinate.
        from_zero_based: If True, input is 0-based half-open and will
            be converted to 1-based inclusive for CLI output.
            If False, input is already 1-based inclusive.

    Returns:
        Region string in 1-based inclusive format for CLI.

    Example:
        >>> format_region_for_cli("chr1", 999, 2000, from_zero_based=True)
        'chr1:1000-2000'
        >>> format_region_for_cli("chr1", 1000, 2000, from_zero_based=False)
        'chr1:1000-2000'
    """
    if from_zero_based:
        # Convert 0-based half-open to 1-based inclusive
        cli_start = start + 1
        cli_end = end
    else:
        cli_start = start
        cli_end = end

    return f"{seqid}:{cli_start}-{cli_end}"
