"""Genomic interval operations.

This module provides utilities for working with genomic intervals:

- Overlap detection
- Interval merging
- Interval subtraction
- Coverage calculation

Example:
    >>> from helixforge.utils.intervals import find_overlaps, merge_intervals
    >>> overlaps = find_overlaps(query, targets)
    >>> merged = merge_intervals(intervals)

TODO:
    - Implement interval operations
    - Add interval tree for efficient queries
    - Support for strand-aware operations
    - Add BED format support
"""

from typing import NamedTuple

# =============================================================================
# Data Structures
# =============================================================================


class Interval(NamedTuple):
    """A simple genomic interval.

    Attributes:
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
    """

    start: int
    end: int

    @property
    def length(self) -> int:
        """Get interval length."""
        return self.end - self.start

    def overlaps(self, other: "Interval") -> bool:
        """Check if this interval overlaps another."""
        return self.start < other.end and other.start < self.end

    def contains(self, position: int) -> bool:
        """Check if this interval contains a position."""
        return self.start <= position < self.end


class GenomicInterval(NamedTuple):
    """A genomic interval with chromosome and strand.

    Attributes:
        seqid: Chromosome/contig identifier.
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
        strand: Strand (+ or -).
    """

    seqid: str
    start: int
    end: int
    strand: str = "+"

    @property
    def length(self) -> int:
        """Get interval length."""
        return self.end - self.start

    def to_interval(self) -> Interval:
        """Convert to simple Interval."""
        return Interval(self.start, self.end)


# =============================================================================
# Overlap Operations
# =============================================================================


def overlaps(a: Interval, b: Interval) -> bool:
    """Check if two intervals overlap.

    Args:
        a: First interval.
        b: Second interval.

    Returns:
        True if intervals overlap.
    """
    return a.start < b.end and b.start < a.end


def overlap_length(a: Interval, b: Interval) -> int:
    """Calculate overlap length between two intervals.

    Args:
        a: First interval.
        b: Second interval.

    Returns:
        Overlap length (0 if no overlap).
    """
    if not overlaps(a, b):
        return 0
    return min(a.end, b.end) - max(a.start, b.start)


def find_overlaps(
    query: Interval,
    targets: list[Interval],
) -> list[tuple[int, Interval]]:
    """Find all intervals that overlap a query.

    Args:
        query: Query interval.
        targets: List of target intervals.

    Returns:
        List of (index, interval) tuples for overlapping intervals.
    """
    result = []
    for i, target in enumerate(targets):
        if overlaps(query, target):
            result.append((i, target))
    return result


# =============================================================================
# Merge Operations
# =============================================================================


def merge_intervals(intervals: list[Interval]) -> list[Interval]:
    """Merge overlapping intervals.

    Args:
        intervals: List of intervals to merge.

    Returns:
        List of merged intervals.
    """
    if not intervals:
        return []

    # Sort by start position
    sorted_intervals = sorted(intervals, key=lambda x: x.start)

    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        if current.start <= last.end:
            # Overlapping, extend the last interval
            merged[-1] = Interval(last.start, max(last.end, current.end))
        else:
            # Non-overlapping, add new interval
            merged.append(current)

    return merged


def merge_genomic_intervals(
    intervals: list[GenomicInterval],
    strand_aware: bool = False,
) -> list[GenomicInterval]:
    """Merge overlapping genomic intervals.

    Args:
        intervals: List of genomic intervals.
        strand_aware: If True, only merge same-strand intervals.

    Returns:
        List of merged genomic intervals.
    """
    # TODO: Implement genomic interval merging
    raise NotImplementedError("merge_genomic_intervals not yet implemented")


# =============================================================================
# Subtraction Operations
# =============================================================================


def subtract_intervals(
    intervals: list[Interval],
    to_remove: list[Interval],
) -> list[Interval]:
    """Subtract intervals from a set of intervals.

    Args:
        intervals: Base intervals.
        to_remove: Intervals to remove.

    Returns:
        Remaining intervals after subtraction.
    """
    # TODO: Implement interval subtraction
    raise NotImplementedError("subtract_intervals not yet implemented")


# =============================================================================
# Coverage Operations
# =============================================================================


def calculate_coverage(
    intervals: list[Interval],
    region: Interval,
) -> float:
    """Calculate fraction of a region covered by intervals.

    Args:
        intervals: Intervals providing coverage.
        region: Region to calculate coverage for.

    Returns:
        Coverage fraction (0.0 to 1.0).
    """
    if region.length == 0:
        return 0.0

    # Merge intervals and intersect with region
    merged = merge_intervals(intervals)
    covered = 0

    for interval in merged:
        # Intersect with region
        start = max(interval.start, region.start)
        end = min(interval.end, region.end)
        if start < end:
            covered += end - start

    return covered / region.length


def interval_gaps(
    intervals: list[Interval],
    region: Interval,
) -> list[Interval]:
    """Find gaps between intervals within a region.

    Args:
        intervals: Intervals to find gaps between.
        region: Bounding region.

    Returns:
        List of gap intervals.
    """
    # TODO: Implement gap finding
    raise NotImplementedError("interval_gaps not yet implemented")
