"""Utility functions for HelixForge.

This module provides common utilities used across HelixForge:

- Interval operations (overlap, merge, subtract)
- Sequence manipulation
- Logging configuration
- Genomic region parsing and validation

Example:
    >>> from helixforge.utils import intervals, sequences
    >>> overlaps = intervals.find_overlaps(region_a, regions_b)
    >>> from helixforge.utils.regions import parse_region, GenomicRegion
    >>> region = parse_region("chr1:1000-2000")
"""

from helixforge.utils.regions import (
    GenomicRegion,
    format_region_for_cli,
    parse_region,
    region_from_scaffold,
    region_to_str,
    validate_region,
)

# TODO: Import and expose main functions once implemented
# from helixforge.utils.intervals import find_overlaps, merge_intervals
# from helixforge.utils.sequences import reverse_complement, translate
# from helixforge.utils.logging import setup_logging

__all__ = [
    "GenomicRegion",
    "parse_region",
    "validate_region",
    "region_to_str",
    "region_from_scaffold",
    "format_region_for_cli",
]
