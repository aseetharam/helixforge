"""Homology search interfaces for BLAST and DIAMOND.

This module provides wrappers for running protein homology searches
using BLAST and DIAMOND.

Features:
    - Run DIAMOND/BLAST searches with optimized parameters
    - Parse search results into structured objects
    - Support for both local and cluster execution
    - Result caching for repeated queries

Example:
    >>> from helixforge.homology.search import run_diamond
    >>> hits = run_diamond(
    ...     query="proteins.fa",
    ...     database="uniprot.dmnd",
    ...     evalue=1e-5,
    ... )

TODO:
    - Implement run_diamond function
    - Implement run_blast function
    - Add result parsing
    - Add database building utilities
    - Support for cluster submission
    - Add result caching
"""

from pathlib import Path
from typing import NamedTuple

# =============================================================================
# Data Structures
# =============================================================================


class HomologyHit(NamedTuple):
    """A single homology search hit.

    Attributes:
        query_id: Query sequence identifier.
        subject_id: Subject/target sequence identifier.
        percent_identity: Percent sequence identity.
        alignment_length: Length of the alignment.
        mismatches: Number of mismatches.
        gap_opens: Number of gap openings.
        query_start: Start position in query (1-based).
        query_end: End position in query (1-based).
        subject_start: Start position in subject (1-based).
        subject_end: End position in subject (1-based).
        evalue: E-value of the hit.
        bitscore: Bit score of the hit.
        query_coverage: Fraction of query covered.
        subject_coverage: Fraction of subject covered.
    """

    query_id: str
    subject_id: str
    percent_identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bitscore: float
    query_coverage: float
    subject_coverage: float


# =============================================================================
# Search Functions
# =============================================================================


def run_diamond(
    query: Path | str,
    database: Path | str,
    output: Path | str | None = None,
    evalue: float = 1e-5,
    max_target_seqs: int = 5,
    threads: int = 1,
    sensitive: bool = False,
    block_size: float = 2.0,
) -> list[HomologyHit]:
    """Run DIAMOND BLASTP search.

    Args:
        query: Path to query FASTA file.
        database: Path to DIAMOND database (.dmnd).
        output: Path for output file. If None, uses temp file.
        evalue: E-value threshold.
        max_target_seqs: Maximum hits per query.
        threads: Number of threads to use.
        sensitive: Use sensitive mode (slower but more sensitive).
        block_size: Memory block size in GB.

    Returns:
        List of HomologyHit objects.

    Raises:
        FileNotFoundError: If query or database doesn't exist.
        RuntimeError: If DIAMOND execution fails.
    """
    # TODO: Implement DIAMOND search
    raise NotImplementedError("run_diamond not yet implemented")


def run_blast(
    query: Path | str,
    database: Path | str,
    output: Path | str | None = None,
    program: str = "blastp",
    evalue: float = 1e-5,
    max_target_seqs: int = 5,
    threads: int = 1,
) -> list[HomologyHit]:
    """Run NCBI BLAST search.

    Args:
        query: Path to query FASTA file.
        database: Path to BLAST database.
        output: Path for output file. If None, uses temp file.
        program: BLAST program (blastp, blastn, blastx, tblastn).
        evalue: E-value threshold.
        max_target_seqs: Maximum hits per query.
        threads: Number of threads to use.

    Returns:
        List of HomologyHit objects.

    Raises:
        FileNotFoundError: If query or database doesn't exist.
        RuntimeError: If BLAST execution fails.
    """
    # TODO: Implement BLAST search
    raise NotImplementedError("run_blast not yet implemented")


def parse_diamond_output(output_path: Path | str) -> list[HomologyHit]:
    """Parse DIAMOND output file (format 6).

    Args:
        output_path: Path to DIAMOND output file.

    Returns:
        List of HomologyHit objects.
    """
    # TODO: Implement parsing
    raise NotImplementedError("parse_diamond_output not yet implemented")


def build_diamond_database(
    fasta: Path | str,
    output: Path | str,
    threads: int = 1,
) -> Path:
    """Build a DIAMOND database from a FASTA file.

    Args:
        fasta: Path to input FASTA file.
        output: Path for output database (without .dmnd extension).
        threads: Number of threads to use.

    Returns:
        Path to the created database.

    Raises:
        RuntimeError: If database building fails.
    """
    # TODO: Implement database building
    raise NotImplementedError("build_diamond_database not yet implemented")
