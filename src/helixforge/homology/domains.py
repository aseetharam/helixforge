"""Protein domain annotation using InterProScan.

This module provides tools for annotating predicted proteins with
protein domains and functional annotations.

Features:
    - Run InterProScan on predicted proteins
    - Parse InterProScan results
    - Map domains to genomic coordinates
    - Use domains for gene model validation

Example:
    >>> from helixforge.homology.domains import run_interproscan
    >>> domains = run_interproscan("proteins.fa")
    >>> for domain in domains:
    ...     print(f"{domain.protein_id}: {domain.name}")

TODO:
    - Implement InterProScan wrapper
    - Add result parsing
    - Map domains to genomic coordinates
    - Support for pre-computed domain databases
    - Integrate with gene validation
"""

from pathlib import Path
from typing import NamedTuple

# =============================================================================
# Data Structures
# =============================================================================


class ProteinDomain(NamedTuple):
    """A protein domain annotation.

    Attributes:
        protein_id: Identifier of the protein.
        source: Source database (Pfam, SMART, etc.).
        domain_id: Domain accession.
        name: Domain name.
        description: Domain description.
        start: Start position in protein (1-based).
        end: End position in protein (1-based).
        score: Domain score/evalue.
        interpro_id: InterPro accession (if available).
        go_terms: Associated GO terms.
    """

    protein_id: str
    source: str
    domain_id: str
    name: str
    description: str
    start: int
    end: int
    score: float
    interpro_id: str | None
    go_terms: list[str]


class GenomicDomain(NamedTuple):
    """A domain mapped to genomic coordinates.

    Attributes:
        domain: The original protein domain.
        seqid: Chromosome/contig identifier.
        genomic_start: Start position in genome (0-based).
        genomic_end: End position in genome (0-based).
        strand: Strand (+ or -).
        exons: Exon coordinates containing this domain.
    """

    domain: ProteinDomain
    seqid: str
    genomic_start: int
    genomic_end: int
    strand: str
    exons: list[tuple[int, int]]


# =============================================================================
# InterProScan Functions
# =============================================================================


def run_interproscan(
    fasta: Path | str,
    output_dir: Path | str | None = None,
    applications: list[str] | None = None,
    cpu: int = 1,
    disable_precalc: bool = False,
) -> list[ProteinDomain]:
    """Run InterProScan on a protein FASTA file.

    Args:
        fasta: Path to protein FASTA file.
        output_dir: Directory for output files.
        applications: List of applications to run (e.g., ["Pfam", "SMART"]).
                     If None, runs all available.
        cpu: Number of CPUs to use.
        disable_precalc: Disable pre-calculated match lookup.

    Returns:
        List of ProteinDomain annotations.

    Raises:
        FileNotFoundError: If FASTA file doesn't exist.
        RuntimeError: If InterProScan execution fails.
    """
    # TODO: Implement InterProScan wrapper
    raise NotImplementedError("run_interproscan not yet implemented")


def parse_interproscan_tsv(path: Path | str) -> list[ProteinDomain]:
    """Parse InterProScan TSV output.

    Args:
        path: Path to InterProScan TSV output file.

    Returns:
        List of ProteinDomain annotations.
    """
    # TODO: Implement parsing
    raise NotImplementedError("parse_interproscan_tsv not yet implemented")


def parse_interproscan_gff(path: Path | str) -> list[ProteinDomain]:
    """Parse InterProScan GFF3 output.

    Args:
        path: Path to InterProScan GFF3 output file.

    Returns:
        List of ProteinDomain annotations.
    """
    # TODO: Implement parsing
    raise NotImplementedError("parse_interproscan_gff not yet implemented")


# =============================================================================
# Mapping Functions
# =============================================================================


def map_domains_to_genome(
    domains: list[ProteinDomain],
    transcript_mapping: dict[str, list[tuple[int, int, int]]],
) -> list[GenomicDomain]:
    """Map protein domains to genomic coordinates.

    Args:
        domains: List of protein domains.
        transcript_mapping: Dict mapping protein_id to list of
                           (genomic_start, genomic_end, phase) for each codon.

    Returns:
        List of GenomicDomain objects.
    """
    # TODO: Implement coordinate mapping
    raise NotImplementedError("map_domains_to_genome not yet implemented")


# =============================================================================
# Validation Functions
# =============================================================================


def validate_domains_consistency(
    gene_domains: list[ProteinDomain],
    expected_domains: list[str] | None = None,
) -> tuple[bool, list[str]]:
    """Check if domains are consistent with expectations.

    Args:
        gene_domains: Domains found in the gene.
        expected_domains: Optional list of expected domain IDs.

    Returns:
        Tuple of (is_consistent, list of issues).
    """
    # TODO: Implement validation
    raise NotImplementedError("validate_domains_consistency not yet implemented")
