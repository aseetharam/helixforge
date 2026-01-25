"""RNA-seq evidence collection for isoform analysis.

This module collects and organizes RNA-seq evidence for isoform
reconstruction, including:

- Splice junction extraction and filtering
- Coverage evidence for exons
- First/last exon evidence from read termini

Example:
    >>> from helixforge.isoforms.evidence import EvidenceCollector
    >>> collector = EvidenceCollector()
    >>> evidence = collector.collect(bam_files, region)

TODO:
    - Implement EvidenceCollector class
    - Add junction filtering and normalization
    - Collect coverage evidence
    - Support for strand-specific libraries
    - Handle multiple replicates/conditions
"""

from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

if TYPE_CHECKING:
    from helixforge.io.bam import BamReader, SpliceJunction

# =============================================================================
# Data Structures
# =============================================================================


class JunctionEvidence(NamedTuple):
    """Consolidated evidence for a splice junction.

    Attributes:
        seqid: Chromosome/contig identifier.
        donor: Donor site position (0-based).
        acceptor: Acceptor site position (0-based).
        strand: Inferred strand.
        total_reads: Total supporting reads.
        unique_reads: Uniquely mapped supporting reads.
        samples: Number of samples with support.
        is_canonical: Whether junction is canonical.
        confidence: Confidence score (0.0 to 1.0).
    """

    seqid: str
    donor: int
    acceptor: int
    strand: str
    total_reads: int
    unique_reads: int
    samples: int
    is_canonical: bool
    confidence: float


class ExonEvidence(NamedTuple):
    """Evidence supporting an exon.

    Attributes:
        seqid: Chromosome/contig identifier.
        start: Exon start (0-based).
        end: Exon end (0-based, exclusive).
        strand: Strand.
        mean_coverage: Mean coverage across exon.
        donor_support: Reads supporting donor site.
        acceptor_support: Reads supporting acceptor site.
        is_first: Evidence this is a first exon.
        is_last: Evidence this is a last exon.
    """

    seqid: str
    start: int
    end: int
    strand: str
    mean_coverage: float
    donor_support: int
    acceptor_support: int
    is_first: bool
    is_last: bool


class EvidenceCollection(NamedTuple):
    """Complete evidence for a genomic region.

    Attributes:
        seqid: Chromosome/contig identifier.
        start: Region start (0-based).
        end: Region end (0-based, exclusive).
        junctions: Splice junction evidence.
        exons: Exon evidence.
        coverage: Per-base coverage array.
    """

    seqid: str
    start: int
    end: int
    junctions: list[JunctionEvidence]
    exons: list[ExonEvidence]
    coverage: np.ndarray


# =============================================================================
# Collector Class
# =============================================================================


class EvidenceCollector:
    """Collects RNA-seq evidence for isoform reconstruction.

    Extracts and consolidates splice junctions, coverage, and
    terminal exon evidence from BAM files.

    Attributes:
        min_junction_reads: Minimum reads for a junction.
        min_coverage: Minimum coverage for an exon.
        strand_specific: Whether libraries are strand-specific.

    Example:
        >>> collector = EvidenceCollector(min_junction_reads=3)
        >>> evidence = collector.collect(bam_readers, "chr1", 0, 100000)
    """

    def __init__(
        self,
        min_junction_reads: int = 3,
        min_coverage: float = 1.0,
        strand_specific: bool = False,
    ) -> None:
        """Initialize the evidence collector.

        Args:
            min_junction_reads: Minimum reads for a junction.
            min_coverage: Minimum coverage for an exon.
            strand_specific: Whether libraries are strand-specific.
        """
        self.min_junction_reads = min_junction_reads
        self.min_coverage = min_coverage
        self.strand_specific = strand_specific

    def collect(
        self,
        bam_readers: list["BamReader"],
        seqid: str,
        start: int,
        end: int,
    ) -> EvidenceCollection:
        """Collect all evidence from a genomic region.

        Args:
            bam_readers: List of BAM readers (one per sample).
            seqid: Chromosome/contig identifier.
            start: Region start (0-based).
            end: Region end (0-based, exclusive).

        Returns:
            EvidenceCollection with all evidence types.
        """
        # TODO: Implement evidence collection
        raise NotImplementedError("collect not yet implemented")

    def collect_junctions(
        self,
        bam_readers: list["BamReader"],
        seqid: str,
        start: int,
        end: int,
    ) -> list[JunctionEvidence]:
        """Collect and consolidate splice junctions.

        Args:
            bam_readers: List of BAM readers.
            seqid: Chromosome/contig identifier.
            start: Region start (0-based).
            end: Region end (0-based, exclusive).

        Returns:
            List of consolidated JunctionEvidence.
        """
        # TODO: Implement junction collection
        raise NotImplementedError("collect_junctions not yet implemented")

    def collect_coverage(
        self,
        bam_readers: list["BamReader"],
        seqid: str,
        start: int,
        end: int,
    ) -> np.ndarray:
        """Collect merged coverage from all samples.

        Args:
            bam_readers: List of BAM readers.
            seqid: Chromosome/contig identifier.
            start: Region start (0-based).
            end: Region end (0-based, exclusive).

        Returns:
            Array of coverage values.
        """
        # TODO: Implement coverage collection
        raise NotImplementedError("collect_coverage not yet implemented")


# =============================================================================
# Utility Functions
# =============================================================================


def consolidate_junctions(
    junctions: list["SpliceJunction"],
    min_reads: int = 1,
) -> list[JunctionEvidence]:
    """Consolidate junctions from multiple samples.

    Args:
        junctions: Raw junctions from all samples.
        min_reads: Minimum total reads to keep.

    Returns:
        Consolidated junction evidence.
    """
    # TODO: Implement consolidation
    raise NotImplementedError("consolidate_junctions not yet implemented")


def filter_junctions(
    junctions: list[JunctionEvidence],
    min_reads: int = 3,
    require_canonical: bool = False,
    min_samples: int = 1,
) -> list[JunctionEvidence]:
    """Filter junctions by evidence criteria.

    Args:
        junctions: Junction evidence to filter.
        min_reads: Minimum total reads.
        require_canonical: Only keep canonical junctions.
        min_samples: Minimum samples with support.

    Returns:
        Filtered junction evidence.
    """
    # TODO: Implement filtering
    raise NotImplementedError("filter_junctions not yet implemented")
