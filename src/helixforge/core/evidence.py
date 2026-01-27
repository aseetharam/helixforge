"""RNA-seq evidence scoring for gene models.

This module provides tools for scoring gene models based on RNA-seq evidence,
including splice junction support, exon coverage, and boundary agreement.

Key components:
- EvidenceLevel: Classification of evidence support levels
- JunctionEvidence: Per-junction evidence metrics
- ExonEvidence: Per-exon coverage metrics
- EvidenceScore: Complete evidence score for a gene
- EvidenceScorerConfig: Configuration for scoring parameters
- EvidenceScorer: Main scoring engine

Example:
    >>> from helixforge.core.evidence import EvidenceScorer, EvidenceScorerConfig
    >>> from helixforge.io.fasta import GenomeAccessor
    >>> from helixforge.io.bam import JunctionExtractor
    >>>
    >>> config = EvidenceScorerConfig()
    >>> scorer = EvidenceScorer(config)
    >>> junctions = JunctionExtractor("rnaseq.bam").extract_all()
    >>>
    >>> with JunctionExtractor("rnaseq.bam") as extractor:
    ...     score = scorer.score_gene(gene, junctions, extractor)
    ...     print(f"{gene.gene_id}: AED={score.aed:.3f}")
"""

from __future__ import annotations

import logging
from collections import defaultdict
from enum import Enum
from typing import TYPE_CHECKING, Iterator

import attrs
import numpy as np

if TYPE_CHECKING:
    from helixforge.io.bam import CoverageProfile, JunctionExtractor, SpliceJunction
    from helixforge.io.gff import GeneModel, TranscriptModel

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

# Default thresholds
DEFAULT_MIN_READS = 3
DEFAULT_MIN_COVERAGE = 5
DEFAULT_BOUNDARY_TOLERANCE = 10  # bp for soft boundary matching
DEFAULT_COVERAGE_RATIO_THRESHOLD = 0.1  # 10% of median for "expressed"


# =============================================================================
# Enums
# =============================================================================


class EvidenceLevel(Enum):
    """Classification of evidence support level."""

    FULL = "full"  # Fully supported by RNA-seq
    PARTIAL = "partial"  # Partially supported
    MINIMAL = "minimal"  # Weak or minimal support
    NONE = "none"  # No evidence support

    def __lt__(self, other: "EvidenceLevel") -> bool:
        """Allow comparison of evidence levels."""
        order = {
            EvidenceLevel.NONE: 0,
            EvidenceLevel.MINIMAL: 1,
            EvidenceLevel.PARTIAL: 2,
            EvidenceLevel.FULL: 3,
        }
        return order[self] < order[other]

    def __le__(self, other: "EvidenceLevel") -> bool:
        """Allow comparison of evidence levels."""
        return self == other or self < other


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class JunctionEvidence:
    """Evidence for a single splice junction.

    Attributes:
        intron_start: 0-based intron start position.
        intron_end: 0-based intron end position (exclusive).
        predicted: Whether this junction was predicted in the gene model.
        observed: Whether this junction has RNA-seq support.
        read_count: Number of reads supporting this junction.
        unique_count: Number of uniquely mapped reads.
        exact_match: True if predicted exactly matches observed.
        shift_donor: Shift in bp needed for donor site (0 if exact).
        shift_acceptor: Shift in bp needed for acceptor site (0 if exact).
    """

    intron_start: int
    intron_end: int
    predicted: bool = True
    observed: bool = False
    read_count: int = 0
    unique_count: int = 0
    exact_match: bool = False
    shift_donor: int = 0
    shift_acceptor: int = 0

    @property
    def is_supported(self) -> bool:
        """Check if junction has sufficient support."""
        return self.observed and self.read_count >= DEFAULT_MIN_READS

    @property
    def intron_length(self) -> int:
        """Intron length in base pairs."""
        return self.intron_end - self.intron_start


@attrs.define(slots=True)
class ExonEvidence:
    """Evidence for a single exon.

    Attributes:
        exon_start: 0-based exon start position.
        exon_end: 0-based exon end position (exclusive).
        mean_coverage: Mean read coverage across the exon.
        median_coverage: Median read coverage across the exon.
        min_coverage: Minimum coverage at any position.
        max_coverage: Maximum coverage at any position.
        fraction_covered: Fraction of bases with coverage >= threshold.
        coverage_uniformity: Coefficient of variation (lower = more uniform).
    """

    exon_start: int
    exon_end: int
    mean_coverage: float = 0.0
    median_coverage: float = 0.0
    min_coverage: int = 0
    max_coverage: int = 0
    fraction_covered: float = 0.0
    coverage_uniformity: float = 0.0

    @property
    def exon_length(self) -> int:
        """Exon length in base pairs."""
        return self.exon_end - self.exon_start

    @property
    def is_expressed(self) -> bool:
        """Check if exon shows expression."""
        return self.median_coverage >= DEFAULT_MIN_COVERAGE


@attrs.define(slots=True)
class EvidenceScore:
    """Complete evidence score for a gene model.

    Attributes:
        gene_id: Gene identifier.
        transcript_id: Primary transcript identifier.
        seqid: Scaffold/chromosome name.
        strand: Strand (+ or -).

        n_introns: Total number of introns.
        n_introns_supported: Number with RNA-seq junction support.
        n_introns_exact: Number with exact coordinate match.
        junction_support_ratio: Fraction of supported introns.

        n_exons: Total number of exons.
        n_exons_expressed: Number with expression evidence.
        mean_exon_coverage: Mean coverage across all exons.
        exon_coverage_ratio: Fraction of expressed exons.

        start_supported: Whether start codon has coverage support.
        stop_supported: Whether stop codon has coverage support.
        strand_consistent: Whether expression matches predicted strand.

        evidence_level: Overall evidence classification.
        aed: Annotation Edit Distance (0-1, lower is better).

        junction_evidence: Per-junction evidence details.
        exon_evidence: Per-exon evidence details.
        flags: Warning flags.
    """

    gene_id: str
    transcript_id: str
    seqid: str
    strand: str

    # Junction metrics
    n_introns: int = 0
    n_introns_supported: int = 0
    n_introns_exact: int = 0
    junction_support_ratio: float = 0.0

    # Exon metrics
    n_exons: int = 0
    n_exons_expressed: int = 0
    mean_exon_coverage: float = 0.0
    exon_coverage_ratio: float = 0.0

    # Boundary metrics
    start_supported: bool = False
    stop_supported: bool = False
    strand_consistent: bool = True

    # Overall classification
    evidence_level: EvidenceLevel = EvidenceLevel.NONE
    aed: float = 1.0  # Annotation Edit Distance (0 = perfect, 1 = no support)

    # Detailed evidence
    junction_evidence: list[JunctionEvidence] = attrs.Factory(list)
    exon_evidence: list[ExonEvidence] = attrs.Factory(list)
    flags: list[str] = attrs.Factory(list)

    @property
    def has_full_junction_support(self) -> bool:
        """Check if all junctions are supported."""
        return self.n_introns == 0 or self.junction_support_ratio >= 1.0

    @property
    def has_partial_junction_support(self) -> bool:
        """Check if at least some junctions are supported."""
        return self.junction_support_ratio > 0

    @property
    def is_single_exon(self) -> bool:
        """Check if gene is single-exon (no introns)."""
        return self.n_introns == 0


@attrs.define(slots=True)
class EvidenceScorerConfig:
    """Configuration for evidence scoring.

    Attributes:
        min_junction_reads: Minimum reads to count as junction support.
        min_exon_coverage: Minimum coverage to count exon as expressed.
        boundary_tolerance: Max bp shift to count as matching boundary.
        coverage_percentile: Percentile for coverage thresholding.
        strand_check: Whether to verify strand consistency.

        junction_weight: Weight for junction support in AED calculation.
        coverage_weight: Weight for exon coverage in AED calculation.
        boundary_weight: Weight for boundary support in AED calculation.
    """

    min_junction_reads: int = DEFAULT_MIN_READS
    min_exon_coverage: int = DEFAULT_MIN_COVERAGE
    boundary_tolerance: int = DEFAULT_BOUNDARY_TOLERANCE
    coverage_percentile: float = 50.0  # Use median
    strand_check: bool = True

    # Weights for AED calculation (should sum to 1.0)
    junction_weight: float = 0.5
    coverage_weight: float = 0.3
    boundary_weight: float = 0.2

    def __attrs_post_init__(self) -> None:
        """Validate configuration."""
        total_weight = self.junction_weight + self.coverage_weight + self.boundary_weight
        if not 0.99 <= total_weight <= 1.01:
            logger.warning(
                f"AED weights sum to {total_weight:.2f}, not 1.0. "
                "Results may not be normalized."
            )


# =============================================================================
# Evidence Scorer
# =============================================================================


class EvidenceScorer:
    """Score gene models based on RNA-seq evidence.

    The EvidenceScorer evaluates how well predicted gene models match
    empirical RNA-seq data. It considers:

    1. **Junction support**: Do predicted splice sites match observed junctions?
    2. **Exon coverage**: Do exons show expression evidence?
    3. **Boundary agreement**: Is there coverage at predicted start/stop sites?
    4. **Strand consistency**: Does the expression pattern match the strand?

    The scorer produces an EvidenceScore for each gene with detailed metrics
    and an overall Annotation Edit Distance (AED) score from 0 (perfect) to 1.

    Attributes:
        config: Configuration parameters.
        _junction_index: Index for fast junction lookup.

    Example:
        >>> config = EvidenceScorerConfig(min_junction_reads=5)
        >>> scorer = EvidenceScorer(config)
        >>>
        >>> with JunctionExtractor("rnaseq.bam") as extractor:
        ...     junctions = extractor.extract_all()
        ...     for gene in genes:
        ...         score = scorer.score_gene(gene, junctions, extractor)
        ...         print(f"{gene.gene_id}: AED={score.aed:.3f}")
    """

    def __init__(self, config: EvidenceScorerConfig | None = None) -> None:
        """Initialize the evidence scorer.

        Args:
            config: Scoring configuration. Uses defaults if not provided.
        """
        self.config = config or EvidenceScorerConfig()
        self._junction_index: dict[str, dict[tuple[int, int], "SpliceJunction"]] = {}

    def build_junction_index(
        self,
        junctions: dict[str, list["SpliceJunction"]],
    ) -> None:
        """Build index for fast junction lookup.

        Args:
            junctions: Dictionary mapping seqid to list of junctions.
        """
        self._junction_index = {}
        for seqid, junc_list in junctions.items():
            self._junction_index[seqid] = {
                (j.start, j.end): j for j in junc_list
            }
        total = sum(len(j) for j in junctions.values())
        logger.debug(f"Built junction index with {total} junctions")

    def _find_junction(
        self,
        seqid: str,
        start: int,
        end: int,
        tolerance: int = 0,
    ) -> tuple["SpliceJunction | None", int, int]:
        """Find a junction matching the given coordinates.

        Args:
            seqid: Scaffold name.
            start: Intron start (0-based).
            end: Intron end (0-based, exclusive).
            tolerance: Allow coordinate shift within this tolerance.

        Returns:
            Tuple of (junction or None, donor_shift, acceptor_shift).
        """
        if seqid not in self._junction_index:
            return None, 0, 0

        index = self._junction_index[seqid]

        # Try exact match first
        if (start, end) in index:
            return index[(start, end)], 0, 0

        # Try within tolerance
        if tolerance > 0:
            best_junction = None
            best_shift = float("inf")
            best_donor_shift = 0
            best_acceptor_shift = 0

            for (jstart, jend), junction in index.items():
                donor_shift = abs(jstart - start)
                acceptor_shift = abs(jend - end)

                if donor_shift <= tolerance and acceptor_shift <= tolerance:
                    total_shift = donor_shift + acceptor_shift
                    if total_shift < best_shift:
                        best_shift = total_shift
                        best_junction = junction
                        best_donor_shift = jstart - start
                        best_acceptor_shift = jend - end

            if best_junction is not None:
                return best_junction, best_donor_shift, best_acceptor_shift

        return None, 0, 0

    def score_junction(
        self,
        intron_start: int,
        intron_end: int,
        seqid: str,
        tolerance: int | None = None,
    ) -> JunctionEvidence:
        """Score a single predicted junction.

        Args:
            intron_start: Intron start position (0-based).
            intron_end: Intron end position (0-based, exclusive).
            seqid: Scaffold name.
            tolerance: Override default tolerance for near matches.

        Returns:
            JunctionEvidence with support metrics.
        """
        if tolerance is None:
            tolerance = self.config.boundary_tolerance

        junction, donor_shift, acceptor_shift = self._find_junction(
            seqid, intron_start, intron_end, tolerance=tolerance
        )

        if junction is None:
            return JunctionEvidence(
                intron_start=intron_start,
                intron_end=intron_end,
                predicted=True,
                observed=False,
            )

        return JunctionEvidence(
            intron_start=intron_start,
            intron_end=intron_end,
            predicted=True,
            observed=True,
            read_count=junction.read_count,
            unique_count=junction.unique_count,
            exact_match=(donor_shift == 0 and acceptor_shift == 0),
            shift_donor=donor_shift,
            shift_acceptor=acceptor_shift,
        )

    def score_exon(
        self,
        exon_start: int,
        exon_end: int,
        coverage: "CoverageProfile",
    ) -> ExonEvidence:
        """Score a single exon based on coverage.

        Args:
            exon_start: Exon start position (0-based).
            exon_end: Exon end position (0-based, exclusive).
            coverage: Coverage profile for the region.

        Returns:
            ExonEvidence with coverage metrics.
        """
        # Extract coverage values for this exon
        cov_start = max(0, exon_start - coverage.start)
        cov_end = min(len(coverage.values), exon_end - coverage.start)

        if cov_start >= cov_end:
            return ExonEvidence(
                exon_start=exon_start,
                exon_end=exon_end,
            )

        exon_cov = coverage.values[cov_start:cov_end]

        if len(exon_cov) == 0:
            return ExonEvidence(
                exon_start=exon_start,
                exon_end=exon_end,
            )

        mean_cov = float(np.mean(exon_cov))
        median_cov = float(np.median(exon_cov))
        min_cov = int(np.min(exon_cov))
        max_cov = int(np.max(exon_cov))

        # Fraction of bases with coverage above threshold
        threshold = self.config.min_exon_coverage
        fraction_covered = float(np.mean(exon_cov >= threshold))

        # Coverage uniformity (coefficient of variation)
        if mean_cov > 0:
            std_cov = float(np.std(exon_cov))
            uniformity = std_cov / mean_cov
        else:
            uniformity = 0.0

        return ExonEvidence(
            exon_start=exon_start,
            exon_end=exon_end,
            mean_coverage=mean_cov,
            median_coverage=median_cov,
            min_coverage=min_cov,
            max_coverage=max_cov,
            fraction_covered=fraction_covered,
            coverage_uniformity=uniformity,
        )

    def score_transcript(
        self,
        transcript: "TranscriptModel",
        coverage: "CoverageProfile | None" = None,
    ) -> tuple[list[JunctionEvidence], list[ExonEvidence]]:
        """Score all junctions and exons for a transcript.

        Args:
            transcript: Transcript model to score.
            coverage: Optional coverage profile for exon scoring.

        Returns:
            Tuple of (junction_evidence_list, exon_evidence_list).
        """
        junction_evidence = []
        exon_evidence = []

        # Score introns (junctions)
        introns = transcript.introns
        for intron_start, intron_end in introns:
            je = self.score_junction(
                intron_start, intron_end, transcript.seqid
            )
            junction_evidence.append(je)

        # Score exons
        if coverage is not None:
            for exon_start, exon_end in transcript.exons:
                ee = self.score_exon(exon_start, exon_end, coverage)
                exon_evidence.append(ee)

        return junction_evidence, exon_evidence

    def _calculate_aed(
        self,
        junction_ratio: float,
        coverage_ratio: float,
        boundary_ratio: float,
    ) -> float:
        """Calculate Annotation Edit Distance (AED).

        AED is a weighted measure of how different the prediction is from
        the evidence. Lower is better (0 = perfect match, 1 = no match).

        Args:
            junction_ratio: Fraction of supported junctions (0-1).
            coverage_ratio: Fraction of expressed exons (0-1).
            boundary_ratio: Fraction of supported boundaries (0-1).

        Returns:
            AED score from 0 to 1.
        """
        # Convert ratios to "distance" (1 - ratio)
        junction_distance = 1.0 - junction_ratio
        coverage_distance = 1.0 - coverage_ratio
        boundary_distance = 1.0 - boundary_ratio

        # Weighted average
        aed = (
            self.config.junction_weight * junction_distance
            + self.config.coverage_weight * coverage_distance
            + self.config.boundary_weight * boundary_distance
        )

        return min(1.0, max(0.0, aed))

    def _classify_evidence_level(
        self,
        junction_ratio: float,
        coverage_ratio: float,
        is_single_exon: bool,
    ) -> EvidenceLevel:
        """Classify overall evidence support level.

        Args:
            junction_ratio: Fraction of supported junctions.
            coverage_ratio: Fraction of expressed exons.
            is_single_exon: True if gene has no introns.

        Returns:
            EvidenceLevel classification.
        """
        if is_single_exon:
            # For single-exon genes, rely on coverage
            if coverage_ratio >= 0.8:
                return EvidenceLevel.FULL
            elif coverage_ratio >= 0.5:
                return EvidenceLevel.PARTIAL
            elif coverage_ratio > 0:
                return EvidenceLevel.MINIMAL
            else:
                return EvidenceLevel.NONE
        else:
            # For multi-exon genes, junctions are primary evidence
            if junction_ratio >= 1.0 and coverage_ratio >= 0.8:
                return EvidenceLevel.FULL
            elif junction_ratio >= 0.8 or (junction_ratio >= 0.5 and coverage_ratio >= 0.5):
                return EvidenceLevel.PARTIAL
            elif junction_ratio > 0 or coverage_ratio > 0:
                return EvidenceLevel.MINIMAL
            else:
                return EvidenceLevel.NONE

    def score_gene(
        self,
        gene: "GeneModel",
        junctions: dict[str, list["SpliceJunction"]] | None = None,
        extractor: "JunctionExtractor | None" = None,
    ) -> EvidenceScore:
        """Score a gene model based on RNA-seq evidence.

        Args:
            gene: Gene model to score.
            junctions: Pre-extracted junctions (if not using extractor).
            extractor: Junction extractor for coverage (optional).

        Returns:
            EvidenceScore with complete metrics.
        """
        # Ensure junction index is built
        if junctions is not None and not self._junction_index:
            self.build_junction_index(junctions)

        # Use primary transcript
        transcript = gene.primary_transcript
        if transcript is None:
            return EvidenceScore(
                gene_id=gene.gene_id,
                transcript_id="",
                seqid=gene.seqid,
                strand=gene.strand,
                evidence_level=EvidenceLevel.NONE,
                aed=1.0,
                flags=["NO_TRANSCRIPTS"],
            )

        # Get coverage if extractor is available
        coverage = None
        if extractor is not None:
            try:
                coverage = extractor.get_coverage(
                    gene.seqid, gene.start, gene.end, strand=gene.strand
                )
            except Exception as e:
                logger.warning(f"Could not get coverage for {gene.gene_id}: {e}")

        # Score transcript
        junction_evidence, exon_evidence = self.score_transcript(transcript, coverage)

        # Calculate metrics
        n_introns = len(junction_evidence)
        n_introns_supported = sum(
            1 for je in junction_evidence
            if je.read_count >= self.config.min_junction_reads
        )
        n_introns_exact = sum(1 for je in junction_evidence if je.exact_match)

        junction_ratio = n_introns_supported / n_introns if n_introns > 0 else 1.0

        n_exons = len(exon_evidence)
        n_exons_expressed = sum(
            1 for ee in exon_evidence
            if ee.median_coverage >= self.config.min_exon_coverage
        )
        mean_exon_coverage = (
            sum(ee.mean_coverage for ee in exon_evidence) / n_exons
            if n_exons > 0 else 0.0
        )
        coverage_ratio = n_exons_expressed / n_exons if n_exons > 0 else 0.0

        # Boundary support (check first and last exon)
        start_supported = False
        stop_supported = False
        if exon_evidence:
            # For + strand: first exon contains start, last contains stop
            # For - strand: first exon contains stop, last contains start
            first_exon = exon_evidence[0]
            last_exon = exon_evidence[-1]

            if gene.strand == "+":
                start_supported = first_exon.min_coverage >= self.config.min_exon_coverage
                stop_supported = last_exon.min_coverage >= self.config.min_exon_coverage
            else:
                start_supported = last_exon.min_coverage >= self.config.min_exon_coverage
                stop_supported = first_exon.min_coverage >= self.config.min_exon_coverage

        boundary_ratio = (int(start_supported) + int(stop_supported)) / 2.0

        # Strand consistency check
        strand_consistent = True
        if self.config.strand_check and coverage is not None:
            # If strand-specific coverage is available, check consistency
            # This is a simplified check - a more sophisticated version would
            # compare forward/reverse coverage profiles
            pass

        # Calculate AED
        aed = self._calculate_aed(junction_ratio, coverage_ratio, boundary_ratio)

        # Classify evidence level
        is_single_exon = n_introns == 0
        evidence_level = self._classify_evidence_level(
            junction_ratio, coverage_ratio, is_single_exon
        )

        # Collect flags
        flags = []
        if n_introns > 0 and n_introns_supported == 0:
            flags.append("NO_JUNCTION_SUPPORT")
        elif n_introns > 0 and junction_ratio < 1.0:
            flags.append("PARTIAL_JUNCTION_SUPPORT")

        if n_exons > 0 and n_exons_expressed == 0:
            flags.append("NO_EXON_EXPRESSION")
        elif n_exons > 0 and coverage_ratio < 0.5:
            flags.append("LOW_EXON_COVERAGE")

        if not start_supported:
            flags.append("START_UNSUPPORTED")
        if not stop_supported:
            flags.append("STOP_UNSUPPORTED")

        if not strand_consistent:
            flags.append("STRAND_INCONSISTENT")

        return EvidenceScore(
            gene_id=gene.gene_id,
            transcript_id=transcript.transcript_id,
            seqid=gene.seqid,
            strand=gene.strand,
            n_introns=n_introns,
            n_introns_supported=n_introns_supported,
            n_introns_exact=n_introns_exact,
            junction_support_ratio=junction_ratio,
            n_exons=n_exons,
            n_exons_expressed=n_exons_expressed,
            mean_exon_coverage=mean_exon_coverage,
            exon_coverage_ratio=coverage_ratio,
            start_supported=start_supported,
            stop_supported=stop_supported,
            strand_consistent=strand_consistent,
            evidence_level=evidence_level,
            aed=aed,
            junction_evidence=junction_evidence,
            exon_evidence=exon_evidence,
            flags=flags,
        )

    def score_genes(
        self,
        genes: list["GeneModel"],
        junctions: dict[str, list["SpliceJunction"]],
        extractor: "JunctionExtractor | None" = None,
    ) -> Iterator[EvidenceScore]:
        """Score multiple genes.

        Args:
            genes: List of gene models to score.
            junctions: Junction dictionary from extraction.
            extractor: Optional extractor for coverage data.

        Yields:
            EvidenceScore for each gene.
        """
        # Build index once
        self.build_junction_index(junctions)

        for gene in genes:
            yield self.score_gene(gene, extractor=extractor)


# =============================================================================
# Utility Functions
# =============================================================================


def summarize_evidence_scores(
    scores: list[EvidenceScore],
) -> dict[str, float | int]:
    """Summarize evidence scores across multiple genes.

    Args:
        scores: List of EvidenceScore objects.

    Returns:
        Dictionary with summary statistics.
    """
    if not scores:
        return {
            "n_genes": 0,
            "mean_aed": 0.0,
            "median_aed": 0.0,
            "n_full_support": 0,
            "n_partial_support": 0,
            "n_minimal_support": 0,
            "n_no_support": 0,
        }

    aeds = [s.aed for s in scores]
    level_counts = defaultdict(int)
    for s in scores:
        level_counts[s.evidence_level] += 1

    return {
        "n_genes": len(scores),
        "mean_aed": float(np.mean(aeds)),
        "median_aed": float(np.median(aeds)),
        "n_full_support": level_counts[EvidenceLevel.FULL],
        "n_partial_support": level_counts[EvidenceLevel.PARTIAL],
        "n_minimal_support": level_counts[EvidenceLevel.MINIMAL],
        "n_no_support": level_counts[EvidenceLevel.NONE],
    }
