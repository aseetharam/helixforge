"""Confidence scoring for gene models using Helixer predictions.

This module implements multi-factor confidence scoring for gene predictions
by aggregating per-base class probabilities from Helixer HDF5 output.

The goal is to identify:
- High-confidence gene models (ready for use)
- Low-confidence models (need manual review or filtering)
- Specific problematic regions within genes (weak exons, uncertain boundaries)

Metrics implemented:
- mean_prob: Mean probability of called class across gene
- min_prob: Minimum probability (weakest point)
- median_prob: Median probability
- entropy: Shannon entropy (prediction uncertainty)
- boundary_sharpness: Sharpness at exon/intron boundaries
- coding_consistency: CDS frame consistency
- exon_min: Worst exon score

Example:
    >>> from helixforge.core.confidence import ConfidenceCalculator
    >>> from helixforge.io.hdf5 import HelixerHDF5Reader
    >>> from helixforge.io.fasta import GenomeAccessor
    >>> reader = HelixerHDF5Reader("predictions.h5", "genome.fa.fai")
    >>> genome = GenomeAccessor("genome.fa")
    >>> calc = ConfidenceCalculator(reader, genome)
    >>> for gene in genes:
    ...     score = calc.score_gene(gene)
    ...     print(f"{gene.gene_id}: {score.overall_score:.3f} ({score.confidence_class})")
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Iterator

import attrs
import numpy as np

if TYPE_CHECKING:
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GeneModel
    from helixforge.io.hdf5 import HelixerHDF5Reader

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

# Helixer class indices
CLASS_INTERGENIC = 0
CLASS_UTR = 1
CLASS_CDS = 2
CLASS_INTRON = 3

# Default thresholds
DEFAULT_LOW_CONF_THRESHOLD = 0.7
DEFAULT_HIGH_CONF_THRESHOLD = 0.85
DEFAULT_MEDIUM_CONF_THRESHOLD = 0.70
DEFAULT_WINDOW_SIZE = 20

# Weights for overall score calculation
DEFAULT_WEIGHTS = {
    "mean_prob": 0.30,
    "min_prob": 0.20,
    "boundary_sharpness": 0.20,
    "coding_consistency": 0.15,
    "worst_exon_score": 0.15,
}


# =============================================================================
# Enums
# =============================================================================


class ConfidenceMetric(Enum):
    """Available confidence metrics."""

    MEAN_PROB = "mean_prob"
    MIN_PROB = "min_prob"
    MEDIAN_PROB = "median_prob"
    ENTROPY = "entropy"
    BOUNDARY_SHARPNESS = "boundary"
    CODING_CONSISTENCY = "coding"
    EXON_MIN = "exon_min"


class ConfidenceClass(Enum):
    """Confidence classification levels."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class GeneConfidence:
    """Confidence scores for a single gene model.

    Attributes:
        gene_id: Gene identifier.
        seqid: Scaffold/chromosome name.
        start: Start position (0-based).
        end: End position (0-based, exclusive).
        strand: Strand (+ or -).
        mean_prob: Mean probability of called class.
        min_prob: Minimum probability (weakest point).
        median_prob: Median probability.
        entropy: Mean entropy across gene.
        boundary_sharpness: Sharpness at exon boundaries.
        coding_consistency: CDS frame consistency score.
        exon_scores: Mean probability per exon.
        intron_scores: Mean probability per intron.
        worst_exon_idx: Index of worst-scoring exon.
        worst_exon_score: Score of worst exon.
        low_confidence_regions: List of (start, end, score) tuples.
        confidence_class: "high", "medium", or "low".
        flags: List of specific issues identified.
    """

    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str

    # Aggregate metrics
    mean_prob: float
    min_prob: float
    median_prob: float
    entropy: float
    boundary_sharpness: float
    coding_consistency: float

    # Per-feature metrics
    exon_scores: list[float] = attrs.Factory(list)
    intron_scores: list[float] = attrs.Factory(list)
    worst_exon_idx: int = 0
    worst_exon_score: float = 0.0

    # Problematic regions (0-based genomic coordinates)
    low_confidence_regions: list[tuple[int, int, float]] = attrs.Factory(list)

    # Overall classification
    confidence_class: str = "medium"
    flags: list[str] = attrs.Factory(list)

    @property
    def overall_score(self) -> float:
        """Calculate weighted composite score (0-1).

        Weighting:
        - mean_prob: 0.3
        - min_prob: 0.2
        - boundary_sharpness: 0.2
        - coding_consistency: 0.15
        - worst_exon_score: 0.15
        """
        return (
            DEFAULT_WEIGHTS["mean_prob"] * self.mean_prob
            + DEFAULT_WEIGHTS["min_prob"] * self.min_prob
            + DEFAULT_WEIGHTS["boundary_sharpness"] * self.boundary_sharpness
            + DEFAULT_WEIGHTS["coding_consistency"] * self.coding_consistency
            + DEFAULT_WEIGHTS["worst_exon_score"] * self.worst_exon_score
        )

    @property
    def n_low_conf_regions(self) -> int:
        """Number of low-confidence regions."""
        return len(self.low_confidence_regions)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for DataFrame/TSV output."""
        return {
            "gene_id": self.gene_id,
            "seqid": self.seqid,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "mean_prob": self.mean_prob,
            "min_prob": self.min_prob,
            "median_prob": self.median_prob,
            "entropy": self.entropy,
            "boundary_sharpness": self.boundary_sharpness,
            "coding_consistency": self.coding_consistency,
            "worst_exon_score": self.worst_exon_score,
            "worst_exon_idx": self.worst_exon_idx,
            "overall_score": self.overall_score,
            "confidence_class": self.confidence_class,
            "flags": ",".join(self.flags) if self.flags else "",
            "n_low_conf_regions": self.n_low_conf_regions,
            "n_exons": len(self.exon_scores),
        }


@attrs.define(slots=True)
class RegionConfidence:
    """Confidence for a genomic region (for visualization).

    Attributes:
        seqid: Scaffold/chromosome name.
        start: Start position (0-based).
        end: End position (0-based, exclusive).
        per_base_entropy: Entropy at each position.
        per_base_max_prob: Maximum probability at each position.
        smoothed_confidence: Rolling window smoothed confidence.
    """

    seqid: str
    start: int
    end: int
    per_base_entropy: np.ndarray
    per_base_max_prob: np.ndarray
    smoothed_confidence: np.ndarray

    @property
    def length(self) -> int:
        """Length of the region."""
        return self.end - self.start

    @property
    def mean_entropy(self) -> float:
        """Mean entropy across region."""
        return float(np.mean(self.per_base_entropy))

    @property
    def mean_confidence(self) -> float:
        """Mean max probability across region."""
        return float(np.mean(self.per_base_max_prob))


# =============================================================================
# Calculator Class
# =============================================================================


class ConfidenceCalculator:
    """Calculate confidence metrics for gene models using Helixer HDF5 predictions.

    Thread-safe for parallel processing across genes.

    Attributes:
        hdf5_reader: Helixer HDF5 reader instance.
        genome: Genome accessor instance (optional, for splice site validation).
        low_conf_threshold: Threshold for low-confidence regions.
        window_size: Window size for smoothing and boundary detection.

    Example:
        >>> reader = HelixerHDF5Reader("predictions.h5", "genome.fa.fai")
        >>> genome = GenomeAccessor("genome.fa")
        >>> calc = ConfidenceCalculator(reader, genome)
        >>> score = calc.score_gene(gene)
    """

    def __init__(
        self,
        hdf5_reader: "HelixerHDF5Reader",
        genome: "GenomeAccessor | None" = None,
        low_conf_threshold: float = DEFAULT_LOW_CONF_THRESHOLD,
        window_size: int = DEFAULT_WINDOW_SIZE,
        weights: dict[str, float] | None = None,
    ) -> None:
        """Initialize the confidence calculator.

        Args:
            hdf5_reader: Helixer HDF5 reader for predictions.
            genome: Genome accessor for sequence (optional).
            low_conf_threshold: Threshold for low-confidence detection.
            window_size: Window size for smoothing and boundary detection.
            weights: Custom weights for overall score (optional).
        """
        self.hdf5_reader = hdf5_reader
        self.genome = genome
        self.low_conf_threshold = low_conf_threshold
        self.window_size = window_size
        self.weights = weights or DEFAULT_WEIGHTS.copy()

    def score_gene(self, gene: "GeneModel") -> GeneConfidence:
        """Calculate all confidence metrics for a gene.

        Steps:
        1. Extract predictions for gene region from HDF5
        2. Map exon/intron/CDS coordinates to array indices
        3. Calculate per-base metrics
        4. Aggregate to feature and gene level
        5. Identify low-confidence regions
        6. Assign confidence class and flags

        Args:
            gene: GeneModel to score.

        Returns:
            GeneConfidence with all calculated metrics.
        """
        # 1. Extract predictions for gene region (strand-aware if available)
        preds = self.hdf5_reader.get_predictions_for_region(
            gene.seqid, gene.start, gene.end, strand=gene.strand
        )

        # 2. Calculate per-base metrics
        max_probs = np.max(preds, axis=1)
        class_calls = np.argmax(preds, axis=1)
        entropy = self._calculate_entropy(preds)

        # 3. Aggregate metrics
        mean_prob = float(np.mean(max_probs))
        min_prob = float(np.min(max_probs))
        median_prob = float(np.median(max_probs))
        mean_entropy = float(np.mean(entropy))

        # 4. Get feature-level scores
        exon_scores = []
        intron_scores = []
        exon_boundaries = []

        # Get transcript features (use primary/first transcript)
        if gene.transcripts:
            transcript = gene.transcripts[0]

            # Score exons
            for exon_start, exon_end in transcript.exons:
                # Convert to array indices (relative to gene start)
                arr_start = max(0, exon_start - gene.start)
                arr_end = min(len(max_probs), exon_end - gene.start)

                if arr_end > arr_start:
                    exon_score = float(np.mean(max_probs[arr_start:arr_end]))
                    exon_scores.append(exon_score)
                    exon_boundaries.extend([arr_start, arr_end])

            # Score introns
            for intron in transcript.introns:
                intron_start, intron_end = intron
                arr_start = max(0, intron_start - gene.start)
                arr_end = min(len(max_probs), intron_end - gene.start)

                if arr_end > arr_start:
                    # For introns, check if INTRON class is predicted
                    intron_probs = preds[arr_start:arr_end, CLASS_INTRON]
                    intron_score = float(np.mean(intron_probs))
                    intron_scores.append(intron_score)

        # Calculate worst exon
        if exon_scores:
            worst_exon_idx = int(np.argmin(exon_scores))
            worst_exon_score = exon_scores[worst_exon_idx]
        else:
            worst_exon_idx = 0
            worst_exon_score = mean_prob

        # 5. Calculate boundary sharpness
        boundary_sharpness = self._calculate_boundary_sharpness(
            preds, exon_boundaries
        )

        # 6. Calculate coding consistency
        coding_consistency = self._calculate_coding_consistency(
            preds, gene, class_calls
        )

        # 7. Identify low-confidence regions
        low_conf_regions = self._identify_low_confidence_regions(
            max_probs, gene.start, self.low_conf_threshold
        )

        # 8. Classify confidence and generate flags
        # Create preliminary score object for classification
        preliminary = GeneConfidence(
            gene_id=gene.gene_id,
            seqid=gene.seqid,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            mean_prob=mean_prob,
            min_prob=min_prob,
            median_prob=median_prob,
            entropy=mean_entropy,
            boundary_sharpness=boundary_sharpness,
            coding_consistency=coding_consistency,
            exon_scores=exon_scores,
            intron_scores=intron_scores,
            worst_exon_idx=worst_exon_idx,
            worst_exon_score=worst_exon_score,
            low_confidence_regions=low_conf_regions,
        )

        confidence_class, flags = self._classify_confidence(preliminary)

        # Return final score with classification
        return GeneConfidence(
            gene_id=gene.gene_id,
            seqid=gene.seqid,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            mean_prob=mean_prob,
            min_prob=min_prob,
            median_prob=median_prob,
            entropy=mean_entropy,
            boundary_sharpness=boundary_sharpness,
            coding_consistency=coding_consistency,
            exon_scores=exon_scores,
            intron_scores=intron_scores,
            worst_exon_idx=worst_exon_idx,
            worst_exon_score=worst_exon_score,
            low_confidence_regions=low_conf_regions,
            confidence_class=confidence_class,
            flags=flags,
        )

    def score_genes_parallel(
        self,
        genes: Iterable["GeneModel"],
        n_workers: int = 1,
        chunk_size: int = 100,
    ) -> Iterator[GeneConfidence]:
        """Score multiple genes, optionally in parallel.

        Uses chunked processing to manage memory.

        Args:
            genes: Iterable of GeneModel objects.
            n_workers: Number of parallel workers.
            chunk_size: Genes per processing chunk.

        Yields:
            GeneConfidence for each gene.
        """
        genes_list = list(genes)

        if n_workers <= 1:
            # Serial processing
            for gene in genes_list:
                yield self.score_gene(gene)
        else:
            # Parallel processing
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                # Submit jobs in chunks
                for i in range(0, len(genes_list), chunk_size):
                    chunk = genes_list[i : i + chunk_size]
                    futures = {
                        executor.submit(self.score_gene, gene): gene
                        for gene in chunk
                    }

                    for future in as_completed(futures):
                        try:
                            yield future.result()
                        except Exception as e:
                            gene = futures[future]
                            logger.error(f"Error scoring gene {gene.gene_id}: {e}")
                            # Yield a low-confidence score on error
                            yield GeneConfidence(
                                gene_id=gene.gene_id,
                                seqid=gene.seqid,
                                start=gene.start,
                                end=gene.end,
                                strand=gene.strand,
                                mean_prob=0.0,
                                min_prob=0.0,
                                median_prob=0.0,
                                entropy=2.0,
                                boundary_sharpness=0.0,
                                coding_consistency=0.0,
                                confidence_class="low",
                                flags=["scoring_error"],
                            )

    def get_region_confidence(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str = "+",
    ) -> RegionConfidence:
        """Get per-base confidence for a region (for visualization).

        Args:
            seqid: Scaffold name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).
            strand: Strand ("+" or "-") for strand-aware prediction retrieval.

        Returns:
            RegionConfidence with per-base metrics.
        """
        preds = self.hdf5_reader.get_predictions_for_region(seqid, start, end, strand=strand)

        per_base_entropy = self._calculate_entropy(preds)
        per_base_max_prob = np.max(preds, axis=1)

        # Smooth with rolling window
        smoothed = self._smooth_array(per_base_max_prob, self.window_size)

        return RegionConfidence(
            seqid=seqid,
            start=start,
            end=end,
            per_base_entropy=per_base_entropy,
            per_base_max_prob=per_base_max_prob,
            smoothed_confidence=smoothed,
        )

    # =========================================================================
    # Internal Methods
    # =========================================================================

    def _calculate_entropy(self, probs: np.ndarray) -> np.ndarray:
        """Calculate Shannon entropy at each position.

        entropy = -sum(p * log2(p)) for each position
        High entropy = uncertain prediction

        Args:
            probs: Array of shape (positions, classes).

        Returns:
            Array of entropy values.
        """
        # Clip to avoid log(0)
        probs_safe = np.clip(probs, 1e-10, 1.0)

        # Shannon entropy: -sum(p * log2(p))
        entropy = -np.sum(probs_safe * np.log2(probs_safe), axis=1)

        return entropy

    def _calculate_boundary_sharpness(
        self,
        probs: np.ndarray,
        exon_boundaries: list[int],
    ) -> float:
        """Measure prediction sharpness at exon/intron boundaries.

        Sharp transitions (high confidence) vs gradual (uncertain).
        Uses gradient magnitude in window around boundaries.

        Args:
            probs: Prediction array (positions, classes).
            exon_boundaries: List of boundary positions (array indices).

        Returns:
            Boundary sharpness score (0-1).
        """
        if not exon_boundaries or len(probs) < self.window_size * 2:
            return 0.5  # Default neutral score

        # Get max prob at each position
        max_probs = np.max(probs, axis=1)

        sharpness_scores = []

        for boundary in exon_boundaries:
            # Skip boundaries too close to edges
            if boundary < self.window_size or boundary >= len(probs) - self.window_size:
                continue

            # Get window around boundary
            window_start = boundary - self.window_size // 2
            window_end = boundary + self.window_size // 2

            window = max_probs[window_start:window_end]

            # Calculate gradient magnitude
            gradient = np.abs(np.diff(window))
            max_gradient = np.max(gradient) if len(gradient) > 0 else 0

            # Sharp transition has high gradient at boundary
            sharpness_scores.append(min(1.0, max_gradient * 2))

        if sharpness_scores:
            return float(np.mean(sharpness_scores))
        return 0.5

    def _calculate_coding_consistency(
        self,
        probs: np.ndarray,
        gene: "GeneModel",
        class_calls: np.ndarray,
    ) -> float:
        """Check if CDS predictions maintain consistent pattern.

        Helixer doesn't explicitly model frame, but consistent high CDS
        probability across a CDS suggests reliable prediction.

        Args:
            probs: Prediction array.
            gene: Gene model.
            class_calls: Argmax class calls.

        Returns:
            Coding consistency score (0-1).
        """
        if not gene.transcripts:
            return 0.5

        transcript = gene.transcripts[0]

        if not transcript.cds:
            return 0.5  # No CDS to evaluate

        consistency_scores = []

        for cds_tuple in transcript.cds:
            cds_start, cds_end = cds_tuple[0], cds_tuple[1]

            # Convert to array indices
            arr_start = max(0, cds_start - gene.start)
            arr_end = min(len(probs), cds_end - gene.start)

            if arr_end <= arr_start:
                continue

            # Get CDS region predictions
            cds_probs = probs[arr_start:arr_end, CLASS_CDS]
            cds_calls = class_calls[arr_start:arr_end]

            # Consistency metrics:
            # 1. What fraction is called as CDS?
            cds_call_fraction = np.mean(cds_calls == CLASS_CDS)

            # 2. Mean CDS probability
            mean_cds_prob = np.mean(cds_probs)

            # 3. Variance (low variance = consistent)
            prob_variance = np.var(cds_probs)
            consistency = 1.0 - min(1.0, prob_variance * 4)

            # Combined score
            score = (cds_call_fraction + mean_cds_prob + consistency) / 3
            consistency_scores.append(score)

        if consistency_scores:
            return float(np.mean(consistency_scores))
        return 0.5

    def _identify_low_confidence_regions(
        self,
        max_probs: np.ndarray,
        gene_start: int,
        threshold: float,
        min_length: int = 10,
    ) -> list[tuple[int, int, float]]:
        """Find contiguous regions below confidence threshold.

        Args:
            max_probs: Array of max probabilities.
            gene_start: Gene start position (for coordinate conversion).
            threshold: Confidence threshold.
            min_length: Minimum region length to report.

        Returns:
            List of (genomic_start, genomic_end, mean_score) tuples.
        """
        regions = []

        # Find positions below threshold
        below_threshold = max_probs < threshold

        # Find contiguous runs
        in_region = False
        region_start = 0

        for i, below in enumerate(below_threshold):
            if below and not in_region:
                # Start of a low-confidence region
                in_region = True
                region_start = i
            elif not below and in_region:
                # End of a low-confidence region
                in_region = False
                if i - region_start >= min_length:
                    mean_score = float(np.mean(max_probs[region_start:i]))
                    regions.append(
                        (gene_start + region_start, gene_start + i, mean_score)
                    )

        # Handle region extending to end
        if in_region and len(max_probs) - region_start >= min_length:
            mean_score = float(np.mean(max_probs[region_start:]))
            regions.append(
                (gene_start + region_start, gene_start + len(max_probs), mean_score)
            )

        return regions

    def _classify_confidence(
        self,
        scores: GeneConfidence,
    ) -> tuple[str, list[str]]:
        """Assign confidence class and identify specific flags.

        Classes:
        - "high": overall_score >= 0.85, no critical flags
        - "medium": overall_score >= 0.70, minor flags ok
        - "low": overall_score < 0.70 or critical flags

        Flags:
        - "weak_exon": any exon score < 0.6
        - "uncertain_boundary": boundary_sharpness < 0.7
        - "high_entropy": mean entropy > 1.0
        - "low_min_prob": min_prob < 0.5
        - "short_low_conf": gene < 300bp with low confidence

        Args:
            scores: Preliminary GeneConfidence object.

        Returns:
            Tuple of (confidence_class, list of flags).
        """
        flags = []
        critical_flag = False

        # Check for weak exons
        if scores.exon_scores and min(scores.exon_scores) < 0.6:
            flags.append("weak_exon")
            critical_flag = True

        # Check boundary sharpness
        if scores.boundary_sharpness < 0.5:
            flags.append("uncertain_boundary")
            critical_flag = True
        elif scores.boundary_sharpness < 0.7:
            flags.append("uncertain_boundary")

        # Check entropy
        if scores.entropy > 1.5:
            flags.append("high_entropy")
            critical_flag = True
        elif scores.entropy > 1.0:
            flags.append("high_entropy")

        # Check minimum probability
        if scores.min_prob < 0.3:
            flags.append("very_low_min_prob")
            critical_flag = True
        elif scores.min_prob < 0.5:
            flags.append("low_min_prob")

        # Check for short low-confidence genes
        gene_length = scores.end - scores.start
        if gene_length < 300 and scores.mean_prob < 0.7:
            flags.append("short_low_conf")

        # Check for low-confidence regions
        if scores.n_low_conf_regions > 3:
            flags.append("many_low_conf_regions")
            critical_flag = True
        elif scores.n_low_conf_regions > 0:
            flags.append("has_low_conf_regions")

        # Determine class
        overall = scores.overall_score

        if critical_flag or overall < DEFAULT_MEDIUM_CONF_THRESHOLD:
            confidence_class = "low"
        elif overall >= DEFAULT_HIGH_CONF_THRESHOLD and not flags:
            confidence_class = "high"
        elif overall >= DEFAULT_HIGH_CONF_THRESHOLD:
            confidence_class = "medium"
        else:
            confidence_class = "medium"

        return confidence_class, flags

    def _smooth_array(
        self,
        arr: np.ndarray,
        window_size: int,
    ) -> np.ndarray:
        """Apply rolling mean smoothing.

        Args:
            arr: Input array.
            window_size: Smoothing window size.

        Returns:
            Smoothed array (same length as input).
        """
        if len(arr) < window_size:
            return arr

        # Use convolution for rolling mean
        kernel = np.ones(window_size) / window_size
        smoothed = np.convolve(arr, kernel, mode="same")

        return smoothed


# =============================================================================
# Output Writers
# =============================================================================


class ConfidenceWriter:
    """Write confidence scores to various formats."""

    @staticmethod
    def to_tsv(
        scores: Iterable[GeneConfidence],
        output_path: Path | str,
    ) -> None:
        """Write tab-separated confidence report.

        Columns:
        gene_id, seqid, start, end, strand, mean_prob, min_prob,
        median_prob, entropy, boundary_sharpness, coding_consistency,
        worst_exon_score, overall_score, confidence_class, flags,
        n_low_conf_regions

        Args:
            scores: Iterable of GeneConfidence objects.
            output_path: Output file path.
        """
        output_path = Path(output_path)

        columns = [
            "gene_id",
            "seqid",
            "start",
            "end",
            "strand",
            "mean_prob",
            "min_prob",
            "median_prob",
            "entropy",
            "boundary_sharpness",
            "coding_consistency",
            "worst_exon_score",
            "overall_score",
            "confidence_class",
            "flags",
            "n_low_conf_regions",
            "n_exons",
        ]

        with open(output_path, "w") as f:
            # Write header
            f.write("\t".join(columns) + "\n")

            # Write data
            for score in scores:
                data = score.to_dict()
                row = [str(data.get(col, "")) for col in columns]
                f.write("\t".join(row) + "\n")

        logger.info(f"Wrote confidence scores to {output_path}")

    @staticmethod
    def to_dataframe(
        scores: Iterable[GeneConfidence],
    ) -> "Any":  # Returns pd.DataFrame
        """Convert to pandas DataFrame.

        Args:
            scores: Iterable of GeneConfidence objects.

        Returns:
            pandas DataFrame with all confidence metrics.
        """
        import pandas as pd

        data = [s.to_dict() for s in scores]
        return pd.DataFrame(data)

    @staticmethod
    def to_bed(
        scores: Iterable[GeneConfidence],
        output_path: Path | str,
        score_column: str = "overall_score",
    ) -> None:
        """Write BED file with scores for genome browser visualization.

        Score field (0-1000) derived from specified metric.
        Color by confidence class (green/yellow/red).

        Args:
            scores: Iterable of GeneConfidence objects.
            output_path: Output file path.
            score_column: Metric to use for BED score field.
        """
        output_path = Path(output_path)

        # RGB colors for confidence classes
        colors = {
            "high": "0,150,0",  # Green
            "medium": "255,200,0",  # Yellow
            "low": "200,0,0",  # Red
        }

        with open(output_path, "w") as f:
            # BED header
            f.write('track name="HelixForge Confidence" itemRgb="On"\n')

            for score in scores:
                # Get the score value
                score_dict = score.to_dict()
                score_val = score_dict.get(score_column, score.overall_score)

                # Convert to 0-1000 range for BED
                bed_score = min(1000, int(score_val * 1000))

                # Get color
                rgb = colors.get(score.confidence_class, "128,128,128")

                # BED12 format: chrom, start, end, name, score, strand,
                #               thickStart, thickEnd, itemRgb
                f.write(
                    f"{score.seqid}\t{score.start}\t{score.end}\t"
                    f"{score.gene_id}\t{bed_score}\t{score.strand}\t"
                    f"{score.start}\t{score.end}\t{rgb}\n"
                )

        logger.info(f"Wrote BED file to {output_path}")

    @staticmethod
    def low_confidence_regions_bed(
        scores: Iterable[GeneConfidence],
        output_path: Path | str,
    ) -> None:
        """Write BED of low-confidence regions within genes.

        Args:
            scores: Iterable of GeneConfidence objects.
            output_path: Output file path.
        """
        output_path = Path(output_path)

        with open(output_path, "w") as f:
            f.write('track name="Low Confidence Regions" color="200,0,0"\n')

            for score in scores:
                for region_start, region_end, region_score in score.low_confidence_regions:
                    bed_score = int(region_score * 1000)
                    name = f"{score.gene_id}_lcr"

                    f.write(
                        f"{score.seqid}\t{region_start}\t{region_end}\t"
                        f"{name}\t{bed_score}\t{score.strand}\n"
                    )

        logger.info(f"Wrote low-confidence regions to {output_path}")


# =============================================================================
# Legacy Compatibility
# =============================================================================

# Aliases for backward compatibility with stub
ConfidenceScore = GeneConfidence
ConfidenceScorer = ConfidenceCalculator


# =============================================================================
# Utility Functions
# =============================================================================


def score_genes_from_files(
    h5_path: Path | str,
    gff_path: Path | str,
    fasta_path: Path | str,
    low_conf_threshold: float = DEFAULT_LOW_CONF_THRESHOLD,
    n_workers: int = 1,
) -> list[GeneConfidence]:
    """Convenience function to score genes from file paths.

    Args:
        h5_path: Path to Helixer HDF5 file.
        gff_path: Path to GFF3 file.
        fasta_path: Path to FASTA file.
        low_conf_threshold: Threshold for low-confidence regions.
        n_workers: Number of parallel workers.

    Returns:
        List of GeneConfidence objects.
    """
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser
    from helixforge.io.hdf5 import HelixerHDF5Reader

    # Get FAI path
    fasta_path = Path(fasta_path)
    fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")

    with HelixerHDF5Reader(h5_path, fai_path) as reader:
        with GenomeAccessor(fasta_path) as genome:
            parser = GFF3Parser(gff_path)
            genes = list(parser.iter_genes())

            calc = ConfidenceCalculator(
                reader, genome, low_conf_threshold=low_conf_threshold
            )

            scores = list(calc.score_genes_parallel(genes, n_workers=n_workers))

    return scores
