"""Main refinement pipeline combining splice correction, boundary adjustment,
confidence scoring, and evidence scoring.

This module provides the primary HelixForge pipeline for refining Helixer
gene predictions using RNA-seq evidence.

Key components:
- RefineConfig: Configuration for the refinement pipeline
- RefinedGene: Result container for a refined gene
- RefinePipeline: Main orchestration class
- RefineReportWriter: Report generation utilities

Example:
    >>> from helixforge.core.refine import RefinePipeline, RefineConfig
    >>> from helixforge.io.fasta import GenomeAccessor
    >>> from helixforge.io.hdf5 import HelixerHDF5Reader
    >>>
    >>> config = RefineConfig()
    >>> pipeline = RefinePipeline(
    ...     genome=genome,
    ...     hdf5_reader=hdf5_reader,
    ...     bam_files=[Path("rnaseq.bam")],
    ...     config=config,
    ... )
    >>> for result in pipeline.refine_genes(genes):
    ...     print(f"{result.gene.gene_id}: conf={result.confidence_score:.3f}")
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator

import attrs

if TYPE_CHECKING:
    from helixforge.io.bam import JunctionExtractor, SpliceJunction
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GeneModel
    from helixforge.io.hdf5 import HelixerHDF5Reader
    from helixforge.core.confidence import GeneConfidence
    from helixforge.core.evidence import EvidenceScore
    from helixforge.core.splice import GeneSpliceReport

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

DEFAULT_MAX_SHIFT = 15
DEFAULT_MIN_JUNCTION_READS = 3
DEFAULT_MIN_TISSUES = 1
DEFAULT_BOUNDARY_WINDOW = 30
DEFAULT_CONFIDENCE_THRESHOLD = 0.7
DEFAULT_VERY_LOW_THRESHOLD = 0.5
DEFAULT_COVERAGE_NORMALIZATION = 50.0


# =============================================================================
# Configuration
# =============================================================================


@attrs.define(slots=True)
class RefineConfig:
    """Configuration for the refine pipeline.

    Attributes:
        max_shift: Maximum splice site correction distance in bp.
        min_junction_reads: Minimum reads to consider junction supported.
        min_tissues: Minimum samples supporting a junction (for multi-BAM).
        adjust_boundaries: Whether to adjust start/stop codon boundaries.
        boundary_search_window: Window size for boundary search.
        confidence_threshold: Threshold for low confidence flagging.
        confidence_very_low_threshold: Threshold for very low confidence.
        coverage_normalization: Expected median coverage for evidence scoring.
        stranded: Whether RNA-seq is strand-specific.
        include_introns: Add intron features to output GFF3.
        include_exon_coverage: Add per-exon coverage attributes.
    """

    # Splice refinement
    max_shift: int = DEFAULT_MAX_SHIFT
    min_junction_reads: int = DEFAULT_MIN_JUNCTION_READS
    min_tissues: int = DEFAULT_MIN_TISSUES

    # Boundary adjustment
    adjust_boundaries: bool = True
    boundary_search_window: int = DEFAULT_BOUNDARY_WINDOW

    # Confidence scoring
    confidence_threshold: float = DEFAULT_CONFIDENCE_THRESHOLD
    confidence_very_low_threshold: float = DEFAULT_VERY_LOW_THRESHOLD

    # Evidence scoring
    coverage_normalization: float = DEFAULT_COVERAGE_NORMALIZATION
    stranded: bool = False
    min_exon_coverage: int = 5
    boundary_tolerance: int = 10
    skip_coverage: bool = False

    # Output options
    include_introns: bool = False
    include_exon_coverage: bool = False


# =============================================================================
# Result Container
# =============================================================================


@attrs.define(slots=True)
class RefinedGene:
    """Result of refining a single gene.

    Attributes:
        gene: The refined GeneModel.
        original_gene: Original GeneModel before refinement.
        confidence: HDF5-based confidence scores.
        evidence: RNA-seq evidence scores.
        splice_report: Splice refinement report.
        splice_corrections: Number of splice site corrections made.
        boundary_adjusted: Whether boundaries were adjusted.
        flags: Combined quality flags.
    """

    gene: "GeneModel"
    original_gene: "GeneModel"

    # Scores
    confidence: "GeneConfidence | None" = None
    evidence: "EvidenceScore | None" = None
    splice_report: "GeneSpliceReport | None" = None

    # Refinement details
    splice_corrections: int = 0
    boundary_adjusted: bool = False

    # Combined flags
    flags: list[str] = attrs.Factory(list)

    @property
    def confidence_score(self) -> float | None:
        """Get HDF5-based confidence score."""
        return self.confidence.overall_score if self.confidence else None

    @property
    def evidence_score(self) -> float | None:
        """Get RNA-seq evidence score (structural score)."""
        if self.evidence is None:
            return None
        # Use (1 - AED) as the evidence score for consistency
        return 1.0 - self.evidence.aed

    @property
    def aed(self) -> float | None:
        """Get Annotation Edit Distance (MAKER-compatible)."""
        return self.evidence.aed if self.evidence else None

    @property
    def junction_support_str(self) -> str:
        """Get junction support as fraction string."""
        if self.evidence is None:
            return "N/A"
        if self.evidence.n_introns == 0:
            return "N/A"
        return f"{self.evidence.n_introns_supported}/{self.evidence.n_introns}"

    def to_report_dict(self) -> dict[str, Any]:
        """Convert to dictionary for report TSV."""
        n_exons = 0
        if self.gene.transcripts:
            n_exons = len(self.gene.transcripts[0].exons)

        return {
            "gene_id": self.gene.gene_id,
            "seqid": self.gene.seqid,
            "start": self.gene.start,
            "end": self.gene.end,
            "strand": self.gene.strand,
            "n_exons": n_exons,
            "confidence_score": round(self.confidence_score, 4) if self.confidence_score else None,
            "evidence_score": round(self.evidence_score, 4) if self.evidence_score else None,
            "aed": round(self.aed, 4) if self.aed else None,
            "junction_support": self.junction_support_str,
            "mean_coverage": round(self.evidence.mean_exon_coverage, 1) if self.evidence else None,
            "splice_corrections": self.splice_corrections,
            "boundary_adjusted": self.boundary_adjusted,
            "flags": ",".join(self.flags) if self.flags else "",
        }


# =============================================================================
# Main Pipeline
# =============================================================================


class RefinePipeline:
    """Main refinement pipeline.

    Orchestrates:
    1. Splice site correction using RNA-seq junctions and PWM scoring
    2. Boundary adjustment for start/stop codons
    3. Confidence scoring from Helixer HDF5 predictions
    4. Evidence scoring from RNA-seq coverage
    5. Flag aggregation

    Attributes:
        genome: Genome accessor.
        hdf5_reader: Helixer HDF5 reader.
        config: Pipeline configuration.

    Example:
        >>> pipeline = RefinePipeline(
        ...     genome=genome,
        ...     hdf5_reader=hdf5_reader,
        ...     bam_files=[Path("rnaseq.bam")],
        ... )
        >>> for result in pipeline.refine_genes(genes):
        ...     print(f"{result.gene.gene_id}: {result.confidence_score:.3f}")
    """

    def __init__(
        self,
        genome: "GenomeAccessor",
        hdf5_reader: "HelixerHDF5Reader",
        bam_files: list[Path] | None = None,
        junction_files: list[Path] | None = None,
        config: RefineConfig | None = None,
        verbose: bool = False,
        quiet: bool = False,
    ) -> None:
        """Initialize pipeline.

        Args:
            genome: Genome accessor.
            hdf5_reader: Helixer HDF5 reader for confidence scoring.
            bam_files: RNA-seq BAM files (sorted, indexed).
            junction_files: Pre-extracted junction files.
            config: Pipeline configuration.
            verbose: Show detailed progress messages.
            quiet: Suppress all progress messages.

        Raises:
            ValueError: If no RNA-seq evidence is provided.
        """
        self.genome = genome
        self.hdf5_reader = hdf5_reader
        self.bam_files = bam_files or []
        self.junction_files = junction_files or []
        self.config = config or RefineConfig()
        self.verbose = verbose
        self.quiet = quiet

        if not self.bam_files and not self.junction_files:
            raise ValueError(
                "RNA-seq evidence required: provide bam_files or junction_files"
            )

        # Components initialized lazily
        self._junctions: dict[str, list["SpliceJunction"]] | None = None
        self._splice_refiner = None
        self._boundary_adjuster = None
        self._confidence_calculator = None
        self._evidence_scorer = None
        self._primary_extractor: "JunctionExtractor | None" = None

        self._init_components()

    def _log(self, message: str, indent: int = 0) -> None:
        """Log a progress message if verbose mode is enabled.

        Args:
            message: Message to display.
            indent: Number of spaces to indent.
        """
        if self.verbose and not self.quiet:
            prefix = "  " * indent
            print(f"{prefix}{message}")

    def _init_components(self) -> None:
        """Initialize pipeline components."""
        from helixforge.core.splice import SpliceRefiner, PositionWeightMatrix
        from helixforge.core.boundaries import BoundaryAdjuster
        from helixforge.core.confidence import ConfidenceCalculator
        from helixforge.core.evidence import EvidenceScorer, EvidenceScorerConfig

        # Load junctions
        self._junctions = self._load_junctions()

        # Splice refiner
        self._log("Initializing splice refiner...", indent=1)
        try:
            self._log("Loading PWM models...", indent=2)
            donor_pwm, acceptor_pwm = PositionWeightMatrix.load_plant_defaults()
        except Exception as e:
            logger.warning(f"Could not load default PWMs: {e}. Using PWM-free mode.")
            donor_pwm, acceptor_pwm = None, None

        self._splice_refiner = SpliceRefiner(
            genome=self.genome,
            junctions=self._junctions,
            donor_pwm=donor_pwm,
            acceptor_pwm=acceptor_pwm,
            max_shift=self.config.max_shift,
            min_junction_reads=self.config.min_junction_reads,
        )

        # Boundary adjuster
        if self.config.adjust_boundaries:
            self._log("Initializing boundary adjuster...", indent=1)
            self._boundary_adjuster = BoundaryAdjuster(genome=self.genome)

        # Confidence calculator
        self._log("Initializing confidence calculator...", indent=1)
        self._confidence_calculator = ConfidenceCalculator(
            hdf5_reader=self.hdf5_reader,
            genome=self.genome,
            low_conf_threshold=self.config.confidence_threshold,
        )

        # Evidence scorer
        self._log("Initializing evidence scorer...", indent=1)
        evidence_config = EvidenceScorerConfig(
            min_junction_reads=self.config.min_junction_reads,
            strand_check=self.config.stranded,
        )
        self._evidence_scorer = EvidenceScorer(evidence_config)
        self._log("Building junction index...", indent=2)
        self._evidence_scorer.build_junction_index(self._junctions)

    def _load_junctions(self) -> dict[str, list["SpliceJunction"]]:
        """Load junctions from BAM or junction files."""
        from helixforge.io.bam import (
            JunctionExtractor,
            aggregate_junctions_multi_sample,
            aggregate_junctions_from_files,
            multi_sample_to_standard_junctions,
        )

        junctions: dict[str, list] = {}

        # From BAM files
        if self.bam_files:
            n_bams = len(self.bam_files)
            if n_bams == 1:
                # Single BAM - use JunctionExtractor directly
                self._log(f"Extracting junctions from BAM: {self.bam_files[0].name}...", indent=1)
                extractor = JunctionExtractor(
                    self.bam_files[0],
                    min_reads=1,  # Get all, filter later
                )
                junctions = extractor.extract_all(min_reads=1)
                self._primary_extractor = extractor
            else:
                # Multiple BAMs - use multi-sample aggregation
                # Show progress for each BAM
                self._log(f"Extracting junctions from {n_bams} BAM files...", indent=1)
                for i, bam_path in enumerate(self.bam_files, 1):
                    self._log(f"Processing BAM {i}/{n_bams}: {bam_path.name}...", indent=2)

                self._log(f"Aggregating junctions (min_tissues={self.config.min_tissues})...", indent=1)
                multi_junctions = aggregate_junctions_multi_sample(
                    bam_paths=self.bam_files,
                    min_reads_per_sample=1,
                    min_samples=self.config.min_tissues,
                )
                junctions = multi_sample_to_standard_junctions(multi_junctions)

                # Keep first extractor for coverage
                self._log("Opening primary BAM for coverage analysis...", indent=1)
                self._primary_extractor = JunctionExtractor(
                    self.bam_files[0],
                    min_reads=1,
                )

        # From junction files
        if self.junction_files:
            n_files = len(self.junction_files)
            if n_files == 1:
                from helixforge.io.bam import load_junctions_auto
                self._log(f"Loading junctions from: {self.junction_files[0].name}...", indent=1)
                file_junctions = load_junctions_auto(
                    self.junction_files[0],
                    min_reads=1,
                )
                for seqid, juncs in file_junctions.items():
                    if seqid not in junctions:
                        junctions[seqid] = []
                    junctions[seqid].extend(juncs)
            else:
                self._log(f"Loading junctions from {n_files} files...", indent=1)
                for i, junc_path in enumerate(self.junction_files, 1):
                    self._log(f"Processing file {i}/{n_files}: {junc_path.name}...", indent=2)

                self._log(f"Aggregating junctions (min_tissues={self.config.min_tissues})...", indent=1)
                multi_junctions = aggregate_junctions_from_files(
                    file_paths=self.junction_files,
                    min_reads_per_sample=1,
                    min_samples=self.config.min_tissues,
                )
                file_junctions = multi_sample_to_standard_junctions(multi_junctions)
                for seqid, juncs in file_junctions.items():
                    if seqid not in junctions:
                        junctions[seqid] = []
                    junctions[seqid].extend(juncs)

        total = sum(len(j) for j in junctions.values())
        self._log(f"Loaded {total:,} junctions from {len(junctions)} scaffold(s)", indent=1)
        logger.info(f"Loaded {total:,} junctions from {len(junctions)} scaffolds")

        return junctions

    def refine_gene(self, gene: "GeneModel") -> RefinedGene:
        """Refine a single gene through all pipeline steps.

        Steps:
        1. Splice site correction
        2. Boundary adjustment
        3. Confidence scoring
        4. Evidence scoring
        5. Flag aggregation

        Args:
            gene: Gene model to refine.

        Returns:
            RefinedGene with all scores and attributes.
        """
        original_gene = deepcopy(gene)

        # Step 1: Splice refinement
        try:
            refined_gene, splice_report = self._splice_refiner.refine_gene(gene)
            splice_corrections = splice_report.n_corrected
        except Exception as e:
            logger.warning(f"Splice refinement failed for {gene.gene_id}: {e}")
            refined_gene = gene
            splice_report = None
            splice_corrections = 0

        # Step 2: Boundary adjustment
        boundary_adjusted = False
        if self._boundary_adjuster and self.config.adjust_boundaries:
            try:
                refined_gene, start_info = self._boundary_adjuster.adjust_start_codon(
                    refined_gene,
                    search_window=self.config.boundary_search_window,
                )
                refined_gene, stop_info = self._boundary_adjuster.adjust_stop_codon(
                    refined_gene,
                    search_window=self.config.boundary_search_window,
                )
                boundary_adjusted = (
                    start_info.get("n_adjusted", 0) > 0
                    or stop_info.get("n_adjusted", 0) > 0
                )
            except Exception as e:
                logger.warning(f"Boundary adjustment failed for {gene.gene_id}: {e}")

        # Step 3: Confidence scoring
        confidence = None
        try:
            confidence = self._confidence_calculator.score_gene(refined_gene)
        except Exception as e:
            logger.warning(f"Confidence scoring failed for {gene.gene_id}: {e}")

        # Step 4: Evidence scoring
        evidence = None
        try:
            evidence = self._evidence_scorer.score_gene(
                refined_gene,
                junctions=self._junctions,
                extractor=self._primary_extractor,
            )
        except Exception as e:
            logger.warning(f"Evidence scoring failed for {gene.gene_id}: {e}")

        # Step 5: Aggregate flags
        flags = self._aggregate_flags(
            confidence, evidence, splice_report, boundary_adjusted
        )

        # Add attributes to gene
        self._add_gene_attributes(refined_gene, confidence, evidence, flags)

        return RefinedGene(
            gene=refined_gene,
            original_gene=original_gene,
            confidence=confidence,
            evidence=evidence,
            splice_report=splice_report,
            splice_corrections=splice_corrections,
            boundary_adjusted=boundary_adjusted,
            flags=flags,
        )

    def refine_genes(
        self,
        genes: Iterator["GeneModel"] | list["GeneModel"],
        n_workers: int = 1,
    ) -> Iterator[RefinedGene]:
        """Refine multiple genes.

        Args:
            genes: Iterator or list of gene models.
            n_workers: Number of parallel workers.

        Yields:
            RefinedGene for each input gene.
        """
        if n_workers == 1:
            for gene in genes:
                yield self.refine_gene(gene)
        else:
            # Parallel execution
            from concurrent.futures import ThreadPoolExecutor, as_completed

            gene_list = list(genes) if not isinstance(genes, list) else genes

            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(self.refine_gene, gene): gene
                    for gene in gene_list
                }

                for future in as_completed(futures):
                    try:
                        yield future.result()
                    except Exception as e:
                        gene = futures[future]
                        logger.error(f"Failed to refine {gene.gene_id}: {e}")

    def refine_region(
        self,
        genes: list["GeneModel"],
        seqid: str,
        start: int,
        end: int,
    ) -> list[RefinedGene]:
        """Refine genes in a specific region.

        Optimized for parallel execution: pre-loads data for region.

        Args:
            genes: List of genes to refine.
            seqid: Scaffold/chromosome name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            List of RefinedGene results.
        """
        results = []
        for gene in genes:
            if gene.seqid == seqid:
                results.append(self.refine_gene(gene))
        return results

    def _aggregate_flags(
        self,
        confidence: "GeneConfidence | None",
        evidence: "EvidenceScore | None",
        splice_report: "GeneSpliceReport | None",
        boundary_adjusted: bool,
    ) -> list[str]:
        """Aggregate flags from all components."""
        flags = []

        # Confidence flags
        if confidence:
            if confidence.overall_score < self.config.confidence_very_low_threshold:
                flags.append("VERY_LOW_CONF")
            elif confidence.overall_score < self.config.confidence_threshold:
                flags.append("LOW_CONF")

            if confidence.worst_exon_score < 0.6:
                flags.append("WEAK_EXON")

            # Add flags from confidence
            flags.extend(confidence.flags)

        # Evidence flags
        if evidence:
            flags.extend(evidence.flags)

        # Splice flags
        if splice_report:
            if splice_report.n_corrected > 0:
                flags.append("SPLICE_CORRECTED")
            if splice_report.n_noncanonical > 0:
                flags.append("NONCANONICAL_SPLICE")
            if splice_report.n_introns > 0 and splice_report.n_supported == 0:
                flags.append("NO_JUNCTION_SUPPORT")

        # Boundary flag
        if boundary_adjusted:
            flags.append("BOUNDARY_ADJUSTED")

        # Deduplicate while preserving order
        seen = set()
        unique_flags = []
        for flag in flags:
            if flag not in seen:
                seen.add(flag)
                unique_flags.append(flag)

        return unique_flags

    def _add_gene_attributes(
        self,
        gene: "GeneModel",
        confidence: "GeneConfidence | None",
        evidence: "EvidenceScore | None",
        flags: list[str],
    ) -> None:
        """Add quality attributes to gene model."""
        if confidence:
            gene.attributes["confidence_score"] = f"{confidence.overall_score:.4f}"

        if evidence:
            # Evidence score = 1 - AED for consistency
            evidence_score = 1.0 - evidence.aed
            gene.attributes["evidence_score"] = f"{evidence_score:.4f}"
            gene.attributes["aed"] = f"{evidence.aed:.4f}"

            if evidence.n_introns > 0:
                gene.attributes["junction_support"] = (
                    f"{evidence.n_introns_supported}/{evidence.n_introns}"
                )
            gene.attributes["mean_coverage"] = f"{evidence.mean_exon_coverage:.1f}"

        if flags:
            gene.attributes["flags"] = ",".join(flags)

    def close(self) -> None:
        """Close any open resources."""
        if self._primary_extractor is not None:
            self._primary_extractor.close()
            self._primary_extractor = None


# =============================================================================
# Report Writer
# =============================================================================


class RefineReportWriter:
    """Write refinement reports."""

    @staticmethod
    def write_report(
        results: list[RefinedGene],
        output_path: Path,
    ) -> None:
        """Write comprehensive refinement report TSV.

        Args:
            results: List of RefinedGene results.
            output_path: Output file path.
        """
        headers = [
            "gene_id",
            "seqid",
            "start",
            "end",
            "strand",
            "n_exons",
            "confidence_score",
            "evidence_score",
            "aed",
            "junction_support",
            "mean_coverage",
            "splice_corrections",
            "boundary_adjusted",
            "flags",
        ]

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
            writer.writeheader()

            for result in results:
                writer.writerow(result.to_report_dict())

        logger.info(f"Wrote refinement report with {len(results)} genes to {output_path}")

    @staticmethod
    def write_summary(results: list[RefinedGene]) -> dict[str, Any]:
        """Generate summary statistics.

        Args:
            results: List of RefinedGene results.

        Returns:
            Dictionary with summary statistics.
        """
        total = len(results)
        if total == 0:
            return {"total_genes": 0}

        # Confidence stats
        conf_scores = [
            r.confidence_score for r in results
            if r.confidence_score is not None
        ]

        # Evidence stats
        ev_scores = [
            r.evidence_score for r in results
            if r.evidence_score is not None
        ]

        # AED stats
        aed_scores = [
            r.aed for r in results
            if r.aed is not None
        ]

        # Splice stats
        splice_corrected = sum(1 for r in results if r.splice_corrections > 0)
        total_corrections = sum(r.splice_corrections for r in results)

        # Junction support
        full_support = sum(
            1 for r in results
            if r.evidence
            and r.evidence.n_introns > 0
            and r.evidence.n_introns_supported == r.evidence.n_introns
        )
        multi_exon = sum(
            1 for r in results
            if r.evidence and r.evidence.n_introns > 0
        )

        return {
            "total_genes": total,
            "mean_confidence": (
                sum(conf_scores) / len(conf_scores)
                if conf_scores else None
            ),
            "mean_evidence_score": (
                sum(ev_scores) / len(ev_scores)
                if ev_scores else None
            ),
            "mean_aed": (
                sum(aed_scores) / len(aed_scores)
                if aed_scores else None
            ),
            "genes_splice_corrected": splice_corrected,
            "total_splice_corrections": total_corrections,
            "multi_exon_genes": multi_exon,
            "full_junction_support": full_support,
            "boundary_adjusted": sum(1 for r in results if r.boundary_adjusted),
        }

    @staticmethod
    def format_summary(summary: dict[str, Any]) -> str:
        """Format summary as text.

        Args:
            summary: Summary dictionary from write_summary().

        Returns:
            Formatted text string.
        """
        lines = [
            "=" * 50,
            "REFINEMENT SUMMARY",
            "=" * 50,
            f"Total genes:             {summary['total_genes']:,}",
        ]

        if summary.get("mean_confidence"):
            lines.append(f"Mean confidence:         {summary['mean_confidence']:.4f}")

        if summary.get("mean_evidence_score"):
            lines.append(f"Mean evidence score:     {summary['mean_evidence_score']:.4f}")

        if summary.get("mean_aed"):
            lines.append(f"Mean AED:                {summary['mean_aed']:.4f}")

        lines.extend([
            f"Splice corrected:        {summary.get('genes_splice_corrected', 0):,} genes",
            f"Total corrections:       {summary.get('total_splice_corrections', 0):,}",
            f"Full junction support:   {summary.get('full_junction_support', 0):,} / "
            f"{summary.get('multi_exon_genes', 0):,} multi-exon genes",
            f"Boundaries adjusted:     {summary.get('boundary_adjusted', 0):,} genes",
            "=" * 50,
        ])

        return "\n".join(lines)

    @staticmethod
    def write_summary_tsv(results: list[RefinedGene], output_path: Path) -> None:
        """Write summary report TSV.

        Args:
            results: List of RefinedGene results.
            output_path: Output file path.
        """
        RefineReportWriter.write_report(results, output_path)

    @staticmethod
    def write_splice_details_tsv(
        results: list[RefinedGene],
        output_path: Path,
    ) -> None:
        """Write detailed splice corrections TSV.

        Args:
            results: List of RefinedGene results.
            output_path: Output file path.
        """
        headers = [
            "gene_id",
            "transcript_id",
            "seqid",
            "intron_idx",
            "site",
            "original_position",
            "corrected_position",
            "shift",
            "reason",
            "confidence",
            "rnaseq_reads",
        ]

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
            writer.writeheader()

            for result in results:
                if result.splice_report is None:
                    continue

                for correction in result.splice_report.corrections:
                    writer.writerow({
                        "gene_id": correction.gene_id,
                        "transcript_id": correction.transcript_id,
                        "seqid": result.gene.seqid,
                        "intron_idx": correction.intron_index,
                        "site": correction.site,
                        "original_position": correction.original_position,
                        "corrected_position": correction.corrected_position,
                        "shift": correction.shift,
                        "reason": correction.reason,
                        "confidence": f"{correction.confidence:.4f}",
                        "rnaseq_reads": correction.rnaseq_reads,
                    })

        logger.info(f"Wrote splice details to {output_path}")

    @staticmethod
    def write_unsupported_introns_bed(
        results: list[RefinedGene],
        output_path: Path,
    ) -> None:
        """Write BED file of unsupported introns.

        Args:
            results: List of RefinedGene results.
            output_path: Output file path.
        """
        with open(output_path, "w") as f:
            for result in results:
                if result.splice_report is None:
                    continue

                # Get intron coordinates from gene's first transcript
                if not result.gene.transcripts:
                    continue

                transcript = result.gene.transcripts[0]
                introns = transcript.introns

                # Write unsupported introns
                for idx in result.splice_report.unsupported_introns:
                    if idx < len(introns):
                        intron_start, intron_end = introns[idx]
                        # BED is 0-based half-open (already in that format)
                        name = f"{result.gene.gene_id}_intron{idx}"
                        strand = result.gene.strand
                        f.write(
                            f"{result.gene.seqid}\t{intron_start}\t{intron_end}\t"
                            f"{name}\t0\t{strand}\n"
                        )

        logger.info(f"Wrote unsupported introns BED to {output_path}")
