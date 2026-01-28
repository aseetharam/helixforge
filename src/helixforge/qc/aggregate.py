"""QC aggregation for combining module results.

This module provides the QCAggregator class that combines results from
confidence, splice, and homology modules into unified GeneQC objects.

Example:
    >>> from helixforge.qc import QCAggregator, QCAggregatorConfig
    >>> config = QCAggregatorConfig()
    >>> aggregator = QCAggregator(config)
    >>> gene_qcs = aggregator.aggregate(
    ...     confidence_results=confidence_scores,
    ...     splice_results=splice_reports,
    ...     homology_results=validation_results,
    ... )
"""

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import attrs

from helixforge.qc.flags import (
    FlagCategory,
    FlagSeverity,
    Flags,
    GeneQC,
    QCFlag,
)


# =============================================================================
# Configuration
# =============================================================================


@attrs.define
class QCAggregatorConfig:
    """Configuration for QC aggregation.

    Attributes:
        confidence_weight: Weight for confidence score in overall QC.
        splice_weight: Weight for splice score in overall QC.
        homology_weight: Weight for homology score in overall QC.
        high_confidence_threshold: Minimum score for high confidence tier.
        medium_confidence_threshold: Minimum score for medium tier.
        low_confidence_threshold: Minimum score for low tier.
        weak_exon_threshold: Threshold for weak exon flag.
        high_entropy_threshold: Threshold for high entropy flag.
        min_splice_support: Minimum RNA-seq reads for splice support.
        min_homology_coverage: Minimum coverage for complete homology.
        min_homology_identity: Minimum identity for homology match.
    """

    # Score weights
    confidence_weight: float = 0.5
    splice_weight: float = 0.25
    homology_weight: float = 0.25

    # Tier thresholds
    high_confidence_threshold: float = 0.85
    medium_confidence_threshold: float = 0.70
    low_confidence_threshold: float = 0.50

    # Confidence flag thresholds
    weak_exon_threshold: float = 0.50
    high_entropy_threshold: float = 1.5
    uncertain_boundary_threshold: float = 0.60

    # Splice thresholds
    min_splice_support: int = 3
    short_intron_threshold: int = 20

    # Homology thresholds
    min_homology_coverage: float = 0.80
    min_homology_identity: float = 0.50

    # Combined AED weights (should sum to 1.0)
    # These weight the contribution of each evidence type to combined_aed
    aed_rnaseq_weight: float = 0.4  # Weight for RNA-seq evidence
    aed_homology_weight: float = 0.4  # Weight for homology evidence
    aed_confidence_weight: float = 0.2  # Weight for Helixer confidence

    def validate(self) -> None:
        """Validate configuration values.

        Raises:
            ValueError: If any configuration value is invalid.
        """
        # Check tier weights sum to 1
        total = self.confidence_weight + self.splice_weight + self.homology_weight
        if abs(total - 1.0) > 0.001:
            raise ValueError(f"Tier weights must sum to 1.0, got {total}")

        # Check AED weights sum to 1
        aed_total = self.aed_rnaseq_weight + self.aed_homology_weight + self.aed_confidence_weight
        if abs(aed_total - 1.0) > 0.001:
            raise ValueError(f"AED weights must sum to 1.0, got {aed_total}")

        # Check thresholds are in order
        if not (0 <= self.low_confidence_threshold <= self.medium_confidence_threshold
                <= self.high_confidence_threshold <= 1):
            raise ValueError("Confidence thresholds must be ordered: low <= medium <= high")


# =============================================================================
# QC Aggregator
# =============================================================================


@attrs.define
class QCAggregator:
    """Aggregates QC results from multiple modules.

    Combines confidence scores, splice refinement reports, and homology
    validation results into unified GeneQC objects with appropriate flags.

    Attributes:
        config: Aggregator configuration.

    Example:
        >>> config = QCAggregatorConfig()
        >>> aggregator = QCAggregator(config)
        >>> gene_qcs = aggregator.aggregate(
        ...     confidence_results=confidence_scores,
        ...     splice_results=splice_reports,
        ... )
    """

    config: QCAggregatorConfig = attrs.Factory(QCAggregatorConfig)

    def aggregate(
        self,
        gene_ids: Iterable[str] | None = None,
        confidence_results: dict[str, Any] | None = None,
        splice_results: dict[str, Any] | None = None,
        homology_results: dict[str, Any] | None = None,
    ) -> dict[str, GeneQC]:
        """Aggregate results from all modules.

        Args:
            gene_ids: Explicit list of gene IDs to process. If None,
                inferred from provided results.
            confidence_results: Dict mapping gene_id to confidence metrics.
                Expected keys: overall_score, mean_prob, min_prob, entropy,
                boundary_sharpness, exon_min, flags (list of flag codes).
            splice_results: Dict mapping gene_id to splice metrics.
                Expected keys: canonical_count, noncanonical_count,
                supported_junctions, unsupported_junctions, corrections,
                flags (list of flag codes).
            homology_results: Dict mapping gene_id to homology validation.
                Expected keys: status, query_coverage, subject_coverage,
                identity, is_chimeric, is_fragmented, te_overlap,
                flags (list of flag codes).

        Returns:
            Dict mapping gene_id to GeneQC objects.
        """
        # Collect all gene IDs
        if gene_ids is None:
            all_ids: set[str] = set()
            if confidence_results:
                all_ids.update(confidence_results.keys())
            if splice_results:
                all_ids.update(splice_results.keys())
            if homology_results:
                all_ids.update(homology_results.keys())
            gene_ids = all_ids

        # Build GeneQC for each gene
        results: dict[str, GeneQC] = {}
        for gene_id in gene_ids:
            gene_qc = GeneQC(gene_id=gene_id)

            # Add confidence data
            if confidence_results and gene_id in confidence_results:
                self._add_confidence(gene_qc, confidence_results[gene_id])

            # Add splice data
            if splice_results and gene_id in splice_results:
                self._add_splice(gene_qc, splice_results[gene_id])

            # Add homology data
            if homology_results and gene_id in homology_results:
                self._add_homology(gene_qc, homology_results[gene_id])

            # Calculate combined AED
            self._calculate_combined_aed(gene_qc)

            # Classify tier
            gene_qc.classify_tier(
                high_threshold=self.config.high_confidence_threshold,
                medium_threshold=self.config.medium_confidence_threshold,
                low_threshold=self.config.low_confidence_threshold,
            )

            results[gene_id] = gene_qc

        return results

    def _calculate_combined_aed(self, gene_qc: GeneQC) -> None:
        """Calculate combined AED incorporating all evidence types.

        The combined AED weights RNA-seq evidence, homology evidence, and
        Helixer confidence to provide a more complete picture than RNA-seq
        alone. This addresses the issue where genes with good homology but
        no RNA-seq expression would incorrectly appear poorly supported.

        Formula:
            combined_evidence = (w_rnaseq * rnaseq_score +
                                 w_homology * homology_score +
                                 w_confidence * confidence_score)
            combined_aed = 1 - combined_evidence

        Args:
            gene_qc: The GeneQC object to update.
        """
        # Extract RNA-seq AED from metrics (stored as 'aed' in refine report)
        rnaseq_aed = None
        if gene_qc.splice_metrics and "aed" in gene_qc.splice_metrics:
            try:
                rnaseq_aed = float(gene_qc.splice_metrics["aed"])
            except (ValueError, TypeError):
                pass
        elif gene_qc.confidence_metrics and "aed" in gene_qc.confidence_metrics:
            try:
                rnaseq_aed = float(gene_qc.confidence_metrics["aed"])
            except (ValueError, TypeError):
                pass

        # Store rnaseq_aed for clarity
        gene_qc.rnaseq_aed = rnaseq_aed

        # Convert RNA-seq AED to evidence score (0-1, higher = better)
        rnaseq_evidence = (1.0 - rnaseq_aed) if rnaseq_aed is not None else 0.0

        # Get homology evidence (0-1, higher = better)
        # homology_score is already query_coverage * identity (normalized)
        homology_evidence = gene_qc.homology_score if gene_qc.homology_score is not None else 0.0

        # Get confidence evidence (0-1, higher = better)
        confidence_evidence = gene_qc.confidence_score if gene_qc.confidence_score is not None else 0.0

        # Calculate weighted combined evidence
        w_rnaseq = self.config.aed_rnaseq_weight
        w_homology = self.config.aed_homology_weight
        w_confidence = self.config.aed_confidence_weight

        # Handle missing evidence types by redistributing weights
        available_weight = 0.0
        combined_evidence = 0.0

        if rnaseq_aed is not None:
            combined_evidence += w_rnaseq * rnaseq_evidence
            available_weight += w_rnaseq
        if gene_qc.homology_score is not None:
            combined_evidence += w_homology * homology_evidence
            available_weight += w_homology
        if gene_qc.confidence_score is not None:
            combined_evidence += w_confidence * confidence_evidence
            available_weight += w_confidence

        # Normalize by available weight to handle missing data
        if available_weight > 0:
            combined_evidence = combined_evidence / available_weight
            gene_qc.combined_aed = 1.0 - combined_evidence
        else:
            gene_qc.combined_aed = None

    def _add_confidence(self, gene_qc: GeneQC, data: dict[str, Any]) -> None:
        """Add confidence module results to GeneQC.

        Args:
            gene_qc: The GeneQC object to update.
            data: Confidence results dict.
        """
        # Store overall score - check multiple possible column names
        if "overall_score" in data:
            gene_qc.confidence_score = data["overall_score"]
        elif "confidence_score" in data:
            gene_qc.confidence_score = data["confidence_score"]
        elif "score" in data:
            gene_qc.confidence_score = data["score"]

        # Store detailed metrics
        gene_qc.confidence_metrics = {
            k: v for k, v in data.items()
            if k not in ("flags", "gene_id")
        }

        # Apply confidence-based flags
        score = gene_qc.confidence_score
        if score is not None:
            if score < self.config.low_confidence_threshold:
                gene_qc.add_flag(Flags.VERY_LOW_CONFIDENCE)
            elif score < self.config.medium_confidence_threshold:
                gene_qc.add_flag(Flags.LOW_CONFIDENCE)

        # Check exon minimum
        if "exon_min" in data:
            if data["exon_min"] < self.config.weak_exon_threshold:
                gene_qc.add_flag(Flags.WEAK_EXON)

        # Check entropy
        if "entropy" in data:
            if data["entropy"] > self.config.high_entropy_threshold:
                gene_qc.add_flag(Flags.HIGH_ENTROPY)

        # Check boundary sharpness
        if "boundary_sharpness" in data:
            if data["boundary_sharpness"] < self.config.uncertain_boundary_threshold:
                gene_qc.add_flag(Flags.UNCERTAIN_BOUNDARY)

        # Add flags from the module itself
        if "flags" in data and data["flags"]:
            flag_data = data["flags"]
            # Handle both string (comma-separated) and list formats
            if isinstance(flag_data, str):
                flag_codes = [f.strip() for f in flag_data.split(",") if f.strip()]
            else:
                flag_codes = flag_data
            for flag_code in flag_codes:
                flag = Flags.get_by_code(flag_code)
                if flag:
                    gene_qc.add_flag(flag)

    def _add_splice(self, gene_qc: GeneQC, data: dict[str, Any]) -> None:
        """Add splice module results to GeneQC.

        Args:
            gene_qc: The GeneQC object to update.
            data: Splice results dict (may come from refine report or splice module).
        """
        # Calculate splice score - handle multiple possible column formats
        total_junctions = data.get("canonical_count", 0) + data.get("noncanonical_count", 0)
        if total_junctions > 0:
            canonical = data.get("canonical_count", 0)
            supported = data.get("supported_junctions", 0)
            # Score based on canonical and supported proportions
            canonical_ratio = canonical / total_junctions
            supported_ratio = supported / total_junctions if total_junctions > 0 else 0
            gene_qc.splice_score = 0.6 * canonical_ratio + 0.4 * supported_ratio
        elif "junction_support" in data:
            # Refine report format - junction_support may be a ratio (0-1) or fraction "5/5"
            junction_support = data.get("junction_support")
            if junction_support is not None:
                # Handle fraction format like "5/5" or "N/A"
                if isinstance(junction_support, str):
                    if "/" in junction_support:
                        try:
                            num, denom = junction_support.split("/")
                            if denom and int(denom) > 0:
                                gene_qc.splice_score = int(num) / int(denom)
                        except (ValueError, ZeroDivisionError):
                            pass  # Leave splice_score as None
                    # "N/A" or other non-parseable strings are ignored
                elif isinstance(junction_support, (int, float)):
                    gene_qc.splice_score = float(junction_support)
        elif "evidence_score" in data:
            # Refine report also has evidence_score which incorporates splice/junction info
            evidence_score = data.get("evidence_score")
            if evidence_score is not None and gene_qc.splice_score is None:
                gene_qc.splice_score = evidence_score

        # Store detailed metrics
        gene_qc.splice_metrics = {
            k: v for k, v in data.items()
            if k not in ("flags", "gene_id")
        }

        # Check for unsupported junctions
        unsupported = data.get("unsupported_junctions", 0)
        if unsupported > 0:
            gene_qc.add_flag(Flags.UNSUPPORTED_JUNCTION)

        # Check for non-canonical splice sites
        noncanonical = data.get("noncanonical_count", 0)
        if noncanonical > 0:
            gene_qc.add_flag(Flags.NON_CANONICAL_SPLICE)

        # Check for corrections made - handle both column names
        corrections = data.get("corrections", data.get("splice_corrections", 0))
        if corrections and corrections > 0:
            gene_qc.add_flag(Flags.SHIFTED_JUNCTION)

        # Check for short introns
        if "min_intron_length" in data:
            if data["min_intron_length"] < self.config.short_intron_threshold:
                gene_qc.add_flag(Flags.SHORT_INTRON)

        # Add flags from the module itself
        if "flags" in data and data["flags"]:
            flag_data = data["flags"]
            # Handle both string (comma-separated) and list formats
            if isinstance(flag_data, str):
                flag_codes = [f.strip() for f in flag_data.split(",") if f.strip()]
            else:
                flag_codes = flag_data
            for flag_code in flag_codes:
                flag = Flags.get_by_code(flag_code)
                if flag:
                    gene_qc.add_flag(flag)

    def _add_homology(self, gene_qc: GeneQC, data: dict[str, Any]) -> None:
        """Add homology module results to GeneQC.

        Args:
            gene_qc: The GeneQC object to update.
            data: Homology results dict.
        """
        # Calculate homology score - ensure numeric values
        try:
            coverage = float(data.get("query_coverage") or 0)
        except (ValueError, TypeError):
            coverage = 0
        try:
            identity = float(data.get("identity") or 0)
        except (ValueError, TypeError):
            identity = 0

        # Get status (case-insensitive)
        status = data.get("status", "")
        status_lower = status.lower() if isinstance(status, str) else ""

        if coverage > 0 and identity > 0:
            # Identity is typically 0-100, coverage is 0-1
            # Normalize identity to 0-1 if needed
            if identity > 1:
                identity = identity / 100.0
            gene_qc.homology_score = coverage * identity
        elif status_lower == "complete":
            gene_qc.homology_score = 1.0
        elif status_lower == "partial":
            gene_qc.homology_score = 0.5
        elif status_lower == "no_hit":
            gene_qc.homology_score = 0.0

        # Store detailed metrics
        gene_qc.homology_metrics = {
            k: v for k, v in data.items()
            if k not in ("flags", "gene_id")
        }

        # Check homology status (case-insensitive)
        status = data.get("status", "")
        if isinstance(status, str):
            status_lower = status.lower()
            if status_lower == "no_hit":
                gene_qc.add_flag(Flags.NO_HOMOLOGY)
            elif status_lower == "partial":
                gene_qc.add_flag(Flags.PARTIAL_HOMOLOGY)

        # Check for coverage below threshold
        if coverage < self.config.min_homology_coverage and coverage > 0:
            gene_qc.add_flag(Flags.PARTIAL_HOMOLOGY)

        # Check for chimeric - handle string "yes"/"no" or boolean
        is_chimeric = data.get("is_chimeric", False)
        if is_chimeric is True or (isinstance(is_chimeric, str) and is_chimeric.lower() in ("yes", "true", "1")):
            gene_qc.add_flag(Flags.CHIMERIC)

        # Check for fragmented - handle string "yes"/"no" or boolean
        is_fragmented = data.get("is_fragmented", False)
        if is_fragmented is True or (isinstance(is_fragmented, str) and is_fragmented.lower() in ("yes", "true", "1")):
            gene_qc.add_flag(Flags.FRAGMENTED)

        # Check for TE overlap
        te_overlap = data.get("te_overlap", 0) or 0
        try:
            te_overlap = float(te_overlap)
        except (ValueError, TypeError):
            te_overlap = 0
        if te_overlap > 0:
            gene_qc.add_flag(Flags.TE_OVERLAP)

        # Add flags from the module itself
        if "flags" in data and data["flags"]:
            flag_data = data["flags"]
            # Handle both string (comma-separated) and list formats
            if isinstance(flag_data, str):
                flag_codes = [f.strip() for f in flag_data.split(",") if f.strip()]
            else:
                flag_codes = flag_data
            for flag_code in flag_codes:
                flag = Flags.get_by_code(flag_code)
                if flag:
                    gene_qc.add_flag(flag)

    def aggregate_from_files(
        self,
        refine_tsv: Path | str | None = None,
        confidence_tsv: Path | str | None = None,
        splice_tsv: Path | str | None = None,
        homology_tsv: Path | str | None = None,
    ) -> dict[str, GeneQC]:
        """Aggregate results from TSV files.

        Args:
            refine_tsv: Path to refine report TSV (contains confidence, splice, evidence).
            confidence_tsv: Path to confidence scores TSV.
            splice_tsv: Path to splice report TSV (deprecated).
            homology_tsv: Path to homology validation TSV.

        Returns:
            Dict mapping gene_id to GeneQC objects.

        Note:
            If refine_tsv is provided, it takes precedence over confidence_tsv
            and splice_tsv for those metrics, as refine output contains both.
        """
        confidence_results = None
        splice_results = None
        homology_results = None

        # Load refine TSV - contains confidence, splice, and evidence data
        if refine_tsv:
            refine_results = self._load_tsv(Path(refine_tsv))
            # Use refine results for both confidence and splice
            confidence_results = refine_results
            splice_results = refine_results

        # Confidence TSV overrides if specified separately
        if confidence_tsv:
            confidence_results = self._load_tsv(Path(confidence_tsv))

        # Splice TSV overrides if specified separately
        if splice_tsv:
            splice_results = self._load_tsv(Path(splice_tsv))

        if homology_tsv:
            homology_results = self._load_tsv(Path(homology_tsv))

        return self.aggregate(
            confidence_results=confidence_results,
            splice_results=splice_results,
            homology_results=homology_results,
        )

    def _load_tsv(self, path: Path) -> dict[str, dict[str, Any]]:
        """Load TSV file into dict keyed by gene_id.

        Args:
            path: Path to TSV file.

        Returns:
            Dict mapping gene_id to row data.
        """
        import csv

        results: dict[str, dict[str, Any]] = {}
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gene_id = row.get("gene_id", row.get("Gene", row.get("gene")))
                if not gene_id:
                    continue

                # Convert numeric values
                parsed: dict[str, Any] = {}
                for key, value in row.items():
                    if value == "":
                        parsed[key] = None
                    else:
                        try:
                            if "." in value:
                                parsed[key] = float(value)
                            else:
                                parsed[key] = int(value)
                        except (ValueError, TypeError):
                            parsed[key] = value

                results[gene_id] = parsed

        return results


# =============================================================================
# Summary Functions
# =============================================================================


def summarize_qc_results(gene_qcs: dict[str, GeneQC]) -> dict[str, Any]:
    """Generate summary statistics for QC results.

    Args:
        gene_qcs: Dict mapping gene_id to GeneQC objects.

    Returns:
        Dict with summary statistics.
    """
    total = len(gene_qcs)
    if total == 0:
        return {
            "total_genes": 0,
            "tier_counts": {},
            "tier_percentages": {},
            "flag_counts": {},
            "severity_counts": {},
            "category_counts": {},
            "avg_confidence": None,
            "avg_splice": None,
            "avg_homology": None,
        }

    # Count by tier
    tier_counts: dict[str, int] = {}
    for qc in gene_qcs.values():
        tier_counts[qc.tier] = tier_counts.get(qc.tier, 0) + 1

    tier_percentages = {tier: count / total * 100 for tier, count in tier_counts.items()}

    # Count flags
    flag_counts: dict[str, int] = {}
    severity_counts: dict[str, int] = {s.value: 0 for s in FlagSeverity}
    category_counts: dict[str, int] = {c.value: 0 for c in FlagCategory}

    for qc in gene_qcs.values():
        for flag in qc.flags:
            flag_counts[flag.code] = flag_counts.get(flag.code, 0) + 1
            severity_counts[flag.severity.value] += 1
            category_counts[flag.category.value] += 1

    # Average scores
    conf_scores = [qc.confidence_score for qc in gene_qcs.values() if qc.confidence_score is not None]
    splice_scores = [qc.splice_score for qc in gene_qcs.values() if qc.splice_score is not None]
    homology_scores = [qc.homology_score for qc in gene_qcs.values() if qc.homology_score is not None]

    return {
        "total_genes": total,
        "tier_counts": tier_counts,
        "tier_percentages": tier_percentages,
        "flag_counts": flag_counts,
        "severity_counts": severity_counts,
        "category_counts": category_counts,
        "avg_confidence": sum(conf_scores) / len(conf_scores) if conf_scores else None,
        "avg_splice": sum(splice_scores) / len(splice_scores) if splice_scores else None,
        "avg_homology": sum(homology_scores) / len(homology_scores) if homology_scores else None,
    }


def get_genes_by_tier(
    gene_qcs: dict[str, GeneQC],
    tier: str,
) -> list[GeneQC]:
    """Get all genes in a specific tier.

    Args:
        gene_qcs: Dict mapping gene_id to GeneQC objects.
        tier: Tier to filter by (high, medium, low, reject).

    Returns:
        List of GeneQC objects in the tier.
    """
    return [qc for qc in gene_qcs.values() if qc.tier == tier]


def get_genes_with_flag(
    gene_qcs: dict[str, GeneQC],
    flag: QCFlag,
) -> list[GeneQC]:
    """Get all genes with a specific flag.

    Args:
        gene_qcs: Dict mapping gene_id to GeneQC objects.
        flag: The flag to filter by.

    Returns:
        List of GeneQC objects with the flag.
    """
    return [qc for qc in gene_qcs.values() if qc.has_flag(flag)]


def export_qc_tsv(
    gene_qcs: dict[str, GeneQC],
    output_path: Path | str,
    include_metrics: bool = True,
) -> None:
    """Export QC results to TSV file.

    Args:
        gene_qcs: Dict mapping gene_id to GeneQC objects.
        output_path: Path to output TSV file.
        include_metrics: Whether to include detailed metrics columns.
    """
    import csv

    # Base columns
    columns = [
        "gene_id",
        "tier",
        "flag_count",
        "flag_codes",
        "max_severity",
        "confidence_score",
        "splice_score",
        "homology_score",
        "rnaseq_aed",
        "combined_aed",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()

        def format_score(score: Any) -> str:
            """Format score as float string, handling non-numeric gracefully."""
            if score is None:
                return ""
            if isinstance(score, (int, float)):
                return f"{score:.4f}"
            # Try to convert string to float
            try:
                return f"{float(score):.4f}"
            except (ValueError, TypeError):
                return ""

        for gene_id, qc in sorted(gene_qcs.items()):
            row = {
                "gene_id": gene_id,
                "tier": qc.tier,
                "flag_count": qc.flag_count,
                "flag_codes": ",".join(qc.flag_codes),
                "max_severity": qc.max_severity.value if qc.max_severity else "",
                "confidence_score": format_score(qc.confidence_score),
                "splice_score": format_score(qc.splice_score),
                "homology_score": format_score(qc.homology_score),
                "rnaseq_aed": format_score(qc.rnaseq_aed),
                "combined_aed": format_score(qc.combined_aed),
            }
            writer.writerow(row)
