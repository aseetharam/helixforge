"""Filtering gene models based on QC criteria.

This module provides tools for filtering gene predictions based on
QC flags, confidence scores, tiers, and other criteria.

Example:
    >>> from helixforge.qc import GeneFilter, FilterCriteria
    >>> from helixforge.qc.flags import Flags, FlagSeverity
    >>>
    >>> # Create filter with custom criteria
    >>> criteria = FilterCriteria(
    ...     min_confidence=0.7,
    ...     max_severity=FlagSeverity.WARNING,
    ...     exclude_flags=[Flags.INTERNAL_STOP],
    ... )
    >>> gene_filter = GeneFilter(criteria)
    >>> result = gene_filter.apply(gene_qcs)
    >>>
    >>> # Or use preset profiles
    >>> result = GeneFilter.high_confidence().apply(gene_qcs)
    >>> result = GeneFilter.publication_ready().apply(gene_qcs)
"""

from collections.abc import Callable, Iterable
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
# Filter Criteria
# =============================================================================


@attrs.define
class FilterCriteria:
    """Criteria for filtering gene models.

    Attributes:
        name: Name of the filter profile.
        description: Description of what this filter selects.
        min_confidence: Minimum confidence score (0-1).
        max_severity: Maximum allowed flag severity.
        allowed_tiers: List of allowed tiers (None = all).
        exclude_flags: Flags that disqualify a gene.
        require_flags: Flags that must be present.
        exclude_categories: Flag categories that disqualify a gene.
        max_flag_count: Maximum number of flags allowed.
        custom_filter: Custom filter function.
    """

    name: str = "custom"
    description: str = ""

    # Score thresholds
    min_confidence: float | None = None
    min_splice_score: float | None = None
    min_homology_score: float | None = None

    # Severity threshold
    max_severity: FlagSeverity | None = None

    # Tier filter
    allowed_tiers: list[str] | None = None

    # Flag-based filters
    exclude_flags: list[QCFlag] = attrs.Factory(list)
    require_flags: list[QCFlag] = attrs.Factory(list)
    exclude_categories: list[FlagCategory] = attrs.Factory(list)

    # Count limits
    max_flag_count: int | None = None

    # Custom filter
    custom_filter: Callable[[GeneQC], bool] | None = None


# =============================================================================
# Filter Result
# =============================================================================


@attrs.define
class FilterResult:
    """Result of a filtering operation.

    Attributes:
        passed: GeneQC objects that passed the filter.
        failed: GeneQC objects that failed the filter.
        criteria: The criteria used for filtering.
        statistics: Statistics about the filtering.
    """

    passed: list[GeneQC]
    failed: list[GeneQC]
    criteria: FilterCriteria
    statistics: dict[str, Any] = attrs.Factory(dict)

    @property
    def total_count(self) -> int:
        """Total number of genes processed."""
        return len(self.passed) + len(self.failed)

    @property
    def pass_count(self) -> int:
        """Number of genes that passed."""
        return len(self.passed)

    @property
    def fail_count(self) -> int:
        """Number of genes that failed."""
        return len(self.failed)

    @property
    def pass_rate(self) -> float:
        """Proportion of genes that passed."""
        if self.total_count == 0:
            return 0.0
        return self.pass_count / self.total_count

    def passed_ids(self) -> list[str]:
        """Get list of gene IDs that passed."""
        return [qc.gene_id for qc in self.passed]

    def failed_ids(self) -> list[str]:
        """Get list of gene IDs that failed."""
        return [qc.gene_id for qc in self.failed]


# =============================================================================
# Gene Filter
# =============================================================================


@attrs.define
class GeneFilter:
    """Filter genes based on QC criteria.

    Provides methods to filter GeneQC objects based on various criteria
    including confidence scores, flag severity, tiers, and custom functions.

    Example:
        >>> criteria = FilterCriteria(min_confidence=0.7)
        >>> gene_filter = GeneFilter(criteria)
        >>> result = gene_filter.apply(gene_qcs)
        >>> print(f"Passed: {result.pass_count}/{result.total_count}")
    """

    criteria: FilterCriteria

    def apply(
        self,
        gene_qcs: dict[str, GeneQC] | Iterable[GeneQC],
    ) -> FilterResult:
        """Apply filter criteria to genes.

        Args:
            gene_qcs: Dict or iterable of GeneQC objects.

        Returns:
            FilterResult with passed/failed genes and statistics.
        """
        # Convert dict to list if needed
        if isinstance(gene_qcs, dict):
            genes = list(gene_qcs.values())
        else:
            genes = list(gene_qcs)

        passed: list[GeneQC] = []
        failed: list[GeneQC] = []
        fail_reasons: dict[str, int] = {}

        for qc in genes:
            is_pass, reason = self._check_gene(qc)
            if is_pass:
                passed.append(qc)
            else:
                failed.append(qc)
                if reason:
                    fail_reasons[reason] = fail_reasons.get(reason, 0) + 1

        statistics = {
            "total": len(genes),
            "passed": len(passed),
            "failed": len(failed),
            "pass_rate": len(passed) / len(genes) if genes else 0,
            "fail_reasons": fail_reasons,
        }

        return FilterResult(
            passed=passed,
            failed=failed,
            criteria=self.criteria,
            statistics=statistics,
        )

    def _check_gene(self, qc: GeneQC) -> tuple[bool, str | None]:
        """Check if a gene passes all criteria.

        Args:
            qc: The GeneQC to check.

        Returns:
            Tuple of (passed, failure_reason).
        """
        criteria = self.criteria

        # Check confidence score
        if criteria.min_confidence is not None:
            if qc.confidence_score is None or qc.confidence_score < criteria.min_confidence:
                return False, "low_confidence"

        # Check splice score
        if criteria.min_splice_score is not None:
            if qc.splice_score is None or qc.splice_score < criteria.min_splice_score:
                return False, "low_splice_score"

        # Check homology score
        if criteria.min_homology_score is not None:
            if qc.homology_score is None or qc.homology_score < criteria.min_homology_score:
                return False, "low_homology_score"

        # Check max severity
        if criteria.max_severity is not None:
            if qc.max_severity is not None and qc.max_severity > criteria.max_severity:
                return False, "high_severity"

        # Check allowed tiers
        if criteria.allowed_tiers is not None:
            if qc.tier not in criteria.allowed_tiers:
                return False, "excluded_tier"

        # Check excluded flags
        if criteria.exclude_flags:
            if qc.has_any_flag(criteria.exclude_flags):
                return False, "excluded_flag"

        # Check required flags
        if criteria.require_flags:
            if not qc.has_all_flags(criteria.require_flags):
                return False, "missing_required_flag"

        # Check excluded categories
        if criteria.exclude_categories:
            for cat in criteria.exclude_categories:
                if qc.flags_by_category(cat):
                    return False, f"excluded_category_{cat.value}"

        # Check flag count
        if criteria.max_flag_count is not None:
            if qc.flag_count > criteria.max_flag_count:
                return False, "too_many_flags"

        # Check custom filter
        if criteria.custom_filter is not None:
            if not criteria.custom_filter(qc):
                return False, "custom_filter"

        return True, None

    # -------------------------------------------------------------------------
    # Preset Filter Profiles
    # -------------------------------------------------------------------------

    @classmethod
    def high_confidence(cls) -> "GeneFilter":
        """Create filter for high-confidence genes only.

        Selects genes with:
        - High tier classification
        - No warning or higher severity flags
        - Confidence >= 0.85
        """
        return cls(FilterCriteria(
            name="high_confidence",
            description="High-confidence genes with no warnings",
            min_confidence=0.85,
            allowed_tiers=["high"],
            max_severity=FlagSeverity.INFO,
        ))

    @classmethod
    def publication_ready(cls) -> "GeneFilter":
        """Create filter for publication-ready genes.

        Selects genes with:
        - High or medium tier
        - No error or critical flags
        - No internal stops or frameshifts
        - Confidence >= 0.70
        """
        return cls(FilterCriteria(
            name="publication_ready",
            description="Publication-ready genes without major issues",
            min_confidence=0.70,
            allowed_tiers=["high", "medium"],
            max_severity=FlagSeverity.WARNING,
            exclude_flags=[
                Flags.INTERNAL_STOP,
                Flags.FRAMESHIFT,
                Flags.CHIMERIC,
            ],
        ))

    @classmethod
    def has_homology(cls) -> "GeneFilter":
        """Create filter for genes with homology support.

        Selects genes with:
        - Homology score > 0
        - No NO_HOMOLOGY flag
        """
        return cls(FilterCriteria(
            name="has_homology",
            description="Genes with protein homology support",
            min_homology_score=0.01,
            exclude_flags=[Flags.NO_HOMOLOGY],
        ))

    @classmethod
    def has_rnaseq_support(cls) -> "GeneFilter":
        """Create filter for genes with RNA-seq support.

        Selects genes without unsupported junctions flag.
        """
        return cls(FilterCriteria(
            name="has_rnaseq_support",
            description="Genes with RNA-seq splice junction support",
            exclude_flags=[Flags.UNSUPPORTED_JUNCTION, Flags.NO_EVIDENCE],
        ))

    @classmethod
    def no_structural_issues(cls) -> "GeneFilter":
        """Create filter for genes without structural issues.

        Selects genes without:
        - Internal stop codons
        - Frameshifts
        - Missing start/stop codons
        """
        return cls(FilterCriteria(
            name="no_structural_issues",
            description="Genes without structural problems",
            exclude_flags=[
                Flags.INTERNAL_STOP,
                Flags.FRAMESHIFT,
                Flags.NO_START_CODON,
                Flags.NO_STOP_CODON,
            ],
        ))

    @classmethod
    def tier_filter(cls, tiers: list[str]) -> "GeneFilter":
        """Create filter for specific tiers.

        Args:
            tiers: List of tiers to include.
        """
        return cls(FilterCriteria(
            name=f"tiers_{'+'.join(tiers)}",
            description=f"Genes in tiers: {', '.join(tiers)}",
            allowed_tiers=tiers,
        ))

    @classmethod
    def custom(
        cls,
        name: str,
        filter_func: Callable[[GeneQC], bool],
        description: str = "",
    ) -> "GeneFilter":
        """Create filter with custom function.

        Args:
            name: Name of the filter.
            filter_func: Function that returns True if gene passes.
            description: Description of the filter.
        """
        return cls(FilterCriteria(
            name=name,
            description=description,
            custom_filter=filter_func,
        ))


# =============================================================================
# Tiered Output Functions
# =============================================================================


def split_by_tier(
    gene_qcs: dict[str, GeneQC] | Iterable[GeneQC],
) -> dict[str, list[GeneQC]]:
    """Split genes into tier groups.

    Args:
        gene_qcs: Dict or iterable of GeneQC objects.

    Returns:
        Dict mapping tier name to list of GeneQC objects.
    """
    if isinstance(gene_qcs, dict):
        genes = list(gene_qcs.values())
    else:
        genes = list(gene_qcs)

    result: dict[str, list[GeneQC]] = {}
    for qc in genes:
        if qc.tier not in result:
            result[qc.tier] = []
        result[qc.tier].append(qc)

    return result


def write_tiered_gene_lists(
    gene_qcs: dict[str, GeneQC],
    output_dir: Path | str,
    prefix: str = "genes",
) -> dict[str, Path]:
    """Write gene ID lists for each tier.

    Args:
        gene_qcs: Dict mapping gene_id to GeneQC objects.
        output_dir: Directory to write files.
        prefix: Prefix for output file names.

    Returns:
        Dict mapping tier to output file path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tiered = split_by_tier(gene_qcs)
    output_files: dict[str, Path] = {}

    for tier, genes in tiered.items():
        output_path = output_dir / f"{prefix}_{tier}.txt"
        with open(output_path, "w") as f:
            for qc in sorted(genes, key=lambda x: x.gene_id):
                f.write(f"{qc.gene_id}\n")
        output_files[tier] = output_path

    return output_files


def write_filtered_gff(
    gene_qcs: dict[str, GeneQC],
    input_gff: Path | str,
    output_gff: Path | str,
    filter_criteria: FilterCriteria | None = None,
    tiers: list[str] | None = None,
) -> int:
    """Write filtered GFF file containing only passing genes.

    Args:
        gene_qcs: Dict mapping gene_id to GeneQC objects.
        input_gff: Path to input GFF3 file.
        output_gff: Path to output GFF3 file.
        filter_criteria: Filter criteria to apply.
        tiers: If specified, only include genes from these tiers.

    Returns:
        Number of genes written.
    """
    # Determine which genes to include
    if filter_criteria:
        gene_filter = GeneFilter(filter_criteria)
        result = gene_filter.apply(gene_qcs)
        include_ids = set(result.passed_ids())
    elif tiers:
        include_ids = {
            qc.gene_id for qc in gene_qcs.values()
            if qc.tier in tiers
        }
    else:
        include_ids = set(gene_qcs.keys())

    # Track current gene and whether to include it
    current_gene: str | None = None
    include_current = False
    genes_written = 0

    with open(input_gff) as fin, open(output_gff, "w") as fout:
        for line in fin:
            # Pass through comments and directives
            if line.startswith("#"):
                fout.write(line)
                continue

            if not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            attributes = parts[8]

            # Parse gene ID from attributes
            gene_id = None
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    gene_id = attr[3:]
                elif attr.startswith("Parent="):
                    # Get the gene ID from the parent
                    parent = attr[7:]
                    if parent.startswith("gene:"):
                        gene_id = parent
                    elif "." in parent:
                        gene_id = parent.split(".")[0]
                    else:
                        gene_id = parent

            if feature_type == "gene":
                if gene_id and gene_id.startswith("ID="):
                    gene_id = gene_id[3:]
                current_gene = gene_id
                include_current = current_gene in include_ids
                if include_current:
                    fout.write(line)
                    genes_written += 1
            else:
                if include_current:
                    fout.write(line)

    return genes_written


# =============================================================================
# Filter Statistics Functions
# =============================================================================


def calculate_filter_statistics(
    original: Iterable[GeneQC],
    filtered: Iterable[GeneQC],
) -> dict[str, Any]:
    """Calculate statistics about filtering.

    Args:
        original: Original gene list.
        filtered: Filtered gene list.

    Returns:
        Dictionary of statistics.
    """
    orig_list = list(original)
    filt_list = list(filtered)

    return {
        "original_count": len(orig_list),
        "filtered_count": len(filt_list),
        "removed_count": len(orig_list) - len(filt_list),
        "retention_rate": len(filt_list) / len(orig_list) if orig_list else 0,
    }


def summarize_filter_results(result: FilterResult) -> str:
    """Generate summary text for filter results.

    Args:
        result: The FilterResult to summarize.

    Returns:
        Multi-line summary string.
    """
    lines = [
        f"Filter: {result.criteria.name}",
        f"Description: {result.criteria.description}",
        "",
        f"Total genes: {result.total_count}",
        f"Passed: {result.pass_count} ({result.pass_rate:.1%})",
        f"Failed: {result.fail_count}",
    ]

    if result.statistics.get("fail_reasons"):
        lines.append("")
        lines.append("Failure reasons:")
        for reason, count in sorted(
            result.statistics["fail_reasons"].items(),
            key=lambda x: -x[1],
        ):
            lines.append(f"  {reason}: {count}")

    return "\n".join(lines)
