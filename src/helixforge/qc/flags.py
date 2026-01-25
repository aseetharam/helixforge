"""QC flags for gene model issues.

This module defines quality control flags that indicate potential
issues with gene predictions. Flags are organized by severity and
category.

Flag Categories:
    - CONFIDENCE: Prediction confidence issues
    - SPLICE: Splice site issues
    - HOMOLOGY: Homology-based issues
    - STRUCTURE: Gene structure issues
    - ANNOTATION: Annotation issues

Example:
    >>> from helixforge.qc.flags import Flags, GeneQC, FlagSeverity
    >>> gene_qc = GeneQC(gene_id="gene001")
    >>> gene_qc.add_flag(Flags.LOW_CONFIDENCE)
    >>> gene_qc.add_flag(Flags.WEAK_EXON)
    >>> print(gene_qc.max_severity)
    FlagSeverity.WARNING
    >>> print(gene_qc.tier)
    medium
"""

from collections.abc import Iterable
from enum import Enum
from typing import Any

import attrs


# =============================================================================
# Enums
# =============================================================================


class FlagSeverity(Enum):
    """Severity levels for QC flags.

    Ordered from least to most severe:
    - INFO: Informational, not necessarily an issue
    - WARNING: Potential issue, may need review
    - ERROR: Definite issue, likely needs correction
    - CRITICAL: Severe issue, gene may be invalid
    """

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"

    def __lt__(self, other: "FlagSeverity") -> bool:
        order = [FlagSeverity.INFO, FlagSeverity.WARNING, FlagSeverity.ERROR, FlagSeverity.CRITICAL]
        return order.index(self) < order.index(other)

    def __le__(self, other: "FlagSeverity") -> bool:
        return self == other or self < other

    def __gt__(self, other: "FlagSeverity") -> bool:
        return not self <= other

    def __ge__(self, other: "FlagSeverity") -> bool:
        return not self < other


class FlagCategory(Enum):
    """Categories of QC flags."""

    CONFIDENCE = "confidence"  # Prediction confidence issues
    SPLICE = "splice"  # Splice site issues
    HOMOLOGY = "homology"  # Homology-based issues
    STRUCTURE = "structure"  # Gene structure issues
    ANNOTATION = "annotation"  # Annotation issues


# =============================================================================
# QCFlag Data Class
# =============================================================================


@attrs.define(frozen=True)
class QCFlag:
    """A quality control flag for a gene.

    Frozen (immutable) attrs class that is hashable for use in sets.

    Attributes:
        code: Short flag code (e.g., "LOW_CONF").
        name: Human-readable flag name.
        description: Detailed description of the issue.
        category: Flag category.
        severity: Severity level.
    """

    code: str
    name: str
    description: str
    category: FlagCategory
    severity: FlagSeverity

    def __hash__(self) -> int:
        return hash(self.code)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, QCFlag):
            return False
        return self.code == other.code


# =============================================================================
# Flags Registry
# =============================================================================


class Flags:
    """Registry of all QC flags.

    All flag definitions are class attributes. Access flags directly:
        >>> Flags.LOW_CONFIDENCE
        >>> Flags.WEAK_EXON
        >>> Flags.get_all()
    """

    # -------------------------------------------------------------------------
    # Confidence flags
    # -------------------------------------------------------------------------

    LOW_CONFIDENCE = QCFlag(
        code="LOW_CONF",
        name="Low Confidence",
        description="Overall confidence score below medium threshold (<0.70)",
        category=FlagCategory.CONFIDENCE,
        severity=FlagSeverity.WARNING,
    )

    VERY_LOW_CONFIDENCE = QCFlag(
        code="VERY_LOW_CONF",
        name="Very Low Confidence",
        description="Overall confidence score very low (<0.50)",
        category=FlagCategory.CONFIDENCE,
        severity=FlagSeverity.ERROR,
    )

    WEAK_EXON = QCFlag(
        code="WEAK_EXON",
        name="Weak Exon",
        description="One or more exons have very low prediction probability",
        category=FlagCategory.CONFIDENCE,
        severity=FlagSeverity.WARNING,
    )

    HIGH_ENTROPY = QCFlag(
        code="HIGH_ENTROPY",
        name="High Entropy",
        description="High prediction uncertainty (entropy) across gene region",
        category=FlagCategory.CONFIDENCE,
        severity=FlagSeverity.WARNING,
    )

    UNCERTAIN_BOUNDARY = QCFlag(
        code="UNCERTAIN_BOUND",
        name="Uncertain Boundary",
        description="Low boundary sharpness at exon/intron transitions",
        category=FlagCategory.CONFIDENCE,
        severity=FlagSeverity.INFO,
    )

    # -------------------------------------------------------------------------
    # Splice flags
    # -------------------------------------------------------------------------

    NON_CANONICAL_SPLICE = QCFlag(
        code="NON_CANON_SPLICE",
        name="Non-canonical Splice Site",
        description="Splice site is not GT-AG, GC-AG, or AT-AC",
        category=FlagCategory.SPLICE,
        severity=FlagSeverity.WARNING,
    )

    UNSUPPORTED_JUNCTION = QCFlag(
        code="UNSUP_JUNCTION",
        name="Unsupported Junction",
        description="Splice junction has no RNA-seq support",
        category=FlagCategory.SPLICE,
        severity=FlagSeverity.WARNING,
    )

    SHIFTED_JUNCTION = QCFlag(
        code="SHIFTED_JUNCTION",
        name="Shifted Junction",
        description="Splice junction was corrected based on RNA-seq evidence",
        category=FlagCategory.SPLICE,
        severity=FlagSeverity.INFO,
    )

    LOW_PWM_SCORE = QCFlag(
        code="LOW_PWM",
        name="Low PWM Score",
        description="Splice site has low position weight matrix score",
        category=FlagCategory.SPLICE,
        severity=FlagSeverity.INFO,
    )

    SHORT_INTRON = QCFlag(
        code="SHORT_INTRON",
        name="Short Intron",
        description="Intron is shorter than typical minimum (< 20bp)",
        category=FlagCategory.SPLICE,
        severity=FlagSeverity.WARNING,
    )

    # -------------------------------------------------------------------------
    # Homology flags
    # -------------------------------------------------------------------------

    NO_HOMOLOGY = QCFlag(
        code="NO_HOMOLOGY",
        name="No Homology",
        description="No homology hits found in protein database",
        category=FlagCategory.HOMOLOGY,
        severity=FlagSeverity.INFO,
    )

    PARTIAL_HOMOLOGY = QCFlag(
        code="PARTIAL_HOMOLOGY",
        name="Partial Homology",
        description="Protein has partial coverage of homology hit",
        category=FlagCategory.HOMOLOGY,
        severity=FlagSeverity.WARNING,
    )

    CHIMERIC = QCFlag(
        code="CHIMERIC",
        name="Chimeric Gene",
        description="Gene appears to be a fusion of multiple genes",
        category=FlagCategory.HOMOLOGY,
        severity=FlagSeverity.ERROR,
    )

    FRAGMENTED = QCFlag(
        code="FRAGMENTED",
        name="Fragmented Gene",
        description="Gene appears to be a fragment of a larger gene",
        category=FlagCategory.HOMOLOGY,
        severity=FlagSeverity.WARNING,
    )

    TE_OVERLAP = QCFlag(
        code="TE_OVERLAP",
        name="TE Overlap",
        description="Gene overlaps transposable element annotation",
        category=FlagCategory.HOMOLOGY,
        severity=FlagSeverity.WARNING,
    )

    # -------------------------------------------------------------------------
    # Structure flags
    # -------------------------------------------------------------------------

    NO_START_CODON = QCFlag(
        code="NO_START",
        name="No Start Codon",
        description="CDS does not begin with a start codon (ATG)",
        category=FlagCategory.STRUCTURE,
        severity=FlagSeverity.ERROR,
    )

    NO_STOP_CODON = QCFlag(
        code="NO_STOP",
        name="No Stop Codon",
        description="CDS does not end with a stop codon",
        category=FlagCategory.STRUCTURE,
        severity=FlagSeverity.ERROR,
    )

    INTERNAL_STOP = QCFlag(
        code="INTERNAL_STOP",
        name="Internal Stop Codon",
        description="CDS contains an internal stop codon",
        category=FlagCategory.STRUCTURE,
        severity=FlagSeverity.CRITICAL,
    )

    FRAMESHIFT = QCFlag(
        code="FRAMESHIFT",
        name="Frameshift",
        description="CDS has a frameshift (length not divisible by 3)",
        category=FlagCategory.STRUCTURE,
        severity=FlagSeverity.ERROR,
    )

    SHORT_CDS = QCFlag(
        code="SHORT_CDS",
        name="Short CDS",
        description="CDS is shorter than typical minimum length",
        category=FlagCategory.STRUCTURE,
        severity=FlagSeverity.WARNING,
    )

    SHORT_EXON = QCFlag(
        code="SHORT_EXON",
        name="Short Exon",
        description="Gene contains a very short exon (< 10bp)",
        category=FlagCategory.STRUCTURE,
        severity=FlagSeverity.INFO,
    )

    LONG_INTRON = QCFlag(
        code="LONG_INTRON",
        name="Long Intron",
        description="Gene contains an unusually long intron",
        category=FlagCategory.STRUCTURE,
        severity=FlagSeverity.INFO,
    )

    # -------------------------------------------------------------------------
    # Annotation flags
    # -------------------------------------------------------------------------

    OVERLAPPING_GENE = QCFlag(
        code="OVERLAP_GENE",
        name="Overlapping Gene",
        description="Gene overlaps another gene on the same strand",
        category=FlagCategory.ANNOTATION,
        severity=FlagSeverity.INFO,
    )

    NO_EVIDENCE = QCFlag(
        code="NO_EVIDENCE",
        name="No Evidence",
        description="No RNA-seq evidence supports this gene",
        category=FlagCategory.ANNOTATION,
        severity=FlagSeverity.WARNING,
    )

    LOW_COVERAGE = QCFlag(
        code="LOW_COVERAGE",
        name="Low Coverage",
        description="Low RNA-seq coverage across gene",
        category=FlagCategory.ANNOTATION,
        severity=FlagSeverity.INFO,
    )

    @classmethod
    def get_all(cls) -> list[QCFlag]:
        """Get all registered flags.

        Returns:
            List of all QCFlag instances.
        """
        return [
            value for name, value in vars(cls).items()
            if isinstance(value, QCFlag)
        ]

    @classmethod
    def get_by_category(cls, category: FlagCategory) -> list[QCFlag]:
        """Get all flags in a category.

        Args:
            category: The category to filter by.

        Returns:
            List of QCFlag instances in the category.
        """
        return [flag for flag in cls.get_all() if flag.category == category]

    @classmethod
    def get_by_severity(cls, severity: FlagSeverity) -> list[QCFlag]:
        """Get all flags with a severity level.

        Args:
            severity: The severity to filter by.

        Returns:
            List of QCFlag instances with the severity.
        """
        return [flag for flag in cls.get_all() if flag.severity == severity]

    @classmethod
    def get_by_code(cls, code: str) -> QCFlag | None:
        """Get a flag by its code.

        Args:
            code: The flag code.

        Returns:
            The QCFlag or None if not found.
        """
        for flag in cls.get_all():
            if flag.code == code:
                return flag
        return None


# =============================================================================
# GeneQC Class
# =============================================================================


@attrs.define
class GeneQC:
    """Aggregated QC information for a single gene.

    This class collects all QC flags and metrics for a gene from
    various modules (confidence, splice, homology) and computes
    an overall tier classification.

    Attributes:
        gene_id: The gene identifier.
        flags: Set of QC flags assigned to this gene.
        tier: Classification tier (high, medium, low, reject).
        confidence_score: Overall confidence score (0-1).
        splice_score: Splice site quality score (0-1).
        homology_score: Homology validation score (0-1).
        confidence_metrics: Detailed confidence metrics dict.
        splice_metrics: Detailed splice metrics dict.
        homology_metrics: Detailed homology metrics dict.
        notes: Additional notes or comments.

    Example:
        >>> gene_qc = GeneQC(gene_id="gene001")
        >>> gene_qc.confidence_score = 0.65
        >>> gene_qc.add_flag(Flags.LOW_CONFIDENCE)
        >>> gene_qc.classify_tier()
        >>> print(gene_qc.tier)
        'medium'
    """

    gene_id: str
    flags: set[QCFlag] = attrs.Factory(set)
    tier: str = "unclassified"

    # Module scores (0-1 scale)
    confidence_score: float | None = None
    splice_score: float | None = None
    homology_score: float | None = None

    # Detailed metrics from each module
    confidence_metrics: dict[str, Any] = attrs.Factory(dict)
    splice_metrics: dict[str, Any] = attrs.Factory(dict)
    homology_metrics: dict[str, Any] = attrs.Factory(dict)

    # Additional metadata
    notes: list[str] = attrs.Factory(list)

    def add_flag(self, flag: QCFlag) -> None:
        """Add a QC flag to this gene.

        Args:
            flag: The QCFlag to add.
        """
        self.flags.add(flag)

    def add_flags(self, flags: Iterable[QCFlag]) -> None:
        """Add multiple QC flags to this gene.

        Args:
            flags: Iterable of QCFlag objects to add.
        """
        for flag in flags:
            self.flags.add(flag)

    def has_flag(self, flag: QCFlag) -> bool:
        """Check if this gene has a specific flag.

        Args:
            flag: The QCFlag to check for.

        Returns:
            True if the gene has the flag.
        """
        return flag in self.flags

    def has_any_flag(self, flags: Iterable[QCFlag]) -> bool:
        """Check if this gene has any of the specified flags.

        Args:
            flags: Iterable of QCFlag objects to check.

        Returns:
            True if the gene has any of the flags.
        """
        return bool(self.flags.intersection(flags))

    def has_all_flags(self, flags: Iterable[QCFlag]) -> bool:
        """Check if this gene has all of the specified flags.

        Args:
            flags: Iterable of QCFlag objects to check.

        Returns:
            True if the gene has all of the flags.
        """
        return set(flags).issubset(self.flags)

    @property
    def max_severity(self) -> FlagSeverity | None:
        """Get the maximum severity among all flags.

        Returns:
            The highest severity level, or None if no flags.
        """
        if not self.flags:
            return None
        return max(flag.severity for flag in self.flags)

    @property
    def has_critical(self) -> bool:
        """Check if gene has any critical severity flags."""
        return any(f.severity == FlagSeverity.CRITICAL for f in self.flags)

    @property
    def has_error(self) -> bool:
        """Check if gene has any error severity flags."""
        return any(f.severity == FlagSeverity.ERROR for f in self.flags)

    @property
    def has_warning(self) -> bool:
        """Check if gene has any warning severity flags."""
        return any(f.severity == FlagSeverity.WARNING for f in self.flags)

    @property
    def flag_count(self) -> int:
        """Get the total number of flags."""
        return len(self.flags)

    @property
    def flag_codes(self) -> list[str]:
        """Get list of flag codes."""
        return sorted([f.code for f in self.flags])

    def flags_by_category(self, category: FlagCategory) -> set[QCFlag]:
        """Get flags in a specific category.

        Args:
            category: The category to filter by.

        Returns:
            Set of flags in the category.
        """
        return {f for f in self.flags if f.category == category}

    def flags_by_severity(self, severity: FlagSeverity) -> set[QCFlag]:
        """Get flags with a specific severity.

        Args:
            severity: The severity to filter by.

        Returns:
            Set of flags with the severity.
        """
        return {f for f in self.flags if f.severity == severity}

    def classify_tier(
        self,
        high_threshold: float = 0.85,
        medium_threshold: float = 0.70,
        low_threshold: float = 0.50,
    ) -> str:
        """Classify the gene into a quality tier.

        Tier classification logic:
        - reject: Has CRITICAL flags or confidence < low_threshold
        - low: Has ERROR flags or confidence < medium_threshold
        - medium: Has WARNING flags or confidence < high_threshold
        - high: No significant flags and confidence >= high_threshold

        Args:
            high_threshold: Minimum score for high tier.
            medium_threshold: Minimum score for medium tier.
            low_threshold: Minimum score for low tier (below = reject).

        Returns:
            The tier classification string.
        """
        # Check for critical issues first
        if self.has_critical:
            self.tier = "reject"
            return self.tier

        # Check confidence score
        conf = self.confidence_score
        if conf is not None and conf < low_threshold:
            self.tier = "reject"
            return self.tier

        # Check for error-level issues
        if self.has_error:
            self.tier = "low"
            return self.tier

        if conf is not None and conf < medium_threshold:
            self.tier = "low"
            return self.tier

        # Check for warnings
        if self.has_warning:
            self.tier = "medium"
            return self.tier

        if conf is not None and conf < high_threshold:
            self.tier = "medium"
            return self.tier

        # No significant issues
        self.tier = "high"
        return self.tier

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary representation.

        Returns:
            Dictionary with all gene QC information.
        """
        return {
            "gene_id": self.gene_id,
            "tier": self.tier,
            "flag_count": self.flag_count,
            "flag_codes": self.flag_codes,
            "max_severity": self.max_severity.value if self.max_severity else None,
            "confidence_score": self.confidence_score,
            "splice_score": self.splice_score,
            "homology_score": self.homology_score,
            "confidence_metrics": self.confidence_metrics,
            "splice_metrics": self.splice_metrics,
            "homology_metrics": self.homology_metrics,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "GeneQC":
        """Create GeneQC from dictionary.

        Args:
            data: Dictionary with gene QC data.

        Returns:
            New GeneQC instance.
        """
        gene_qc = cls(
            gene_id=data["gene_id"],
            tier=data.get("tier", "unclassified"),
            confidence_score=data.get("confidence_score"),
            splice_score=data.get("splice_score"),
            homology_score=data.get("homology_score"),
            confidence_metrics=data.get("confidence_metrics", {}),
            splice_metrics=data.get("splice_metrics", {}),
            homology_metrics=data.get("homology_metrics", {}),
            notes=data.get("notes", []),
        )

        # Restore flags from codes
        for code in data.get("flag_codes", []):
            flag = Flags.get_by_code(code)
            if flag:
                gene_qc.add_flag(flag)

        return gene_qc


# =============================================================================
# Utility Functions
# =============================================================================


def summarize_flags(flags: Iterable[QCFlag]) -> dict[str, int]:
    """Summarize flags by category.

    Args:
        flags: Iterable of QC flags.

    Returns:
        Dict mapping category name to count.
    """
    summary: dict[str, int] = {}
    for flag in flags:
        cat = flag.category.value
        summary[cat] = summary.get(cat, 0) + 1
    return summary


def max_severity(flags: Iterable[QCFlag]) -> FlagSeverity | None:
    """Get the maximum severity from a collection of flags.

    Args:
        flags: Iterable of QC flags.

    Returns:
        Maximum severity, or None if empty.
    """
    flag_list = list(flags)
    if not flag_list:
        return None
    return max(f.severity for f in flag_list)


def filter_by_severity(
    flags: Iterable[QCFlag],
    min_severity: FlagSeverity,
) -> list[QCFlag]:
    """Filter flags by minimum severity.

    Args:
        flags: Iterable of QC flags.
        min_severity: Minimum severity to include.

    Returns:
        List of flags at or above the minimum severity.
    """
    return [f for f in flags if f.severity >= min_severity]


def filter_by_category(
    flags: Iterable[QCFlag],
    categories: Iterable[FlagCategory],
) -> list[QCFlag]:
    """Filter flags by categories.

    Args:
        flags: Iterable of QC flags.
        categories: Categories to include.

    Returns:
        List of flags in the specified categories.
    """
    cat_set = set(categories)
    return [f for f in flags if f.category in cat_set]
