"""QC flags for gene model issues.

This module defines quality control flags that indicate potential
issues with gene predictions. Flags are organized by severity and
category.

Flag Categories:
    - Structure: Gene structure issues (splice sites, frameshifts)
    - Sequence: Sequence-based issues (start/stop codons)
    - Evidence: Evidence support issues
    - Homology: Homology-based issues

Example:
    >>> from helixforge.qc.flags import apply_qc_flags, QCFlag
    >>> flags = apply_qc_flags(gene, genome)
    >>> for flag in flags:
    ...     print(f"{flag.code}: {flag.description}")

TODO:
    - Define all QC flags
    - Implement flag detection functions
    - Add severity levels
    - Create flag aggregation utilities
"""

from enum import Enum
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from helixforge.core.models import GeneModel
    from helixforge.io.fasta import FastaReader

# =============================================================================
# Enums
# =============================================================================


class FlagSeverity(Enum):
    """Severity levels for QC flags."""

    INFO = "info"  # Informational, not necessarily an issue
    WARNING = "warning"  # Potential issue, may need review
    ERROR = "error"  # Definite issue, likely needs correction
    CRITICAL = "critical"  # Severe issue, gene may be invalid


class FlagCategory(Enum):
    """Categories of QC flags."""

    STRUCTURE = "structure"  # Gene structure issues
    SEQUENCE = "sequence"  # Sequence-based issues
    EVIDENCE = "evidence"  # Evidence support issues
    HOMOLOGY = "homology"  # Homology-based issues
    LENGTH = "length"  # Length-related issues


# =============================================================================
# Data Structures
# =============================================================================


class QCFlag(NamedTuple):
    """A quality control flag for a gene.

    Attributes:
        code: Short flag code (e.g., "NO_START").
        category: Flag category.
        severity: Severity level.
        description: Human-readable description.
        details: Additional details.
        position: Genomic position of issue (if applicable).
    """

    code: str
    category: FlagCategory
    severity: FlagSeverity
    description: str
    details: str = ""
    position: int | None = None


# =============================================================================
# Flag Definitions
# =============================================================================

# Structure flags
FLAG_NON_CANONICAL_SPLICE = QCFlag(
    code="NON_CANONICAL_SPLICE",
    category=FlagCategory.STRUCTURE,
    severity=FlagSeverity.WARNING,
    description="Non-canonical splice site detected",
)

FLAG_SHORT_INTRON = QCFlag(
    code="SHORT_INTRON",
    category=FlagCategory.STRUCTURE,
    severity=FlagSeverity.WARNING,
    description="Intron shorter than minimum length",
)

FLAG_LONG_INTRON = QCFlag(
    code="LONG_INTRON",
    category=FlagCategory.STRUCTURE,
    severity=FlagSeverity.INFO,
    description="Intron longer than typical maximum",
)

# Sequence flags
FLAG_NO_START_CODON = QCFlag(
    code="NO_START_CODON",
    category=FlagCategory.SEQUENCE,
    severity=FlagSeverity.ERROR,
    description="No start codon found at CDS start",
)

FLAG_NO_STOP_CODON = QCFlag(
    code="NO_STOP_CODON",
    category=FlagCategory.SEQUENCE,
    severity=FlagSeverity.ERROR,
    description="No stop codon found at CDS end",
)

FLAG_INTERNAL_STOP = QCFlag(
    code="INTERNAL_STOP",
    category=FlagCategory.SEQUENCE,
    severity=FlagSeverity.CRITICAL,
    description="Internal stop codon in CDS",
)

FLAG_FRAMESHIFT = QCFlag(
    code="FRAMESHIFT",
    category=FlagCategory.SEQUENCE,
    severity=FlagSeverity.ERROR,
    description="Frameshift detected in CDS",
)

# Evidence flags
FLAG_NO_EVIDENCE = QCFlag(
    code="NO_EVIDENCE",
    category=FlagCategory.EVIDENCE,
    severity=FlagSeverity.WARNING,
    description="No RNA-seq evidence supports this gene",
)

FLAG_LOW_COVERAGE = QCFlag(
    code="LOW_COVERAGE",
    category=FlagCategory.EVIDENCE,
    severity=FlagSeverity.INFO,
    description="Low RNA-seq coverage",
)

FLAG_UNSUPPORTED_JUNCTION = QCFlag(
    code="UNSUPPORTED_JUNCTION",
    category=FlagCategory.EVIDENCE,
    severity=FlagSeverity.WARNING,
    description="Splice junction has no RNA-seq support",
)

# Homology flags
FLAG_NO_HOMOLOGY = QCFlag(
    code="NO_HOMOLOGY",
    category=FlagCategory.HOMOLOGY,
    severity=FlagSeverity.INFO,
    description="No homology to known proteins",
)

FLAG_PARTIAL_HOMOLOGY = QCFlag(
    code="PARTIAL_HOMOLOGY",
    category=FlagCategory.HOMOLOGY,
    severity=FlagSeverity.WARNING,
    description="Partial homology coverage",
)

# Length flags
FLAG_SHORT_CDS = QCFlag(
    code="SHORT_CDS",
    category=FlagCategory.LENGTH,
    severity=FlagSeverity.WARNING,
    description="CDS shorter than minimum length",
)

FLAG_SHORT_EXON = QCFlag(
    code="SHORT_EXON",
    category=FlagCategory.LENGTH,
    severity=FlagSeverity.INFO,
    description="Exon shorter than typical minimum",
)


# =============================================================================
# Flag Application Functions
# =============================================================================


def apply_qc_flags(
    gene: "GeneModel",
    genome: "FastaReader",
    check_evidence: bool = True,
) -> list[QCFlag]:
    """Apply all QC checks and return flags.

    Args:
        gene: Gene model to check.
        genome: Genome sequence reader.
        check_evidence: Whether to check evidence flags.

    Returns:
        List of QCFlag objects.
    """
    # TODO: Implement flag application
    raise NotImplementedError("apply_qc_flags not yet implemented")


def check_start_codon(
    gene: "GeneModel",
    genome: "FastaReader",
) -> QCFlag | None:
    """Check for valid start codon.

    Args:
        gene: Gene model to check.
        genome: Genome sequence reader.

    Returns:
        QCFlag if issue found, None otherwise.
    """
    # TODO: Implement check
    raise NotImplementedError("check_start_codon not yet implemented")


def check_stop_codon(
    gene: "GeneModel",
    genome: "FastaReader",
) -> QCFlag | None:
    """Check for valid stop codon.

    Args:
        gene: Gene model to check.
        genome: Genome sequence reader.

    Returns:
        QCFlag if issue found, None otherwise.
    """
    # TODO: Implement check
    raise NotImplementedError("check_stop_codon not yet implemented")


def check_splice_sites(
    gene: "GeneModel",
    genome: "FastaReader",
) -> list[QCFlag]:
    """Check all splice sites for canonical dinucleotides.

    Args:
        gene: Gene model to check.
        genome: Genome sequence reader.

    Returns:
        List of QCFlags for non-canonical sites.
    """
    # TODO: Implement check
    raise NotImplementedError("check_splice_sites not yet implemented")


def check_internal_stops(
    gene: "GeneModel",
    genome: "FastaReader",
) -> list[QCFlag]:
    """Check for internal stop codons.

    Args:
        gene: Gene model to check.
        genome: Genome sequence reader.

    Returns:
        List of QCFlags for internal stops.
    """
    # TODO: Implement check
    raise NotImplementedError("check_internal_stops not yet implemented")


# =============================================================================
# Flag Summary Functions
# =============================================================================


def summarize_flags(flags: list[QCFlag]) -> dict[str, int]:
    """Summarize flags by category.

    Args:
        flags: List of QC flags.

    Returns:
        Dict mapping category to count.
    """
    summary: dict[str, int] = {}
    for flag in flags:
        cat = flag.category.value
        summary[cat] = summary.get(cat, 0) + 1
    return summary


def max_severity(flags: list[QCFlag]) -> FlagSeverity | None:
    """Get the maximum severity from a list of flags.

    Args:
        flags: List of QC flags.

    Returns:
        Maximum severity, or None if list is empty.
    """
    if not flags:
        return None

    severity_order = [
        FlagSeverity.INFO,
        FlagSeverity.WARNING,
        FlagSeverity.ERROR,
        FlagSeverity.CRITICAL,
    ]

    max_idx = 0
    for flag in flags:
        idx = severity_order.index(flag.severity)
        max_idx = max(max_idx, idx)

    return severity_order[max_idx]
