"""Isoform selection and prioritization.

This module provides tools for selecting representative transcripts
and prioritizing isoforms based on various criteria:

- Protein-coding potential
- Expression level
- Evolutionary conservation
- Functional annotation

Example:
    >>> from helixforge.isoforms.select import select_representative
    >>> primary = select_representative(gene.transcripts)

TODO:
    - Implement representative selection
    - Add scoring criteria
    - Support for expression-based selection
    - APPRIS-like principal isoform selection
    - Handle non-coding transcripts
"""

from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from helixforge.core.models import GeneModel, Transcript

# =============================================================================
# Data Structures
# =============================================================================


class IsoformScore(NamedTuple):
    """Composite score for isoform ranking.

    Attributes:
        transcript_id: Transcript identifier.
        overall_score: Combined ranking score.
        coding_score: Protein-coding potential score.
        expression_score: Expression-based score.
        conservation_score: Conservation-based score.
        structure_score: Gene structure score.
        is_principal: Whether this is the principal isoform.
    """

    transcript_id: str
    overall_score: float
    coding_score: float
    expression_score: float
    conservation_score: float
    structure_score: float
    is_principal: bool


class SelectionCriteria(NamedTuple):
    """Criteria weights for isoform selection.

    Attributes:
        coding_weight: Weight for coding potential.
        expression_weight: Weight for expression level.
        conservation_weight: Weight for conservation.
        structure_weight: Weight for gene structure quality.
        prefer_longest_cds: Prefer longer CDS when tied.
    """

    coding_weight: float = 0.3
    expression_weight: float = 0.3
    conservation_weight: float = 0.2
    structure_weight: float = 0.2
    prefer_longest_cds: bool = True


# =============================================================================
# Selection Functions
# =============================================================================


def select_representative(
    transcripts: list["Transcript"],
    criteria: SelectionCriteria | None = None,
) -> "Transcript":
    """Select the representative/principal transcript.

    Args:
        transcripts: List of transcripts for a gene.
        criteria: Selection criteria weights.

    Returns:
        The selected representative transcript.

    Raises:
        ValueError: If transcripts list is empty.
    """
    # TODO: Implement selection
    raise NotImplementedError("select_representative not yet implemented")


def rank_isoforms(
    transcripts: list["Transcript"],
    criteria: SelectionCriteria | None = None,
) -> list[IsoformScore]:
    """Rank all isoforms by composite score.

    Args:
        transcripts: List of transcripts to rank.
        criteria: Selection criteria weights.

    Returns:
        List of IsoformScore, sorted by overall_score descending.
    """
    # TODO: Implement ranking
    raise NotImplementedError("rank_isoforms not yet implemented")


def filter_isoforms(
    transcripts: list["Transcript"],
    min_cds_length: int = 100,
    require_start_codon: bool = True,
    require_stop_codon: bool = True,
) -> list["Transcript"]:
    """Filter isoforms by basic quality criteria.

    Args:
        transcripts: Transcripts to filter.
        min_cds_length: Minimum CDS length in nucleotides.
        require_start_codon: Require canonical start codon.
        require_stop_codon: Require stop codon.

    Returns:
        Filtered list of transcripts.
    """
    # TODO: Implement filtering
    raise NotImplementedError("filter_isoforms not yet implemented")


# =============================================================================
# Scoring Functions
# =============================================================================


def score_coding_potential(transcript: "Transcript") -> float:
    """Score protein-coding potential of a transcript.

    Args:
        transcript: Transcript to score.

    Returns:
        Coding potential score (0.0 to 1.0).
    """
    # TODO: Implement scoring
    raise NotImplementedError("score_coding_potential not yet implemented")


def score_structure(transcript: "Transcript") -> float:
    """Score gene structure quality.

    Considers:
    - Canonical splice sites
    - Exon/intron lengths
    - Start/stop codons

    Args:
        transcript: Transcript to score.

    Returns:
        Structure score (0.0 to 1.0).
    """
    # TODO: Implement scoring
    raise NotImplementedError("score_structure not yet implemented")


def collapse_redundant(
    transcripts: list["Transcript"],
    identity_threshold: float = 0.99,
) -> list["Transcript"]:
    """Collapse nearly identical transcripts.

    Args:
        transcripts: Transcripts to collapse.
        identity_threshold: Minimum identity to collapse.

    Returns:
        List of unique transcripts.
    """
    # TODO: Implement collapsing
    raise NotImplementedError("collapse_redundant not yet implemented")
