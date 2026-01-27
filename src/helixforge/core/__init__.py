"""Core refinement logic for HelixForge.

This module contains the fundamental algorithms and data structures
for refining gene predictions:

- Confidence scoring
- Splice site analysis
- Gene boundary refinement
- RNA-seq evidence scoring
- Merge/split detection

Example:
    >>> from helixforge.core.evidence import EvidenceScorer, EvidenceScorerConfig
    >>> from helixforge.core.splice import SpliceRefiner
"""

from helixforge.core.evidence import (
    EvidenceLevel,
    EvidenceScore,
    EvidenceScorer,
    EvidenceScorerConfig,
    ExonEvidence,
    JunctionEvidence,
    summarize_evidence_scores,
)

__all__: list[str] = [
    # Evidence scoring
    "EvidenceLevel",
    "EvidenceScore",
    "EvidenceScorer",
    "EvidenceScorerConfig",
    "ExonEvidence",
    "JunctionEvidence",
    "summarize_evidence_scores",
]
