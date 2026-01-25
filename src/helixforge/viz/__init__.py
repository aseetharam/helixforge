"""Visualization tools for HelixForge.

This module provides visualization tools for gene models and evidence:

- Locus-level views with evidence tracks
- Confidence visualization (gene-level and genome-wide)
- Genome-wide statistics plots
- Interactive visualization (Dash/Panel)

Example:
    >>> from helixforge.viz import plot_locus, plot_gene_confidence
    >>> fig = plot_locus(gene, evidence)
    >>> fig.savefig("locus.png")
    >>>
    >>> # Confidence visualization
    >>> fig = plot_gene_confidence(gene, confidence, region_conf)
    >>> fig = plot_confidence_distribution(scores)
"""

from helixforge.viz.genome import (
    plot_confidence_by_scaffold,
    plot_confidence_distribution,
    plot_confidence_metrics_comparison,
    plot_flag_distribution,
)
from helixforge.viz.locus import (
    LocusPlotter,
    plot_gene_confidence,
    plot_gene_confidence_batch,
    plot_locus,
)

__all__ = [
    # Locus-level
    "LocusPlotter",
    "plot_locus",
    "plot_gene_confidence",
    "plot_gene_confidence_batch",
    # Genome-wide
    "plot_confidence_distribution",
    "plot_confidence_by_scaffold",
    "plot_confidence_metrics_comparison",
    "plot_flag_distribution",
]
