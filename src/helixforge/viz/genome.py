"""Genome-wide visualization and statistics.

This module provides genome-wide visualization and statistics plots:

- Gene density plots
- Confidence score distributions
- QC flag summaries
- Chromosome ideograms

Example:
    >>> from helixforge.viz.genome import plot_statistics, plot_confidence_distribution
    >>> fig = plot_statistics(genes)
    >>> fig.savefig("genome_stats.png")
    >>>
    >>> # Confidence distribution
    >>> from helixforge.core.confidence import GeneConfidence
    >>> fig = plot_confidence_distribution(confidence_scores, "distribution.html")
"""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import numpy as np

if TYPE_CHECKING:
    import matplotlib.figure
    import plotly.graph_objects as go

    from helixforge.core.confidence import GeneConfidence
    from helixforge.core.models import GeneModel

logger = logging.getLogger(__name__)

# Color schemes
CONFIDENCE_COLORS = {
    "high": "#2ecc71",  # Green
    "medium": "#f39c12",  # Orange/Yellow
    "low": "#e74c3c",  # Red
}

# =============================================================================
# Statistics Plots
# =============================================================================


def plot_statistics(
    genes: list["GeneModel"],
    output: Path | str | None = None,
) -> "matplotlib.figure.Figure":
    """Generate summary statistics plots.

    Creates a multi-panel figure with:
    - Confidence score histogram
    - Gene length distribution
    - Exon count distribution
    - QC flag summary

    Args:
        genes: Gene models to summarize.
        output: If provided, save to this path.

    Returns:
        Matplotlib figure.
    """
    # TODO: Implement statistics plotting
    raise NotImplementedError("plot_statistics not yet implemented")


def plot_confidence_distribution(
    scores: Iterable["GeneConfidence"],
    output_path: Path | str | None = None,
    format: str = "html",
    title: str = "Confidence Score Distribution",
) -> Any:
    """Plot confidence score distribution as histogram and pie chart.

    Creates a two-panel figure:
    - Left: Histogram of overall confidence scores with density curve
    - Right: Pie chart of confidence class distribution

    Args:
        scores: Iterable of GeneConfidence objects.
        output_path: If provided, save figure to this path.
        format: Output format ("html" for Plotly, "png"/"pdf" for Matplotlib).
        title: Plot title.

    Returns:
        Plotly Figure (if format="html") or Matplotlib Figure.

    Example:
        >>> from helixforge.core.confidence import ConfidenceCalculator
        >>> scores = list(calc.score_genes_parallel(genes))
        >>> fig = plot_confidence_distribution(scores, "dist.html")
    """
    # Convert to list to allow multiple passes
    scores_list = list(scores)

    if format == "html":
        return _plot_confidence_distribution_plotly(scores_list, output_path, title)
    else:
        return _plot_confidence_distribution_matplotlib(
            scores_list, output_path, format, title
        )


def _plot_confidence_distribution_plotly(
    scores: list["GeneConfidence"],
    output_path: Path | str | None,
    title: str,
) -> "go.Figure":
    """Create interactive Plotly distribution plot."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Extract data
    overall_scores = [s.overall_score for s in scores]
    class_counts = Counter(s.confidence_class for s in scores)

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=["Score Distribution", "Confidence Classes"],
        specs=[[{"type": "histogram"}, {"type": "pie"}]],
    )

    # Left panel: Histogram
    fig.add_trace(
        go.Histogram(
            x=overall_scores,
            nbinsx=30,
            name="Overall Score",
            marker_color="#3498db",
            opacity=0.7,
        ),
        row=1,
        col=1,
    )

    # Add vertical lines for thresholds
    fig.add_vline(
        x=0.70,
        line_dash="dash",
        line_color="orange",
        annotation_text="Medium",
        row=1,
        col=1,
    )
    fig.add_vline(
        x=0.85,
        line_dash="dash",
        line_color="green",
        annotation_text="High",
        row=1,
        col=1,
    )

    # Right panel: Pie chart
    labels = ["High", "Medium", "Low"]
    values = [
        class_counts.get("high", 0),
        class_counts.get("medium", 0),
        class_counts.get("low", 0),
    ]
    colors = [
        CONFIDENCE_COLORS["high"],
        CONFIDENCE_COLORS["medium"],
        CONFIDENCE_COLORS["low"],
    ]

    fig.add_trace(
        go.Pie(
            labels=labels,
            values=values,
            marker_colors=colors,
            textinfo="percent+value",
            name="Classes",
        ),
        row=1,
        col=2,
    )

    # Calculate statistics
    mean_score = np.mean(overall_scores)
    median_score = np.median(overall_scores)
    std_score = np.std(overall_scores)

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{title}</b><br>"
                f"<sup>N={len(scores)} | Mean={mean_score:.3f} | "
                f"Median={median_score:.3f} | Std={std_score:.3f}</sup>"
            ),
            x=0.5,
        ),
        height=500,
        showlegend=False,
    )

    fig.update_xaxes(title_text="Overall Score", range=[0, 1], row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)

    if output_path:
        output_path = Path(output_path)
        fig.write_html(str(output_path))
        logger.info(f"Saved confidence distribution to {output_path}")

    return fig


def _plot_confidence_distribution_matplotlib(
    scores: list["GeneConfidence"],
    output_path: Path | str | None,
    format: str,
    title: str,
) -> "matplotlib.figure.Figure":
    """Create static Matplotlib distribution plot."""
    import matplotlib.pyplot as plt

    overall_scores = [s.overall_score for s in scores]
    class_counts = Counter(s.confidence_class for s in scores)

    fig, (ax_hist, ax_pie) = plt.subplots(1, 2, figsize=(12, 5))

    # Histogram
    ax_hist.hist(overall_scores, bins=30, color="#3498db", alpha=0.7, edgecolor="black")
    ax_hist.axvline(x=0.70, color="orange", linestyle="--", label="Medium threshold")
    ax_hist.axvline(x=0.85, color="green", linestyle="--", label="High threshold")
    ax_hist.set_xlabel("Overall Score")
    ax_hist.set_ylabel("Count")
    ax_hist.set_xlim(0, 1)
    ax_hist.set_title("Score Distribution")
    ax_hist.legend()

    # Pie chart
    labels = ["High", "Medium", "Low"]
    values = [
        class_counts.get("high", 0),
        class_counts.get("medium", 0),
        class_counts.get("low", 0),
    ]
    colors = [
        CONFIDENCE_COLORS["high"],
        CONFIDENCE_COLORS["medium"],
        CONFIDENCE_COLORS["low"],
    ]

    ax_pie.pie(
        values,
        labels=labels,
        colors=colors,
        autopct="%1.1f%%",
        startangle=90,
    )
    ax_pie.set_title("Confidence Classes")

    # Statistics
    mean_score = np.mean(overall_scores)
    median_score = np.median(overall_scores)
    std_score = np.std(overall_scores)

    fig.suptitle(
        f"{title}\nN={len(scores)} | Mean={mean_score:.3f} | "
        f"Median={median_score:.3f} | Std={std_score:.3f}",
        fontsize=12,
        fontweight="bold",
    )

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
        logger.info(f"Saved confidence distribution to {output_path}")

    return fig


def plot_length_distributions(
    genes: list["GeneModel"],
    ax: Any | None = None,
) -> Any:
    """Plot gene/exon/intron length distributions.

    Args:
        genes: Gene models to analyze.
        ax: Optional matplotlib axes.

    Returns:
        Matplotlib axes.
    """
    # TODO: Implement length distributions
    raise NotImplementedError("plot_length_distributions not yet implemented")


def plot_flag_summary(
    genes: list["GeneModel"],
    ax: Any | None = None,
) -> Any:
    """Plot QC flag summary.

    Args:
        genes: Gene models with QC flags.
        ax: Optional matplotlib axes.

    Returns:
        Matplotlib axes.
    """
    # TODO: Implement flag summary
    raise NotImplementedError("plot_flag_summary not yet implemented")


# =============================================================================
# Genome View
# =============================================================================


def plot_ideogram(
    genes: list["GeneModel"],
    chromosome_sizes: dict[str, int],
    output: Path | str | None = None,
) -> "matplotlib.figure.Figure":
    """Plot chromosome ideogram with gene positions.

    Args:
        genes: Gene models to display.
        chromosome_sizes: Dict mapping chromosome to size.
        output: If provided, save to this path.

    Returns:
        Matplotlib figure.
    """
    # TODO: Implement ideogram
    raise NotImplementedError("plot_ideogram not yet implemented")


def plot_gene_density(
    genes: list["GeneModel"],
    chromosome: str,
    chromosome_size: int,
    window_size: int = 100000,
    ax: Any | None = None,
) -> Any:
    """Plot gene density along a chromosome.

    Args:
        genes: Gene models to count.
        chromosome: Chromosome to plot.
        chromosome_size: Size of the chromosome.
        window_size: Window size for density calculation.
        ax: Optional matplotlib axes.

    Returns:
        Matplotlib axes.
    """
    # TODO: Implement gene density
    raise NotImplementedError("plot_gene_density not yet implemented")


# =============================================================================
# Comparative Plots
# =============================================================================


def plot_comparison_stats(
    gene_sets: dict[str, list["GeneModel"]],
    output: Path | str | None = None,
) -> "matplotlib.figure.Figure":
    """Plot comparison statistics for multiple gene sets.

    Args:
        gene_sets: Dict mapping set name to genes.
        output: If provided, save to this path.

    Returns:
        Matplotlib figure.
    """
    # TODO: Implement comparison stats
    raise NotImplementedError("plot_comparison_stats not yet implemented")


# =============================================================================
# Confidence-specific Genome Plots
# =============================================================================


def plot_confidence_by_scaffold(
    scores: Iterable["GeneConfidence"],
    output_path: Path | str | None = None,
    format: str = "html",
    max_scaffolds: int = 30,
    title: str = "Confidence by Scaffold",
) -> Any:
    """Plot confidence score distribution per scaffold as box plots.

    Creates box plots showing confidence score distribution for each scaffold,
    sorted by median score (lowest to highest).

    Args:
        scores: Iterable of GeneConfidence objects.
        output_path: If provided, save figure to this path.
        format: Output format ("html" for Plotly, "png"/"pdf" for Matplotlib).
        max_scaffolds: Maximum number of scaffolds to display.
        title: Plot title.

    Returns:
        Plotly Figure (if format="html") or Matplotlib Figure.

    Example:
        >>> scores = list(calc.score_genes_parallel(genes))
        >>> fig = plot_confidence_by_scaffold(scores, "by_scaffold.html")
    """
    # Convert to list and group by scaffold
    scores_list = list(scores)

    if format == "html":
        return _plot_confidence_by_scaffold_plotly(
            scores_list, output_path, max_scaffolds, title
        )
    else:
        return _plot_confidence_by_scaffold_matplotlib(
            scores_list, output_path, format, max_scaffolds, title
        )


def _plot_confidence_by_scaffold_plotly(
    scores: list["GeneConfidence"],
    output_path: Path | str | None,
    max_scaffolds: int,
    title: str,
) -> "go.Figure":
    """Create interactive Plotly box plot."""
    import plotly.graph_objects as go

    # Group scores by scaffold
    scaffold_scores: dict[str, list[float]] = {}
    for s in scores:
        if s.seqid not in scaffold_scores:
            scaffold_scores[s.seqid] = []
        scaffold_scores[s.seqid].append(s.overall_score)

    # Sort scaffolds by median score
    scaffold_medians = {
        k: np.median(v) for k, v in scaffold_scores.items()
    }
    sorted_scaffolds = sorted(scaffold_medians.keys(), key=lambda x: scaffold_medians[x])

    # Limit to max_scaffolds
    if len(sorted_scaffolds) > max_scaffolds:
        sorted_scaffolds = sorted_scaffolds[:max_scaffolds]
        note = f" (showing {max_scaffolds} lowest-confidence scaffolds)"
    else:
        note = ""

    fig = go.Figure()

    for scaffold in sorted_scaffolds:
        scaffold_data = scaffold_scores[scaffold]
        median = scaffold_medians[scaffold]

        # Color based on median
        if median >= 0.85:
            color = CONFIDENCE_COLORS["high"]
        elif median >= 0.70:
            color = CONFIDENCE_COLORS["medium"]
        else:
            color = CONFIDENCE_COLORS["low"]

        fig.add_trace(
            go.Box(
                y=scaffold_data,
                name=scaffold,
                marker_color=color,
                boxmean=True,
                hoverinfo="y+name",
            )
        )

    # Add threshold lines
    fig.add_hline(
        y=0.70,
        line_dash="dash",
        line_color="orange",
        annotation_text="Medium threshold",
    )
    fig.add_hline(
        y=0.85,
        line_dash="dash",
        line_color="green",
        annotation_text="High threshold",
    )

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b>{note}",
            x=0.5,
        ),
        xaxis_title="Scaffold",
        yaxis_title="Overall Score",
        yaxis=dict(range=[0, 1]),
        showlegend=False,
        height=500,
    )

    if output_path:
        output_path = Path(output_path)
        fig.write_html(str(output_path))
        logger.info(f"Saved scaffold confidence plot to {output_path}")

    return fig


def _plot_confidence_by_scaffold_matplotlib(
    scores: list["GeneConfidence"],
    output_path: Path | str | None,
    format: str,
    max_scaffolds: int,
    title: str,
) -> "matplotlib.figure.Figure":
    """Create static Matplotlib box plot."""
    import matplotlib.pyplot as plt

    # Group scores by scaffold
    scaffold_scores: dict[str, list[float]] = {}
    for s in scores:
        if s.seqid not in scaffold_scores:
            scaffold_scores[s.seqid] = []
        scaffold_scores[s.seqid].append(s.overall_score)

    # Sort scaffolds by median score
    scaffold_medians = {
        k: np.median(v) for k, v in scaffold_scores.items()
    }
    sorted_scaffolds = sorted(scaffold_medians.keys(), key=lambda x: scaffold_medians[x])

    # Limit to max_scaffolds
    if len(sorted_scaffolds) > max_scaffolds:
        sorted_scaffolds = sorted_scaffolds[:max_scaffolds]
        note = f" (showing {max_scaffolds} lowest-confidence scaffolds)"
    else:
        note = ""

    # Prepare data for plotting
    data = [scaffold_scores[s] for s in sorted_scaffolds]

    # Create figure
    fig_width = max(10, len(sorted_scaffolds) * 0.5)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    # Create box plot
    bp = ax.boxplot(data, patch_artist=True)

    # Color boxes based on median
    for i, (patch, scaffold) in enumerate(zip(bp["boxes"], sorted_scaffolds)):
        median = scaffold_medians[scaffold]
        if median >= 0.85:
            color = CONFIDENCE_COLORS["high"]
        elif median >= 0.70:
            color = CONFIDENCE_COLORS["medium"]
        else:
            color = CONFIDENCE_COLORS["low"]
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Add threshold lines
    ax.axhline(y=0.70, color="orange", linestyle="--", label="Medium threshold")
    ax.axhline(y=0.85, color="green", linestyle="--", label="High threshold")

    ax.set_xticks(range(1, len(sorted_scaffolds) + 1))
    ax.set_xticklabels(sorted_scaffolds, rotation=45, ha="right")
    ax.set_xlabel("Scaffold")
    ax.set_ylabel("Overall Score")
    ax.set_ylim(0, 1)
    ax.set_title(f"{title}{note}")
    ax.legend(loc="upper right")

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
        logger.info(f"Saved scaffold confidence plot to {output_path}")

    return fig


def plot_confidence_metrics_comparison(
    scores: Iterable["GeneConfidence"],
    output_path: Path | str | None = None,
    format: str = "html",
    title: str = "Confidence Metrics Comparison",
) -> Any:
    """Plot comparison of different confidence metrics.

    Creates a multi-panel figure showing distributions and correlations
    between different confidence metrics.

    Args:
        scores: Iterable of GeneConfidence objects.
        output_path: If provided, save figure to this path.
        format: Output format ("html" for Plotly, "png"/"pdf" for Matplotlib).
        title: Plot title.

    Returns:
        Plotly Figure (if format="html") or Matplotlib Figure.
    """
    scores_list = list(scores)

    if format == "html":
        return _plot_metrics_comparison_plotly(scores_list, output_path, title)
    else:
        return _plot_metrics_comparison_matplotlib(scores_list, output_path, format, title)


def _plot_metrics_comparison_plotly(
    scores: list["GeneConfidence"],
    output_path: Path | str | None,
    title: str,
) -> "go.Figure":
    """Create interactive Plotly metrics comparison."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Extract metrics
    metrics = {
        "Mean Prob": [s.mean_prob for s in scores],
        "Min Prob": [s.min_prob for s in scores],
        "Median Prob": [s.median_prob for s in scores],
        "Entropy": [s.entropy for s in scores],
        "Boundary": [s.boundary_sharpness for s in scores],
        "Coding": [s.coding_consistency for s in scores],
        "Overall": [s.overall_score for s in scores],
    }

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            "Metric Distributions",
            "Mean vs Min Probability",
            "Entropy vs Overall Score",
            "Boundary vs Coding Consistency",
        ],
    )

    # Panel 1: Box plots of all metrics
    for name, values in metrics.items():
        fig.add_trace(
            go.Box(y=values, name=name, boxmean=True),
            row=1,
            col=1,
        )

    # Panel 2: Mean vs Min scatter
    colors = [CONFIDENCE_COLORS.get(s.confidence_class, "#95a5a6") for s in scores]
    fig.add_trace(
        go.Scatter(
            x=metrics["Mean Prob"],
            y=metrics["Min Prob"],
            mode="markers",
            marker=dict(color=colors, size=5, opacity=0.6),
            name="Genes",
        ),
        row=1,
        col=2,
    )
    fig.update_xaxes(title_text="Mean Prob", row=1, col=2)
    fig.update_yaxes(title_text="Min Prob", row=1, col=2)

    # Panel 3: Entropy vs Overall
    fig.add_trace(
        go.Scatter(
            x=metrics["Entropy"],
            y=metrics["Overall"],
            mode="markers",
            marker=dict(color=colors, size=5, opacity=0.6),
            name="Genes",
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    fig.update_xaxes(title_text="Entropy", row=2, col=1)
    fig.update_yaxes(title_text="Overall Score", row=2, col=1)

    # Panel 4: Boundary vs Coding
    fig.add_trace(
        go.Scatter(
            x=metrics["Boundary"],
            y=metrics["Coding"],
            mode="markers",
            marker=dict(color=colors, size=5, opacity=0.6),
            name="Genes",
            showlegend=False,
        ),
        row=2,
        col=2,
    )
    fig.update_xaxes(title_text="Boundary Sharpness", row=2, col=2)
    fig.update_yaxes(title_text="Coding Consistency", row=2, col=2)

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0.5),
        height=800,
        showlegend=False,
    )

    if output_path:
        output_path = Path(output_path)
        fig.write_html(str(output_path))
        logger.info(f"Saved metrics comparison to {output_path}")

    return fig


def _plot_metrics_comparison_matplotlib(
    scores: list["GeneConfidence"],
    output_path: Path | str | None,
    format: str,
    title: str,
) -> "matplotlib.figure.Figure":
    """Create static Matplotlib metrics comparison."""
    import matplotlib.pyplot as plt

    # Extract metrics
    metrics = {
        "Mean Prob": [s.mean_prob for s in scores],
        "Min Prob": [s.min_prob for s in scores],
        "Median Prob": [s.median_prob for s in scores],
        "Entropy": [s.entropy for s in scores],
        "Boundary": [s.boundary_sharpness for s in scores],
        "Coding": [s.coding_consistency for s in scores],
        "Overall": [s.overall_score for s in scores],
    }

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Colors for scatter plots
    colors = [CONFIDENCE_COLORS.get(s.confidence_class, "#95a5a6") for s in scores]

    # Panel 1: Box plots
    ax1 = axes[0, 0]
    ax1.boxplot(
        list(metrics.values()),
        labels=list(metrics.keys()),
    )
    ax1.set_ylabel("Value")
    ax1.set_title("Metric Distributions")
    ax1.tick_params(axis="x", rotation=45)

    # Panel 2: Mean vs Min
    ax2 = axes[0, 1]
    ax2.scatter(metrics["Mean Prob"], metrics["Min Prob"], c=colors, alpha=0.6, s=10)
    ax2.set_xlabel("Mean Prob")
    ax2.set_ylabel("Min Prob")
    ax2.set_title("Mean vs Min Probability")

    # Panel 3: Entropy vs Overall
    ax3 = axes[1, 0]
    ax3.scatter(metrics["Entropy"], metrics["Overall"], c=colors, alpha=0.6, s=10)
    ax3.set_xlabel("Entropy")
    ax3.set_ylabel("Overall Score")
    ax3.set_title("Entropy vs Overall Score")

    # Panel 4: Boundary vs Coding
    ax4 = axes[1, 1]
    ax4.scatter(metrics["Boundary"], metrics["Coding"], c=colors, alpha=0.6, s=10)
    ax4.set_xlabel("Boundary Sharpness")
    ax4.set_ylabel("Coding Consistency")
    ax4.set_title("Boundary vs Coding Consistency")

    fig.suptitle(title, fontsize=14, fontweight="bold")
    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
        logger.info(f"Saved metrics comparison to {output_path}")

    return fig


def plot_flag_distribution(
    scores: Iterable["GeneConfidence"],
    output_path: Path | str | None = None,
    format: str = "html",
    title: str = "Confidence Flag Distribution",
) -> Any:
    """Plot distribution of confidence flags.

    Creates a bar chart showing the frequency of different confidence flags.

    Args:
        scores: Iterable of GeneConfidence objects.
        output_path: If provided, save figure to this path.
        format: Output format ("html" for Plotly, "png"/"pdf" for Matplotlib).
        title: Plot title.

    Returns:
        Plotly Figure (if format="html") or Matplotlib Figure.
    """
    scores_list = list(scores)

    # Count flags
    flag_counts: Counter[str] = Counter()
    for s in scores_list:
        for flag in s.flags:
            flag_counts[flag] += 1

    if not flag_counts:
        logger.warning("No flags found in scores")

    if format == "html":
        import plotly.graph_objects as go

        # Sort by count
        sorted_flags = flag_counts.most_common()
        flags = [f[0] for f in sorted_flags]
        counts = [f[1] for f in sorted_flags]

        fig = go.Figure(
            go.Bar(
                x=flags,
                y=counts,
                marker_color="#e74c3c",
            )
        )

        fig.update_layout(
            title=dict(text=f"<b>{title}</b>", x=0.5),
            xaxis_title="Flag",
            yaxis_title="Count",
            height=400,
        )

        if output_path:
            output_path = Path(output_path)
            fig.write_html(str(output_path))
            logger.info(f"Saved flag distribution to {output_path}")

        return fig
    else:
        import matplotlib.pyplot as plt

        sorted_flags = flag_counts.most_common()
        flags = [f[0] for f in sorted_flags]
        counts = [f[1] for f in sorted_flags]

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(flags, counts, color="#e74c3c")
        ax.set_xlabel("Flag")
        ax.set_ylabel("Count")
        ax.set_title(title)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        if output_path:
            output_path = Path(output_path)
            fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
            logger.info(f"Saved flag distribution to {output_path}")

        return fig


# =============================================================================
# Splice Refinement Summary Plots
# =============================================================================


def plot_splice_summary(
    reports: list[Any],  # list[GeneSpliceReport]
    output_path: Path | str | None = None,
    format: str = "html",
    title: str = "Splice Refinement Summary",
) -> Any:
    """Genome-wide splice refinement summary.

    Creates a multi-panel figure showing:
    1. Support ratio histogram
    2. Corrections per gene histogram
    3. Splice type pie chart
    4. Correction distance distribution

    Args:
        reports: List of GeneSpliceReport objects.
        output_path: If provided, save figure to this path.
        format: Output format ("html" for Plotly, "png"/"pdf" for Matplotlib).
        title: Plot title.

    Returns:
        Plotly Figure (if format="html") or Matplotlib Figure.

    Example:
        >>> reports = [refiner.refine_gene(g)[1] for g in genes]
        >>> fig = plot_splice_summary(reports, "splice_summary.html")
    """
    if format == "html":
        return _plot_splice_summary_plotly(reports, output_path, title)
    else:
        return _plot_splice_summary_matplotlib(reports, output_path, format, title)


def _plot_splice_summary_plotly(
    reports: list[Any],
    output_path: Path | str | None,
    title: str,
) -> Any:
    """Create interactive Plotly splice summary."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Extract data
    support_ratios = [r.support_ratio for r in reports]
    n_corrected = [r.n_corrected for r in reports]
    n_canonical = sum(r.n_canonical for r in reports)
    n_noncanonical = sum(r.n_noncanonical for r in reports)

    # Extract correction distances
    correction_distances = []
    for r in reports:
        for c in r.corrections:
            correction_distances.append(abs(c.shift))

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=[
            "Support Ratio Distribution",
            "Corrections per Gene",
            "Splice Site Types",
            "Correction Distance Distribution",
        ],
        specs=[
            [{"type": "histogram"}, {"type": "histogram"}],
            [{"type": "pie"}, {"type": "histogram"}],
        ],
    )

    # Panel 1: Support ratio histogram
    fig.add_trace(
        go.Histogram(
            x=support_ratios,
            nbinsx=20,
            name="Support Ratio",
            marker_color="#3498db",
        ),
        row=1,
        col=1,
    )
    fig.update_xaxes(title_text="Support Ratio", range=[0, 1], row=1, col=1)
    fig.update_yaxes(title_text="Gene Count", row=1, col=1)

    # Panel 2: Corrections per gene histogram
    fig.add_trace(
        go.Histogram(
            x=n_corrected,
            nbinsx=max(1, max(n_corrected) if n_corrected else 1),
            name="Corrections",
            marker_color="#27ae60",
        ),
        row=1,
        col=2,
    )
    fig.update_xaxes(title_text="Corrections per Gene", row=1, col=2)
    fig.update_yaxes(title_text="Gene Count", row=1, col=2)

    # Panel 3: Splice type pie chart
    fig.add_trace(
        go.Pie(
            labels=["Canonical (GT-AG/GC-AG)", "Non-canonical"],
            values=[n_canonical, n_noncanonical],
            marker_colors=["#27ae60", "#e74c3c"],
            textinfo="percent+value",
            name="Splice Types",
        ),
        row=2,
        col=1,
    )

    # Panel 4: Correction distance histogram
    if correction_distances:
        fig.add_trace(
            go.Histogram(
                x=correction_distances,
                nbinsx=15,
                name="Distances",
                marker_color="#f39c12",
            ),
            row=2,
            col=2,
        )
    fig.update_xaxes(title_text="Correction Distance (bp)", row=2, col=2)
    fig.update_yaxes(title_text="Count", row=2, col=2)

    # Calculate summary stats
    total_genes = len(reports)
    total_introns = sum(r.n_introns for r in reports)
    total_supported = sum(r.n_supported for r in reports)
    total_corrected = sum(r.n_corrected for r in reports)
    fully_supported = sum(1 for r in reports if r.fully_supported)

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{title}</b><br>"
                f"<sup>Genes: {total_genes} | Introns: {total_introns} | "
                f"Supported: {total_supported} ({total_supported/total_introns*100:.1f}%) | "
                f"Corrected: {total_corrected} | Fully Supported Genes: {fully_supported}</sup>"
            ),
            x=0.5,
        ),
        height=700,
        showlegend=False,
    )

    if output_path:
        output_path = Path(output_path)
        fig.write_html(str(output_path))
        logger.info(f"Saved splice summary to {output_path}")

    return fig


def _plot_splice_summary_matplotlib(
    reports: list[Any],
    output_path: Path | str | None,
    format: str,
    title: str,
) -> Any:
    """Create static Matplotlib splice summary."""
    import matplotlib.pyplot as plt

    # Extract data
    support_ratios = [r.support_ratio for r in reports]
    n_corrected = [r.n_corrected for r in reports]
    n_canonical = sum(r.n_canonical for r in reports)
    n_noncanonical = sum(r.n_noncanonical for r in reports)

    correction_distances = []
    for r in reports:
        for c in r.corrections:
            correction_distances.append(abs(c.shift))

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel 1: Support ratio
    ax1 = axes[0, 0]
    ax1.hist(support_ratios, bins=20, color="#3498db", edgecolor="black")
    ax1.set_xlabel("Support Ratio")
    ax1.set_ylabel("Gene Count")
    ax1.set_title("Support Ratio Distribution")
    ax1.set_xlim(0, 1)

    # Panel 2: Corrections per gene
    ax2 = axes[0, 1]
    max_corr = max(n_corrected) if n_corrected else 1
    ax2.hist(n_corrected, bins=range(max_corr + 2), color="#27ae60", edgecolor="black")
    ax2.set_xlabel("Corrections per Gene")
    ax2.set_ylabel("Gene Count")
    ax2.set_title("Corrections per Gene")

    # Panel 3: Splice types
    ax3 = axes[1, 0]
    ax3.pie(
        [n_canonical, n_noncanonical],
        labels=["Canonical", "Non-canonical"],
        colors=["#27ae60", "#e74c3c"],
        autopct="%1.1f%%",
        startangle=90,
    )
    ax3.set_title("Splice Site Types")

    # Panel 4: Correction distances
    ax4 = axes[1, 1]
    if correction_distances:
        ax4.hist(correction_distances, bins=15, color="#f39c12", edgecolor="black")
    ax4.set_xlabel("Correction Distance (bp)")
    ax4.set_ylabel("Count")
    ax4.set_title("Correction Distance Distribution")

    # Summary
    total_genes = len(reports)
    total_introns = sum(r.n_introns for r in reports)
    total_supported = sum(r.n_supported for r in reports)
    total_corrected = sum(r.n_corrected for r in reports)

    fig.suptitle(
        f"{title}\n"
        f"Genes: {total_genes} | Introns: {total_introns} | "
        f"Supported: {total_supported} ({total_supported/total_introns*100:.1f}%) | "
        f"Corrected: {total_corrected}",
        fontsize=12,
        fontweight="bold",
    )

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
        logger.info(f"Saved splice summary to {output_path}")

    return fig
