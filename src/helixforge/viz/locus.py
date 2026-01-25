"""Locus-level visualization.

This module provides detailed visualization of individual gene loci,
including gene models, RNA-seq evidence, and annotations.

Features:
    - Gene structure diagrams
    - Coverage tracks
    - Splice junction arcs
    - Domain annotations
    - Confidence heatmaps

Example:
    >>> from helixforge.viz.locus import LocusPlotter, plot_gene_confidence
    >>> plotter = LocusPlotter()
    >>> fig = plotter.plot(gene, evidence)
    >>> fig.savefig("locus.png")
    >>>
    >>> # Confidence visualization
    >>> from helixforge.core.confidence import GeneConfidence, RegionConfidence
    >>> fig = plot_gene_confidence(gene, confidence, region_conf)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    import matplotlib.figure
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go

    from helixforge.core.confidence import GeneConfidence, RegionConfidence
    from helixforge.core.models import GeneModel
    from helixforge.io.gff import GeneModel as GFFGeneModel
    from helixforge.isoforms.evidence import EvidenceCollection

logger = logging.getLogger(__name__)

# =============================================================================
# Plotter Class
# =============================================================================


class LocusPlotter:
    """Creates detailed locus visualizations.

    Generates multi-track visualizations showing gene structure,
    RNA-seq evidence, and annotations.

    Attributes:
        figsize: Figure size (width, height) in inches.
        dpi: Resolution for raster output.
        style: Matplotlib style to use.

    Example:
        >>> plotter = LocusPlotter(figsize=(12, 6))
        >>> fig = plotter.plot(gene, evidence)
    """

    def __init__(
        self,
        figsize: tuple[float, float] = (12, 6),
        dpi: int = 150,
        style: str = "default",
    ) -> None:
        """Initialize the locus plotter.

        Args:
            figsize: Figure size in inches.
            dpi: Resolution for raster output.
            style: Matplotlib style name.
        """
        self.figsize = figsize
        self.dpi = dpi
        self.style = style

    def plot(
        self,
        gene: "GeneModel",
        evidence: "EvidenceCollection | None" = None,
        show_domains: bool = True,
        show_coverage: bool = True,
        show_junctions: bool = True,
    ) -> "matplotlib.figure.Figure":
        """Create a locus visualization.

        Args:
            gene: Gene model to visualize.
            evidence: Optional evidence to display.
            show_domains: Show domain annotations.
            show_coverage: Show coverage track.
            show_junctions: Show junction arcs.

        Returns:
            Matplotlib figure.
        """
        # TODO: Implement plotting
        raise NotImplementedError("plot not yet implemented")

    def plot_gene_structure(
        self,
        ax: Any,  # matplotlib.axes.Axes
        gene: "GeneModel",
    ) -> None:
        """Plot gene structure on an axes.

        Args:
            ax: Matplotlib axes.
            gene: Gene model to plot.
        """
        # TODO: Implement gene structure plotting
        raise NotImplementedError("plot_gene_structure not yet implemented")

    def plot_coverage(
        self,
        ax: Any,  # matplotlib.axes.Axes
        evidence: "EvidenceCollection",
    ) -> None:
        """Plot coverage track on an axes.

        Args:
            ax: Matplotlib axes.
            evidence: Evidence with coverage data.
        """
        # TODO: Implement coverage plotting
        raise NotImplementedError("plot_coverage not yet implemented")

    def plot_junctions(
        self,
        ax: Any,  # matplotlib.axes.Axes
        evidence: "EvidenceCollection",
    ) -> None:
        """Plot splice junction arcs on an axes.

        Args:
            ax: Matplotlib axes.
            evidence: Evidence with junction data.
        """
        # TODO: Implement junction plotting
        raise NotImplementedError("plot_junctions not yet implemented")

    def save(
        self,
        fig: "matplotlib.figure.Figure",
        path: Path | str,
    ) -> None:
        """Save figure to file.

        Args:
            fig: Figure to save.
            path: Output file path.
        """
        # TODO: Implement saving
        raise NotImplementedError("save not yet implemented")


# =============================================================================
# Convenience Functions
# =============================================================================


def plot_locus(
    gene: "GeneModel",
    evidence: "EvidenceCollection | None" = None,
    output: Path | str | None = None,
    **kwargs: Any,
) -> "matplotlib.figure.Figure":
    """Quick function to plot a gene locus.

    Args:
        gene: Gene model to visualize.
        evidence: Optional evidence to display.
        output: If provided, save to this path.
        **kwargs: Additional arguments for LocusPlotter.

    Returns:
        Matplotlib figure.
    """
    plotter = LocusPlotter(**kwargs)
    fig = plotter.plot(gene, evidence)

    if output is not None:
        plotter.save(fig, output)

    return fig


def plot_comparison(
    genes: list["GeneModel"],
    labels: list[str] | None = None,
    output: Path | str | None = None,
) -> "matplotlib.figure.Figure":
    """Plot multiple gene models for comparison.

    Args:
        genes: Gene models to compare.
        labels: Labels for each gene model.
        output: If provided, save to this path.

    Returns:
        Matplotlib figure.
    """
    # TODO: Implement comparison plotting
    raise NotImplementedError("plot_comparison not yet implemented")


# =============================================================================
# Confidence Visualization Functions
# =============================================================================


# Color schemes
CONFIDENCE_COLORS = {
    "high": "#2ecc71",  # Green
    "medium": "#f39c12",  # Orange/Yellow
    "low": "#e74c3c",  # Red
}

HELIXER_CLASS_COLORS = {
    0: "#3498db",  # Intergenic - Blue
    1: "#9b59b6",  # UTR - Purple
    2: "#2ecc71",  # CDS - Green
    3: "#f39c12",  # Intron - Orange
}

HELIXER_CLASS_NAMES = {
    0: "Intergenic",
    1: "UTR",
    2: "CDS",
    3: "Intron",
}


def plot_gene_confidence(
    gene: "GFFGeneModel",
    confidence: "GeneConfidence",
    region_conf: "RegionConfidence | None" = None,
    output_path: Path | str | None = None,
    format: str = "html",
    show_legend: bool = True,
) -> Any:
    """Create detailed confidence visualization for a gene.

    Generates a multi-panel plot showing:
    - Gene structure (exons, introns, CDS)
    - Per-base probability heatmap (if region_conf provided)
    - Entropy track (if region_conf provided)
    - Smoothed confidence track
    - Low-confidence regions highlighted
    - Summary metrics panel

    Args:
        gene: Gene model to visualize.
        confidence: GeneConfidence scores for the gene.
        region_conf: Optional RegionConfidence with per-base data.
        output_path: If provided, save figure to this path.
        format: Output format ("html" for Plotly, "png"/"pdf" for Matplotlib).
        show_legend: Whether to show legend.

    Returns:
        Plotly Figure (if format="html") or Matplotlib Figure.

    Example:
        >>> calc = ConfidenceCalculator(reader, genome)
        >>> score = calc.score_gene(gene)
        >>> region = calc.get_region_confidence(gene.seqid, gene.start, gene.end)
        >>> fig = plot_gene_confidence(gene, score, region, "gene1.html")
    """
    if format == "html":
        return _plot_gene_confidence_plotly(
            gene, confidence, region_conf, output_path, show_legend
        )
    else:
        return _plot_gene_confidence_matplotlib(
            gene, confidence, region_conf, output_path, format, show_legend
        )


def _plot_gene_confidence_plotly(
    gene: "GFFGeneModel",
    confidence: "GeneConfidence",
    region_conf: "RegionConfidence | None",
    output_path: Path | str | None,
    show_legend: bool,
) -> "go.Figure":
    """Create interactive Plotly visualization."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Determine number of rows based on available data
    n_rows = 3  # Gene structure, confidence track, metrics
    if region_conf is not None:
        n_rows = 5  # Add probability heatmap and entropy track

    row_heights = (
        [0.15, 0.25, 0.20, 0.20, 0.20]
        if n_rows == 5
        else [0.25, 0.35, 0.40]
    )

    subplot_titles = (
        ["Gene Structure", "Probability Heatmap", "Entropy", "Confidence", ""]
        if n_rows == 5
        else ["Gene Structure", "Confidence", ""]
    )

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        row_heights=row_heights[:n_rows],
        shared_xaxes=True,
        vertical_spacing=0.05,
        subplot_titles=subplot_titles[:n_rows],
    )

    gene_start = confidence.start
    gene_end = confidence.end
    x_positions = np.arange(gene_start, gene_end)

    current_row = 1

    # Row 1: Gene Structure
    _add_gene_structure_trace(fig, gene, confidence, current_row)
    current_row += 1

    # Rows 2-3: Probability heatmap and entropy (if region_conf provided)
    if region_conf is not None:
        # Probability heatmap
        _add_probability_heatmap(fig, region_conf, gene_start, current_row)
        current_row += 1

        # Entropy track
        fig.add_trace(
            go.Scatter(
                x=x_positions,
                y=region_conf.per_base_entropy,
                mode="lines",
                name="Entropy",
                line=dict(color="#8e44ad", width=1),
                fill="tozeroy",
                fillcolor="rgba(142, 68, 173, 0.3)",
            ),
            row=current_row,
            col=1,
        )
        fig.update_yaxes(title_text="Entropy", row=current_row, col=1)
        current_row += 1

    # Confidence track
    if region_conf is not None:
        y_data = region_conf.smoothed_confidence
    else:
        # Create simple confidence bar
        y_data = np.full(gene_end - gene_start, confidence.mean_prob)

    fig.add_trace(
        go.Scatter(
            x=x_positions,
            y=y_data,
            mode="lines",
            name="Confidence",
            line=dict(color="#27ae60", width=2),
            fill="tozeroy",
            fillcolor="rgba(39, 174, 96, 0.3)",
        ),
        row=current_row,
        col=1,
    )

    # Add low-confidence regions as red highlights
    for lcr_start, lcr_end, lcr_score in confidence.low_confidence_regions:
        fig.add_vrect(
            x0=lcr_start,
            x1=lcr_end,
            fillcolor="rgba(231, 76, 60, 0.3)",
            layer="below",
            line_width=0,
            row=current_row,
            col=1,
        )

    # Add threshold line
    fig.add_hline(
        y=0.7,
        line_dash="dash",
        line_color="red",
        annotation_text="Threshold",
        row=current_row,
        col=1,
    )

    fig.update_yaxes(
        title_text="Confidence",
        range=[0, 1],
        row=current_row,
        col=1,
    )

    # Add metrics annotation in the last row space
    metrics_text = _format_metrics_text(confidence)
    fig.add_annotation(
        x=0.5,
        y=0.02,
        xref="paper",
        yref="paper",
        text=metrics_text,
        showarrow=False,
        font=dict(family="monospace", size=10),
        align="left",
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="#cccccc",
        borderwidth=1,
    )

    # Update layout
    conf_color = CONFIDENCE_COLORS.get(confidence.confidence_class, "#95a5a6")
    fig.update_layout(
        title=dict(
            text=(
                f"<b>{confidence.gene_id}</b> | "
                f"{confidence.seqid}:{confidence.start:,}-{confidence.end:,} "
                f"({confidence.strand}) | "
                f'<span style="color:{conf_color}">{confidence.confidence_class.upper()}</span> '
                f"(Score: {confidence.overall_score:.3f})"
            ),
            x=0.5,
        ),
        height=150 * n_rows + 100,
        showlegend=show_legend,
        xaxis_title="Genomic Position",
        hovermode="x unified",
    )

    # Set consistent x-axis range
    fig.update_xaxes(range=[gene_start, gene_end])

    if output_path:
        output_path = Path(output_path)
        fig.write_html(str(output_path))
        logger.info(f"Saved confidence plot to {output_path}")

    return fig


def _add_gene_structure_trace(
    fig: "go.Figure",
    gene: "GFFGeneModel",
    confidence: "GeneConfidence",
    row: int,
) -> None:
    """Add gene structure visualization to figure."""
    import plotly.graph_objects as go

    gene_start = confidence.start
    gene_end = confidence.end

    # Gene line (thin line across gene)
    fig.add_trace(
        go.Scatter(
            x=[gene_start, gene_end],
            y=[0.5, 0.5],
            mode="lines",
            line=dict(color="#2c3e50", width=2),
            name="Gene",
            showlegend=False,
        ),
        row=row,
        col=1,
    )

    # Add exons as blocks
    if gene.transcripts:
        transcript = gene.transcripts[0]

        for i, (exon_start, exon_end) in enumerate(transcript.exons):
            # Get exon score for color
            if i < len(confidence.exon_scores):
                score = confidence.exon_scores[i]
                color = _score_to_color(score)
            else:
                color = "#3498db"

            # Exon rectangle
            fig.add_shape(
                type="rect",
                x0=exon_start,
                x1=exon_end,
                y0=0.25,
                y1=0.75,
                fillcolor=color,
                line=dict(color="#2c3e50", width=1),
                row=row,
                col=1,
            )

            # Add exon label
            fig.add_annotation(
                x=(exon_start + exon_end) / 2,
                y=0.9,
                text=f"E{i+1}",
                showarrow=False,
                font=dict(size=8),
                row=row,
                col=1,
            )

        # Highlight CDS if present
        for cds_tuple in transcript.cds:
            cds_start, cds_end = cds_tuple[0], cds_tuple[1]
            fig.add_shape(
                type="rect",
                x0=cds_start,
                x1=cds_end,
                y0=0.3,
                y1=0.7,
                fillcolor="rgba(46, 204, 113, 0.8)",
                line=dict(color="#27ae60", width=2),
                row=row,
                col=1,
            )

    # Add strand arrow
    arrow_x = gene_start + (gene_end - gene_start) * 0.1
    arrow_text = "→" if confidence.strand == "+" else "←"
    fig.add_annotation(
        x=arrow_x,
        y=0.5,
        text=arrow_text,
        showarrow=False,
        font=dict(size=16),
        row=row,
        col=1,
    )

    fig.update_yaxes(
        range=[0, 1],
        showticklabels=False,
        row=row,
        col=1,
    )


def _add_probability_heatmap(
    fig: "go.Figure",
    region_conf: "RegionConfidence",
    gene_start: int,
    row: int,
) -> None:
    """Add probability heatmap to figure."""
    import plotly.graph_objects as go

    # Use max prob for simple heatmap
    x_positions = np.arange(gene_start, gene_start + len(region_conf.per_base_max_prob))

    fig.add_trace(
        go.Heatmap(
            x=x_positions,
            y=["Confidence"],
            z=[region_conf.per_base_max_prob],
            colorscale=[
                [0, "#e74c3c"],  # Low - Red
                [0.5, "#f39c12"],  # Medium - Yellow
                [1, "#27ae60"],  # High - Green
            ],
            zmin=0,
            zmax=1,
            showscale=True,
            colorbar=dict(
                title="Prob",
                len=0.3,
                y=0.8,
            ),
        ),
        row=row,
        col=1,
    )

    fig.update_yaxes(showticklabels=False, row=row, col=1)


def _score_to_color(score: float) -> str:
    """Convert confidence score to color."""
    if score >= 0.85:
        return "#27ae60"  # Green
    elif score >= 0.70:
        return "#f39c12"  # Yellow/Orange
    else:
        return "#e74c3c"  # Red


def _format_metrics_text(confidence: "GeneConfidence") -> str:
    """Format confidence metrics as text."""
    flags_str = ", ".join(confidence.flags) if confidence.flags else "None"

    return (
        f"<b>Confidence Metrics</b><br>"
        f"Overall Score: {confidence.overall_score:.3f} ({confidence.confidence_class})<br>"
        f"Mean Prob: {confidence.mean_prob:.3f} | "
        f"Min Prob: {confidence.min_prob:.3f} | "
        f"Median Prob: {confidence.median_prob:.3f}<br>"
        f"Entropy: {confidence.entropy:.3f} | "
        f"Boundary Sharpness: {confidence.boundary_sharpness:.3f} | "
        f"Coding Consistency: {confidence.coding_consistency:.3f}<br>"
        f"Exons: {len(confidence.exon_scores)} | "
        f"Worst Exon: #{confidence.worst_exon_idx + 1} ({confidence.worst_exon_score:.3f}) | "
        f"Low Conf Regions: {confidence.n_low_conf_regions}<br>"
        f"Flags: {flags_str}"
    )


def _plot_gene_confidence_matplotlib(
    gene: "GFFGeneModel",
    confidence: "GeneConfidence",
    region_conf: "RegionConfidence | None",
    output_path: Path | str | None,
    format: str,
    show_legend: bool,
) -> "matplotlib.figure.Figure":
    """Create static Matplotlib visualization."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # Determine number of rows
    n_rows = 3 if region_conf is None else 5
    fig, axes = plt.subplots(
        n_rows, 1,
        figsize=(14, 2 * n_rows),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 2, 2, 2, 1][:n_rows]},
    )

    if n_rows == 3:
        axes = [axes[0], None, None, axes[1], axes[2]]

    gene_start = confidence.start
    gene_end = confidence.end
    x_positions = np.arange(gene_start, gene_end)

    # Row 0: Gene structure
    ax_struct = axes[0]
    ax_struct.set_ylabel("Gene")
    ax_struct.set_ylim(0, 1)
    ax_struct.set_xlim(gene_start, gene_end)

    # Gene line
    ax_struct.plot([gene_start, gene_end], [0.5, 0.5], "k-", linewidth=2)

    # Exons
    if gene.transcripts:
        transcript = gene.transcripts[0]
        for i, (exon_start, exon_end) in enumerate(transcript.exons):
            if i < len(confidence.exon_scores):
                color = _score_to_color(confidence.exon_scores[i])
            else:
                color = "#3498db"

            rect = Rectangle(
                (exon_start, 0.25),
                exon_end - exon_start,
                0.5,
                facecolor=color,
                edgecolor="black",
            )
            ax_struct.add_patch(rect)

    ax_struct.set_yticks([])
    ax_struct.set_title(
        f"{confidence.gene_id} | {confidence.seqid}:{confidence.start:,}-{confidence.end:,} "
        f"({confidence.strand}) | {confidence.confidence_class.upper()} "
        f"(Score: {confidence.overall_score:.3f})",
        fontsize=12,
        fontweight="bold",
    )

    # Rows 1-2: Probability heatmap and entropy (if available)
    if region_conf is not None:
        ax_heat = axes[1]
        ax_entropy = axes[2]

        # Heatmap
        ax_heat.imshow(
            [region_conf.per_base_max_prob],
            aspect="auto",
            cmap="RdYlGn",
            vmin=0,
            vmax=1,
            extent=[gene_start, gene_end, 0, 1],
        )
        ax_heat.set_ylabel("Prob")
        ax_heat.set_yticks([])

        # Entropy
        ax_entropy.fill_between(
            x_positions,
            region_conf.per_base_entropy,
            alpha=0.5,
            color="#8e44ad",
        )
        ax_entropy.plot(x_positions, region_conf.per_base_entropy, color="#8e44ad", linewidth=0.5)
        ax_entropy.set_ylabel("Entropy")

    # Confidence track
    ax_conf = axes[3]
    if region_conf is not None:
        y_data = region_conf.smoothed_confidence
    else:
        y_data = np.full(gene_end - gene_start, confidence.mean_prob)

    ax_conf.fill_between(x_positions, y_data, alpha=0.5, color="#27ae60")
    ax_conf.plot(x_positions, y_data, color="#27ae60", linewidth=1)

    # Low-confidence regions
    for lcr_start, lcr_end, _ in confidence.low_confidence_regions:
        ax_conf.axvspan(lcr_start, lcr_end, alpha=0.3, color="red")

    ax_conf.axhline(y=0.7, color="red", linestyle="--", label="Threshold")
    ax_conf.set_ylabel("Confidence")
    ax_conf.set_ylim(0, 1)
    if show_legend:
        ax_conf.legend(loc="upper right")

    # Metrics panel
    ax_metrics = axes[4] if n_rows == 5 else axes[2]
    ax_metrics.axis("off")

    metrics_text = (
        f"Overall Score: {confidence.overall_score:.3f} | "
        f"Mean: {confidence.mean_prob:.3f} | "
        f"Min: {confidence.min_prob:.3f} | "
        f"Entropy: {confidence.entropy:.3f}\n"
        f"Boundary: {confidence.boundary_sharpness:.3f} | "
        f"Coding: {confidence.coding_consistency:.3f} | "
        f"Worst Exon: #{confidence.worst_exon_idx + 1} ({confidence.worst_exon_score:.3f})\n"
        f"Flags: {', '.join(confidence.flags) if confidence.flags else 'None'}"
    )
    ax_metrics.text(
        0.5, 0.5, metrics_text,
        transform=ax_metrics.transAxes,
        ha="center", va="center",
        fontfamily="monospace",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"),
    )

    # Set x-axis label on bottom plot
    axes[n_rows - 1].set_xlabel("Genomic Position")

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
        logger.info(f"Saved confidence plot to {output_path}")

    return fig


# =============================================================================
# Batch Plotting Functions
# =============================================================================


def plot_gene_confidence_batch(
    genes: list["GFFGeneModel"],
    confidences: list["GeneConfidence"],
    output_dir: Path | str,
    format: str = "html",
    calc: Any = None,  # ConfidenceCalculator for region data
) -> list[Path]:
    """Generate confidence plots for multiple genes.

    Args:
        genes: List of gene models.
        confidences: List of corresponding GeneConfidence objects.
        output_dir: Directory to save plots.
        format: Output format ("html", "png", "pdf").
        calc: Optional ConfidenceCalculator for getting RegionConfidence.

    Returns:
        List of paths to generated plots.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = []
    for gene, conf in zip(genes, confidences):
        # Get region confidence if calculator provided
        region_conf = None
        if calc is not None:
            region_conf = calc.get_region_confidence(conf.seqid, conf.start, conf.end)

        # Generate filename
        ext = "html" if format == "html" else format
        filename = f"{conf.gene_id}_confidence.{ext}"
        output_path = output_dir / filename

        plot_gene_confidence(
            gene, conf, region_conf,
            output_path=output_path,
            format=format,
        )
        paths.append(output_path)

    logger.info(f"Generated {len(paths)} confidence plots in {output_dir}")
    return paths


# =============================================================================
# Splice Refinement Visualization
# =============================================================================


def plot_splice_refinement(
    gene: "GFFGeneModel",
    original_gene: "GFFGeneModel",
    report: Any,  # GeneSpliceReport
    junctions: list[Any] | None = None,  # list[SpliceJunction]
    output_path: Path | str | None = None,
    format: str = "html",
) -> Any:
    """Visualize splice site corrections for a gene.

    Creates a multi-track visualization showing:
    1. Original gene model (gray)
    2. Refined gene model (colored by correction status)
    3. RNA-seq junction arcs with read counts
    4. Correction annotations (arrows showing shifts)

    Args:
        gene: Refined gene model.
        original_gene: Original gene model before refinement.
        report: GeneSpliceReport with correction details.
        junctions: Optional list of SpliceJunction objects.
        output_path: If provided, save figure to this path.
        format: Output format ("html" for Plotly, "png"/"pdf" for Matplotlib).

    Returns:
        Plotly Figure (if format="html") or Matplotlib Figure.

    Example:
        >>> refiner = SpliceRefiner(genome, junctions)
        >>> refined, report = refiner.refine_gene(gene)
        >>> fig = plot_splice_refinement(refined, gene, report, junctions)
    """
    if format == "html":
        return _plot_splice_refinement_plotly(
            gene, original_gene, report, junctions, output_path
        )
    else:
        return _plot_splice_refinement_matplotlib(
            gene, original_gene, report, junctions, output_path, format
        )


def _plot_splice_refinement_plotly(
    gene: "GFFGeneModel",
    original_gene: "GFFGeneModel",
    report: Any,
    junctions: list[Any] | None,
    output_path: Path | str | None,
) -> Any:
    """Create interactive Plotly splice refinement plot."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    n_rows = 4 if junctions else 3
    row_heights = [0.25, 0.25, 0.25, 0.25][:n_rows]

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        row_heights=row_heights,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=["Original Gene", "Refined Gene", "Junctions", "Corrections"][:n_rows],
    )

    gene_start = min(gene.start, original_gene.start)
    gene_end = max(gene.end, original_gene.end)

    # Row 1: Original gene (gray)
    _add_gene_model_trace(fig, original_gene, 1, color="#95a5a6", name="Original")

    # Row 2: Refined gene (colored)
    _add_gene_model_trace(fig, gene, 2, color="#27ae60", name="Refined")

    # Row 3: Junctions (if provided)
    current_row = 3
    if junctions:
        for j in junctions:
            # Draw arc for junction
            mid_x = (j.start + j.end) / 2
            arc_height = 0.5 + np.log1p(j.read_count) / 10

            # Create arc path
            x_arc = np.linspace(j.start, j.end, 50)
            y_arc = arc_height * np.sin(np.linspace(0, np.pi, 50))

            fig.add_trace(
                go.Scatter(
                    x=x_arc,
                    y=y_arc,
                    mode="lines",
                    line=dict(color="#3498db", width=1 + np.log1p(j.read_count) / 3),
                    name=f"Junction ({j.read_count} reads)",
                    showlegend=False,
                    hoverinfo="text",
                    hovertext=f"{j.start}-{j.end}: {j.read_count} reads",
                ),
                row=current_row,
                col=1,
            )

        fig.update_yaxes(range=[0, 2], showticklabels=False, row=current_row, col=1)
        current_row += 1

    # Row 4 (or 3): Corrections
    if report.corrections:
        for corr in report.corrections:
            # Draw arrow from original to corrected position
            fig.add_annotation(
                x=corr.corrected_position,
                y=0.5,
                ax=corr.original_position,
                ay=0.5,
                xref=f"x{current_row}",
                yref=f"y{current_row}",
                axref=f"x{current_row}",
                ayref=f"y{current_row}",
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=2,
                arrowcolor="#e74c3c" if corr.shift != 0 else "#27ae60",
            )

            # Add label
            fig.add_trace(
                go.Scatter(
                    x=[corr.corrected_position],
                    y=[0.7],
                    mode="text",
                    text=[f"{corr.site}: {corr.shift:+d}"],
                    textposition="top center",
                    textfont=dict(size=9),
                    showlegend=False,
                ),
                row=current_row,
                col=1,
            )

    fig.update_yaxes(range=[0, 1], showticklabels=False, row=current_row, col=1)

    # Layout
    conf_class_color = "#27ae60" if report.fully_supported else "#f39c12"
    fig.update_layout(
        title=dict(
            text=(
                f"<b>{report.gene_id}</b> Splice Refinement | "
                f"Introns: {report.n_introns} | "
                f"Supported: {report.n_supported} | "
                f"Corrected: {report.n_corrected} | "
                f'<span style="color:{conf_class_color}">{"Fully Supported" if report.fully_supported else "Partial Support"}</span>'
            ),
            x=0.5,
        ),
        height=150 * n_rows + 100,
        showlegend=False,
        xaxis_title="Genomic Position",
    )

    fig.update_xaxes(range=[gene_start - 100, gene_end + 100])

    if output_path:
        output_path = Path(output_path)
        fig.write_html(str(output_path))
        logger.info(f"Saved splice refinement plot to {output_path}")

    return fig


def _add_gene_model_trace(
    fig: Any,
    gene: "GFFGeneModel",
    row: int,
    color: str,
    name: str,
) -> None:
    """Add gene model visualization to a subplot."""
    import plotly.graph_objects as go

    if not gene.transcripts:
        return

    transcript = gene.transcripts[0]

    # Gene line
    fig.add_trace(
        go.Scatter(
            x=[gene.start, gene.end],
            y=[0.5, 0.5],
            mode="lines",
            line=dict(color=color, width=2),
            name=name,
            showlegend=False,
        ),
        row=row,
        col=1,
    )

    # Exons as rectangles
    for exon_start, exon_end in transcript.exons:
        fig.add_shape(
            type="rect",
            x0=exon_start,
            x1=exon_end,
            y0=0.25,
            y1=0.75,
            fillcolor=color,
            line=dict(color="black", width=1),
            row=row,
            col=1,
        )

    fig.update_yaxes(range=[0, 1], showticklabels=False, row=row, col=1)


def _plot_splice_refinement_matplotlib(
    gene: "GFFGeneModel",
    original_gene: "GFFGeneModel",
    report: Any,
    junctions: list[Any] | None,
    output_path: Path | str | None,
    format: str,
) -> Any:
    """Create static Matplotlib splice refinement plot."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Rectangle

    n_rows = 4 if junctions else 3
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 2 * n_rows), sharex=True)

    gene_start = min(gene.start, original_gene.start)
    gene_end = max(gene.end, original_gene.end)

    # Row 0: Original gene
    ax_orig = axes[0]
    _add_gene_model_matplotlib(ax_orig, original_gene, "#95a5a6")
    ax_orig.set_ylabel("Original")
    ax_orig.set_title(f"{report.gene_id} Splice Refinement")

    # Row 1: Refined gene
    ax_ref = axes[1]
    _add_gene_model_matplotlib(ax_ref, gene, "#27ae60")
    ax_ref.set_ylabel("Refined")

    current_row = 2

    # Row 2: Junctions
    if junctions:
        ax_junc = axes[current_row]
        for j in junctions:
            # Draw arc
            mid = (j.start + j.end) / 2
            width = j.end - j.start
            height = 0.3 + np.log1p(j.read_count) / 20

            arc_x = np.linspace(j.start, j.end, 50)
            arc_y = height * np.sin(np.linspace(0, np.pi, 50))

            ax_junc.plot(
                arc_x, arc_y,
                color="#3498db",
                linewidth=1 + np.log1p(j.read_count) / 5,
            )

        ax_junc.set_ylabel("Junctions")
        ax_junc.set_ylim(0, 1)
        ax_junc.set_yticks([])
        current_row += 1

    # Row 3: Corrections
    ax_corr = axes[current_row]
    for corr in report.corrections:
        arrow = FancyArrowPatch(
            (corr.original_position, 0.5),
            (corr.corrected_position, 0.5),
            arrowstyle="->",
            color="#e74c3c" if corr.shift != 0 else "#27ae60",
            mutation_scale=15,
        )
        ax_corr.add_patch(arrow)
        ax_corr.annotate(
            f"{corr.site}: {corr.shift:+d}",
            (corr.corrected_position, 0.7),
            fontsize=8,
            ha="center",
        )

    ax_corr.set_ylabel("Corrections")
    ax_corr.set_ylim(0, 1)
    ax_corr.set_yticks([])
    ax_corr.set_xlabel("Genomic Position")

    for ax in axes:
        ax.set_xlim(gene_start - 100, gene_end + 100)

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
        logger.info(f"Saved splice refinement plot to {output_path}")

    return fig


def _add_gene_model_matplotlib(ax: Any, gene: "GFFGeneModel", color: str) -> None:
    """Add gene model to matplotlib axes."""
    from matplotlib.patches import Rectangle

    if not gene.transcripts:
        return

    transcript = gene.transcripts[0]

    # Gene line
    ax.plot([gene.start, gene.end], [0.5, 0.5], color=color, linewidth=2)

    # Exons
    for exon_start, exon_end in transcript.exons:
        rect = Rectangle(
            (exon_start, 0.25),
            exon_end - exon_start,
            0.5,
            facecolor=color,
            edgecolor="black",
        )
        ax.add_patch(rect)

    ax.set_ylim(0, 1)
    ax.set_yticks([])


def plot_splice_site_logo(
    sequences: list[str],
    site_type: str,  # "donor" or "acceptor"
    output_path: Path | str | None = None,
) -> Any:
    """Generate sequence logo for splice sites.

    Creates a sequence logo visualization showing nucleotide frequencies
    at each position around splice sites.

    Args:
        sequences: List of aligned splice site sequences.
        site_type: "donor" or "acceptor".
        output_path: If provided, save figure to this path.

    Returns:
        Matplotlib Figure.

    Example:
        >>> donor_seqs = ["CAGGTAAGT", "AAGGTGAGT", "GAGGTAAGA"]
        >>> fig = plot_splice_site_logo(donor_seqs, "donor", "donor_logo.png")
    """
    import matplotlib.pyplot as plt

    if not sequences:
        raise ValueError("No sequences provided")

    seq_len = len(sequences[0])
    if not all(len(s) == seq_len for s in sequences):
        raise ValueError("All sequences must have same length")

    # Calculate frequencies
    freqs = np.zeros((seq_len, 4))
    for seq in sequences:
        seq = seq.upper()
        for i, base in enumerate(seq):
            if base == "A":
                freqs[i, 0] += 1
            elif base == "C":
                freqs[i, 1] += 1
            elif base == "G":
                freqs[i, 2] += 1
            elif base == "T":
                freqs[i, 3] += 1

    # Normalize
    freqs = freqs / len(sequences)

    # Calculate information content (bits)
    background = 0.25  # Equal background
    info_content = np.zeros((seq_len, 4))
    for i in range(seq_len):
        for j in range(4):
            if freqs[i, j] > 0:
                info_content[i, j] = freqs[i, j] * np.log2(freqs[i, j] / background)
            else:
                info_content[i, j] = 0

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 3))

    bases = ["A", "C", "G", "T"]
    colors = {"A": "#2ecc71", "C": "#3498db", "G": "#f39c12", "T": "#e74c3c"}

    positions = (
        list(range(-3, 7)) if site_type == "donor" else list(range(-14, 2))
    )
    if len(positions) != seq_len:
        positions = list(range(seq_len))

    for i in range(seq_len):
        # Sort bases by information content at this position
        sorted_indices = np.argsort(info_content[i])
        y_offset = 0

        for idx in sorted_indices:
            height = max(0, info_content[i, idx])
            if height > 0.01:
                ax.bar(
                    i,
                    height,
                    bottom=y_offset,
                    color=colors[bases[idx]],
                    width=0.9,
                    edgecolor="black",
                    linewidth=0.5,
                )
                if height > 0.2:
                    ax.text(
                        i,
                        y_offset + height / 2,
                        bases[idx],
                        ha="center",
                        va="center",
                        fontsize=10,
                        fontweight="bold",
                    )
                y_offset += height

    ax.set_xlim(-0.5, seq_len - 0.5)
    ax.set_ylim(0, 2)
    ax.set_xticks(range(seq_len))
    ax.set_xticklabels(positions)
    ax.set_xlabel("Position relative to splice site")
    ax.set_ylabel("Information (bits)")
    ax.set_title(f"{site_type.capitalize()} Splice Site Logo (n={len(sequences)})")

    # Add vertical line at splice site
    if site_type == "donor":
        ax.axvline(x=2.5, color="red", linestyle="--", alpha=0.5)
    else:
        ax.axvline(x=13.5, color="red", linestyle="--", alpha=0.5)

    plt.tight_layout()

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved splice site logo to {output_path}")

    return fig
