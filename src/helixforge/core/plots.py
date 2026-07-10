"""Confidence-score plotting for the standalone ``helixforge confidence`` command."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import numpy as np

if TYPE_CHECKING:
    import matplotlib.figure
    import plotly.graph_objects as go

    from helixforge.core.confidence import GeneConfidence, RegionConfidence
    from helixforge.core.gff import GeneModel as GFFGeneModel

logger = logging.getLogger(__name__)

# Color scheme (from v1 viz/genome.py).
CONFIDENCE_COLORS = {
    "high": "#2ecc71",  # Green
    "medium": "#f39c12",  # Orange/Yellow
    "low": "#e74c3c",  # Red
}


# =============================================================================
# Genome-wide confidence distribution  (v1 viz/genome.py)
# =============================================================================


def plot_confidence_distribution(
    scores: Iterable["GeneConfidence"],
    output_path: Path | str | None = None,
    format: str = "html",
    title: str = "Confidence Score Distribution",
) -> Any:
    """Plot confidence score distribution as a histogram and an empirical CDF.

    Creates a two-panel figure:
    - Left: histogram of overall confidence scores, with mean/median marked and
      the 5th–95th percentile band shaded; class thresholds are subtle markers.
    - Right: empirical CDF (cumulative fraction vs score) with median and
      quartiles marked, answering "what fraction of genes score below X".

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
    """Create interactive Plotly distribution plot (histogram + empirical CDF).

    Left: histogram with the mean/median marked and the 5th–95th percentile band
    shaded, so the eye lands on where the data actually is. The high/medium class
    thresholds stay but are rendered as subtle secondary markers. Right: the
    empirical CDF (cumulative fraction vs score) with median and quartiles
    marked, so a user reads "what fraction of genes score below X" directly.
    """
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    overall = np.asarray([s.overall_score for s in scores], dtype=float)

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=["Score Distribution", "Cumulative Distribution (CDF)"],
        specs=[[{"type": "xy"}, {"type": "xy"}]],
    )

    # Distribution statistics.
    mean_score = float(np.mean(overall))
    median_score = float(np.median(overall))
    std_score = float(np.std(overall))
    p5, p25, p75, p95 = (float(np.percentile(overall, q)) for q in (5, 25, 75, 95))

    # ---- Left panel: histogram --------------------------------------------
    fig.add_trace(
        go.Histogram(
            x=overall,
            nbinsx=30,
            name="Overall Score",
            marker_color="#3498db",
            opacity=0.7,
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Shade the 5th–95th percentile band — where the bulk of the data lives.
    fig.add_vrect(
        x0=p5,
        x1=p95,
        fillcolor="#3498db",
        opacity=0.10,
        line_width=0,
        annotation_text="5–95%",
        annotation_position="top left",
        row=1,
        col=1,
    )
    # Mean (solid) and median (dashed): the eye-catching markers.
    fig.add_vline(
        x=mean_score,
        line_color="#c0392b",
        line_width=2,
        annotation_text=f"mean {mean_score:.2f}",
        annotation_position="top",
        row=1,
        col=1,
    )
    fig.add_vline(
        x=median_score,
        line_dash="dash",
        line_color="#16a085",
        line_width=2,
        annotation_text=f"median {median_score:.2f}",
        annotation_position="bottom",
        row=1,
        col=1,
    )
    # Class thresholds: subtle secondary markers (faint dotted grey, no label).
    for thr in (0.70, 0.85):
        fig.add_vline(
            x=thr,
            line_dash="dot",
            line_color="grey",
            line_width=1,
            opacity=0.5,
            row=1,
            col=1,
        )

    # ---- Right panel: empirical CDF ---------------------------------------
    order = np.sort(overall)
    cdf = np.arange(1, len(order) + 1, dtype=float) / len(order)
    fig.add_trace(
        go.Scatter(
            x=order,
            y=cdf,
            mode="lines",
            line=dict(color="#3498db", width=2, shape="hv"),
            name="CDF",
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    # Median + quartiles as subtle horizontal/vertical guides.
    for frac, score, label in (
        (0.25, p25, "Q1"),
        (0.50, median_score, "median"),
        (0.75, p75, "Q3"),
    ):
        fig.add_vline(
            x=score,
            line_dash="dot",
            line_color="grey",
            line_width=1,
            opacity=0.6,
            annotation_text=f"{label} {score:.2f}",
            annotation_position="top",
            row=1,
            col=2,
        )
        fig.add_hline(
            y=frac,
            line_dash="dot",
            line_color="grey",
            line_width=1,
            opacity=0.6,
            row=1,
            col=2,
        )

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{title}</b><br>"
                f"<sup>N={len(scores)} | Mean={mean_score:.3f} | "
                f"Median={median_score:.3f} | Std={std_score:.3f}</sup>"
            ),
            x=0.5,
        ),
        height=520,
        margin=dict(b=90),
        showlegend=False,
        annotations=[
            *fig.layout.annotations,
            dict(
                text=(
                    "Scores are genome-relative — a lower mean (e.g. ~0.64 in maize) "
                    "is a baseline for that genome, not a defect. "
                    "Pick a cutoff from this distribution, not the fixed 0.85/0.70 lines."
                ),
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0.5,
                y=-0.18,
                xanchor="center",
                font=dict(size=10, color="grey"),
            ),
        ],
    )

    fig.update_xaxes(title_text="Overall Score", range=[0, 1], row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_xaxes(title_text="Overall Score", range=[0, 1], row=1, col=2)
    fig.update_yaxes(title_text="Cumulative fraction", range=[0, 1], row=1, col=2)

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
    """Create static Matplotlib distribution plot (histogram + empirical CDF).

    Left: histogram with mean/median marked and the 5th–95th percentile band
    shaded; the high/medium class thresholds stay as subtle secondary markers.
    Right: the empirical CDF with median and quartiles marked, so a user can read
    "what fraction of genes score below X" directly.
    """
    import matplotlib.pyplot as plt

    overall = np.asarray([s.overall_score for s in scores], dtype=float)

    mean_score = float(np.mean(overall))
    median_score = float(np.median(overall))
    std_score = float(np.std(overall))
    p5, p25, p75, p95 = (float(np.percentile(overall, q)) for q in (5, 25, 75, 95))

    fig, (ax_hist, ax_cdf) = plt.subplots(1, 2, figsize=(12, 5))

    # ---- Left panel: histogram --------------------------------------------
    ax_hist.hist(overall, bins=30, color="#3498db", alpha=0.7, edgecolor="black")
    # Shade where the bulk of the data lives (5th–95th percentile band).
    ax_hist.axvspan(p5, p95, color="#3498db", alpha=0.12, label="5–95% band")
    # Mean (solid) and median (dashed) — the markers the eye should find.
    ax_hist.axvline(
        mean_score, color="#c0392b", linewidth=2, label=f"mean {mean_score:.2f}"
    )
    ax_hist.axvline(
        median_score,
        color="#16a085",
        linewidth=2,
        linestyle="--",
        label=f"median {median_score:.2f}",
    )
    # Class thresholds: subtle secondary markers (faint dotted grey).
    ax_hist.axvline(0.70, color="grey", linewidth=1, linestyle=":", alpha=0.5)
    ax_hist.axvline(0.85, color="grey", linewidth=1, linestyle=":", alpha=0.5)
    ax_hist.set_xlabel("Overall Score")
    ax_hist.set_ylabel("Count")
    ax_hist.set_xlim(0, 1)
    ax_hist.set_title("Score Distribution")
    ax_hist.legend(fontsize=8)

    # ---- Right panel: empirical CDF ---------------------------------------
    order = np.sort(overall)
    cdf = np.arange(1, len(order) + 1, dtype=float) / len(order)
    ax_cdf.step(order, cdf, where="post", color="#3498db", linewidth=2)
    # Median + quartiles as subtle guides.
    for frac, score in ((0.25, p25), (0.50, median_score), (0.75, p75)):
        ax_cdf.axvline(score, color="grey", linewidth=1, linestyle=":", alpha=0.6)
        ax_cdf.axhline(frac, color="grey", linewidth=1, linestyle=":", alpha=0.6)
    ax_cdf.annotate(
        f"median {median_score:.2f}", xy=(median_score, 0.5), fontsize=8, color="grey"
    )
    ax_cdf.set_xlabel("Overall Score")
    ax_cdf.set_ylabel("Cumulative fraction")
    ax_cdf.set_xlim(0, 1)
    ax_cdf.set_ylim(0, 1)
    ax_cdf.set_title("Cumulative Distribution (CDF)")

    fig.suptitle(
        f"{title}\nN={len(scores)} | Mean={mean_score:.3f} | "
        f"Median={median_score:.3f} | Std={std_score:.3f}",
        fontsize=12,
        fontweight="bold",
    )
    # Caption: scores are genome-relative.
    fig.text(
        0.5,
        0.005,
        "Scores are genome-relative — a lower mean (e.g. ~0.64 in maize) is a "
        "baseline for that genome, not a defect. Pick a cutoff from this "
        "distribution, not the fixed 0.85/0.70 lines.",
        ha="center",
        fontsize=8,
        color="grey",
        wrap=True,
    )

    plt.tight_layout(rect=(0, 0.04, 1, 1))

    if output_path:
        output_path = Path(output_path)
        fig.savefig(output_path, format=format, dpi=150, bbox_inches="tight")
        logger.info(f"Saved confidence distribution to {output_path}")

    return fig


# =============================================================================
# Per-gene confidence plots  (v1 viz/locus.py)
# =============================================================================


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

    row_heights = [0.15, 0.25, 0.20, 0.20, 0.20] if n_rows == 5 else [0.25, 0.35, 0.40]

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
                text=f"E{i + 1}",
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
        n_rows,
        1,
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
        ax_entropy.plot(
            x_positions, region_conf.per_base_entropy, color="#8e44ad", linewidth=0.5
        )
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
        0.5,
        0.5,
        metrics_text,
        transform=ax_metrics.transAxes,
        ha="center",
        va="center",
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
            gene,
            conf,
            region_conf,
            output_path=output_path,
            format=format,
        )
        paths.append(output_path)

    logger.info(f"Generated {len(paths)} confidence plots in {output_dir}")
    return paths


# =============================================================================
# Splice Refinement Visualization
# =============================================================================
