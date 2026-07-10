"""Static per-locus figures."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

_VIZ_HINT = "matplotlib is required for plotting: pip install 'helixforge[viz]'"

if TYPE_CHECKING:
    from helixforge.io.hdf5 import HDF5ConfidenceReader
    from helixforge.reconcile.models import ReconciledGene, SpliceJunction


def _mpl() -> tuple[Any, Any, Any]:
    """Lazily import matplotlib with a non-interactive backend."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.patches import Rectangle
    except ImportError as exc:  # pragma: no cover - exercised via monkeypatch
        raise ImportError(_VIZ_HINT) from exc
    return plt, Rectangle, LineCollection


# Glyph heights (in lane units).
_EXON_H = 0.30
_CDS_H = 0.55
_TIER_COLORS = {1: "#1b7837", 2: "#7fbf7b", 3: "#d9b365", 4: "#b35806"}


# ---------------------------------------------------------------------------
# Model normalisation: accept ReconciledGene, HelixerLocus, or transcripts
# ---------------------------------------------------------------------------


def _models_from(obj: Any) -> list[tuple[str, list[Any], list[Any], str, int, int]]:
    """Return ``[(label, exons, cds, strand, start, end), ...]`` for a gene-like.

    ``ReconciledGene`` → one entry per transcript; ``HelixerLocus`` (single
    model) → one entry. ``None`` → empty list.
    """
    if obj is None:
        return []
    if hasattr(obj, "transcripts"):  # ReconciledGene
        primary = getattr(obj, "primary_transcript_id", None)
        out = []
        for t in obj.transcripts:
            label = t.transcript_id + (" *" if t.transcript_id == primary else "")
            out.append(
                (label, list(t.exons), list(t.cds or []), t.strand, t.start, t.end)
            )
        return out
    # HelixerLocus / single-model object.
    return [
        (
            obj.gene_id,
            list(obj.exons),
            list(obj.cds or []),
            obj.strand,
            obj.start,
            obj.end,
        )
    ]


# ---------------------------------------------------------------------------
# Drawing primitives
# ---------------------------------------------------------------------------


def _draw_models(
    ax: Any,
    models: list[tuple[str, list[Any], list[Any], str, int, int]],
    Rectangle: Any,
    x0: int,
    x1: int,
) -> int:
    """Draw stacked isoform glyphs; return the number of lanes used."""
    n = len(models)
    for lane, (label, exons, cds, strand, start, end) in enumerate(models):
        y = n - 1 - lane  # first model on top
        # intron backbone
        ax.plot([start, end], [y, y], color="#555555", lw=0.8, zorder=1)
        # strand arrow at the coding 5' end
        ax.annotate(
            "▶" if strand == "+" else "◀",
            xy=(start if strand == "+" else end, y),
            color="#555555",
            fontsize=6,
            va="center",
            ha="left" if strand == "+" else "right",
            zorder=3,
        )
        cds_spans = [(c.start, c.end) for c in cds]
        for ex in exons:
            ax.add_patch(
                Rectangle(
                    (ex.start, y - _EXON_H / 2),
                    ex.end - ex.start,
                    _EXON_H,
                    facecolor="#cccccc",
                    edgecolor="#888888",
                    lw=0.4,
                    zorder=2,
                )
            )
        for cs, ce in cds_spans:
            ax.add_patch(
                Rectangle(
                    (cs, y - _CDS_H / 2),
                    ce - cs,
                    _CDS_H,
                    facecolor="#2166ac",
                    edgecolor="#1a4f86",
                    lw=0.4,
                    zorder=4,
                )
            )
        ax.text(x0, y + _CDS_H / 2 + 0.05, label, fontsize=6, va="bottom", ha="left")
    ax.set_ylim(-0.7, max(n, 1) - 0.3)
    ax.set_xlim(x0, x1)
    ax.set_yticks([])
    return n


def _draw_junction_arcs(
    ax: Any,
    junctions: list[SpliceJunction] | None,
    gene: Any,
    y_base: float,
) -> None:
    """Semicircular arcs for junctions overlapping the gene; lw ∝ read_count."""
    import numpy as np

    if not junctions:
        return
    relevant = [
        j
        for j in junctions
        if j.seqid == gene.seqid
        and j.strand == gene.strand
        and j.donor >= gene.start
        and j.acceptor <= gene.end
    ]
    if not relevant:
        return
    max_reads = max(j.read_count for j in relevant) or 1
    for j in relevant:
        mid = (j.donor + j.acceptor) / 2.0
        rad = (j.acceptor - j.donor) / 2.0
        theta = np.linspace(0, np.pi, 40)
        xs = mid + rad * np.cos(theta)
        ys = y_base + 0.4 * np.sin(theta)
        lw = 0.5 + 3.0 * (j.read_count / max_reads)
        ax.plot(xs, ys, color="#762a83", lw=lw, alpha=0.8, zorder=2)


def _draw_confidence_strip(
    ax: Any,
    h5_reader: HDF5ConfidenceReader,
    seqid: str,
    x0: int,
    x1: int,
    y_base: float,
    height: float = 0.4,
) -> None:
    """Helixer max(CDS,UTR) per-base confidence as a heat strip."""
    import numpy as np

    try:
        preds = h5_reader.get_per_base_predictions(seqid, x0, x1)
    except (KeyError, IndexError, ValueError):
        return
    conf = np.maximum(preds[:, 2], preds[:, 1])  # CDS, UTR channels
    ax.imshow(
        conf[np.newaxis, :],
        aspect="auto",
        cmap="magma",
        vmin=0.0,
        vmax=1.0,
        extent=(x0, x1, y_base, y_base + height),
        zorder=1,
    )
    ax.text(x0, y_base + height + 0.02, "Helixer conf.", fontsize=5, va="bottom")


def _coverage_array(
    coverage: Any,
    seqid: str,
    x0: int,
    x1: int,
) -> list[float] | None:
    """Return a per-base coverage list, or ``None``.

    ``coverage`` may be a ``CoverageCalculator``-like object (has
    ``region_coverage_array``) or an already-computed sequence.
    """
    if coverage is None:
        return None
    if hasattr(coverage, "region_coverage_array"):
        try:
            arr: list[float] = coverage.region_coverage_array(seqid, x0, x1)
            return arr
        except (KeyError, IndexError, ValueError):
            return None
    return list(coverage)


def _draw_coverage(
    ax: Any,
    arr: list[float],
    x0: int,
    y_base: float,
    height: float = 0.8,
) -> None:
    if not arr:
        return
    import numpy as np

    xs = np.arange(x0, x0 + len(arr))
    mx = max(arr) or 1
    ys = y_base + height * (np.asarray(arr, dtype=float) / mx)
    ax.fill_between(xs, y_base, ys, color="#4393c3", alpha=0.6, zorder=1)
    ax.text(x0, y_base + height + 0.02, "coverage", fontsize=5, va="bottom")


def _draw_as_labels(ax: Any, gene: ReconciledGene, y: float) -> None:
    for e in gene.as_events:
        mid = (e.start + e.end) / 2.0
        ax.annotate(
            e.kind,
            xy=(mid, y),
            fontsize=5,
            color="#b2182b",
            ha="center",
            va="top",
        )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def plot_locus(
    gene: ReconciledGene,
    genome: Any = None,
    junctions: list[SpliceJunction] | None = None,
    coverage: Any = None,
    h5_reader: HDF5ConfidenceReader | None = None,
    out_path: str | Path | None = None,
    before: Any = None,
    figsize: tuple[float, float] | None = None,
) -> Any:
    """Render a per-locus figure; return the matplotlib ``Figure``.

    Evidence arguments are all optional and degrade gracefully when absent. When
    ``before`` (a ``ReconciledGene`` or ``HelixerLocus``) is given, a two-panel
    before/after figure is produced. If ``out_path`` is set the figure is saved
    (format inferred from the suffix, SVG/PDF for the manuscript).
    """
    plt, Rectangle, _ = _mpl()

    models = _models_from(gene)
    before_models = _models_from(before)
    x0 = min([gene.start, *[m[4] for m in before_models]])
    x1 = max([gene.end, *[m[5] for m in before_models]])
    pad = max(1, int((x1 - x0) * 0.02))
    x0, x1 = x0 - pad, x1 + pad

    if before is not None:
        fig, (ax_before, ax) = plt.subplots(
            2,
            1,
            figsize=figsize or (10, 6),
            gridspec_kw={"height_ratios": [1, 3]},
        )
        _draw_models(ax_before, before_models, Rectangle, x0, x1)
        ax_before.set_title(
            f"Helixer (before), {getattr(before, 'gene_id', '')}",
            fontsize=8,
            loc="left",
        )
    else:
        fig, ax = plt.subplots(figsize=figsize or (10, 4))

    n_lanes = _draw_models(ax, models, Rectangle, x0, x1)
    ax.set_title(
        f"{gene.gene_id}  (tier {gene.tier}, {gene.origin}, "
        f"{len(gene.transcripts)} iso)",
        fontsize=8,
        loc="left",
    )
    _draw_as_labels(ax, gene, -0.55)

    # Evidence tracks stacked below the isoforms (negative y region).
    y_track = -0.7
    _draw_junction_arcs(ax, junctions, gene, y_track)
    if h5_reader is not None:
        y_track -= 0.6
        _draw_confidence_strip(ax, h5_reader, gene.seqid, x0, x1, y_track)
    arr = _coverage_array(coverage, gene.seqid, x0, x1)
    if arr:
        y_track -= 1.0
        _draw_coverage(ax, arr, x0, y_track)

    ax.set_ylim(y_track - 0.2, max(n_lanes, 1) - 0.2)
    ax.set_xlabel(f"{gene.seqid} position (0-based)", fontsize=7)
    fig.tight_layout()

    if out_path is not None:
        fig.savefig(str(out_path), bbox_inches="tight")
    return fig


def _as_complexity(gene: Any) -> int:
    return len(gene.as_events)


def plot_loci(
    genes: list[ReconciledGene],
    out_dir: str | Path,
    top_n_by_as_complexity: int | None = None,
    fmt: str = "svg",
    **kw: Any,
) -> list[Path]:
    """Batch-render loci to ``out_dir/<gene_id>.<fmt>``; return the list of paths.

    With ``top_n_by_as_complexity`` only the N most AS-complex genes (by
    ``len(gene.as_events)``) are drawn. Figures are closed after saving to bound
    memory across large gene sets.
    """
    plt, _, _ = _mpl()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = list(genes)
    if top_n_by_as_complexity is not None:
        selected = sorted(selected, key=_as_complexity, reverse=True)[
            :top_n_by_as_complexity
        ]

    paths = []
    for gene in selected:
        out_path = out_dir / f"{gene.gene_id}.{fmt}"
        fig = plot_locus(gene, out_path=out_path, **kw)
        plt.close(fig)
        paths.append(out_path)
    return paths
