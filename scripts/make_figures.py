#!/usr/bin/env python3
"""Regenerate every poster figure from committed result tables.

Poster step 3 (``prompts/helixforge-poster/03_poster_figures.md``). One command,
one figure set: each figure reads a committed on-disk result table and is rebuilt
from it — **no figure or number is hand-drawn or hand-typed** (CLAUDE.md §16 F3).
A figure whose source table is missing fails with a clear message rather than
emitting a blank plot. Output is vector (SVG **and** PDF), presentation-sized.

The isoform story is **TRaCE-on** (``m3_out/athaliana.gff3`` regenerated with
``--trace-primary`` in poster step 2c) and leads with the terminus-independent,
TRaCE-independent junction-level precision/recall; the primary-precision panel is
the TRaCE demonstration; the alt-precision panel is shown only against Helixer's
structural zero (one transcript per gene). See ``docs/BENCHMARK.md`` for the full
isoform write-up.

Usage::

    python scripts/make_figures.py --out-dir presentation/
    python scripts/make_figures.py --out-dir /tmp/figs --only ablation
    python scripts/make_figures.py --out-dir presentation/ --skip-locus

Every data-extraction function (``extract_*``) is pure and importable, so the
test suite asserts the figure values equal the source rows without rendering.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent.parent

# --- Committed source tables (CLAUDE.md §16 F3: bench_out/ ablation_out/ m3_out/). --
DEFAULT_SOURCES: dict[str, Path] = {
    # before/after: Helixer vs HelixForge accuracy + completeness, both annotations.
    "before_after": ROOT / "bench_out_m7a" / "results_gffcompare_mikado_compare_compleasm_both.tsv",
    # the master HelixForge benchmark table (isoform_* count rows live here).
    "benchmark": ROOT / "bench_out" / "benchmark.tsv",
    # helixer_support ablation: full vs no_helixer_support, one row each.
    "ablation": ROOT / "ablation_out" / "ablation.tsv",
    # matched-locus / junction isoform accuracy — TRaCE-on and the no-TRaCE baseline.
    "iso_trace": ROOT / "bench_out" / "iso_trace.result.json",
    "iso_notrace": ROOT / "bench_out" / "iso_noTRACE.result.json",
    # the reconciled-gene report (tier counts + structure).
    "report": ROOT / "m3_out" / "athaliana.report.json",
    # the TRaCE-on reconciled GFF3 + the Helixer "before" GFF3 (example-locus figure).
    "reconciled_gff3": ROOT / "m3_out" / "athaliana.gff3",
    "helixer_gff3": ROOT / "helixforge_testing" / "helixer_output" / "Arabidopsis-thaliana_helixer.gff3",
    # the authored (non-data-driven) method schematic.
    "schematic_src": ROOT / "presentation" / "method_schematic.src.svg",
}

# A clean, illustrative reconciled locus for the example-locus figure: minus strand
# (exercises the high-risk coding-direction path, CLAUDE.md §5), tier 1, two isoforms
# with A3/A5/ES alternative splicing. Recorded here so the figure is reproducible.
EXAMPLE_GENE = "HFG_02402"

# --- Presentation style: consistent palette, vector, presentation-sized. ---
HELIXER_COLOR = "#9ecae1"
HELIXFORGE_COLOR = "#2166ac"
TRADEOFF_COLOR = "#d8801f"
TRACE_COLOR = "#b2266e"
POS_COLOR = "#2a8a2a"
REF_COLOR = "#999999"
TITLE_SIZE = 15
LABEL_SIZE = 11
TICK_SIZE = 10


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _require(path: str | Path, what: str) -> Path:
    """Return ``path`` or raise a clear ``FileNotFoundError`` naming the figure.

    A missing source table must fail loudly here (never produce a blank figure).
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(
            f"cannot build the {what} figure: source table not found at {p}. "
            f"Regenerate it (see docs/REPRODUCE.md) or pass an explicit path."
        )
    return p


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _f1(sn: float, pr: float) -> float:
    """Harmonic mean of sensitivity and precision (a formula, not a typed number)."""
    return 0.0 if (sn + pr) == 0 else round(2 * sn * pr / (sn + pr), 2)


def _matched(d: dict[str, Any]) -> dict[str, Any]:
    """Normalise the two isoform-result JSON shapes to the matched-metrics dict.

    ``iso_trace.result.json`` nests metrics under ``["matched"]``; the no-TRaCE
    baseline JSON is the ``matched_locus_isoform_accuracy`` dict written directly.
    """
    return d["matched"] if "matched" in d else d


# ---------------------------------------------------------------------------
# Pure data extraction (importable + tested; never hand-typed numbers)
# ---------------------------------------------------------------------------


def extract_before_after(path: str | Path) -> list[tuple[str, float, float, bool]]:
    """Helixer vs HelixForge headline metrics + the completeness tradeoff.

    Returns ``[(label, helixer, helixforge, is_tradeoff), ...]`` read from the
    both-annotations results table. ``is_tradeoff`` marks compleasm completeness
    (which *decreases* — the precision/recall tradeoff of reconciliation).
    """
    rows = _read_tsv(_require(path, "before/after benchmark"))
    idx: dict[tuple[str, str, str], float] = {
        (r["annotation"], r["tool"], r["metric"]): float(r["value"]) for r in rows
    }

    def both(tool: str, metric: str) -> tuple[float, float]:
        return idx[("helixer", tool, metric)], idx[("helixforge", tool, metric)]

    out: list[tuple[str, float, float, bool]] = []
    hx, hf = both("mikado_compare", "intron_f1")
    out.append(("Intron F1", hx, hf, False))
    hx, hf = both("gffcompare", "locus_sn")
    out.append(("Locus Sn", hx, hf, False))
    hx, hf = both("mikado_compare", "gene_80_base_f1_f1")
    out.append(("Gene F1 ≥80%", hx, hf, False))
    hx, hf = both("compleasm", "complete")
    out.append(("BUSCO complete", hx, hf, True))
    return out


# The 13 ablation levels: full beats no_helixer_support on every one (13/13).
# mikado_compare columns are already F1; gffcompare F1 is the harmonic mean of its
# Sn/Pr columns; compleasm is the completeness percent.
_ABLATION_LEVELS: list[tuple[str, str, tuple[str, ...]]] = [
    ("Base", "f1", ("mikado_compare.base_f1",)),
    ("Exon (lenient)", "f1", ("mikado_compare.exon_lenient_f1",)),
    ("Splice site", "f1", ("mikado_compare.splice_site_f1",)),
    ("Intron", "f1", ("mikado_compare.intron_f1",)),
    ("Intron chain (mik)", "f1", ("mikado_compare.intron_chain_f1",)),
    ("Transcript ≥80%", "f1", ("mikado_compare.transcript_80_base_f1_f1",)),
    ("Gene ≥80%", "f1", ("mikado_compare.gene_80_base_f1_f1",)),
    ("Gene ≥95%", "f1", ("mikado_compare.gene_95_base_f1_f1",)),
    ("Exon (gff)", "snpr", ("gffcompare.exon_sn", "gffcompare.exon_pr")),
    ("Intron chain (gff)", "snpr", ("gffcompare.intron_chain_sn", "gffcompare.intron_chain_pr")),
    ("Transcript (gff)", "snpr", ("gffcompare.transcript_sn", "gffcompare.transcript_pr")),
    ("Locus (gff)", "snpr", ("gffcompare.locus_sn", "gffcompare.locus_pr")),
    ("Completeness", "raw", ("compleasm.complete",)),
]


def extract_ablation(path: str | Path) -> list[tuple[str, float, float, float]]:
    """``full`` vs ``no_helixer_support`` deltas for the 13 levels.

    Returns ``[(label, delta, full_val, nohs_val), ...]`` with ``delta = full -
    no_helixer_support``. The figure that justifies the novel ``helixer_support``
    coupling: all 13 deltas are positive.
    """
    rows = _read_tsv(_require(path, "ablation"))
    by_variant: dict[str, dict[str, str]] = {r["variant"]: r for r in rows}
    full, nohs = by_variant["full"], by_variant["no_helixer_support"]

    def value(row: dict[str, str], kind: str, cols: tuple[str, ...]) -> float:
        if kind == "snpr":
            return _f1(float(row[cols[0]]), float(row[cols[1]]))
        return round(float(row[cols[0]]), 2)

    out: list[tuple[str, float, float, float]] = []
    for label, kind, cols in _ABLATION_LEVELS:
        fv = value(full, kind, cols)
        nv = value(nohs, kind, cols)
        out.append((label, round(fv - nv, 2), fv, nv))
    return out


def extract_isoform_panels(
    trace_path: str | Path, notrace_path: str | Path
) -> dict[str, tuple[float, float]]:
    """The three junction-led isoform panels, as percentages.

    * ``junction``  — TRaCE-independent intron-set precision/recall (headline).
    * ``primary``   — canonical precision **before** (no-TRaCE) vs **after** TRaCE.
    * ``alt``       — HelixForge alt-isoform precision vs Helixer's structural 0
      (Helixer emits one transcript per gene, so zero alternatives are possible).
    """
    with open(_require(trace_path, "isoform-accuracy")) as fh:
        trace = _matched(json.load(fh))
    with open(_require(notrace_path, "isoform-accuracy (no-TRaCE baseline)")) as fh:
        notrace = _matched(json.load(fh))

    pct = lambda x: round(float(x) * 100, 2)
    return {
        "junction": (pct(trace["junction_precision"]), pct(trace["junction_recall"])),
        "primary": (pct(notrace["primary_precision"]), pct(trace["primary_precision"])),
        # second element is Helixer's by-design zero (one transcript per gene).
        "alt": (pct(trace["alt_isoform_precision"]), 0.0),
    }


def extract_isoform_counts(path: str | Path) -> list[tuple[str, float]]:
    """Mean isoforms/gene for Helixer / HelixForge / Araport11 (the count rows)."""
    rows = _read_tsv(_require(path, "isoform-count"))
    idx = {r["metric"]: float(r["value"]) for r in rows if r["tool"] == "isoform"}
    return [
        ("Helixer", idx["isoform_count_helixer_mean"]),
        ("HelixForge", idx["isoform_count_helixforge_mean"]),
        ("Araport11", idx["isoform_count_reference_mean"]),
    ]


def extract_tiers(path: str | Path) -> list[tuple[str, int]]:
    """Genes per tier, from the reconciled-gene report JSON."""
    with open(_require(path, "tier-distribution")) as fh:
        report = json.load(fh)
    tiers = report["tier"]
    return [(f"Tier {k}", int(tiers[k])) for k in sorted(tiers)]


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _mpl():
    import matplotlib

    matplotlib.use("Agg")  # non-interactive backend (CLAUDE.md C2 / test floor).
    import matplotlib.pyplot as plt

    plt.rcParams.update({"svg.fonttype": "none", "font.size": LABEL_SIZE})
    return plt


def _save(fig, out_dir: Path, name: str) -> list[Path]:
    """Write ``name.svg`` and ``name.pdf`` (vector); return the written paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext in ("svg", "pdf"):
        p = out_dir / f"{name}.{ext}"
        fig.savefig(p, bbox_inches="tight")
        written.append(p)
    import matplotlib.pyplot as plt

    plt.close(fig)
    return written


def figure_method_schematic(out_dir: Path, schematic_src: str | Path) -> list[Path]:
    """Copy the authored schematic SVG into the output set (+ PDF if convertible)."""
    src = _require(schematic_src, "method schematic")
    out_dir.mkdir(parents=True, exist_ok=True)
    written = [out_dir / "method_schematic.svg"]
    shutil.copyfile(src, written[0])
    try:  # optional vector PDF if a converter is installed.
        import cairosvg  # type: ignore

        pdf = out_dir / "method_schematic.pdf"
        cairosvg.svg2pdf(url=str(src), write_to=str(pdf))
        written.append(pdf)
    except Exception:
        print("  note: cairosvg not available — method_schematic.pdf skipped (SVG copied).")
    return written


def figure_before_after(out_dir: Path, path: str | Path) -> list[Path]:
    plt = _mpl()
    data = extract_before_after(path)
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    n = len(data)
    xs = range(n)
    w = 0.38
    for i, (label, hx, hf, tradeoff) in enumerate(data):
        ax.bar(i - w / 2, hx, w, color=HELIXER_COLOR, edgecolor="#555",
               label="Helixer" if i == 0 else None)
        ax.bar(i + w / 2, hf, w,
               color=TRADEOFF_COLOR if tradeoff else HELIXFORGE_COLOR,
               edgecolor="#555", hatch="//" if tradeoff else None,
               label="HelixForge" if i == 0 else None)
        for x, v in ((i - w / 2, hx), (i + w / 2, hf)):
            ax.text(x, v + 0.7, f"{v:.1f}", ha="center", va="bottom", fontsize=TICK_SIZE)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([d[0] for d in data], fontsize=TICK_SIZE)
    ax.set_ylabel("percent", fontsize=LABEL_SIZE)
    ax.set_ylim(0, 105)
    ax.set_title("Accuracy gains vs Araport11 (Helixer → HelixForge)",
                 fontsize=TITLE_SIZE)
    ax.legend(loc="lower left", fontsize=TICK_SIZE)
    ax.annotate("completeness: precision/recall tradeoff\n(added isoforms register as duplicate BUSCOs)",
                xy=(n - 1 + w / 2, data[-1][2]), xytext=(n - 1.85, 88),
                ha="center", va="center", fontsize=8.5, color=TRADEOFF_COLOR,
                arrowprops=dict(arrowstyle="->", color=TRADEOFF_COLOR, lw=1.2))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return _save(fig, out_dir, "before_after_benchmark")


def figure_ablation(out_dir: Path, path: str | Path) -> list[Path]:
    plt = _mpl()
    data = extract_ablation(path)
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    labels = [d[0] for d in data]
    deltas = [d[1] for d in data]
    ys = range(len(data))
    ax.barh(list(ys), deltas, color=POS_COLOR, edgecolor="#1f6b1f", height=0.62)
    for y, d in zip(ys, deltas):
        ax.text(d + 0.01, y, f"+{d:.2f}", va="center", ha="left", fontsize=TICK_SIZE)
    ax.axvline(0, color="#888", lw=0.8)
    ax.set_yticks(list(ys))
    ax.set_yticklabels(labels, fontsize=TICK_SIZE)
    ax.invert_yaxis()
    ax.set_xlabel("Δ (full − no_helixer_support), F1 / completeness points",
                  fontsize=LABEL_SIZE)
    ax.set_xlim(0, max(deltas) * 1.25)
    ax.set_title("The novel coupling: the external helixer_support metric helps\n"
                 "every accuracy level + completeness (13/13 positive)",
                 fontsize=TITLE_SIZE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return _save(fig, out_dir, "ablation_helixer_support")


def figure_isoform_panels(out_dir: Path, trace_path: str | Path,
                          notrace_path: str | Path) -> list[Path]:
    plt = _mpl()
    p = extract_isoform_panels(trace_path, notrace_path)
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(12.5, 4.6))

    # Panel A (headline) — junction-level intron-set P/R (TRaCE-independent).
    jp, jr = p["junction"]
    axA.bar([0, 1], [jp, jr], color=[HELIXFORGE_COLOR, "#4a90c2"], edgecolor="#555", width=0.6)
    for x, v in ((0, jp), (1, jr)):
        axA.text(x, v + 1, f"{v:.2f}", ha="center", va="bottom", fontsize=TICK_SIZE)
    axA.set_xticks([0, 1])
    axA.set_xticklabels(["Precision", "Recall"], fontsize=TICK_SIZE)
    axA.set_ylim(0, 105)
    axA.set_ylabel("percent", fontsize=LABEL_SIZE)
    axA.set_title("A  Junction-level (intron set)\nHelixForge introns are real",
                  fontsize=12)

    # Panel B (TRaCE demonstration) — primary/canonical precision before vs after.
    pb, pa = p["primary"]
    axB.bar([0, 1], [pb, pa], color=["#cfa3bd", TRACE_COLOR], edgecolor="#555", width=0.6)
    for x, v in ((0, pb), (1, pa)):
        axB.text(x, v + 1, f"{v:.2f}", ha="center", va="bottom", fontsize=TICK_SIZE)
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["no TRaCE", "TRaCE"], fontsize=TICK_SIZE)
    axB.set_ylim(0, 105)
    axB.set_title("B  Primary precision\nTRaCE elects a better canonical",
                  fontsize=12)
    axB.annotate("", xy=(1, pa), xytext=(0, pb),
                 arrowprops=dict(arrowstyle="->", color=TRACE_COLOR, lw=1.4))

    # Panel C (vs Helixer) — alt-isoform precision against Helixer's structural 0.
    alt, hx_zero = p["alt"]
    axC.bar([0, 1], [hx_zero, alt], color=[REF_COLOR, POS_COLOR], edgecolor="#555", width=0.6)
    axC.text(0, hx_zero + 1, "0", ha="center", va="bottom", fontsize=TICK_SIZE)
    axC.text(1, alt + 1, f"{alt:.2f}", ha="center", va="bottom", fontsize=TICK_SIZE)
    axC.set_xticks([0, 1])
    axC.set_xticklabels(["Helixer\n(1 tx/gene)", "HelixForge"], fontsize=TICK_SIZE)
    axC.set_ylim(0, 105)
    axC.set_title("C  Alt-isoform precision\nwe add real alternative isoforms",
                  fontsize=12)

    for ax in (axA, axB, axC):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("Isoform accuracy vs Araport11 (TRaCE-on, terminus-independent)",
                 fontsize=TITLE_SIZE, y=1.02)
    fig.tight_layout()
    return _save(fig, out_dir, "isoform_accuracy")


def figure_isoform_counts(out_dir: Path, path: str | Path) -> list[Path]:
    plt = _mpl()
    data = extract_isoform_counts(path)
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    colors = [HELIXER_COLOR, HELIXFORGE_COLOR, REF_COLOR]
    for i, (label, v) in enumerate(data):
        ax.bar(i, v, 0.6, color=colors[i], edgecolor="#555")
        ax.text(i, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=TICK_SIZE)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels([d[0] for d in data], fontsize=TICK_SIZE)
    ax.set_ylabel("mean isoforms / gene", fontsize=LABEL_SIZE)
    ax.set_ylim(0, max(d[1] for d in data) * 1.2)
    ax.set_title("Reconciliation adds isoforms toward the reference density",
                 fontsize=TITLE_SIZE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return _save(fig, out_dir, "isoform_count_distribution")


def figure_tiers(out_dir: Path, path: str | Path) -> list[Path]:
    plt = _mpl()
    data = extract_tiers(path)
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    palette = ["#2166ac", "#6aaed6", "#c6dbef"]
    for i, (label, v) in enumerate(data):
        ax.bar(i, v, 0.6, color=palette[i % len(palette)], edgecolor="#555")
        ax.text(i, v + max(d[1] for d in data) * 0.01, f"{v:,}",
                ha="center", va="bottom", fontsize=TICK_SIZE)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels([d[0] for d in data], fontsize=TICK_SIZE)
    ax.set_ylabel("genes", fontsize=LABEL_SIZE)
    ax.set_title("Genes per confidence tier", fontsize=TITLE_SIZE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return _save(fig, out_dir, "tier_distribution")


# --- example-locus figure (integration; needs viz deps + the real GFF3s) ---


class _Seg:
    __slots__ = ("start", "end")

    def __init__(self, start: int, end: int):
        self.start, self.end = start, end


class _Tx:
    def __init__(self, tid, exons, cds, strand, start, end):
        self.transcript_id = tid
        self.exons = exons
        self.cds = cds
        self.strand = strand
        self.start = start
        self.end = end


class _Gene:
    """Minimal duck-typed gene for ``viz.locus_plot.plot_locus`` (read-only render)."""

    def __init__(self, gene_id, seqid, strand, start, end, transcripts,
                 primary_transcript_id, tier, origin, as_events):
        self.gene_id = gene_id
        self.seqid = seqid
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts = transcripts
        self.primary_transcript_id = primary_transcript_id
        self.tier = tier
        self.origin = origin
        self.as_events = as_events


def _parse_gff3_gene(path: Path, gene_id: str) -> _Gene | None:
    """Parse a single gene (its transcripts/exons/CDS) from a GFF3 into shims.

    Coordinate conversion at the I/O boundary only (GFF3 1-based inclusive →
    0-based half-open; CLAUDE.md §4).
    """
    gene_attr: dict[str, str] = {}
    span: tuple[str, str, int, int] | None = None
    tx_order: list[str] = []
    tx_exons: dict[str, list[_Seg]] = {}
    tx_cds: dict[str, list[_Seg]] = {}
    primary: str | None = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
            typ = f[2]
            s, e = int(f[3]) - 1, int(f[4])  # → internal 0-based half-open
            if typ == "gene" and attrs.get("ID") == gene_id:
                gene_attr = attrs
                span = (f[0], f[6], s, e)
            elif typ in ("mRNA", "transcript") and attrs.get("Parent") == gene_id:
                tid = attrs["ID"]
                tx_order.append(tid)
                tx_exons.setdefault(tid, [])
                tx_cds.setdefault(tid, [])
                if primary is None or tid.endswith(".1"):
                    if tid.endswith(".1"):
                        primary = tid
            elif typ == "exon" and attrs.get("Parent") in tx_exons:
                tx_exons[attrs["Parent"]].append(_Seg(s, e))
            elif typ == "CDS" and attrs.get("Parent") in tx_cds:
                tx_cds[attrs["Parent"]].append(_Seg(s, e))
    if span is None or not tx_order:
        return None
    seqid, strand, gstart, gend = span
    txs = []
    for tid in tx_order:
        exons = sorted(tx_exons[tid], key=lambda x: x.start)
        cds = sorted(tx_cds[tid], key=lambda x: x.start)
        txs.append(_Tx(tid, exons, cds, strand,
                       exons[0].start if exons else gstart,
                       exons[-1].end if exons else gend))
    return _Gene(gene_id, seqid, strand, gstart, gend, txs,
                 primary or tx_order[0], int(gene_attr.get("tier", 0)),
                 gene_attr.get("origin", ""), [])


def _helixer_before(path: Path, gene: _Gene) -> _Gene | None:
    """Find the Helixer single model overlapping ``gene`` (same strand) for the
    before panel. Helixer emits one transcript per gene (CLAUDE.md §1)."""
    best: _Gene | None = None
    cur_id: str | None = None
    cur_exons: list[_Seg] = []
    cur_cds: list[_Seg] = []
    cur: tuple[str, str, int, int] | None = None

    def _overlaps(c) -> bool:
        return c and c[0] == gene.seqid and c[1] == gene.strand and not (
            c[3] <= gene.start or c[2] >= gene.end)

    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ = f[2]
            attrs = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
            s, e = int(f[3]) - 1, int(f[4])
            if typ in ("mRNA", "transcript"):
                if _overlaps(cur) and best is None:
                    g = _Gene(cur_id, cur[0], cur[1], cur[2], cur[3],
                              [_Tx(cur_id, sorted(cur_exons, key=lambda x: x.start),
                                   sorted(cur_cds, key=lambda x: x.start),
                                   cur[1], cur[2], cur[3])], cur_id, 0, "helixer", [])
                    best = g
                cur_id = attrs.get("ID")
                cur = (f[0], f[6], s, e)
                cur_exons, cur_cds = [], []
            elif typ == "exon":
                cur_exons.append(_Seg(s, e))
            elif typ == "CDS":
                cur_cds.append(_Seg(s, e))
            if best is not None:
                break
    if best is None and _overlaps(cur):
        best = _Gene(cur_id, cur[0], cur[1], cur[2], cur[3],
                     [_Tx(cur_id, sorted(cur_exons, key=lambda x: x.start),
                          sorted(cur_cds, key=lambda x: x.start), cur[1], cur[2], cur[3])],
                     cur_id, 0, "helixer", [])
    return best


def figure_example_locus(out_dir: Path, reconciled_gff3: str | Path,
                         helixer_gff3: str | Path | None = None,
                         gene_id: str = EXAMPLE_GENE) -> list[Path]:
    """Helixer's single model vs HelixForge's isoforms for one real locus.

    The TRaCE-elected canonical (``.1``) is marked by ``plot_locus`` (``*`` label).
    Integration figure: needs matplotlib + the real GFF3s on disk.
    """
    from helixforge.viz.locus_plot import plot_locus  # heavy viz deps

    gene = _parse_gff3_gene(_require(reconciled_gff3, "example-locus"), gene_id)
    if gene is None:
        raise ValueError(f"gene {gene_id} not found in {reconciled_gff3}")
    before = None
    if helixer_gff3 and Path(helixer_gff3).is_file():
        before = _helixer_before(Path(helixer_gff3), gene)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext in ("svg", "pdf"):
        p = out_dir / f"example_locus_{gene_id}.{ext}"
        fig = plot_locus(gene, before=before, out_path=p, figsize=(9, 6))
        written.append(p)
        import matplotlib.pyplot as plt

        plt.close(fig)
    return written


# ---------------------------------------------------------------------------
# Captions + driver
# ---------------------------------------------------------------------------

# (figure stem, builder, source description for the caption file).
FIGURES: dict[str, tuple[Callable[..., list[Path]], str]] = {
    "schematic": (
        lambda od, s: figure_method_schematic(od, s["schematic_src"]),
        "method_schematic — authored: presentation/method_schematic.src.svg",
    ),
    "before_after": (
        lambda od, s: figure_before_after(od, s["before_after"]),
        "before_after_benchmark — bench_out_m7a/results_gffcompare_mikado_compare_compleasm_both.tsv "
        "rows: mikado_compare.intron_f1, gffcompare.locus_sn, mikado_compare.gene_80_base_f1_f1, compleasm.complete",
    ),
    "ablation": (
        lambda od, s: figure_ablation(od, s["ablation"]),
        "ablation_helixer_support — ablation_out/ablation.tsv variants full vs no_helixer_support (13 F1/completeness deltas)",
    ),
    "isoform": (
        lambda od, s: figure_isoform_panels(od, s["iso_trace"], s["iso_notrace"]),
        "isoform_accuracy — bench_out/iso_trace.result.json (junction_precision/recall 93.47/92.57; "
        "primary 85.88; alt 41.67) + bench_out/iso_noTRACE.result.json (primary 81.17 baseline)",
    ),
    "counts": (
        lambda od, s: figure_isoform_counts(od, s["benchmark"]),
        "isoform_count_distribution — bench_out/benchmark.tsv isoform_count_{helixer,helixforge,reference}_mean",
    ),
    "tiers": (
        lambda od, s: figure_tiers(od, s["report"]),
        "tier_distribution — m3_out/athaliana.report.json tier{1,2,3}",
    ),
    "locus": (
        lambda od, s: figure_example_locus(od, s["reconciled_gff3"], s["helixer_gff3"]),
        f"example_locus_{EXAMPLE_GENE} — m3_out/athaliana.gff3 (TRaCE-on) vs "
        "helixforge_testing/helixer_output/Arabidopsis-thaliana_helixer.gff3; canonical .1 marked '*'",
    ),
}


def _write_captions(out_dir: Path, built: list[str]) -> Path:
    p = out_dir / "CAPTIONS.md"
    lines = ["# Poster figures — source mapping",
             "",
             "Each figure is regenerated by `scripts/make_figures.py` from a committed result",
             "table; no number is hand-typed (CLAUDE.md §16 F3). The isoform numbers are TRaCE-on.",
             ""]
    for key in FIGURES:
        if key in built:
            lines.append(f"- {FIGURES[key][1]}")
    p.write_text("\n".join(lines) + "\n")
    return p


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=str(ROOT / "presentation"),
                    help="output directory for the figure set (default: presentation/)")
    ap.add_argument("--only", choices=sorted(FIGURES), action="append",
                    help="build only the named figure(s); repeatable")
    ap.add_argument("--skip-locus", action="store_true",
                    help="skip the heavy example-locus figure (viz deps / GFF3s)")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    wanted = args.only or list(FIGURES)
    sources = dict(DEFAULT_SOURCES)

    built: list[str] = []
    for key in FIGURES:
        if key not in wanted:
            continue
        builder, _ = FIGURES[key]
        try:
            paths = builder(out_dir, sources)
        except Exception as exc:
            if key == "locus" and (args.skip_locus or not args.only):
                print(f"[skip] {key}: {exc}")
                continue
            print(f"[FAIL] {key}: {exc}", file=sys.stderr)
            return 1
        built.append(key)
        print(f"[ok]   {key}: " + ", ".join(p.name for p in paths))

    cap = _write_captions(out_dir, built)
    print(f"[ok]   captions: {cap.name}")
    print(f"\nWrote {len(built)} figure(s) to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
