"""Genome-wide QC report from a GFF3 annotation file."""

from __future__ import annotations

import base64
import io
import json
import math
from collections import Counter
from pathlib import Path
from statistics import mean, median
from typing import Any

from helixforge.io.gff import build_gffutils_db
from helixforge.utils.regions import gff3_to_internal


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_TRANSCRIPT_TYPES = ("mRNA", "transcript")


def _attr_first(feature: Any, key: str, default: str | None = None) -> str | None:
    """First value of a gffutils multi-valued attribute, or *default*."""
    vals = feature.attributes.get(key, [])
    if vals:
        return str(vals[0])
    return default


def _attr_list(feature: Any, key: str) -> list[str]:
    """All values of a gffutils multi-valued attribute."""
    return [str(v) for v in feature.attributes.get(key, [])]


def _n50(lengths: list[int] | list[float]) -> int | float:
    """N50 of a list of lengths."""
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2
    cumul: int | float = 0
    for length in s:
        cumul += length
        if cumul >= half:
            return length
    return s[-1]


# ---------------------------------------------------------------------------
# Core: collect statistics from a GFF3
# ---------------------------------------------------------------------------


def _collect_stats(
    gff3_path: str | Path,
    genome_path: str | Path | None = None,
    helixer_h5_path: str | Path | None = None,
) -> dict[str, Any]:
    """Parse a GFF3 and return the statistics dict."""
    db = build_gffutils_db(str(gff3_path))

    # Per-gene accumulators
    tier_counts: Counter[str] = Counter()
    origin_counts: Counter[str] = Counter()
    biotype_counts: Counter[str] = Counter()
    flag_counts: Counter[str] = Counter()
    as_event_counts: Counter[str] = Counter()

    total_genes = 0
    total_transcripts = 0
    multi_isoform_genes = 0
    coding_genes = 0

    cds_lengths: list[int] = []
    exon_counts: list[int] = []
    intron_lengths: list[int] = []

    for gene in db.features_of_type("gene", order_by=("seqid", "start")):
        total_genes += 1

        tier = _attr_first(gene, "tier", "unknown") or "unknown"
        origin = _attr_first(gene, "origin", "unknown") or "unknown"
        biotype = (
            _attr_first(gene, "gene_biotype")
            or _attr_first(gene, "biotype")
            or "unknown"
        )
        tier_counts[tier] += 1
        origin_counts[origin] += 1
        biotype_counts[biotype] += 1

        # Collect transcripts under this gene
        txs: list[Any] = []
        for ttype in _TRANSCRIPT_TYPES:
            txs.extend(db.children(gene, featuretype=ttype, order_by="start"))
        total_transcripts += len(txs)
        if len(txs) > 1:
            multi_isoform_genes += 1

        gene_has_cds = False
        for tx in txs:
            # Flags (on transcript features in HelixForge output)
            for flag in _attr_list(tx, "flags"):
                flag_counts[flag] += 1

            # AS events
            for ev in _attr_list(tx, "as_events"):
                kind = ev.split("@")[0] if "@" in ev else ev
                as_event_counts[kind] += 1

            # Exons
            exons_list = sorted(
                db.children(tx, featuretype="exon", order_by="start"),
                key=lambda e: e.start,
            )
            n_exons = len(exons_list)
            exon_counts.append(n_exons)

            # Intron lengths (gaps between consecutive exons, 0-based half-open)
            for i in range(len(exons_list) - 1):
                _, prev_end = gff3_to_internal(exons_list[i].start, exons_list[i].end)
                next_start, _ = gff3_to_internal(
                    exons_list[i + 1].start, exons_list[i + 1].end
                )
                intron_len = next_start - prev_end
                if intron_len > 0:
                    intron_lengths.append(intron_len)

            # CDS total length per transcript
            cds_feats = sorted(
                db.children(tx, featuretype="CDS", order_by="start"),
                key=lambda c: c.start,
            )
            if cds_feats:
                gene_has_cds = True
                total_cds_len = sum(c.end - c.start + 1 for c in cds_feats)
                cds_lengths.append(total_cds_len)

        if gene_has_cds:
            coding_genes += 1

    # Confidence scores from HDF5 (optional)
    confidence_scores: list[float] | None = None
    if helixer_h5_path is not None:
        confidence_scores = _read_confidence_scores(str(helixer_h5_path), db)

    # Contig stats from genome FASTA (optional)
    genome_stats: dict[str, Any] | None = None
    if genome_path is not None:
        genome_stats = _genome_stats(str(genome_path))

    stats: dict[str, Any] = {
        "total_genes": total_genes,
        "total_transcripts": total_transcripts,
        "multi_isoform_genes": multi_isoform_genes,
        "coding_genes": coding_genes,
        "tier": dict(sorted(tier_counts.items())),
        "origin": dict(sorted(origin_counts.items())),
        "biotype": dict(sorted(biotype_counts.items())),
        "flags": dict(sorted(flag_counts.items(), key=lambda kv: -kv[1])),
        "as_events": dict(sorted(as_event_counts.items(), key=lambda kv: -kv[1])),
        "cds_length": _distribution_stats(cds_lengths),
        "exon_count": _distribution_stats(exon_counts),
        "intron_length": _distribution_stats(intron_lengths),
    }
    if confidence_scores:
        stats["confidence"] = _distribution_stats(
            [int(round(s * 1000)) / 1000 for s in confidence_scores]
        )
    if genome_stats:
        stats["genome"] = genome_stats

    return stats


def _distribution_stats(values: list[int] | list[float]) -> dict[str, Any]:
    """Min, median, mean, max, N50 for a list of numeric values."""
    if not values:
        return {"count": 0, "min": 0, "median": 0, "mean": 0.0, "max": 0, "n50": 0}
    sorted_vals = sorted(values)
    return {
        "count": len(sorted_vals),
        "min": sorted_vals[0],
        "median": median(sorted_vals),
        "mean": round(mean(sorted_vals), 2),
        "max": sorted_vals[-1],
        "n50": _n50(sorted_vals),
    }


def _read_confidence_scores(h5_path: str, db: Any) -> list[float]:
    """Per-gene mean confidence from the Helixer HDF5."""
    try:
        from helixforge.io.hdf5 import HDF5ConfidenceReader
    except ImportError:
        return []
    try:
        reader = HDF5ConfidenceReader(h5_path)
    except (FileNotFoundError, ValueError):
        return []
    scores: list[float] = []
    try:
        for gene in db.features_of_type("gene", order_by=("seqid", "start")):
            start, end = gff3_to_internal(gene.start, gene.end)
            try:
                conf = reader.get_region_confidence(gene.seqid, start, end)
                if not math.isnan(conf):
                    scores.append(conf)
            except (KeyError, ValueError):
                pass
    finally:
        reader.close()
    return scores


def _genome_stats(fasta_path: str) -> dict[str, Any]:
    """Basic genome stats from a FASTA index."""
    try:
        from pyfaidx import Fasta
    except ImportError:
        return {}
    try:
        fa = Fasta(fasta_path)
    except Exception:
        return {}
    lengths = [len(fa[k]) for k in fa.keys()]
    return {
        "contigs": len(lengths),
        "total_bp": sum(lengths),
        "n50": _n50(lengths),
        "longest": max(lengths) if lengths else 0,
    }


# ---------------------------------------------------------------------------
# Output formatters
# ---------------------------------------------------------------------------


def _write_json(stats: dict[str, Any], output_path: str | Path) -> None:
    with open(str(output_path), "w") as fh:
        json.dump(stats, fh, indent=2, default=str)
        fh.write("\n")


def _write_tsv(stats: dict[str, Any], output_path: str | Path) -> None:
    with open(str(output_path), "w") as fh:
        _flatten_to_tsv(stats, fh, prefix="")


def _flatten_to_tsv(d: dict[str, Any], fh: Any, prefix: str) -> None:
    for key, value in d.items():
        full_key = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, dict):
            _flatten_to_tsv(value, fh, full_key)
        else:
            fh.write(f"{full_key}\t{value}\n")


def _make_chart_base64(title: str, labels: list[str], values: list[int | float]) -> str:
    """Render a matplotlib bar chart to a base64-encoded PNG string."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return ""

    fig, ax = plt.subplots(figsize=(6, max(2.5, 0.4 * len(labels))))
    bars = ax.barh(range(len(labels)), values, color="#4a90d9")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.invert_yaxis()
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_width() + max(values) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            str(val) if isinstance(val, int) else f"{val:.2f}",
            va="center",
            fontsize=8,
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def _make_histogram_base64(
    title: str, values: list[int] | list[float], xlabel: str, bins: int = 30
) -> str:
    """Render a matplotlib histogram to a base64-encoded PNG string."""
    if not values:
        return ""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return ""

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.hist(values, bins=bins, color="#4a90d9", edgecolor="#2c5f8a")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def _write_html(
    stats: dict[str, Any], output_path: str | Path, genome_name: str
) -> None:
    sections: list[str] = []

    # Header
    sections.append(
        f"<h1>HelixForge QC Report</h1>\n"
        f"<p class='sub'>Genome: <b>{genome_name}</b> &mdash; "
        f"{stats['total_genes']} genes</p>"
    )

    # Summary statistics table
    summary_rows = [
        ("Total genes", stats["total_genes"]),
        ("Total transcripts", stats["total_transcripts"]),
        ("Multi-isoform genes", stats["multi_isoform_genes"]),
        ("Coding genes", stats["coding_genes"]),
    ]
    if "genome" in stats:
        gs = stats["genome"]
        summary_rows.extend(
            [
                ("Genome contigs", gs.get("contigs", "—")),
                ("Genome total bp", f"{gs.get('total_bp', 0):,}"),
                ("Genome N50", f"{gs.get('n50', 0):,}"),
            ]
        )
    sections.append(_html_table("Summary Statistics", summary_rows))

    # Tier distribution chart
    tier = stats.get("tier", {})
    if tier:
        chart = _make_chart_base64(
            "Tier Distribution",
            [f"Tier {k}" for k in tier],
            list(tier.values()),
        )
        sections.append(_html_chart_section("Tier Distribution", chart, tier))

    # Biotype distribution chart
    biotype = stats.get("biotype", {})
    if biotype:
        chart = _make_chart_base64(
            "Biotype Distribution", list(biotype.keys()), list(biotype.values())
        )
        sections.append(_html_chart_section("Biotype Distribution", chart, biotype))

    # CDS length histogram
    cds_info = stats.get("cds_length", {})
    if cds_info.get("count", 0) > 0:
        sections.append(_html_dist_table("CDS Length Distribution", cds_info))

    # Intron length distribution
    intron_info = stats.get("intron_length", {})
    if intron_info.get("count", 0) > 0:
        sections.append(_html_dist_table("Intron Length Distribution", intron_info))

    # Exon count distribution
    exon_info = stats.get("exon_count", {})
    if exon_info.get("count", 0) > 0:
        sections.append(_html_dist_table("Exon Count Distribution", exon_info))

    # Confidence score histogram
    if "confidence" in stats:
        conf_info = stats["confidence"]
        sections.append(_html_dist_table("Confidence Score Distribution", conf_info))

    # QC flags table
    flags = stats.get("flags", {})
    if flags:
        flag_rows = [(k, v) for k, v in flags.items()]
        sections.append(_html_table("QC Flags", flag_rows))

    # AS event types
    as_events = stats.get("as_events", {})
    if as_events:
        as_rows = [(k, v) for k, v in as_events.items()]
        sections.append(_html_table("Alternative Splicing Events", as_rows))

    html = _HTML_TEMPLATE.format(
        genome_name=genome_name,
        body="\n".join(sections),
    )
    with open(str(output_path), "w") as fh:
        fh.write(html)


def _html_table(title: str, rows: list[tuple[str, Any]]) -> str:
    body = "\n".join(f"<tr><td>{k}</td><td class='v'>{v}</td></tr>" for k, v in rows)
    return f"<h2>{title}</h2>\n<table>\n{body}\n</table>"


def _html_dist_table(title: str, info: dict[str, Any]) -> str:
    rows = [
        ("Count", info.get("count", 0)),
        ("Min", info.get("min", 0)),
        ("Median", info.get("median", 0)),
        ("Mean", info.get("mean", 0)),
        ("Max", info.get("max", 0)),
        ("N50", info.get("n50", 0)),
    ]
    return _html_table(title, rows)


def _html_chart_section(title: str, chart_b64: str, data: dict[str, Any]) -> str:
    parts = [f"<h2>{title}</h2>"]
    if chart_b64:
        parts.append(f'<img src="data:image/png;base64,{chart_b64}" alt="{title}">')
    rows = [(k, v) for k, v in data.items()]
    body = "\n".join(f"<tr><td>{k}</td><td class='v'>{v}</td></tr>" for k, v in rows)
    parts.append(f"<table>\n{body}\n</table>")
    return "\n".join(parts)


_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>HelixForge QC Report — {genome_name}</title>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #222; max-width: 900px; }}
 h1 {{ margin-bottom: 0; }}
 .sub {{ color: #777; margin-top: 0.2rem; }}
 table {{ border-collapse: collapse; margin: 0.5rem 0 1.5rem; min-width: 22rem; }}
 td {{ border: 1px solid #ddd; padding: 0.3rem 0.7rem; }}
 td.v {{ text-align: right; font-variant-numeric: tabular-nums; }}
 h2 {{ margin-top: 1.5rem; }}
 img {{ max-width: 100%; margin: 0.5rem 0; }}
</style>
</head>
<body>
{body}
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_qc_report(
    gff3_path: str | Path,
    output_path: str | Path,
    genome_path: str | Path | None = None,
    helixer_h5_path: str | Path | None = None,
    fmt: str = "html",
) -> dict[str, Any]:
    """Generate genome-wide QC report from a GFF3 annotation.

    Collects: gene count by tier/origin/biotype; transcript and multi-isoform
    counts; CDS length distribution (min, median, mean, max, N50); exon count
    distribution; intron length distribution; QC flag frequency; AS event type
    counts. Optionally adds confidence score distribution from Helixer HDF5.

    Returns:
        Statistics dict (same data regardless of output format).
    """
    stats = _collect_stats(gff3_path, genome_path, helixer_h5_path)

    genome_name = Path(gff3_path).stem

    if fmt == "json":
        _write_json(stats, output_path)
    elif fmt == "tsv":
        _write_tsv(stats, output_path)
    else:
        _write_html(stats, output_path, genome_name)

    return stats
