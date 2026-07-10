"""Genome-level annotation + evidence-support report."""

from __future__ import annotations

import json
import math
from statistics import mean, median
from typing import TYPE_CHECKING, Any

from helixforge.stats.before_after import (
    _junction_support_from_set,
    annotation_summary,
)
from helixforge.stats.evidence_concordance import compute_aed
from helixforge.stats.summary import summarize_genes
from helixforge.utils.atomic import atomic_write
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene, SpliceJunction

_log = get_logger(__name__)

NAN = float("nan")


def _package_version() -> str | None:
    try:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("helixforge")
        except PackageNotFoundError:
            return None
    except Exception:  # noqa: BLE001 - version is decorative; never sink the report
        return None


def _has_utr(tx: Any) -> bool:
    """True if the transcript has UTR (exon span extends beyond the CDS span)."""
    if not tx.cds:
        return False
    cds_lo = min(c.start for c in tx.cds)
    cds_hi = max(c.end for c in tx.cds)
    ex_lo = min(e.start for e in tx.exons)
    ex_hi = max(e.end for e in tx.exons)
    return bool(ex_lo < cds_lo or ex_hi > cds_hi)


def _quantile(values: list[float], q: float) -> float:
    """Linear-interpolation quantile of a *sorted* list (q ∈ [0, 1])."""
    if not values:
        return NAN
    if len(values) == 1:
        return values[0]
    pos = q * (len(values) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return values[lo]
    frac = pos - lo
    return values[lo] * (1 - frac) + values[hi] * frac


# ---------------------------------------------------------------------------
# D4 — evidence-support / mapping-rate summary
# ---------------------------------------------------------------------------


def support_summary(
    genes: list[ReconciledGene],
    *,
    junctions: list[SpliceJunction] | None = None,
    min_reads: int = 3,
    tpm_threshold: float = 0.0,
    bam_stats: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """The "how trustworthy is this annotation" summary.

    - ``mapping_rate`` — overall RNA-seq mapping rate from ``bam_stats``
      (``{mapped, unmapped}`` or a precomputed ``mapping_rate``); ``None`` if no
      BAM stats are given.
    - ``frac_genes_junction_supported`` — fraction of genes whose primary
      transcript has ≥1 junction-supported intron. Uses the model's stored
      ``junction_support_fraction``; if that is ``None`` and a ``junctions`` list
      is supplied, it is computed against the read-count-filtered junction set.
      Single-exon genes (no introns) count as *not* junction-supported.
    - ``frac_genes_tpm_pass`` — fraction of genes whose primary transcript has
      ``tpm >= tpm_threshold`` (genes with no TPM are not counted as passing).
    - ``aed`` — distribution (n / mean / median / q1 / q3 / min / max) of the
      per-gene representative AED (:func:`compute_aed`).

    All fractions are over the full gene set so they are directly hand-checkable.
    """
    n = len(genes)

    # mapping rate
    mapping_rate: float | None = None
    if bam_stats is not None:
        if bam_stats.get("mapping_rate") is not None:
            mapping_rate = float(bam_stats["mapping_rate"])
        else:
            mapped = float(bam_stats.get("mapped", 0))
            unmapped = float(bam_stats.get("unmapped", 0))
            total = mapped + unmapped
            mapping_rate = (mapped / total) if total else 0.0

    junction_set = None
    if junctions:
        junction_set = {
            (j.seqid, j.donor, j.acceptor, j.strand)
            for j in junctions
            if j.read_count >= min_reads
        }

    junction_supported = 0
    tpm_pass = 0
    aeds: list[float] = []
    for g in genes:
        p = _primary(g)
        # junction support
        jsf = p.junction_support_fraction
        if jsf is None and junction_set is not None:
            jsf = _junction_support_from_set(g.seqid, g.strand, p.exons, junction_set)
        if jsf is not None and jsf > 0:
            junction_supported += 1
        # TPM pass
        if p.tpm is not None and p.tpm >= tpm_threshold:
            tpm_pass += 1
        # AED
        a = _gene_aed(g, p, junction_set)
        if a is not None:
            aeds.append(a)

    aeds_sorted = sorted(aeds)
    aed_dist = {
        "n": len(aeds_sorted),
        "mean": mean(aeds_sorted) if aeds_sorted else NAN,
        "median": median(aeds_sorted) if aeds_sorted else NAN,
        "q1": _quantile(aeds_sorted, 0.25),
        "q3": _quantile(aeds_sorted, 0.75),
        "min": aeds_sorted[0] if aeds_sorted else NAN,
        "max": aeds_sorted[-1] if aeds_sorted else NAN,
    }

    return {
        "gene_count": n,
        "mapping_rate": mapping_rate,
        "genes_junction_supported": junction_supported,
        "frac_genes_junction_supported": (junction_supported / n) if n else NAN,
        "genes_tpm_pass": tpm_pass,
        "tpm_threshold": tpm_threshold,
        "frac_genes_tpm_pass": (tpm_pass / n) if n else NAN,
        "aed": aed_dist,
    }


def _primary(gene: ReconciledGene) -> Any:
    for t in gene.transcripts:
        if t.transcript_id == gene.primary_transcript_id:
            return t
    return gene.transcripts[0]


def _gene_aed(
    gene: ReconciledGene,
    primary: Any,
    junction_set: set[tuple[str, int, int, str]] | None,
) -> float | None:
    """Representative AED for a gene from its stored evidence fields."""
    junction = primary.junction_support_fraction
    if junction is None and junction_set is not None:
        junction = _junction_support_from_set(
            gene.seqid, gene.strand, primary.exons, junction_set
        )
    expression = None if primary.tpm is None else (1.0 if primary.tpm > 0 else 0.0)
    return compute_aed(junction, expression, None)


# ---------------------------------------------------------------------------
# D3 — genome-level report
# ---------------------------------------------------------------------------


def build_report(
    genes: list[ReconciledGene],
    *,
    genome: Any = None,
    h5_reader: Any = None,
    junctions: list[SpliceJunction] | None = None,
    min_reads: int = 3,
    tpm_threshold: float = 0.0,
    bam_stats: dict[str, Any] | None = None,
    completeness: dict[str, Any] | None = None,
    run_stats: dict[str, Any] | None = None,
    protein_lengths: dict[str, int] | None = None,
) -> dict[str, Any]:
    """Assemble the genome-level annotation + support report dict (D3/D4).

    ``genome`` / ``h5_reader`` / ``junctions`` enrich the structural section via
    :func:`annotation_summary` (codon-accurate completeness, AED, Helixer
    support) when present; with none of them the structural numbers fall back to
    structural proxies. ``completeness`` is a pre-computed
    ``{busco, compleasm, omark}`` dict (the bench tools are run by the caller, not
    here — they are off the default path). ``run_stats`` is ``RunStats.as_dict()``.

    The returned dict is JSON-serialisable and deterministic (sorted keys on
    serialisation). It also carries a ``multiqc`` block (custom-content schema) so
    the file ingests directly into a MultiQC dashboard.
    """
    junction_set = None
    if junctions:
        junction_set = {
            (j.seqid, j.donor, j.acceptor, j.strand)
            for j in junctions
            if j.read_count >= min_reads
        }

    summ = annotation_summary(
        genes,
        genome=genome,
        h5_reader=h5_reader,
        junction_set=junction_set,
        protein_lengths=protein_lengths,
    )
    counts = summarize_genes(genes)

    # exons-per-transcript + % with UTR (not in annotation_summary)
    all_exon_counts = [len(t.exons) for g in genes for t in g.transcripts]
    coding = counts["coding"]
    utr_genes = sum(1 for g in genes if _primary(g).cds and _has_utr(_primary(g)))
    total_exons = sum(all_exon_counts)

    structure = {
        "gene_count": summ["gene_count"],
        "transcript_count": summ["transcript_count"],
        "exon_count": total_exons,
        "coding_gene_count": summ["coding_gene_count"],
        "mono_exon_genes": summ["mono_exon_genes"],
        "multi_exon_genes": summ["multi_exon_genes"],
        "frac_multi_exonic": (
            summ["multi_exon_genes"] / summ["gene_count"] if summ["gene_count"] else NAN
        ),
        "mean_exons_per_transcript": mean(all_exon_counts) if all_exon_counts else NAN,
        "mean_cds_length": summ["mean_cds_length"],
        "median_cds_length": summ["median_cds_length"],
        "isoforms_per_gene_mean": summ["isoforms_per_gene_mean"],
        "isoforms_per_gene_median": summ["isoforms_per_gene_median"],
        "isoforms_per_gene_max": summ["isoforms_per_gene_max"],
        "isoforms_per_gene_distribution": {
            str(k): v for k, v in sorted(summ["isoforms_per_gene_distribution"].items())
        },
        "pct_with_utr": (100.0 * utr_genes / coding) if coding else NAN,
        "pct_complete_orfs": summ["pct_complete_orfs"],
        "mean_aed": summ["mean_aed"],
        "mean_helixer_support": summ["mean_helixer_support"],
    }

    # EXPRESSED / LOW / SILENT
    from collections import Counter

    status = Counter(g.classification.status for g in genes)
    n = len(genes)
    classification = {
        "counts": dict(sorted(status.items())),
        "fractions": {
            k: (status.get(k, 0) / n if n else NAN)
            for k in ("EXPRESSED", "LOW", "SILENT")
        },
    }

    support = support_summary(
        genes,
        junctions=junctions,
        min_reads=min_reads,
        tpm_threshold=tpm_threshold,
        bam_stats=bam_stats,
    )

    report: dict[str, Any] = {
        "tool": "HelixForge",
        "version": _package_version(),
        "structure": structure,
        "biotype": counts["biotype"],
        "tier": {str(k): v for k, v in counts["tier"].items()},
        "origin": counts["origin"],
        "flags": counts["flags"],
        "as_kinds": counts["as_kinds"],
        "classification": classification,
        "completeness": completeness,
        "support": support,
        "run_stats": run_stats,
    }
    report["multiqc"] = _multiqc_block(report)
    return report


def _multiqc_block(report: dict[str, Any]) -> dict[str, Any]:
    """A MultiQC custom-content block summarising the headline numbers.

    MultiQC ingests a JSON whose ``data`` maps a sample id → ``{metric: value}``;
    ``plot_type='generalstats'`` surfaces these in the general-stats table. We
    expose the scalar headline metrics (NaNs dropped — MultiQC dislikes them).
    """
    s = report["structure"]
    sup = report["support"]
    flat = {
        "genes": s["gene_count"],
        "transcripts": s["transcript_count"],
        "coding_genes": s["coding_gene_count"],
        "mean_cds_length": s["mean_cds_length"],
        "isoforms_per_gene_mean": s["isoforms_per_gene_mean"],
        "frac_multi_exonic": s["frac_multi_exonic"],
        "pct_with_utr": s["pct_with_utr"],
        "pct_complete_orfs": s["pct_complete_orfs"],
        "mean_aed": s["mean_aed"],
        "mapping_rate": sup["mapping_rate"],
        "frac_genes_junction_supported": sup["frac_genes_junction_supported"],
    }
    data = {k: v for k, v in flat.items() if v is not None and not _is_nan(v)}
    return {
        "id": "helixforge_annotation",
        "section_name": "HelixForge Annotation",
        "plot_type": "generalstats",
        "data": {"HelixForge": data},
    }


def _is_nan(v: Any) -> bool:
    return isinstance(v, float) and math.isnan(v)


# ---------------------------------------------------------------------------
# Serialisation — JSON (MultiQC-compatible) + self-contained HTML
# ---------------------------------------------------------------------------


def _json_safe(value: Any) -> Any:
    """Recursively replace NaN floats with None so the JSON is valid + diffable."""
    if isinstance(value, float):
        return None if math.isnan(value) else value
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def write_report_json(report: dict[str, Any], path: str) -> str:
    """Write ``report`` as deterministic JSON (sorted keys, NaN→null). Returns path."""
    with atomic_write(path) as fh:
        json.dump(_json_safe(report), fh, indent=2, sort_keys=True)
        fh.write("\n")
    return str(path)


def _fmt(v: Any) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        if math.isnan(v):
            return "—"
        return f"{v:.4g}"
    return str(v)


def _table(title: str, rows: list[tuple[str, Any]]) -> str:
    body = "\n".join(
        f"<tr><td>{k}</td><td class='v'>{_fmt(v)}</td></tr>" for k, v in rows
    )
    return f"<h2>{title}</h2>\n<table>\n{body}\n</table>"


def write_report_html(report: dict[str, Any], path: str) -> str:
    """Render ``report`` to a self-contained ``report.html`` (deterministic). Returns path."""
    s = report["structure"]
    sup = report["support"]
    cls = report["classification"]
    sections: list[str] = []
    sections.append(
        _table(
            "Structure",
            [
                ("Genes", s["gene_count"]),
                ("Transcripts (mRNA)", s["transcript_count"]),
                ("Exons", s["exon_count"]),
                ("Coding genes", s["coding_gene_count"]),
                ("Mono-exonic genes", s["mono_exon_genes"]),
                ("Multi-exonic genes", s["multi_exon_genes"]),
                ("Fraction multi-exonic", s["frac_multi_exonic"]),
                ("Mean exons / transcript", s["mean_exons_per_transcript"]),
                ("Mean CDS length", s["mean_cds_length"]),
                ("Median CDS length", s["median_cds_length"]),
                ("Isoforms / gene (mean)", s["isoforms_per_gene_mean"]),
                ("Isoforms / gene (max)", s["isoforms_per_gene_max"]),
                ("% genes with UTR", s["pct_with_utr"]),
                ("% complete ORFs", s["pct_complete_orfs"]),
                ("Mean AED", s["mean_aed"]),
                ("Mean Helixer support", s["mean_helixer_support"]),
            ],
        )
    )
    sections.append(_table("Biotype", sorted(report["biotype"].items())))
    sections.append(_table("Tier", sorted(report["tier"].items())))
    sections.append(_table("Origin", sorted(report["origin"].items())))
    sections.append(
        _table(
            "Expression class",
            [
                (k, f"{cls['counts'].get(k, 0)} ({_fmt(100.0 * cls['fractions'][k])}%)")
                for k in ("EXPRESSED", "LOW", "SILENT")
            ],
        )
    )
    sections.append(
        _table(
            "Evidence support",
            [
                ("RNA-seq mapping rate", sup["mapping_rate"]),
                ("Genes junction-supported", sup["genes_junction_supported"]),
                ("Fraction junction-supported", sup["frac_genes_junction_supported"]),
                (f"Genes TPM ≥ {sup['tpm_threshold']}", sup["genes_tpm_pass"]),
                ("Fraction TPM-pass", sup["frac_genes_tpm_pass"]),
                ("AED mean", sup["aed"]["mean"]),
                ("AED median", sup["aed"]["median"]),
            ],
        )
    )
    if report.get("completeness"):
        comp_rows = sorted(_flatten_scalars(report["completeness"]).items())
        sections.append(_table("Completeness (BUSCO / compleasm / OMArk)", comp_rows))
    if report.get("run_stats"):
        rs_rows = sorted(_flatten_scalars(report["run_stats"]).items())
        sections.append(_table("Run decision telemetry", rs_rows))

    html = _HTML_TEMPLATE.format(
        version=_fmt(report.get("version")), body="\n".join(sections)
    )
    with atomic_write(path) as fh:
        fh.write(html)
    return str(path)


def _flatten_scalars(d: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    """Flatten one level of nested dicts into ``{key: scalar}`` for HTML tables."""
    out: dict[str, Any] = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten_scalars(v, prefix=f"{key}."))
        else:
            out[key] = v
    return out


_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>HelixForge annotation report</title>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #222; }}
 h1 {{ margin-bottom: 0; }}
 .sub {{ color: #777; margin-top: 0.2rem; }}
 table {{ border-collapse: collapse; margin: 0.5rem 0 1.5rem; min-width: 22rem; }}
 td {{ border: 1px solid #ddd; padding: 0.3rem 0.7rem; }}
 td.v {{ text-align: right; font-variant-numeric: tabular-nums; }}
 h2 {{ margin-top: 1.5rem; }}
</style>
</head>
<body>
<h1>HelixForge annotation report</h1>
<p class="sub">HelixForge {version}</p>
{body}
</body>
</html>
"""
