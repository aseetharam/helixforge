"""Every headline number in the poster docs traces to a committed result table.

No hand-typed numbers: each literal asserted to appear in ``README.md`` /
``docs/POSTER_NARRATIVE.md`` / ``docs/BENCHMARK.md`` is checked against the value
in a committed ``bench_out/`` / ``ablation_out/`` table (parse both, assert equal
within rounding). CLAUDE.md §16 F3/F5.

The committed source tables (small TSV/JSON, un-ignored in ``.gitignore``):
- ``bench_out_m7a/results_gffcompare_mikado_compare_compleasm_both.tsv`` — Helixer
  AND HelixForge accuracy + completeness (the before/after structural numbers).
- ``bench_out/iso_trace.result.json`` — TRaCE-on isoform table (junction, primary,
  alt precisions).
- ``bench_out/iso_noTRACE.result.json`` — the no-TRaCE primary-precision baseline.
- ``bench_out/benchmark.tsv`` — HelixForge master table (isoform_* rows).
- ``ablation_out/ablation.tsv`` — full vs no_helixer_support (the 13/13 figure).
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
README = ROOT / "README.md"
NARRATIVE = ROOT / "docs" / "POSTER_NARRATIVE.md"
BENCHMARK = ROOT / "docs" / "BENCHMARK.md"

# Tolerance: a doc literal may round a committed value to <=2 decimals.
TOL = 0.06


# --------------------------------------------------------------------------- #
# Committed-table parsers.
# --------------------------------------------------------------------------- #
def _require(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f"committed result table missing: {path.relative_to(ROOT)}")
    return path


def _both_tsv() -> dict[tuple[str, str, str], float]:
    """(annotation, tool, metric) -> value for the before/after table."""
    p = _require(ROOT / "bench_out_m7a" / "results_gffcompare_mikado_compare_compleasm_both.tsv")
    out: dict[tuple[str, str, str], float] = {}
    with p.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[(row["annotation"], row["tool"], row["metric"])] = float(row["value"])
            except (TypeError, ValueError):
                pass
    return out


def _benchmark_tsv() -> dict[tuple[str, str], float]:
    """(tool, metric) -> value for the HelixForge master table."""
    p = _require(ROOT / "bench_out" / "benchmark.tsv")
    out: dict[tuple[str, str], float] = {}
    with p.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[(row["tool"], row["metric"])] = float(row["value"])
            except (TypeError, ValueError):
                pass
    return out


def _iso(name: str) -> dict:
    p = _require(ROOT / "bench_out" / name)
    return json.loads(p.read_text())


def _ablation_rows() -> dict[str, dict[str, str]]:
    p = _require(ROOT / "ablation_out" / "ablation.tsv")
    with p.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["variant"]: r for r in rows}


# --------------------------------------------------------------------------- #
# Helpers.
# --------------------------------------------------------------------------- #
def _has_number(text: str, literal: str) -> bool:
    """``literal`` appears as a standalone number token (not a substring of a
    longer number, so '93.0' does not match '93.02')."""
    pat = r"(?<![\d.])" + re.escape(literal) + r"(?![\d])"
    return re.search(pat, text) is not None


def _check(doc_text: str, literal: str, source_value: float) -> None:
    assert _has_number(doc_text, literal), f"{literal!r} not found as a number in doc"
    assert abs(float(literal) - source_value) <= TOL, (
        f"doc literal {literal} does not match committed value {source_value}"
    )


# --------------------------------------------------------------------------- #
# Tests.
# --------------------------------------------------------------------------- #
def test_structural_before_after_trace_to_both_tsv():
    """intron F1 / locus Sn / gene F1>=80% before->after match the committed
    Helixer-and-HelixForge results table, in README and the narrative."""
    both = _both_tsv()
    readme, narr = README.read_text(), NARRATIVE.read_text()

    checks = [
        ("89.0", both[("helixer", "mikado_compare", "intron_f1")]),
        ("93.0", both[("helixforge", "mikado_compare", "intron_f1")]),
        ("66.9", both[("helixer", "gffcompare", "locus_sn")]),
        ("75.0", both[("helixforge", "gffcompare", "locus_sn")]),
        ("67.0", both[("helixer", "mikado_compare", "gene_80_base_f1_f1")]),
        ("74.7", both[("helixforge", "mikado_compare", "gene_80_base_f1_f1")]),
    ]
    for literal, value in checks:
        _check(readme, literal, value)
        _check(narr, literal, value)


def test_isoform_headline_numbers_trace():
    """Junction 93.47/92.57, primary 81.17->85.88, alt 41.67 — the required
    isoform headline numbers — match the committed iso_*.result.json tables."""
    trace = _iso("iso_trace.result.json")["matched"]
    notrace = _iso("iso_noTRACE.result.json")
    readme, narr, bench = README.read_text(), NARRATIVE.read_text(), BENCHMARK.read_text()

    # Junction P/R (terminus-independent, TRaCE-independent): README + narrative + bench.
    for doc in (readme, narr, bench):
        _check(doc, "93.47", trace["junction_precision"] * 100)
        _check(doc, "92.57", trace["junction_recall"] * 100)

    # Primary precision baseline -> TRaCE-on: narrative + bench.
    for doc in (narr, bench):
        _check(doc, "81.17", notrace["primary_precision"] * 100)
        _check(doc, "85.88", trace["primary_precision"] * 100)

    # Alt-isoform precision (vs Helixer's structural 0): all three.
    for doc in (readme, narr, bench):
        _check(doc, "41.67", trace["alt_isoform_precision"] * 100)


def test_benchmark_md_isoform_section_traces():
    """The BENCHMARK.md isoform section's sensitivity-headroom and count numbers
    trace to benchmark.tsv isoform_* rows."""
    bt = _benchmark_tsv()
    bench = BENCHMARK.read_text()
    _check(bench, "18.13", bt[("isoform", "isoform_intron_chain_sn_multiiso")])
    _check(bench, "44.97", bt[("isoform", "isoform_intron_chain_pr_multiiso")])
    # genes that gained isoforms (4,536) appears with the thousands comma stripped.
    gained = int(bt[("isoform", "isoform_count_genes_gained_isoforms")])
    assert gained == 4536
    assert "4,536" in bench and "4,536" in NARRATIVE.read_text()


def test_completeness_tradeoff_traces():
    """The honest completeness tradeoff (99.09 -> 97.38, Missing 0.67 -> 2.25)
    matches the committed compleasm rows and is stated in LIMITATIONS + BENCHMARK."""
    both = _both_tsv()
    bt = _benchmark_tsv()
    limits = (ROOT / "docs" / "LIMITATIONS.md").read_text()
    bench = BENCHMARK.read_text()
    for doc in (limits, bench):
        _check(doc, "99.09", both[("helixer", "compleasm", "complete")])
        _check(doc, "97.38", bt[("compleasm", "complete")])
        _check(doc, "0.67", both[("helixer", "compleasm", "missing")])
        _check(doc, "2.25", bt[("compleasm", "missing")])


def test_ablation_is_13_of_13_positive():
    """The '13/13' claim recomputes from ablation.tsv: full beats
    no_helixer_support on all 13 F1/completeness levels (0 negative)."""
    rows = _ablation_rows()
    full, nohs = rows["full"], rows["no_helixer_support"]

    # The 13 levels mirror scripts/make_figures.py::_ABLATION_LEVELS.
    f1_cols = [
        "mikado_compare.base_f1",
        "mikado_compare.exon_lenient_f1",
        "mikado_compare.splice_site_f1",
        "mikado_compare.intron_f1",
        "mikado_compare.intron_chain_f1",
        "mikado_compare.transcript_80_base_f1_f1",
        "mikado_compare.gene_80_base_f1_f1",
        "mikado_compare.gene_95_base_f1_f1",
        "compleasm.complete",
    ]
    snpr_pairs = [
        ("gffcompare.exon_sn", "gffcompare.exon_pr"),
        ("gffcompare.intron_chain_sn", "gffcompare.intron_chain_pr"),
        ("gffcompare.transcript_sn", "gffcompare.transcript_pr"),
        ("gffcompare.locus_sn", "gffcompare.locus_pr"),
    ]

    def _f1(sn: float, pr: float) -> float:
        return 0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr)

    deltas = []
    for col in f1_cols:
        deltas.append(float(full[col]) - float(nohs[col]))
    for sn_col, pr_col in snpr_pairs:
        full_f1 = _f1(float(full[sn_col]), float(full[pr_col]))
        nohs_f1 = _f1(float(nohs[sn_col]), float(nohs[pr_col]))
        deltas.append(full_f1 - nohs_f1)

    assert len(deltas) == 13
    n_positive = sum(1 for d in deltas if d > 0)
    assert n_positive == 13, f"expected 13/13 positive, got {n_positive}/13"

    # The docs state exactly "13/13".
    assert "13/13" in README.read_text()
    assert "13/13" in NARRATIVE.read_text()
