"""Tests for ``scripts/make_figures.py`` (poster step 3).

Hermetic: every test builds a tiny **fixture** result table holding the exact
poster numbers and asserts the figure-extraction reads those values back and the
figure renders to vector (SVG/PDF). No dependency on the large gitignored
``bench_out/`` / ``m3_out/`` tables, so the suite is green anywhere.

CLAUDE.md §12: concrete literal coordinates/numbers; non-interactive backend.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")  # non-interactive backend (test floor).

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "make_figures.py"


def _load():
    spec = importlib.util.spec_from_file_location("make_figures", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mf = _load()


# ---------------------------------------------------------------------------
# Fixture builders — tiny tables with the literal poster numbers.
# ---------------------------------------------------------------------------


def _before_after_tsv(p: Path) -> Path:
    rows = [
        ("helixer", "mikado_compare", "intron_f1", "89.04"),
        ("helixforge", "mikado_compare", "intron_f1", "92.97"),
        ("helixer", "gffcompare", "locus_sn", "66.9"),
        ("helixforge", "gffcompare", "locus_sn", "75.0"),
        ("helixer", "mikado_compare", "gene_80_base_f1_f1", "67.05"),
        ("helixforge", "mikado_compare", "gene_80_base_f1_f1", "74.67"),
        ("helixer", "compleasm", "complete", "99.09"),
        ("helixforge", "compleasm", "complete", "97.38"),
    ]
    lines = ["annotation\ttool\tmetric\tvalue\tstatus"]
    lines += [f"{a}\t{t}\t{m}\t{v}\tok" for a, t, m, v in rows]
    p.write_text("\n".join(lines) + "\n")
    return p


def _ablation_tsv(p: Path) -> Path:
    # Only the columns extract_ablation reads; full beats no_helixer_support on all.
    cols = [
        "variant",
        "mikado_compare.base_f1", "mikado_compare.exon_lenient_f1",
        "mikado_compare.splice_site_f1", "mikado_compare.intron_f1",
        "mikado_compare.intron_chain_f1", "mikado_compare.transcript_80_base_f1_f1",
        "mikado_compare.gene_80_base_f1_f1", "mikado_compare.gene_95_base_f1_f1",
        "gffcompare.exon_sn", "gffcompare.exon_pr",
        "gffcompare.intron_chain_sn", "gffcompare.intron_chain_pr",
        "gffcompare.transcript_sn", "gffcompare.transcript_pr",
        "gffcompare.locus_sn", "gffcompare.locus_pr",
        "compleasm.complete",
    ]
    full = ["full", "86.86", "84.56", "92.05", "92.97", "55.56", "55.05",
            "74.67", "44.03", "74.1", "81.0", "47.5", "65.9", "47.7", "65.7",
            "75.1", "78.0", "97.38"]
    nohs = ["no_helixer_support", "86.74", "84.42", "92.0", "92.9", "55.0",
            "54.56", "74.01", "43.74", "74.0", "80.8", "47.0", "65.2", "47.2",
            "65.3", "74.4", "77.6", "96.96"]
    p.write_text("\t".join(cols) + "\n" + "\t".join(full) + "\n" + "\t".join(nohs) + "\n")
    return p


def _iso_trace_json(p: Path) -> Path:
    p.write_text(json.dumps({"matched": {
        "junction_precision": 0.9347330948845998,
        "junction_recall": 0.9257446008605205,
        "primary_precision": 0.8588481050432707,
        "alt_isoform_precision": 0.4166666666666667,
    }}))
    return p


def _iso_notrace_json(p: Path) -> Path:
    # Written directly (no "matched" wrapper), as scripts compute the baseline.
    p.write_text(json.dumps({
        "primary_precision": 0.8116980005968367,
        "junction_precision": 0.9347330948845998,
        "junction_recall": 0.9257446008605205,
        "alt_isoform_precision": 0.4434826883910387,
    }))
    return p


def _benchmark_tsv(p: Path) -> Path:
    rows = [
        ("isoform", "isoform_count_helixer_mean", "1.0"),
        ("isoform", "isoform_count_helixforge_mean", "1.2788567645111455"),
        ("isoform", "isoform_count_reference_mean", "1.7521605496293617"),
        ("agat", "number_of_gene", "27184.0"),  # ignored by the extractor
    ]
    lines = ["tool\tmetric\tvalue\tstatus"]
    lines += [f"{t}\t{m}\t{v}\tok" for t, m, v in rows]
    p.write_text("\n".join(lines) + "\n")
    return p


def _report_json(p: Path) -> Path:
    p.write_text(json.dumps({"tier": {"1": 20631, "2": 5565, "3": 990}}))
    return p


def _is_vector(path: Path) -> bool:
    head = path.read_bytes()[:8]
    if path.suffix == ".svg":
        return b"<svg" in path.read_bytes()[:600] or head.startswith(b"<?xml")
    return head.startswith(b"%PDF")


# ---------------------------------------------------------------------------
# 1. before/after values == source rows.
# ---------------------------------------------------------------------------


def test_before_after_values_match_source(tmp_path):
    data = mf.extract_before_after(_before_after_tsv(tmp_path / "ba.tsv"))
    assert data == [
        ("Intron F1", 89.04, 92.97, False),
        ("Locus Sn", 66.9, 75.0, False),
        ("Gene F1 ≥80%", 67.05, 74.67, False),
        ("BUSCO complete", 99.09, 97.38, True),  # tradeoff flagged
    ]


# ---------------------------------------------------------------------------
# 2. ablation: 13 levels, all positive deltas (the 13/13 result).
# ---------------------------------------------------------------------------


def test_ablation_13_of_13_positive(tmp_path):
    data = mf.extract_ablation(_ablation_tsv(tmp_path / "abl.tsv"))
    assert len(data) == 13
    assert all(delta > 0 for _, delta, _, _ in data)
    by = {label: (delta, fv, nv) for label, delta, fv, nv in data}
    # mikado F1 read directly:
    assert by["Gene ≥80%"] == (0.66, 74.67, 74.01)
    # gffcompare F1 computed as harmonic mean of its Sn/Pr columns:
    assert by["Locus (gff)"] == (0.55, 76.52, 75.97)
    assert by["Completeness"][1] == 97.38  # full completeness


# ---------------------------------------------------------------------------
# 3. the three isoform panels (junction 93.47/92.57, primary 81.17/85.88,
#    alt 41.67 vs 0).
# ---------------------------------------------------------------------------


def test_isoform_panels_values(tmp_path):
    panels = mf.extract_isoform_panels(
        _iso_trace_json(tmp_path / "t.json"),
        _iso_notrace_json(tmp_path / "n.json"),
    )
    assert panels["junction"] == (93.47, 92.57)
    assert panels["primary"] == (81.17, 85.88)  # before (no-TRaCE) -> after (TRaCE)
    assert panels["alt"] == (41.67, 0.0)  # HelixForge vs Helixer's structural zero


def test_isoform_count_and_tier_values(tmp_path):
    counts = mf.extract_isoform_counts(_benchmark_tsv(tmp_path / "b.tsv"))
    assert counts[0] == ("Helixer", 1.0)
    assert counts[1][0] == "HelixForge" and round(counts[1][1], 2) == 1.28
    assert counts[2][0] == "Araport11" and round(counts[2][1], 2) == 1.75
    tiers = mf.extract_tiers(_report_json(tmp_path / "r.json"))
    assert tiers == [("Tier 1", 20631), ("Tier 2", 5565), ("Tier 3", 990)]


# ---------------------------------------------------------------------------
# 4. each data-driven figure writes vector SVG + PDF from its table.
# ---------------------------------------------------------------------------


def test_figures_write_vector_svg_and_pdf(tmp_path):
    out = tmp_path / "figs"
    written: list[Path] = []
    written += mf.figure_before_after(out, _before_after_tsv(tmp_path / "ba.tsv"))
    written += mf.figure_ablation(out, _ablation_tsv(tmp_path / "abl.tsv"))
    written += mf.figure_isoform_panels(
        out, _iso_trace_json(tmp_path / "t.json"), _iso_notrace_json(tmp_path / "n.json"))
    written += mf.figure_isoform_counts(out, _benchmark_tsv(tmp_path / "b.tsv"))
    written += mf.figure_tiers(out, _report_json(tmp_path / "r.json"))

    suffixes = {p.suffix for p in written}
    assert suffixes == {".svg", ".pdf"}
    for p in written:
        assert p.is_file() and p.stat().st_size > 0
        assert _is_vector(p), f"{p} is not vector"


# ---------------------------------------------------------------------------
# 5. a missing source table fails with a clear message — not a blank figure.
# ---------------------------------------------------------------------------


def test_missing_source_table_raises_clear_error(tmp_path):
    out = tmp_path / "figs"
    missing = tmp_path / "does_not_exist.tsv"
    with pytest.raises(FileNotFoundError) as exc:
        mf.figure_before_after(out, missing)
    msg = str(exc.value)
    assert "before/after" in msg and str(missing) in msg
    # No blank figure left behind.
    assert not out.exists() or not any(out.glob("before_after_benchmark.*"))


def test_missing_isoform_json_raises_clear_error(tmp_path):
    with pytest.raises(FileNotFoundError) as exc:
        mf.extract_isoform_panels(tmp_path / "nope_trace.json", tmp_path / "nope_nt.json")
    assert "isoform-accuracy" in str(exc.value)


# ---------------------------------------------------------------------------
# 6. method schematic copies the authored SVG (vector).
# ---------------------------------------------------------------------------


def test_method_schematic_copies_authored_svg(tmp_path):
    src = tmp_path / "schem.src.svg"
    src.write_text('<?xml version="1.0"?>\n<svg xmlns="http://www.w3.org/2000/svg"></svg>\n')
    out = tmp_path / "figs"
    written = mf.figure_method_schematic(out, src)
    svg = out / "method_schematic.svg"
    assert svg in written and svg.is_file()
    assert _is_vector(svg)


# ---------------------------------------------------------------------------
# 7. example-locus figure (integration — needs viz deps + real GFF3s).
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_example_locus_figure_marks_canonical(tmp_path):
    pytest.importorskip("matplotlib")
    gff3 = mf.DEFAULT_SOURCES["reconciled_gff3"]
    if not Path(gff3).is_file():
        pytest.skip("real reconciled GFF3 not present")
    # Parse the chosen example gene and confirm the primary (.1) is the canonical.
    gene = mf._parse_gff3_gene(Path(gff3), mf.EXAMPLE_GENE)
    assert gene is not None
    assert gene.primary_transcript_id == f"{mf.EXAMPLE_GENE}.1"
    assert len(gene.transcripts) >= 2  # a multi-isoform locus
    written = mf.figure_example_locus(out_dir=tmp_path, reconciled_gff3=gff3,
                                      helixer_gff3=mf.DEFAULT_SOURCES["helixer_gff3"])
    assert any(p.suffix == ".svg" and p.is_file() for p in written)
    for p in written:
        assert _is_vector(p)
