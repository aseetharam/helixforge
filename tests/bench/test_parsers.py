"""Phase 17 D2 — benchmark stat parsers against captured real-format samples.

The samples under ``tests/bench/data/`` are tiny, redacted captures of the
*actual* tool output formats (not mocks-of-mocks). Each parser must pull the
right numbers from them and survive malformed/empty input. Floor: 12.
"""

from pathlib import Path

import pytest

from helixforge.bench.wrappers import (
    parse_agat_stats,
    parse_busco_summary,
    parse_compleasm_summary,
    parse_gffcompare_stats,
    parse_mikado_compare_stats,
    parse_omark_summary,
)

DATA = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# mikado compare
# ---------------------------------------------------------------------------

def test_mikado_compare_real_levels():
    out = parse_mikado_compare_stats(DATA / "mikado_compare.stats")
    assert out["base"] == {"sn": 78.34, "pr": 85.21, "f1": 81.63}
    assert out["intron_chain"] == {"sn": 55.20, "pr": 60.13, "f1": 57.56}
    assert out["gene"]["f1"] == 71.29


def test_mikado_compare_transcript_variants_distinct_keys():
    out = parse_mikado_compare_stats(DATA / "mikado_compare.stats")
    # Real mikado has several "Transcript level (...)" rows → distinct slugs.
    assert out["transcript_stringent"]["sn"] == 40.11
    assert "transcript_80_base_f1" in out
    assert out["exon_stringent"]["pr"] == 70.45


def test_mikado_compare_empty(tmp_path):
    p = tmp_path / "empty.stats"
    p.write_text("")
    assert parse_mikado_compare_stats(p) == {}


# ---------------------------------------------------------------------------
# gffcompare
# ---------------------------------------------------------------------------

def test_gffcompare_real_levels():
    out = parse_gffcompare_stats(DATA / "gffcompare.stats")
    assert out["base"] == {"sn": 78.3, "pr": 85.2}
    assert out["intron_chain"] == {"sn": 55.2, "pr": 60.1}
    assert out["transcript"]["pr"] == 64.4
    assert out["locus"]["sn"] == 70.5


def test_gffcompare_ignores_non_level_lines():
    out = parse_gffcompare_stats(DATA / "gffcompare.stats")
    # "Matching intron chains:" etc. must not become level rows.
    assert "matching_intron_chains" not in out
    assert set(out) == {"base", "exon", "intron", "intron_chain", "transcript", "locus"}


# ---------------------------------------------------------------------------
# compleasm
# ---------------------------------------------------------------------------

def test_compleasm_real():
    out = parse_compleasm_summary(DATA / "compleasm_summary.txt")
    assert out["single"] == 96.86
    assert out["duplicated"] == 1.42
    assert out["missing"] == 1.10
    assert out["n"] == 5004
    assert out["complete"] == pytest.approx(98.28)


def test_compleasm_malformed_lines_skipped(tmp_path):
    p = tmp_path / "summary.txt"
    p.write_text("garbage line\nS:50.0%, 1\nnot a metric\n")
    out = parse_compleasm_summary(p)
    assert out["single"] == 50.0
    assert "complete" not in out  # no duplicated line → not derivable


# ---------------------------------------------------------------------------
# BUSCO
# ---------------------------------------------------------------------------

def test_busco_real():
    out = parse_busco_summary(DATA / "busco_short_summary.txt")
    assert out["complete"] == 98.1
    assert out["single"] == 96.4
    assert out["duplicated"] == 1.7
    assert out["fragmented"] == 0.5
    assert out["missing"] == 1.4
    assert out["n"] == 4596


def test_busco_empty_returns_empty(tmp_path):
    p = tmp_path / "short_summary.txt"
    p.write_text("# no results line here\n")
    assert parse_busco_summary(p) == {}


# ---------------------------------------------------------------------------
# OMArk (bracketed real format) — completeness + consistency + spurious signal
# ---------------------------------------------------------------------------

def test_omark_real_bracketed_completeness():
    out = parse_omark_summary(DATA / "omark.sum")
    assert out["single"] == 90.40
    assert out["duplicated"] == 5.10
    assert out["missing"] == 4.50
    assert out["complete"] == pytest.approx(95.50)


def test_omark_real_consistency_and_fragmentation():
    out = parse_omark_summary(DATA / "omark.sum")
    assert out["consistent"] == 85.00
    assert out["inconsistent"] == 10.00
    # spurious-isoform-relevant signals (directly relevant to the isoform claim).
    assert out["fragmented"] == 3.00
    assert out["partial_hits"] == 6.00


def test_omark_flat_format_still_parses(tmp_path):
    # Backward-compat with the older flat "Single:90.00% (n)" form.
    p = tmp_path / "proteins.sum"
    p.write_text("Single:90.00% (13042)\nConsistent:85.00% (29750)\n")
    out = parse_omark_summary(p)
    assert out["single"] == 90.0
    assert out["consistent"] == 85.0


# ---------------------------------------------------------------------------
# AGAT
# ---------------------------------------------------------------------------

def test_agat_real_counts():
    out = parse_agat_stats(DATA / "agat_statistics.txt")
    assert out["number_of_gene"] == 27184.0
    assert out["number_of_mrna"] == 34765.0
    assert out["number_of_single_exon_gene"] == 3000.0


def test_agat_real_means_and_labels_with_units():
    out = parse_agat_stats(DATA / "agat_statistics.txt")
    assert out["mean_exons_per_mrna"] == 5.2
    assert out["mean_cds_length_bp"] == 1235.0  # label "(bp)" preserved through slug


def test_agat_skips_section_headers(tmp_path):
    p = tmp_path / "agat.txt"
    p.write_text("----------\nCompute mrna with isoforms\nNumber of gene    100\n")
    out = parse_agat_stats(p)
    assert out == {"number_of_gene": 100.0}
