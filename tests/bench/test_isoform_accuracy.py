"""Poster step 2 — isoform-level accuracy metric (floor: 8).

Synthetic GFF3 fixtures + the captured ``mikado compare`` ``.stats`` sample under
``tests/bench/data/`` (real format). The external tool is mocked everywhere
except the ``integration``-marked end-to-end, which self-skips without ``mikado``.
Both strands are exercised for the multi-isoform restriction and the subset
writer (CLAUDE.md §12 — minus-strand tests are mandatory).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from helixforge.bench import isoform as iso
from helixforge.bench.isoform import (
    _multiiso_genes,
    _parse_genes,
    _summary,
    _transcript_counts,
    _write_genes_gff3,
    isoform_accuracy,
    isoform_result_to_rows,
)
from helixforge.bench.wrappers import BenchmarkError, parse_mikado_compare_stats

DATA = Path(__file__).parent / "data"


# ---------------------------------------------------------------------------
# Synthetic GFF3 builders (1-based inclusive, the on-disk convention)
# ---------------------------------------------------------------------------


def _gene_block(gid: str, seqid: str, strand: str, transcripts: list[tuple[str, list[tuple[int, int]]]]) -> list[str]:
    """One gene with the given transcripts (each: tid, list of 1-based exon spans)."""
    lo = min(s for _, exs in transcripts for s, _ in exs)
    hi = max(e for _, exs in transcripts for _, e in exs)
    out = [f"{seqid}\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}"]
    for tid, exons in transcripts:
        t_lo = min(s for s, _ in exons)
        t_hi = max(e for _, e in exons)
        out.append(f"{seqid}\ttest\tmRNA\t{t_lo}\t{t_hi}\t.\t{strand}\t.\tID={tid};Parent={gid}")
        for i, (s, e) in enumerate(sorted(exons), start=1):
            out.append(f"{seqid}\ttest\texon\t{s}\t{e}\t.\t{strand}\t.\tID={tid}.e{i};Parent={tid}")
    return out


def _write_gff3(path: Path, blocks: list[list[str]]) -> Path:
    path.write_text("##gff-version 3\n" + "\n".join(ln for b in blocks for ln in b) + "\n")
    return path


def _plus_mixed(tmp_path: Path) -> Path:
    """Two single-isoform genes + two multi-isoform genes, all on + strand."""
    return _write_gff3(tmp_path / "plus.gff3", [
        _gene_block("g1", "chr1", "+", [("g1.1", [(1000, 1200)])]),
        _gene_block("g2", "chr1", "+", [
            ("g2.1", [(2000, 2200), (2400, 2600)]),
            ("g2.2", [(2000, 2200), (2450, 2600)]),
        ]),
        _gene_block("g3", "chr1", "+", [("g3.1", [(3000, 3300)])]),
        _gene_block("g4", "chr1", "+", [
            ("g4.1", [(4000, 4200), (4400, 4600), (4800, 5000)]),
            ("g4.2", [(4000, 4200), (4800, 5000)]),
            ("g4.3", [(4000, 4250), (4400, 4600)]),
        ]),
    ])


def _minus_mixed(tmp_path: Path) -> Path:
    """One single-isoform gene + one multi-isoform gene on the - strand."""
    return _write_gff3(tmp_path / "minus.gff3", [
        _gene_block("m1", "chr2", "-", [("m1.1", [(6000, 6300)])]),
        _gene_block("m2", "chr2", "-", [
            ("m2.1", [(7000, 7200), (7400, 7600)]),
            ("m2.2", [(7000, 7200), (7450, 7600)]),
        ]),
    ])


# ---------------------------------------------------------------------------
# Multi-isoform restriction — both strands
# ---------------------------------------------------------------------------


def test_multiiso_restriction_selects_right_genes_plus(tmp_path):
    genes = _parse_genes(_plus_mixed(tmp_path))
    multi = _multiiso_genes(genes, min_isoforms=2)
    assert sorted(g["gene_id"] for g in multi) == ["g2", "g4"]


def test_multiiso_restriction_minus_strand(tmp_path):
    genes = _parse_genes(_minus_mixed(tmp_path))
    multi = _multiiso_genes(genes, min_isoforms=2)
    assert [g["gene_id"] for g in multi] == ["m2"]
    assert all(g["strand"] == "-" for g in multi)


def test_multiiso_threshold_three(tmp_path):
    genes = _parse_genes(_plus_mixed(tmp_path))
    # Only g4 has >= 3 isoforms.
    assert [g["gene_id"] for g in _multiiso_genes(genes, min_isoforms=3)] == ["g4"]


# ---------------------------------------------------------------------------
# Subset GFF3 writer — round-trips only the multi-iso genes, coords preserved
# ---------------------------------------------------------------------------


def test_subset_gff3_roundtrips_plus(tmp_path):
    genes = _parse_genes(_plus_mixed(tmp_path))
    multi = _multiiso_genes(genes, min_isoforms=2)
    sub = tmp_path / "sub_plus.gff3"
    assert _write_genes_gff3(multi, sub) == 2
    reparsed = {g["gene_id"]: g for g in _parse_genes(sub)}
    assert set(reparsed) == {"g2", "g4"}
    # Concrete literal coordinates (internal 0-based half-open: 2000-1=1999 .. 2200).
    g2_exons = reparsed["g2"]["transcripts"][0]["exons"]
    assert (g2_exons[0].start, g2_exons[0].end) == (1999, 2200)
    assert reparsed["g2"]["strand"] == "+"


def test_subset_gff3_roundtrips_minus(tmp_path):
    genes = _parse_genes(_minus_mixed(tmp_path))
    multi = _multiiso_genes(genes, min_isoforms=2)
    sub = tmp_path / "sub_minus.gff3"
    assert _write_genes_gff3(multi, sub) == 1
    reparsed = {g["gene_id"]: g for g in _parse_genes(sub)}
    assert set(reparsed) == {"m2"}
    assert reparsed["m2"]["strand"] == "-"
    m2_exons = reparsed["m2"]["transcripts"][0]["exons"]
    assert (m2_exons[0].start, m2_exons[0].end) == (6999, 7200)


# ---------------------------------------------------------------------------
# Isoform-count summary on a known set
# ---------------------------------------------------------------------------


def test_isoform_count_summary_known_set():
    s = _summary([1, 1, 2, 3], min_isoforms=2)
    assert s["mean"] == 1.75
    assert s["median"] == 1.5
    assert s["multiiso_genes"] == 2
    assert s["total_genes"] == 4


def test_transcript_counts_from_gff(tmp_path):
    counts = _transcript_counts(_parse_genes(_plus_mixed(tmp_path)))
    assert sorted(counts) == [1, 1, 2, 3]


# ---------------------------------------------------------------------------
# isoform_accuracy — mocked mikado compare (real .stats format parsed)
# ---------------------------------------------------------------------------


def _patch_mikado(monkeypatch):
    """Make run_mikado_compare parse the captured real-format sample (no subprocess)."""
    levels = parse_mikado_compare_stats(DATA / "mikado_compare.stats")
    monkeypatch.setattr(iso, "run_mikado_compare",
                        lambda *a, **k: dict(levels))
    return levels


def test_intron_chain_and_transcript_f1_parsed(tmp_path, monkeypatch):
    _patch_mikado(monkeypatch)
    pred = _plus_mixed(tmp_path)
    ref = _write_gff3(tmp_path / "ref.gff3", [
        _gene_block("r1", "chr1", "+", [
            ("r1.1", [(2000, 2200), (2400, 2600)]),
            ("r1.2", [(2000, 2200), (2450, 2600)]),
        ]),
    ])
    res = isoform_accuracy(pred, ref, tmp_path / "iso", min_isoforms=2)
    # From the sample .stats: intron_chain f1 = 57.56, transcript (>=80% base) f1 = 62.26.
    assert res["whole"]["intron_chain"]["f1"] == 57.56
    assert res["whole"]["transcript"]["f1"] == 62.26
    assert res["multiiso"]["intron_chain"]["f1"] == 57.56
    assert res["multiiso"]["n_pred_genes"] == 2  # g2, g4
    assert res["multiiso"]["n_ref_genes"] == 1


def test_helixer_counts_default_one_per_gene(tmp_path, monkeypatch):
    _patch_mikado(monkeypatch)
    pred = _plus_mixed(tmp_path)
    ref = _minus_mixed(tmp_path)
    res = isoform_accuracy(pred, ref, tmp_path / "iso")
    c = res["counts"]
    # No helixer_gff3 → Helixer is 1 transcript/gene by design (CLAUDE.md §1).
    assert c["helixer_mean"] == 1.0 and c["helixer_median"] == 1.0
    # HelixForge output: counts [1,1,2,3] → mean 1.75; 2 genes gained isoforms.
    assert c["helixforge_mean"] == 1.75
    assert c["genes_gained_isoforms"] == 2


def test_helixer_counts_from_real_input(tmp_path, monkeypatch):
    _patch_mikado(monkeypatch)
    pred = _plus_mixed(tmp_path)
    ref = _minus_mixed(tmp_path)
    helixer = _write_gff3(tmp_path / "helixer.gff3", [
        _gene_block("h1", "chr1", "+", [("h1.1", [(1000, 1200)])]),
        _gene_block("h2", "chr1", "+", [("h2.1", [(2000, 2200)])]),
    ])
    res = isoform_accuracy(pred, ref, tmp_path / "iso", helixer_gff3=helixer)
    assert res["counts"]["helixer_mean"] == 1.0  # genuinely 1/gene here too


# ---------------------------------------------------------------------------
# OMArk attach / omit
# ---------------------------------------------------------------------------


def test_omark_attached_when_db_given(tmp_path, monkeypatch):
    _patch_mikado(monkeypatch)
    from helixforge.bench.wrappers import parse_omark_summary
    monkeypatch.setattr(iso, "run_omark",
                        lambda *a, **k: parse_omark_summary(DATA / "omark.sum"))
    res = isoform_accuracy(
        _plus_mixed(tmp_path), _minus_mixed(tmp_path), tmp_path / "iso",
        proteins_fa=tmp_path / "p.fa", omadb=tmp_path / "oma",
    )
    assert "omark" in res
    assert res["omark"]["fragmented"] == 3.0


def test_omark_omitted_without_db(tmp_path, monkeypatch):
    _patch_mikado(monkeypatch)
    res = isoform_accuracy(_plus_mixed(tmp_path), _minus_mixed(tmp_path), tmp_path / "iso")
    assert "omark" not in res


# ---------------------------------------------------------------------------
# Flatten → benchmark rows
# ---------------------------------------------------------------------------


def test_result_to_rows_metric_names(tmp_path, monkeypatch):
    _patch_mikado(monkeypatch)
    monkeypatch.setattr(iso, "run_omark",
                        lambda *a, **k: {"fragmented": 3.0, "consistent": 85.0})
    res = isoform_accuracy(
        _plus_mixed(tmp_path), _minus_mixed(tmp_path), tmp_path / "iso",
        proteins_fa=tmp_path / "p.fa", omadb=tmp_path / "oma",
    )
    rows = isoform_result_to_rows(res)
    metrics = {r["metric"]: r["value"] for r in rows}
    assert all(r["tool"] == "isoform" and r["status"] == "ok" for r in rows)
    assert metrics["isoform_intron_chain_f1"] == 57.56
    assert metrics["isoform_transcript_f1"] == 62.26
    assert metrics["isoform_intron_chain_f1_multiiso"] == 57.56
    assert metrics["isoform_count_genes_gained_isoforms"] == 2.0
    assert metrics["isoform_count_helixer_mean"] == 1.0
    assert metrics["isoform_omark_fragmented"] == 3.0


# ---------------------------------------------------------------------------
# Integration — real mikado (self-skips when absent)
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.skipif(shutil.which("mikado") is None, reason="mikado not installed")
def test_isoform_accuracy_end_to_end(tmp_path):
    pred = _plus_mixed(tmp_path)
    ref = _plus_mixed(tmp_path)  # identical → high self-comparison F1
    res = isoform_accuracy(pred, ref, tmp_path / "iso")
    assert "intron_chain" in res["whole"]
    assert res["multiiso"]["n_pred_genes"] == 2
