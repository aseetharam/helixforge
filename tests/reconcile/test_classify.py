"""Tests for reconcile/classify.py (Phase 3). Floor: 27.

CoverageCalculator is mocked (no real BAM); aggregated dicts are built from the
real StringTie aggregator or by hand.
"""

import pytest

from helixforge.io.stringtie import StringTieParser
from helixforge.reconcile import classify as classify_mod
from helixforge.reconcile.classify import classify_loci
from helixforge.reconcile.models import Exon, HelixerLocus, StringTieTranscript


# --- builders ---

def _locus(gene_id, seqid, start, end, strand, exons=None):
    return HelixerLocus(gene_id, seqid, start, end, strand, exons=exons or [])


def _st(tid, seqid, strand, exons, tpm, sample_id):
    start = min(e.start for e in exons)
    end = max(e.end for e in exons)
    return StringTieTranscript(
        tid, "STRG", seqid, start, end, strand, exons, tpm, sample_id
    )


# --- a fake CoverageCalculator factory ---

def _make_fake_cov(table_by_path):
    class _Fake:
        def __init__(self, path):
            self.path = path

        @classmethod
        def from_bam(cls, path):
            return cls(path)

        @classmethod
        def from_bigwig(cls, path):
            return cls(path)

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def mean_coverage(self, seqid, start, end):
            return table_by_path[self.path].get((seqid, start, end), 0.0)

    return _Fake


# --------------------------------------------------------------------------
# StringTie classification
# --------------------------------------------------------------------------

def _stringtie_setup():
    locusA = _locus("A", "chr1", 100, 300, "+", [Exon(100, 150), Exon(200, 300)])
    locusB = _locus("B", "chr1", 400, 500, "+", [Exon(400, 500)])
    locusC = _locus("C", "chr1", 600, 700, "-", [Exon(600, 700)])
    locusD = _locus("D", "chr1", 800, 900, "+", [Exon(800, 900)])
    txs = [
        _st("tA1", "chr1", "+", [Exon(100, 150), Exon(200, 300)], 5.0, "s1"),
        _st("tA2", "chr1", "+", [Exon(100, 150), Exon(200, 300)], 6.0, "s2"),
        _st("tB", "chr1", "+", [Exon(400, 500)], 0.3, "s1"),
        _st("tDanti", "chr1", "-", [Exon(800, 900)], 10.0, "s1"),
    ]
    agg = StringTieParser().aggregate_across_samples(txs, min_tpm=0.0, min_samples=1)
    return [locusA, locusB, locusC, locusD], txs, agg


def test_stringtie_expressed():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    a = res[0]
    assert a.status == "EXPRESSED"
    assert a.evidence_source == "stringtie"


def test_stringtie_expressed_max_tpm():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    assert res[0].max_tpm == pytest.approx(6.0)


def test_stringtie_expressed_num_samples():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    assert res[0].num_samples_expressed == 2


def test_stringtie_low_below_threshold():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    b = res[1]
    assert b.status == "LOW"
    assert b.evidence_source == "stringtie"
    assert b.max_tpm == pytest.approx(0.3)


def test_stringtie_no_overlap_silent_none():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    # locusC (chr1 -, 600-700) has no overlapping StringTie
    assert res[2].status == "SILENT"
    assert res[2].evidence_source == "none"


def test_stringtie_antisense_ignored():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    # locusD (+) overlaps only an antisense (-) transcript -> not classified
    assert res[3].status == "SILENT"
    assert res[3].evidence_source == "none"


def test_stringtie_order_preserved():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    assert [r.locus_id for r in res] == ["A", "B", "C", "D"]


def test_stringtie_min_samples_filter():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(
        loci, stringtie_aggregated=agg, stringtie_transcripts=txs, min_samples=3
    )
    # A has only 2 samples -> drops below threshold -> LOW
    assert res[0].status == "LOW"


def test_stringtie_min_tpm_relaxed():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(
        loci, stringtie_aggregated=agg, stringtie_transcripts=txs, min_tpm=0.1
    )
    # B at 0.3 now qualifies
    assert res[1].status == "EXPRESSED"


def test_stringtie_without_aggregated_uses_own_tpm():
    loci, txs, _ = _stringtie_setup()
    res = classify_loci(loci, stringtie_transcripts=txs)
    a = res[0]
    assert a.status == "EXPRESSED"
    assert a.max_tpm == pytest.approx(6.0)
    # no aggregation -> per-transcript samples = 1
    assert a.num_samples_expressed == 1


def test_stringtie_precedence_over_coverage(monkeypatch):
    loci, txs, agg = _stringtie_setup()
    fake = _make_fake_cov({"bam1": {}})
    monkeypatch.setattr(classify_mod, "CoverageCalculator", fake)
    res = classify_loci(
        loci, stringtie_aggregated=agg, stringtie_transcripts=txs, bam_paths=["bam1"]
    )
    # A is classified by StringTie, not coverage
    assert res[0].evidence_source == "stringtie"


# --------------------------------------------------------------------------
# Coverage fallback
# --------------------------------------------------------------------------

def test_coverage_expressed_weighted(monkeypatch):
    locus = _locus("A", "chr1", 100, 300, "+", [Exon(100, 150), Exon(200, 300)])
    table = {"bam1": {("chr1", 100, 150): 10.0, ("chr1", 200, 300): 4.0}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"])
    # weighted: (10*50 + 4*100)/150 = 6.0
    assert res[0].status == "EXPRESSED"
    assert res[0].evidence_source == "bam_coverage"
    assert res[0].mean_coverage == pytest.approx(6.0)


def test_coverage_weighted_prompt_example(monkeypatch):
    locus = _locus("E", "chr1", 0, 400, "+", [Exon(0, 200), Exon(300, 400)])
    table = {"bam1": {("chr1", 0, 200): 50.0, ("chr1", 300, 400): 10.0}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"])
    # (200*50 + 100*10)/300 = 36.666...
    assert res[0].mean_coverage == pytest.approx(36.6667, abs=1e-3)


def test_coverage_low(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {"bam1": {("chr1", 0, 100): 0.5}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"])
    assert res[0].status == "LOW"


def test_coverage_silent(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {"bam1": {("chr1", 0, 100): 0.0}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"])
    assert res[0].status == "SILENT"
    assert res[0].evidence_source == "bam_coverage"


def test_coverage_no_exons_uses_span(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", exons=[])
    table = {"bam1": {("chr1", 0, 100): 5.0}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"])
    assert res[0].mean_coverage == pytest.approx(5.0)


def test_coverage_multisample_mean(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {
        "b1": {("chr1", 0, 100): 4.0},
        "b2": {("chr1", 0, 100): 2.0},
    }
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["b1", "b2"])
    # mean of per-source exonic coverage: (4 + 2)/2 = 3.0
    assert res[0].mean_coverage == pytest.approx(3.0)


def test_coverage_threshold_custom(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {"bam1": {("chr1", 0, 100): 5.0}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"], coverage_threshold=10.0)
    # 5.0 < 10.0 threshold but >= near_zero -> LOW
    assert res[0].status == "LOW"


def test_bigwig_preferred_over_bam(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {
        "bam1": {("chr1", 0, 100): 0.0},
        "bw1": {("chr1", 0, 100): 9.0},
    }
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"], bigwig_paths=["bw1"])
    assert res[0].evidence_source == "bigwig"
    assert res[0].mean_coverage == pytest.approx(9.0)
    assert res[0].status == "EXPRESSED"


def test_coverage_fallback_when_no_stringtie_overlap(monkeypatch):
    # locusC has no StringTie overlap -> falls back to coverage
    loci, txs, agg = _stringtie_setup()
    table = {"bam1": {("chr1", 600, 700): 8.0}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci(
        loci, stringtie_aggregated=agg, stringtie_transcripts=txs, bam_paths=["bam1"]
    )
    assert res[2].status == "EXPRESSED"
    assert res[2].evidence_source == "bam_coverage"


# --------------------------------------------------------------------------
# No-data mode
# --------------------------------------------------------------------------

def test_no_data_all_silent():
    loci = [_locus("A", "chr1", 0, 100, "+"), _locus("B", "chr1", 200, 300, "-")]
    res = classify_loci(loci)
    assert all(r.status == "SILENT" for r in res)
    assert all(r.evidence_source == "none" for r in res)


def test_no_data_order_preserved():
    loci = [_locus("A", "chr1", 0, 100, "+"), _locus("B", "chr1", 200, 300, "-")]
    res = classify_loci(loci)
    assert [r.locus_id for r in res] == ["A", "B"]


def test_empty_loci_returns_empty():
    assert classify_loci([]) == []


def test_aggregated_only_without_transcripts_silent():
    # aggregated provided but no transcripts -> cannot determine overlap -> SILENT
    loci = [_locus("A", "chr1", 100, 300, "+", [Exon(100, 300)])]
    res = classify_loci(loci, stringtie_aggregated={"chr1:+:((100, 300),)": {
        "max_tpm": 9.0, "num_samples": 2}})
    assert res[0].status == "SILENT"
    assert res[0].evidence_source == "none"


# --------------------------------------------------------------------------
# Boundary + additional coverage cases
# --------------------------------------------------------------------------

def test_coverage_at_threshold_is_expressed(monkeypatch):
    # exactly coverage_threshold (>= comparison) -> EXPRESSED
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {"bam1": {("chr1", 0, 100): 2.0}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"])
    assert res[0].status == "EXPRESSED"


def test_coverage_at_near_zero_is_low(monkeypatch):
    # exactly near_zero_coverage (>= comparison) -> LOW
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {"bam1": {("chr1", 0, 100): 0.1}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"])
    assert res[0].status == "LOW"


def test_coverage_custom_near_zero(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {"bam1": {("chr1", 0, 100): 0.05}}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bam_paths=["bam1"], near_zero_coverage=0.01)
    # 0.05 >= custom near_zero 0.01 -> LOW (not SILENT)
    assert res[0].status == "LOW"


def test_bigwig_multisample_mean(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {
        "bw1": {("chr1", 0, 100): 6.0},
        "bw2": {("chr1", 0, 100): 2.0},
    }
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    res = classify_loci([locus], bigwig_paths=["bw1", "bw2"])
    assert res[0].evidence_source == "bigwig"
    assert res[0].mean_coverage == pytest.approx(4.0)


def test_stringtie_low_records_num_samples():
    loci, txs, agg = _stringtie_setup()
    res = classify_loci(loci, stringtie_aggregated=agg, stringtie_transcripts=txs)
    # locusB LOW but still records the supporting sample count (1)
    assert res[1].num_samples_expressed == 1


# --------------------------------------------------------------------------
# Phase 21 D3 — pooled coverage opens each source once (assesment §2.4)
# --------------------------------------------------------------------------

def _make_counting_cov(table_by_path, opens):
    """Like ``_make_fake_cov`` but records every ``from_bam``/``from_bigwig`` open."""
    class _Fake:
        def __init__(self, path):
            self.path = path

        @classmethod
        def from_bam(cls, path):
            opens.append(path)
            return cls(path)

        @classmethod
        def from_bigwig(cls, path):
            opens.append(path)
            return cls(path)

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def mean_coverage(self, seqid, start, end):
            return table_by_path[self.path].get((seqid, start, end), 0.0)

    return _Fake


def test_pool_opens_each_source_once_across_many_loci(monkeypatch):
    loci = [_locus(f"g{i}", "chr1", 0, 100, "+", [Exon(0, 100)]) for i in range(25)]
    table = {"b1": {("chr1", 0, 100): 5.0}, "b2": {("chr1", 0, 100): 3.0}}
    opens = []
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_counting_cov(table, opens))
    res = classify_loci(loci, bam_paths=["b1", "b2"])
    # 25 loci x 2 BAMs == 50 opens under the old per-locus path; pooled => 2.
    assert opens == ["b1", "b2"]
    assert all(r.mean_coverage == pytest.approx(4.0) for r in res)  # (5+3)/2


def test_pooled_coverage_equals_per_locus_reference(monkeypatch):
    # Mixed strands; coverage is strand-agnostic so the pooled result must equal
    # a hand-computed per-locus weighted mean (the old per-locus-open behavior).
    loci = [
        _locus("A", "chr1", 100, 300, "+", [Exon(100, 150), Exon(200, 300)]),
        _locus("B", "chr1", 0, 100, "-", [Exon(0, 100)]),
    ]
    table = {"b1": {
        ("chr1", 100, 150): 10.0, ("chr1", 200, 300): 4.0, ("chr1", 0, 100): 1.5,
    }}
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_fake_cov(table))
    pooled = classify_loci(loci, bam_paths=["b1"])

    ref = []
    for locus in loci:
        regions = [(e.start, e.end) for e in locus.exons]
        total = sum(e - s for s, e in regions)
        w = sum(table["b1"][(locus.seqid, s, e)] * (e - s) for s, e in regions) / total
        ref.append(w)
    assert [r.mean_coverage for r in pooled] == pytest.approx(ref)


def test_pool_bigwig_preferred_never_opens_bam(monkeypatch):
    locus = _locus("A", "chr1", 0, 100, "+", [Exon(0, 100)])
    table = {"bw1": {("chr1", 0, 100): 7.0}}
    opens = []
    monkeypatch.setattr(classify_mod, "CoverageCalculator", _make_counting_cov(table, opens))
    res = classify_loci([locus], bam_paths=["b1"], bigwig_paths=["bw1"])
    # bigWig is preferred, so only it is opened — the BAM handle is never created.
    assert opens == ["bw1"]
    assert res[0].evidence_source == "bigwig"
    assert res[0].mean_coverage == pytest.approx(7.0)
