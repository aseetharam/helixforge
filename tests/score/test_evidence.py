"""Tests for the standalone RNA-seq ``evidence`` scorer (Phase 15). Floor: 20.

Junction union/merge from BAM + STAR SJ, intron precision/recall/F1 on concrete
supported / contradicted / novel cases (both strands), single-exon handling, TPM
linkage by structure hash, optional coverage, and the region filter. Coordinates
are literal internal 0-based half-open (CLAUDE.md §12).
"""

import pytest

from helixforge.io.bam import CoverageCalculator
from helixforge.reconcile.models import Exon, SpliceJunction
from helixforge.score.evidence import (
    TSV_COLUMNS,
    _transcript_coverage_metrics,
    collect_coverage,
    collect_junctions,
    rna_aed_from_ratios,
    score_annotation,
    score_transcript_evidence,
    summarize_evidence,
    write_evidence_tsv,
)

# Literal structures matching the GFF3 fixture.
SUPPORTED_PLUS = [Exon(100, 150), Exon(250, 300)]   # intron (150, 250)
CONTRADICTED_PLUS = [Exon(100, 160), Exon(250, 300)]  # intron (160, 250)
NOVEL_PLUS = [Exon(600, 650), Exon(760, 800)]       # intron (650, 760)
SINGLE_EXON = [Exon(100, 150)]
SUPPORTED_MINUS = [Exon(100, 150), Exon(250, 300)]  # chr2 -, intron (150, 250)


def _tx(seqid, strand, exons, tid="t"):
    return {
        "transcript_id": tid, "seqid": seqid, "strand": strand, "exons": exons,
        "start": min(e.start for e in exons), "end": max(e.end for e in exons),
    }


def _plus_junc(read_count=5):
    return [SpliceJunction("chr1", 150, 250, "+", read_count=read_count)]


# ---------------------------------------------------------------------------
# collect_junctions
# ---------------------------------------------------------------------------

def test_collect_requires_a_source():
    with pytest.raises(ValueError):
        collect_junctions()


def test_collect_bam_only(ev_bam_path):
    juncs = collect_junctions(bam_paths=[ev_bam_path])
    keys = {(j.seqid, j.donor, j.acceptor, j.strand): j for j in juncs}
    assert (("chr1", 150, 250, "+")) in keys
    assert (("chr2", 150, 250, "-")) in keys
    assert keys[("chr1", 150, 250, "+")].read_count == 5
    assert keys[("chr2", 150, 250, "-")].read_count == 3


def test_collect_star_only(ev_star_sj_path):
    juncs = collect_junctions(star_sj_paths=[ev_star_sj_path])
    keys = {(j.donor, j.acceptor) for j in juncs}
    # (500, 600) had only 2 reads -> below default min_reads=3, filtered
    assert (150, 250) in keys
    assert (350, 450) in keys
    assert (500, 600) not in keys


def test_collect_star_min_reads_relaxed(ev_star_sj_path):
    juncs = collect_junctions(star_sj_paths=[ev_star_sj_path], min_reads=1)
    assert (500, 600) in {(j.donor, j.acceptor) for j in juncs}


def test_collect_union_merges_bam_and_sj(ev_bam_path, ev_star_sj_path):
    juncs = collect_junctions(bam_paths=[ev_bam_path], star_sj_paths=[ev_star_sj_path])
    merged = {(j.seqid, j.donor, j.acceptor, j.strand): j for j in juncs}
    j = merged[("chr1", 150, 250, "+")]
    assert j.read_count == 13          # 5 (BAM) + 8 (STAR)
    assert j.samples == 2              # two distinct source files
    assert merged[("chr1", 350, 450, "+")].read_count == 6


def test_collect_multisample_bam_merge(ev_bam_path, ev_bam_path2):
    juncs = collect_junctions(bam_paths=[ev_bam_path, ev_bam_path2])
    merged = {(j.donor, j.acceptor, j.strand): j for j in juncs}
    j = merged[(150, 250, "+")]
    assert j.read_count == 7           # 5 + 2
    assert j.samples == 2


def test_collect_region_filter(ev_bam_path, ev_star_sj_path):
    juncs = collect_junctions(
        bam_paths=[ev_bam_path], star_sj_paths=[ev_star_sj_path], region="chr2"
    )
    assert all(j.seqid == "chr2" for j in juncs)
    assert len(juncs) == 1


def test_collect_region_span_filter(ev_star_sj_path):
    juncs = collect_junctions(star_sj_paths=[ev_star_sj_path], region="chr1:340-460")
    assert {(j.donor, j.acceptor) for j in juncs} == {(350, 450)}


def test_collect_sorted(ev_bam_path, ev_star_sj_path):
    juncs = collect_junctions(bam_paths=[ev_bam_path], star_sj_paths=[ev_star_sj_path])
    keys = [(j.seqid, j.donor, j.acceptor) for j in juncs]
    assert keys == sorted(keys)


# ---------------------------------------------------------------------------
# score_transcript_evidence — concordance classes, both strands
# ---------------------------------------------------------------------------

def test_transcript_supported_plus():
    m = score_transcript_evidence(_tx("chr1", "+", SUPPORTED_PLUS), _plus_junc())
    assert m["num_introns"] == 1
    assert m["supported"] == 1
    assert m["contradicted"] == 0
    assert m["novel_in_data"] == 0
    assert m["junction_support_fraction"] == 1.0
    assert m["intron_precision"] == 1.0
    assert m["intron_recall"] == 1.0
    assert m["intron_f1"] == 1.0


def test_transcript_contradicted_plus():
    m = score_transcript_evidence(_tx("chr1", "+", CONTRADICTED_PLUS), _plus_junc())
    assert m["supported"] == 0
    assert m["contradicted"] == 1
    assert m["novel_in_data"] == 0
    assert m["junction_support_fraction"] == 0.0
    assert m["intron_precision"] == 0.0


def test_transcript_novel_plus():
    m = score_transcript_evidence(_tx("chr1", "+", NOVEL_PLUS), _plus_junc())
    assert m["supported"] == 0
    assert m["contradicted"] == 0
    assert m["novel_in_data"] == 1
    assert m["junction_support_fraction"] == 0.0
    # no qualifying junction in span -> recall (and hence f1) undefined
    assert m["intron_recall"] is None
    assert m["intron_f1"] is None


def test_transcript_supported_minus():
    juncs = [SpliceJunction("chr2", 150, 250, "-", read_count=3)]
    m = score_transcript_evidence(_tx("chr2", "-", SUPPORTED_MINUS), juncs)
    assert m["supported"] == 1
    assert m["junction_support_fraction"] == 1.0
    assert m["intron_f1"] == 1.0


def test_transcript_single_exon_no_introns():
    m = score_transcript_evidence(_tx("chr1", "+", SINGLE_EXON), _plus_junc())
    assert m["num_introns"] == 0
    assert m["supported"] == 0
    assert m["junction_support_fraction"] is None
    assert m["intron_precision"] is None
    assert m["intron_recall"] is None
    assert m["intron_f1"] is None


def test_transcript_min_reads_gating():
    weak = [SpliceJunction("chr1", 150, 250, "+", read_count=2)]
    novel = score_transcript_evidence(_tx("chr1", "+", SUPPORTED_PLUS), weak, min_reads=3)
    assert novel["supported"] == 0 and novel["novel_in_data"] == 1
    ok = score_transcript_evidence(_tx("chr1", "+", SUPPORTED_PLUS), weak, min_reads=1)
    assert ok["supported"] == 1


def test_transcript_coverage_passthrough():
    m = score_transcript_evidence(_tx("chr1", "+", SUPPORTED_PLUS), _plus_junc(), coverage=7.0)
    assert m["mean_coverage"] == 7.0


def test_transcript_no_coverage_key_when_absent():
    m = score_transcript_evidence(_tx("chr1", "+", SUPPORTED_PLUS), _plus_junc())
    assert "mean_coverage" not in m


def test_transcript_tpm_linkage_by_overlap_plus():
    # A StringTie transcript whose exons overlap the model (here exact, but
    # overlap — not structural identity — is what matters) → its TPM is assigned.
    idx = {"chr1": [("+", 100, 300, ((100, 150), (250, 300)), 9.0)]}
    m = score_transcript_evidence(
        _tx("chr1", "+", SUPPORTED_PLUS), _plus_junc(), tpm_index=idx
    )
    assert m["tpm"] == 9.0


def test_transcript_tpm_linkage_inexact_overlap_plus():
    # Different exon boundaries than the model, but still overlapping exonically:
    # the old exact-structure-hash match missed this (the poster bug); overlap
    # matching assigns the TPM.
    idx = {"chr1": [("+", 100, 320, ((100, 160), (250, 320)), 7.5)]}
    m = score_transcript_evidence(
        _tx("chr1", "+", SUPPORTED_PLUS), _plus_junc(), tpm_index=idx
    )
    assert m["tpm"] == 7.5


def test_transcript_tpm_linkage_by_overlap_minus():
    idx = {"chr2": [("-", 100, 300, ((100, 150), (250, 300)), 4.0)]}
    m = score_transcript_evidence(
        _tx("chr2", "-", SUPPORTED_MINUS), _plus_junc(), tpm_index=idx
    )
    assert m["tpm"] == 4.0


def test_transcript_tpm_absent_when_no_exonic_overlap():
    # StringTie transcript is on chr1 + but lies entirely in the model's intron
    # (no exonic base shared) → no TPM.
    idx = {"chr1": [("+", 160, 240, ((160, 240),), 5.0)]}
    m = score_transcript_evidence(
        _tx("chr1", "+", SUPPORTED_PLUS), _plus_junc(), tpm_index=idx
    )
    assert "tpm" not in m


def test_transcript_tpm_absent_on_strand_mismatch():
    # Same coordinates but opposite strand → not a match.
    idx = {"chr1": [("-", 100, 300, ((100, 150), (250, 300)), 9.0)]}
    m = score_transcript_evidence(
        _tx("chr1", "+", SUPPORTED_PLUS), _plus_junc(), tpm_index=idx
    )
    assert "tpm" not in m


def test_transcript_tpm_best_overlap_wins():
    # Two candidates: the one with greater exonic overlap is chosen (not the
    # higher TPM with less overlap).
    idx = {"chr1": [
        ("+", 100, 300, ((100, 150), (250, 300)), 2.0),  # full overlap (100 bp)
        ("+", 100, 130, ((100, 130),), 99.0),            # partial overlap (30 bp)
    ]}
    m = score_transcript_evidence(
        _tx("chr1", "+", SUPPORTED_PLUS), _plus_junc(), tpm_index=idx
    )
    assert m["tpm"] == 2.0


def test_transcript_no_exons_raises():
    with pytest.raises(ValueError):
        score_transcript_evidence(
            {"transcript_id": "x", "seqid": "chr1", "strand": "+", "exons": []},
            _plus_junc(),
        )


# ---------------------------------------------------------------------------
# score_annotation — end to end on the synthetic GFF3
# ---------------------------------------------------------------------------

@pytest.fixture
def scored(ev_gff3_path, ev_bam_path, ev_star_sj_path, ev_stringtie_gtfs):
    df = score_annotation(
        ev_gff3_path,
        bam_paths=[ev_bam_path], star_sj_paths=[ev_star_sj_path],
        stringtie_gtfs=ev_stringtie_gtfs,
    )
    return {row["transcript_id"]: row for _, row in df.iterrows()}, df


def test_score_annotation_columns(scored):
    _, df = scored
    assert list(df.columns) == TSV_COLUMNS


def test_score_annotation_row_count(scored):
    rows, _ = scored
    assert set(rows) == {"gA.t1", "gA.t2", "gB.t1", "gS.t1", "gM.t1"}


def test_score_annotation_supported_row(scored):
    rows, _ = scored
    r = rows["gA.t1"]
    assert r["supported"] == 1
    assert r["junction_support_fraction"] == 1.0
    assert r["intron_f1"] == 1.0


def test_score_annotation_contradicted_row(scored):
    rows, _ = scored
    assert rows["gA.t2"]["contradicted"] == 1
    assert rows["gA.t2"]["supported"] == 0


def test_score_annotation_novel_row(scored):
    rows, _ = scored
    assert rows["gB.t1"]["novel_in_data"] == 1


def test_score_annotation_minus_supported(scored):
    rows, _ = scored
    assert rows["gM.t1"]["strand"] == "-"
    assert rows["gM.t1"]["supported"] == 1


def test_score_annotation_single_exon(scored):
    rows, _ = scored
    import pandas as pd
    r = rows["gS.t1"]
    assert r["num_introns"] == 0
    assert pd.isna(r["junction_support_fraction"])


def test_score_annotation_tpm_linked(scored):
    import pandas as pd
    rows, _ = scored
    # Overlap-based (not exact-structure): the aggregated max TPM (9.0 across
    # samples) is assigned to every model overlapping the StringTie assembly on
    # the same strand, even with different exon boundaries (the poster bug fix).
    assert rows["gA.t1"]["tpm"] == 9.0          # exact structure
    assert rows["gA.t2"]["tpm"] == 9.0          # different boundary, still overlaps
    assert rows["gS.t1"]["tpm"] == 9.0          # single-exon overlap of the assembly
    assert pd.isna(rows["gB.t1"]["tpm"])        # no overlapping StringTie transcript
    assert pd.isna(rows["gM.t1"]["tpm"])        # chr2 minus — no StringTie there


def test_score_annotation_coverage(scored):
    rows, _ = scored
    assert rows["gA.t1"]["mean_coverage"] == 5.0    # 5 reads cover both exons
    assert rows["gM.t1"]["mean_coverage"] == 3.0
    assert rows["gB.t1"]["mean_coverage"] == 0.0


def test_score_annotation_coverage_none_without_bam(ev_gff3_path, ev_star_sj_path):
    df = score_annotation(ev_gff3_path, star_sj_paths=[ev_star_sj_path])
    assert df["mean_coverage"].isna().all()


def test_score_annotation_tpm_none_without_stringtie(ev_gff3_path, ev_bam_path):
    df = score_annotation(ev_gff3_path, bam_paths=[ev_bam_path])
    assert df["tpm"].isna().all()


def test_score_annotation_region_filter(ev_gff3_path, ev_bam_path):
    df = score_annotation(ev_gff3_path, bam_paths=[ev_bam_path], region="chr2")
    assert set(df["transcript_id"]) == {"gM.t1"}


def test_score_annotation_region_span_filter(ev_gff3_path, ev_bam_path):
    df = score_annotation(ev_gff3_path, bam_paths=[ev_bam_path], region="chr1:601-800")
    assert set(df["transcript_id"]) == {"gB.t1"}


# ---------------------------------------------------------------------------
# boundary_ratio — supported when the *terminal exon is expressed* (per-base
# median >= threshold), NOT when the per-exon minimum reaches it. Regression for
# the poster bug where a deeply covered multi-exon gene got boundary_ratio=0: a
# single low base anywhere in a terminal exon (a coverage dip, or a Helixer UTR
# annotated past where reads reach) zeroed the minimum and pinned rna_aed at 0.2.
# Cache value = (mean, median, min).
# ---------------------------------------------------------------------------

# A deeply expressed 3-exon gene whose terminal exons have a 0-depth base inside
# (min=0) but a high median — exactly the #000007 pathology.
_DIP_EXONS = [Exon(0, 100), Exon(200, 300), Exon(400, 500)]
_DIP_CACHE = {
    ("chr1", 0, 100): (50.0, 131.0, 0.0),    # median 131, min 0 (interior/UTR-tip dip)
    ("chr1", 200, 300): (50.0, 50.0, 40.0),  # interior exon
    ("chr1", 400, 500): (50.0, 131.0, 0.0),  # median 131, min 0
}


def test_boundary_ratio_uses_median_not_exon_min_plus():
    # Old per-exon-min logic: first/last min == 0 -> boundary 0. Median logic:
    # both terminal exons are expressed (median 131) -> boundary 1.0.
    _, cov_ratio, boundary = _transcript_coverage_metrics("chr1", "+", _DIP_EXONS, _DIP_CACHE, 5)
    assert boundary == 1.0
    assert cov_ratio == 1.0          # all three medians >= 5


def test_boundary_ratio_uses_median_not_exon_min_minus():
    # Same cache, minus strand: ratio is symmetric -> still 1.0.
    _, _, boundary = _transcript_coverage_metrics("chr1", "-", _DIP_EXONS, _DIP_CACHE, 5)
    assert boundary == 1.0


def test_well_supported_multiexon_reaches_zero_aed():
    # #000007-style: junctions all supported (jr=1), all exons covered (cr=1),
    # both termini expressed (br=1) -> rna_aed 0, not the old 0.2 (br pinned to 0).
    _, cov_ratio, boundary = _transcript_coverage_metrics("chr1", "+", _DIP_EXONS, _DIP_CACHE, 5)
    assert rna_aed_from_ratios(1.0, cov_ratio, boundary) == 0.0
    # The bug value for contrast: a perfectly junction/coverage-supported gene
    # with boundary wrongly 0 was pinned here.
    assert rna_aed_from_ratios(1.0, 1.0, 0.0) == pytest.approx(0.2)


def test_boundary_ratio_half_when_one_terminal_exon_unexpressed():
    # The high-coordinate terminal exon is not expressed (median 1 < 5) -> 0.5.
    cache = dict(_DIP_CACHE)
    cache[("chr1", 400, 500)] = (1.0, 1.0, 0.0)
    _, _, boundary_plus = _transcript_coverage_metrics("chr1", "+", _DIP_EXONS, cache, 5)
    assert boundary_plus == 0.5
    # On the minus strand the same unexpressed high-coordinate exon is the *start*;
    # the ratio value is identical (symmetric) but the uncovered terminus is the 5' one.
    _, _, boundary_minus = _transcript_coverage_metrics("chr1", "-", _DIP_EXONS, cache, 5)
    assert boundary_minus == 0.5


def test_boundary_ratio_zero_when_both_terminal_exons_unexpressed():
    cache = {
        ("chr1", 0, 100): (2.0, 2.0, 0.0),       # median 2 < 5
        ("chr1", 200, 300): (50.0, 50.0, 40.0),  # interior covered
        ("chr1", 400, 500): (1.0, 1.0, 0.0),     # median 1 < 5
    }
    _, _, boundary = _transcript_coverage_metrics("chr1", "+", _DIP_EXONS, cache, 5)
    assert boundary == 0.0


# ---------------------------------------------------------------------------
# summarize_evidence + write_evidence_tsv
# ---------------------------------------------------------------------------

def test_summarize_evidence(scored):
    _, df = scored
    s = summarize_evidence(df)
    assert s["n_transcripts"] == 5
    assert s["n_multi_exon"] == 4
    assert s["n_single_exon"] == 1
    assert s["n_fully_supported"] == 2          # gA.t1, gM.t1
    assert s["fraction_fully_supported"] == 0.5
    assert s["mean_junction_support_fraction"] == 0.5


def test_summarize_evidence_empty():
    import pandas as pd
    s = summarize_evidence(pd.DataFrame(columns=TSV_COLUMNS))
    assert s["n_transcripts"] == 0
    assert s["fraction_fully_supported"] is None


def test_summarize_evidence_reports_aed_denominators():
    # RNA AED on all 3 transcripts; protein AED on only 1 -> the means are over
    # different denominators and the summary must surface each count so they are
    # never silently compared on unlike sets.
    import pandas as pd
    df = pd.DataFrame(
        [
            {"transcript_id": "t1", "num_introns": 1, "rna_aed": 0.2, "protein_aed": 0.1},
            {"transcript_id": "t2", "num_introns": 1, "rna_aed": 0.4, "protein_aed": None},
            {"transcript_id": "t3", "num_introns": 0, "rna_aed": 0.6, "protein_aed": None},
        ],
        columns=TSV_COLUMNS,
    )
    s = summarize_evidence(df)
    assert s["n_rna_aed"] == 3
    assert s["n_protein_aed"] == 1
    assert s["mean_rna_aed"] == pytest.approx((0.2 + 0.4 + 0.6) / 3)
    assert s["mean_protein_aed"] == pytest.approx(0.1)  # over the 1 hit only


def test_write_evidence_tsv(scored, tmp_path):
    _, df = scored
    out = write_evidence_tsv(df, tmp_path / "evidence.tsv")
    text = out.read_text()
    assert "\t".join(TSV_COLUMNS) in text.splitlines()[0]
    assert "gA.t1" in text


# ---------------------------------------------------------------------------
# Extract-once + parallelism (the performance fix): each BAM is opened once and
# scanned once per gene locus — never per transcript — and ``--threads`` engages
# real cores via process pools while keeping output identical to ``--threads 1``.
# ---------------------------------------------------------------------------

def test_collect_coverage_opens_each_bam_once(ev_bam_path, ev_bam_path2, monkeypatch):
    # Many exon regions, two BAMs: ``from_bam`` is called exactly once per BAM,
    # never once per region/transcript (the old per-transcript-BAM-access bug).
    opens = []
    orig = CoverageCalculator.from_bam.__func__

    def spy(cls, path, **kw):
        opens.append(str(path))
        return orig(cls, path, **kw)

    monkeypatch.setattr(CoverageCalculator, "from_bam", classmethod(spy))
    regions = [("chr1", 100, 150), ("chr1", 250, 300), ("chr1", 120, 200)]
    collect_coverage([ev_bam_path, ev_bam_path2], regions, threads=1)
    assert sorted(opens) == sorted([ev_bam_path, ev_bam_path2])  # one open per BAM


def test_collect_coverage_queries_scale_with_loci_not_transcripts(
    ev_bam_path, monkeypatch
):
    # Five exon regions collapse to two merged spans ([100,300) and [400,450)).
    # The number of pileup queries is O(merged spans) — independent of how many
    # transcripts share those exons — so it stays 2, not 5.
    calls = []
    orig = CoverageCalculator.region_coverage_array

    def spy(self, seqid, start, end, **kw):
        calls.append((seqid, start, end))
        return orig(self, seqid, start, end, **kw)

    monkeypatch.setattr(CoverageCalculator, "region_coverage_array", spy)
    regions = [
        ("chr1", 100, 150),
        ("chr1", 150, 250),
        ("chr1", 250, 300),
        ("chr1", 400, 420),
        ("chr1", 420, 450),
    ]
    collect_coverage([ev_bam_path], regions, threads=1)
    assert calls == [("chr1", 100, 300), ("chr1", 400, 450)]


def _two_bam_scored(gff3, bams, threads):
    df = score_annotation(gff3, bam_paths=bams, threads=threads)
    return df.reset_index(drop=True)


def test_score_annotation_threads_match_serial(
    ev_gff3_path, ev_bam_path, ev_bam_path2, ev_star_sj_path
):
    serial = _two_bam_scored(ev_gff3_path, [ev_bam_path, ev_bam_path2], 1)
    parallel = _two_bam_scored(ev_gff3_path, [ev_bam_path, ev_bam_path2], 4)
    # Byte-identical output regardless of thread count (rows pre-sorted).
    assert serial.equals(parallel)


def test_score_annotation_threads_use_process_pool(
    ev_gff3_path, ev_bam_path, ev_bam_path2, monkeypatch
):
    # With ``threads > 1`` the scoring (and multi-BAM extraction) runs in a
    # ProcessPoolExecutor — Python CPU work the GIL would otherwise serialize.
    import helixforge.score.evidence as ev

    seen = []
    real_pool = ev.ProcessPoolExecutor

    def recording_pool(*args, **kwargs):
        seen.append(kwargs.get("max_workers"))
        return real_pool(*args, **kwargs)

    monkeypatch.setattr(ev, "ProcessPoolExecutor", recording_pool)
    score_annotation(ev_gff3_path, bam_paths=[ev_bam_path, ev_bam_path2], threads=2)
    # at least one process pool was created with the requested worker count
    assert 2 in seen
