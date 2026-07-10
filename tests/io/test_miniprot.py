"""Tests for MiniprotParser (Phase 2). Floor: 20. Both strands covered."""

import pytest

from helixforge.io.miniprot import MiniprotParser


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        MiniprotParser("/no/such.gff3")


def test_parse_count(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse()
    assert len(alns) == 4


def test_parse_sorted_by_seqid_start(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse()
    assert [a.protein_id for a in alns] == ["protA", "protC", "protA", "protB"]


def test_mrna_coordinate_conversion(miniprot_gff_path):
    # MP1 mRNA 101-300 -> internal (100, 300)
    a = MiniprotParser(miniprot_gff_path).parse()[0]
    assert a.start == 100
    assert a.end == 300


def test_cds_segments_converted_and_sorted(miniprot_gff_path):
    a = MiniprotParser(miniprot_gff_path).parse()[0]
    assert [(c.start, c.end) for c in a.cds_segments] == [(100, 150), (200, 300)]


def test_cds_phase(miniprot_gff_path):
    a = MiniprotParser(miniprot_gff_path).parse()[0]
    assert all(c.phase == 0 for c in a.cds_segments)


def test_identity_fraction_already_normalized(miniprot_gff_path):
    a = MiniprotParser(miniprot_gff_path).parse()[0]
    assert a.identity == pytest.approx(0.95)


def test_identity_percent_normalized(miniprot_gff_path):
    # MP2 has Identity=88.0% -> 0.88
    alns = MiniprotParser(miniprot_gff_path).parse()
    mp2 = next(a for a in alns if a.protein_id == "protB")
    assert mp2.identity == pytest.approx(0.88)


def test_query_coverage_full(miniprot_gff_path):
    # MP1 covers protA 1-100, denom 100 -> coverage 1.0
    a = MiniprotParser(miniprot_gff_path).parse()[0]
    assert a.query_coverage == pytest.approx(1.0)


def test_query_coverage_partial(miniprot_gff_path):
    # MP3 covers protA 40-100, denom 100 -> 61/100 = 0.61
    alns = MiniprotParser(miniprot_gff_path).parse()
    mp3 = next(a for a in alns if a.protein_id == "protA" and a.start == 400)
    assert mp3.query_coverage == pytest.approx(0.61)


def test_rank_parsed(miniprot_gff_path):
    a = MiniprotParser(miniprot_gff_path).parse()[0]
    assert a.rank == 1


def test_score_from_column6(miniprot_gff_path):
    a = MiniprotParser(miniprot_gff_path).parse()[0]
    assert a.score == 200.0


def test_minus_strand_alignment(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse()
    mp2 = next(a for a in alns if a.protein_id == "protB")
    assert mp2.strand == "-"
    assert (mp2.start, mp2.end) == (100, 250)


def test_filter_min_identity(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse(min_identity=0.8)
    # MP3 (0.60) dropped; MP1(0.95), MP4(0.90), MP2(0.88) remain
    assert len(alns) == 3
    assert all(a.identity >= 0.8 for a in alns)


def test_filter_min_coverage(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse(min_coverage=0.7)
    # MP3 (0.61) dropped
    assert all(a.query_coverage >= 0.7 for a in alns)
    assert len(alns) == 3


def test_parse_for_region(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse_for_region("chr1", 100, 300)
    pids = sorted(a.protein_id for a in alns)
    # MP1 (protA) and MP4 (protC) overlap [100,300) on chr1
    assert pids == ["protA", "protC"]


def test_parse_for_region_excludes_non_overlap(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse_for_region("chr1", 100, 110)
    # only MP1 overlaps [100,110); MP4 starts at 119
    assert [a.protein_id for a in alns] == ["protA"]


def test_get_best_per_locus_sorted_by_rank_then_score(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse()
    best = MiniprotParser.get_best_per_locus(alns, "chr1", 100, 300, "+")
    # MP1 (rank1) before MP4 (rank2); MP3 doesn't overlap; MP2 wrong strand/seqid
    assert [a.rank for a in best] == [1, 2]
    assert best[0].score == 200.0


def test_get_best_per_locus_excludes_wrong_strand(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse()
    best = MiniprotParser.get_best_per_locus(alns, "chr1", 100, 300, "-")
    assert best == []


def test_get_best_per_locus_minus_strand(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse()
    best = MiniprotParser.get_best_per_locus(alns, "chr2", 100, 250, "-")
    assert [a.protein_id for a in best] == ["protB"]


def test_all_coverage_in_unit_range(miniprot_gff_path):
    alns = MiniprotParser(miniprot_gff_path).parse()
    assert all(0.0 <= a.query_coverage <= 1.0 for a in alns)
    assert all(0.0 <= a.identity <= 1.0 for a in alns)


# --------------------------------------------------------------------------
# Phase 19 D3: gffutils DB reuse via dbfn / keep_db. Default (None) preserves
# the historical in-memory behavior; a persistent dbfn parses identically.
# --------------------------------------------------------------------------

def _aln_tuples(alns):
    return [
        (a.protein_id, a.seqid, a.start, a.end, a.strand,
         [(c.start, c.end, c.phase) for c in a.cds_segments])
        for a in alns
    ]


def test_persistent_dbfn_matches_memory(miniprot_gff_path, tmp_path):
    in_memory = MiniprotParser(miniprot_gff_path).parse()
    dbfn = tmp_path / "miniprot.gffutils.db"
    on_disk = MiniprotParser(miniprot_gff_path, dbfn=str(dbfn), keep_db=True).parse()
    assert dbfn.exists()
    assert _aln_tuples(on_disk) == _aln_tuples(in_memory)
    # reuse path: second construction opens the on-disk DB without rebuilding
    mtime = dbfn.stat().st_mtime_ns
    reused = MiniprotParser(miniprot_gff_path, dbfn=str(dbfn), keep_db=True).parse()
    assert dbfn.stat().st_mtime_ns == mtime
    assert _aln_tuples(reused) == _aln_tuples(in_memory)
