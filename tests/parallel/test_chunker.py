"""Tests for parallel/chunker.py — v1 strategies, gene-respecting (realign).

Floor: 8. Concrete literal coordinates only; both strands present via the
``synthetic_genome`` fixture (chr1 plus-strand loci, chr2 a minus-strand locus).
Every cut must land in an inter-locus gap so a gene is never split.
"""

from helixforge.parallel.chunker import ChunkStrategy, plan_chunks
from helixforge.parallel.plan import _region_bounds
from helixforge.reconcile.locus import load_helixer_loci


def _loci(synthetic_genome):
    return load_helixer_loci(synthetic_genome["gff3"])


def _lengths(synthetic_genome):
    return dict(synthetic_genome["scaffold_len"])


def test_scaffold_strategy_one_chunk_per_scaffold(synthetic_genome):
    specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                        strategy="scaffold", min_boundary_gap=1000)
    # chr1 (5 loci) + chr2 (1 minus-strand locus) → exactly two whole-scaffold chunks.
    assert len(specs) == 2
    assert specs[0].region == "chr1:0-20000"
    assert specs[1].region == "chr2:0-5000"
    assert specs[0].num_loci == 5
    assert specs[1].locus_ids == ["G6"]  # minus-strand locus on its own scaffold


def test_genes_strategy_two_loci_per_chunk(synthetic_genome):
    specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                        strategy="genes", chunk_size=2, min_boundary_gap=1000)
    # chr1 cuts in the 5000 bp gap after G2 and the 6000 bp gap after G4.
    assert [s.region for s in specs] == [
        "chr1:0-3500", "chr1:3500-10000", "chr1:10000-20000", "chr2:0-5000",
    ]
    assert [s.locus_ids for s in specs] == [
        ["G1", "G2"], ["G3", "G4"], ["G5"], ["G6"],
    ]


def test_genes_strategy_cut_falls_in_gap_not_in_a_locus(synthetic_genome):
    specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                        strategy="genes", chunk_size=2, min_boundary_gap=1000)
    # The cut between chunk 0 and 1 (3500) sits inside the G2→G3 gap [1000, 6000).
    _seqid, _lo, hi = _region_bounds(specs[0].region)
    assert hi == 3500
    assert 1000 < hi < 6000


def test_size_strategy_snaps_to_gap(synthetic_genome):
    # 5 kb windows: span exceeds 5000 after G3 (0..6500) → cut in the 6000 bp gap
    # after G4? No — first eligible gap once span>=5000 is the G4→G5 gap [7000,13000).
    specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                        strategy="size", chunk_size=5000, min_chunk_size=1000,
                        min_boundary_gap=1000)
    # chr1 → one cut at the G4→G5 gap midpoint (10000); chr2 whole.
    assert [s.region for s in specs] == [
        "chr1:0-10000", "chr1:10000-20000", "chr2:0-5000",
    ]


def test_adaptive_strategy_targets_chunk_count(synthetic_genome):
    specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                        strategy="adaptive", target_chunks=3, min_boundary_gap=1000)
    # 6 loci / 3 → 2 loci per chunk; same cuts as the genes(2) case.
    assert [s.num_loci for s in specs] == [2, 2, 1, 1]
    assert specs[3].locus_ids == ["G6"]


def test_max_chunk_size_splits_long_scaffold(synthetic_genome):
    specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                        strategy="scaffold", max_chunk_size=5000,
                        min_chunk_size=1000, min_boundary_gap=1000)
    chr1 = [s for s in specs if s.region.startswith("chr1")]
    # chr1 is split (more than one chunk); chr2 stays whole.
    assert len(chr1) >= 2
    assert any(s.region == "chr2:0-5000" for s in specs)


def test_every_locus_owned_exactly_once(synthetic_genome):
    for strat, kw in [("scaffold", {}), ("genes", {"chunk_size": 2}),
                      ("size", {"chunk_size": 5000}),
                      ("adaptive", {"target_chunks": 3})]:
        specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                            strategy=strat, min_boundary_gap=1000, **kw)
        owned = [lid for s in specs for lid in s.locus_ids]
        assert sorted(owned) == ["G1", "G2", "G3", "G4", "G5", "G6"]
        assert len(owned) == len(set(owned))  # no locus owned twice


def test_minus_strand_locus_handled(synthetic_genome):
    # chr2's only locus is minus strand; it must still get its own whole-scaffold
    # chunk regardless of strategy (strand never changes coordinate storage).
    specs = plan_chunks(_loci(synthetic_genome), _lengths(synthetic_genome),
                        strategy="genes", chunk_size=10, min_boundary_gap=1000)
    chr2 = [s for s in specs if s.region.startswith("chr2")]
    assert len(chr2) == 1
    assert chr2[0].locus_ids == ["G6"]
    assert chr2[0].region == "chr2:0-5000"


def test_strategy_enum_accepts_strings(synthetic_genome):
    assert ChunkStrategy("scaffold") is ChunkStrategy.BY_SCAFFOLD
    assert ChunkStrategy("adaptive") is ChunkStrategy.ADAPTIVE
