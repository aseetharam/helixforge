"""Tests for parallel/plan.py — partition + HFG-range reservation (Phase 16 D1).

Floor: 14. Concrete literal coordinates; both strands present in the fixture.
The core invariants: boundaries only in gaps >= threshold, no gene ever split,
reserved ranges disjoint and sum-cover all loci, reseed keeps existing HFGs.
"""

import pytest

from helixforge.parallel.plan import (
    Chunk,
    Plan,
    _region_bounds,
    partition_by_strategy,
    partition_genome,
    read_plan,
    reserve_id_ranges,
    write_plan,
)
from helixforge.reconcile.locus import load_helixer_loci


def _locus_in_one_region(locus, plan):
    """How many chunk regions contain ``locus`` wholly (must be exactly 1)."""
    owners = 0
    for chunk in plan.chunks:
        for region in chunk.regions:
            seqid, lo, hi = _region_bounds(region)
            if seqid != locus.seqid:
                continue
            if lo is None or (lo <= locus.start and locus.end <= hi):
                owners += 1
    return owners


# --------------------------------------------------------------------------
# partition_genome
# --------------------------------------------------------------------------

def test_partition_cut_points_in_gaps(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    # chr1 → 3 chunks at midpoints of the two big gaps; chr2 → 1 whole-scaffold.
    regions = [c.regions[0] for c in plan.chunks]
    assert regions == ["chr1:0-3500", "chr1:3500-10000", "chr1:10000-20000", "chr2"]


def test_partition_chunk_count(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    assert len(plan) == 4


def test_partition_num_loci_per_chunk(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    assert [c.num_loci for c in plan.chunks] == [2, 2, 1, 1]


def test_partition_no_gene_split(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    loci = load_helixer_loci(synthetic_genome["gff3"])
    for locus in loci:
        assert _locus_in_one_region(locus, plan) == 1


def test_partition_covers_all_loci(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    assert sum(c.num_loci for c in plan.chunks) == synthetic_genome["total_loci"]
    assert plan.total_loci == 6


def test_partition_small_gap_not_cut(synthetic_genome):
    # With a huge min_boundary_gap no chr1 gap qualifies → one chunk per scaffold.
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=1, min_boundary_gap=100_000,
    )
    assert len(plan) == 2  # chr1 (whole) + chr2 (whole)
    assert plan.chunks[0].regions == ["chr1"]
    assert plan.chunks[0].num_loci == 5


def test_partition_minus_strand_locus_kept(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    chr2 = [c for c in plan.chunks if c.regions == ["chr2"]]
    assert len(chr2) == 1
    assert chr2[0].locus_ids == ["G6"]


def test_partition_target_chunks(synthetic_genome):
    # target_chunks=3 over 6 loci → target_loci_per_chunk = ceil(6/3) = 2.
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_chunks=3, min_boundary_gap=1000,
    )
    assert [c.num_loci for c in plan.chunks] == [2, 2, 1, 1]


def test_partition_default_one_chunk_per_scaffold(synthetic_genome):
    # No target → never cut a scaffold.
    plan = partition_genome(synthetic_genome["fasta"], synthetic_genome["gff3"])
    assert len(plan) == 2
    assert [c.num_loci for c in plan.chunks] == [5, 1]


def test_partition_accepts_fai(synthetic_genome):
    # Build a .fai-style index and partition off it (no FASTA read needed).
    fai = synthetic_genome["tmp_path"] / "genome.fasta.fai"
    fai.write_text("chr1\t20000\t6\t20000\t20001\nchr2\t5000\t6\t5000\t5001\n")
    plan = partition_genome(
        str(fai), synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    assert plan.chunks[-1].regions == ["chr2"]


# --------------------------------------------------------------------------
# reserve_id_ranges
# --------------------------------------------------------------------------

def test_reserve_ranges_contiguous_disjoint(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(plan)
    assert [c.id_range for c in plan.chunks] == [(1, 3), (3, 5), (5, 6), (6, 7)]
    assert [c.id_base for c in plan.chunks] == [1, 3, 5, 6]


def test_reserve_ranges_sum_cover_all_loci(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(plan)
    span = sum(hi - lo for lo, hi in (c.id_range for c in plan.chunks))
    assert span == synthetic_genome["total_loci"]


def test_reserve_novel_ranges_disjoint(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(plan)
    assert [c.novel_range for c in plan.chunks] == [
        (90000, 90002), (90002, 90004), (90004, 90005), (90005, 90006),
    ]
    assert all(c.novel_base >= 90000 for c in plan.chunks)


def test_reserve_reseed_from_master_keeps_above_existing(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    # Master already used HFG_00010 (and a novel HFG_90007); new ranges start above.
    reserve_id_ranges(plan, id_map={"G1": "HFG_00010", "NOVEL_x": "HFG_90007"})
    assert plan.chunks[0].id_base == 11
    assert plan.chunks[0].id_range == (11, 13)
    assert plan.chunks[0].novel_base == 90008


def test_reserve_ranges_disjoint_after_reseed(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(plan, id_map={"G1": "HFG_00010"})
    bases = [c.id_range for c in plan.chunks]
    for (a_lo, a_hi), (b_lo, b_hi) in zip(bases, bases[1:]):
        assert b_lo >= a_hi  # contiguous, non-overlapping


# --------------------------------------------------------------------------
# write_plan / read_plan round-trip
# --------------------------------------------------------------------------

def test_plan_roundtrip(synthetic_genome):
    plan = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(plan)
    path = synthetic_genome["tmp_path"] / "plan.json"
    write_plan(plan, path)
    back = read_plan(path)
    assert len(back) == len(plan)
    assert [c.regions for c in back.chunks] == [c.regions for c in plan.chunks]
    assert [c.id_range for c in back.chunks] == [c.id_range for c in plan.chunks]
    assert back.chunks[0].novel_range == (90000, 90002)


def test_partition_invalid_target_raises(synthetic_genome):
    with pytest.raises(ValueError):
        partition_genome(
            synthetic_genome["fasta"], synthetic_genome["gff3"],
            target_loci_per_chunk=0,
        )


# --------------------------------------------------------------------------
# Phase 32 D3: small-scaffold bin-packing (fragmented draft assemblies)
# --------------------------------------------------------------------------

def _write_fragmented(tmp_path, n_scaffolds=10, scaf_len=800):
    """N tiny single-locus scaffolds (alternating strand) — a fragmented draft."""
    fasta = tmp_path / "frag.fasta"
    gff3 = tmp_path / "frag.gff3"
    fa_lines, gff_lines = [], ["##gff-version 3"]
    for i in range(n_scaffolds):
        seqid = f"ctg{i:03d}"
        strand = "+" if i % 2 == 0 else "-"
        fa_lines.append(f">{seqid}")
        fa_lines.append("A" * scaf_len)
        # one gene 101-300 (1-based) -> internal (100, 300)
        gff_lines.append(f"{seqid}\tTest\tgene\t101\t300\t.\t{strand}\t.\tID={seqid}_g")
        gff_lines.append(f"{seqid}\tTest\tmRNA\t101\t300\t.\t{strand}\t.\tID={seqid}_g.t1;Parent={seqid}_g")
        gff_lines.append(f"{seqid}\tTest\texon\t101\t300\t.\t{strand}\t.\tID={seqid}_g.t1.exon1;Parent={seqid}_g.t1")
    fasta.write_text("\n".join(fa_lines) + "\n")
    gff3.write_text("\n".join(gff_lines) + "\n")
    return str(fasta), str(gff3)


def test_pack_many_small_scaffolds_into_few_chunks(tmp_path):
    fasta, gff3 = _write_fragmented(tmp_path, n_scaffolds=10)
    plan = partition_genome(
        fasta, gff3, pack_small_scaffolds=True, max_scaffolds_per_chunk=3,
    )
    # 10 scaffolds packed 3+3+3+1 -> 4 chunks, not 10.
    assert len(plan) == 4
    assert [len(c.regions) for c in plan.chunks] == [3, 3, 3, 1]
    # bare-seqid (whole-scaffold) regions
    assert plan.chunks[0].regions == ["ctg000", "ctg001", "ctg002"]


def test_pack_no_gene_split_and_covers_all(tmp_path):
    fasta, gff3 = _write_fragmented(tmp_path, n_scaffolds=10)
    plan = partition_genome(
        fasta, gff3, pack_small_scaffolds=True, max_scaffolds_per_chunk=4,
    )
    loci = load_helixer_loci(gff3)
    for locus in loci:
        assert _locus_in_one_region(locus, plan) == 1
    assert sum(c.num_loci for c in plan.chunks) == len(loci) == 10


def test_pack_ranges_disjoint_and_cover_all_loci(tmp_path):
    fasta, gff3 = _write_fragmented(tmp_path, n_scaffolds=10)
    plan = partition_genome(
        fasta, gff3, pack_small_scaffolds=True, max_scaffolds_per_chunk=3,
    )
    reserve_id_ranges(plan)  # asserts disjointness internally
    # the union of reserved HFG numbers must cover all 10 loci exactly
    total_range = sum(hi - lo for (lo, hi) in (c.id_range for c in plan.chunks))
    assert total_range == 10
    # ranges are contiguous + disjoint
    ordered = sorted((c.id_range for c in plan.chunks))
    for (a_lo, a_hi), (b_lo, b_hi) in zip(ordered, ordered[1:]):
        assert b_lo >= a_hi


def test_pack_target_loci_flushes_bin(tmp_path):
    fasta, gff3 = _write_fragmented(tmp_path, n_scaffolds=10)
    # target 2 loci/chunk -> each packed chunk holds 2 single-locus scaffolds.
    plan = partition_genome(
        fasta, gff3, target_loci_per_chunk=2, pack_small_scaffolds=True,
        max_scaffolds_per_chunk=100,
    )
    assert len(plan) == 5
    assert all(c.num_loci == 2 for c in plan.chunks)
    assert all(len(c.regions) == 2 for c in plan.chunks)


def test_pack_default_off_is_one_chunk_per_scaffold(tmp_path):
    fasta, gff3 = _write_fragmented(tmp_path, n_scaffolds=10)
    plan = partition_genome(fasta, gff3)  # packing off by default
    assert len(plan) == 10
    assert all(len(c.regions) == 1 for c in plan.chunks)


def test_pack_large_scaffold_not_packed(tmp_path):
    # One big chromosome + several small contigs: the big one keeps its own
    # chunk; the small contigs are packed together.
    fasta = tmp_path / "mix.fasta"
    gff3 = tmp_path / "mix.gff3"
    fa, gf = [], ["##gff-version 3"]
    # big chromosome 2 Mb (> 1 Mb small_scaffold_bp), one gene
    fa.append(">chrBig"); fa.append("A" * 2_000_000)
    gf.append("chrBig\tTest\tgene\t101\t300\t.\t+\t.\tID=big_g")
    gf.append("chrBig\tTest\tmRNA\t101\t300\t.\t+\t.\tID=big_g.t1;Parent=big_g")
    gf.append("chrBig\tTest\texon\t101\t300\t.\t+\t.\tID=big_g.t1.exon1;Parent=big_g.t1")
    for i in range(4):
        seqid = f"ctg{i:03d}"
        fa.append(f">{seqid}"); fa.append("A" * 800)
        gf.append(f"{seqid}\tTest\tgene\t101\t300\t.\t-\t.\tID={seqid}_g")
        gf.append(f"{seqid}\tTest\tmRNA\t101\t300\t.\t-\t.\tID={seqid}_g.t1;Parent={seqid}_g")
        gf.append(f"{seqid}\tTest\texon\t101\t300\t.\t-\t.\tID={seqid}_g.t1.exon1;Parent={seqid}_g.t1")
    fasta.write_text("\n".join(fa) + "\n")
    gff3.write_text("\n".join(gf) + "\n")
    plan = partition_genome(
        str(fasta), str(gff3), pack_small_scaffolds=True, max_scaffolds_per_chunk=10,
    )
    # chrBig is its own single-region chunk; the 4 small contigs share one chunk.
    assert len(plan) == 2
    big = [c for c in plan.chunks if c.regions == ["chrBig"]]
    packed = [c for c in plan.chunks if len(c.regions) == 4]
    assert len(big) == 1
    assert len(packed) == 1
    assert packed[0].regions == ["ctg000", "ctg001", "ctg002", "ctg003"]


def test_pack_invalid_max_scaffolds_raises(tmp_path):
    fasta, gff3 = _write_fragmented(tmp_path, n_scaffolds=3)
    with pytest.raises(ValueError):
        partition_genome(
            fasta, gff3, pack_small_scaffolds=True, max_scaffolds_per_chunk=0,
        )


# --------------------------------------------------------------------------
# partition_by_strategy — v1 strategies + id reservation (realign)
# --------------------------------------------------------------------------

def test_strategy_plan_scaffold_default(synthetic_genome):
    plan = partition_by_strategy(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        strategy="scaffold", min_boundary_gap=1000,
    )
    assert plan.strategy == "scaffold"
    assert [c.regions for c in plan.chunks] == [["chr1:0-20000"], ["chr2:0-5000"]]
    assert plan.total_loci == 6


def test_strategy_plan_genes_reserves_disjoint_ranges(synthetic_genome):
    plan = partition_by_strategy(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        strategy="genes", chunk_size=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(plan)
    # Every locus covered exactly once across the strategy chunks.
    loci = load_helixer_loci(synthetic_genome["gff3"])
    for locus in loci:
        assert _locus_in_one_region(locus, plan) == 1
    # Ranges laid end to end and pairwise disjoint.
    bases = [c.id_base for c in plan.chunks]
    assert bases == [1, 3, 5, 6]


def test_strategy_plan_round_trips_through_json(synthetic_genome, tmp_path):
    plan = partition_by_strategy(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        strategy="adaptive", target_chunks=3, min_boundary_gap=1000,
    )
    reserve_id_ranges(plan)
    out = tmp_path / "plan.json"
    write_plan(plan, out)
    loaded = read_plan(out)
    assert loaded.strategy == "adaptive"
    assert [c.regions for c in loaded.chunks] == [c.regions for c in plan.chunks]
    assert [c.id_base for c in loaded.chunks] == [c.id_base for c in plan.chunks]
