"""Tests for reconcile/mikado_integrate.py — CORE IP (Phase 6). Floor: 30."""

import pytest

from helixforge.reconcile.mikado_integrate import (
    NOVEL_ID_BASE,
    IdAllocator,
    assign_gene_id,
    assign_novel_gene_id,
    assign_tier,
    build_reconciled_gene,
    intron_chain,
    match_loci,
    reciprocal_cds_overlap,
    reconcile,
)
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    HelixerLocus,
    LocusClassification,
    MikadoLocus,
    TranscriptCandidate,
)


# --- builders ---

def _h(gid, seqid, start, end, strand, exons, cds=None):
    return HelixerLocus(gid, seqid, start, end, strand, exons=exons, cds=cds)


def _mt(tid, seqid, strand, exons, cds=None, score=None, protein=None, blast=None):
    start = min(e.start for e in exons)
    end = max(e.end for e in exons)
    return TranscriptCandidate(
        tid, "L", "mikado", seqid, start, end, strand, exons,
        cds=cds, combined_score=score, protein_id=protein, blast_score=blast,
    )


def _m(lid, seqid, start, end, strand, transcripts):
    return MikadoLocus(lid, seqid, start, end, strand, transcripts)


@pytest.fixture
def scenario():
    """A scenario exercising all five correspondence cases."""
    helixer = [
        _h("H1", "chr1", 1000, 2000, "+", [Exon(1000, 1200), Exon(1500, 2000)]),
        _h("H2", "chr1", 3000, 5000, "+", [Exon(3000, 3500), Exon(4500, 5000)]),
        _h("H3a", "chr2", 1000, 1500, "+", [Exon(1000, 1500)]),
        _h("H3b", "chr2", 2000, 2500, "+", [Exon(2000, 2500)]),
        _h("H4a", "chr2", 5000, 5500, "+", [Exon(5000, 5500)]),
        _h("H4b", "chr2", 6000, 6500, "+", [Exon(6000, 6500)]),
        _h("H5", "chr3", 5000, 5500, "-", [Exon(5000, 5500)]),
    ]
    mikado = [
        _m("M1", "chr1", 1000, 2000, "+", [
            _mt("M1.1", "chr1", "+", [Exon(1000, 1200), Exon(1500, 2000)],
                cds=[CDSSegment(1050, 1200, 0), CDSSegment(1500, 1800, 0)],
                score=10.0, protein="P1")]),
        _m("M2a", "chr1", 3000, 3500, "+", [
            _mt("M2a.1", "chr1", "+", [Exon(3000, 3500)],
                cds=[CDSSegment(3000, 3300, 0)], score=8.0, protein="Pa")]),
        _m("M2b", "chr1", 4500, 5000, "+", [
            _mt("M2b.1", "chr1", "+", [Exon(4500, 5000)],
                cds=[CDSSegment(4500, 4800, 0)], score=9.0, protein="Pb")]),
        _m("M3", "chr2", 1000, 2500, "+", [
            _mt("M3.1", "chr2", "+", [Exon(1000, 1500), Exon(2000, 2500)],
                cds=[CDSSegment(1000, 1500, 0), CDSSegment(2000, 2400, 0)],
                score=12.0, protein="P3")]),
        _m("M4", "chr2", 5000, 6500, "+", [
            _mt("M4.1", "chr2", "+", [Exon(5000, 6500)],
                cds=[CDSSegment(5000, 6200, 0)], score=11.0, protein="P4")]),
        _m("M5", "chr3", 8000, 8500, "+", [
            _mt("M5.1", "chr3", "+", [Exon(8000, 8500)],
                cds=[CDSSegment(8000, 8300, 0)], score=7.0, protein="P5")]),
    ]
    classifications = [
        LocusClassification("H1", "EXPRESSED"),
        LocusClassification("H2", "EXPRESSED"),
        LocusClassification("H3a", "EXPRESSED"),
        LocusClassification("H3b", "EXPRESSED"),
        LocusClassification("H4a", "EXPRESSED"),
        LocusClassification("H4b", "EXPRESSED"),
        LocusClassification("H5", "SILENT"),
    ]
    return helixer, classifications, mikado


# --------------------------------------------------------------------------
# match_loci
# --------------------------------------------------------------------------

def test_match_one_to_one(scenario):
    helixer, _, mikado = scenario
    corr = match_loci(helixer, mikado)
    assert any(h.gene_id == "H1" and m.locus_id == "M1" for h, m in corr.one_to_one)


def test_match_split(scenario):
    helixer, _, mikado = scenario
    corr = match_loci(helixer, mikado)
    h, ms = next((h, ms) for h, ms in corr.splits if h.gene_id == "H2")
    assert {m.locus_id for m in ms} == {"M2a", "M2b"}


def test_match_merge(scenario):
    helixer, _, mikado = scenario
    corr = match_loci(helixer, mikado)
    merge_groups = {tuple(sorted(h.gene_id for h in hs)) for hs, m in corr.merges}
    assert ("H3a", "H3b") in merge_groups
    assert ("H4a", "H4b") in merge_groups


def test_match_helixer_only(scenario):
    helixer, _, mikado = scenario
    corr = match_loci(helixer, mikado)
    assert any(h.gene_id == "H5" for h in corr.helixer_only)


def test_match_novel(scenario):
    helixer, _, mikado = scenario
    corr = match_loci(helixer, mikado)
    assert any(m.locus_id == "M5" for m in corr.novel)


def test_match_strand_blocks_edge():
    helixer = [_h("H", "chr1", 100, 200, "+", [Exon(100, 200)])]
    mikado = [_m("M", "chr1", 100, 200, "-", [_mt("M.1", "chr1", "-", [Exon(100, 200)])])]
    corr = match_loci(helixer, mikado)
    assert [h.gene_id for h in corr.helixer_only] == ["H"]
    assert [m.locus_id for m in corr.novel] == ["M"]


# --------------------------------------------------------------------------
# reconcile: counts, origins, flags
# --------------------------------------------------------------------------

def test_reconcile_gene_count(scenario):
    genes, _, _ = reconcile(*scenario)
    # 1 (1:1) + 2 (split) + 1 (merge) + 2 (merge rejected) + 1 (helixer-only) = 7
    assert len(genes) == 7


def test_reconcile_origins(scenario):
    genes, _, _ = reconcile(*scenario)
    from collections import Counter
    counts = Counter(g.origin for g in genes)
    assert counts["mikado_1to1"] == 1
    assert counts["split"] == 2
    assert counts["merge"] == 1
    assert counts["helixer_backstop"] == 3  # H4a, H4b (rejected) + H5


def test_reconcile_no_helixer_lost(scenario):
    helixer, _, _ = scenario
    _, id_map, _ = reconcile(*scenario)
    for h in helixer:
        assert h.gene_id in id_map


def test_reconcile_one_to_one_gene(scenario):
    genes, _, _ = reconcile(*scenario)
    g = next(g for g in genes if g.origin == "mikado_1to1")
    assert g.gene_id == "HFG_00001"
    assert g.tier == 1  # CDS + protein


def test_reconcile_split_flags_and_ids(scenario):
    genes, _, _ = reconcile(*scenario)
    splits = sorted((g for g in genes if g.origin == "split"), key=lambda g: g.gene_id)
    assert [g.gene_id for g in splits] == ["HFG_00002_a", "HFG_00002_b"]
    assert all(any(f.name == "LOCUS_SPLIT" for f in g.flags) for g in splits)


def test_reconcile_merge_accepted(scenario):
    genes, _, _ = reconcile(*scenario)
    merge = next(g for g in genes if g.origin == "merge")
    assert "H3b" in merge.merged_from
    assert any(f.name == "LOCUS_MERGE" for f in merge.flags)


def test_reconcile_merge_keeps_lowest_hfg(scenario):
    # H3a/H3b merge happens after H1(1) and H2(2) -> lowest fresh HFG is 00003
    genes, _, _ = reconcile(*scenario)
    merge = next(g for g in genes if g.origin == "merge")
    assert merge.gene_id == "HFG_00003"


def test_reconcile_merge_rejected_no_bridge(scenario):
    genes, _, _ = reconcile(*scenario)
    rejected = [g for g in genes if any(f.name == "MERGE_REJECTED" for f in g.flags)]
    # H4a + H4b kept separate (M4 has no bridging intron)
    assert len(rejected) == 2
    assert all(g.origin == "helixer_backstop" for g in rejected)


def test_reconcile_helixer_only_flags(scenario):
    genes, _, _ = reconcile(*scenario)
    h5 = next(g for g in genes if g.seqid == "chr3")
    names = {f.name for f in h5.flags}
    assert "HELIXER_ONLY" in names
    assert "NO_EXPRESSION" in names  # SILENT
    assert h5.tier == 4


def test_reconcile_novel_not_admitted_by_default(scenario):
    genes, _, admissions = reconcile(*scenario)
    assert all(g.origin != "novel" for g in genes)
    assert any(a.reason == "novel_not_admitted" and a.transcript_id == "M5.1"
               for a in admissions)


def test_reconcile_novel_admitted_when_enabled(scenario):
    genes, _, _ = reconcile(*scenario, admit_novel=True)
    novel = next(g for g in genes if g.origin == "novel")
    assert any(f.name == "NOVEL_LOCUS" for f in novel.flags)
    from helixforge.reconcile.mikado_integrate import _hfg_num
    assert _hfg_num(novel.gene_id) >= NOVEL_ID_BASE


def test_reconcile_novel_evidence_floor_blocks(scenario):
    # M5 score is 7.0; floor 8.0 should block admission
    genes, _, _ = reconcile(*scenario, admit_novel=True, novel_evidence_floor=8.0)
    assert all(g.origin != "novel" for g in genes)


def test_reconcile_novel_evidence_floor_passes(scenario):
    genes, _, _ = reconcile(*scenario, admit_novel=True, novel_evidence_floor=5.0)
    assert any(g.origin == "novel" for g in genes)


# --------------------------------------------------------------------------
# Phase 29 D1 — paralog/tandem-array-aware merge guards
# --------------------------------------------------------------------------

from helixforge.reconcile.models import SpliceJunction  # noqa: E402


def _merge_case(strand, *, with_cds=False):
    """Two adjacent loci (1000-1500, 2000-2500) bridged by a Mikado intron.

    The bridging gap is (1500, 2000). Returns (helixer, classifications, mikado).
    Concrete literal coords; works on either strand.
    """
    if with_cds:
        ha = _h("Ha", "chrZ", 1000, 1500, strand, [Exon(1000, 1500)],
                cds=[CDSSegment(1000, 1300, 0)])
        hb = _h("Hb", "chrZ", 2000, 2500, strand, [Exon(2000, 2500)],
                cds=[CDSSegment(2000, 2300, 0)])
    else:
        ha = _h("Ha", "chrZ", 1000, 1500, strand, [Exon(1000, 1500)])
        hb = _h("Hb", "chrZ", 2000, 2500, strand, [Exon(2000, 2500)])
    helixer = [ha, hb]
    classifications = [
        LocusClassification("Ha", "EXPRESSED"),
        LocusClassification("Hb", "EXPRESSED"),
    ]
    mikado = [
        _m("MZ", "chrZ", 1000, 2500, strand, [
            _mt("MZ.1", "chrZ", strand,
                [Exon(1000, 1500), Exon(2000, 2500)],
                cds=[CDSSegment(1000, 1500, 0), CDSSegment(2000, 2400, 0)],
                score=10.0, protein="PZ")]),
    ]
    return helixer, classifications, mikado


def _bridge_junction(strand, canonical, reads=5):
    return SpliceJunction("chrZ", 1500, 2000, strand, reads, canonical=canonical)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_merge_accepts_canonical_bridge_with_reads(strand):
    helixer, cls, mikado = _merge_case(strand)
    junctions = [_bridge_junction(strand, "GT-AG", reads=5)]
    genes, id_map, _ = reconcile(helixer, cls, mikado, junctions=junctions,
                                 merge_min_gap_reads=3)
    assert any(g.origin == "merge" for g in genes)
    assert "Ha" in id_map and "Hb" in id_map


@pytest.mark.parametrize("strand", ["+", "-"])
def test_merge_rejects_noncanonical_bridge(strand):
    helixer, cls, mikado = _merge_case(strand)
    junctions = [_bridge_junction(strand, "non-canonical", reads=50)]
    genes, id_map, _ = reconcile(helixer, cls, mikado, junctions=junctions,
                                 merge_min_gap_reads=3)
    # no fusion; both loci kept as flagged backstops (no gene lost)
    assert all(g.origin != "merge" for g in genes)
    rejected = [g for g in genes if any(f.name == "MERGE_REJECTED" for f in g.flags)]
    assert len(rejected) == 2
    assert all(g.origin == "helixer_backstop" for g in rejected)
    assert "Ha" in id_map and "Hb" in id_map


def test_merge_rejects_insufficient_gap_reads():
    helixer, cls, mikado = _merge_case("+")
    junctions = [_bridge_junction("+", "GT-AG", reads=2)]  # below threshold
    genes, _, _ = reconcile(helixer, cls, mikado, junctions=junctions,
                            merge_min_gap_reads=5)
    assert all(g.origin != "merge" for g in genes)
    assert sum(any(f.name == "MERGE_REJECTED" for f in g.flags) for g in genes) == 2


def test_merge_rejects_when_no_junction_spans_gap():
    helixer, cls, mikado = _merge_case("+")
    # A junction elsewhere on the contig, not spanning the inter-genic gap.
    junctions = [SpliceJunction("chrZ", 100, 400, "+", 99, canonical="GT-AG")]
    genes, _, _ = reconcile(helixer, cls, mikado, junctions=junctions)
    assert all(g.origin != "merge" for g in genes)


def test_merge_unevaluated_motif_accepts_with_reads():
    # canonical=None (motif never evaluated) is treated as "not known non-canonical"
    helixer, cls, mikado = _merge_case("+")
    junctions = [SpliceJunction("chrZ", 1500, 2000, "+", 5, canonical=None)]
    genes, _, _ = reconcile(helixer, cls, mikado, junctions=junctions,
                            merge_min_gap_reads=3)
    assert any(g.origin == "merge" for g in genes)


def test_merge_legacy_when_no_junctions_supplied():
    # junctions=None preserves the legacy any-Mikado-bridging-intron behavior.
    helixer, cls, mikado = _merge_case("+")
    genes, _, _ = reconcile(helixer, cls, mikado)
    assert any(g.origin == "merge" for g in genes)


class _DictGenome:
    """Minimal genome accessor returning preset sequences by (seqid, start, end)."""

    def __init__(self, mapping):
        self._m = mapping

    def get_sequence(self, seqid, start, end):
        return self._m[(seqid, start, end)]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_merge_rejects_high_identity_paralogs(strand):
    helixer, cls, mikado = _merge_case(strand, with_cds=True)
    same = "ACGTGCAATTGC" * 25  # 300 bp, identical for both loci → identity 1.0
    genome = _DictGenome({
        ("chrZ", 1000, 1300): same,
        ("chrZ", 2000, 2300): same,
    })
    junctions = [_bridge_junction(strand, "GT-AG", reads=50)]  # good bridge present
    genes, id_map, _ = reconcile(
        helixer, cls, mikado, junctions=junctions,
        genome=genome, paralog_identity_threshold=0.9,
    )
    # high CDS identity ⇒ recent paralogs ⇒ NOT fused, despite a valid bridge
    assert all(g.origin != "merge" for g in genes)
    assert sum(any(f.name == "MERGE_REJECTED" for f in g.flags) for g in genes) == 2
    assert "Ha" in id_map and "Hb" in id_map


def test_merge_low_identity_paralogs_still_fuse():
    helixer, cls, mikado = _merge_case("+", with_cds=True)
    genome = _DictGenome({
        ("chrZ", 1000, 1300): "ACGTACGTACGT" * 25,
        ("chrZ", 2000, 2300): "TTGGCCAATTGG" * 25,  # disjoint k-mers → identity ~0
    })
    junctions = [_bridge_junction("+", "GT-AG", reads=50)]
    genes, _, _ = reconcile(
        helixer, cls, mikado, junctions=junctions,
        genome=genome, paralog_identity_threshold=0.9,
    )
    assert any(g.origin == "merge" for g in genes)


# --------------------------------------------------------------------------
# Regression — a rejected merge must emit UNIQUE ids even when a prior run
# accepted it and coupled the constituent loci to one shared HFG in id_map.
# (assessment-v4.md §2.3; the scorer:stats IndexError crash on TAIR10.)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_merge_rejected_decollides_shared_hfg(strand):
    """A prior accepted merge maps Ha+Hb → one HFG; on rejection they must split."""
    helixer, cls, mikado = _merge_case(strand)
    # Simulate the persisted id_map from a run where the merge was ACCEPTED:
    # both constituent loci share the lowest HFG (CLAUDE.md §11).
    id_map = {"Ha": "HFG_00010", "Hb": "HFG_00010"}
    junctions = [_bridge_junction(strand, "non-canonical", reads=50)]  # force reject
    genes, out_map, _ = reconcile(
        helixer, cls, mikado, id_map=id_map, junctions=junctions,
        merge_min_gap_reads=3,
    )
    rejected = [g for g in genes if any(f.name == "MERGE_REJECTED" for f in g.flags)]
    assert len(rejected) == 2
    # Distinct gene ids AND distinct transcript ids (the actual file-level bug).
    gene_ids = [g.gene_id for g in rejected]
    tx_ids = [t.transcript_id for g in rejected for t in g.transcripts]
    assert len(set(gene_ids)) == 2
    assert len(set(tx_ids)) == 2
    # The representative (lowest-start) locus keeps the canonical historical id.
    assert out_map["Ha"] == "HFG_00010"
    assert out_map["Hb"] != "HFG_00010"


def test_merge_rejected_fresh_run_ids_unchanged():
    """With a fresh id_map each released locus already gets its own number."""
    helixer, cls, mikado = _merge_case("+")
    junctions = [_bridge_junction("+", "non-canonical", reads=50)]
    genes, _, _ = reconcile(helixer, cls, mikado, junctions=junctions,
                            merge_min_gap_reads=3)
    rejected = [g for g in genes if any(f.name == "MERGE_REJECTED" for f in g.flags)]
    assert len({g.gene_id for g in rejected}) == 2


# --------------------------------------------------------------------------
# Phase 29 D3 — merge-rejected loci tier sensibly (no gene lost)
# --------------------------------------------------------------------------

def test_merge_rejected_locus_tiers_as_backstop_expressed():
    helixer, cls, mikado = _merge_case("+")  # both EXPRESSED, no CDS
    junctions = [_bridge_junction("+", "non-canonical", reads=50)]
    genes, _, _ = reconcile(helixer, cls, mikado, junctions=junctions)
    rejected = [g for g in genes if any(f.name == "MERGE_REJECTED" for f in g.flags)]
    assert len(rejected) == 2
    # CDS-less + EXPRESSED ⇒ Tier 3 (not a coding Tier 1/2)
    assert all(g.tier == 3 for g in rejected)


def test_merge_rejected_locus_tiers_silent_is_tier4():
    helixer, _cls, mikado = _merge_case("+")
    cls = [LocusClassification("Ha", "SILENT"), LocusClassification("Hb", "SILENT")]
    junctions = [_bridge_junction("+", "non-canonical", reads=50)]
    genes, _, _ = reconcile(helixer, cls, mikado, junctions=junctions)
    rejected = [g for g in genes if any(f.name == "MERGE_REJECTED" for f in g.flags)]
    assert len(rejected) == 2
    assert all(g.tier == 4 for g in rejected)


# --------------------------------------------------------------------------
# ID stability
# --------------------------------------------------------------------------

def test_id_stability_across_runs(scenario):
    genes1, id_map1, _ = reconcile(*scenario)
    ids1 = {g.gene_id for g in genes1}
    genes2, id_map2, _ = reconcile(*scenario, id_map=id_map1)
    ids2 = {g.gene_id for g in genes2}
    assert ids1 == ids2


def test_id_stability_same_mapping(scenario):
    _, id_map1, _ = reconcile(*scenario)
    _, id_map2, _ = reconcile(*scenario, id_map=id_map1)
    for key, value in id_map1.items():
        assert id_map2[key] == value


def test_id_map_not_mutated_in_place(scenario):
    helixer, cls, mikado = scenario
    original = {"H1": "HFG_00001"}
    reconcile(helixer, cls, mikado, id_map=original)
    assert original == {"H1": "HFG_00001"}  # caller's dict untouched


def test_existing_id_reused(scenario):
    helixer, cls, mikado = scenario
    genes, _, _ = reconcile(helixer, cls, mikado, id_map={"H1": "HFG_00042"})
    g = next(g for g in genes if g.origin == "mikado_1to1")
    assert g.gene_id == "HFG_00042"


# --------------------------------------------------------------------------
# assign_gene_id / assign_tier / helpers
# --------------------------------------------------------------------------

def test_assign_gene_id_new():
    id_map = {}
    h = _h("Hx", "chr1", 0, 100, "+", [Exon(0, 100)])
    assert assign_gene_id(h, id_map) == "HFG_00001"
    assert id_map["Hx"] == "HFG_00001"


def test_assign_gene_id_stable():
    id_map = {"Hx": "HFG_00007"}
    h = _h("Hx", "chr1", 0, 100, "+", [Exon(0, 100)])
    assert assign_gene_id(h, id_map) == "HFG_00007"


def test_assign_gene_id_merge_keeps_lowest():
    id_map = {"Ha": "HFG_00005", "Hb": "HFG_00002"}
    ha = _h("Ha", "chr1", 0, 100, "+", [Exon(0, 100)])
    chosen = assign_gene_id(ha, id_map, merged_from=["Hb"])
    assert chosen == "HFG_00002"
    assert id_map["Ha"] == "HFG_00002"


def test_assign_tier_matrix():
    cds = [CDSSegment(0, 99, 0)]
    p_full = _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds, protein="P")
    p_nocds = _mt("t", "chr1", "+", [Exon(0, 99)])
    expr = LocusClassification("x", "EXPRESSED")
    silent = LocusClassification("x", "SILENT")
    assert assign_tier(p_full, expr, "mikado_1to1") == 1
    assert assign_tier(p_nocds, expr, "split") == 2
    cds_noprot = _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds)
    assert assign_tier(cds_noprot, expr, "merge") == 2
    assert assign_tier(p_nocds, expr, "helixer_backstop") == 3
    assert assign_tier(p_nocds, silent, "helixer_backstop") == 4


def test_assign_tier_blast_score_grants_tier1():
    # Tier 1 via BLAST homology (no protein_id accession, as in real Mikado output)
    cds = [CDSSegment(0, 99, 0)]
    expr = LocusClassification("x", "EXPRESSED")
    blast_hit = _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds, blast=816.0)
    assert assign_tier(blast_hit, expr, "mikado_1to1") == 1
    # CDS but zero BLAST and no accession -> Tier 2
    no_homol = _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds, blast=0.0)
    assert assign_tier(no_homol, expr, "split") == 2
    # homology but no CDS -> Tier 2 (Tier 1 needs both)
    no_cds = _mt("t", "chr1", "+", [Exon(0, 99)], blast=500.0)
    assert assign_tier(no_cds, expr, "mikado_1to1") == 2


def test_has_homology_property():
    cds = [CDSSegment(0, 99, 0)]
    assert _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds, protein="P").has_homology
    assert _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds, blast=12.0).has_homology
    assert not _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds).has_homology
    assert not _mt("t", "chr1", "+", [Exon(0, 99)], cds=cds, blast=0.0).has_homology


def test_intron_chain():
    t = _mt("t", "chr1", "+", [Exon(0, 100), Exon(200, 300), Exon(400, 500)])
    assert intron_chain(t) == ((100, 200), (300, 400))


def test_intron_chain_single_exon():
    t = _mt("t", "chr1", "+", [Exon(0, 100)])
    assert intron_chain(t) == ()


def test_reciprocal_cds_overlap_no_cds():
    a = _mt("a", "chr1", "+", [Exon(0, 100)])
    b = _mt("b", "chr1", "+", [Exon(0, 100)])
    assert reciprocal_cds_overlap(a, b) == 0.0


def test_reciprocal_cds_overlap_full():
    a = _mt("a", "chr1", "+", [Exon(0, 300)], cds=[CDSSegment(0, 300, 0)])
    b = _mt("b", "chr1", "+", [Exon(0, 300)], cds=[CDSSegment(0, 300, 0)])
    assert reciprocal_cds_overlap(a, b) == 1.0


# --------------------------------------------------------------------------
# admissions + primary
# --------------------------------------------------------------------------

def test_admissions_recorded(scenario):
    _, _, admissions = reconcile(*scenario)
    # at least one admitted primary per built gene (7 genes)
    primaries = [a for a in admissions if a.reason == "primary"]
    assert len(primaries) == 7


def test_primary_is_top_score():
    # two isoforms in one mikado locus; higher score becomes .1
    helixer = [_h("H", "chr1", 1000, 2000, "+", [Exon(1000, 1200), Exon(1500, 2000)])]
    mikado = [_m("M", "chr1", 1000, 2000, "+", [
        _mt("lo", "chr1", "+", [Exon(1000, 2000)], score=3.0),
        _mt("hi", "chr1", "+", [Exon(1000, 1200), Exon(1500, 2000)], score=9.0),
    ])]
    cls = [LocusClassification("H", "EXPRESSED")]
    genes, _, _ = reconcile(helixer, cls, mikado)
    g = genes[0]
    primary = next(t for t in g.transcripts if t.is_primary)
    assert primary.transcript_id == "HFG_00001.1"
    assert primary.combined_score == 9.0


def test_isoform_ids_ordered_by_score():
    helixer = [_h("H", "chr1", 1000, 2000, "+", [Exon(1000, 1200), Exon(1500, 2000)])]
    mikado = [_m("M", "chr1", 1000, 2000, "+", [
        _mt("lo", "chr1", "+", [Exon(1000, 2000)], score=3.0),
        _mt("hi", "chr1", "+", [Exon(1000, 1200), Exon(1500, 2000)], score=9.0),
    ])]
    cls = [LocusClassification("H", "EXPRESSED")]
    genes, _, _ = reconcile(helixer, cls, mikado)
    ids = sorted(t.transcript_id for t in genes[0].transcripts)
    assert ids == ["HFG_00001.1", "HFG_00001.2"]


def test_build_reconciled_gene_backstop_single_isoform():
    h = _h("H", "chr1", 1000, 1500, "+", [Exon(1000, 1500)])
    cls = LocusClassification("H", "LOW")
    gene = build_reconciled_gene(h, None, cls, "helixer_backstop", {})
    assert len(gene.transcripts) == 1
    assert gene.transcripts[0].source == "helixer"
    assert gene.transcripts[0].cds is None  # CDS comes in Phase 7
    assert gene.tier == 3


# --------------------------------------------------------------------------
# IdAllocator — disjoint-range allocation (Phase 16 §C3)
# --------------------------------------------------------------------------

def test_id_allocator_default_reproduces_current_ids(scenario):
    """Passing the default IdAllocator must yield byte-identical ids."""
    genes_default, idmap_default, _ = reconcile(*scenario)
    genes_alloc, idmap_alloc, _ = reconcile(*scenario, allocator=IdAllocator())
    assert idmap_default == idmap_alloc
    assert [g.gene_id for g in genes_default] == [g.gene_id for g in genes_alloc]


def test_id_allocator_default_first_number_is_one():
    h = _h("H", "chr1", 1000, 1500, "+", [Exon(1000, 1500)])
    assert assign_gene_id(h, {}, allocator=IdAllocator()) == "HFG_00001"


def test_id_allocator_nondefault_base_in_range():
    h = _h("H", "chr1", 1000, 1500, "+", [Exon(1000, 1500)])
    assert assign_gene_id(h, {}, allocator=IdAllocator(base=500)) == "HFG_00500"


def test_id_allocator_nondefault_novel_base():
    m = _m("M9", "chr1", 1000, 1500, "+", [_mt("M9.1", "chr1", "+", [Exon(1000, 1500)])])
    hfg = assign_novel_gene_id(m, {}, allocator=IdAllocator(novel_base=90500))
    assert hfg == "HFG_90500"


def test_id_allocator_next_number_skips_used():
    alloc = IdAllocator(base=100)
    # 100 already taken in the seed map → next free in the range is 101.
    assert alloc.next_number({"X": "HFG_00100"}) == 101


def test_reconcile_allocator_base_keeps_ids_in_range(scenario):
    """A chunk-style allocator confines all Helixer-anchored ids to its range."""
    genes, _, _ = reconcile(*scenario, allocator=IdAllocator(base=1000))
    anchored = [g for g in genes if g.origin != "novel"]
    # gene_id is "HFG_01000" or a split child "HFG_01000_a"; the number is part[1].
    nums = [int(g.gene_id.split("_")[1]) for g in anchored]
    assert all(n >= 1000 for n in nums)


# --------------------------------------------------------------------------
# IdAllocator — stateful O(1) cursor (Phase 20 §1.1). Must reproduce the old
# "lowest free number at/above base" scanner exactly.
# --------------------------------------------------------------------------

def _old_scanner(id_map, base=1, novel_base=NOVEL_ID_BASE, novel=False):
    """The pre-Phase-20 stateless scanner (rebuilds `used` every call)."""
    used = set()
    for value in id_map.values():
        try:
            used.add(int(value.split("_")[1]))
        except (IndexError, ValueError):
            continue
    n = novel_base if novel else base
    while n in used:
        n += 1
    return n


def test_stateful_allocator_matches_old_scanner_with_foreign_seed():
    # Seed carries an in-range used id (3), a far out-of-range/foreign HFG
    # (99999), and a non-HFG garbage value that must be ignored.
    seed = {"pre1": "HFG_00003", "pre2": "HFG_99999", "junk": "not_an_hfg"}
    alloc = IdAllocator()
    new_map, old_map = dict(seed), dict(seed)
    new_seq, old_seq = [], []
    for i in range(20):
        n_new = alloc.next_number(new_map)
        new_map[f"new{i}"] = f"HFG_{n_new:05d}"
        n_old = _old_scanner(old_map)
        old_map[f"new{i}"] = f"HFG_{n_old:05d}"
        new_seq.append(n_new)
        old_seq.append(n_old)
    assert new_seq == old_seq
    # concrete: base=1, used {3, 99999} -> 1, 2, skip 3, 4, 5, ...
    assert new_seq[:5] == [1, 2, 4, 5, 6]


def test_stateful_allocator_skips_in_range_seed():
    alloc = IdAllocator(base=1)
    # 2 is already taken in the seed -> 1 then skip 2 -> 3.
    m = {"x": "HFG_00002"}
    a = alloc.next_number(m)
    m[f"k{a}"] = f"HFG_{a:05d}"
    b = alloc.next_number(m)
    assert (a, b) == (1, 3)


def test_stateful_allocator_novel_range_honored():
    alloc = IdAllocator()
    m = {}
    a = alloc.next_number(m, novel=True)
    m["k1"] = f"HFG_{a:05d}"
    b = alloc.next_number(m, novel=True)
    assert (a, b) == (NOVEL_ID_BASE, NOVEL_ID_BASE + 1)


def test_stateful_allocator_default_bases_unchanged():
    alloc = IdAllocator()
    assert alloc.base == 1
    assert alloc.novel_base == NOVEL_ID_BASE == 90000
    assert alloc.next_number({}) == 1
    assert IdAllocator().next_number({}, novel=True) == 90000


def test_stateful_allocator_reserved_range_floor():
    alloc = IdAllocator(base=500)
    m = {}
    nums = []
    for i in range(10):
        n = alloc.next_number(m)
        m[f"k{i}"] = f"HFG_{n:05d}"
        nums.append(n)
        assert n >= 500
    assert nums == list(range(500, 510))


def test_stateful_allocators_disjoint_ranges_never_collide():
    # Two chunk allocators with disjoint reserved bases, seeded from a shared
    # master id_map, must never hand out colliding numbers (CLAUDE.md §C3).
    master = {"existing": "HFG_00001"}
    chunk_a = IdAllocator(base=10)
    chunk_b = IdAllocator(base=100)
    a_map, b_map = dict(master), dict(master)
    a_nums, b_nums = [], []
    for i in range(5):
        na = chunk_a.next_number(a_map)
        a_map[f"a{i}"] = f"HFG_{na:05d}"
        a_nums.append(na)
        nb = chunk_b.next_number(b_map)
        b_map[f"b{i}"] = f"HFG_{nb:05d}"
        b_nums.append(nb)
    assert a_nums == [10, 11, 12, 13, 14]
    assert b_nums == [100, 101, 102, 103, 104]
    assert set(a_nums).isdisjoint(b_nums)


def test_stateful_allocator_composes_with_reserve_id_ranges():
    from helixforge.parallel.plan import Chunk, Plan, reserve_id_ranges

    plan = Plan(chunks=[
        Chunk(chunk_id="c0", regions=["chr1:0-1000"], num_loci=3),
        Chunk(chunk_id="c1", regions=["chr1:2000-3000"], num_loci=2),
    ])
    reserve_id_ranges(plan, id_map={})
    c0, c1 = plan.chunks
    # Each chunk drives an allocator from its reserved base; numbers stay in range.
    alloc0 = IdAllocator(base=c0.id_base)
    alloc1 = IdAllocator(base=c1.id_base)
    m0, m1 = {}, {}
    n0 = [alloc0.next_number(m0) for _ in range(c0.num_loci)]
    n1 = [alloc1.next_number(m1) for _ in range(c1.num_loci)]
    assert all(c0.id_range[0] <= n < c0.id_range[1] for n in n0)
    assert all(c1.id_range[0] <= n < c1.id_range[1] for n in n1)
    assert set(n0).isdisjoint(n1)


def test_allocator_no_full_rescan_per_call(monkeypatch):
    """Probe: _hfg_num is called only while seeding, never per allocation."""
    import helixforge.reconcile.mikado_integrate as mi

    calls = {"n": 0}
    orig = mi._hfg_num

    def counting(hfg):
        calls["n"] += 1
        return orig(hfg)

    monkeypatch.setattr(mi, "_hfg_num", counting)
    seed = {f"k{i}": f"HFG_{i + 1:05d}" for i in range(50)}
    alloc = mi.IdAllocator()
    m = dict(seed)
    for i in range(1000):
        n = alloc.next_number(m)
        m[f"new{i}"] = f"HFG_{n:05d}"
    # Seeded once over the 50-entry map; the 1000 allocations add zero rescans.
    # (The old scanner would have called _hfg_num ~550k times.)
    assert calls["n"] == 50


def test_allocator_empty_seed_zero_hfg_num_calls(monkeypatch):
    import helixforge.reconcile.mikado_integrate as mi

    calls = {"n": 0}
    orig = mi._hfg_num
    monkeypatch.setattr(mi, "_hfg_num", lambda h: (calls.__setitem__("n", calls["n"] + 1), orig(h))[1])

    alloc = mi.IdAllocator()
    m = {}
    for i in range(500):
        n = alloc.next_number(m)
        m[f"k{i}"] = f"HFG_{n:05d}"
    assert calls["n"] == 0
    assert m["k0"] == "HFG_00001"


# --------------------------------------------------------------------------
# Intrinsic Helixer ORF carried into the backstop transcript (TPM/biotype fix)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_backstop_inherits_helixer_cds(strand):
    # A Helixer-only gene must carry the model's own CDS so it is coding on the
    # strength of its intrinsic ORF (450 nt = mod-3 complete). Tier 2 (CDS, no
    # homology) instead of the silent-backstop Tier 3/4.
    h = _h("H1", "chr1", 1000, 1600, strand,
           [Exon(1000, 1300), Exon(1400, 1600)],
           cds=[CDSSegment(1000, 1300, 0), CDSSegment(1400, 1550, 0)])
    gene = build_reconciled_gene(
        h, None, LocusClassification("H1", "SILENT"), "helixer_backstop", {})
    primary = gene.transcripts[0]
    assert primary.cds is not None
    assert [(c.start, c.end) for c in primary.cds] == [(1000, 1300), (1400, 1550)]
    assert primary.total_cds_length == 450
    assert gene.tier == 2  # CDS, no homology
    assert not primary.cds_partial


@pytest.mark.parametrize("strand", ["+", "-"])
def test_backstop_partial_helixer_cds_flags_partial(strand):
    # A non-mod-3 Helixer CDS (200 nt) is a partial ORF: CDS still carried,
    # cds_partial True, PARTIAL_ORF flag emitted.
    h = _h("H1", "chr1", 1000, 1600, strand,
           [Exon(1000, 1200)], cds=[CDSSegment(1000, 1200, 0)])
    gene = build_reconciled_gene(
        h, None, LocusClassification("H1", "SILENT"), "helixer_backstop", {})
    primary = gene.transcripts[0]
    assert primary.cds is not None
    assert primary.cds_partial
    assert any(f.name == "PARTIAL_ORF" for f in gene.flags)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_backstop_no_helixer_cds_stays_cdsless(strand):
    # A genuinely CDS-less Helixer locus stays CDS-less (true non-coding).
    h = _h("H1", "chr1", 1000, 1600, strand,
           [Exon(1000, 1200), Exon(1400, 1600)], cds=None)
    gene = build_reconciled_gene(
        h, None, LocusClassification("H1", "EXPRESSED"), "helixer_backstop", {})
    assert gene.transcripts[0].cds is None
