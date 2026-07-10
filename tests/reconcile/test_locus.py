"""Tests for reconcile/locus.py (Phase 3). Floor: 30."""

import pytest

from helixforge.io.hdf5 import HDF5ConfidenceReader
from helixforge.reconcile.locus import (
    build_locus_index,
    build_locus_index_by_scaffold,
    enrich_with_confidence,
    filter_loci_by_region,
    load_helixer_loci,
    merge_overlapping_loci,
    _merge_exon_lists,
)
from helixforge.reconcile.models import Exon, HelixerLocus


def _locus(gene_id, seqid, start, end, strand, confidence=None, exons=None):
    return HelixerLocus(
        gene_id=gene_id, seqid=seqid, start=start, end=end, strand=strand,
        confidence=confidence, exons=exons or [],
    )


# --------------------------------------------------------------------------
# load_helixer_loci
# --------------------------------------------------------------------------

def test_load_count_after_merge(load_gff_path):
    loci = load_helixer_loci(load_gff_path)
    # gene1+gene2 merge -> 4 loci total
    assert len(loci) == 4


def test_load_sorted_by_seqid_start(load_gff_path):
    loci = load_helixer_loci(load_gff_path)
    keys = [(g.seqid, g.start) for g in loci]
    assert keys == sorted(keys)


def test_load_merged_locus_coords(load_gff_path):
    loci = load_helixer_loci(load_gff_path)
    merged = next(g for g in loci if g.gene_id == "gene1")
    assert (merged.start, merged.end) == (100, 300)


def test_load_keeps_antisense_separate(load_gff_path):
    loci = load_helixer_loci(load_gff_path)
    ids = {g.gene_id for g in loci}
    # gene4 is antisense to gene1/gene2 and must survive un-merged
    assert "gene4" in ids
    g4 = next(g for g in loci if g.gene_id == "gene4")
    assert g4.strand == "-"
    assert (g4.start, g4.end) == (120, 250)


def test_load_without_h5_confidence_none(load_gff_path):
    loci = load_helixer_loci(load_gff_path)
    assert all(g.confidence is None for g in loci)


def test_load_with_h5_sets_confidence(load_gff_path, load_h5_path):
    loci = load_helixer_loci(load_gff_path, h5_path=load_h5_path)
    for g in loci:
        assert g.confidence is not None
        assert 0.0 <= g.confidence <= 1.0


def test_load_with_h5_confidence_value(load_gff_path, load_h5_path):
    loci = load_helixer_loci(load_gff_path, h5_path=load_h5_path)
    g3 = next(g for g in loci if g.gene_id == "gene3")
    assert g3.confidence == pytest.approx(0.90, abs=1e-5)


def test_load_region_filter_tuple(load_gff_path):
    loci = load_helixer_loci(load_gff_path, region=("chr1", 0, 260))
    ids = {g.gene_id for g in loci}
    # gene1+gene2 (overlap [0,260)) merge; gene4 kept; gene3/gene5 excluded
    assert ids == {"gene1", "gene4"}


def test_load_region_filter_string(load_gff_path):
    loci = load_helixer_loci(load_gff_path, region="chr2:0-300")
    assert [g.gene_id for g in loci] == ["gene5"]


# --------------------------------------------------------------------------
# enrich_with_confidence
# --------------------------------------------------------------------------

def test_enrich_sets_confidence(load_h5_path):
    loci = [_locus("g", "chr1", 100, 300, "+")]
    with HDF5ConfidenceReader(load_h5_path) as reader:
        enriched = enrich_with_confidence(loci, reader)
    assert enriched[0].confidence == pytest.approx(0.90, abs=1e-5)


def test_enrich_missing_scaffold_leaves_none(load_h5_path):
    loci = [_locus("g", "chrZ", 100, 300, "+")]
    with HDF5ConfidenceReader(load_h5_path) as reader:
        enriched = enrich_with_confidence(loci, reader)
    assert enriched[0].confidence is None


def test_enrich_out_of_bounds_leaves_none(load_h5_path):
    # chr2 length is 256; this region runs past the end
    loci = [_locus("g", "chr2", 100, 9999, "+")]
    with HDF5ConfidenceReader(load_h5_path) as reader:
        enriched = enrich_with_confidence(loci, reader)
    assert enriched[0].confidence is None


def test_enrich_preserves_order_and_count(load_h5_path):
    loci = [
        _locus("a", "chr1", 0, 100, "+"),
        _locus("b", "chrZ", 0, 100, "+"),
        _locus("c", "chr2", 0, 100, "-"),
    ]
    with HDF5ConfidenceReader(load_h5_path) as reader:
        enriched = enrich_with_confidence(loci, reader)
    assert [g.gene_id for g in enriched] == ["a", "b", "c"]


# --------------------------------------------------------------------------
# merge_overlapping_loci
# --------------------------------------------------------------------------

def test_merge_same_strand_overlap():
    loci = [_locus("a", "chr1", 100, 200, "+"), _locus("b", "chr1", 150, 300, "+")]
    merged = merge_overlapping_loci(loci)
    assert len(merged) == 1
    assert (merged[0].start, merged[0].end) == (100, 300)


def test_merge_keeps_first_gene_id():
    loci = [_locus("a", "chr1", 100, 200, "+"), _locus("b", "chr1", 150, 300, "+")]
    merged = merge_overlapping_loci(loci)
    assert merged[0].gene_id == "a"


def test_merge_opposite_strand_not_merged():
    loci = [_locus("a", "chr1", 100, 250, "+"), _locus("b", "chr1", 120, 300, "-")]
    merged = merge_overlapping_loci(loci)
    assert len(merged) == 2


def test_merge_touching_not_merged():
    # [100,200) and [200,300) touch but do not overlap
    loci = [_locus("a", "chr1", 100, 200, "+"), _locus("b", "chr1", 200, 300, "+")]
    merged = merge_overlapping_loci(loci)
    assert len(merged) == 2


def test_merge_non_overlapping_separate():
    loci = [_locus("a", "chr1", 100, 200, "+"), _locus("b", "chr1", 400, 500, "+")]
    merged = merge_overlapping_loci(loci)
    assert len(merged) == 2


def test_merge_three_way_chain():
    loci = [
        _locus("a", "chr1", 100, 200, "+"),
        _locus("b", "chr1", 150, 280, "+"),
        _locus("c", "chr1", 260, 400, "+"),
    ]
    merged = merge_overlapping_loci(loci)
    assert len(merged) == 1
    assert (merged[0].start, merged[0].end) == (100, 400)


def test_merge_confidence_mean():
    loci = [
        _locus("a", "chr1", 100, 200, "+", confidence=0.8),
        _locus("b", "chr1", 150, 300, "+", confidence=0.4),
    ]
    merged = merge_overlapping_loci(loci)
    assert merged[0].confidence == pytest.approx(0.6)


def test_merge_confidence_all_none():
    loci = [_locus("a", "chr1", 100, 200, "+"), _locus("b", "chr1", 150, 300, "+")]
    merged = merge_overlapping_loci(loci)
    assert merged[0].confidence is None


def test_merge_union_exons():
    a = _locus("a", "chr1", 100, 200, "+", exons=[Exon(100, 150), Exon(170, 200)])
    b = _locus("b", "chr1", 150, 300, "+", exons=[Exon(180, 220), Exon(260, 300)])
    merged = merge_overlapping_loci([a, b])
    # (170,200) and (180,220) overlap -> (170,220)
    assert [(e.start, e.end) for e in merged[0].exons] == [
        (100, 150), (170, 220), (260, 300)
    ]


def test_merge_drops_cds():
    a = HelixerLocus("a", "chr1", 100, 200, "+", exons=[Exon(100, 200)],
                     cds=None)
    b = HelixerLocus("b", "chr1", 150, 300, "+", exons=[Exon(150, 300)])
    merged = merge_overlapping_loci([a, b])
    assert merged[0].cds is None


def test_merge_provenance_singleton():
    loci = [_locus("a", "chr1", 100, 200, "+")]
    merged, prov = merge_overlapping_loci(loci, return_provenance=True)
    assert prov == {"a": ["a"]}


def test_merge_provenance_records_constituents():
    loci = [
        _locus("a", "chr1", 100, 200, "+"),
        _locus("b", "chr1", 150, 280, "+"),
        _locus("c", "chr1", 260, 400, "+"),
    ]
    merged, prov = merge_overlapping_loci(loci, return_provenance=True)
    assert prov == {"a": ["a", "b", "c"]}


def test_merge_provenance_mixed():
    loci = [
        _locus("a", "chr1", 100, 200, "+"),
        _locus("b", "chr1", 150, 300, "+"),
        _locus("c", "chr1", 500, 600, "+"),
    ]
    merged, prov = merge_overlapping_loci(loci, return_provenance=True)
    assert prov == {"a": ["a", "b"], "c": ["c"]}


def test_merge_preserves_singleton_object():
    single = _locus("a", "chr1", 100, 200, "+", confidence=0.7)
    merged = merge_overlapping_loci([single])
    assert merged[0] is single


# --------------------------------------------------------------------------
# _merge_exon_lists
# --------------------------------------------------------------------------

def test_merge_exon_lists_overlap():
    a = [Exon(100, 200)]
    b = [Exon(150, 250)]
    assert [(e.start, e.end) for e in _merge_exon_lists(a, b)] == [(100, 250)]


def test_merge_exon_lists_touching_not_merged():
    a = [Exon(100, 200)]
    b = [Exon(200, 300)]
    assert [(e.start, e.end) for e in _merge_exon_lists(a, b)] == [(100, 200), (200, 300)]


def test_merge_exon_lists_disjoint():
    a = [Exon(100, 150)]
    b = [Exon(300, 400)]
    assert [(e.start, e.end) for e in _merge_exon_lists(a, b)] == [(100, 150), (300, 400)]


def test_merge_exon_lists_empty():
    assert _merge_exon_lists([], []) == []


def test_merge_exon_lists_nested():
    a = [Exon(100, 400)]
    b = [Exon(150, 250)]
    assert [(e.start, e.end) for e in _merge_exon_lists(a, b)] == [(100, 400)]


# --------------------------------------------------------------------------
# filter_loci_by_region
# --------------------------------------------------------------------------

def test_filter_region_overlap():
    loci = [
        _locus("a", "chr1", 100, 200, "+"),
        _locus("b", "chr1", 400, 500, "+"),
        _locus("c", "chr2", 100, 200, "+"),
    ]
    kept = filter_loci_by_region(loci, "chr1", 0, 250)
    assert [g.gene_id for g in kept] == ["a"]


def test_filter_region_boundary_half_open():
    loci = [_locus("a", "chr1", 100, 200, "+")]
    # query [200,300) does not overlap [100,200)
    assert filter_loci_by_region(loci, "chr1", 200, 300) == []


def test_filter_region_wrong_scaffold():
    loci = [_locus("a", "chr1", 100, 200, "+")]
    assert filter_loci_by_region(loci, "chr2", 0, 1000) == []


# --------------------------------------------------------------------------
# index builders
# --------------------------------------------------------------------------

def test_build_locus_index_query():
    loci = [_locus("a", "chr1", 100, 200, "+"), _locus("b", "chr1", 400, 500, "+")]
    idx = build_locus_index(loci)
    assert idx.query(150, 160) == [0]


def test_build_locus_index_by_scaffold_keys():
    loci = [
        _locus("a", "chr1", 100, 200, "+"),
        _locus("b", "chr2", 100, 200, "+"),
    ]
    indices = build_locus_index_by_scaffold(loci)
    assert set(indices.keys()) == {"chr1", "chr2"}


def test_build_locus_index_by_scaffold_stores_original_indices():
    loci = [
        _locus("a", "chr1", 100, 200, "+"),
        _locus("b", "chr2", 100, 200, "+"),
        _locus("c", "chr1", 400, 500, "+"),
    ]
    indices = build_locus_index_by_scaffold(loci)
    # data payload carries the original-list index (2 for locus c)
    hits = indices["chr1"].query_with_data(420, 430)
    assert hits[0][2] == 2


def test_build_locus_index_by_scaffold_no_cross_scaffold_hits():
    loci = [
        _locus("a", "chr1", 100, 200, "+"),
        _locus("b", "chr2", 100, 200, "+"),
    ]
    indices = build_locus_index_by_scaffold(loci)
    # querying chr1 index must not return the chr2 locus
    assert indices["chr1"].query(100, 200) == [0]
    assert len(indices["chr1"]) == 1
