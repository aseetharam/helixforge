"""Phase 24 D4 — boundary-stitch pass. Floor: 3.

A cross-chunk merge split by a partition boundary is recovered when a verified
junction bridges the gap and both genes sit within ``flank`` of the boundary;
otherwise the pass is a no-op. Both strands; concrete literal coordinates
(CLAUDE.md §12). The bridging predicate is the reconciler's own
``_has_bridging_introns``.
"""

import pytest

from helixforge.parallel.plan import Chunk, Plan
from helixforge.parallel.stitch import boundary_stitch
from helixforge.qc.flags import LOCUS_MERGE
from helixforge.reconcile.mikado_integrate import IdAllocator, build_reconciled_gene
from helixforge.reconcile.models import Exon, HelixerLocus, LocusClassification, SpliceJunction


def _backstop(seed, seqid, start, end, strand, id_map, alloc):
    h = HelixerLocus(gene_id=seed, seqid=seqid, start=start, end=end,
                     strand=strand, exons=[Exon(start, end)])
    cls = LocusClassification(seed, "SILENT")
    return build_reconciled_gene(h, None, cls, "helixer_backstop", id_map,
                                 allocator=alloc)


def _two_genes(strand):
    """Genes A (2300-2400) and B (2600-2700) straddling a boundary at 2500."""
    id_map, alloc = {}, IdAllocator()
    a = _backstop("hA", "chr1", 2300, 2400, strand, id_map, alloc)
    b = _backstop("hB", "chr1", 2600, 2700, strand, id_map, alloc)
    return [a, b]


def _plan(flank=200):
    # Two chr1 chunks cut at 2500 → a single interior boundary there.
    return Plan(
        chunks=[
            Chunk("chunk_0000", ["chr1:0-2500"], 1),
            Chunk("chunk_0001", ["chr1:2500-10000"], 1),
        ],
        flank=flank,
    )


# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def test_cross_boundary_merge_recovered_plus():
    genes = _two_genes("+")
    # A junction spanning the gap [2400, 2600): donor <= 2400, acceptor >= 2600.
    junctions = [SpliceJunction("chr1", 2350, 2650, "+", read_count=12)]
    stitched, recovered = boundary_stitch(genes, _plan(), junctions)
    assert recovered == 1
    assert len(stitched) == 1
    merged = stitched[0]
    # Lowest HFG kept (A), the other recorded in merged_from, span unioned.
    assert merged.gene_id == genes[0].gene_id
    assert genes[1].gene_id in merged.merged_from
    assert (merged.start, merged.end) == (2300, 2700)
    assert len(merged.transcripts) == 2
    assert LOCUS_MERGE in merged.flags


def test_minus_strand_merge_recovered():
    genes = _two_genes("-")
    junctions = [SpliceJunction("chr1", 2350, 2650, "-", read_count=9)]
    stitched, recovered = boundary_stitch(genes, _plan(), junctions)
    assert recovered == 1
    assert stitched[0].strand == "-"
    assert (stitched[0].start, stitched[0].end) == (2300, 2700)


# ---------------------------------------------------------------------------
# No-ops
# ---------------------------------------------------------------------------

def test_no_bridging_junction_is_noop():
    genes = _two_genes("+")
    # Junction does NOT span the full gap (acceptor 2450 < right gene start 2600).
    junctions = [SpliceJunction("chr1", 2350, 2450, "+", read_count=12)]
    stitched, recovered = boundary_stitch(genes, _plan(), junctions)
    assert recovered == 0
    assert {g.gene_id for g in stitched} == {g.gene_id for g in genes}


def test_strand_mismatch_not_merged():
    genes = _two_genes("-")               # genes on minus strand
    junctions = [SpliceJunction("chr1", 2350, 2650, "+", read_count=12)]  # plus
    _stitched, recovered = boundary_stitch(genes, _plan(), junctions)
    assert recovered == 0


def test_genes_beyond_flank_not_merged():
    genes = _two_genes("+")
    junctions = [SpliceJunction("chr1", 2350, 2650, "+", read_count=12)]
    # flank=50: gene A ends at 2400 (boundary 2500) → 100 bp away > flank → skip.
    _stitched, recovered = boundary_stitch(genes, _plan(flank=50), junctions)
    assert recovered == 0
