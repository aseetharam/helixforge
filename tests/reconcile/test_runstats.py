"""Phase 23 D3 — RunStats decision telemetry. Floor: 4 (both strands).

The counters must match hand-computed values on a synthetic reconcile, and
threading the stats object must not change any gene/tier/origin count (it only
observes). Both strands are exercised for the merge/split/backstop paths.
"""

from helixforge.qc.flags import HELIXER_ONLY, NO_EXPRESSION, PARTIAL_ORF
from helixforge.reconcile.cds import assign_backstop_cds
from helixforge.reconcile.fallback import refine_backstop_gene
from helixforge.reconcile.mikado_integrate import reconcile
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    HelixerLocus,
    LocusClassification,
    MikadoLocus,
    MiniprotAlignment,
    ReconciledGene,
    SpliceJunction,
    TranscriptCandidate,
)
from helixforge.reconcile.runstats import RunStats


# --- builders ---------------------------------------------------------------

def _h(gid, seqid, start, end, strand, exons):
    return HelixerLocus(gid, seqid, start, end, strand, exons=exons)


def _mt(tid, seqid, strand, exons, cds=None, score=10.0):
    return TranscriptCandidate(
        tid, "L", "mikado", seqid, min(e.start for e in exons),
        max(e.end for e in exons), strand, exons, cds=cds, combined_score=score,
    )


def _m(lid, seqid, start, end, strand, transcripts):
    return MikadoLocus(lid, seqid, start, end, strand, transcripts)


# --- basic counter mechanics ------------------------------------------------

def test_bump_headline_and_as_dict():
    s = RunStats()
    s.bump("merges_accepted")
    s.bump("splits", 3)
    h = s.headline()
    assert h["merges_accepted"] == 1 and h["splits"] == 3
    assert "tier_counts" not in h            # headline is the flat numbers only
    assert s.as_dict()["splits"] == 3
    assert s.as_dict()["tier_counts"] == {}  # derived, empty until populated


def test_populate_from_genes_tier_origin_partial():
    def gene(gid, tier, origin, flags):
        tx = TranscriptCandidate(f"{gid}.1", gid, "helixer", "chr1", 100, 200,
                                 "+", [Exon(100, 200)], is_primary=True)
        return ReconciledGene(gid, "chr1", 100, 200, "+", tier, [tx], f"{gid}.1",
                              LocusClassification(gid, "EXPRESSED"), origin,
                              flags=list(flags))

    genes = [
        gene("HFG_00001", 1, "mikado_1to1", []),
        gene("HFG_00002", 1, "mikado_1to1", [PARTIAL_ORF]),
        gene("HFG_00003", 3, "helixer_backstop", [HELIXER_ONLY]),
    ]
    s = RunStats()
    s.populate_from_genes(genes)
    assert s.tier_counts == {1: 2, 3: 1}
    assert s.origin_counts == {"helixer_backstop": 1, "mikado_1to1": 2}
    assert s.partial_orfs == 1
    # idempotent
    s.populate_from_genes(genes)
    assert s.tier_counts == {1: 2, 3: 1}


# --- reconcile decision counts (both strands) -------------------------------

def _merge_split_scenario(strand):
    """Build a set that yields 1 accepted merge, 1 rejected merge, 1 split."""
    # Merge-accepted: two Helixer loci bridged by a Mikado intron (1500..2000).
    ha = _h("Ha", "chr1", 1000, 1500, strand, [Exon(1000, 1500)])
    hb = _h("Hb", "chr1", 2000, 2500, strand, [Exon(2000, 2500)])
    m_accept = _m("Macc", "chr1", 1000, 2500, strand, [
        _mt("Macc.1", "chr1", strand, [Exon(1000, 1500), Exon(2000, 2500)],
            cds=[CDSSegment(1000, 1300, 0), CDSSegment(2000, 2300, 0)])])
    # Merge-rejected: two loci over a single-exon Mikado locus (no bridge).
    hc = _h("Hc", "chr2", 5000, 5500, strand, [Exon(5000, 5500)])
    hd = _h("Hd", "chr2", 6000, 6500, strand, [Exon(6000, 6500)])
    m_reject = _m("Mrej", "chr2", 5000, 6500, strand, [
        _mt("Mrej.1", "chr2", strand, [Exon(5000, 6500)],
            cds=[CDSSegment(5000, 6200, 0)])])
    # Split: one Helixer locus over two Mikado loci.
    he = _h("He", "chr3", 3000, 5000, strand, [Exon(3000, 3500), Exon(4500, 5000)])
    m_s1 = _m("Ms1", "chr3", 3000, 3500, strand, [
        _mt("Ms1.1", "chr3", strand, [Exon(3000, 3500)],
            cds=[CDSSegment(3000, 3300, 0)])])
    m_s2 = _m("Ms2", "chr3", 4500, 5000, strand, [
        _mt("Ms2.1", "chr3", strand, [Exon(4500, 5000)],
            cds=[CDSSegment(4500, 4800, 0)])])
    helixer = [ha, hb, hc, hd, he]
    mikado = [m_accept, m_reject, m_s1, m_s2]
    cls = [LocusClassification(h.gene_id, "EXPRESSED") for h in helixer]
    return helixer, cls, mikado


def _assert_merge_split_counts(strand):
    helixer, cls, mikado = _merge_split_scenario(strand)
    stats = RunStats()
    genes, _id, _adm = reconcile(helixer, cls, mikado, stats=stats)
    assert stats.merges_accepted == 1
    assert stats.merges_rejected == 1
    assert stats.splits == 1
    # observation only: a no-stats run yields the identical gene set.
    genes_b, _, _ = reconcile(helixer, cls, mikado)
    assert [g.gene_id for g in genes] == [g.gene_id for g in genes_b]
    assert [g.origin for g in genes] == [g.origin for g in genes_b]


def test_reconcile_counts_plus_strand():
    _assert_merge_split_counts("+")


def test_reconcile_counts_minus_strand():
    _assert_merge_split_counts("-")


def test_reconcile_counts_dropped_redundant():
    # Two identical isoforms in one Mikado locus → one dropped as redundant.
    exons = [Exon(1000, 1200), Exon(1500, 2000)]
    cds = [CDSSegment(1000, 1200, 0), CDSSegment(1500, 1900, 0)]
    h = _h("H1", "chr1", 1000, 2000, "+", exons)
    m = _m("M1", "chr1", 1000, 2000, "+", [
        _mt("M1.1", "chr1", "+", exons, cds=cds, score=10.0),
        _mt("M1.2", "chr1", "+", exons, cds=cds, score=8.0),
    ])
    stats = RunStats()
    reconcile([h], [LocusClassification("H1", "EXPRESSED")], [m], stats=stats)
    assert stats.isoforms_dropped_redundant == 1


# --- backstop junction correction counts (both strands) ---------------------

def _backstop_gene(exon_bounds, strand):
    tx = TranscriptCandidate(
        "HFG_00001.1", "HFG_00001", "helixer", "chr1",
        exon_bounds[0][0], exon_bounds[-1][1], strand,
        [Exon(*b) for b in exon_bounds], is_primary=True,
    )
    return ReconciledGene(
        "HFG_00001", "chr1", exon_bounds[0][0], exon_bounds[-1][1], strand, 4,
        [tx], "HFG_00001.1", LocusClassification("HFG_00001", "SILENT"),
        "helixer_backstop", flags=[HELIXER_ONLY, NO_EXPRESSION],
    )


def _jn(donor, acceptor, strand, reads=10):
    return SpliceJunction("chr1", donor, acceptor, strand, read_count=reads)


def test_refine_backstop_counts_applied_and_reverted():
    bounds = [(1000, 1200), (1300, 1500), (1600, 1800)]
    for strand in ("+", "-"):
        # intron0 contradicted by a junction that lands a valid structure → applied.
        applied_js = [_jn(1200, 1320, strand), _jn(1500, 1600, strand)]
        s = RunStats()
        refine_backstop_gene(_backstop_gene(bounds, strand), applied_js, stats=s)
        assert s.junction_corrections_applied == 1
        assert s.junction_corrections_reverted == 0

        # contradicting junction would make a <20 bp intron → transaction reverts.
        revert_js = [_jn(1200, 1215, strand), _jn(1500, 1600, strand)]
        s2 = RunStats()
        refine_backstop_gene(_backstop_gene(bounds, strand), revert_js, stats=s2)
        assert s2.junction_corrections_applied == 0
        assert s2.junction_corrections_reverted == 1


def test_assign_backstop_cds_counts_rescue_source():
    for strand in ("+", "-"):
        gene = _backstop_gene([(1000, 1200), (1300, 1500)], strand)
        aln = MiniprotAlignment(
            protein_id="sp|TEST", seqid="chr1", start=1000, end=1500, strand=strand,
            cds_segments=[CDSSegment(1000, 1150, 0), CDSSegment(1300, 1450, 0)],
            query_coverage=0.9, identity=0.85, score=500.0, rank=0,
        )
        s = RunStats()
        assign_backstop_cds(gene, [aln], stats=s)
        assert s.backstop_rescued_miniprot == 1
        assert s.backstop_rescued_none == 0

        # No alignment → no rescue source available → counted as "none".
        s2 = RunStats()
        assign_backstop_cds(gene, [], stats=s2)
        assert s2.backstop_rescued_none == 1
        assert s2.backstop_rescued_miniprot == 0
