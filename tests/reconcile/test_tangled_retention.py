"""Phase 18 D3 — retain Mikado loci in many-to-many components.

In a tangled tandem array, a Mikado locus can lose all its Helixer claimants
during many-to-many decomposition. It must NOT be silently dropped: it is
retained as a novel-style gene (flagged NOVEL_LOCUS + FROM_TANGLED_LOCUS) when
``admit_novel``, or recorded as a not-admitted admission otherwise. No Helixer
gene is ever lost.

Concrete literal coordinates only (CLAUDE.md §12).
"""

from helixforge.reconcile.mikado_integrate import (
    Correspondence,
    _decompose_many_to_many,
    match_loci,
    reconcile,
)
from helixforge.reconcile.models import (
    Exon,
    HelixerLocus,
    LocusClassification,
    MikadoLocus,
    TranscriptCandidate,
)


def helixer(gid, start, end, strand="+"):
    return HelixerLocus(
        gene_id=gid, seqid="chr1", start=start, end=end, strand=strand,
        exons=[Exon(start, end)],
    )


def mikado(lid, start, end, strand="+"):
    tx = TranscriptCandidate(
        transcript_id=f"{lid}.1", locus_id=lid, source="mikado", seqid="chr1",
        start=start, end=end, strand=strand, exons=[Exon(start, end)],
        combined_score=10.0,
    )
    return MikadoLocus(
        locus_id=lid, seqid="chr1", start=start, end=end, strand=strand,
        transcripts=[tx],
    )


# A 2-Helixer × 2-Mikado tangle where both Helixer loci prefer m1, leaving m2
# connected (via h1) but claimant-less.
#   h1 [1000,2000)  overlaps m1 fully (1.0) and m2 partially (0.43)
#   h2 [1050,1950)  overlaps m1 (0.94); below-threshold to m2
#   m1 [1000,1900)  the claimed locus
#   m2 [1850,2200)  loses all claimants
def tangle():
    return (
        [helixer("AT1", 1000, 2000), helixer("AT2", 1050, 1950)],
        [mikado("m1", 1000, 1900), mikado("m2", 1850, 2200)],
    )


def classifications(*gids):
    return [LocusClassification(g, "EXPRESSED") for g in gids]


# --- match_loci / _decompose_many_to_many: retained, not dropped ---

def test_decompose_retains_claimantless_locus():
    corr = Correspondence()
    hs = [helixer("AT1", 1000, 2000), helixer("AT2", 1050, 1950)]
    ms = [mikado("m1", 1000, 1900), mikado("m2", 1850, 2200)]
    _decompose_many_to_many(corr, hs, ms)
    assert any(m.locus_id == "m2" for m in corr.novel)
    assert "m2" in corr.tangled_loci
    # the claimed locus becomes a merge (both Helixer kept)
    assert corr.merges and corr.merges[0][1].locus_id == "m1"


def test_match_loci_marks_tangled_locus_novel():
    hloci, mloci = tangle()
    corr = match_loci(hloci, mloci, reciprocal_overlap=0.3)
    assert {m.locus_id for m in corr.novel} == {"m2"}
    assert corr.tangled_loci == {"m2"}


# --- reconcile, admit_novel=True: emitted with both flags, no Helixer lost ---

def test_reconcile_admit_novel_emits_tangled_gene():
    hloci, mloci = tangle()
    genes, id_map, admissions = reconcile(
        hloci, classifications("AT1", "AT2"), mloci,
        admit_novel=True, reciprocal_overlap=0.3,
    )
    novel_genes = [g for g in genes if g.origin == "novel"]
    assert len(novel_genes) == 1
    flag_names = {f.name for f in novel_genes[0].flags}
    assert "NOVEL_LOCUS" in flag_names
    assert "FROM_TANGLED_LOCUS" in flag_names
    # both Helixer loci retained (merge keeps both ids in the map)
    assert "AT1" in id_map and "AT2" in id_map


# --- reconcile, admit_novel=False: audited, not dropped, no Helixer lost ---

def test_reconcile_no_admit_audits_tangled_locus():
    hloci, mloci = tangle()
    genes, id_map, admissions = reconcile(
        hloci, classifications("AT1", "AT2"), mloci,
        admit_novel=False, reciprocal_overlap=0.3,
    )
    # m2 is not emitted as a gene...
    assert all(g.origin != "novel" for g in genes)
    # ...but it left an audit record (not silently dropped).
    m2_records = [a for a in admissions if a.gene_id == "m2"]
    assert m2_records and all(a.admitted is False for a in m2_records)
    # no Helixer gene lost
    assert "AT1" in id_map and "AT2" in id_map


def test_reconcile_no_helixer_lost_in_tangle():
    hloci, mloci = tangle()
    for admit in (True, False):
        _, id_map, _ = reconcile(
            hloci, classifications("AT1", "AT2"), mloci,
            admit_novel=admit, reciprocal_overlap=0.3,
        )
        assert {"AT1", "AT2"} <= set(id_map)
