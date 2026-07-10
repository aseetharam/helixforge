"""Tests for reconcile/pseudogene.py — Phase 29 D2/D3 (assessment §2.3). Floor: 4.

A homology-backed disabled ORF (premature stop / mod-3 frameshift) is a
pseudogene candidate, not a failed coding gene. A broken ORF *without* homology
stays a normal flagged gene; a clean gene is untouched. Both strands. Concrete
literal coordinates.
"""

import pytest

from helixforge.qc.flags import INTERNAL_STOP, NO_STOP
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.reconcile.pseudogene import (
    apply_pseudogene_typing,
    is_pseudogene_candidate,
)


def _coding_gene(strand, *, protein=True, blast=None, tier=1, biotype=None,
                 cds=((1000, 1300, 0),)):
    """A coding primary transcript: exon (1000,1500), CDS 300 bp (mod-3)."""
    tid = "HFG_00001.1"
    t = TranscriptCandidate(
        tid, "HFG_00001", "mikado", "chr1", 1000, 1500, strand,
        [Exon(1000, 1500)],
        cds=[CDSSegment(*c) for c in cds],
        protein_id=("P1" if protein else None),
        blast_score=blast,
        is_primary=True,
    )
    return ReconciledGene(
        "HFG_00001", "chr1", 1000, 1500, strand, tier, [t], tid,
        LocusClassification("HFG_00001", "EXPRESSED"), "mikado_1to1",
        biotype=biotype,
    )


def _names(gene):
    return {f.name for f in gene.flags}


@pytest.mark.parametrize("strand", ["+", "-"])
def test_homology_plus_internal_stop_is_pseudogene(strand):
    gene = _coding_gene(strand)  # protein-backed, Tier 1
    out = apply_pseudogene_typing(gene, [INTERNAL_STOP])
    assert out.biotype == "pseudogene"
    assert "PSEUDOGENE_CANDIDATE" in _names(out)
    assert out.tier == 2  # demoted out of Tier 1


def test_no_homology_broken_orf_stays_normal():
    gene = _coding_gene("+", protein=False, blast=None)
    out = apply_pseudogene_typing(gene, [INTERNAL_STOP])
    assert out is gene  # unchanged identity
    assert out.biotype is None
    assert "PSEUDOGENE_CANDIDATE" not in _names(out)


def test_clean_homology_gene_unaffected():
    gene = _coding_gene("+")
    out = apply_pseudogene_typing(gene, [])  # no lesion
    assert out is gene
    assert out.biotype is None


def test_partial_orf_no_stop_is_not_a_lesion():
    # A 5'/3'-partial ORF (NO_STOP from a truncated terminus) is truncated, not
    # disabled — it must NOT be typed pseudogene on a NO_STOP flag alone.
    gene = _coding_gene("+")
    assert not is_pseudogene_candidate(gene, [NO_STOP])
    out = apply_pseudogene_typing(gene, [NO_STOP])
    assert out.biotype is None


def test_blast_score_threshold_gates_pseudogene():
    # blast-only homology: below the floor → not a pseudogene; at/above → yes.
    weak = _coding_gene("+", protein=False, blast=2.0)
    assert not is_pseudogene_candidate(weak, [INTERNAL_STOP], min_homology=5.0)
    strong = _coding_gene("+", protein=False, blast=9.0)
    assert is_pseudogene_candidate(strong, [INTERNAL_STOP], min_homology=5.0)


def test_existing_biotype_not_overridden():
    gene = _coding_gene("+", biotype="protein_coding")
    out = apply_pseudogene_typing(gene, [INTERNAL_STOP])
    assert out is gene
    assert out.biotype == "protein_coding"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_pseudogene_tier2_stays_tier2(strand):
    # A Tier-2 (CDS-only) homology gene with a lesion is still typed pseudogene;
    # tier demotion only fires for Tier 1, so Tier 2 is left as-is.
    gene = _coding_gene(strand, tier=2)
    out = apply_pseudogene_typing(gene, [INTERNAL_STOP])
    assert out.biotype == "pseudogene"
    assert out.tier == 2
