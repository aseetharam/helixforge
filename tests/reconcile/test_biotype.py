"""Tests for reconcile/biotype.py — Phase 30 D1/D3 (assessment §2.1). Floor: 8.

Turns expressed ORF-less loci into first-class non-coding annotations:
protein_coding (credible ORF) / lncRNA (expressed + spliced + >=200 nt + quiet
CDS channel) / ncRNA_undetermined (ORF-less but short/mono-exonic/low-expr/CDS
still loud) / pseudogene (set upstream, never overridden). Both strands, concrete
literal coordinates.
"""

import pytest

from helixforge.reconcile.biotype import assign_biotype, classify_biotype
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)


def _gene(
    strand,
    *,
    exons,
    cds=None,
    status="EXPRESSED",
    tier=2,
    origin="mikado_1to1",
    biotype=None,
):
    """Build a one-transcript gene from literal exon tuples."""
    start = min(s for s, _ in exons)
    end = max(e for _, e in exons)
    tid = "HFG_00001.1"
    cds_segs = [CDSSegment(*c) for c in cds] if cds else None
    t = TranscriptCandidate(
        tid, "HFG_00001", "stringtie", "chr1", start, end, strand,
        [Exon(s, e) for s, e in exons],
        cds=cds_segs,
        is_primary=True,
    )
    return ReconciledGene(
        "HFG_00001", "chr1", start, end, strand, tier, [t], tid,
        LocusClassification("HFG_00001", status), origin, biotype=biotype,
    )


# --- protein_coding ---------------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_credible_orf_is_protein_coding(strand):
    gene = _gene(strand, exons=[(1000, 1500)], cds=[(1000, 1300, 0)], tier=1)
    out = assign_biotype(gene, cds_channel_conf=0.9)
    assert out.biotype == "protein_coding"
    assert out.transcripts[0].biotype == "protein_coding"
    assert out.tier == 1  # coding tier unchanged


# --- putative coding (good ORF, no evidence) --------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_good_orf_without_evidence_is_putative_coding(strand):
    # A complete ORF with neither expression nor homology is still protein_coding
    # (absence of corroboration is neutral) but flagged PUTATIVE_CODING.
    gene = _gene(
        strand, exons=[(1000, 1500)], cds=[(1000, 1300, 0)],
        status="SILENT", tier=2, origin="helixer_backstop",
    )
    out = assign_biotype(gene)
    assert out.biotype == "protein_coding"
    assert any(f.name == "PUTATIVE_CODING" for f in out.flags)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_expressed_coding_is_not_putative(strand):
    # Expression corroborates the ORF -> coding, no PUTATIVE_CODING flag.
    gene = _gene(
        strand, exons=[(1000, 1500)], cds=[(1000, 1300, 0)],
        status="EXPRESSED", tier=2,
    )
    out = assign_biotype(gene)
    assert out.biotype == "protein_coding"
    assert not any(f.name == "PUTATIVE_CODING" for f in out.flags)


def test_homology_backed_coding_is_not_putative():
    # protein_id (homology) corroborates the ORF even when silent -> not putative.
    start, end = 1000, 1500
    t = TranscriptCandidate(
        "HFG_00001.1", "HFG_00001", "stringtie", "chr1", start, end, "+",
        [Exon(start, end)], cds=[CDSSegment(1000, 1300, 0)],
        protein_id="sp|P12345", is_primary=True,
    )
    gene = ReconciledGene(
        "HFG_00001", "chr1", start, end, "+", 1, [t], "HFG_00001.1",
        LocusClassification("HFG_00001", "SILENT"), "helixer_backstop",
    )
    out = assign_biotype(gene)
    assert out.biotype == "protein_coding"
    assert not any(f.name == "PUTATIVE_CODING" for f in out.flags)


# --- lncRNA -----------------------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_expressed_spliced_long_quiet_is_lncrna(strand):
    # 2 exons -> 1 intron (spliced); 200 + 200 = 400 nt (>= 200); CDS channel quiet.
    gene = _gene(strand, exons=[(1000, 1200), (1400, 1600)], status="EXPRESSED", tier=2)
    out = assign_biotype(gene, cds_channel_conf=0.1)
    assert out.biotype == "lncRNA"
    assert out.transcripts[0].biotype == "lncRNA"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_lncrna_without_cds_channel_signal(strand):
    # No HDF5 signal (None) must not block a lncRNA call (other 4 criteria hold).
    gene = _gene(strand, exons=[(1000, 1200), (1400, 1600)], status="EXPRESSED")
    out = assign_biotype(gene, cds_channel_conf=None)
    assert out.biotype == "lncRNA"


# --- ncRNA_undetermined -----------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_loud_cds_channel_demotes_to_undetermined(strand):
    # Expressed + spliced + long, but Helixer CDS channel still loud -> undetermined.
    gene = _gene(strand, exons=[(1000, 1200), (1400, 1600)], status="EXPRESSED")
    out = assign_biotype(gene, cds_channel_conf=0.8)
    assert out.biotype == "ncRNA_undetermined"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_monoexonic_is_undetermined(strand):
    # Single exon (not spliced) -> not lncRNA even if expressed + long + quiet.
    gene = _gene(strand, exons=[(1000, 1600)], status="EXPRESSED")
    out = assign_biotype(gene, cds_channel_conf=0.1)
    assert out.biotype == "ncRNA_undetermined"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_short_spliced_is_undetermined(strand):
    # Spliced + expressed + quiet but only 50 + 50 = 100 nt (< 200) -> undetermined.
    gene = _gene(strand, exons=[(1000, 1050), (1200, 1250)], status="EXPRESSED")
    out = assign_biotype(gene, cds_channel_conf=0.1)
    assert out.biotype == "ncRNA_undetermined"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_low_expression_spliced_is_undetermined(strand):
    # Spliced + long + quiet but only LOW expression -> undetermined (not lncRNA).
    gene = _gene(strand, exons=[(1000, 1200), (1400, 1600)], status="LOW", tier=3)
    out = assign_biotype(gene, cds_channel_conf=0.1)
    assert out.biotype == "ncRNA_undetermined"


# --- pseudogene precedence --------------------------------------------------

def test_existing_pseudogene_not_overridden():
    gene = _gene("+", exons=[(1000, 1200), (1400, 1600)], biotype="pseudogene", tier=2)
    out = assign_biotype(gene, cds_channel_conf=0.1)
    assert out.biotype == "pseudogene"  # upstream biotype wins
    assert out is gene  # returned unchanged


def test_classify_biotype_pure_decision():
    # classify_biotype returns the string without mutating the gene.
    gene = _gene("+", exons=[(1000, 1200), (1400, 1600)], status="EXPRESSED")
    assert classify_biotype(gene, cds_channel_conf=0.1) == "lncRNA"
    assert gene.biotype is None  # not mutated


# --- D3 tier coherence ------------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_noncoding_demotes_out_of_cds_tier(strand):
    # A Mikado-origin ORF-less gene sits at Tier 2; as lncRNA it must leave the
    # CDS+homology tier -> EXPRESSED => Tier 3.
    gene = _gene(strand, exons=[(1000, 1200), (1400, 1600)], status="EXPRESSED", tier=2)
    out = assign_biotype(gene, cds_channel_conf=0.1)
    assert out.biotype == "lncRNA"
    assert out.tier == 3


def test_silent_noncoding_tiers_to_4():
    # SILENT ORF-less mono-exonic gene at Tier 2 -> ncRNA_undetermined, Tier 4.
    gene = _gene("+", exons=[(1000, 1600)], status="SILENT", tier=2)
    out = assign_biotype(gene)
    assert out.biotype == "ncRNA_undetermined"
    assert out.tier == 4


def test_backstop_tier3_preserved():
    # A CDS-less backstop already at Tier 3 keeps its tier when typed lncRNA.
    gene = _gene(
        "+", exons=[(1000, 1200), (1400, 1600)], status="EXPRESSED",
        tier=3, origin="helixer_backstop",
    )
    out = assign_biotype(gene, cds_channel_conf=0.1)
    assert out.biotype == "lncRNA"
    assert out.tier == 3
