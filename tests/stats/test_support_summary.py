"""Phase 31 D4 — evidence-support / mapping-rate summary.

Hand-computed fractions only (CLAUDE.md §12); both strands. Mapping-rate
arithmetic is checked from explicit counts; the BAM-stats reader is checked
against the fully-mapped synthetic BAM fixture.
"""

from __future__ import annotations

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.stats.report import support_summary


def _gene(gid, strand="+", tpm=5.0, jsf=1.0, cds=True):
    exons = [Exon(1000, 1200), Exon(1300, 1500)]
    cdsseg = [CDSSegment(1050, 1200, 0), CDSSegment(1300, 1450, 0)] if cds else None
    tx = TranscriptCandidate(
        transcript_id=f"{gid}.1", locus_id=gid, source="mikado", seqid="chr1",
        start=1000, end=1500, strand=strand, exons=exons, cds=cdsseg,
        tpm=tpm, junction_support_fraction=jsf, confidence=0.9,
        combined_score=10.0, is_primary=True,
    )
    return ReconciledGene(
        gene_id=gid, seqid="chr1", start=1000, end=1500, strand=strand, tier=1,
        transcripts=[tx], primary_transcript_id=f"{gid}.1",
        classification=LocusClassification(locus_id=gid, status="EXPRESSED", max_tpm=tpm),
        origin="mikado_1to1",
    )


def test_support_fractions_hand_computed():
    # 4 genes (both strands): junction support {1.0, 0.5, 0.0, None→0}; tpm {5,0.4,2,0}.
    genes = [
        _gene("HFG_00001", strand="+", tpm=5.0, jsf=1.0),
        _gene("HFG_00002", strand="-", tpm=0.4, jsf=0.5),
        _gene("HFG_00003", strand="+", tpm=2.0, jsf=0.0),
        _gene("HFG_00004", strand="-", tpm=0.0, jsf=0.0),
    ]
    ss = support_summary(genes, tpm_threshold=0.5, bam_stats={"mapped": 75, "unmapped": 25})
    # junction-supported: genes 1 (1.0) and 2 (0.5) → 2/4.
    assert ss["genes_junction_supported"] == 2
    assert ss["frac_genes_junction_supported"] == pytest.approx(0.5)
    # tpm >= 0.5: genes 1 (5.0) and 3 (2.0) → 2/4.
    assert ss["genes_tpm_pass"] == 2
    assert ss["frac_genes_tpm_pass"] == pytest.approx(0.5)
    # mapping rate 75 / (75+25) = 0.75.
    assert ss["mapping_rate"] == pytest.approx(0.75)


def test_support_mapping_rate_precomputed_and_none():
    genes = [_gene("HFG_00001")]
    # explicit mapping_rate wins over counts
    ss = support_summary(genes, bam_stats={"mapping_rate": 0.42})
    assert ss["mapping_rate"] == pytest.approx(0.42)
    # no bam stats → None (not an error)
    assert support_summary(genes)["mapping_rate"] is None


def test_support_aed_distribution():
    genes = [_gene("HFG_00001", tpm=5.0, jsf=1.0), _gene("HFG_00002", tpm=0.0, jsf=0.0)]
    ss = support_summary(genes)
    # gene1: mean(1.0, 1.0)=1 → AED 0; gene2: mean(0.0, 0.0)=0 → AED 1.
    assert ss["aed"]["n"] == 2
    assert ss["aed"]["min"] == pytest.approx(0.0)
    assert ss["aed"]["max"] == pytest.approx(1.0)
    assert ss["aed"]["mean"] == pytest.approx(0.5)
