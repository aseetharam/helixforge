"""Unit tests for ``helixforge.stats.summary.summarize_genes`` (Phase 12 prep).

Synthetic genes only, concrete literal coordinates, both strands (CLAUDE.md
§12). This is the fast, data-independent counterpart to the heavy
``tests/test_m3_golden.py`` gate: it proves the count helper itself is correct
regardless of whether the real M3 data is present.
"""

from __future__ import annotations

from helixforge.qc.flags import HELIXER_ONLY, NO_STOP
from helixforge.reconcile.models import (
    ASEvent,
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.stats.summary import summarize_genes


def _classification(locus_id, status="EXPRESSED"):
    return LocusClassification(
        locus_id=locus_id,
        status=status,
        max_tpm=10.0,
        mean_coverage=5.0,
        evidence_source="stringtie",
        num_samples_expressed=1,
    )


def _plus_coding_gene():
    """Plus-strand gene, 2 transcripts (multi-isoform), CDS present, 1 AS event."""
    exons = [Exon(1000, 1200), Exon(1300, 1500), Exon(1600, 1800)]
    cds = [CDSSegment(1050, 1200, 0), CDSSegment(1300, 1500, 0), CDSSegment(1600, 1649, 1)]
    t1 = TranscriptCandidate(
        transcript_id="HFG_00001.1", locus_id="HFG_00001", source="mikado",
        seqid="chr1", start=1000, end=1800, strand="+", exons=exons, cds=cds,
        tpm=10.0, junction_support_fraction=1.0, confidence=0.9,
        combined_score=15.0, is_primary=True,
    )
    t2 = TranscriptCandidate(
        transcript_id="HFG_00001.2", locus_id="HFG_00001", source="mikado",
        seqid="chr1", start=1000, end=1800, strand="+",
        exons=[Exon(1000, 1200), Exon(1600, 1800)], cds=None,
        tpm=3.0, junction_support_fraction=1.0, confidence=0.7,
        combined_score=8.0, is_primary=False,
    )
    return ReconciledGene(
        gene_id="HFG_00001", seqid="chr1", start=1000, end=1800, strand="+",
        tier=1, transcripts=[t1, t2], primary_transcript_id="HFG_00001.1",
        classification=_classification("HFG_00001"), origin="mikado_1to1",
        as_events=[ASEvent("ES", "chr1", 1300, 1500, "+", support_read_count=20)],
        flags=[],
    )


def _minus_backstop_gene():
    """Minus-strand backstop gene, single transcript, no CDS, two flags."""
    exons = [Exon(2000, 2200), Exon(2300, 2500), Exon(2600, 2800)]
    t1 = TranscriptCandidate(
        transcript_id="HFG_00002.1", locus_id="HFG_00002", source="helixer",
        seqid="chr1", start=2000, end=2800, strand="-", exons=exons, cds=None,
        tpm=0.0, junction_support_fraction=0.0, confidence=0.5,
        combined_score=1.0, is_primary=True,
    )
    return ReconciledGene(
        gene_id="HFG_00002", seqid="chr1", start=2000, end=2800, strand="-",
        tier=3, transcripts=[t1], primary_transcript_id="HFG_00002.1",
        classification=_classification("HFG_00002", status="SILENT"),
        origin="helixer_backstop", as_events=[],
        flags=[HELIXER_ONLY, NO_STOP],
    )


def test_summarize_empty():
    s = summarize_genes([])
    assert s == {
        "genes": 0, "total_tx": 0, "multi": 0, "coding": 0,
        "tier": {}, "origin": {}, "biotype": {}, "as_kinds": {}, "flags": {},
    }


def test_summarize_both_strands():
    s = summarize_genes([_plus_coding_gene(), _minus_backstop_gene()])
    assert s["genes"] == 2
    assert s["total_tx"] == 3          # 2 + 1
    assert s["multi"] == 1             # only the plus gene has >1 transcript
    assert s["coding"] == 1            # only the plus gene has a CDS
    assert s["tier"] == {1: 1, 3: 1}
    assert s["origin"] == {"helixer_backstop": 1, "mikado_1to1": 1}
    assert s["as_kinds"] == {"ES": 1}
    assert s["flags"] == {"HELIXER_ONLY": 1, "NO_STOP": 1}


def test_summarize_count_dicts_are_sorted():
    """Nested count dicts have deterministic (sorted) key order."""
    s = summarize_genes([_minus_backstop_gene(), _plus_coding_gene()])
    assert list(s["flags"].keys()) == ["HELIXER_ONLY", "NO_STOP"]
    assert list(s["tier"].keys()) == [1, 3]
