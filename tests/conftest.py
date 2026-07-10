"""Shared test fixtures (Phase 0).

Concrete literal coordinates only (CLAUDE.md §12). Synthetic fixtures only —
no real biological data. Both strands are represented.
"""

import pytest

from helixforge.reconcile.models import (
    ASEvent,
    CDSSegment,
    Exon,
    HelixerLocus,
    LocusClassification,
    MikadoLocus,
    MiniprotAlignment,
    ReconciledGene,
    StringTieTranscript,
    TranscriptCandidate,
)
from helixforge.utils.sequences import reverse_complement


class MockGenome:
    """A tiny in-memory genome. ``get_sequence`` reverse-complements on minus."""

    def __init__(self, sequences):
        self.sequences = dict(sequences)

    def get_sequence(self, seqid, start, end, strand="+"):
        seq = self.sequences[seqid][start:end]
        if strand == "-":
            return reverse_complement(seq)
        return seq


# A 60 bp contig holding a plus-strand ORF: ATG AAA CCC GGG TTT TAA at [0, 18).
_CHR1 = "ATGAAACCCGGGTTTTAA" + "A" * 42


@pytest.fixture
def mock_genome():
    return MockGenome({"chr1": _CHR1})


@pytest.fixture
def sample_exons_plus():
    """Three plus-strand exons; introns at 1200-1300 and 1500-1600."""
    return [Exon(1000, 1200), Exon(1300, 1500), Exon(1600, 1800)]


@pytest.fixture
def sample_exons_minus():
    """Three minus-strand exons in a distinct region."""
    return [Exon(2000, 2200), Exon(2300, 2500), Exon(2600, 2800)]


@pytest.fixture
def sample_cds_segments():
    """CDS within ``sample_exons_plus``; total length 399 (divisible by 3)."""
    return [
        CDSSegment(1050, 1200, 0),
        CDSSegment(1300, 1500, 0),
        CDSSegment(1600, 1649, 1),
    ]


@pytest.fixture
def sample_helixer_locus(sample_exons_plus, sample_cds_segments):
    return HelixerLocus(
        gene_id="HELIXER_chr1_1",
        seqid="chr1",
        start=1000,
        end=1800,
        strand="+",
        confidence=0.9,
        exons=sample_exons_plus,
        cds=sample_cds_segments,
    )


@pytest.fixture
def sample_stringtie_transcript(sample_exons_plus):
    return StringTieTranscript(
        transcript_id="ST.1",
        gene_id="STRG.1",
        seqid="chr1",
        start=1000,
        end=1800,
        strand="+",
        exons=sample_exons_plus,
        tpm=12.5,
        sample_id="sample_a",
        coverage=8.0,
    )


@pytest.fixture
def sample_miniprot_alignment(sample_cds_segments):
    return MiniprotAlignment(
        protein_id="sp|P12345",
        seqid="chr1",
        start=1000,
        end=1800,
        strand="+",
        cds_segments=sample_cds_segments,
        query_coverage=0.95,
        identity=0.88,
        score=500.0,
        rank=0,
    )


@pytest.fixture
def sample_transcript_candidate(sample_exons_plus, sample_cds_segments):
    return TranscriptCandidate(
        transcript_id="HFG_00001.1",
        locus_id="HFG_00001",
        source="mikado",
        seqid="chr1",
        start=1000,
        end=1800,
        strand="+",
        exons=sample_exons_plus,
        cds=sample_cds_segments,
        tpm=12.5,
        junction_support_fraction=1.0,
        confidence=0.9,
        combined_score=15.0,
        is_primary=True,
    )


@pytest.fixture
def sample_as_event():
    return ASEvent("ES", "chr1", 1300, 1500, "+", support_read_count=20)


@pytest.fixture
def sample_classification():
    return LocusClassification(
        locus_id="HFG_00001",
        status="EXPRESSED",
        max_tpm=12.5,
        mean_coverage=8.0,
        evidence_source="stringtie",
        num_samples_expressed=1,
    )


@pytest.fixture
def sample_mikado_locus(sample_transcript_candidate):
    return MikadoLocus(
        locus_id="mikado.1",
        seqid="chr1",
        start=1000,
        end=1800,
        strand="+",
        transcripts=[sample_transcript_candidate],
        metrics={"cdna_length": 600},
        scores={"score": 15.0},
    )


@pytest.fixture
def sample_reconciled_gene(
    sample_transcript_candidate, sample_classification, sample_as_event
):
    return ReconciledGene(
        gene_id="HFG_00001",
        seqid="chr1",
        start=1000,
        end=1800,
        strand="+",
        tier=1,
        transcripts=[sample_transcript_candidate],
        primary_transcript_id="HFG_00001.1",
        classification=sample_classification,
        origin="mikado_1to1",
        as_events=[sample_as_event],
        flags=[],
        merged_from=[],
    )
