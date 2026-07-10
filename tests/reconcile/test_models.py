"""Tests for all data models (Phase 0). Floor: 75 tests.

Concrete literal coordinates only. Both strands exercised throughout.
"""

import attrs
import pytest

from helixforge.reconcile.models import (
    ASEvent,
    CDSSegment,
    Exon,
    HelixerLocus,
    Interval,
    IsoformAdmission,
    LocusClassification,
    MikadoLocus,
    MiniprotAlignment,
    QCFlag,
    ReconciledGene,
    SpliceJunction,
    StringTieTranscript,
    TranscriptCandidate,
)


# --------------------------------------------------------------------------
# Interval
# --------------------------------------------------------------------------

def test_interval_valid():
    iv = Interval(1000, 1200)
    assert iv.start == 1000
    assert iv.end == 1200


def test_interval_len():
    assert len(Interval(1000, 1200)) == 200


def test_interval_rejects_end_le_start():
    with pytest.raises(ValueError):
        Interval(1200, 1000)


def test_interval_rejects_zero_width():
    with pytest.raises(ValueError):
        Interval(1000, 1000)


def test_interval_rejects_negative_start():
    with pytest.raises(ValueError):
        Interval(-1, 100)


def test_interval_rejects_non_int():
    with pytest.raises(ValueError):
        Interval(1000.5, 1200)


def test_interval_rejects_bool():
    with pytest.raises(ValueError):
        Interval(True, 1200)


# --------------------------------------------------------------------------
# Exon
# --------------------------------------------------------------------------

def test_exon_valid():
    e = Exon(1000, 1200)
    assert len(e) == 200


def test_exon_single_base_min_width():
    assert len(Exon(1000, 1001)) == 1


def test_exon_rejects_inverted():
    with pytest.raises(ValueError):
        Exon(1200, 1000)


def test_exon_rejects_negative():
    with pytest.raises(ValueError):
        Exon(-5, 10)


# --------------------------------------------------------------------------
# CDSSegment
# --------------------------------------------------------------------------

def test_cds_segment_valid_phase0():
    c = CDSSegment(1000, 1150, 0)
    assert len(c) == 150
    assert c.phase == 0


@pytest.mark.parametrize("phase", [0, 1, 2])
def test_cds_segment_all_valid_phases(phase):
    assert CDSSegment(1000, 1150, phase).phase == phase


def test_cds_segment_default_phase_zero():
    assert CDSSegment(1000, 1150).phase == 0


def test_cds_segment_rejects_phase_3():
    with pytest.raises(ValueError):
        CDSSegment(1000, 1150, 3)


def test_cds_segment_rejects_negative_phase():
    with pytest.raises(ValueError):
        CDSSegment(1000, 1150, -1)


def test_cds_segment_rejects_inverted():
    with pytest.raises(ValueError):
        CDSSegment(1150, 1000, 0)


# --------------------------------------------------------------------------
# SpliceJunction
# --------------------------------------------------------------------------

def test_splice_junction_valid_plus():
    j = SpliceJunction("chr1", 1200, 1300, "+", read_count=10)
    assert j.intron_length == 100
    assert j.samples == 1


def test_splice_junction_valid_minus():
    j = SpliceJunction("chr1", 2200, 2300, "-", read_count=5, samples=3)
    assert j.intron_length == 100
    assert j.samples == 3


def test_splice_junction_rejects_donor_ge_acceptor():
    with pytest.raises(ValueError):
        SpliceJunction("chr1", 1300, 1200, "+", read_count=10)


def test_splice_junction_rejects_bad_strand():
    with pytest.raises(ValueError):
        SpliceJunction("chr1", 1200, 1300, ".", read_count=10)


def test_splice_junction_rejects_negative_read_count():
    with pytest.raises(ValueError):
        SpliceJunction("chr1", 1200, 1300, "+", read_count=-1)


def test_splice_junction_rejects_zero_samples():
    with pytest.raises(ValueError):
        SpliceJunction("chr1", 1200, 1300, "+", read_count=10, samples=0)


def test_splice_junction_zero_read_count_allowed():
    assert SpliceJunction("chr1", 1200, 1300, "+", read_count=0).read_count == 0


# --------------------------------------------------------------------------
# HelixerLocus
# --------------------------------------------------------------------------

def test_helixer_locus_valid_plus(sample_helixer_locus):
    assert sample_helixer_locus.strand == "+"
    assert sample_helixer_locus.span.start == 1000
    assert sample_helixer_locus.span.end == 1800


def test_helixer_locus_span_is_interval(sample_helixer_locus):
    assert isinstance(sample_helixer_locus.span, Interval)
    assert len(sample_helixer_locus.span) == 800


def test_helixer_locus_valid_minus(sample_exons_minus):
    loc = HelixerLocus("g", "chr1", 2000, 2800, "-", exons=sample_exons_minus)
    assert loc.strand == "-"


def test_helixer_locus_empty_exons_allowed():
    loc = HelixerLocus("g", "chr1", 1000, 1800, "+")
    assert loc.exons == []


def test_helixer_locus_confidence_none_ok():
    loc = HelixerLocus("g", "chr1", 1000, 1800, "+")
    assert loc.confidence is None


def test_helixer_locus_confidence_in_range():
    assert HelixerLocus("g", "chr1", 1000, 1800, "+", confidence=0.0).confidence == 0.0
    assert HelixerLocus("g", "chr1", 1000, 1800, "+", confidence=1.0).confidence == 1.0


def test_helixer_locus_rejects_confidence_above_one():
    with pytest.raises(ValueError):
        HelixerLocus("g", "chr1", 1000, 1800, "+", confidence=1.5)


def test_helixer_locus_rejects_confidence_negative():
    with pytest.raises(ValueError):
        HelixerLocus("g", "chr1", 1000, 1800, "+", confidence=-0.1)


def test_helixer_locus_rejects_bad_strand():
    with pytest.raises(ValueError):
        HelixerLocus("g", "chr1", 1000, 1800, "*")


def test_helixer_locus_rejects_exon_outside_span(sample_exons_plus):
    with pytest.raises(ValueError):
        HelixerLocus("g", "chr1", 1000, 1700, "+", exons=sample_exons_plus)


def test_helixer_locus_rejects_unsorted_exons():
    bad = [Exon(1300, 1500), Exon(1000, 1200)]
    with pytest.raises(ValueError):
        HelixerLocus("g", "chr1", 1000, 1800, "+", exons=bad)


def test_helixer_locus_rejects_overlapping_exons():
    bad = [Exon(1000, 1300), Exon(1200, 1500)]
    with pytest.raises(ValueError):
        HelixerLocus("g", "chr1", 1000, 1800, "+", exons=bad)


def test_helixer_locus_rejects_cds_outside_exons(sample_exons_plus):
    bad_cds = [CDSSegment(1210, 1290, 0)]  # in intron, no exon contains it
    with pytest.raises(ValueError):
        HelixerLocus("g", "chr1", 1000, 1800, "+", exons=sample_exons_plus, cds=bad_cds)


# --------------------------------------------------------------------------
# StringTieTranscript
# --------------------------------------------------------------------------

def test_stringtie_valid(sample_stringtie_transcript):
    assert sample_stringtie_transcript.tpm == 12.5
    assert sample_stringtie_transcript.sample_id == "sample_a"


def test_stringtie_minus_strand(sample_exons_minus):
    t = StringTieTranscript("t", "g", "chr1", 2000, 2800, "-", sample_exons_minus, 3.0, "s")
    assert t.strand == "-"


def test_stringtie_rejects_negative_tpm(sample_exons_plus):
    with pytest.raises(ValueError):
        StringTieTranscript("t", "g", "chr1", 1000, 1800, "+", sample_exons_plus, -1.0, "s")


def test_stringtie_rejects_empty_exons():
    with pytest.raises(ValueError):
        StringTieTranscript("t", "g", "chr1", 1000, 1800, "+", [], 5.0, "s")


def test_stringtie_rejects_negative_coverage(sample_exons_plus):
    with pytest.raises(ValueError):
        StringTieTranscript(
            "t", "g", "chr1", 1000, 1800, "+", sample_exons_plus, 5.0, "s", coverage=-2.0
        )


def test_stringtie_rejects_exon_outside_span(sample_exons_plus):
    with pytest.raises(ValueError):
        StringTieTranscript("t", "g", "chr1", 1100, 1800, "+", sample_exons_plus, 5.0, "s")


# --------------------------------------------------------------------------
# MiniprotAlignment
# --------------------------------------------------------------------------

def test_miniprot_valid(sample_miniprot_alignment):
    assert sample_miniprot_alignment.identity == 0.88
    assert sample_miniprot_alignment.rank == 0


def test_miniprot_minus_strand(sample_cds_segments):
    a = MiniprotAlignment("p", "chr1", 1000, 1800, "-", sample_cds_segments, 0.9, 0.8, 100.0)
    assert a.strand == "-"


def test_miniprot_rejects_coverage_above_one(sample_cds_segments):
    with pytest.raises(ValueError):
        MiniprotAlignment("p", "chr1", 1000, 1800, "+", sample_cds_segments, 1.2, 0.8, 100.0)


def test_miniprot_rejects_identity_above_one(sample_cds_segments):
    with pytest.raises(ValueError):
        MiniprotAlignment("p", "chr1", 1000, 1800, "+", sample_cds_segments, 0.9, 1.1, 100.0)


def test_miniprot_rejects_empty_cds():
    with pytest.raises(ValueError):
        MiniprotAlignment("p", "chr1", 1000, 1800, "+", [], 0.9, 0.8, 100.0)


def test_miniprot_rejects_cds_outside_span(sample_cds_segments):
    with pytest.raises(ValueError):
        MiniprotAlignment("p", "chr1", 1100, 1800, "+", sample_cds_segments, 0.9, 0.8, 100.0)


def test_miniprot_coverage_boundary_values(sample_cds_segments):
    a = MiniprotAlignment("p", "chr1", 1000, 1800, "+", sample_cds_segments, 0.0, 1.0, 1.0)
    assert a.query_coverage == 0.0
    assert a.identity == 1.0


# --------------------------------------------------------------------------
# LocusClassification
# --------------------------------------------------------------------------

@pytest.mark.parametrize("status", ["EXPRESSED", "LOW", "SILENT"])
def test_classification_valid_statuses(status):
    assert LocusClassification("loc", status).status == status


def test_classification_rejects_bad_status():
    with pytest.raises(ValueError):
        LocusClassification("loc", "MAYBE")


@pytest.mark.parametrize("src", ["stringtie", "bam_coverage", "bigwig", "none"])
def test_classification_valid_sources(src):
    assert LocusClassification("loc", "LOW", evidence_source=src).evidence_source == src


def test_classification_rejects_bad_source():
    with pytest.raises(ValueError):
        LocusClassification("loc", "LOW", evidence_source="rnaseq")


def test_classification_rejects_negative_num_samples():
    with pytest.raises(ValueError):
        LocusClassification("loc", "LOW", num_samples_expressed=-1)


def test_classification_defaults():
    c = LocusClassification("loc", "SILENT")
    assert c.evidence_source == "none"
    assert c.num_samples_expressed == 0
    assert c.max_tpm is None


# --------------------------------------------------------------------------
# TranscriptCandidate
# --------------------------------------------------------------------------

def test_candidate_valid(sample_transcript_candidate):
    t = sample_transcript_candidate
    assert t.is_primary is True
    assert t.source == "mikado"


@pytest.mark.parametrize("source", ["stringtie", "helixer", "miniprot", "mikado"])
def test_candidate_valid_sources(source, sample_exons_plus):
    t = TranscriptCandidate("t", "loc", source, "chr1", 1000, 1800, "+", sample_exons_plus)
    assert t.source == source


def test_candidate_rejects_bad_source(sample_exons_plus):
    with pytest.raises(ValueError):
        TranscriptCandidate("t", "loc", "augustus", "chr1", 1000, 1800, "+", sample_exons_plus)


def test_candidate_num_exons(sample_transcript_candidate):
    assert sample_transcript_candidate.num_exons == 3


def test_candidate_num_introns(sample_transcript_candidate):
    assert sample_transcript_candidate.num_introns == 2


def test_candidate_total_exon_length(sample_transcript_candidate):
    # 200 + 200 + 200
    assert sample_transcript_candidate.total_exon_length == 600


def test_candidate_total_cds_length(sample_transcript_candidate):
    # 150 + 200 + 49
    assert sample_transcript_candidate.total_cds_length == 399


def test_candidate_total_cds_length_none(sample_exons_plus):
    t = TranscriptCandidate("t", "loc", "stringtie", "chr1", 1000, 1800, "+", sample_exons_plus)
    assert t.total_cds_length == 0


def test_candidate_introns_coordinates(sample_transcript_candidate):
    introns = sample_transcript_candidate.introns
    assert len(introns) == 2
    assert introns[0].start == 1200
    assert introns[0].end == 1300
    assert introns[1].start == 1500
    assert introns[1].end == 1600


def test_candidate_single_exon_no_introns():
    t = TranscriptCandidate("t", "loc", "helixer", "chr1", 1000, 1200, "+", [Exon(1000, 1200)])
    assert t.num_exons == 1
    assert t.num_introns == 0
    assert t.introns == []


def test_candidate_minus_strand(sample_exons_minus):
    t = TranscriptCandidate("t", "loc", "mikado", "chr1", 2000, 2800, "-", sample_exons_minus)
    assert t.strand == "-"
    assert t.num_introns == 2


def test_candidate_minus_strand_with_cds(sample_exons_minus):
    cds = [CDSSegment(2050, 2200, 0), CDSSegment(2300, 2500, 0), CDSSegment(2600, 2649, 1)]
    t = TranscriptCandidate(
        "t", "loc", "mikado", "chr1", 2000, 2800, "-", sample_exons_minus, cds=cds
    )
    assert t.total_cds_length == 399  # divisible by 3


def test_candidate_rejects_cds_not_mod3(sample_exons_plus):
    bad = [CDSSegment(1050, 1200, 0)]  # length 150 ok... make it not mod3
    bad = [CDSSegment(1050, 1201, 0)]  # length 151
    # 1201 is within exon (1000,1200)? no -> would also fail within-exon; choose within exon
    bad = [CDSSegment(1050, 1100, 0), CDSSegment(1300, 1302, 0)]  # 50 + 2 = 52, not mod3
    with pytest.raises(ValueError):
        TranscriptCandidate("t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus, cds=bad)


def test_candidate_allows_partial_cds_plus(sample_exons_plus):
    # 50 + 2 = 52 bp, not divisible by 3, but a partial ORF -> accepted
    partial = [CDSSegment(1050, 1100, 0), CDSSegment(1300, 1302, 0)]
    t = TranscriptCandidate(
        "t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus,
        cds=partial, cds_partial=True,
    )
    assert t.cds_partial is True
    assert t.total_cds_length == 52  # not mod3, allowed because partial


def test_candidate_allows_partial_cds_minus(sample_exons_minus):
    # minus strand partial ORF: 50 + 2 = 52 bp within sample_exons_minus
    partial = [CDSSegment(2050, 2100, 0), CDSSegment(2300, 2302, 0)]
    t = TranscriptCandidate(
        "t", "loc", "mikado", "chr1", 2000, 2800, "-", sample_exons_minus,
        cds=partial, cds_partial=True,
    )
    assert t.cds_partial is True
    assert t.total_cds_length == 52


def test_candidate_partial_default_false(sample_exons_plus):
    t = TranscriptCandidate("t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus)
    assert t.cds_partial is False


def test_candidate_complete_cds_still_enforces_mod3(sample_exons_plus):
    # same non-mod3 CDS WITHOUT the partial flag must still be rejected
    bad = [CDSSegment(1050, 1100, 0), CDSSegment(1300, 1302, 0)]
    with pytest.raises(ValueError):
        TranscriptCandidate(
            "t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus,
            cds=bad, cds_partial=False,
        )


def test_candidate_rejects_cds_outside_exons(sample_exons_plus):
    bad = [CDSSegment(1210, 1290, 0)]  # in intron
    with pytest.raises(ValueError):
        TranscriptCandidate("t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus, cds=bad)


def test_candidate_rejects_fraction_above_one(sample_exons_plus):
    with pytest.raises(ValueError):
        TranscriptCandidate(
            "t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus,
            junction_support_fraction=1.5,
        )


def test_candidate_rejects_confidence_negative(sample_exons_plus):
    with pytest.raises(ValueError):
        TranscriptCandidate(
            "t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus, confidence=-0.5
        )


def test_candidate_rejects_empty_exons():
    with pytest.raises(ValueError):
        TranscriptCandidate("t", "loc", "mikado", "chr1", 1000, 1800, "+", [])


def test_candidate_rejects_bad_strand(sample_exons_plus):
    with pytest.raises(ValueError):
        TranscriptCandidate("t", "loc", "mikado", "chr1", 1000, 1800, "x", sample_exons_plus)


def test_candidate_defaults(sample_exons_plus):
    t = TranscriptCandidate("t", "loc", "stringtie", "chr1", 1000, 1800, "+", sample_exons_plus)
    assert t.cds is None
    assert t.is_primary is False
    assert t.tpm is None


# --------------------------------------------------------------------------
# QCFlag
# --------------------------------------------------------------------------

def test_qcflag_valid():
    f = QCFlag("X", "structure", "ERROR", "desc")
    assert f.name == "X"
    assert f.description == "desc"


def test_qcflag_is_frozen():
    f = QCFlag("X", "structure", "ERROR")
    with pytest.raises(attrs.exceptions.FrozenInstanceError):
        f.name = "Y"


@pytest.mark.parametrize("cat", ["evidence", "confidence", "structure", "homology", "splice", "locus"])
def test_qcflag_valid_categories(cat):
    assert QCFlag("X", cat, "INFO").category == cat


@pytest.mark.parametrize("sev", ["INFO", "WARNING", "ERROR", "CRITICAL"])
def test_qcflag_valid_severities(sev):
    assert QCFlag("X", "locus", sev).severity == sev


def test_qcflag_rejects_bad_category():
    with pytest.raises(ValueError):
        QCFlag("X", "weird", "INFO")


def test_qcflag_rejects_bad_severity():
    with pytest.raises(ValueError):
        QCFlag("X", "locus", "FATAL")


def test_qcflag_default_description():
    assert QCFlag("X", "locus", "INFO").description == ""


# --------------------------------------------------------------------------
# ASEvent (v3)
# --------------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["ES", "IR", "A5", "A3", "ALT_TSS", "ALT_TES", "MX"])
def test_as_event_valid_kinds(kind):
    e = ASEvent(kind, "chr1", 1300, 1500, "+")
    assert e.kind == kind


def test_as_event_valid_minus(sample_exons_minus):
    e = ASEvent("IR", "chr1", 2300, 2500, "-", support_read_count=15)
    assert e.strand == "-"
    assert e.support_read_count == 15


def test_as_event_is_frozen():
    e = ASEvent("ES", "chr1", 1300, 1500, "+")
    with pytest.raises(attrs.exceptions.FrozenInstanceError):
        e.kind = "IR"


def test_as_event_rejects_bad_kind():
    with pytest.raises(ValueError):
        ASEvent("RI", "chr1", 1300, 1500, "+")


def test_as_event_rejects_bad_coords():
    with pytest.raises(ValueError):
        ASEvent("ES", "chr1", 1500, 1300, "+")


def test_as_event_rejects_bad_strand():
    with pytest.raises(ValueError):
        ASEvent("ES", "chr1", 1300, 1500, ".")


def test_as_event_rejects_negative_support():
    with pytest.raises(ValueError):
        ASEvent("ES", "chr1", 1300, 1500, "+", support_read_count=-1)


def test_as_event_default_support_zero():
    assert ASEvent("ES", "chr1", 1300, 1500, "+").support_read_count == 0


# --------------------------------------------------------------------------
# MikadoLocus (v3)
# --------------------------------------------------------------------------

def test_mikado_locus_valid(sample_mikado_locus):
    assert sample_mikado_locus.locus_id == "mikado.1"
    assert len(sample_mikado_locus.transcripts) == 1
    assert sample_mikado_locus.metrics["cdna_length"] == 600


def test_mikado_locus_default_metrics_scores():
    loc = MikadoLocus("m", "chr1", 1000, 1800, "+", transcripts=[])
    assert loc.metrics == {}
    assert loc.scores == {}


def test_mikado_locus_rejects_non_mikado_source(sample_exons_plus):
    t = TranscriptCandidate("t", "loc", "stringtie", "chr1", 1000, 1800, "+", sample_exons_plus)
    with pytest.raises(ValueError):
        MikadoLocus("m", "chr1", 1000, 1800, "+", transcripts=[t])


def test_mikado_locus_rejects_seqid_mismatch(sample_exons_plus):
    t = TranscriptCandidate("t", "loc", "mikado", "chr2", 1000, 1800, "+", sample_exons_plus)
    with pytest.raises(ValueError):
        MikadoLocus("m", "chr1", 1000, 1800, "+", transcripts=[t])


def test_mikado_locus_rejects_strand_mismatch(sample_exons_plus):
    t = TranscriptCandidate("t", "loc", "mikado", "chr1", 1000, 1800, "+", sample_exons_plus)
    with pytest.raises(ValueError):
        MikadoLocus("m", "chr1", 1000, 1800, "-", transcripts=[t])


def test_mikado_locus_minus_strand(sample_exons_minus):
    t = TranscriptCandidate("t", "loc", "mikado", "chr1", 2000, 2800, "-", sample_exons_minus)
    loc = MikadoLocus("m", "chr1", 2000, 2800, "-", transcripts=[t])
    assert loc.strand == "-"


def test_mikado_locus_rejects_bad_coords():
    with pytest.raises(ValueError):
        MikadoLocus("m", "chr1", 1800, 1000, "+", transcripts=[])


# --------------------------------------------------------------------------
# ReconciledGene (v3)
# --------------------------------------------------------------------------

def test_reconciled_gene_valid(sample_reconciled_gene):
    g = sample_reconciled_gene
    assert g.gene_id == "HFG_00001"
    assert g.tier == 1
    assert g.origin == "mikado_1to1"
    assert len(g.as_events) == 1


@pytest.mark.parametrize("tier", [1, 2, 3, 4])
def test_reconciled_gene_valid_tiers(tier, sample_transcript_candidate, sample_classification):
    g = ReconciledGene(
        "g", "chr1", 1000, 1800, "+", tier, [sample_transcript_candidate],
        "HFG_00001.1", sample_classification, "mikado_1to1",
    )
    assert g.tier == tier


def test_reconciled_gene_rejects_bad_tier(sample_transcript_candidate, sample_classification):
    with pytest.raises(ValueError):
        ReconciledGene(
            "g", "chr1", 1000, 1800, "+", 5, [sample_transcript_candidate],
            "HFG_00001.1", sample_classification, "mikado_1to1",
        )


@pytest.mark.parametrize(
    "origin", ["mikado_1to1", "split", "merge", "helixer_backstop", "novel"]
)
def test_reconciled_gene_valid_origins(origin, sample_transcript_candidate, sample_classification):
    g = ReconciledGene(
        "g", "chr1", 1000, 1800, "+", 1, [sample_transcript_candidate],
        "HFG_00001.1", sample_classification, origin,
    )
    assert g.origin == origin


def test_reconciled_gene_rejects_bad_origin(sample_transcript_candidate, sample_classification):
    with pytest.raises(ValueError):
        ReconciledGene(
            "g", "chr1", 1000, 1800, "+", 1, [sample_transcript_candidate],
            "HFG_00001.1", sample_classification, "fusion",
        )


def test_reconciled_gene_rejects_no_transcripts(sample_classification):
    with pytest.raises(ValueError):
        ReconciledGene(
            "g", "chr1", 1000, 1800, "+", 1, [], "x", sample_classification, "novel"
        )


def test_reconciled_gene_rejects_primary_not_in_transcripts(
    sample_transcript_candidate, sample_classification
):
    with pytest.raises(ValueError):
        ReconciledGene(
            "g", "chr1", 1000, 1800, "+", 1, [sample_transcript_candidate],
            "NOT_THERE", sample_classification, "mikado_1to1",
        )


def test_reconciled_gene_rejects_transcript_strand_mismatch(
    sample_exons_minus, sample_classification
):
    t = TranscriptCandidate("t", "loc", "mikado", "chr1", 2000, 2800, "-", sample_exons_minus)
    with pytest.raises(ValueError):
        ReconciledGene(
            "g", "chr1", 2000, 2800, "+", 1, [t], "t", sample_classification, "novel"
        )


def test_reconciled_gene_rejects_transcript_seqid_mismatch(
    sample_exons_plus, sample_classification
):
    t = TranscriptCandidate("t", "loc", "mikado", "chr2", 1000, 1800, "+", sample_exons_plus)
    with pytest.raises(ValueError):
        ReconciledGene(
            "g", "chr1", 1000, 1800, "+", 1, [t], "t", sample_classification, "novel"
        )


def test_reconciled_gene_minus_strand(sample_exons_minus, sample_classification):
    t = TranscriptCandidate("t", "loc", "mikado", "chr1", 2000, 2800, "-", sample_exons_minus)
    g = ReconciledGene(
        "g", "chr1", 2000, 2800, "-", 2, [t], "t", sample_classification, "helixer_backstop"
    )
    assert g.strand == "-"


def test_reconciled_gene_defaults(sample_transcript_candidate, sample_classification):
    g = ReconciledGene(
        "g", "chr1", 1000, 1800, "+", 1, [sample_transcript_candidate],
        "HFG_00001.1", sample_classification, "mikado_1to1",
    )
    assert g.as_events == []
    assert g.flags == []
    assert g.merged_from == []


def test_reconciled_gene_merged_from(sample_transcript_candidate, sample_classification):
    g = ReconciledGene(
        "g", "chr1", 1000, 1800, "+", 1, [sample_transcript_candidate],
        "HFG_00001.1", sample_classification, "merge", merged_from=["HELIXER_1", "HELIXER_2"],
    )
    assert g.merged_from == ["HELIXER_1", "HELIXER_2"]


def test_reconciled_gene_evolve_revalidates(sample_reconciled_gene):
    with pytest.raises(ValueError):
        attrs.evolve(sample_reconciled_gene, tier=9)


# --------------------------------------------------------------------------
# IsoformAdmission (v3)
# --------------------------------------------------------------------------

def test_isoform_admission_valid():
    a = IsoformAdmission("t.1", "g1", True, "primary")
    assert a.admitted is True
    assert a.reason == "primary"


def test_isoform_admission_rejected_record():
    a = IsoformAdmission("t.2", "g1", False, "redundant", redundant_with="t.1")
    assert a.admitted is False
    assert a.redundant_with == "t.1"


def test_isoform_admission_with_novel_event(sample_as_event):
    a = IsoformAdmission("t.3", "g1", True, "novel AS", novel_as_event=sample_as_event)
    assert a.novel_as_event.kind == "ES"


def test_isoform_admission_rejects_non_bool_admitted():
    with pytest.raises(ValueError):
        IsoformAdmission("t.1", "g1", "yes", "reason")


def test_isoform_admission_rejects_bad_novel_event():
    with pytest.raises(ValueError):
        IsoformAdmission("t.1", "g1", True, "reason", novel_as_event="ES")


def test_isoform_admission_defaults():
    a = IsoformAdmission("t.1", "g1", True, "ok")
    assert a.novel_as_event is None
    assert a.score is None
    assert a.redundant_with is None
