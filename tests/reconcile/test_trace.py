"""Tests for reconcile/trace.py — TRaCE canonical-transcript election (Phase 33b).

Floor: 20. Concrete literal coordinates; both strands. Synthetic candidates +
synthetic per-sample StringTie transcripts with known TPM / overlap.
"""

from __future__ import annotations

import pytest

from helixforge.reconcile.as_events import derive_as_events
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    StringTieTranscript,
    TranscriptCandidate,
)
from helixforge.reconcile.trace import (
    Ballot,
    TraceParams,
    aed_ballot,
    domain_coverage_by_transcript,
    elect,
    length_ballot,
    proportion_overlap,
    structural_aed,
    trace_order,
)


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def _tc(tid, exons, *, strand="+", cds=None, score=None, seqid="chr1"):
    start = min(e.start for e in exons)
    end = max(e.end for e in exons)
    return TranscriptCandidate(
        tid, "G", "mikado", seqid, start, end, strand, exons,
        cds=cds, combined_score=score,
    )


def _st(tid, exons, *, tpm, strand="+", seqid="chr1", sample="s1"):
    start = min(e.start for e in exons)
    end = max(e.end for e in exons)
    return StringTieTranscript(tid, "g", seqid, start, end, strand, exons, tpm, sample)


def _ballot(ranks, voter="x"):
    return Ballot(ranks=ranks, kind="length", voter=voter)


# ---------------------------------------------------------------------------
# Geometry: proportion overlap + structural AED
# ---------------------------------------------------------------------------

def test_proportion_overlap_half():
    # candidate 100bp, assembled 100bp, overlap 50 → 50 / min(100,100) = 0.5
    assert proportion_overlap([(100, 200)], [(150, 250)]) == 0.5


def test_structural_aed_identical_is_zero():
    a = [(100, 200), (300, 400)]
    assert structural_aed(a, a) == 0.0


def test_structural_aed_intron_retention_partial():
    # spliced candidate vs single-exon (intron-retaining) assembled transcript.
    cand = [(100, 200), (300, 400)]   # 200 bp, on overlap region [100,400]
    asm = [(100, 400)]                # 300 bp
    # ov = 200; sn = 200/300, sp = 200/200=1 → AED = 1 - (0.6667+1)/2
    assert structural_aed(cand, asm) == pytest.approx(0.166667, abs=1e-5)


def test_structural_aed_disjoint_is_one():
    assert structural_aed([(100, 200)], [(5000, 5100)]) == 1.0


# ---------------------------------------------------------------------------
# AED ballot (sample voter) — cutoffs honored, both strands
# ---------------------------------------------------------------------------

def test_aed_ballot_ranks_by_overlap_region_aed():
    a = _tc("A", [Exon(100, 200), Exon(300, 400)])
    b = _tc("B", [Exon(100, 250)])
    sample = [_st("s.t1", [Exon(100, 200), Exon(300, 400)], tpm=10.0)]
    ranks = aed_ballot([a, b], sample)
    assert ranks == {"A": 1, "B": 2}


def test_aed_ballot_max_aed_excludes_distant_candidate():
    a = _tc("A", [Exon(100, 200), Exon(300, 400)])
    b = _tc("B", [Exon(100, 250)])          # AED ~0.167 to the sample
    sample = [_st("s.t1", [Exon(100, 200), Exon(300, 400)], tpm=10.0)]
    ranks = aed_ballot([a, b], sample, max_aed=0.1)
    assert ranks == {"A": 1, "B": None}


def test_aed_ballot_min_tpm_excludes_low_expression():
    a = _tc("A", [Exon(100, 200), Exon(300, 400)])
    sample = [_st("s.t1", [Exon(100, 200), Exon(300, 400)], tpm=0.3)]
    ranks = aed_ballot([a], sample, min_tpm=0.5)
    assert ranks == {"A": None}


def test_aed_ballot_min_overlap_excludes_nonoverlapping():
    a = _tc("A", [Exon(100, 200)])
    far = _tc("F", [Exon(5000, 5100)])
    sample = [_st("s.t1", [Exon(100, 200)], tpm=10.0)]
    ranks = aed_ballot([a, far], sample)
    assert ranks == {"A": 1, "F": None}


def test_aed_ballot_minus_strand():
    a = _tc("A", [Exon(100, 200), Exon(300, 400)], strand="-")
    b = _tc("B", [Exon(100, 400)], strand="-")
    sample = [_st("s.t1", [Exon(100, 200), Exon(300, 400)], tpm=8.0, strand="-")]
    ranks = aed_ballot([a, b], sample)
    assert ranks == {"A": 1, "B": 2}


def test_aed_ballot_wrong_strand_does_not_vote():
    a = _tc("A", [Exon(100, 200)], strand="+")
    sample = [_st("s.t1", [Exon(100, 200)], tpm=10.0, strand="-")]
    assert aed_ballot([a], sample) == {"A": None}


# ---------------------------------------------------------------------------
# Length ballots
# ---------------------------------------------------------------------------

def test_length_ballot_protein_ranks_by_cds():
    a = _tc("A", [Exon(100, 400)], cds=[CDSSegment(100, 400, 0)])   # 300
    b = _tc("B", [Exon(100, 250)], cds=[CDSSegment(100, 250, 0)])   # 150
    assert length_ballot([a, b], "protein_length") == {"A": 1, "B": 2}


def test_length_ballot_cdna_ranks_by_exon_length():
    a = _tc("A", [Exon(100, 300)])   # 200
    b = _tc("B", [Exon(100, 400)])   # 300
    assert length_ballot([a, b], "cdna_length") == {"A": 2, "B": 1}


def test_length_ballot_domain_uses_values():
    a = _tc("A", [Exon(100, 200)])
    b = _tc("B", [Exon(100, 200)])
    c = _tc("C", [Exon(100, 200)])
    ranks = length_ballot([a, b, c], "domain_coverage",
                          {"A": 1.0, "B": 0.0, "C": 0.5})
    assert ranks == {"A": 1, "C": 2, "B": 3}


def test_length_ballot_ties_share_rank():
    a = _tc("A", [Exon(100, 400)])   # 300
    b = _tc("B", [Exon(100, 400)])   # 300
    c = _tc("C", [Exon(100, 250)])   # 150
    assert length_ballot([a, b, c], "cdna_length") == {"A": 1, "B": 1, "C": 3}


def test_length_ballot_bad_key_raises():
    a = _tc("A", [Exon(100, 200)])
    with pytest.raises(ValueError):
        length_ballot([a], "nonsense")


# ---------------------------------------------------------------------------
# RCV election — worked example, determinism, tie-breaks
# ---------------------------------------------------------------------------

def _worked_example_candidates():
    return [
        _tc("A", [Exon(100, 200)], score=1.0),
        _tc("B", [Exon(100, 200)], score=2.0),
        _tc("C", [Exon(100, 200)], score=3.0),
    ]


def _worked_example_ballots():
    # 3 ballots A>B>C, 2 ballots B>A>C → A wins seat 1 (3 vs 2 rank-1);
    # after removing A, B wins seat 2 on the transferred rank-2 votes.
    return [
        _ballot({"A": 1, "B": 2, "C": 3}),
        _ballot({"A": 1, "B": 2, "C": 3}),
        _ballot({"A": 1, "C": 2, "B": 3}),
        _ballot({"B": 1, "A": 2, "C": 3}),
        _ballot({"B": 1, "A": 2, "C": 3}),
    ]


def test_elect_worked_example_rank1_plurality():
    order = elect(_worked_example_candidates(), _worked_example_ballots(), {"x": 1.0})
    assert order == ["A", "B", "C"]


def test_elect_second_seat_decided_by_rank2():
    # Remove the seat-1 winner up front: the same ballots must seat B over C on
    # the rank-2 votes that A's ballots transfer.
    cands = [
        _tc("B", [Exon(100, 200)], score=2.0),
        _tc("C", [Exon(100, 200)], score=3.0),
    ]
    ballots = [
        _ballot({"B": 1, "C": 2}),
        _ballot({"B": 1, "C": 2}),
        _ballot({"C": 1, "B": 2}),
        _ballot({"B": 1, "C": 2}),
    ]
    assert elect(cands, ballots, {"x": 1.0}) == ["B", "C"]


def test_elect_is_deterministic_across_calls():
    cands = _worked_example_candidates()
    ballots = _worked_example_ballots()
    first = elect(cands, ballots, {"x": 1.0})
    for _ in range(5):
        assert elect(cands, ballots, {"x": 1.0}) == first


def test_elect_tiebreak_prefers_higher_combined_score():
    # No ballots → purely the deterministic tie-break (combined_score desc).
    cands = [_tc("A", [Exon(100, 200)], score=1.0),
             _tc("B", [Exon(100, 200)], score=2.0)]
    assert elect(cands, [], {}) == ["B", "A"]


def test_elect_tiebreak_then_transcript_id():
    # Equal scores → lowest transcript_id wins the tie.
    cands = [_tc("B", [Exon(100, 200)], score=1.0),
             _tc("A", [Exon(100, 200)], score=1.0)]
    assert elect(cands, [], {}) == ["A", "B"]


def test_electorate_balancing_changes_outcome():
    cands = [_tc("A", [Exon(100, 200)], score=1.0),
             _tc("B", [Exon(100, 200)], score=5.0)]
    # 4 sample ballots favour A; one length ballot (weight 1) favours B.
    sample_ballots = [
        Ballot({"A": 1, "B": 2}, "sample", f"s{i}") for i in range(4)
    ]
    length_ballots = [Ballot({"B": 1, "A": 2}, "length", "protein")]
    ballots = sample_ballots + length_ballots
    weights = {"protein": 1.0}
    # Unbalanced: 4 sample rank-1 votes swamp the single length voter → A.
    assert elect(cands, ballots, weights, balance_samples=False)[0] == "A"
    # Balanced: the 4 samples collectively carry the length-voter total (1.0),
    # the bloc ties, and B's higher combined_score breaks it → B.
    assert elect(cands, ballots, weights, balance_samples=True)[0] == "B"


# ---------------------------------------------------------------------------
# trace_order — fallbacks, end-to-end promotion, domain voter
# ---------------------------------------------------------------------------

def test_trace_order_single_transcript_returns_input():
    a = _tc("A", [Exon(100, 200)])
    assert trace_order([a], {"s1": []}) == [a]


def test_trace_order_no_evidence_returns_input_order():
    a = _tc("A", [Exon(100, 200), Exon(300, 400)], score=1.0)
    b = _tc("B", [Exon(100, 400)], score=5.0)
    # No samples, no domain coverage → fall back to the input (combined_score) order.
    assert trace_order([a, b], {}) == [a, b]


def test_trace_order_promotes_spliced_over_intron_retention():
    # B has the higher Mikado combined_score but is intron-retention-inflated;
    # the RNA-seq sample matches the spliced A, so TRaCE elects A as primary.
    a = _tc("A", [Exon(100, 200), Exon(300, 400)], score=1.0)
    b = _tc("B", [Exon(100, 400)], score=5.0)
    samples = {"s1": [_st("s1.t", [Exon(100, 200), Exon(300, 400)], tpm=10.0)]}
    ordered = trace_order([a, b], samples)
    assert [t.transcript_id for t in ordered] == ["A", "B"]


def test_trace_order_minus_strand_promotion():
    a = _tc("A", [Exon(100, 200), Exon(300, 400)], strand="-", score=1.0)
    b = _tc("B", [Exon(100, 400)], strand="-", score=5.0)
    samples = {"s1": [_st("s1.t", [Exon(100, 200), Exon(300, 400)],
                          tpm=10.0, strand="-")]}
    ordered = trace_order([a, b], samples)
    assert [t.transcript_id for t in ordered] == ["A", "B"]


def test_trace_order_domain_voter_present_vs_absent():
    # Two equal-length candidates, no RNA-seq match difference: the domain voter
    # (highest weight) decides when coverage is supplied.
    a = _tc("A", [Exon(100, 400)], score=2.0)
    b = _tc("B", [Exon(100, 400)], score=1.0)
    samples = {"s1": [_st("s1.t", [Exon(100, 400)], tpm=10.0)]}
    cov = {"A": 0.0, "B": 1.0}
    with_domain = trace_order(
        [a, b], samples, domain_coverage=cov,
        params=TraceParams(use_domain=True),
    )
    assert with_domain[0].transcript_id == "B"   # domain coverage promotes B
    # Without the domain voter, the tie falls to combined_score → A.
    without = trace_order(
        [a, b], samples, domain_coverage=cov,
        params=TraceParams(use_domain=False),
    )
    assert without[0].transcript_id == "A"


def test_domain_coverage_by_transcript_proxy():
    class _Rec:
        def __init__(self, complete):
            self.domain_complete = complete

    a = _tc("HFG_1.1", [Exon(100, 200)])
    b = _tc("HFG_1.2", [Exon(100, 200)])
    gene = _gene("HFG_1", [a, b])
    cov = domain_coverage_by_transcript(
        [gene], {"HFG_1.1": _Rec(True), "HFG_1.2": _Rec(False)}
    )
    assert cov == {"HFG_1.1": 1.0, "HFG_1.2": 0.0}
    assert domain_coverage_by_transcript([gene], {}) == {}


# ---------------------------------------------------------------------------
# Pipeline wiring: _maybe_trace_reorder
# ---------------------------------------------------------------------------

def _gene(gene_id, transcripts):
    primary = next((t for t in transcripts if t.is_primary), transcripts[0])
    start = min(t.start for t in transcripts)
    end = max(t.end for t in transcripts)
    return ReconciledGene(
        gene_id=gene_id, seqid=primary.seqid, start=start, end=end,
        strand=primary.strand, tier=2, transcripts=transcripts,
        primary_transcript_id=primary.transcript_id,
        classification=LocusClassification(gene_id, "EXPRESSED"),
        origin="mikado_1to1",
        as_events=derive_as_events(transcripts),
    )


def _two_isoform_gene():
    # .1 is the combined_score primary (spliced is NOT — IR isoform scored higher);
    # the sample favours the spliced .2, so TRaCE should flip the primary.
    t1 = TranscriptCandidate(
        "HFG_00001.1", "HFG_00001", "mikado", "chr1", 100, 400, "+",
        [Exon(100, 400)], combined_score=5.0, is_primary=True,
    )
    t2 = TranscriptCandidate(
        "HFG_00001.2", "HFG_00001", "mikado", "chr1", 100, 400, "+",
        [Exon(100, 200), Exon(300, 400)], combined_score=1.0,
    )
    return _gene("HFG_00001", [t1, t2])


def _trace_inputs():
    return {"per_sample": {"s1": [
        _st("s1.t", [Exon(100, 200), Exon(300, 400)], tpm=10.0),
    ]}}


def test_maybe_trace_reorder_off_is_noop():
    from helixforge.reconcile.pipeline import PipelineConfig, _maybe_trace_reorder

    gene = _two_isoform_gene()
    cfg = PipelineConfig(genome_fasta="g", helixer_gff3="h", trace_primary=False)
    out, func = _maybe_trace_reorder([gene], cfg, _trace_inputs(), {})
    assert out[0] is gene
    assert out[0].primary_transcript_id == "HFG_00001.1"


def test_maybe_trace_reorder_flips_primary_and_sets_rank():
    from helixforge.reconcile.pipeline import PipelineConfig, _maybe_trace_reorder

    gene = _two_isoform_gene()
    before_events = sorted((e.kind, e.start, e.end) for e in gene.as_events)
    cfg = PipelineConfig(genome_fasta="g", helixer_gff3="h", trace_primary=True)
    out, _ = _maybe_trace_reorder([gene], cfg, _trace_inputs(), {})
    g = out[0]
    # The new primary (.1) is the spliced structure (two exons).
    primary = next(t for t in g.transcripts if t.transcript_id == "HFG_00001.1")
    assert primary.num_exons == 2
    assert g.primary_transcript_id == "HFG_00001.1"
    # trace_rank persisted 1..N in elected order.
    assert [t.trace_rank for t in g.transcripts] == [1, 2]
    # AS-event set is preserved exactly (reorder-only, count-neutral).
    assert sorted((e.kind, e.start, e.end) for e in g.as_events) == before_events


def test_maybe_trace_reorder_rekeys_functional_records():
    from helixforge.reconcile.pipeline import PipelineConfig, _maybe_trace_reorder

    gene = _two_isoform_gene()
    # A functional record attached to the old primary (.1, the IR isoform).
    functional = {"HFG_00001.1": object()}
    rec = functional["HFG_00001.1"]
    cfg = PipelineConfig(genome_fasta="g", helixer_gff3="h", trace_primary=True)
    _, new_func = _maybe_trace_reorder([gene], cfg, _trace_inputs(), functional)
    # The IR isoform moved to .2, so its record follows it to the new id.
    assert new_func == {"HFG_00001.2": rec}
