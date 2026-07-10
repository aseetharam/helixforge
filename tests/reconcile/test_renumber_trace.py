"""_renumber explicit-order support (Phase 33b D2). Floor: 4. Both strands.

The TRaCE election supplies an explicit isoform order to
``mikado_integrate._renumber``; the default (no order) must reproduce the
historical ``combined_score`` sort exactly so every existing caller is unchanged.
"""

from __future__ import annotations

from helixforge.reconcile.mikado_integrate import _renumber
from helixforge.reconcile.models import Exon, TranscriptCandidate


def _tc(tid, *, score, strand="+", start=100, end=200):
    return TranscriptCandidate(
        tid, "G", "mikado", "chr1", start, end, strand,
        [Exon(start, end)], combined_score=score,
    )


def test_default_order_is_combined_score_desc():
    txs = [_tc("t1", score=1.0), _tc("t2", score=5.0), _tc("t3", score=3.0)]
    out = _renumber(txs, "HFG_00007")
    assert [t.transcript_id for t in out] == [
        "HFG_00007.1", "HFG_00007.2", "HFG_00007.3"
    ]
    # .1 is the highest score (t2); primary flag set on it alone.
    assert out[0].is_primary and out[0].num_exons == 1
    assert out[0].start == 100  # t2's structure
    assert [t.is_primary for t in out] == [True, False, False]


def test_default_order_tiebreak_by_transcript_id():
    # Equal scores → ascending transcript_id decides (historical behavior).
    txs = [_tc("t_b", score=2.0), _tc("t_a", score=2.0)]
    out = _renumber(txs, "HFG_1")
    # t_a sorts before t_b → it becomes .1.
    assert out[0].locus_id == "HFG_1"
    assert [t.transcript_id for t in out] == ["HFG_1.1", "HFG_1.2"]


def test_explicit_order_numbers_by_supplied_order():
    txs = [_tc("t1", score=5.0), _tc("t2", score=1.0), _tc("t3", score=3.0)]
    # Force the lowest-scored t2 to be primary via an explicit election order.
    out = _renumber(txs, "HFG_00009", order=["t2", "t3", "t1"])
    assert [t.transcript_id for t in out] == [
        "HFG_00009.1", "HFG_00009.2", "HFG_00009.3"
    ]
    assert out[0].is_primary
    # .1 must be the structure that was 't2'; identify it by being the only one
    # whose pre-renumber score was 1.0.
    assert out[0].combined_score == 1.0
    assert out[1].combined_score == 3.0
    assert out[2].combined_score == 5.0


def test_explicit_order_minus_strand_and_unnamed_appended():
    txs = [
        _tc("t1", score=5.0, strand="-"),
        _tc("t2", score=1.0, strand="-"),
        _tc("t3", score=9.0, strand="-"),
    ]
    # Only name two; the unnamed remainder falls back to the score order, last.
    out = _renumber(txs, "HFG_2", order=["t2", "t1"])
    assert [t.combined_score for t in out] == [1.0, 5.0, 9.0]
    assert out[0].is_primary and out[0].strand == "-"
    assert [t.transcript_id for t in out] == ["HFG_2.1", "HFG_2.2", "HFG_2.3"]
