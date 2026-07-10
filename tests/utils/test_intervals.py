"""Tests for the NCLS-backed IntervalIndex (Phase 0). Floor: 12."""

import pytest

from helixforge.reconcile.models import CDSSegment, Exon
from helixforge.utils.intervals import IntervalIndex, bounds


def test_empty_index_is_empty():
    idx = IntervalIndex()
    assert idx.is_empty
    assert len(idx) == 0


def test_empty_index_query_returns_empty():
    assert IntervalIndex().query(100, 200) == []


def test_add_and_len():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200), (1300, 1500)])
    assert len(idx) == 2
    assert not idx.is_empty


def test_query_single_overlap():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200), (1300, 1500)])
    assert idx.query(1100, 1150) == [0]


def test_query_no_overlap():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200), (1300, 1500)])
    assert idx.query(1250, 1290) == []


def test_query_multiple_overlap():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200), (1100, 1300), (1250, 1400)])
    assert idx.query(1150, 1260) == [0, 1, 2]


def test_query_half_open_boundary_excludes_touching():
    # interval [1000,1200) and query [1200,1300) do not overlap (half-open)
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200)])
    assert idx.query(1200, 1300) == []


def test_query_half_open_boundary_includes_overlap():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200)])
    assert idx.query(1199, 1300) == [0]


def test_query_with_data_returns_payload():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200, "exon_a"), (1300, 1500, "exon_b")])
    assert idx.query_with_data(1100, 1150) == [(1000, 1200, "exon_a")]


def test_query_with_data_multiple_payloads():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1300, "a"), (1200, 1400, "b")])
    assert idx.query_with_data(1250, 1260) == [(1000, 1300, "a"), (1200, 1400, "b")]


def test_add_intervals_incremental():
    idx = IntervalIndex()
    idx.add_intervals([(1000, 1200)])
    idx.add_intervals([(1300, 1500)])
    assert len(idx) == 2
    assert idx.query(1350, 1400) == [1]


def test_add_rejects_zero_width():
    idx = IntervalIndex()
    with pytest.raises(ValueError):
        idx.add_intervals([(1000, 1000)])


def test_add_rejects_inverted():
    idx = IntervalIndex()
    with pytest.raises(ValueError):
        idx.add_intervals([(1200, 1000)])


def test_add_rejects_too_short_tuple():
    idx = IntervalIndex()
    with pytest.raises(ValueError):
        idx.add_intervals([(1000,)])


def test_query_ids_sorted():
    idx = IntervalIndex()
    idx.add_intervals([(1300, 1500), (1000, 1200), (1100, 1400)])
    # ids are insertion order; overlapping query returns sorted ids
    assert idx.query(1050, 1450) == [0, 1, 2]


# ---------------------------------------------------------------------------
# Phase 19 D1: bounds() unified here (was duplicated in mikado_integrate +
# as_events). Coordinates are stored (low, high) regardless of strand
# (CLAUDE.md §4), so bounds is strand-agnostic by design.
# ---------------------------------------------------------------------------

def test_bounds_on_tuple_input():
    # plus strand and minus strand both store (low, high) — same tuple shape
    assert bounds((1000, 1200)) == (1000, 1200)
    # extra payload elements (start, end, *data) are ignored
    assert bounds((1300, 1500, "exon_a", "+")) == (1300, 1500)


def test_bounds_on_object_input_both_strands():
    plus_exon = Exon(1000, 1200)            # interval from a '+' gene
    minus_cds = CDSSegment(2000, 2300, 0)   # interval from a '-' gene
    assert bounds(plus_exon) == (1000, 1200)
    assert bounds(minus_cds) == (2000, 2300)
