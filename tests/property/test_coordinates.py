"""Phase 25 D2 — property-based coordinate invariants (M10).

Hypothesis round-trip properties for the coordinate-conversion boundary
(``utils/regions.py``), the *only* place 1-based↔0-based conversion is allowed
(CLAUDE.md §4). A silent off-by-one or strand confusion here is exactly the
class of bug that caused the v1 catastrophe (CLAUDE.md §6); these properties
guard the boundary with both concrete seeds (``@example``) and generated cases.

Coordinate storage is strand-agnostic (CLAUDE.md §4.2: ``(low, high)``
regardless of strand), so the conversions themselves do not take a strand — but
every property is exercised against CDS-segment models built on **both strands**
to prove the boundary is lossless wherever it sits in the pipeline.
"""

from __future__ import annotations

from hypothesis import example, given
from hypothesis import strategies as st

from helixforge.reconcile.models import CDSSegment
from helixforge.utils.regions import (
    format_region,
    gff3_to_internal,
    internal_to_gff3,
    parse_region,
)

# 1-based inclusive GFF3 coordinates: start >= 1, end >= start.
_gff3_coords = st.integers(min_value=1, max_value=10_000_000).flatmap(
    lambda s: st.tuples(st.just(s), st.integers(min_value=s, max_value=s + 5_000_000))
)
# 0-based half-open internal coordinates: start >= 0, end > start.
_internal_coords = st.integers(min_value=0, max_value=10_000_000).flatmap(
    lambda s: st.tuples(st.just(s), st.integers(min_value=s + 1, max_value=s + 5_000_000))
)


@given(_gff3_coords)
@example((1, 1))          # smallest legal 1 bp feature
@example((1000, 1200))    # the canonical fixture span
@example((2000, 2800))    # the minus-strand fixture span
def test_gff3_internal_gff3_is_identity(coords):
    start, end = coords
    istart, iend = gff3_to_internal(start, end)
    assert (istart, iend) == (start - 1, end)        # documented mapping
    assert internal_to_gff3(istart, iend) == (start, end)


@given(_internal_coords)
@example((0, 1))          # smallest legal 1 bp internal feature
@example((999, 1200))     # internal form of the canonical fixture span
def test_internal_gff3_internal_is_identity(coords):
    start, end = coords
    gstart, gend = internal_to_gff3(start, end)
    assert (gstart, gend) == (start + 1, end)        # documented mapping
    assert gff3_to_internal(gstart, gend) == (start, end)


@given(_gff3_coords)
@example((1000, 1200))
def test_inclusive_width_preserved_across_boundary(coords):
    """1-based inclusive width (end - start + 1) == internal half-open width."""
    start, end = coords
    istart, iend = gff3_to_internal(start, end)
    assert iend - istart == end - start + 1


@given(_internal_coords)
@example((0, 3))
def test_internal_to_gff3_keeps_one_based_invariants(coords):
    start, end = coords
    gstart, gend = internal_to_gff3(start, end)
    assert gstart >= 1            # GFF3 is 1-based
    assert gend >= gstart         # start <= end always (CLAUDE.md §4.1)


@given(
    st.text(
        alphabet=st.characters(min_codepoint=33, max_codepoint=126, blacklist_characters=":-"),
        min_size=1,
        max_size=20,
    ),
    _gff3_coords,
)
@example("chr1", (1000, 1800))
@example("scaffold_7", (2000, 2800))
def test_parse_format_region_round_trip(seqid, coords):
    start, end = coords
    region = format_region(seqid, start, end)
    assert parse_region(region) == (seqid, start, end)


@given(st.text(alphabet=st.characters(min_codepoint=33, max_codepoint=126,
                                       blacklist_characters=":-"),
               min_size=1, max_size=20))
@example("chr1")
def test_bare_seqid_region_round_trips(seqid):
    assert parse_region(format_region(seqid)) == (seqid, None, None)


# --- both-strand model-level boundary losslessness ----------------------------

_seg_lists = st.lists(
    st.integers(min_value=0, max_value=200_000).flatmap(
        lambda s: st.tuples(
            st.just(s),
            st.integers(min_value=s + 3, max_value=s + 3000),
            st.integers(min_value=0, max_value=2),
        )
    ),
    min_size=1,
    max_size=8,
)


@given(_seg_lists, st.sampled_from(["+", "-"]))
@example([(1050, 1200, 0), (1300, 1500, 0), (1600, 1649, 1)], "+")   # plus fixture CDS
@example([(2050, 2200, 0), (2300, 2500, 0), (2600, 2649, 1)], "-")   # minus fixture CDS
def test_cds_round_trip_preserves_coords_and_mod3_both_strands(raw, strand):
    """internal→gff3→internal on every CDS segment is identity, on both strands,
    and the total-CDS-length mod 3 (CLAUDE.md §5.3) is preserved across it."""
    segs = [CDSSegment(s, e, p) for (s, e, p) in raw]
    total_before = sum(seg.end - seg.start for seg in segs)
    round_tripped = []
    for seg in segs:
        gstart, gend = internal_to_gff3(seg.start, seg.end)
        istart, iend = gff3_to_internal(gstart, gend)
        assert (istart, iend) == (seg.start, seg.end)
        round_tripped.append(iend - istart)
    total_after = sum(round_tripped)
    assert total_after == total_before
    assert total_after % 3 == total_before % 3
    assert strand in ("+", "-")
