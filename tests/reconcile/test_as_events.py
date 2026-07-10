"""Tests for reconcile/as_events.py (Phase 6). Floor: 18. Both strands."""

from helixforge.reconcile.as_events import (
    derive_as_events,
    is_redundant,
    overlap_bases,
    reciprocal_overlap,
)
from helixforge.reconcile.models import CDSSegment, Exon, TranscriptCandidate


def _tc(tid, exons, strand="+", cds=None, primary=False, score=None, seqid="chr1"):
    start = min(e.start for e in exons)
    end = max(e.end for e in exons)
    return TranscriptCandidate(
        tid, "L", "mikado", seqid, start, end, strand, exons,
        cds=cds, combined_score=score, is_primary=primary,
    )


def _kinds_with_coords(events):
    return {(e.kind, e.start, e.end) for e in events}


# --- overlap helpers ---

def test_overlap_bases_simple():
    assert overlap_bases([Exon(100, 200)], [Exon(150, 250)]) == 50


def test_overlap_bases_disjoint():
    assert overlap_bases([Exon(100, 150)], [Exon(300, 400)]) == 0


def test_overlap_bases_multi():
    a = [Exon(100, 200), Exon(300, 400)]
    b = [Exon(150, 350)]
    assert overlap_bases(a, b) == 50 + 50


def test_reciprocal_overlap_identical():
    assert reciprocal_overlap([Exon(100, 200)], [Exon(100, 200)]) == 1.0


def test_reciprocal_overlap_partial():
    # ov=50; la=100, lb=200 -> min(0.5, 0.25)=0.25
    assert reciprocal_overlap([Exon(100, 200)], [Exon(150, 350)]) == 0.25


def test_reciprocal_overlap_empty():
    assert reciprocal_overlap([], [Exon(100, 200)]) == 0.0


# --- derive_as_events: single isoform ---

def test_single_isoform_no_events():
    p = _tc("t1", [Exon(100, 200)], primary=True)
    assert derive_as_events([p]) == []


# --- exon skipping (plus) ---

def test_exon_skipping_plus():
    primary = _tc("p", [Exon(100, 150), Exon(200, 250), Exon(300, 350)], primary=True)
    alt = _tc("a", [Exon(100, 150), Exon(300, 350)])
    events = derive_as_events([primary, alt])
    assert ("ES", 200, 250) in _kinds_with_coords(events)


# --- intron retention (plus) ---

def test_intron_retention_plus():
    primary = _tc("p", [Exon(100, 150), Exon(200, 250)], primary=True)
    alt = _tc("a", [Exon(100, 250)])
    events = derive_as_events([primary, alt])
    assert ("IR", 150, 200) in _kinds_with_coords(events)


def test_intron_retention_no_false_es():
    primary = _tc("p", [Exon(100, 150), Exon(200, 250)], primary=True)
    alt = _tc("a", [Exon(100, 250)])
    kinds = {e.kind for e in derive_as_events([primary, alt])}
    assert "ES" not in kinds  # retention must not be misread as skipping


# --- alt 3' splice site (plus) ---

def test_alt_acceptor_a3_plus():
    primary = _tc("p", [Exon(100, 150), Exon(200, 300)], primary=True)
    alt = _tc("a", [Exon(100, 150), Exon(230, 300)])
    events = derive_as_events([primary, alt])
    # share donor 150, differ acceptor 200/230 -> A3 on + strand
    assert ("A3", 200, 230) in _kinds_with_coords(events)


# --- alt 5' splice site (plus) ---

def test_alt_donor_a5_plus():
    primary = _tc("p", [Exon(100, 150), Exon(200, 300)], primary=True)
    alt = _tc("a", [Exon(100, 170), Exon(200, 300)])
    events = derive_as_events([primary, alt])
    # share acceptor 200, differ donor 150/170 -> A5 on + strand
    assert ("A5", 150, 170) in _kinds_with_coords(events)


# --- ALT_TSS / ALT_TES (plus) ---

def test_alt_tss_plus():
    primary = _tc("p", [Exon(100, 150), Exon(200, 300)], primary=True)
    alt = _tc("a", [Exon(120, 150), Exon(200, 300)])
    events = derive_as_events([primary, alt])
    assert ("ALT_TSS", 100, 120) in _kinds_with_coords(events)


def test_alt_tes_plus():
    primary = _tc("p", [Exon(100, 150), Exon(200, 300)], primary=True)
    alt = _tc("a", [Exon(100, 150), Exon(200, 280)])
    events = derive_as_events([primary, alt])
    assert ("ALT_TES", 280, 300) in _kinds_with_coords(events)


# --- minus strand ---

def test_alt_tss_minus_uses_high_coord():
    primary = _tc("p", [Exon(100, 150), Exon(200, 300)], strand="-", primary=True)
    alt = _tc("a", [Exon(100, 150), Exon(200, 280)], strand="-")
    events = derive_as_events([primary, alt])
    # minus TSS is the high coordinate (transcript end differs 300 vs 280)
    assert ("ALT_TSS", 280, 300) in _kinds_with_coords(events)


def test_alt_tes_minus_uses_low_coord():
    primary = _tc("p", [Exon(100, 150), Exon(200, 300)], strand="-", primary=True)
    alt = _tc("a", [Exon(120, 150), Exon(200, 300)], strand="-")
    events = derive_as_events([primary, alt])
    # minus TES is the low coordinate (transcript start differs 100 vs 120)
    assert ("ALT_TES", 100, 120) in _kinds_with_coords(events)


def test_alt_site_minus_flips_class():
    primary = _tc("p", [Exon(100, 150), Exon(200, 300)], strand="-", primary=True)
    alt = _tc("a", [Exon(100, 150), Exon(230, 300)], strand="-")
    events = derive_as_events([primary, alt])
    # share donor (low) 150, differ high -> on minus, high boundary is 5' -> A5
    assert ("A5", 200, 230) in _kinds_with_coords(events)


def test_event_strand_recorded():
    primary = _tc("p", [Exon(100, 150), Exon(200, 280)], strand="-", primary=True)
    alt = _tc("a", [Exon(100, 150), Exon(200, 300)], strand="-")
    events = derive_as_events([primary, alt])
    assert all(e.strand == "-" for e in events)


# --- mutually exclusive exons ---

def test_mutually_exclusive_exons():
    # shared flanking exons (100,150) and (400,450); competing middle exons
    primary = _tc("p", [Exon(100, 150), Exon(200, 250), Exon(400, 450)], primary=True)
    alt = _tc("a", [Exon(100, 150), Exon(300, 350), Exon(400, 450)])
    events = derive_as_events([primary, alt])
    assert ("MX", 200, 350) in _kinds_with_coords(events)


# --- redundancy ---

def test_is_redundant_identical():
    a = _tc("a", [Exon(100, 300)])
    b = _tc("b", [Exon(100, 300)])
    assert is_redundant(a, b, 0.6, 0.6)


def test_is_redundant_disjoint_false():
    a = _tc("a", [Exon(100, 200)])
    b = _tc("b", [Exon(500, 600)])
    assert not is_redundant(a, b, 0.6, 0.6)


def test_is_redundant_cds_disagreement_false():
    a = _tc("a", [Exon(100, 400)], cds=[CDSSegment(100, 250, 0)])
    b = _tc("b", [Exon(100, 400)], cds=[CDSSegment(300, 399, 0)])
    # cDNA overlap high, but CDS do not overlap -> not redundant
    assert not is_redundant(a, b, 0.6, 0.6)


def test_is_redundant_cds_agreement_true():
    a = _tc("a", [Exon(100, 400)], cds=[CDSSegment(100, 250, 0)])
    b = _tc("b", [Exon(100, 400)], cds=[CDSSegment(100, 250, 0)])
    assert is_redundant(a, b, 0.6, 0.6)


def test_is_redundant_no_cds_uses_cdna():
    a = _tc("a", [Exon(100, 300)])
    b = _tc("b", [Exon(100, 300)])
    assert is_redundant(a, b, 0.6, 0.6)


# --- Phase 20 §2.5: linear (two-pointer) overlap_bases ---

def _brute_overlap(a_intervals, b_intervals):
    """Reference O(Ea×Eb) overlap (the pre-Phase-20 nested double loop)."""
    total = 0
    for a in a_intervals:
        for b in b_intervals:
            ov = min(a.end, b.end) - max(a.start, b.start)
            if ov > 0:
                total += ov
    return total


def _disjoint_intervals(rng, n, base):
    """A sorted list of n disjoint Exons starting at/after ``base``."""
    out = []
    cursor = base
    for _ in range(n):
        cursor += rng.randint(0, 40)        # inter-exon gap (>=0)
        length = rng.randint(1, 60)         # exon length (>=1)
        out.append(Exon(cursor, cursor + length))
        cursor += length
    return out


def test_overlap_bases_matches_bruteforce_randomized():
    import random
    rng = random.Random(20)
    for _ in range(200):
        a = _disjoint_intervals(rng, rng.randint(0, 8), 1000)
        b = _disjoint_intervals(rng, rng.randint(0, 8), 1000)
        assert overlap_bases(a, b) == _brute_overlap(a, b)


def test_overlap_bases_handbuilt_cases():
    # touching, not overlapping
    assert overlap_bases([Exon(100, 200)], [Exon(200, 300)]) == 0
    # fully nested: smaller inside larger
    assert overlap_bases([Exon(100, 500)], [Exon(200, 300)]) == 100
    # multi-interval partial overlaps, order-independent
    a = [Exon(100, 200), Exon(300, 400), Exon(500, 600)]
    b = [Exon(150, 350), Exon(550, 650)]
    expected = 50 + 50 + 50  # 150-200, 300-350, 550-600
    assert overlap_bases(a, b) == expected
    assert overlap_bases(b, a) == expected
    # identical chains -> full length
    assert overlap_bases(a, a) == 100 + 100 + 100
    # empty lists
    assert overlap_bases([], [Exon(100, 200)]) == 0
    assert overlap_bases([Exon(100, 200)], []) == 0


def test_overlap_bases_unsorted_input_defensive_sort():
    # Two-pointer sorts defensively, so an out-of-order list is still correct.
    a = [Exon(500, 600), Exon(100, 200)]
    b = [Exon(150, 550)]
    assert overlap_bases(a, b) == 50 + 50  # 150-200 and 500-550


def test_overlap_bases_minus_strand_intervals_match_bruteforce():
    # Coordinates are stored (low, high) regardless of strand; build the exon
    # lists of a minus-strand pair and confirm the linear result still matches.
    minus_a = [Exon(1000, 1200), Exon(1500, 1800), Exon(2000, 2100)]
    minus_b = [Exon(1100, 1550), Exon(1900, 2050)]
    assert overlap_bases(minus_a, minus_b) == _brute_overlap(minus_a, minus_b)
    # concrete: 1100-1200 (100) + 1500-1550 (50) + 2000-2050 (50)
    assert overlap_bases(minus_a, minus_b) == 200


def test_is_redundant_decisions_unchanged_both_strands():
    for strand in ("+", "-"):
        # identical structure + CDS -> redundant
        a = _tc("a", [Exon(100, 250), Exon(400, 600)], strand=strand,
                cds=[CDSSegment(151, 250, 0), CDSSegment(400, 499, 0)])
        b = _tc("b", [Exon(100, 250), Exon(400, 600)], strand=strand,
                cds=[CDSSegment(151, 250, 0), CDSSegment(400, 499, 0)])
        assert is_redundant(a, b, 0.6, 0.6)
        # different intron chain (extra exon) -> not redundant even on minus
        c = _tc("c", [Exon(100, 250), Exon(300, 350), Exon(400, 600)], strand=strand)
        d = _tc("d", [Exon(100, 250), Exon(400, 600)], strand=strand)
        assert not is_redundant(c, d, 0.6, 0.6)


def test_is_redundant_precomputed_chain_matches_inline():
    from helixforge.reconcile.as_events import _intron_chain
    a = _tc("a", [Exon(100, 250), Exon(400, 600)])
    b = _tc("b", [Exon(100, 250), Exon(400, 600)])
    inline = is_redundant(a, b, 0.6, 0.6)
    cached = is_redundant(a, b, 0.6, 0.6,
                          a_chain=_intron_chain(a), b_chain=_intron_chain(b))
    assert inline == cached is True
