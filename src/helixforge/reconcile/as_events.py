"""Alternative-splicing event taxonomy."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from helixforge.reconcile.models import ASEvent, TranscriptCandidate
from helixforge.utils.intervals import bounds as _bounds


def overlap_bases(a_intervals: Sequence[Any], b_intervals: Sequence[Any]) -> int:
    """Total overlapping bases between two lists of disjoint intervals.

    Two-pointer sort-merge over the (disjoint) interval lists: O(Ea+Eb) after an
    O(E log E) defensive sort, replacing the old O(Ea×Eb) nested loop (Phase 20
    §2.5). The integer result is identical for disjoint inputs (the real domain,
    exon/CDS lists are sorted, non-overlapping by model invariant), which the
    locus-matching inner loop and isoform-redundancy filter both consume.
    """
    a = sorted(_bounds(x) for x in a_intervals)
    b = sorted(_bounds(x) for x in b_intervals)
    i = j = 0
    total = 0
    while i < len(a) and j < len(b):
        sa, ea = a[i]
        sb, eb = b[j]
        ov = min(ea, eb) - max(sa, sb)
        if ov > 0:
            total += ov
        # Advance the interval that ends first; ties advance b (then a's tail
        # cannot overlap b's next interval, so a advances next iteration).
        if ea < eb:
            i += 1
        else:
            j += 1
    return total


def _total_length(intervals: Sequence[Any]) -> int:
    return sum(e - s for s, e in (_bounds(i) for i in intervals))


def reciprocal_overlap(a_intervals: Sequence[Any], b_intervals: Sequence[Any]) -> float:
    """Reciprocal overlap fraction ``min(ov/len_a, ov/len_b)`` (0 if either empty)."""
    la = _total_length(a_intervals)
    lb = _total_length(b_intervals)
    if la == 0 or lb == 0:
        return 0.0
    ov = overlap_bases(a_intervals, b_intervals)
    return min(ov / la, ov / lb)


def _intron_chain(t: TranscriptCandidate) -> frozenset[tuple[int, int]]:
    return frozenset((i.start, i.end) for i in t.introns)


def is_redundant(
    a: TranscriptCandidate,
    b: TranscriptCandidate,
    min_cds_overlap: float,
    min_cdna_overlap: float,
    *,
    a_chain: frozenset[tuple[int, int]] | None = None,
    b_chain: frozenset[tuple[int, int]] | None = None,
) -> bool:
    """True if isoforms ``a`` and ``b`` are structurally near-duplicate.

    Requires high cDNA (exon) reciprocal overlap AND **identical intron chains**
    AND, when both carry CDS, CDS reciprocal overlap ``>= min_cds_overlap``. The
    intron-chain check is what keeps genuine alternatives (ES/IR/A5/A3/MX), which
    can share a lot of cDNA, from being dropped. Defense in depth: Mikado pick
    should already remove most redundancy.

    ``a_chain``/``b_chain`` accept a precomputed ``frozenset`` intron chain so a
    pairwise sweep (``_filter_redundant``) computes each chain once instead of
    per pair (Phase 20 §2.5); omitted, they fall back to computing it here.
    """
    if reciprocal_overlap(a.exons, b.exons) < min_cdna_overlap:
        return False
    ac = a_chain if a_chain is not None else _intron_chain(a)
    bc = b_chain if b_chain is not None else _intron_chain(b)
    if ac != bc:
        return False
    if a.cds and b.cds:
        return reciprocal_overlap(a.cds, b.cds) >= min_cds_overlap
    return True


def derive_as_events(transcripts: list[TranscriptCandidate]) -> list[ASEvent]:
    """Compute the AS-event set for a gene's isoforms (relative to the primary)."""
    if len(transcripts) < 2:
        return []

    primary = next((t for t in transcripts if t.is_primary), transcripts[0])
    seqid, strand = primary.seqid, primary.strand
    events: list[ASEvent] = []
    seen: set[tuple[str, int, int]] = set()

    def add(kind: str, start: int, end: int) -> None:
        if start >= end:
            return
        key = (kind, start, end)
        if key in seen:
            return
        seen.add(key)
        events.append(ASEvent(kind, seqid, start, end, strand))

    p_exons = [_bounds(e) for e in primary.exons]
    p_introns = [_bounds(i) for i in primary.introns]
    p_exon_set = set(p_exons)
    p_intron_set = set(p_introns)

    for q in transcripts:
        if q is primary:
            continue
        q_exons = [_bounds(e) for e in q.exons]
        q_introns = [_bounds(i) for i in q.introns]
        q_exon_set = set(q_exons)
        q_intron_set = set(q_introns)

        _detect_skips(add, p_exons, q_exons, q_exon_set, q_introns)
        _detect_skips(add, q_exons, p_exons, p_exon_set, p_introns)
        _detect_retention(add, p_introns, q_intron_set, q_exons)
        _detect_retention(add, q_introns, p_intron_set, p_exons)
        _detect_alt_sites(add, p_introns, q_introns, strand)
        _detect_termini(add, primary, q, strand)
        _detect_mutually_exclusive(add, p_exons, q_exons, p_exon_set, q_exon_set)

    return events


def _detect_skips(
    add: Callable[[str, int, int], None],
    exons: list[tuple[int, int]],
    other_exons: list[tuple[int, int]],
    other_exon_set: set[tuple[int, int]],
    other_introns: list[tuple[int, int]],
) -> None:
    """Exon in ``exons`` spliced over by ``other`` (an intron spans it, no exon)."""
    for ex in exons:
        if ex in other_exon_set:
            continue
        spanned = any(i[0] <= ex[0] and i[1] >= ex[1] for i in other_introns)
        has_exon_here = any(oe[0] < ex[1] and ex[0] < oe[1] for oe in other_exons)
        if spanned and not has_exon_here:
            add("ES", ex[0], ex[1])


def _detect_retention(
    add: Callable[[str, int, int], None],
    introns: list[tuple[int, int]],
    other_intron_set: set[tuple[int, int]],
    other_exons: list[Any],
) -> None:
    """Intron retained as a single exon in the other isoform."""
    for i in introns:
        if i in other_intron_set:
            continue
        if any(
            oe[0] <= i[0] and oe[1] >= i[1] for oe in (_bounds(e) for e in other_exons)
        ):
            add("IR", i[0], i[1])


def _detect_alt_sites(
    add: Callable[[str, int, int], None],
    p_introns: list[tuple[int, int]],
    q_introns: list[tuple[int, int]],
    strand: str,
) -> None:
    """Introns sharing exactly one boundary -> alt 5'/3' splice site."""
    for ip in p_introns:
        for iq in q_introns:
            if ip == iq:
                continue
            share_low = ip[0] == iq[0]
            share_high = ip[1] == iq[1]
            if share_low and not share_high:
                lo, hi = sorted((ip[1], iq[1]))
                add("A3" if strand == "+" else "A5", lo, hi)
            elif share_high and not share_low:
                lo, hi = sorted((ip[0], iq[0]))
                add("A5" if strand == "+" else "A3", lo, hi)


def _detect_termini(
    add: Callable[[str, int, int], None],
    primary: TranscriptCandidate,
    q: TranscriptCandidate,
    strand: str,
) -> None:
    """Differing transcript termini -> ALT_TSS / ALT_TES (strand-aware)."""
    p_lo, p_hi = primary.start, primary.end
    q_lo, q_hi = q.start, q.end
    if strand == "+":
        if p_lo != q_lo:
            add("ALT_TSS", *sorted((p_lo, q_lo)))
        if p_hi != q_hi:
            add("ALT_TES", *sorted((p_hi, q_hi)))
    else:
        if p_hi != q_hi:
            add("ALT_TSS", *sorted((p_hi, q_hi)))
        if p_lo != q_lo:
            add("ALT_TES", *sorted((p_lo, q_lo)))


def _detect_mutually_exclusive(
    add: Callable[[str, int, int], None],
    p_exons: list[tuple[int, int]],
    q_exons: list[tuple[int, int]],
    p_exon_set: set[tuple[int, int]],
    q_exon_set: set[tuple[int, int]],
) -> None:
    """Two isoforms use different, non-overlapping exons in the same slot."""
    shared = sorted(p_exon_set & q_exon_set)
    for k in range(len(shared) - 1):
        gap_lo, gap_hi = shared[k][1], shared[k + 1][0]
        p_in = [
            e
            for e in p_exons
            if e not in q_exon_set and e[0] >= gap_lo and e[1] <= gap_hi
        ]
        q_in = [
            e
            for e in q_exons
            if e not in p_exon_set and e[0] >= gap_lo and e[1] <= gap_hi
        ]
        for ep in p_in:
            for eq in q_in:
                if not (
                    ep[0] < eq[1] and eq[0] < ep[1]
                ):  # non-overlapping alternatives
                    add("MX", min(ep[0], eq[0]), max(ep[1], eq[1]))
