"""Reconciliation onto the Helixer gene set."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import attrs

from helixforge.constants import (
    MERGE_MIN_GAP_READS,
    NOVEL_ID_BASE,
    PARALOG_IDENTITY_K,
)
from helixforge.qc.flags import (
    FROM_TANGLED_LOCUS,
    HELIXER_ONLY,
    LOCUS_MERGE,
    LOCUS_SPLIT,
    MERGE_REJECTED,
    NO_EXPRESSION,
    NOVEL_LOCUS,
    PARTIAL_ORF,
    as_event_flag,
    dedup_flags,
)
from helixforge.reconcile.as_events import (
    _intron_chain,
    derive_as_events,
    is_redundant,
    overlap_bases,
    reciprocal_overlap,
)
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    HelixerLocus,
    Interval,
    IsoformAdmission,
    LocusClassification,
    MikadoLocus,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.utils.intervals import IntervalIndex, bounds
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)

# NOVEL_ID_BASE re-exported from helixforge.constants (Phase 19, §1.6): novel
# (Mikado-only) genes get HFGs from a high range, outside the Helixer-anchored
# block, so they never collide with Helixer-derived ids.
MIKADO_ORIGINS = ("mikado_1to1", "split", "merge", "novel")
_SUFFIXES = "abcdefghijklmnopqrstuvwxyz"


# ---------------------------------------------------------------------------
# ID allocation
# ---------------------------------------------------------------------------


def _hfg_num(hfg: str) -> int:
    return int(hfg.split("_")[1])


def _fmt_hfg(n: int) -> str:
    return f"HFG_{n:05d}"


@dataclass
class IdAllocator:
    """Stateful allocator for the next free HFG number (Phase 16 §C3, Phase 20 §1.1).

    Lets a parallel chunk draw new Helixer-anchored numbers from a reserved,
    disjoint ``[base, ...)`` range and novel numbers from ``[novel_base, ...)``
    instead of the single global block, so chunks never need a run-time lock to
    stay collision-free. ``base``/``novel_base`` default to the historical
    ``1``/``NOVEL_ID_BASE``, so a single-process run is byte-for-byte unchanged
    and every existing test + the M3 golden counts hold.

    **Allocation is amortized O(1)** (Phase 20): the historical scanner rebuilt a
    ``used`` set over the *entire* ``id_map`` on every call, Θ(N²) over a fresh
    whole-genome run. This version seeds ``_used`` and the two cursors **once**,
    lazily, from the first ``id_map`` it sees, then advances a monotonically
    rising cursor past ``_used`` and records each number it hands out. The
    "lowest free number at/above base" semantics are preserved exactly: because
    ``id_map`` only ever grows (numbers are never freed) and every fresh number
    is recorded, the cursor reproduces the scanner's strictly-increasing
    sequence. So HFG assignments are identical and M3 IDs do not move.

    The allocator owns mutable state, so a *fresh* instance is used per
    reconcile run / per standalone ``assign_*`` call (callers do
    ``allocator or IdAllocator()``); it must not be shared as a long-lived
    singleton across independent runs (the cursor would not reset).
    """

    base: int = 1
    novel_base: int = NOVEL_ID_BASE
    # Allocation state, seeded lazily on the first next_number call (so the
    # public constructor signature stays base/novel_base only). None => unseeded.
    _used: set[int] | None = field(default=None, init=False, repr=False)
    _cursor: int = field(default=0, init=False, repr=False)
    _novel_cursor: int = field(default=0, init=False, repr=False)

    def _ensure_seeded(self, id_map: dict[str, str]) -> None:
        """Seed ``_used`` + cursors once from ``id_map`` (idempotent after first)."""
        if self._used is not None:
            return
        used: set[int] = set()
        for value in id_map.values():
            try:
                used.add(_hfg_num(value))
            except (IndexError, ValueError):
                continue
        self._used = used
        self._cursor = self.base
        self._novel_cursor = self.novel_base

    def next_number(self, id_map: dict[str, str], novel: bool = False) -> int:
        """Lowest free HFG number at/above the policy's (novel) base, in O(1) amortized.

        On the first call ``id_map`` seeds the allocation state (so a chunk whose
        seed map carries a foreign/out-of-range HFG, e.g. a locus that moved
        chunks across a re-partition, still skips it, exactly like the historical
        scanner). Subsequent calls do not rescan ``id_map``; they advance the
        cursor past the recorded ``_used`` set and record the returned number.
        """
        self._ensure_seeded(id_map)
        assert self._used is not None  # set by _ensure_seeded
        if novel:
            n = self._novel_cursor
            while n in self._used:
                n += 1
            self._used.add(n)
            self._novel_cursor = n + 1
        else:
            n = self._cursor
            while n in self._used:
                n += 1
            self._used.add(n)
            self._cursor = n + 1
        return n


# The historical global policy template (base=1, novel_base=90000). Kept as a
# defined symbol for back-compat; it is NOT used as a shared mutable singleton,
# callers allocate a fresh IdAllocator() per run/call so state never leaks
# across independent runs (Phase 20 §1.1).
_DEFAULT_ALLOCATOR = IdAllocator()


def _next_number(id_map: dict[str, str], novel: bool = False) -> int:
    """Back-compat shim: lowest free HFG number from the default global block.

    Uses a fresh allocator each call so the result depends only on ``id_map``
    (matching the historical stateless scanner for one-off callers)."""
    return IdAllocator().next_number(id_map, novel=novel)


def assign_gene_id(
    helixer_locus: HelixerLocus,
    id_map: dict[str, str],
    merged_from: Iterable[str] = (),
    *,
    allocator: IdAllocator | None = None,
) -> str:
    """Stable HFG for a Helixer locus; merges keep the lowest existing HFG.

    Every constituent Helixer id (the locus + ``merged_from``) is mapped to the
    chosen HFG in ``id_map`` so future runs reuse it. Never renumbers. New
    numbers come from ``allocator`` (default: the global policy).
    """
    alloc = allocator or IdAllocator()
    ids = [helixer_locus.gene_id, *merged_from]
    existing = [id_map[i] for i in ids if i in id_map]
    chosen = (
        min(existing, key=_hfg_num) if existing else _fmt_hfg(alloc.next_number(id_map))
    )
    for i in ids:
        id_map[i] = chosen
    return chosen


def assign_novel_gene_id(
    mikado_locus: MikadoLocus,
    id_map: dict[str, str],
    *,
    allocator: IdAllocator | None = None,
) -> str:
    """Fresh HFG (novel range) for a Mikado-only locus, keyed by its locus id."""
    alloc = allocator or IdAllocator()
    key = mikado_locus.locus_id
    if key in id_map:
        return id_map[key]
    hfg = _fmt_hfg(alloc.next_number(id_map, novel=True))
    id_map[key] = hfg
    return hfg


# ---------------------------------------------------------------------------
# Structural comparison + correspondence
# ---------------------------------------------------------------------------


def intron_chain(t: TranscriptCandidate) -> tuple[tuple[int, int], ...]:
    """The transcript's intron chain as a tuple of ``(start, end)`` pairs."""
    return tuple((i.start, i.end) for i in t.introns)


def reciprocal_cds_overlap(a: TranscriptCandidate, b: TranscriptCandidate) -> float:
    """Reciprocal CDS overlap fraction (0.0 if either lacks CDS)."""
    if a.cds and b.cds:
        return reciprocal_overlap(a.cds, b.cds)
    return 0.0


def _containment_overlap(a_intervals: list[Any], b_intervals: list[Any]) -> float:
    """Fraction of the *smaller* structure covered by the other (``ov/min(la,lb)``).

    Used for locus-correspondence edges: containment (not strict reciprocal min)
    is what lets a large Mikado locus form an edge with each of the small Helixer
    loci it splits from / merges, while still ignoring raw span (CDS/exon based,
    so UTRs do not drive false merges). The strict reciprocal metric lives in
    :func:`reciprocal_cds_overlap` for 1:1 confidence and in ``is_redundant``.
    """
    la = sum(e - s for s, e in (bounds(i) for i in a_intervals))
    lb = sum(e - s for s, e in (bounds(i) for i in b_intervals))
    if la == 0 or lb == 0:
        return 0.0
    return overlap_bases(a_intervals, b_intervals) / min(la, lb)


def _locus_transcript_overlap(
    helixer_locus: HelixerLocus, transcript: TranscriptCandidate
) -> float:
    """Best structural overlap: CDS if both have it, else exon overlap."""
    if helixer_locus.cds and transcript.cds:
        return _containment_overlap(helixer_locus.cds, transcript.cds)
    return _containment_overlap(helixer_locus.exons, transcript.exons)


def _structural_overlap(
    helixer_locus: HelixerLocus, mikado_locus: MikadoLocus
) -> float:
    return max(
        (_locus_transcript_overlap(helixer_locus, t) for t in mikado_locus.transcripts),
        default=0.0,
    )


@dataclass
class Correspondence:
    """The locus-correspondence partition produced by :func:`match_loci`."""

    one_to_one: list[tuple[HelixerLocus, MikadoLocus]] = field(default_factory=list)
    splits: list[tuple[HelixerLocus, list[MikadoLocus]]] = field(default_factory=list)
    merges: list[tuple[list[HelixerLocus], MikadoLocus]] = field(default_factory=list)
    helixer_only: list[HelixerLocus] = field(default_factory=list)
    novel: list[MikadoLocus] = field(default_factory=list)
    # locus_ids of novel loci that came from a tangled many-to-many component
    # (a Mikado locus that lost all its Helixer claimants). Flagged
    # FROM_TANGLED_LOCUS when admitted; never silently dropped.
    tangled_loci: set[str] = field(default_factory=set)


class _UnionFind:
    def __init__(self, nodes: list[tuple[str, int]]) -> None:
        self.parent: dict[tuple[str, int], tuple[str, int]] = {n: n for n in nodes}

    def find(self, x: tuple[str, int]) -> tuple[str, int]:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: tuple[str, int], b: tuple[str, int]) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def match_loci(
    helixer_loci: list[HelixerLocus],
    mikado_loci: list[MikadoLocus],
    reciprocal_overlap: float = 0.5,
) -> Correspondence:
    """Partition Helixer↔Mikado loci into the five correspondence cases.

    An edge is drawn when loci share seqid + strand and their structural
    reciprocal overlap (CDS if available, else exons) is ``>= reciprocal_overlap``.
    Connected components are then classified by their (helixer, mikado) counts.
    """
    # per-scaffold index over mikado loci spans
    mikado_by_scaffold: dict[str, list[tuple[int, int, int]]] = {}
    for mi, m in enumerate(mikado_loci):
        mikado_by_scaffold.setdefault(m.seqid, []).append((m.start, m.end, mi))
    indices: dict[str, IntervalIndex] = {}
    for seqid, intervals in mikado_by_scaffold.items():
        idx = IntervalIndex()
        idx.add_intervals(intervals)
        indices[seqid] = idx

    nodes = [("h", i) for i in range(len(helixer_loci))]
    nodes += [("m", j) for j in range(len(mikado_loci))]
    uf = _UnionFind(nodes)
    h_has_edge = [False] * len(helixer_loci)
    m_has_edge = [False] * len(mikado_loci)

    for hi, h in enumerate(helixer_loci):
        hidx = indices.get(h.seqid)
        if hidx is None:
            continue
        for _s, _e, mi in hidx.query_with_data(h.start, h.end):
            m = mikado_loci[mi]
            if m.strand != h.strand:
                continue
            if _structural_overlap(h, m) >= reciprocal_overlap:
                uf.union(("h", hi), ("m", mi))
                h_has_edge[hi] = True
                m_has_edge[mi] = True

    # group nodes by component
    components: dict[tuple[str, int], list[tuple[str, int]]] = {}
    for node in nodes:
        components.setdefault(uf.find(node), []).append(node)

    corr = Correspondence()
    for members in components.values():
        hs = sorted(
            (helixer_loci[i] for kind, i in members if kind == "h"),
            key=lambda g: g.start,
        )
        ms = sorted(
            (mikado_loci[j] for kind, j in members if kind == "m"),
            key=lambda g: g.start,
        )
        if not ms:
            corr.helixer_only.extend(hs)
        elif not hs:
            corr.novel.extend(ms)
        elif len(hs) == 1 and len(ms) == 1:
            corr.one_to_one.append((hs[0], ms[0]))
        elif len(hs) == 1 and len(ms) >= 2:
            corr.splits.append((hs[0], ms))
        elif len(hs) >= 2 and len(ms) == 1:
            corr.merges.append((hs, ms[0]))
        else:
            _decompose_many_to_many(corr, hs, ms)

    return corr


def _decompose_many_to_many(
    corr: Correspondence, hs: list[HelixerLocus], ms: list[MikadoLocus]
) -> None:
    """Best-effort split of a many-to-many component (assign each Helixer to its
    best-overlapping Mikado; group into 1:1 or merge). Keeps every Helixer.

    A Mikado locus that loses *all* its Helixer claimants is **not** dropped:
    in tangled plant tandem arrays it carries real structure, so it
    is retained as a novel-style locus (subject to the ``admit_novel`` gate and
    flagged FROM_TANGLED_LOCUS). When not admitted it still earns an audit record
    via the normal novel path, never silently lost."""
    # Memoize structural overlap per (helixer, mikado) pair so it is computed
    # once, not recomputed in the max(...) and any later regrouping (Phase 20 §2.5).
    overlap_cache: dict[tuple[int, int], float] = {}

    def s_overlap(h: HelixerLocus, m: MikadoLocus) -> float:
        key = (id(h), id(m))
        if key not in overlap_cache:
            overlap_cache[key] = _structural_overlap(h, m)
        return overlap_cache[key]

    by_mikado: dict[int, list[HelixerLocus]] = {}
    for h in hs:
        best = max(ms, key=lambda m: s_overlap(h, m))
        by_mikado.setdefault(id(best), []).append(h)
    for m in ms:
        group = by_mikado.get(id(m), [])
        if len(group) == 1:
            corr.one_to_one.append((group[0], m))
        elif len(group) >= 2:
            corr.merges.append((sorted(group, key=lambda g: g.start), m))
        else:
            _log.info(
                "many-to-many: Mikado locus %s received no Helixer; retained as "
                "novel-style",
                m.locus_id,
            )
            corr.novel.append(m)
            corr.tangled_loci.add(m.locus_id)


# ---------------------------------------------------------------------------
# Bridging-junction check for merges
# ---------------------------------------------------------------------------


def _has_bridging_introns(
    helixer_loci_sorted: list[HelixerLocus],
    mikado_locus: MikadoLocus,
    junctions: Any = None,
    *,
    min_gap_reads: int = MERGE_MIN_GAP_READS,
    require_canonical: bool = True,
    tolerance: int = 0,
) -> bool:
    """True if every inter-locus gap is spanned by a *verified* bridging intron.

    The base requirement is unchanged: a Mikado intron (Mikado uses
    ``only_confirmed_introns``) must cross the gap between each pair of adjacent
    Helixer loci. The guard is strengthened against fusing adjacent paralogs in
    tandem arrays when a ``junctions`` set/index is supplied: the bridging intron
    must also correspond to a splice junction that

    * is **canonical**, a non-canonical bridge is rejected, and
    * is supported by **>= ``min_gap_reads`` reads** spanning the inter-genic gap
      specifically (not merely *some* intron elsewhere in the Mikado locus).

    ``junctions`` may be a ``list[SpliceJunction]`` or a :class:`JunctionIndex`.
    When it is ``None`` the function is byte-for-byte the legacy any-Mikado-intron
    check, so existing callers (and the model tests) are unchanged. A junction
    whose motif was never evaluated (``canonical is None``) is treated as "not
    known non-canonical" and accepted (the Phase 28 convention), so a junction set
    built without motifs degrades to the read-count guard alone.
    """
    # Local import keeps the module-load graph acyclic (fallback imports models).
    from helixforge.reconcile.fallback import find_matching_junction

    all_introns = [
        (i.start, i.end) for t in mikado_locus.transcripts for i in t.introns
    ]
    for a, b in zip(helixer_loci_sorted, helixer_loci_sorted[1:]):
        gap_lo, gap_hi = a.end, b.start
        if gap_lo >= gap_hi:  # overlapping loci shouldn't reach here
            continue
        bridging = [(s, e) for (s, e) in all_introns if s <= gap_lo and e >= gap_hi]
        if not bridging:
            return False
        if junctions is None:
            continue
        # Require at least one bridging intron backed by a sufficiently-supported,
        # canonical junction spanning the gap.
        ok = False
        for s, e in bridging:
            j = find_matching_junction(
                Interval(s, e),
                junctions,
                mikado_locus.seqid,
                mikado_locus.strand,
                tolerance,
                min_gap_reads,
            )
            if j is None:
                continue
            if require_canonical and j.canonical is not None and not j.is_canonical:
                continue
            ok = True
            break
        if not ok:
            return False
    return True


# ---------------------------------------------------------------------------
# Paralog / recent-duplication identity proxy (alignment-free)
# ---------------------------------------------------------------------------


def _kmer_set(seq: str, k: int) -> frozenset[str]:
    """The set of length-``k`` substrings of ``seq`` (upper-cased)."""
    s = seq.upper()
    if len(s) < k:
        return frozenset()
    return frozenset(s[i : i + k] for i in range(len(s) - k + 1))


def kmer_identity(seq_a: str, seq_b: str, k: int = PARALOG_IDENTITY_K) -> float:
    """Alignment-free pairwise identity proxy: k-mer Jaccard of two sequences.

    A cheap recent-duplication / homeolog signal, near-identical
    adjacent loci yield a Jaccard near 1.0. Returns 0.0 when either sequence is
    shorter than ``k`` (no shared k-mers possible). Order-independent and
    strand-naive: the caller passes coding-direction CDS/exon sequence.
    """
    ka, kb = _kmer_set(seq_a, k), _kmer_set(seq_b, k)
    if not ka or not kb:
        return 0.0
    inter = len(ka & kb)
    union = len(ka | kb)
    return inter / union if union else 0.0


def _locus_cds_sequence(locus: HelixerLocus, genome: Any) -> str | None:
    """Coding-direction CDS sequence for a Helixer locus (exon fallback).

    Used only by the optional paralog-identity guard. Returns ``None`` when no
    genome accessor is available or the sequence cannot be fetched.
    """
    if genome is None:
        return None
    from helixforge.utils.sequences import reverse_complement

    segs = locus.cds if locus.cds else locus.exons
    if not segs:
        return None
    try:
        parts = [genome.get_sequence(locus.seqid, s.start, s.end) for s in segs]
    except Exception:  # noqa: BLE001 - a fetch failure disables the guard, never fatal
        return None
    seq = "".join(parts)
    return reverse_complement(seq) if locus.strand == "-" else seq


def _high_identity_pair(
    helixer_loci_sorted: list[HelixerLocus],
    genome: Any,
    threshold: float,
    k: int = PARALOG_IDENTITY_K,
) -> bool:
    """True if any adjacent pair of loci has k-mer identity >= ``threshold``.

    A positive result marks the component as a recent-duplication / homeolog
    cluster whose merge should be rejected (the two loci are paralogs, not one
    split gene). Disabled (returns False) when ``genome`` is None or a sequence
    cannot be built.
    """
    seqs = [_locus_cds_sequence(h, genome) for h in helixer_loci_sorted]
    for sa, sb in zip(seqs, seqs[1:]):
        if sa is None or sb is None:
            continue
        if kmer_identity(sa, sb, k) >= threshold:
            return True
    return False


# ---------------------------------------------------------------------------
# Tier assignment
# ---------------------------------------------------------------------------


def assign_tier(
    primary: TranscriptCandidate,
    classification: LocusClassification | None,
    origin: str,
) -> int:
    """Assign quality tier 1–4.

    Tier 1 = coherent CDS + protein-homology support. Homology is a hit accession
    (``protein_id``) OR a positive BLAST/DIAMOND score (``blast_score`` from the
    metrics TSV), Mikado's loci GFF3 carries no accession, so the BLAST score is
    the real signal.

    A ``helixer_backstop`` gene that has been **rescued** with a valid projected
    CDS tiers like a Mikado-origin gene: Tier 1 if homology-backed
    (a real miniprot accession), Tier 2 if CDS-only. A still-silent backstop with
    no CDS stays at Tier 3 (expressed/low) or Tier 4 (silent). M3 is unchanged
    because its 832 backstop genes remain CDS-less (no miniprot GFF)."""
    if origin in MIKADO_ORIGINS:
        if primary.cds is not None and primary.has_homology:
            return 1
        return 2
    # helixer_backstop
    if primary.cds is not None:
        return 1 if primary.has_homology else 2
    if classification is not None and classification.status in ("EXPRESSED", "LOW"):
        return 3
    return 4


# ---------------------------------------------------------------------------
# Gene construction
# ---------------------------------------------------------------------------


def _renumber(
    transcripts: list[TranscriptCandidate],
    gene_id: str,
    order: list[str] | None = None,
) -> list[TranscriptCandidate]:
    """Order isoforms and re-id as ``gene_id.N`` (primary = .1).

    Default ordering is by Mikado ``combined_score`` (desc), tie-broken by
    ``transcript_id``, the historical behavior every existing caller and test
    relies on, so omitting ``order`` is byte-for-byte unchanged. When ``order`` is
    supplied (a list of ``transcript_id`` strings, e.g. the TRaCE election result;
    Phase 33b D2), the isoforms are numbered in *that* order instead of re-sorting
    by score; any transcript not named in ``order`` is appended afterward in the
    default score order so nothing is dropped.
    """

    def sort_key(t: TranscriptCandidate) -> tuple[float, str]:
        score = t.combined_score if t.combined_score is not None else float("-inf")
        return (-score, t.transcript_id)

    if order is None:
        ordered = sorted(transcripts, key=sort_key)
    else:
        rank = {tid: i for i, tid in enumerate(order)}
        # Named transcripts first in the supplied order; any unnamed remainder in
        # the default score order (rank == len(order) keeps them stable + last).
        ordered = sorted(
            transcripts,
            key=lambda t: (rank.get(t.transcript_id, len(order)), *sort_key(t)),
        )

    out: list[TranscriptCandidate] = []
    for i, t in enumerate(ordered):
        out.append(
            attrs.evolve(
                t,
                transcript_id=f"{gene_id}.{i + 1}",
                locus_id=gene_id,
                is_primary=(i == 0),
            )
        )
    return out


def _intrinsic_helixer_cds(
    h: HelixerLocus, exons: list[Exon]
) -> tuple[list[CDSSegment] | None, bool]:
    """The Helixer model's own CDS for a backstop transcript, plus a partial flag.

    A Helixer-only (backstop) gene carries the CDS Helixer itself predicted, this
    is the model's *intrinsic* ORF, and a gene must be called coding on the
    strength of that ORF, not only when an external miniprot/TransDecoder backstop
    re-derives one. The CDS is taken verbatim from the Helixer GFF3 (already
    sorted, within-exon, GFF3-phased), never adjusted (CDS boundaries come from
    the model, we only validate). Returns ``(None, False)`` when Helixer predicted
    no CDS (a genuine non-coding prediction) or the CDS is not contained in the
    transcript's exons (e.g. after a locus merge dropped the CDS). The partial
    flag is set when the total CDS length is not a multiple of 3 (a 5'/3'-partial
    ORF at a contig edge); it relaxes the mod-3 model invariant and surfaces as
    ``PARTIAL_ORF``.
    """
    if not h.cds:
        return None, False
    # Every CDS segment must be contained within a single exon, else the model
    # validator would reject it; leave the gene CDS-less in that (rare) case.
    for seg in h.cds:
        if not any(seg.start >= ex.start and seg.end <= ex.end for ex in exons):
            return None, False
    cds = [CDSSegment(c.start, c.end, c.phase) for c in h.cds]
    total = sum(len(c) for c in cds)
    return cds, (total % 3 != 0)


def build_reconciled_gene(
    helixer_locus_or_ids: HelixerLocus | None,
    mikado_locus_or_none: MikadoLocus | None,
    classification: LocusClassification,
    origin: str,
    id_map: dict[str, str],
    *,
    gene_id: str | None = None,
    merged_from: Iterable[str] = (),
    extra_flags: Iterable[Any] = (),
    allocator: IdAllocator | None = None,
) -> ReconciledGene:
    """Build one ``ReconciledGene``: assign id, order isoforms, derive AS events,
    assign tier, attach AS-event + locus flags."""
    if gene_id is None:
        if origin == "novel":
            assert (
                mikado_locus_or_none is not None
            )  # novel origin always has a mikado locus
            gene_id = assign_novel_gene_id(
                mikado_locus_or_none, id_map, allocator=allocator
            )
        else:
            assert (
                helixer_locus_or_ids is not None
            )  # non-novel origin always has a helixer locus
            gene_id = assign_gene_id(
                helixer_locus_or_ids, id_map, merged_from, allocator=allocator
            )

    if mikado_locus_or_none is not None:
        transcripts = _renumber(mikado_locus_or_none.transcripts, gene_id)
    else:
        assert (
            helixer_locus_or_ids is not None
        )  # no mikado locus means helixer locus is present
        h = helixer_locus_or_ids
        exons = list(h.exons) if h.exons else [Exon(h.start, h.end)]
        # Carry the Helixer model's own CDS so the gene is coding on the strength
        # of its intrinsic ORF: not only when a miniprot/TransDecoder backstop
        # re-derives one. Without this every Helixer-only gene loses its ORF here
        # and is mis-called non-coding downstream.
        cds, cds_partial = _intrinsic_helixer_cds(h, exons)
        transcripts = [
            TranscriptCandidate(
                transcript_id=f"{gene_id}.1",
                locus_id=gene_id,
                source="helixer",
                seqid=h.seqid,
                start=h.start,
                end=h.end,
                strand=h.strand,
                exons=exons,
                cds=cds,
                cds_partial=cds_partial,
                is_primary=True,
            )
        ]

    primary = transcripts[0]
    seqid = primary.seqid
    strand = primary.strand
    start = min(t.start for t in transcripts)
    end = max(t.end for t in transcripts)

    as_events = derive_as_events(transcripts)
    tier = assign_tier(primary, classification, origin)

    flags = list(extra_flags)
    if origin == "helixer_backstop":
        flags.append(HELIXER_ONLY)
        if classification is not None and classification.status == "SILENT":
            flags.append(NO_EXPRESSION)
    if getattr(primary, "cds_partial", False):
        flags.append(PARTIAL_ORF)
    for kind in sorted({e.kind for e in as_events}):
        flags.append(as_event_flag(kind))

    return ReconciledGene(
        gene_id=gene_id,
        seqid=seqid,
        start=start,
        end=end,
        strand=strand,
        tier=tier,
        transcripts=transcripts,
        primary_transcript_id=primary.transcript_id,
        classification=classification,
        origin=origin,
        as_events=as_events,
        flags=dedup_flags(flags),
        merged_from=list(merged_from),
    )


def record_admissions(
    gene: ReconciledGene,
    dropped: Iterable[tuple[str, str, float | None]] = (),
) -> list[IsoformAdmission]:
    """One ``IsoformAdmission`` per retained isoform (+ any dropped redundants)."""
    primary = next((t for t in gene.transcripts if t.is_primary), gene.transcripts[0])
    out = []
    for t in gene.transcripts:
        if t.is_primary:
            out.append(
                IsoformAdmission(
                    t.transcript_id,
                    gene.gene_id,
                    True,
                    "primary",
                    score=t.combined_score,
                )
            )
        else:
            events = derive_as_events([primary, t])
            out.append(
                IsoformAdmission(
                    t.transcript_id,
                    gene.gene_id,
                    True,
                    "alternative_splicing",
                    novel_as_event=events[0] if events else None,
                    score=t.combined_score,
                )
            )
    for tid, redundant_with, score in dropped:
        out.append(
            IsoformAdmission(
                tid,
                gene.gene_id,
                False,
                "redundant",
                redundant_with=redundant_with,
                score=score,
            )
        )
    return out


def _filter_redundant(
    transcripts: list[TranscriptCandidate],
    min_cds_overlap: float,
    min_cdna_overlap: float,
) -> tuple[list[TranscriptCandidate], list[tuple[str, str, float | None]]]:
    """Drop near-duplicate isoforms (keep higher-scored); return (kept, dropped)."""

    def score(t: TranscriptCandidate) -> float:
        return t.combined_score if t.combined_score is not None else float("-inf")

    ordered = sorted(transcripts, key=lambda t: (-score(t), t.transcript_id))
    # Memoize each transcript's frozenset intron chain once, not per pair
    # (Phase 20 §2.5), the O(k²) sweep below would otherwise recompute it.
    chain_cache: dict[int, frozenset[tuple[int, int]]] = {}

    def chain_of(t: TranscriptCandidate) -> frozenset[tuple[int, int]]:
        c = chain_cache.get(id(t))
        if c is None:
            c = _intron_chain(t)
            chain_cache[id(t)] = c
        return c

    kept: list[TranscriptCandidate] = []
    dropped: list[tuple[str, str, float | None]] = []
    for t in ordered:
        twin = next(
            (
                k
                for k in kept
                if is_redundant(
                    t,
                    k,
                    min_cds_overlap,
                    min_cdna_overlap,
                    a_chain=chain_of(t),
                    b_chain=chain_of(k),
                )
            ),
            None,
        )
        if twin is None:
            kept.append(t)
        else:
            dropped.append((t.transcript_id, twin.transcript_id, t.combined_score))
    return kept, dropped


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------


def reconcile(
    helixer_loci: list[HelixerLocus],
    classifications: list[LocusClassification],
    mikado_loci: list[MikadoLocus],
    id_map: dict[str, str] | None = None,
    admit_novel: bool = False,
    novel_evidence_floor: float | None = None,
    reciprocal_overlap: float = 0.5,
    min_cds_overlap: float = 0.6,
    min_cdna_overlap: float = 0.6,
    allocator: IdAllocator | None = None,
    stats: Any = None,
    junctions: Any = None,
    merge_min_gap_reads: int = MERGE_MIN_GAP_READS,
    merge_require_canonical: bool = True,
    merge_tolerance: int = 0,
    genome: Any = None,
    paralog_identity_threshold: float | None = None,
    paralog_kmer_k: int = PARALOG_IDENTITY_K,
) -> tuple[list[ReconciledGene], dict[str, str], list[IsoformAdmission]]:
    """Reconcile Mikado loci onto the Helixer gene set.

    Returns ``(genes, id_map, admissions)``, the reconciled genes (sorted by
    seqid/start), the updated Helixer-id → HFG map (persist for stable IDs), and
    the isoform-admission audit.

    ``allocator`` (an :class:`IdAllocator`) decides where new HFG numbers come
    from; the default global policy (``base=1``, ``novel_base=NOVEL_ID_BASE``)
    reproduces the historical single-process numbering. A parallel chunk passes
    its reserved-range allocator so genome-wide IDs stay collision-free (§C3).

    ``stats`` (an optional :class:`~helixforge.reconcile.runstats.RunStats`) is
    bumped at each decision site, splits, merge accept/reject, redundant-isoform
    drops, when supplied. It only **observes**: increments never alter control
    flow, so the gene/tier/origin/AS counts and HFG ids are unchanged.

    Paralog/tandem-array merge guards: when a
    ``junctions`` set/index is supplied, an accepted merge additionally requires
    the bridging intron to be canonical and backed by >= ``merge_min_gap_reads``
    gap-spanning reads (see :func:`_has_bridging_introns`). When ``genome`` is
    given **and** ``paralog_identity_threshold`` is set, an adjacent pair with
    k-mer CDS identity at/above the threshold is treated as a recent duplication
    and its merge is rejected. Both guards default off (``junctions=None`` /
    ``paralog_identity_threshold=None``), preserving the legacy result; a rejected
    merge keeps each Helixer locus as a ``MERGE_REJECTED`` backstop (no gene lost).
    """
    alloc = allocator or IdAllocator()
    id_map = dict(id_map) if id_map else {}
    cls_by_id = {c.locus_id: c for c in classifications}

    def classification_for(locus: HelixerLocus) -> LocusClassification:
        return cls_by_id.get(
            locus.gene_id, LocusClassification(locus.gene_id, "SILENT")
        )

    corr = match_loci(helixer_loci, mikado_loci, reciprocal_overlap)
    genes = []
    admissions = []

    def _count_dropped(dropped: list[tuple[str, str, float | None]]) -> None:
        if stats is not None and dropped:
            stats.bump("isoforms_dropped_redundant", len(dropped))

    # --- 1:1 ---
    for h, m in corr.one_to_one:
        m2, dropped = _adopt(m, min_cds_overlap, min_cdna_overlap)
        _count_dropped(dropped)
        gene = build_reconciled_gene(
            h, m2, classification_for(h), "mikado_1to1", id_map, allocator=alloc
        )
        genes.append(gene)
        admissions.extend(record_admissions(gene, dropped))

    # --- splits ---
    for h, ms in corr.splits:
        if stats is not None:
            stats.bump("splits")
        base = assign_gene_id(h, id_map, allocator=alloc)
        for offset, m in enumerate(sorted(ms, key=lambda g: g.start)):
            child_id = f"{base}_{_SUFFIXES[offset]}"
            m2, dropped = _adopt(m, min_cds_overlap, min_cdna_overlap)
            _count_dropped(dropped)
            gene = build_reconciled_gene(
                h,
                m2,
                classification_for(h),
                "split",
                id_map,
                gene_id=child_id,
                extra_flags=[LOCUS_SPLIT],
                allocator=alloc,
            )
            genes.append(gene)
            admissions.extend(record_admissions(gene, dropped))

    # --- merges ---
    for hs, m in corr.merges:
        hs_sorted = sorted(hs, key=lambda g: g.start)
        paralogs = paralog_identity_threshold is not None and _high_identity_pair(
            hs_sorted, genome, paralog_identity_threshold, paralog_kmer_k
        )
        bridged = not paralogs and _has_bridging_introns(
            hs_sorted,
            m,
            junctions,
            min_gap_reads=merge_min_gap_reads,
            require_canonical=merge_require_canonical,
            tolerance=merge_tolerance,
        )
        if bridged:
            if stats is not None:
                stats.bump("merges_accepted")
            rep = hs_sorted[0]
            others = [g.gene_id for g in hs_sorted[1:]]
            m2, dropped = _adopt(m, min_cds_overlap, min_cdna_overlap)
            _count_dropped(dropped)
            gene = build_reconciled_gene(
                rep,
                m2,
                classification_for(rep),
                "merge",
                id_map,
                merged_from=others,
                extra_flags=[LOCUS_MERGE],
                allocator=alloc,
            )
            genes.append(gene)
            admissions.extend(record_admissions(gene, dropped))
        else:
            if stats is not None:
                stats.bump("merges_rejected")
            # reject: keep each Helixer locus separate as a flagged backstop.
            # A *prior* run may have accepted this merge and mapped every
            # constituent Helixer id to one shared HFG in ``id_map`` (an
            # accepted merge keeps the lowest id for all members).
            # Releasing the loci again must NOT hand that single id to more than
            # one gene: that emits duplicate gene/transcript IDs, violating the
            # stable-unique-ID guarantee (the scorer:stats IndexError crash; see
            # assessment-v4.md §2.3). The lowest-start (representative) locus
            # keeps the canonical id; any later locus whose id would collide is
            # re-allocated a fresh distinct HFG and re-keyed in ``id_map`` so the
            # new id is itself stable across reruns. A fresh run never collides
            # (each locus gets its own number), so this is a no-op there.
            seen_ids: set[str] = set()
            for h in hs_sorted:
                gid = assign_gene_id(h, id_map, allocator=alloc)
                if gid in seen_ids:
                    gid = _fmt_hfg(alloc.next_number(id_map))
                    id_map[h.gene_id] = gid
                seen_ids.add(gid)
                gene = build_reconciled_gene(
                    h,
                    None,
                    classification_for(h),
                    "helixer_backstop",
                    id_map,
                    gene_id=gid,
                    extra_flags=[MERGE_REJECTED],
                    allocator=alloc,
                )
                genes.append(gene)
                admissions.extend(record_admissions(gene))

    # --- helixer-only (completeness backstop) ---
    for h in corr.helixer_only:
        gene = build_reconciled_gene(
            h, None, classification_for(h), "helixer_backstop", id_map, allocator=alloc
        )
        genes.append(gene)
        admissions.extend(record_admissions(gene))

    # --- novel ---
    for m in corr.novel:
        if admit_novel and _meets_evidence_floor(m, novel_evidence_floor):
            m2, dropped = _adopt(m, min_cds_overlap, min_cdna_overlap)
            _count_dropped(dropped)
            novel_flags = [NOVEL_LOCUS]
            if m.locus_id in corr.tangled_loci:
                novel_flags.append(FROM_TANGLED_LOCUS)
            gene = build_reconciled_gene(
                None,
                m2,
                LocusClassification(m.locus_id, "LOW", evidence_source="none"),
                "novel",
                id_map,
                extra_flags=novel_flags,
                allocator=alloc,
            )
            genes.append(gene)
            admissions.extend(record_admissions(gene, dropped))
        else:
            for t in m.transcripts:
                admissions.append(
                    IsoformAdmission(
                        t.transcript_id,
                        m.locus_id,
                        False,
                        "novel_not_admitted",
                        score=t.combined_score,
                    )
                )

    genes.sort(key=lambda g: (g.seqid, g.start))
    return genes, id_map, admissions


def _adopt(
    mikado_locus: MikadoLocus,
    min_cds_overlap: float,
    min_cdna_overlap: float,
) -> tuple[MikadoLocus, list[tuple[str, str, float | None]]]:
    """Filter redundant isoforms; return (mikado_locus_with_kept, dropped)."""
    kept, dropped = _filter_redundant(
        mikado_locus.transcripts, min_cds_overlap, min_cdna_overlap
    )
    return attrs.evolve(mikado_locus, transcripts=kept), dropped


def _meets_evidence_floor(mikado_locus: MikadoLocus, floor: float | None) -> bool:
    if floor is None:
        return True
    best = max(
        (
            t.combined_score
            for t in mikado_locus.transcripts
            if t.combined_score is not None
        ),
        default=None,
    )
    return best is not None and best >= floor
