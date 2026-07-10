"""Junction correction of backstop-gene introns."""

from __future__ import annotations

from typing import Any

import attrs

from helixforge.constants import MIN_EXON_BP, MIN_INTRON_BP
from helixforge.qc.flags import (
    ALL_JUNCTIONS_SUPPORTED,
    NON_CANONICAL_SPLICE,
    NO_JUNCTION_SUPPORT,
    PARTIAL_JUNCTION_SUPPORT,
    dedup_flags,
)
from helixforge.reconcile.models import (
    CANONICAL_SPLICE_MOTIFS,
    Exon,
    Interval,
    ReconciledGene,
    SpliceJunction,
    TranscriptCandidate,
)
from helixforge.utils.intervals import IntervalIndex
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)

# Structural floors verified by the transaction; re-exported from
# helixforge.constants. MIN_INTRON_BP is biology-bearing.

SUPPORTED = "SUPPORTED"
CONTRADICTED = "CONTRADICTED"
NO_EVIDENCE = "NO_EVIDENCE"


# ---------------------------------------------------------------------------
# Junction index
# ---------------------------------------------------------------------------
#
# The list-based matchers below scan the *entire* junction set per intron,
# O(introns × all_junctions), called once per backstop gene. At maize scale
# (J ~ 10⁶) that is the second O(N²)-flavored cliff after the ID allocator.
# JunctionIndex builds the set once into, per (seqid, strand):
#   * an exact ``(donor, acceptor) -> junction`` dict for the tolerance==0 fast
#     path (O(1) lookup);
#   * point IntervalIndexes (NCLS) on donor and on acceptor positions for the
#     tolerance/contradiction queries (O(log J) candidate retrieval).
# Match/contradiction decisions are **identical** to the linear scan: the dict
# stores the highest-read-count junction per (donor, acceptor) with first-on-tie
# (mirroring the scan's strict ``>``), and candidate ids come back in insertion
# (== bucket == original) order so ties resolve to the same junction.


class _Bucket:
    """Per-``(seqid, strand)`` lookup structures."""

    __slots__ = ("juncs", "exact", "donor_idx", "acceptor_idx")

    def __init__(self) -> None:
        self.juncs: list[SpliceJunction] = []
        self.exact: dict[tuple[int, int], SpliceJunction] = {}
        self.donor_idx: IntervalIndex = IntervalIndex()
        self.acceptor_idx: IntervalIndex = IntervalIndex()


class JunctionIndex:
    """Indexed junction set for O(1)/O(log J) match & contradiction queries.

    Build once with :meth:`from_junctions` (or the constructor), then thread it
    wherever a junction *list* would go: :func:`find_matching_junction`,
    :func:`find_contradicting_junction`, :func:`classify_introns` and
    :func:`refine_backstop_gene` all accept either a list (linear scan, unchanged)
    or a ``JunctionIndex`` (indexed path, identical results).
    """

    def __init__(self, junctions: Any = ()) -> None:
        self._buckets: dict[tuple[str, str], _Bucket] = {}
        donor_ivals: dict[tuple[str, str], list[tuple[int, int, int]]] = {}
        acc_ivals: dict[tuple[str, str], list[tuple[int, int, int]]] = {}
        for j in junctions:
            key = (j.seqid, j.strand)
            bucket = self._buckets.get(key)
            if bucket is None:
                bucket = self._buckets[key] = _Bucket()
                donor_ivals[key] = []
                acc_ivals[key] = []
            i = len(bucket.juncs)
            bucket.juncs.append(j)
            # First-on-tie, max-read-count otherwise, mirrors the scan's strict `>`.
            ek = (j.donor, j.acceptor)
            cur = bucket.exact.get(ek)
            if cur is None or j.read_count > cur.read_count:
                bucket.exact[ek] = j
            donor_ivals[key].append((j.donor, j.donor + 1, i))
            acc_ivals[key].append((j.acceptor, j.acceptor + 1, i))
        for key, bucket in self._buckets.items():
            bucket.donor_idx.add_intervals(donor_ivals[key])
            bucket.acceptor_idx.add_intervals(acc_ivals[key])

    @classmethod
    def from_junctions(cls, junctions: Any) -> JunctionIndex:
        """Build a :class:`JunctionIndex` from an iterable of ``SpliceJunction``."""
        return cls(junctions)

    def __len__(self) -> int:
        return sum(len(b.juncs) for b in self._buckets.values())

    def find_matching(
        self,
        intron: Interval,
        seqid: str,
        strand: str,
        tolerance: int = 0,
        min_reads: int = 3,
    ) -> SpliceJunction | None:
        """Indexed equivalent of :func:`find_matching_junction`."""
        bucket = self._buckets.get((seqid, strand))
        if bucket is None:
            return None
        if tolerance == 0:
            j = bucket.exact.get((intron.start, intron.end))
            if j is not None and j.read_count >= min_reads:
                return j
            return None
        cand = bucket.donor_idx.query(
            intron.start - tolerance, intron.start + tolerance + 1
        )
        best = None
        for i in cand:
            j = bucket.juncs[i]
            if j.read_count < min_reads:
                continue
            if (
                abs(j.donor - intron.start) <= tolerance
                and abs(j.acceptor - intron.end) <= tolerance
            ):
                if best is None or j.read_count > best.read_count:
                    best = j
        return best

    def find_contradicting(
        self,
        intron: Interval,
        seqid: str,
        strand: str,
        tolerance: int = 0,
        min_reads: int = 3,
    ) -> SpliceJunction | None:
        """Indexed equivalent of :func:`find_contradicting_junction`."""
        bucket = self._buckets.get((seqid, strand))
        if bucket is None:
            return None
        donor = bucket.donor_idx.query(
            intron.start - tolerance, intron.start + tolerance + 1
        )
        acc = bucket.acceptor_idx.query(
            intron.end - tolerance, intron.end + tolerance + 1
        )
        best = None
        for i in sorted(set(donor) | set(acc)):
            j = bucket.juncs[i]
            if j.read_count < min_reads:
                continue
            donor_match = abs(j.donor - intron.start) <= tolerance
            acceptor_match = abs(j.acceptor - intron.end) <= tolerance
            if donor_match == acceptor_match:
                continue
            if best is None or j.read_count > best.read_count:
                best = j
        return best


# ---------------------------------------------------------------------------
# Junction matching
# ---------------------------------------------------------------------------


def find_matching_junction(
    intron: Interval,
    junctions: list[SpliceJunction] | JunctionIndex,
    seqid: str,
    strand: str,
    tolerance: int = 0,
    min_reads: int = 3,
) -> SpliceJunction | None:
    """Best same-strand junction matching **both** intron boundaries (≥ ``min_reads``).

    A match has ``|donor - intron.start| <= tolerance`` **and**
    ``|acceptor - intron.end| <= tolerance``. Returns the highest-read-count
    match, or ``None``. ``junctions`` may be a list (linear scan) or a
    :class:`JunctionIndex` (O(1)/O(log J), identical result).
    """
    if isinstance(junctions, JunctionIndex):
        return junctions.find_matching(intron, seqid, strand, tolerance, min_reads)
    best = None
    for j in junctions:
        if j.seqid != seqid or j.strand != strand or j.read_count < min_reads:
            continue
        if (
            abs(j.donor - intron.start) <= tolerance
            and abs(j.acceptor - intron.end) <= tolerance
        ):
            if best is None or j.read_count > best.read_count:
                best = j
    return best


def find_contradicting_junction(
    intron: Interval,
    junctions: list[SpliceJunction] | JunctionIndex,
    seqid: str,
    strand: str,
    tolerance: int = 0,
    min_reads: int = 3,
) -> SpliceJunction | None:
    """Best same-strand junction that shares **exactly one** intron boundary.

    Represents a well-supported alternative splice site: one end agrees (within
    ``tolerance``), the other disagrees (beyond ``tolerance``). The disagreeing
    end is what gets corrected. Returns the highest-read-count such junction, or
    ``None``. Junctions sharing *neither* boundary are intentionally ignored,
    relocating both ends to an unrelated junction is too aggressive (it stays
    ``NO_EVIDENCE`` and is flagged, not moved). ``junctions`` may be a list or a
    :class:`JunctionIndex` (identical result).
    """
    if isinstance(junctions, JunctionIndex):
        return junctions.find_contradicting(intron, seqid, strand, tolerance, min_reads)
    best = None
    for j in junctions:
        if j.seqid != seqid or j.strand != strand or j.read_count < min_reads:
            continue
        donor_match = abs(j.donor - intron.start) <= tolerance
        acceptor_match = abs(j.acceptor - intron.end) <= tolerance
        if donor_match == acceptor_match:
            continue  # both match (that's a match) or neither (ignored)
        if best is None or j.read_count > best.read_count:
            best = j
    return best


def classify_introns(
    transcript: TranscriptCandidate,
    junctions: list[SpliceJunction] | JunctionIndex,
    tolerance: int = 0,
    min_reads: int = 3,
) -> list[tuple[str, SpliceJunction | None]]:
    """Classify each intron of ``transcript`` against the junction set.

    Returns a list (genomic-ascending intron order) of ``(status, junction)``:
    ``SUPPORTED`` (matched), ``CONTRADICTED`` (one boundary off, junction is the
    correction target), or ``NO_EVIDENCE`` (junction ``None``).
    """
    seqid, strand = transcript.seqid, transcript.strand
    out: list[tuple[str, SpliceJunction | None]] = []
    for intron in transcript.introns:
        match = find_matching_junction(
            intron, junctions, seqid, strand, tolerance, min_reads
        )
        if match is not None:
            out.append((SUPPORTED, match))
            continue
        contra = find_contradicting_junction(
            intron, junctions, seqid, strand, tolerance, min_reads
        )
        if contra is not None:
            out.append((CONTRADICTED, contra))
        else:
            out.append((NO_EVIDENCE, None))
    return out


# ---------------------------------------------------------------------------
# Transaction: apply corrections, verify once, revert on any failure
# ---------------------------------------------------------------------------


def apply_intron_corrections(
    transcript: TranscriptCandidate,
    corrections: list[tuple[int, int, int]],
) -> TranscriptCandidate | None:
    """Apply intron-boundary corrections atomically; return a new transcript or None.

    ``corrections`` is a list of ``(intron_index, new_donor, new_acceptor)``.
    Intron ``i`` sits between exon ``i`` and exon ``i+1``: ``new_donor`` becomes
    exon ``i``'s end, ``new_acceptor`` becomes exon ``i+1``'s start. **All**
    corrections are applied to a working copy, then the structure is verified
    **once**:

    * every exon ``start < end`` (no 0-width / inverted exons),
    * exons sorted ascending and non-overlapping,
    * the span ``[transcript.start, transcript.end)`` is preserved,
    * every exon ``>= MIN_EXON_BP``,
    * every intron ``>= MIN_INTRON_BP``.

    Any failure, including a ``ValueError`` from model construction, reverts the
    whole transaction by returning ``None``. The original
    transcript object is never modified.
    """
    if not corrections:
        return transcript

    n = len(transcript.exons)
    bounds: list[tuple[int, int]] = [(ex.start, ex.end) for ex in transcript.exons]

    for idx, new_donor, new_acceptor in corrections:
        if idx < 0 or idx + 1 >= n:
            return None  # no such intron, revert
        bounds[idx] = (bounds[idx][0], new_donor)
        bounds[idx + 1] = (new_acceptor, bounds[idx + 1][1])

    if not _verify_bounds(bounds, transcript.start, transcript.end):
        return None

    try:
        new_exons = [Exon(s, e) for (s, e) in bounds]
        return attrs.evolve(transcript, exons=new_exons)
    except ValueError as exc:
        # Construction-time invariant (e.g. CDS no longer within exons) failed.
        _log.info(
            "intron correction for %s reverted: %s", transcript.transcript_id, exc
        )
        return None


def _verify_bounds(
    bounds: list[tuple[int, int]], span_start: int, span_end: int
) -> bool:
    """Verify the full set of structural invariants for ``bounds`` (once)."""
    prev_end = None
    for s, e in bounds:
        if e <= s:
            return False  # 0-width or inverted exon
        if (e - s) < MIN_EXON_BP:
            return False
        if prev_end is not None:
            if s < prev_end:
                return False  # overlap / out of order
            if (s - prev_end) < MIN_INTRON_BP:
                return False  # intron too short
        prev_end = e
    if bounds[0][0] < span_start or bounds[-1][1] > span_end:
        return False  # exon outside the gene span
    return True


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def _cds_introns_coherent(cds_segments: Any, exons: Any) -> bool:
    """True iff every CDS-derived intron coincides with an exon intron.

    CDS segments are stored one-per-exon (split at introns). The intron implied
    between two consecutive CDS segments, ``(cds[i].end, cds[i+1].start)`` in the
    0-based half-open convention, must equal a real exon intron
    ``(exon[j].end, exon[j+1].start)``. When it does not, an *internal* CDS
    boundary (one facing an intron) sits inside an exon instead of on the splice
    site, the CDS is not a valid spliced ORF: Mikado's transcript finalizer
    derives a CDS intron that matches no exon intron and asserts
    ``len(cds_introns) > 0``. Terminal CDS boundaries (5'/3' UTR interior to the
    first/last coding exon) are unconstrained and never checked.
    """
    if cds_segments is None or len(cds_segments) < 2:
        return True
    ex = sorted((e.start, e.end) for e in exons)
    exon_introns = {(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)}
    seg = sorted((c.start, c.end) for c in cds_segments)
    return all(
        (seg[i][1], seg[i + 1][0]) in exon_introns for i in range(len(seg) - 1)
    )


def _reattach_cds(
    corrected: TranscriptCandidate,
    saved_cds: list[Any],
    saved_partial: bool,
) -> TranscriptCandidate:
    """Re-attach a junction-corrected transcript's saved CDS to its new exons.

    The structure is corrected on a CDS-stripped copy (so the model's own ORF
    never blocks a valid junction fix); the CDS is then re-attached here. Fast
    path: the saved CDS still sits within the corrected exons **and** is coherent
    with the corrected introns (every internal CDS boundary lands on a splice
    site), so it is restored verbatim. Otherwise it is re-clipped to the corrected
    exons (complete-ORF phase recompute) and coherence-checked again.

    A correction that *grows* an exon outward past a flush CDS boundary leaves the
    CDS contained but no longer flush, its CDS-derived intron no longer matches
    any exon intron. Re-clipping cannot fix that (clipping only removes intronic
    overhang, it never extends a segment to a moved boundary), so such a CDS is
    rejected as a structurally-invalid spliced ORF and the transcript is left
    CDS-less. A downstream backstop source (miniprot) may still rescue it; else
    ``validate`` flags the CDS-less gene. A flagged gene beats an incoherent one
    (CLAUDE.md §6, reject the source rather than emit CDS Mikado cannot finalize).
    """
    try:
        verbatim = attrs.evolve(
            corrected, cds=list(saved_cds), cds_partial=saved_partial
        )
    except ValueError:
        verbatim = None
    if verbatim is not None and _cds_introns_coherent(verbatim.cds, verbatim.exons):
        return verbatim

    from helixforge.reconcile.cds import project_cds_to_exons

    reproj = project_cds_to_exons(saved_cds, corrected.exons, corrected.strand)
    if reproj is not None and _cds_introns_coherent(reproj, corrected.exons):
        try:
            return attrs.evolve(corrected, cds=reproj)
        except ValueError:
            pass
    return corrected


def refine_backstop_gene(
    gene: ReconciledGene,
    junctions: list[SpliceJunction] | JunctionIndex,
    tolerance: int = 0,
    min_reads: int = 3,
    stats: Any = None,
) -> ReconciledGene:
    """Junction-correct a backstop gene's introns; return the (possibly) new gene.

    Acts only when ``origin == 'helixer_backstop'`` and the gene's transcript has
    introns. Matched introns are kept; contradicted introns are corrected via the
    transaction (reverting all on any invariant failure); no-evidence introns are
    kept and the gene is flagged. Mikado-origin and single-exon genes are
    returned unchanged. Sets a splice-support flag and the transcript's
    ``junction_support_fraction``.

    ``stats`` (optional :class:`~helixforge.reconcile.runstats.RunStats`) records
    how many individual intron corrections were applied vs. reverted (per the
    all-or-nothing transaction); observation only, it never changes the result.
    """
    if gene.origin != "helixer_backstop":
        return gene

    transcript = gene.transcripts[0]
    n_introns = transcript.num_introns
    if n_introns == 0:
        return gene  # single-exon: never insert introns

    # Junction-correct on a CDS-stripped copy: the model's own (Helixer-intrinsic)
    # CDS must not block a valid intron correction on CDS-containment (correct the
    # structure first, then re-project the CDS, the documented order). The CDS is
    # re-attached to the corrected exons at the end.
    saved_cds = transcript.cds
    saved_partial = transcript.cds_partial
    work_tx = attrs.evolve(transcript, cds=None) if saved_cds is not None else transcript

    classified = classify_introns(work_tx, junctions, tolerance, min_reads)
    corrections: list[tuple[int, int, int]] = []
    non_canonical_rejected = False
    for idx, (status, junction) in enumerate(classified):
        if status == CONTRADICTED:
            assert (
                junction is not None
            )  # classify_introns always pairs CONTRADICTED with a junction
            # Canonicity hard gate: this is the one place coordinates are
            # mutated, so a correction that
            # would relocate a splice site to a **non-canonical** motif is
            # rejected outright: the intron stays uncorrected and is flagged. A
            # junction whose motif was never evaluated (``canonical is None``,
            # e.g. no STAR motif / no genome) is accepted as before, keeping the
            # legacy path byte-identical.
            if _is_non_canonical(junction):
                non_canonical_rejected = True
                if stats is not None:
                    stats.bump("junction_corrections_rejected_noncanonical", 1)
                _log.info(
                    "backstop %s: rejected non-canonical junction correction "
                    "(%d-%d, motif=%s)",
                    gene.gene_id,
                    junction.donor,
                    junction.acceptor,
                    junction.canonical,
                )
                continue
            corrections.append((idx, junction.donor, junction.acceptor))

    new_transcript = work_tx
    corrected_ok = False
    if corrections:
        candidate = apply_intron_corrections(work_tx, corrections)
        if candidate is not None:
            new_transcript = candidate
            corrected_ok = True
            if stats is not None:
                stats.bump("junction_corrections_applied", len(corrections))
        else:
            if stats is not None:
                stats.bump("junction_corrections_reverted", len(corrections))
            _log.info(
                "backstop %s: intron correction reverted; keeping original structure",
                gene.gene_id,
            )

    # Final support: matched always counts; a corrected intron counts only if it
    # was both accepted by the canonicity gate (in ``corrections``) and the whole
    # transaction applied. A non-canonical-rejected intron is therefore *not*
    # supported (it was neither matched nor corrected).
    corrected_indices = {idx for idx, _, _ in corrections} if corrected_ok else set()
    supported = sum(
        1
        for i, (status, _) in enumerate(classified)
        if status == SUPPORTED or (status == CONTRADICTED and i in corrected_indices)
    )
    if supported == n_introns:
        splice_flag = ALL_JUNCTIONS_SUPPORTED
    elif supported == 0:
        splice_flag = NO_JUNCTION_SUPPORT
    else:
        splice_flag = PARTIAL_JUNCTION_SUPPORT

    # INFO flag if any intron touches a non-canonical motif (rejected correction
    # target, or a matched junction that is itself non-canonical).
    non_canonical_seen = non_canonical_rejected or any(
        status in (SUPPORTED, CONTRADICTED) and _is_non_canonical(junction)
        for status, junction in classified
    )
    extra_flags = [splice_flag]
    if non_canonical_seen:
        extra_flags.append(NON_CANONICAL_SPLICE)

    # Re-attach the model's CDS to the (possibly corrected) exons.
    if saved_cds is not None:
        new_transcript = _reattach_cds(new_transcript, saved_cds, saved_partial)

    new_transcript = attrs.evolve(
        new_transcript, junction_support_fraction=supported / n_introns
    )
    new_gene = attrs.evolve(
        gene,
        transcripts=[new_transcript],
        flags=dedup_flags([*gene.flags, *extra_flags]),
    )
    return new_gene


def _is_non_canonical(junction: SpliceJunction | None) -> bool:
    """True iff ``junction`` has an evaluated, non-canonical motif.

    ``None`` junction or ``canonical is None`` (never evaluated) -> False, so the
    canonicity gate is inert wherever motif info is unavailable.
    """
    return (
        junction is not None
        and junction.canonical is not None
        and junction.canonical not in CANONICAL_SPLICE_MOTIFS
    )
