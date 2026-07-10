"""TRaCE — Transcript Ranking and Canonical Election."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from helixforge.reconcile.as_events import overlap_bases

if TYPE_CHECKING:
    from helixforge.reconcile.models import StringTieTranscript, TranscriptCandidate


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


@dataclass
class TraceParams:
    """TRaCE cutoffs + vote weights (defaults are the paper's).

    Sample-voter cutoffs (``max_aed`` / ``min_tpm`` / ``min_overlap``) gate which
    candidates a sample ballot may rank. The length-voter weights prioritise
    ``domain`` > ``protein`` > ``cdna`` (the TRaCE figure uses 9 / 6 / 3).
    ``use_domain`` enables the optional domain voter (only when per-isoform domain
    coverage is supplied). ``balance_samples`` scales the *collective* sample-vote
    weight to ``balance_target`` (default = total length-voter weight) so a large
    sample count does not swamp the length voters; with it off each sample carries
    unit weight.
    """

    max_aed: float = 0.5
    min_tpm: float = 0.5
    min_overlap: float = 0.5
    weight_domain: float = 9.0
    weight_protein: float = 6.0
    weight_cdna: float = 3.0
    use_domain: bool = True
    balance_samples: bool = True
    balance_target: float | None = None


@dataclass
class Ballot:
    """One voter's ranking of the candidates.

    ``ranks`` maps a candidate ``transcript_id`` to its 1-based rank (1 = best);
    a candidate the voter does not rank maps to ``None``. Ties share a rank.
    ``kind`` is ``"sample"`` or ``"length"``; ``voter`` is the sample id (sample
    ballots) or the metric name ``"domain"`` / ``"protein"`` / ``"cdna"`` (length
    ballots) — used for weighting + electorate balancing.
    """

    ranks: dict[str, int | None]
    kind: str
    voter: str


# ---------------------------------------------------------------------------
# Geometry helpers (reuse overlap_bases; never reimplement exon overlap)
# ---------------------------------------------------------------------------


def _exon_tuples(t: Any) -> list[tuple[int, int]]:
    return [(e.start, e.end) for e in t.exons]


def _total(intervals: Sequence[tuple[int, int]]) -> int:
    return sum(e - s for s, e in intervals)


def _clip(
    intervals: Sequence[tuple[int, int]], lo: int, hi: int
) -> list[tuple[int, int]]:
    """Clip ``intervals`` to the genomic window ``[lo, hi)``."""
    out: list[tuple[int, int]] = []
    for s, e in intervals:
        ns, ne = max(s, lo), min(e, hi)
        if ne > ns:
            out.append((ns, ne))
    return out


def proportion_overlap(
    a_exons: Sequence[tuple[int, int]], b_exons: Sequence[tuple[int, int]]
) -> float:
    """Exon-overlap fraction of the *smaller* structure (``ov / min(la, lb)``)."""
    la, lb = _total(a_exons), _total(b_exons)
    if la == 0 or lb == 0:
        return 0.0
    return overlap_bases(a_exons, b_exons) / min(la, lb)


def structural_aed(
    candidate_exons: Sequence[tuple[int, int]],
    assembled_exons: Sequence[tuple[int, int]],
) -> float:
    """Annotation edit distance between two exon structures, on the overlap region.

    AED = ``1 - (SN + SP) / 2`` where, restricted to the genomic span both share,
    ``SN`` = fraction of the assembled (evidence) bases recovered and ``SP`` =
    fraction of the candidate bases matching the evidence. Restricting to the
    overlap region means a partially-assembled low-expression transcript is not
    penalised for stopping short (TRaCE convention). Returns ``1.0`` (maximally
    distant) when the structures do not overlap.
    """
    ca0, aa0 = list(candidate_exons), list(assembled_exons)
    if not ca0 or not aa0:
        return 1.0
    lo = max(min(s for s, _ in ca0), min(s for s, _ in aa0))
    hi = min(max(e for _, e in ca0), max(e for _, e in aa0))
    if hi <= lo:
        return 1.0
    ca = _clip(ca0, lo, hi)
    aa = _clip(aa0, lo, hi)
    la, lb = _total(ca), _total(aa)
    if la == 0 or lb == 0:
        return 1.0
    ov = overlap_bases(ca, aa)
    sn = ov / lb
    sp = ov / la
    return round(1.0 - (sn + sp) / 2.0, 6)


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------


def _ranks_from_scores(
    scores: dict[str, float],
    *,
    ascending: bool,
    all_ids: Sequence[str],
) -> dict[str, int | None]:
    """Competition-rank ``scores`` (1 = best); ties share a rank, unscored → None.

    ``ascending`` ranks the smallest value best (AED); otherwise the largest value
    best (length / coverage). Ids in ``all_ids`` absent from ``scores`` map to
    ``None`` (unranked).
    """
    ordered = sorted(scores.items(), key=lambda kv: kv[1] if ascending else -kv[1])
    out: dict[str, int | None] = {}
    prev_val: float | None = None
    prev_rank = 0
    for i, (cid, val) in enumerate(ordered):
        if prev_val is not None and val == prev_val:
            out[cid] = prev_rank
        else:
            out[cid] = i + 1
            prev_rank = i + 1
            prev_val = val
    for cid in all_ids:
        out.setdefault(cid, None)
    return out


def aed_ballot(
    candidates: list[TranscriptCandidate],
    sample_transcripts: list[StringTieTranscript],
    *,
    max_aed: float = 0.5,
    min_tpm: float = 0.5,
    min_overlap: float = 0.5,
) -> dict[str, int | None]:
    """One sample's ballot: rank each candidate by AED to its best assembled tx.

    For each candidate, the sample's *most highly expressed* assembled transcript
    that overlaps it (same seqid+strand, ``tpm >= min_tpm``, proportion-overlap
    ``>= min_overlap``) is chosen; the candidate's AED to it (overlap region only)
    is its score. A candidate whose best-expressed overlapping transcript has AED
    ``> max_aed`` — or has no qualifying transcript — is left unranked (``None``).
    Lower AED ⇒ better rank; ties share a rank.
    """
    usable = [s for s in sample_transcripts if (s.tpm or 0.0) >= min_tpm]
    aeds: dict[str, float] = {}
    for c in candidates:
        c_exons = _exon_tuples(c)
        best_tpm: float | None = None
        best_aed = 1.0
        for s in usable:
            if s.seqid != c.seqid or s.strand != c.strand:
                continue
            s_exons = _exon_tuples(s)
            if proportion_overlap(c_exons, s_exons) < min_overlap:
                continue
            if best_tpm is None or s.tpm > best_tpm:
                best_tpm = s.tpm
                best_aed = structural_aed(c_exons, s_exons)
        if best_tpm is not None and best_aed <= max_aed:
            aeds[c.transcript_id] = best_aed
    return _ranks_from_scores(
        aeds, ascending=True, all_ids=[c.transcript_id for c in candidates]
    )


def length_ballot(
    candidates: list[TranscriptCandidate],
    key: str,
    values: dict[str, float] | None = None,
) -> dict[str, int | None]:
    """Rank candidates by a length-style metric (higher = better rank).

    ``key`` is ``"protein_length"`` (total CDS length), ``"cdna_length"`` (total
    exon length), or ``"domain_coverage"`` (read from ``values``; missing → 0.0).
    Every candidate is ranked (no ``None``); ties share a rank.
    """
    if key == "protein_length":
        scores = {c.transcript_id: float(c.total_cds_length) for c in candidates}
    elif key == "cdna_length":
        scores = {c.transcript_id: float(c.total_exon_length) for c in candidates}
    elif key == "domain_coverage":
        vals = values or {}
        scores = {
            c.transcript_id: float(vals.get(c.transcript_id, 0.0)) for c in candidates
        }
    else:
        raise ValueError(
            "key must be 'protein_length', 'cdna_length', or 'domain_coverage', "
            f"got {key!r}"
        )
    return _ranks_from_scores(
        scores, ascending=False, all_ids=[c.transcript_id for c in candidates]
    )


# ---------------------------------------------------------------------------
# The RCV election
# ---------------------------------------------------------------------------


def _ballot_weights(
    ballots: list[Ballot],
    weights: dict[str, float],
    *,
    balance_samples: bool,
    balance_target: float | None,
) -> list[tuple[Ballot, float]]:
    """Resolve each ballot's weight (length: per-metric; sample: balanced share)."""
    length_total = sum(
        float(weights.get(b.voter, 0.0)) for b in ballots if b.kind == "length"
    )
    sample_ballots = [b for b in ballots if b.kind == "sample"]
    n_samples = len(sample_ballots)
    target = balance_target if balance_target is not None else length_total

    out: list[tuple[Ballot, float]] = []
    for b in ballots:
        if b.kind == "length":
            out.append((b, float(weights.get(b.voter, 0.0))))
        elif balance_samples and n_samples > 0 and target > 0:
            out.append((b, target / n_samples))
        else:
            out.append((b, 1.0))
    return out


def _elect_one(
    remaining: list[str],
    weighted: list[tuple[Ballot, float]],
    tiebreak: dict[str, tuple[float, str]],
) -> str:
    """Elect one seat among ``remaining`` by lexicographic per-level vote vectors.

    Each ballot contributes its weight to the candidate(s) at its best *remaining*
    rank (level 0), then its next rank (level 1), and so on. Candidates are
    compared by their vote vector descending (rank-1 votes, then rank-2, …); a
    full tie is broken by ``tiebreak`` (``combined_score`` desc, ``transcript_id``
    asc).
    """
    vectors: dict[str, list[float]] = {c: [] for c in remaining}
    for b, w in weighted:
        ranked = [(c, b.ranks.get(c)) for c in remaining if b.ranks.get(c) is not None]
        if not ranked:
            continue
        ranked.sort(key=lambda cv: cv[1])  # type: ignore[arg-type,return-value]
        level = -1
        prev: int | None = None
        for c, rv in ranked:
            if rv != prev:
                level += 1
                prev = rv
            vec = vectors[c]
            while len(vec) <= level:
                vec.append(0.0)
            vec[level] += w

    maxlen = max((len(v) for v in vectors.values()), default=0)

    def key(c: str) -> tuple[Any, ...]:
        vec = vectors[c]
        padded = tuple(-(vec[i] if i < len(vec) else 0.0) for i in range(maxlen))
        return padded + tiebreak[c]

    return min(remaining, key=key)


def elect(
    candidates: list[TranscriptCandidate],
    ballots: list[Ballot],
    weights: dict[str, float],
    *,
    balance_samples: bool = True,
    balance_target: float | None = None,
) -> list[str]:
    """Run the multi-round RCV election; return the full ordering of candidate ids.

    Elects one seat at a time (winner first), removing the winner and re-tallying
    among the remaining until all candidates are ordered. Deterministic: the only
    randomness-free tie-break is ``(combined_score desc, transcript_id asc)``.
    """
    weighted = _ballot_weights(
        ballots,
        weights,
        balance_samples=balance_samples,
        balance_target=balance_target,
    )
    tiebreak = {
        c.transcript_id: (
            -(c.combined_score if c.combined_score is not None else float("-inf")),
            c.transcript_id,
        )
        for c in candidates
    }
    remaining = [c.transcript_id for c in candidates]
    order: list[str] = []
    while remaining:
        winner = _elect_one(remaining, weighted, tiebreak)
        order.append(winner)
        remaining.remove(winner)
    return order


# ---------------------------------------------------------------------------
# Domain coverage (optional voter) + top-level ordering
# ---------------------------------------------------------------------------


def domain_coverage_by_transcript(
    genes: list[Any],
    functional_records: dict[str, Any] | None,
) -> dict[str, float]:
    """Per-candidate domain-coverage proxy from functional records (optional path).

    The domain voter needs coverage for **every** candidate, but functional
    annotation currently attaches to the **primary** transcript only
    (``prep/function.py``), and the parsed :class:`~prep.function.FunctionalRecord`
    keeps a boolean ``domain_complete`` rather than a covered-fraction. So this is
    a deliberately coarse proxy: ``1.0`` when a transcript has a record flagged
    ``domain_complete``, else ``0.0``; transcripts with no record are simply
    absent (the voter treats them as ``0.0``). A true per-isoform coverage
    fraction requires running functional annotation on **all** isoforms with
    coordinate-aware domain spans — documented future work. With no records (the
    default pipeline path) this returns ``{}`` and the domain voter is disabled.
    """
    out: dict[str, float] = {}
    if not functional_records:
        return out
    for g in genes:
        for t in g.transcripts:
            rec = functional_records.get(t.transcript_id)
            if rec is None:
                continue
            out[t.transcript_id] = (
                1.0 if getattr(rec, "domain_complete", False) else 0.0
            )
    return out


def trace_order(
    gene_transcripts: list[TranscriptCandidate],
    sample_transcripts_by_sample: dict[str, list[StringTieTranscript]],
    domain_coverage: dict[str, float] | None = None,
    params: TraceParams | None = None,
) -> list[TranscriptCandidate]:
    """Reorder one gene's transcripts by the TRaCE election (winner first).

    Builds the sample + length (+ optional domain) ballots, holds the RCV
    election, and returns ``gene_transcripts`` permuted into the elected order.
    A single-transcript gene, or a gene with **no** expression/domain evidence,
    is returned in its input order unchanged (the caller then keeps the
    ``combined_score`` primary) — length voters alone are intentionally *not*
    allowed to reorder, since that would promote intron-retention-inflated
    longest isoforms (the very failure TRaCE exists to avoid).
    """
    params = params or TraceParams()
    cands = list(gene_transcripts)
    if len(cands) < 2:
        return cands

    sample_ballots: list[Ballot] = []
    for sid in sorted(sample_transcripts_by_sample):
        ranks = aed_ballot(
            cands,
            sample_transcripts_by_sample[sid],
            max_aed=params.max_aed,
            min_tpm=params.min_tpm,
            min_overlap=params.min_overlap,
        )
        if any(v is not None for v in ranks.values()):
            sample_ballots.append(Ballot(ranks=ranks, kind="sample", voter=sid))

    has_domain = params.use_domain and bool(domain_coverage)
    if not sample_ballots and not has_domain:
        return cands  # no-evidence fallback (caller keeps combined_score primary)

    length_ballots: list[Ballot] = []
    if has_domain:
        length_ballots.append(
            Ballot(
                length_ballot(cands, "domain_coverage", domain_coverage),
                "length",
                "domain",
            )
        )
    length_ballots.append(
        Ballot(length_ballot(cands, "protein_length"), "length", "protein")
    )
    length_ballots.append(Ballot(length_ballot(cands, "cdna_length"), "length", "cdna"))

    weights = {
        "domain": params.weight_domain,
        "protein": params.weight_protein,
        "cdna": params.weight_cdna,
    }
    order_ids = elect(
        cands,
        sample_ballots + length_ballots,
        weights,
        balance_samples=params.balance_samples,
        balance_target=params.balance_target,
    )
    by_id = {c.transcript_id: c for c in cands}
    return [by_id[i] for i in order_ids]
