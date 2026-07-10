"""Boundary-stitch pass: recover cross-chunk merges."""

from __future__ import annotations

from collections import namedtuple
from typing import TYPE_CHECKING, Any

import attrs

from helixforge.parallel.plan import _region_bounds
from helixforge.qc.flags import LOCUS_MERGE, dedup_flags
from helixforge.reconcile.mikado_integrate import (
    _has_bridging_introns,
    _hfg_num,
    _renumber,
)
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.parallel.plan import Plan
    from helixforge.reconcile.models import ReconciledGene

_log = get_logger(__name__)

# A lightweight gene span for callers that only have coordinates (e.g. the driver
# reloading an aggregated GFF3); real ``ReconciledGene``s expose the same four
# attributes, so both flow through :func:`find_boundary_merge_candidates`.
GeneSpan = namedtuple("GeneSpan", "gene_id seqid start end strand")

_PseudoIntron = namedtuple("_PseudoIntron", "start end")


class _PseudoLocus:
    """Minimal stand-in so candidate junctions feed ``_has_bridging_introns``."""

    def __init__(self, introns: list[tuple[int, int]]) -> None:
        self.transcripts = [
            namedtuple("_T", "introns")([_PseudoIntron(s, e) for (s, e) in introns])
        ]


# ---------------------------------------------------------------------------
# Boundary geometry
# ---------------------------------------------------------------------------


def _plan_boundaries(plan: Plan) -> dict[str, list[int]]:
    """``{seqid: [interior cut positions]}`` from a plan's chunk regions.

    An interior boundary on a scaffold is a coordinate that is both the ``hi`` of
    one chunk region and the ``lo`` of the next (the shared cut point). Whole-
    scaffold (bare-seqid) chunks contribute none, they were never cut.
    """
    his: dict[str, set[int]] = {}
    los: dict[str, set[int]] = {}
    for chunk in plan.chunks:
        for region in chunk.regions:
            seqid, lo, hi = _region_bounds(region)
            if lo is None or hi is None:
                continue
            los.setdefault(seqid, set()).add(lo)
            his.setdefault(seqid, set()).add(hi)
    out: dict[str, list[int]] = {}
    for seqid, his_set in his.items():
        interior = sorted(his_set & los.get(seqid, set()))
        if interior:
            out[seqid] = interior
    return out


def _bridged(left: Any, right: Any, introns: list[tuple[int, int]]) -> bool:
    """True if a verified intron spans the gap ``[left.end, right.start)``.

    Reuses the reconciler's exact merge predicate
    (``_has_bridging_introns``) so a stitched merge is decided identically to an
    in-chunk one.
    """
    return _has_bridging_introns([left, right], _PseudoLocus(introns))  # type: ignore[arg-type]  # _PseudoLocus duck-types MikadoLocus


# ---------------------------------------------------------------------------
# Candidate detection
# ---------------------------------------------------------------------------


def find_boundary_merge_candidates(
    genes: list[Any],
    plan: Plan,
    junctions: list[Any],
    *,
    flank: int | None = None,
) -> list[tuple[str, str]]:
    """Return ``[(left_gene_id, right_gene_id)]`` recoverable across boundaries.

    A pair qualifies when both genes sit within ``flank`` of the same chunk
    boundary on opposite sides, share strand, are disjoint and ordered, and a
    same-strand verified junction bridges the gap between them. ``flank`` defaults
    to the plan's ``flank``. Pure detection (no mutation); reused by the driver to
    *count* recoverable merges off a reloaded GFF3 and by :func:`boundary_stitch`
    to *apply* them.
    """
    flank = plan.flank if flank is None else flank
    boundaries = _plan_boundaries(plan)
    if not boundaries:
        return []

    jx: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for j in junctions:
        jx.setdefault((j.seqid, j.strand), []).append((j.donor, j.acceptor))

    by_seqid: dict[str, list[Any]] = {}
    for g in genes:
        by_seqid.setdefault(g.seqid, []).append(g)

    pairs: list[tuple[str, str]] = []
    consumed: set[str] = set()
    for seqid, cuts in boundaries.items():
        scaf_genes = sorted(by_seqid.get(seqid, []), key=lambda g: g.start)
        introns = jx.get
        for p in cuts:
            left = [g for g in scaf_genes if p - flank <= g.end <= p]
            right = [g for g in scaf_genes if p <= g.start <= p + flank]
            for gl in left:
                if gl.gene_id in consumed:
                    continue
                for gr in right:
                    if gr.gene_id in consumed or gr.strand != gl.strand:
                        continue
                    if gl.end >= gr.start:  # need a real gap to bridge
                        continue
                    strand_introns = introns((seqid, gl.strand), [])
                    if strand_introns and _bridged(gl, gr, strand_introns):
                        pairs.append((gl.gene_id, gr.gene_id))
                        consumed.update((gl.gene_id, gr.gene_id))
                        break
    return pairs


# ---------------------------------------------------------------------------
# Merge application
# ---------------------------------------------------------------------------


def _order_by_id(
    a: ReconciledGene,
    b: ReconciledGene,
) -> tuple[ReconciledGene, ReconciledGene]:
    """``(rep, other)`` keeping the lowest HFG number as the merge rep (§11)."""
    try:
        na: Any
        nb: Any
        na, nb = _hfg_num(a.gene_id), _hfg_num(b.gene_id)
    except (IndexError, ValueError):
        na, nb = a.gene_id, b.gene_id
    return (a, b) if na <= nb else (b, a)


def _merge_two(gl: ReconciledGene, gr: ReconciledGene) -> ReconciledGene:
    """Fuse two straddling genes into one (lowest HFG kept; ``LOCUS_MERGE``)."""
    rep, other = _order_by_id(gl, gr)
    combined = [*rep.transcripts, *other.transcripts]
    # Preserve a TRaCE-elected ordering (Phase 33b D2): when the chunk runs elected
    # the canonical transcript (``trace_rank`` set), re-sorting the fused gene by
    # ``combined_score`` would discard that election and disagree with a whole-
    # genome run. Concatenating rep-then-other keeps each gene's elected order and
    # the representative's elected primary as the merged primary; with no TRaCE
    # ranks (the default) ``order=None`` reproduces the historical score sort.
    order = (
        [t.transcript_id for t in combined]
        if any(t.trace_rank is not None for t in combined)
        else None
    )
    transcripts = _renumber(combined, rep.gene_id, order=order)
    primary = transcripts[0]
    merged_from = list(rep.merged_from)
    for gid in (other.gene_id, *other.merged_from):
        if gid not in merged_from:
            merged_from.append(gid)
    return attrs.evolve(
        rep,
        start=min(gl.start, gr.start),
        end=max(gl.end, gr.end),
        transcripts=transcripts,
        primary_transcript_id=primary.transcript_id,
        merged_from=merged_from,
        flags=dedup_flags([*rep.flags, *other.flags, LOCUS_MERGE]),
    )


def boundary_stitch(
    genes: list[ReconciledGene],
    plan: Plan,
    junctions: list[Any],
    *,
    flank: int | None = None,
) -> tuple[list[ReconciledGene], int]:
    """Merge cross-boundary pairs in ``genes``; return ``(stitched, n_recovered)``.

    A no-op (returns the list unchanged, ``0``) when no boundary pair is bridged.
    The output is re-sorted by ``(seqid, start)`` so it stays in genomic order.
    """
    pairs = find_boundary_merge_candidates(genes, plan, junctions, flank=flank)
    if not pairs:
        return list(genes), 0

    by_id: dict[str, ReconciledGene] = {g.gene_id: g for g in genes}
    consumed: set[str] = set()
    merged: list[ReconciledGene] = []
    for left_id, right_id in pairs:
        if left_id in consumed or right_id in consumed:
            continue
        merged.append(_merge_two(by_id[left_id], by_id[right_id]))
        consumed.update((left_id, right_id))

    out = [g for g in genes if g.gene_id not in consumed]
    out.extend(merged)
    out.sort(key=lambda g: (g.seqid, g.start))
    _log.info("boundary-stitch recovered %d cross-chunk merge(s)", len(merged))
    return out, len(merged)
