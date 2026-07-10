"""Strategy-based genome chunking."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Callable

from helixforge.constants import PLAN_DEFAULT_MIN_BOUNDARY_GAP
from helixforge.utils.logging import get_logger
from helixforge.utils.regions import format_region

if TYPE_CHECKING:
    from helixforge.reconcile.models import HelixerLocus

_log = get_logger(__name__)

# Defaults mirror v1's chunker (``create_plan``): 10 Mb windows for "size", 100
# loci for "genes", 100 kb minimum chunk for the bp-based strategies.
DEFAULT_SIZE_CHUNK_BP = 10_000_000
DEFAULT_GENES_PER_CHUNK = 100
DEFAULT_MIN_CHUNK_BP = 100_000


class ChunkStrategy(str, Enum):
    """The four v1 chunking strategies (value = the ``--strategy`` flag word)."""

    BY_SCAFFOLD = (
        "scaffold"  # one chunk per scaffold (split big scaffolds via --max-chunk-size)
    )
    BY_SIZE = "size"  # fixed base-pair windows (snapped to inter-locus gaps)
    BY_GENES = "genes"  # fixed number of loci per chunk
    ADAPTIVE = "adaptive"  # balance loci across ~target_chunks chunks


@dataclass
class ChunkSpec:
    """One chunk's region + the Helixer loci it owns (pre-id-reservation)."""

    region: str
    locus_ids: list[str] = field(default_factory=list)

    @property
    def num_loci(self) -> int:
        return len(self.locus_ids)


def _resolve_target(
    strategy: ChunkStrategy,
    total_loci: int,
    chunk_size: int | None,
    target_chunks: int | None,
) -> tuple[int | None, int | None]:
    """Return ``(genes_per_chunk, bp_per_chunk)`` for the chosen strategy.

    Exactly one of the two is set (the other ``None``); ``scaffold`` with no
    ``max_chunk_size`` leaves both ``None`` (never cut a scaffold).
    """
    if strategy is ChunkStrategy.BY_GENES:
        return (chunk_size or DEFAULT_GENES_PER_CHUNK), None
    if strategy is ChunkStrategy.BY_SIZE:
        return None, (chunk_size or DEFAULT_SIZE_CHUNK_BP)
    if strategy is ChunkStrategy.ADAPTIVE:
        if target_chunks is None or target_chunks < 1:
            raise ValueError("adaptive strategy needs --target-chunks >= 1")
        # ceil division: aim for the requested chunk count without overshooting.
        return max(1, -(-total_loci // target_chunks)), None
    return None, None  # BY_SCAFFOLD: cut only on --max-chunk-size (handled below)


def _segment_scaffold(
    loci_sorted: list[HelixerLocus],
    scaf_len: int,
    *,
    min_boundary_gap: int,
    min_span_bp: int,
    should_cut: Callable[[list[HelixerLocus], int], bool],
) -> list[ChunkSpec]:
    """Split one scaffold's sorted loci into gene-respecting segments.

    A cut is taken **after** the current locus when the strategy says the running
    segment is "full" (``should_cut``), the segment already spans at least
    ``min_span_bp`` bp, **and** the gap to the next locus is ``>= min_boundary_gap``
    (so the cut lands in intergenic space). The genomic boundary is the gap
    midpoint. A scaffold with no eligible gap stays one chunk, a gene is never
    split.
    """
    specs: list[ChunkSpec] = []
    current: list[HelixerLocus] = [loci_sorted[0]]
    seg_lo = 0
    for prev, nxt in zip(loci_sorted, loci_sorted[1:]):
        gap = nxt.start - prev.end
        seg_span = prev.end - seg_lo
        if (
            gap >= min_boundary_gap
            and seg_span >= min_span_bp
            and should_cut(current, seg_span)
        ):
            hi = prev.end + gap // 2
            specs.append(
                ChunkSpec(
                    format_region(loci_sorted[0].seqid, seg_lo, hi),
                    [g.gene_id for g in current],
                )
            )
            current = [nxt]
            seg_lo = hi
        else:
            current.append(nxt)
    specs.append(
        ChunkSpec(
            format_region(loci_sorted[0].seqid, seg_lo, scaf_len),
            [g.gene_id for g in current],
        )
    )
    return specs


def plan_chunks(
    loci: list[HelixerLocus],
    scaffold_lengths: dict[str, int],
    *,
    strategy: ChunkStrategy | str = ChunkStrategy.BY_SCAFFOLD,
    chunk_size: int | None = None,
    min_chunk_size: int = DEFAULT_MIN_CHUNK_BP,
    max_chunk_size: int | None = None,
    target_chunks: int | None = None,
    min_boundary_gap: int = PLAN_DEFAULT_MIN_BOUNDARY_GAP,
) -> list[ChunkSpec]:
    """Partition ``loci`` into gene-respecting :class:`ChunkSpec`\\ s by strategy.

    ``loci`` are v3 master Helixer loci; ``scaffold_lengths`` bounds the last
    region per scaffold. The four strategies match v1's vocabulary:

    * ``scaffold``, one chunk per scaffold; with ``max_chunk_size`` a long
      scaffold is split into ``~max_chunk_size`` bp pieces (cut in gaps).
    * ``size``, ``~chunk_size`` bp windows (default 10 Mb), cut in gaps.
    * ``genes``, ``chunk_size`` loci per chunk (default 100).
    * ``adaptive``, balance loci across ``~target_chunks`` chunks.

    Every cut respects ``min_boundary_gap`` so no gene is split; chunks are
    emitted in scaffold order. Scaffolds with no loci produce no chunk.
    """
    strategy = ChunkStrategy(strategy)
    if min_boundary_gap < 1:
        raise ValueError("min_boundary_gap must be >= 1")
    total = len(loci)
    genes_per_chunk, bp_per_chunk = _resolve_target(
        strategy, total, chunk_size, target_chunks
    )

    # The "should we cut after this segment" predicate + the bp floor that keeps
    # the bp-based strategies from making tiny chunks (count-based strategies
    # ignore min_chunk_size, so a short scaffold can still be cut by locus count).
    def _never_cut(_l: list[HelixerLocus], _s: int) -> bool:
        return False

    def _cut_by_span(_loci: list[HelixerLocus], span: int) -> bool:
        return span >= (max_chunk_size or 0)

    def _cut_by_bp(_loci: list[HelixerLocus], span: int) -> bool:
        return span >= (bp_per_chunk or 0)

    def _cut_by_genes(seg: list[HelixerLocus], _span: int) -> bool:
        return len(seg) >= (genes_per_chunk or 0)

    if strategy is ChunkStrategy.BY_SCAFFOLD:
        if max_chunk_size is None:
            should_cut: Callable[[list[HelixerLocus], int], bool] = _never_cut
            min_span = 0
        else:
            should_cut = _cut_by_span
            min_span = min(min_chunk_size, max_chunk_size)
    elif strategy is ChunkStrategy.BY_SIZE:
        should_cut = _cut_by_bp
        min_span = min(min_chunk_size, bp_per_chunk or min_chunk_size)
    else:  # BY_GENES / ADAPTIVE, count-based, no bp floor
        should_cut = _cut_by_genes
        min_span = 0

    by_scaffold: dict[str, list[HelixerLocus]] = {}
    for locus in loci:
        by_scaffold.setdefault(locus.seqid, []).append(locus)

    specs: list[ChunkSpec] = []
    for seqid in sorted(by_scaffold):
        scaf_loci = sorted(by_scaffold[seqid], key=lambda g: g.start)
        scaf_len = scaffold_lengths.get(seqid)
        if scaf_len is None:
            scaf_len = max(g.end for g in scaf_loci)
            _log.warning(
                "scaffold %s absent from index; bounding region at %d", seqid, scaf_len
            )
        specs.extend(
            _segment_scaffold(
                scaf_loci,
                scaf_len,
                min_boundary_gap=min_boundary_gap,
                min_span_bp=min_span,
                should_cut=should_cut,
            )
        )

    _log.info(
        "strategy=%s partitioned %d loci over %d scaffolds into %d chunks "
        "(chunk_size=%s, max_chunk_size=%s, target_chunks=%s, "
        "min_boundary_gap=%d)",
        strategy.value,
        total,
        len(by_scaffold),
        len(specs),
        chunk_size,
        max_chunk_size,
        target_chunks,
        min_boundary_gap,
    )
    return specs
