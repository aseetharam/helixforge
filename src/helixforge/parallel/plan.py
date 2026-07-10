"""Genome partitioning + disjoint HFG-range reservation."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from helixforge.reconcile.models import HelixerLocus

from helixforge.constants import NOVEL_ID_BASE
from helixforge.constants import (
    PLAN_DEFAULT_MIN_BOUNDARY_GAP as DEFAULT_MIN_BOUNDARY_GAP,
)
from helixforge.constants import PLAN_MAX_SCAFFOLDS_PER_CHUNK, PLAN_SMALL_SCAFFOLD_BP
from helixforge.io.fasta import GenomeAccessor
from helixforge.reconcile.locus import load_helixer_loci
from helixforge.reconcile.mikado_integrate import _hfg_num
from helixforge.utils.logging import get_logger
from helixforge.utils.regions import format_region

_log = get_logger(__name__)

# DEFAULT_MIN_BOUNDARY_GAP re-exported from helixforge.constants (Phase 19,
# §1.6): heuristic floor for a cut-eligible inter-locus gap when the caller gives
# none. Wide enough to clear a typical gene/Mikado-merge bridging span; raise it
# to trade parallelism for fewer split merges (see module docstring).


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class Chunk:
    """One independent pipeline unit: a region (or few) + its reserved IDs."""

    chunk_id: str
    regions: list[str]  # region strings (internal 0-based coords)
    num_loci: int
    locus_ids: list[str] = field(default_factory=list)  # Helixer gene_ids in this chunk
    id_base: int = 1
    id_range: tuple[int, int] = (1, 1)  # [lo, hi) HFG numbers reserved (disjoint)
    novel_base: int = NOVEL_ID_BASE
    novel_range: tuple[int, int] = (NOVEL_ID_BASE, NOVEL_ID_BASE)
    est_resources: dict[str, int] = field(default_factory=dict)


@dataclass
class Plan:
    """A whole-genome partition: the chunk list + how it was built."""

    chunks: list[Chunk]
    min_boundary_gap: int = DEFAULT_MIN_BOUNDARY_GAP
    flank: int = 200
    total_loci: int = 0
    strategy: str = "scaffold"  # v1 chunking strategy used (provenance)

    def __len__(self) -> int:
        return len(self.chunks)


# ---------------------------------------------------------------------------
# Scaffold lengths
# ---------------------------------------------------------------------------


def _read_scaffold_lengths(genome_fai: str | Path) -> dict[str, int]:
    """``{seqid: length}`` from a samtools ``.fai`` index or a FASTA path.

    A ``.fai`` (or any two+-column ``name<TAB>length...`` table) is read directly;
    anything else is opened via :class:`GenomeAccessor` (which builds/uses the
    ``.fai``). Lengths only bound the last region per scaffold, so exactness
    beyond the locus span is irrelevant.
    """
    path = Path(genome_fai)
    name = path.name
    if name.endswith(".fai") or name.endswith(".fasta.fai") or name.endswith(".fa.fai"):
        lengths: dict[str, int] = {}
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                cols = line.split("\t")
                lengths[cols[0]] = int(cols[1])
        return lengths
    with GenomeAccessor(str(path)) as genome:
        return genome.get_scaffold_lengths()


# ---------------------------------------------------------------------------
# Partitioning
# ---------------------------------------------------------------------------


def _target_loci_per_chunk(
    total_loci: int,
    target_chunks: int | None,
    target_loci_per_chunk: int | None,
) -> int | None:
    """Resolve the desired loci-per-chunk; ``None`` means do not cut scaffolds."""
    if target_loci_per_chunk is not None:
        if target_loci_per_chunk < 1:
            raise ValueError("target_loci_per_chunk must be >= 1")
        return target_loci_per_chunk
    if target_chunks is not None:
        if target_chunks < 1:
            raise ValueError("target_chunks must be >= 1")
        # ceil division so we never *over*-shoot the requested chunk count.
        return max(1, -(-total_loci // target_chunks))
    return None  # one chunk per scaffold


def _segment_scaffold(
    loci_sorted: list[HelixerLocus],
    target: int | None,
    min_gap: int,
) -> list[tuple[list[HelixerLocus], int | None, int | None]]:
    """Split one scaffold's sorted loci into segments, cutting only in big gaps.

    A cut is made after the current locus when (a) the running segment already
    holds ``target`` loci and (b) the gap to the next locus is ``>= min_gap``. A
    scaffold with no eligible gap (or ``target is None``) stays a single segment,
    so a gene is never split.
    """
    segments: list[tuple[list[HelixerLocus], int | None, int | None]] = []
    current: list[HelixerLocus] = [loci_sorted[0]]
    for prev, nxt in zip(loci_sorted, loci_sorted[1:]):
        gap = nxt.start - prev.end
        eligible = gap >= min_gap
        if target is not None and len(current) >= target and eligible:
            segments.append((current, prev.end, nxt.start))  # carry the cut gap
            current = [nxt]
        else:
            current.append(nxt)
    segments.append((current, None, None))
    return segments


def partition_genome(
    genome_fai: str | Path,
    helixer_gff3: str | Path,
    *,
    target_chunks: int | None = None,
    target_loci_per_chunk: int | None = None,
    min_boundary_gap: int | None = None,
    flank: int = 200,
    pack_small_scaffolds: bool = False,
    small_scaffold_bp: int = PLAN_SMALL_SCAFFOLD_BP,
    max_scaffolds_per_chunk: int = PLAN_MAX_SCAFFOLDS_PER_CHUNK,
) -> Plan:
    """Partition a genome into chunks that never split a gene.

    Loci are the **master** Helixer loci (after same-strand merge, Phase 3), so a
    fused over-split locus stays whole. Each scaffold is cut independently, only
    in inter-locus gaps ``>= min_boundary_gap`` (default
    ``max(flank, DEFAULT_MIN_BOUNDARY_GAP)``), targeting roughly
    ``target_loci_per_chunk`` loci per chunk; a scaffold with no large-enough gap
    becomes one chunk. Returns a :class:`Plan`; call :func:`reserve_id_ranges`
    next to fill the HFG ranges.

    **Small-scaffold bin-packing** — when
    ``pack_small_scaffolds`` is set, every *whole, uncut* scaffold at or below
    ``small_scaffold_bp`` is packed together with other small scaffolds into
    combined chunks (each carrying several bare-seqid ``regions``), so a
    fragmented draft assembly of 10^5-10^6 tiny contigs does not explode into one
    chunk per scaffold (which would blow past Slurm ``MaxArraySize``). A bin is
    flushed when it would exceed ``target_loci_per_chunk`` loci (when a target is
    set) or reaches ``max_scaffolds_per_chunk`` scaffolds. A gene is **never**
    split: only single-segment scaffolds are packed, and each packed scaffold's
    every locus still lies wholly in exactly one chunk. Large/cut scaffolds keep
    one region per chunk as before. Default ``False`` reproduces the historical
    one-chunk-per-scaffold behaviour exactly.

    See the module docstring for the residual cross-boundary-merge limitation.
    """
    if min_boundary_gap is None:
        min_boundary_gap = max(flank, DEFAULT_MIN_BOUNDARY_GAP)
    if min_boundary_gap < 1:
        raise ValueError("min_boundary_gap must be >= 1")
    if max_scaffolds_per_chunk < 1:
        raise ValueError("max_scaffolds_per_chunk must be >= 1")

    lengths = _read_scaffold_lengths(genome_fai)
    loci: list[HelixerLocus] = load_helixer_loci(str(helixer_gff3))
    total = len(loci)
    target = _target_loci_per_chunk(total, target_chunks, target_loci_per_chunk)

    by_scaffold: dict[str, list[HelixerLocus]] = {}
    for locus in loci:
        by_scaffold.setdefault(locus.seqid, []).append(locus)

    chunks: list[Chunk] = []
    # Small, single-segment scaffolds deferred for bin-packing.
    packable: list[tuple[str, list[HelixerLocus]]] = []
    for seqid in sorted(by_scaffold):
        scaf_loci = sorted(by_scaffold[seqid], key=lambda g: g.start)
        scaf_len = lengths.get(seqid)
        if scaf_len is None:
            # No length: fall back to the locus extent (region still covers them).
            scaf_len = max(g.end for g in scaf_loci)
            _log.warning(
                "scaffold %s absent from fai; bounding region at %d", seqid, scaf_len
            )
        segments = _segment_scaffold(scaf_loci, target, min_boundary_gap)
        n_seg = len(segments)

        # A whole, uncut, small scaffold is eligible for packing (deferred so the
        # bins fill deterministically after the cut scaffolds are emitted).
        if pack_small_scaffolds and n_seg == 1 and scaf_len <= small_scaffold_bp:
            packable.append((seqid, scaf_loci))
            continue

        # Cut points: midpoint of each carried gap; first lo=0, last hi=scaf_len.
        lo = 0
        for si, (seg_loci, gap_lo, gap_hi) in enumerate(segments):
            if si < n_seg - 1:
                # gap_lo/gap_hi are always ints for non-last segments
                hi = cast(int, gap_lo) + (cast(int, gap_hi) - cast(int, gap_lo)) // 2
            else:
                hi = scaf_len
            # A single uncut scaffold is addressed by its bare seqid (whole
            # scaffold); a true sub-range carries explicit internal coordinates.
            if n_seg == 1:
                region = seqid
            else:
                region = format_region(seqid, lo, hi)
            chunk_id = f"chunk_{len(chunks):04d}"
            chunks.append(
                Chunk(
                    chunk_id=chunk_id,
                    regions=[region],
                    num_loci=len(seg_loci),
                    locus_ids=[g.gene_id for g in seg_loci],
                    est_resources=_est_resources(len(seg_loci)),
                )
            )
            lo = hi

    # Flush the deferred small scaffolds into combined multi-region chunks.
    _pack_into_chunks(chunks, packable, target, max_scaffolds_per_chunk)

    plan = Plan(
        chunks=chunks, min_boundary_gap=min_boundary_gap, flank=flank, total_loci=total
    )
    _validate_partition(plan, loci)
    _log.info(
        "partitioned %d loci over %d scaffolds into %d chunks "
        "(target_loci_per_chunk=%s, min_boundary_gap=%d, packed=%d small "
        "scaffolds)",
        total,
        len(by_scaffold),
        len(chunks),
        target,
        min_boundary_gap,
        len(packable),
    )
    return plan


def partition_by_strategy(
    genome_fai: str | Path,
    helixer_gff3: str | Path,
    *,
    strategy: str = "scaffold",
    chunk_size: int | None = None,
    min_chunk_size: int = 100_000,
    max_chunk_size: int | None = None,
    target_chunks: int | None = None,
    min_boundary_gap: int | None = None,
    flank: int = 200,
) -> Plan:
    """Partition a genome with v1's chunking strategies (gene-respecting).

    This is the CLI ``parallel plan`` entry point: it adopts v1's transparent
    ``--strategy {scaffold,size,genes,adaptive}`` vocabulary
    (:mod:`helixforge.parallel.chunker`) while keeping v3's two guarantees — no
    gene is ever split across a chunk (every cut lands in an inter-locus gap
    ``>= min_boundary_gap``), and each chunk owns a disjoint set of master loci so
    :func:`reserve_id_ranges` can hand out globally-unique HFG ranges.

    Loci are the **master** Helixer loci (after same-strand merge), so a fused
    over-split locus stays whole. Returns a :class:`Plan`; call
    :func:`reserve_id_ranges` next to fill the HFG ranges.
    """
    from helixforge.parallel.chunker import ChunkStrategy, plan_chunks

    if min_boundary_gap is None:
        min_boundary_gap = max(flank, DEFAULT_MIN_BOUNDARY_GAP)
    if min_boundary_gap < 1:
        raise ValueError("min_boundary_gap must be >= 1")

    lengths = _read_scaffold_lengths(genome_fai)
    loci: list[HelixerLocus] = load_helixer_loci(str(helixer_gff3))
    specs = plan_chunks(
        loci,
        lengths,
        strategy=ChunkStrategy(strategy),
        chunk_size=chunk_size,
        min_chunk_size=min_chunk_size,
        max_chunk_size=max_chunk_size,
        target_chunks=target_chunks,
        min_boundary_gap=min_boundary_gap,
    )
    chunks = [
        Chunk(
            chunk_id=f"chunk_{i:04d}",
            regions=[spec.region],
            num_loci=spec.num_loci,
            locus_ids=list(spec.locus_ids),
            est_resources=_est_resources(spec.num_loci),
        )
        for i, spec in enumerate(specs)
    ]
    plan = Plan(
        chunks=chunks,
        min_boundary_gap=min_boundary_gap,
        flank=flank,
        total_loci=len(loci),
        strategy=ChunkStrategy(strategy).value,
    )
    _validate_partition(plan, loci)
    _log.info(
        "strategy=%s → %d chunks over %d loci", plan.strategy, len(chunks), len(loci)
    )
    return plan


def _pack_into_chunks(
    chunks: list[Chunk],
    packable: list[tuple[str, list[HelixerLocus]]],
    target: int | None,
    max_scaffolds_per_chunk: int,
) -> None:
    """Greedily bin-pack whole small scaffolds into combined chunks.

    Appends combined chunks (each with several bare-seqid ``regions``) to
    ``chunks`` in scaffold order. A bin is flushed before adding a scaffold when
    the bin already holds at least one scaffold **and** either adding the next
    would push it past ``target`` loci (when a target is set) or it already holds
    ``max_scaffolds_per_chunk`` scaffolds. Each whole scaffold is added intact, so
    a gene is never split and every locus stays in exactly one chunk.
    """
    cur_regions: list[str] = []
    cur_ids: list[str] = []
    cur_loci = 0

    def _flush() -> None:
        nonlocal cur_regions, cur_ids, cur_loci
        if not cur_regions:
            return
        chunks.append(
            Chunk(
                chunk_id=f"chunk_{len(chunks):04d}",
                regions=list(cur_regions),
                num_loci=cur_loci,
                locus_ids=list(cur_ids),
                est_resources=_est_resources(cur_loci),
            )
        )
        cur_regions, cur_ids, cur_loci = [], [], 0

    for seqid, scaf_loci in packable:
        n = len(scaf_loci)
        over_target = target is not None and cur_loci + n > target
        over_count = len(cur_regions) >= max_scaffolds_per_chunk
        if cur_regions and (over_target or over_count):
            _flush()
        cur_regions.append(seqid)  # bare seqid = whole scaffold
        cur_ids.extend(g.gene_id for g in scaf_loci)
        cur_loci += n
    _flush()


def _est_resources(num_loci: int) -> dict[str, int]:
    """Light per-chunk resource hint (heuristic; refine with ``suggest``)."""
    return {
        "num_loci": num_loci,
        "mem_mb": max(2048, num_loci * 4),  # ~4 MB/locus floor, min 2 GB
    }


def _validate_partition(plan: Plan, loci: list[HelixerLocus]) -> None:
    """Assert no boundary splits a locus and every locus is covered exactly once.

    Each locus must fall wholly within exactly one chunk region. Region strings
    are internal 0-based half-open; a bare seqid means the whole scaffold. Chunks
    may now hold **multiple regions** (small-scaffold bin-packing):
    the per-region coverage bookkeeping handles that transparently, and a
    bin-packed whole scaffold must still be owned by exactly one chunk (a scaffold
    appearing as a bare seqid in two chunks would surface as a locus with two
    owners below).
    """
    covers: dict[str, list[tuple[int | None, int | None, str]]] = {}
    # A whole-scaffold (bare seqid) region must be the *only* region for its
    # seqid — a packed scaffold can never also be cut/packed elsewhere.
    whole_scaffold_owner: dict[str, str] = {}
    for chunk in plan.chunks:
        for region in chunk.regions:
            seqid, lo, hi = _region_bounds(region)
            covers.setdefault(seqid, []).append((lo, hi, chunk.chunk_id))
            if lo is None:
                if seqid in whole_scaffold_owner:
                    raise ValueError(
                        f"scaffold {seqid} appears as a whole-scaffold region in "
                        f"two chunks ({whole_scaffold_owner[seqid]}, "
                        f"{chunk.chunk_id}) — a packed scaffold must be owned once"
                    )
                whole_scaffold_owner[seqid] = chunk.chunk_id

    seen: dict[str, str] = {}
    for locus in loci:
        owners = []
        for lo, hi, cid in covers.get(locus.seqid, []):
            if lo is None or hi is None:
                within = True
                overlaps = True
            else:
                within = lo <= locus.start and locus.end <= hi
                overlaps = locus.start < hi and locus.end > lo
            if overlaps and not within:
                raise ValueError(
                    f"partition boundary splits locus {locus.gene_id} "
                    f"({locus.seqid}:{locus.start}-{locus.end}) in chunk {cid}"
                )
            if within:
                owners.append(cid)
        if len(owners) != 1:
            raise ValueError(
                f"locus {locus.gene_id} covered by {len(owners)} chunks "
                f"(expected exactly 1): {owners}"
            )
        seen[locus.gene_id] = owners[0]
    if len(seen) != len(loci):
        raise ValueError("internal: locus coverage count mismatch")


def _region_bounds(region: str) -> tuple[str, int | None, int | None]:
    """``(seqid, lo, hi)`` for a region string; bare seqid → ``(seqid, None, None)``."""
    if ":" not in region:
        return region, None, None
    seqid, _, span = region.rpartition(":")
    lo_s, _, hi_s = span.partition("-")
    return seqid, int(lo_s), int(hi_s)


# ---------------------------------------------------------------------------
# HFG range reservation
# ---------------------------------------------------------------------------


def reserve_id_ranges(
    plan: Plan,
    id_map: dict[str, str] | None = None,
) -> Plan:
    """Assign each chunk a disjoint, contiguous HFG range (mutates + returns plan).

    Ranges are laid end to end and sized to each chunk's Helixer-locus count, so
    a chunk can never run out of numbers (one number per locus; splits use letter
    suffixes, merges reuse the lowest existing id). They start **above** any
    number already in the master ``id_map`` (separately for the Helixer-anchored
    block and the ``NOVEL_ID_BASE`` novel block), so a chunk's fresh allocations
    can never collide with an HFG the master already handed to a locus now sitting
    in another chunk. Already-mapped loci keep their HFG regardless (the allocator
    consults the seed map first), giving rerun stability.

    Asserts the reserved ranges (and novel sub-ranges) are pairwise disjoint.
    """
    id_map = id_map or {}
    nums: list[int] = []
    for value in id_map.values():
        try:
            nums.append(_hfg_num(value))
        except (IndexError, ValueError):
            continue
    base_used = [n for n in nums if n < NOVEL_ID_BASE]
    novel_used = [n for n in nums if n >= NOVEL_ID_BASE]
    base_cursor = (max(base_used) + 1) if base_used else 1
    novel_cursor = (max(novel_used) + 1) if novel_used else NOVEL_ID_BASE

    for chunk in plan.chunks:
        size = max(chunk.num_loci, 1)
        chunk.id_base = base_cursor
        chunk.id_range = (base_cursor, base_cursor + size)
        base_cursor += size

        chunk.novel_base = novel_cursor
        chunk.novel_range = (novel_cursor, novel_cursor + size)
        novel_cursor += size

    if base_cursor > NOVEL_ID_BASE:
        _log.warning(
            "reserved Helixer-anchored range reaches %d, at/over the novel base "
            "%d — raise NOVEL_ID_BASE or reduce the genome's locus count",
            base_cursor,
            NOVEL_ID_BASE,
        )
    _assert_ranges_disjoint(plan)
    _log.info(
        "reserved HFG ranges for %d chunks: base [%d, %d), novel from %d",
        len(plan.chunks),
        plan.chunks[0].id_base if plan.chunks else 0,
        base_cursor,
        NOVEL_ID_BASE,
    )
    return plan


def _assert_ranges_disjoint(plan: Plan) -> None:
    """Pairwise-disjointness check over both the base and novel ranges."""
    for label, attr in (("HFG", "id_range"), ("novel", "novel_range")):
        ordered = sorted(plan.chunks, key=lambda c: getattr(c, attr)[0])
        for a, b in zip(ordered, ordered[1:]):
            a_hi = getattr(a, attr)[1]
            b_lo = getattr(b, attr)[0]
            if b_lo < a_hi:
                raise ValueError(
                    f"{label} ranges overlap: {a.chunk_id} {getattr(a, attr)} "
                    f"vs {b.chunk_id} {getattr(b, attr)}"
                )


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------


def write_plan(plan: Plan, path: str | Path) -> Path:
    """Write a :class:`Plan` to JSON (chunks as dicts; tuples become lists)."""
    payload = {
        "strategy": plan.strategy,
        "min_boundary_gap": plan.min_boundary_gap,
        "flank": plan.flank,
        "total_loci": plan.total_loci,
        "chunks": [asdict(c) for c in plan.chunks],
    }
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=False))
    return Path(path)


def read_plan(path: str | Path) -> Plan:
    """Read a :class:`Plan` back from JSON (lists become tuples for ranges)."""
    payload = json.loads(Path(path).read_text())
    chunks: list[Chunk] = []
    for c in payload["chunks"]:
        chunks.append(
            Chunk(
                chunk_id=c["chunk_id"],
                regions=list(c["regions"]),
                num_loci=c["num_loci"],
                locus_ids=list(c.get("locus_ids", [])),
                id_base=c["id_base"],
                id_range=cast(tuple[int, int], tuple(c["id_range"])),
                novel_base=c["novel_base"],
                novel_range=cast(tuple[int, int], tuple(c["novel_range"])),
                est_resources=dict(c.get("est_resources", {})),
            )
        )
    return Plan(
        chunks=chunks,
        min_boundary_gap=payload.get("min_boundary_gap", DEFAULT_MIN_BOUNDARY_GAP),
        flank=payload.get("flank", 200),
        total_loci=payload.get("total_loci", len(chunks)),
        strategy=payload.get("strategy", "scaffold"),
    )
