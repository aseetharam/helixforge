"""Granularity + resource recommendations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# Heuristic knobs (labelled; tune for your cluster / genome) re-exported from
# helixforge.constants (Phase 19, §1.6). The boundary-gap default is suggest's
# own (conservative; clears most Mikado merges), distinct from plan's floor.
from helixforge.constants import (
    MEM_GB_PER_MB_SEQUENCE,
    MIN_CHUNK_MEM_GB,
    MIN_WALLTIME_MIN,
    SUGGEST_DEFAULT_MIN_BOUNDARY_GAP as DEFAULT_MIN_BOUNDARY_GAP,
    TARGET_CHUNK_BP,
    WALLTIME_MIN_PER_MB,
)
from helixforge.parallel.plan import _read_scaffold_lengths
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Genome stats
# ---------------------------------------------------------------------------


@dataclass
class GenomeStats:
    genome_size: int
    scaffold_count: int
    n50: int
    largest: int


def _genome_stats(genome_fai: str) -> GenomeStats:
    lengths = list(_read_scaffold_lengths(genome_fai).values())
    if not lengths:
        raise ValueError("no scaffolds found in the genome index")
    total = sum(lengths)
    lengths_desc = sorted(lengths, reverse=True)
    half = total / 2
    running = 0
    n50 = lengths_desc[-1]
    for length in lengths_desc:
        running += length
        if running >= half:
            n50 = length
            break
    return GenomeStats(
        genome_size=total,
        scaffold_count=len(lengths),
        n50=n50,
        largest=lengths_desc[0],
    )


# ---------------------------------------------------------------------------
# Suggestion
# ---------------------------------------------------------------------------


@dataclass
class Suggestion:
    stats: GenomeStats
    target_chunks: int
    target_loci_per_chunk: int | None
    min_boundary_gap: int
    procs: int
    threads: int
    mem_gb: int
    walltime_min: int
    rationale: list[str] = field(default_factory=list)

    @property
    def walltime(self) -> str:
        """Walltime as a Slurm ``HH:MM:00`` string."""
        hours, minutes = divmod(self.walltime_min, 60)
        return f"{hours:02d}:{minutes:02d}:00"

    def render(self) -> str:
        """A human-readable recommendation block (heuristics, with trade-offs)."""
        s = self.stats
        lines = [
            "HelixForge parallel — suggested plan (HEURISTIC, not a guarantee)",
            f"  genome size      : {s.genome_size:,} bp",
            f"  scaffolds        : {s.scaffold_count:,} (N50 {s.n50:,}, largest {s.largest:,})",
            "",
            f"  target_chunks    : {self.target_chunks}",
            f"  min_boundary_gap : {self.min_boundary_gap}",
            f"  per-chunk procs  : {self.procs}",
            f"  per-chunk threads: {self.threads}",
            f"  per-chunk mem    : {self.mem_gb} GB",
            f"  per-chunk walltime: {self.walltime}",
            "",
            "  paste into `parallel plan` (numbers are directly usable as flags):",
            "    helixforge parallel plan --genome GENOME --gff HELIXER.gff3 \\",
            f"        --strategy adaptive --target-chunks {self.target_chunks} \\",
            f"        --min-boundary-gap {self.min_boundary_gap} -o plan.json",
            "",
            "  trade-offs:",
        ]
        lines += [f"    - {note}" for note in self.rationale]
        return "\n".join(lines)


def suggest_plan(
    genome_fai: str, *, hpc_profile: dict[str, Any] | None = None
) -> Suggestion:
    """Recommend chunking + resources from a genome index and HPC profile.

    ``hpc_profile`` is a dict (all keys optional): ``cores_per_node``,
    ``mem_gb_per_node``, ``max_array_size``, ``walltime_cap_hours``. Returns a
    :class:`Suggestion`; print ``.render()`` for the reasoning.
    """
    hpc = dict(hpc_profile or {})
    cores_per_node = int(hpc.get("cores_per_node", 16))
    mem_gb_per_node = float(hpc.get("mem_gb_per_node", 64))
    max_array_size = int(hpc.get("max_array_size", 1000))
    walltime_cap_hours = float(hpc.get("walltime_cap_hours", 24))

    stats = _genome_stats(genome_fai)
    rationale = []

    # --- chunk count: balance sequence-per-chunk against scheduler limits ---
    by_size = max(1, -(-stats.genome_size // TARGET_CHUNK_BP))  # ceil
    # A scaffold is never split below itself; you cannot have fewer chunks than
    # scaffolds when every scaffold is its own chunk, but you *can* have more by
    # cutting large scaffolds. We recommend by sequence, then clamp to the array.
    target_chunks = min(by_size, max_array_size)
    if by_size > max_array_size:
        rationale.append(
            f"sequence target wanted {by_size} chunks but the array cap is "
            f"{max_array_size}; clamped — chunks will be larger / slower"
        )
    else:
        rationale.append(
            f"~{TARGET_CHUNK_BP // 1_000_000} Mb/chunk → {by_size} chunks "
            f"(under the {max_array_size}-task array cap)"
        )
    if stats.scaffold_count > target_chunks:
        rationale.append(
            f"genome is fragmented ({stats.scaffold_count} scaffolds > "
            f"{target_chunks} target chunks): each scaffold is at least one chunk, "
            "so expect more chunks than the target (and small tail chunks)"
        )

    bp_per_chunk = stats.genome_size / max(target_chunks, 1)
    mb_per_chunk = bp_per_chunk / 1_000_000

    # --- per-chunk resources ---
    threads = min(cores_per_node, 8)
    procs = max(1, min(cores_per_node // 2, 4))
    mem_gb = max(MIN_CHUNK_MEM_GB, int(mb_per_chunk * MEM_GB_PER_MB_SEQUENCE))
    mem_gb = min(mem_gb, int(mem_gb_per_node))
    walltime_min = max(MIN_WALLTIME_MIN, int(mb_per_chunk * WALLTIME_MIN_PER_MB))
    walltime_min = min(walltime_min, int(walltime_cap_hours * 60))

    rationale.append(
        "more chunks ⇒ more parallelism but more inter-chunk boundaries (each a "
        "potential lost cross-boundary Mikado merge) and more aggregate overhead"
    )
    rationale.append(
        f"min_boundary_gap={DEFAULT_MIN_BOUNDARY_GAP} keeps boundaries clear of "
        "typical merges; raise it if you see split merges, lower it for more cut sites"
    )
    rationale.append(
        f"resources are per chunk: serialise's SQLite DB must sit on fast local "
        f"scratch (never network FS); mem scales ~{mem_gb} GB at ~{mb_per_chunk:.0f} Mb/chunk"
    )

    suggestion = Suggestion(
        stats=stats,
        target_chunks=target_chunks,
        target_loci_per_chunk=None,
        min_boundary_gap=DEFAULT_MIN_BOUNDARY_GAP,
        procs=procs,
        threads=threads,
        mem_gb=mem_gb,
        walltime_min=walltime_min,
        rationale=rationale,
    )
    _log.info(
        "suggested %d chunks (%.0f Mb/chunk), %d GB / %s per chunk",
        target_chunks,
        mb_per_chunk,
        mem_gb,
        suggestion.walltime,
    )
    return suggestion
