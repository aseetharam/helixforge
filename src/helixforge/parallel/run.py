"""The ``helixforge run`` whole-genome driver."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from helixforge.parallel.plan import Plan
    from helixforge.reconcile.pipeline import PipelineConfig

from helixforge.parallel.aggregate import aggregate
from helixforge.parallel.plan import (
    partition_genome,
    reserve_id_ranges,
    write_plan,
)
from helixforge.parallel.stitch import GeneSpan, find_boundary_merge_candidates
from helixforge.parallel.suggest import suggest_plan
from helixforge.parallel.tasks import (
    run_local,
    write_slurm_array,
)
from helixforge.reconcile.pipeline import run_pipeline
from helixforge.utils.atomic import atomic_write
from helixforge.utils.logging import get_logger
from helixforge.utils.regions import gff3_to_internal

_log = get_logger(__name__)


@dataclass
class RunResult:
    """Outcome of :func:`run_genome`, what ran, where it landed."""

    mode: str  # "single" | "scatter-local" | "scatter-hypershell" | "scatter-slurm"
    manifest_path: Path
    out_prefix: str
    genes: list[Any] | None = None  # single mode, list[ReconciledGene] at runtime
    aggregate: Any | None = None  # AggregateResult (scatter-local)
    plan: Any | None = None  # Plan
    chunk_prefixes: list[str] = field(default_factory=list)
    script_path: Path | None = None  # scatter-slurm
    boundary_merges_recovered: int = 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_target_chunks(
    scatter: Any,
    hpc_profile: Any,
) -> tuple[int | None, str | None]:
    """Map ``scatter`` (``"auto"`` | int) to a ``target_chunks`` (None ⇒ auto)."""
    if isinstance(scatter, int):
        if scatter < 1:
            raise ValueError(f"scatter chunk count must be >= 1, got {scatter}")
        return scatter, None
    if scatter == "auto":
        return None, "auto"
    raise ValueError(f"scatter must be 'off', 'auto', or an int, got {scatter!r}")


def _read_gene_spans(gff3_path: str) -> list[GeneSpan]:
    """Lightweight ``[GeneSpan]`` from an aggregated HelixForge GFF3 (gene lines).

    Used by the optional boundary-stitch detection (D4): it needs only each gene's
    span + strand + id, so a full reconciled-gene reload is unnecessary.
    """
    spans: list[GeneSpan] = []
    with open(gff3_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            start, end = gff3_to_internal(int(cols[3]), int(cols[4]))
            gene_id = None
            for part in cols[8].split(";"):
                part = part.strip()
                if part.startswith("ID="):
                    gene_id = part[3:]
                    break
            if gene_id is not None:
                spans.append(GeneSpan(gene_id, cols[0], start, end, cols[6]))
    return spans


def _stitch_count(out_prefix: str, plan: Plan, base_config: PipelineConfig) -> int:
    """Count boundary-recoverable cross-chunk merges off the aggregated GFF3 (D4).

    Reloads gene spans + the genome-wide STAR-SJ junctions and reuses the
    reconciler's bridging predicate. Reports the count (the structural merge is
    applied by :func:`~helixforge.parallel.stitch.boundary_stitch` on in-memory
    gene sets). Best-effort: never sinks a run.
    """
    from helixforge.reconcile.pipeline import _merge_star_junctions

    try:
        spans = _read_gene_spans(f"{out_prefix}.gff3")
        junctions = (
            _merge_star_junctions(
                base_config.star_sj_paths, base_config.junction_min_reads
            )
            if base_config.star_sj_paths
            else []
        )
        pairs = find_boundary_merge_candidates(spans, plan, junctions)
        if pairs:
            _log.info(
                "boundary-stitch: %d cross-chunk merge(s) recoverable: %s",
                len(pairs),
                pairs,
            )
        else:
            _log.info("boundary-stitch: no cross-chunk merges to recover")
        return len(pairs)
    except Exception as exc:  # noqa: BLE001 - accuracy pass never sinks a run
        _log.warning("boundary-stitch detection failed: %s", exc)
        return 0


def _write_manifest(path: str | Path, payload: dict[str, Any]) -> Path:
    with atomic_write(str(path)) as fh:
        fh.write(json.dumps(payload, indent=2, sort_keys=False, default=str))
    _log.info("wrote run manifest → %s", path)
    return Path(path)


def _dump_resolved_config(base_config: PipelineConfig, out_prefix: str) -> str | None:
    """Best-effort resolved-config YAML next to the manifest (D2 reuse)."""
    try:
        return base_config.to_yaml(f"{out_prefix}.run.resolved.yaml")
    except Exception as exc:  # noqa: BLE001
        _log.warning("could not write resolved config: %s", exc)
        return None


# ---------------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------------


def run_genome(
    base_config: PipelineConfig,
    *,
    scatter: str | int = "auto",
    hpc: str = "local",
    workers: int = 4,
    target_loci_per_chunk: int | None = None,
    min_boundary_gap: int | None = None,
    hpc_profile: dict[str, Any] | None = None,
    out_prefix: str | None = None,
    master_id_map_path: str | None = None,
    stitch: bool = False,
    script_path: str | None = None,
    hs_bin: str = "hs",
    num_tasks: int | None = None,
) -> RunResult:
    """Drive a whole genome end to end → :class:`RunResult`.

    ``scatter``: ``"off"`` runs the single-process pipeline unchanged; ``"auto"``
    sizes the partition from the genome index (``suggest``); an int targets that
    many chunks. ``hpc``: ``"local"`` runs the chunks now over a bounded
    ``ProcessPoolExecutor`` (``workers``) then aggregates; ``"hypershell"`` runs
    them now via ``hs cluster`` (``num_tasks`` concurrency, defaults to
    ``workers``) then aggregates; ``"slurm"`` emits a Slurm array (one chunk/task)
    and stops, the user submits it and re-invokes ``aggregate`` afterward.
    ``stitch`` (local/hypershell) reports boundary-recoverable
    cross-chunk merges (D4). Honors per-chunk + global checkpoint/resume via
    ``base_config.resume`` (carried into every chunk config).
    """
    out_prefix = out_prefix or base_config.output_prefix
    master_id_map_path = master_id_map_path or base_config.id_map_path

    # --- single-process path: byte-identical to run_pipeline ---
    if scatter == "off":
        _log.info("[run] scatter=off, single-process pipeline")
        genes = run_pipeline(base_config)
        resolved = _dump_resolved_config(base_config, out_prefix)
        manifest = _write_manifest(
            f"{out_prefix}.manifest.json",
            {
                "mode": "single",
                "scatter": "off",
                "out_prefix": out_prefix,
                "num_genes": len(genes),
                "id_map_path": base_config.id_map_path,
                "resolved_config": resolved,
            },
        )
        return RunResult(
            mode="single", manifest_path=manifest, out_prefix=out_prefix, genes=genes
        )

    # --- scatter path: plan → tasks → execute → aggregate ---
    target_chunks, auto = _resolve_target_chunks(scatter, hpc_profile)
    if auto == "auto":
        suggestion = suggest_plan(base_config.genome_fasta, hpc_profile=hpc_profile)
        target_chunks = suggestion.target_chunks
        if min_boundary_gap is None:
            min_boundary_gap = suggestion.min_boundary_gap
        _log.info("[run] scatter=auto → %d target chunks", target_chunks)

    master_id_map = None
    if master_id_map_path and Path(master_id_map_path).exists():
        master_id_map = json.loads(Path(master_id_map_path).read_text())

    plan = partition_genome(
        base_config.genome_fasta,
        base_config.helixer_gff3,
        target_chunks=target_chunks,
        target_loci_per_chunk=target_loci_per_chunk,
        min_boundary_gap=min_boundary_gap,
        flank=base_config.flank,
    )
    reserve_id_ranges(plan, id_map=master_id_map)
    plan_path = write_plan(plan, f"{out_prefix}.plan.json")
    resolved = _dump_resolved_config(base_config, out_prefix)

    # --- emit a Slurm array and stop (chunks run under the scheduler) ---
    if hpc == "slurm":
        script = script_path or f"{out_prefix}.array.sbatch"
        write_slurm_array(plan, base_config, script)
        manifest = _write_manifest(
            f"{out_prefix}.manifest.json",
            {
                "mode": "scatter-slurm",
                "scatter": scatter,
                "hpc": "slurm",
                "out_prefix": out_prefix,
                "num_chunks": len(plan),
                "total_loci": plan.total_loci,
                "plan_path": str(plan_path),
                "slurm_array": str(script),
                "master_id_map_path": master_id_map_path,
                "resolved_config": resolved,
                "next": (
                    f"submit {script}; then `helixforge parallel aggregate` "
                    f"over the per-chunk {out_prefix}.chunk_* prefixes"
                ),
            },
        )
        return RunResult(
            mode="scatter-slurm",
            manifest_path=manifest,
            out_prefix=out_prefix,
            plan=plan,
            script_path=Path(script),
        )

    if hpc not in ("local", "hypershell"):
        raise ValueError(
            f"hpc must be 'local', 'slurm', or 'hypershell', got {hpc!r}"
        )

    # --- run every chunk (local pool or HyperShell), then gather ---
    if hpc == "hypershell":
        from helixforge.parallel.hypershell_backend import run_hypershell

        prefixes = run_hypershell(
            plan, base_config, out_prefix,
            num_tasks=num_tasks or workers, hs_bin=hs_bin,
        )
    else:
        prefixes = run_local(plan, base_config, workers=workers)

    result = aggregate(
        [str(p) for p in prefixes], out_prefix, master_id_map_path=master_id_map_path
    )

    recovered = _stitch_count(out_prefix, plan, base_config) if stitch else 0

    manifest = _write_manifest(
        f"{out_prefix}.manifest.json",
        {
            "mode": "scatter-local" if hpc == "local" else "scatter-hypershell",
            "scatter": scatter,
            "hpc": hpc,
            "workers": workers,
            "out_prefix": out_prefix,
            "num_chunks": len(plan),
            "num_genes": result.num_genes,
            "num_loci": result.num_loci,
            "tier_counts": result.tier_counts,
            "origin_counts": result.origin_counts,
            "plan_path": str(plan_path),
            "gff3_path": str(result.gff3_path),
            "report_path": str(result.report_path),
            "master_id_map_path": master_id_map_path,
            "chunk_prefixes": [str(p) for p in prefixes],
            "boundary_merges_recovered": recovered,
            "stitch": stitch,
            "resolved_config": resolved,
        },
    )
    return RunResult(
        mode="scatter-local" if hpc == "local" else "scatter-hypershell",
        manifest_path=manifest,
        out_prefix=out_prefix,
        aggregate=result,
        plan=plan,
        chunk_prefixes=[str(p) for p in prefixes],
        boundary_merges_recovered=recovered,
    )
