"""Per-chunk config materialisation + dispatch."""

from __future__ import annotations

import dataclasses
import shlex
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from helixforge.reconcile.pipeline import PipelineConfig, run_pipeline
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.parallel.plan import Plan

_log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Per-chunk PipelineConfig
# ---------------------------------------------------------------------------


def _chunk_prefix(base_prefix: str, chunk_id: str) -> str:
    return f"{base_prefix}.{chunk_id}"


def build_chunk_configs(
    plan: Plan, base_config: PipelineConfig
) -> list[PipelineConfig]:
    """Clone ``base_config`` once per chunk with region + reserved ids + paths.

    Each chunk gets its own ``output_prefix`` (``<base>.<chunk_id>``), a per-chunk
    ``id_map_path`` (so chunks never write the same map), and the chunk's
    ``id_base``/``novel_base`` from :func:`reserve_id_ranges`. ``report_path`` and
    ``work_dir`` are reset to ``None`` so ``PipelineConfig.__post_init__``
    re-derives them from the new prefix. Cheap per-region evidence filtering
    (StringTie/SJ by scaffold) is left to the pipeline's ``region`` handling;
    BAMs are queried by region at run time.
    """
    configs = []
    for chunk in plan.chunks:
        if len(chunk.regions) != 1:
            raise NotImplementedError(
                f"chunk {chunk.chunk_id} has {len(chunk.regions)} regions; the "
                "single-region runner supports exactly one region per chunk"
            )
        prefix = _chunk_prefix(base_config.output_prefix, chunk.chunk_id)
        configs.append(
            dataclasses.replace(
                base_config,
                region=chunk.regions[0],
                chunk_id=chunk.chunk_id,
                id_base=chunk.id_base,
                novel_base=chunk.novel_base,
                output_prefix=prefix,
                id_map_path=f"{prefix}.id_map.json",
                report_path=None,
                work_dir=None,
            )
        )
    return configs


# ---------------------------------------------------------------------------
# PipelineConfig → `helixforge reconcile` argv
# ---------------------------------------------------------------------------

# (config_attr, flag, kind). kind ∈ {value, list, bool}; ``bool`` carries the
# (on, off) flag pair matching the click `--x/--no-x` options. ``None`` scalars
# and empty lists are skipped. Mirrors cli.__init__._PIPELINE_OPTIONS.
_ARG_SPEC: list[tuple[str, str | tuple[str, str], str]] = [
    ("genome_fasta", "--genome", "value"),
    ("helixer_gff3", "--helixer", "value"),
    ("helixer_h5", "--helixer-h5", "value"),
    ("stringtie_list", "--stringtie", "list"),
    ("bam_paths", "--bam", "list"),
    ("bigwig_paths", "--bigwig", "list"),
    ("star_sj_paths", "--star-sj", "list"),
    ("miniprot_gff", "--miniprot", "value"),
    ("protein_db", "--protein-db", "value"),
    ("transdecoder_bin_dir", "--transdecoder-bin-dir", "value"),
    (
        "backstop_transdecoder",
        ("--backstop-transdecoder", "--no-backstop-transdecoder"),
        "bool",
    ),
    ("portcullis_bin", "--portcullis-bin", "value"),
    ("mikado_bin", "--mikado-bin", "value"),
    ("diamond_bin", "--diamond-bin", "value"),
    ("use_mikado_configure", ("--mikado-configure", "--no-mikado-configure"), "bool"),
    ("scoring_profile", "--scoring-profile", "value"),
    ("helixer_is_reference", ("--helixer-reference", "--no-helixer-reference"), "bool"),
    ("helixer_support_weight", "--helixer-support-weight", "value"),
    ("output_prefix", "--output-prefix", "value"),
    ("report_path", "--report-path", "value"),
    ("id_map_path", "--id-map-path", "value"),
    ("work_dir", "--work-dir", "value"),
    ("region", "--region", "value"),
    ("chunk_id", "--chunk-id", "value"),
    ("id_base", "--id-base", "value"),
    ("novel_base", "--novel-base", "value"),
    ("procs", "--procs", "value"),
    ("threads", "--threads", "value"),
    ("min_tpm", "--min-tpm", "value"),
    ("min_samples", "--min-samples", "value"),
    ("coverage_threshold", "--coverage-threshold", "value"),
    ("near_zero_coverage", "--near-zero-coverage", "value"),
    ("as_report", ("--as-report", "--no-as-report"), "bool"),
    (
        "only_confirmed_introns",
        ("--only-confirmed-introns", "--allow-unconfirmed-introns"),
        "bool",
    ),
    ("max_isoforms", "--max-isoforms", "value"),
    (
        "keep_retained_introns",
        ("--keep-retained-introns", "--drop-retained-introns"),
        "bool",
    ),
    ("pad", ("--pad", "--no-pad"), "bool"),
    ("chimera_split", ("--chimera-split", "--no-chimera-split"), "bool"),
    ("flank", "--flank", "value"),
    ("reciprocal_overlap", "--reciprocal-overlap", "value"),
    ("min_cds_overlap", "--min-cds-overlap", "value"),
    ("min_cdna_overlap", "--min-cdna-overlap", "value"),
    ("admit_novel", ("--admit-novel", "--no-admit-novel"), "bool"),
    ("novel_evidence_floor", "--novel-evidence-floor", "value"),
    ("short_cds_threshold", "--short-cds-threshold", "value"),
    ("short_exon_threshold", "--short-exon-threshold", "value"),
    ("long_intron_threshold", "--long-intron-threshold", "value"),
    ("junction_tolerance", "--junction-tolerance", "value"),
    ("junction_min_reads", "--junction-min-reads", "value"),
]


def pipeline_config_to_reconcile_argv(
    config: PipelineConfig,
    helixforge_bin: str = "helixforge",
) -> list[str]:
    """Render a :class:`PipelineConfig` as a ``helixforge reconcile`` argv list."""
    argv: list[str] = [helixforge_bin, "reconcile"]
    for attr, flag, kind in _ARG_SPEC:
        value = getattr(config, attr)
        if kind == "value":
            if value is None:
                continue
            argv += [cast(str, flag), str(value)]
        elif kind == "list":
            for item in value or []:
                argv += [cast(str, flag), str(item)]
        elif kind == "bool":
            bool_flags = cast(tuple[str, str], flag)
            argv.append(bool_flags[0] if value else bool_flags[1])
    return argv


def _quote(argv: list[str]) -> str:
    return " ".join(shlex.quote(str(a)) for a in argv)


# Per-chunk varying flags: rendered from placeholders in the default template,
# so they are stripped from the base-config rendering before the suffix is added.
_PER_CHUNK_FLAGS = frozenset(
    {
        "--region",
        "--chunk-id",
        "--id-base",
        "--novel-base",
        "--output-prefix",
        "--report-path",
        "--id-map-path",
        "--work-dir",
    }
)


def default_reconcile_template(
    base_config: PipelineConfig,
    *,
    helixforge_bin: str = "helixforge",
) -> str:
    """Build the out-of-the-box ``--command`` template for ``parallel tasks``.

    Renders ``base_config`` (the genome / Helixer / evidence / scoring inputs the
    user attached) as a ``helixforge reconcile`` command, drops the flags that
    vary per chunk, and appends them back as placeholders so each chunk runs over
    its own region, reserved id range, and output prefix. The result is fully
    overridable: pass ``--command`` to replace it wholesale (v1 model).
    """
    argv = pipeline_config_to_reconcile_argv(base_config, helixforge_bin)
    kept: list[str] = []
    skip_value = False
    for tok in argv:
        if skip_value:
            skip_value = False
            continue
        if tok in _PER_CHUNK_FLAGS:
            skip_value = True
            continue
        kept.append(tok)
    base = _quote(kept)
    suffix = (
        " --region {region} --id-base {id_start} --novel-base {novel_start}"
        " --chunk-id {chunk_id} --output-prefix {output_dir}/{chunk_id}"
        " --id-map-path {output_dir}/{chunk_id}.id_map.json"
    )
    return base + suffix


# ---------------------------------------------------------------------------
# Slurm array + shell driver
# ---------------------------------------------------------------------------


def write_slurm_array(
    plan: Plan,
    base_config: PipelineConfig,
    script_path: str | Path,
    *,
    helixforge_bin: str = "helixforge",
    job_name: str = "helixforge",
    cpus_per_task: int | None = None,
    mem: str | None = None,
    time: str | None = None,
    partition: str | None = None,
    account: str | None = None,
    max_concurrent: int | None = None,
    log_dir: str = "logs",
    extra_sbatch: tuple[str, ...] | list[str] = (),
) -> Path:
    """Write a Slurm array script: one chunk per task (``--array=0-(N-1)``).

    Each task selects its chunk by ``SLURM_ARRAY_TASK_ID`` and runs the rendered
    ``helixforge reconcile`` command. Resource directives default to the chunk
    configs' ``procs`` (cpus-per-task) when not given explicitly.
    """
    configs = build_chunk_configs(plan, base_config)
    n = len(configs)
    if n == 0:
        raise ValueError("plan has no chunks")
    if cpus_per_task is None:
        cpus_per_task = max(c.procs for c in configs)

    array = f"0-{n - 1}"
    if max_concurrent:
        array += f"%{max_concurrent}"

    lines = ["#!/bin/bash"]
    sbatch = [
        ("job-name", job_name),
        ("array", array),
        ("cpus-per-task", cpus_per_task),
        ("mem", mem),
        ("time", time),
        ("partition", partition),
        ("account", account),
        ("output", f"{log_dir}/%x_%a.out"),
        ("error", f"{log_dir}/%x_%a.err"),
    ]
    for key, value in sbatch:
        if value is not None:
            lines.append(f"#SBATCH --{key}={value}")
    for raw in extra_sbatch:
        lines.append(f"#SBATCH {raw}")

    lines += ["", "set -euo pipefail", f"mkdir -p {shlex.quote(log_dir)}", ""]
    lines.append('case "$SLURM_ARRAY_TASK_ID" in')
    for i, config in enumerate(configs):
        cmd = _quote(pipeline_config_to_reconcile_argv(config, helixforge_bin))
        lines.append(f"  {i}) {cmd} ;;")
    lines += [
        '  *) echo "no such chunk: $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;',
        "esac",
        "",
    ]

    Path(script_path).write_text("\n".join(lines))
    _log.info("wrote Slurm array (%d tasks) → %s", n, script_path)
    return Path(script_path)


def write_shell_driver(
    plan: Plan,
    base_config: PipelineConfig,
    script_path: str | Path,
    *,
    helixforge_bin: str = "helixforge",
    workers: int = 1,
) -> Path:
    """Write a portable (scheduler-free) shell driver running every chunk.

    With ``workers == 1`` chunks run sequentially under ``set -e``; with
    ``workers > 1`` each chunk is backgrounded and a simple gate caps concurrency
    to ``workers`` (no scheduler, just ``wait -n``).
    """
    configs = build_chunk_configs(plan, base_config)
    lines = ["#!/bin/bash", "set -euo pipefail", ""]
    if workers and workers > 1:
        lines += [f"MAX_JOBS={int(workers)}", ""]
        for config in configs:
            cmd = _quote(pipeline_config_to_reconcile_argv(config, helixforge_bin))
            lines.append(
                'while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n; done'
            )
            lines.append(f"{cmd} &")
        lines.append("wait")
    else:
        for config in configs:
            cmd = _quote(pipeline_config_to_reconcile_argv(config, helixforge_bin))
            lines.append(cmd)
    lines.append("")

    Path(script_path).write_text("\n".join(lines))
    _log.info(
        "wrote shell driver (%d chunks, workers=%d) → %s",
        len(configs),
        workers,
        script_path,
    )
    return Path(script_path)


# ---------------------------------------------------------------------------
# Local process pool
# ---------------------------------------------------------------------------


def run_chunk(config: PipelineConfig) -> Path:
    """Run one chunk's pipeline; return its ``output_prefix`` (a ``Path``)."""
    run_pipeline(config)
    return Path(config.output_prefix)


def run_local(
    plan: Plan,
    base_config: PipelineConfig,
    workers: int = 4,
    *,
    _executor_cls: Any = ProcessPoolExecutor,
) -> list[Path]:
    """Run every chunk through a bounded pool; return per-chunk output prefixes.

    ``_executor_cls`` is a test seam — the unit suite injects a
    ``ThreadPoolExecutor`` so a mocked :func:`run_pipeline` is observed in-process
    (a real process pool would re-import the unpatched module). Production uses
    the default ``ProcessPoolExecutor``.
    """
    configs = build_chunk_configs(plan, base_config)
    prefixes: list[Path] = []
    with _executor_cls(max_workers=workers) as executor:
        futures: dict[Future[Path], PipelineConfig] = {
            executor.submit(run_chunk, c): c for c in configs
        }
        for future in as_completed(futures):
            chunk_cfg = futures[future]
            try:
                prefixes.append(future.result())
            except Exception:
                _log.exception("chunk %s failed", chunk_cfg.chunk_id)
                raise
    # Stable order (as_completed is arrival order) for reproducible callers.
    prefixes.sort(key=lambda p: str(p))
    _log.info("ran %d chunks locally (workers=%d)", len(configs), workers)
    return prefixes
