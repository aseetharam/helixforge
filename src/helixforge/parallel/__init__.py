"""Large-genome scatter-gather.

Partition a genome into gene-safe region chunks (a gene is never split across a
boundary), reserve each chunk a disjoint HFG id range so the gathered annotation
has globally-unique, rerun-stable ids, run the chunks (local pool or a Slurm
array), then aggregate and optionally boundary-stitch back into one annotation.

Two entry surfaces onto the same primitives:

- ``reconcile --scatter`` (one command, single host or Slurm) → :func:`run_genome`.
- the ``parallel`` command group (scheduler-agnostic four steps):
  ``suggest`` (:func:`suggest_plan`) → ``plan`` (:func:`partition_by_strategy` /
  :func:`partition_genome` + :func:`reserve_id_ranges`) → ``tasks``
  (:func:`generate_task_file` / :func:`default_reconcile_template`) → run with any
  executor (GNU parallel / xargs / Slurm array / HyperShell — see
  :func:`generate_hypershell_command`) → ``aggregate`` (:func:`aggregate`).

This module re-exports the public surface of its submodules for convenience; the
submodules remain the canonical import path.
"""

from helixforge.constants import (
    PLAN_DEFAULT_MIN_BOUNDARY_GAP as DEFAULT_MIN_BOUNDARY_GAP,
)
from helixforge.parallel.aggregate import (
    AggregateResult,
    aggregate,
    prefixes_from_pattern,
)
from helixforge.parallel.chunker import (
    DEFAULT_GENES_PER_CHUNK,
    DEFAULT_MIN_CHUNK_BP,
    DEFAULT_SIZE_CHUNK_BP,
    ChunkSpec,
    ChunkStrategy,
    plan_chunks,
)
from helixforge.parallel.plan import (
    Chunk,
    Plan,
    partition_by_strategy,
    partition_genome,
    read_plan,
    reserve_id_ranges,
    write_plan,
)
from helixforge.parallel.hypershell_backend import (
    run_hypershell,
    write_hypershell_plan,
)
from helixforge.parallel.run import RunResult, run_genome
from helixforge.parallel.stitch import (
    GeneSpan,
    boundary_stitch,
    find_boundary_merge_candidates,
)
from helixforge.parallel.suggest import GenomeStats, Suggestion, suggest_plan
from helixforge.parallel.taskgen import (
    TaskFile,
    estimate_parallelism,
    format_command,
    generate_hypershell_command,
    generate_task_file,
    get_optimal_workers,
    write_example_sbatch,
)
from helixforge.parallel.tasks import (
    build_chunk_configs,
    default_reconcile_template,
    pipeline_config_to_reconcile_argv,
    run_chunk,
    run_local,
    write_shell_driver,
    write_slurm_array,
)

__all__ = [
    # Planning / partitioning
    "Chunk",
    "Plan",
    "partition_genome",
    "partition_by_strategy",
    "reserve_id_ranges",
    "write_plan",
    "read_plan",
    "DEFAULT_MIN_BOUNDARY_GAP",
    # Chunking strategies
    "ChunkStrategy",
    "ChunkSpec",
    "plan_chunks",
    "DEFAULT_SIZE_CHUNK_BP",
    "DEFAULT_GENES_PER_CHUNK",
    "DEFAULT_MIN_CHUNK_BP",
    # Per-chunk config + dispatch
    "build_chunk_configs",
    "pipeline_config_to_reconcile_argv",
    "default_reconcile_template",
    "write_slurm_array",
    "write_shell_driver",
    "run_chunk",
    "run_local",
    # Executor-agnostic task files
    "TaskFile",
    "format_command",
    "generate_task_file",
    "generate_hypershell_command",
    "estimate_parallelism",
    "get_optimal_workers",
    "write_example_sbatch",
    # Whole-genome driver
    "RunResult",
    "run_genome",
    # HyperShell execution backend
    "run_hypershell",
    "write_hypershell_plan",
    # Aggregation + boundary stitch
    "AggregateResult",
    "aggregate",
    "prefixes_from_pattern",
    "GeneSpan",
    "find_boundary_merge_candidates",
    "boundary_stitch",
    # Suggestion
    "GenomeStats",
    "Suggestion",
    "suggest_plan",
]
