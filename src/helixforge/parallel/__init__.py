"""Parallelization utilities for HelixForge.

This module provides tools for parallel execution of HelixForge operations:

- Genome chunking for parallel processing
- Local multiprocessing execution
- Task file generation for HyperShell/GNU Parallel
- Memory monitoring

For HPC clusters, the recommended approach is:
1. Create a chunk plan with GenomeChunker
2. Generate a task file with TaskGenerator
3. Execute with HyperShell or GNU Parallel

Example:
    >>> from helixforge.parallel import GenomeChunker, ChunkStrategy, TaskGenerator
    >>> chunker = GenomeChunker(genome)
    >>> plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)
    >>> gen = TaskGenerator(plan)
    >>> task_file = gen.generate(
    ...     command_template="helixforge confidence --chunk-id {chunk_id} -o out/{chunk_id}.tsv",
    ...     output_path="tasks.txt",
    ... )
    >>> # Execute with: hs launch --parallelism 32 < tasks.txt
"""

from helixforge.parallel.chunker import (
    ChunkPlan,
    ChunkStrategy,
    GenomicChunk,
    GenomeChunker,
    Chunker,  # Legacy alias
    merge_overlapping_results,
    suggest_chunk_parameters,
)

from helixforge.parallel.executor import (
    ExecutorBackend,
    ExecutionStats,
    MemoryMonitor,
    MemoryStats,
    ParallelExecutor,
    TaskResult,
    ChunkProcessor,
    Executor,  # Legacy alias
    get_memory_stats,
    get_optimal_workers,
    track_memory,
)

from helixforge.parallel.taskgen import (
    TaskFile,
    TaskGenerator,
    format_command,
    generate_simple_task_file,
    estimate_parallelism,
    generate_hypershell_command,
)

from helixforge.parallel.slurm import (
    detect_slurm_environment,
    is_slurm_job,
    get_slurm_task_id,
    get_slurm_resources,
    get_chunk_for_task,
    write_example_sbatch,
)

__all__ = [
    # Chunking
    "ChunkPlan",
    "ChunkStrategy",
    "GenomicChunk",
    "GenomeChunker",
    "Chunker",
    "merge_overlapping_results",
    "suggest_chunk_parameters",
    # Execution
    "ExecutorBackend",
    "ExecutionStats",
    "MemoryMonitor",
    "MemoryStats",
    "ParallelExecutor",
    "TaskResult",
    "ChunkProcessor",
    "Executor",
    "get_memory_stats",
    "get_optimal_workers",
    "track_memory",
    # Task Generation (HyperShell/GNU Parallel)
    "TaskFile",
    "TaskGenerator",
    "format_command",
    "generate_simple_task_file",
    "estimate_parallelism",
    "generate_hypershell_command",
    # SLURM Utilities (minimal)
    "detect_slurm_environment",
    "is_slurm_job",
    "get_slurm_task_id",
    "get_slurm_resources",
    "get_chunk_for_task",
    "write_example_sbatch",
]
