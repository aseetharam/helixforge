"""Minimal SLURM utilities.

For users who prefer direct SLURM integration over HyperShell.
This module provides helpers for detecting SLURM environments
and accessing allocated resources.

Recommended approach: Use HyperShell with task files from
parallel/taskgen.py for scheduler-agnostic parallelization.

Example:
    >>> from helixforge.parallel.slurm import detect_slurm_environment
    >>> env = detect_slurm_environment()
    >>> if env:
    ...     print(f"Running in SLURM job {env['SLURM_JOB_ID']}")
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)


# =============================================================================
# Environment Detection
# =============================================================================


def detect_slurm_environment() -> dict | None:
    """Detect if running under SLURM.

    Returns:
        Dict with SLURM env vars if in SLURM, None otherwise.

    Example:
        >>> env = detect_slurm_environment()
        >>> if env:
        ...     print(f"Job ID: {env['SLURM_JOB_ID']}")
    """
    if "SLURM_JOB_ID" not in os.environ:
        return None

    env_vars = [
        "SLURM_JOB_ID",
        "SLURM_ARRAY_JOB_ID",
        "SLURM_ARRAY_TASK_ID",
        "SLURM_ARRAY_TASK_COUNT",
        "SLURM_CPUS_PER_TASK",
        "SLURM_MEM_PER_NODE",
        "SLURM_JOB_NAME",
        "SLURM_SUBMIT_DIR",
        "SLURM_NODELIST",
    ]
    return {k: os.environ[k] for k in env_vars if k in os.environ}


def is_slurm_job() -> bool:
    """Check if running under SLURM.

    Returns:
        True if SLURM environment detected.
    """
    return "SLURM_JOB_ID" in os.environ


def get_slurm_task_id() -> int | None:
    """Get current SLURM array task ID.

    Returns:
        Task ID (0-indexed) or None if not an array job.
    """
    task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
    if task_id is not None:
        return int(task_id)
    return None


# =============================================================================
# Resource Access
# =============================================================================


def get_slurm_resources() -> dict:
    """Get allocated resources from SLURM environment.

    Returns:
        Dict with:
            - cpus: Number of CPUs per task
            - memory_mb: Memory in MB (if available)
            - task_id: Array task ID (if array job)
            - node: Node name (if available)

    Example:
        >>> res = get_slurm_resources()
        >>> print(f"Using {res['cpus']} CPUs")
    """
    env = detect_slurm_environment()
    if not env:
        return {
            "cpus": 1,
            "memory_mb": None,
            "task_id": None,
            "node": None,
        }

    return {
        "cpus": int(env.get("SLURM_CPUS_PER_TASK", 1)),
        "memory_mb": _parse_slurm_memory(env.get("SLURM_MEM_PER_NODE")),
        "task_id": (
            int(env["SLURM_ARRAY_TASK_ID"])
            if "SLURM_ARRAY_TASK_ID" in env
            else None
        ),
        "node": env.get("SLURM_NODELIST"),
    }


def _parse_slurm_memory(mem_str: str | None) -> int | None:
    """Parse SLURM memory string to MB.

    Args:
        mem_str: Memory string (e.g., '16G', '16000M', '16000').

    Returns:
        Memory in MB, or None if parsing fails.
    """
    if not mem_str:
        return None

    mem_str = mem_str.strip().upper()

    try:
        if mem_str.endswith("G"):
            return int(float(mem_str[:-1]) * 1024)
        elif mem_str.endswith("M"):
            return int(mem_str[:-1])
        elif mem_str.endswith("K"):
            return int(float(mem_str[:-1]) / 1024)
        else:
            # Assume MB
            return int(mem_str)
    except ValueError:
        return None


# =============================================================================
# Example SBATCH Templates
# =============================================================================


EXAMPLE_SBATCH_HYPERSHELL = """\
#!/bin/bash
#SBATCH --job-name=helixforge
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

# Example SLURM script for running HyperShell
# Modify partition, time, and resources for your cluster

# Load required modules (site-specific)
# module load python/3.10
# module load conda
# conda activate helixforge

# Ensure output directories exist
mkdir -p outputs logs

# Run HyperShell with allocated CPUs
hs launch --parallelism ${SLURM_CPUS_PER_TASK} < tasks.txt
"""

EXAMPLE_SBATCH_PARALLEL = """\
#!/bin/bash
#SBATCH --job-name=helixforge
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

# Example SLURM script for running GNU Parallel
# Modify partition, time, and resources for your cluster

# Load required modules (site-specific)
# module load python/3.10
# module load parallel
# module load conda
# conda activate helixforge

# Ensure output directories exist
mkdir -p outputs logs

# Run GNU Parallel with allocated CPUs
parallel -j ${SLURM_CPUS_PER_TASK} < tasks.txt
"""


def write_example_sbatch(
    output_path: Path | str,
    executor: str = "hypershell",
) -> Path:
    """Write example SBATCH script for reference.

    Args:
        output_path: Where to write the script.
        executor: Which executor to use ("hypershell" or "parallel").

    Returns:
        Path to the written script.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if executor == "hypershell":
        content = EXAMPLE_SBATCH_HYPERSHELL
    elif executor == "parallel":
        content = EXAMPLE_SBATCH_PARALLEL
    else:
        raise ValueError(f"Unknown executor: {executor}")

    output_path.write_text(content)
    logger.info(f"Wrote example SBATCH script to {output_path}")

    return output_path


# =============================================================================
# Chunk Access (for use within SLURM array jobs)
# =============================================================================


def get_chunk_for_task(
    chunk_plan_path: Path | str,
    task_id: int | None = None,
) -> "GenomicChunk":
    """Get chunk assignment for current SLURM array task.

    Args:
        chunk_plan_path: Path to chunk plan JSON.
        task_id: Task ID (uses SLURM_ARRAY_TASK_ID if None).

    Returns:
        GenomicChunk for this task.

    Raises:
        RuntimeError: If task ID not available.
        IndexError: If task ID exceeds chunk count.
    """
    from helixforge.parallel.chunker import ChunkPlan

    if task_id is None:
        task_id = get_slurm_task_id()
        if task_id is None:
            raise RuntimeError(
                "Task ID not provided and SLURM_ARRAY_TASK_ID not set"
            )

    plan = ChunkPlan.load(chunk_plan_path)

    if task_id >= len(plan.chunks):
        raise IndexError(
            f"Task ID {task_id} exceeds chunk count {len(plan.chunks)}"
        )

    return plan.chunks[task_id]


# =============================================================================
# Exports
# =============================================================================


__all__ = [
    "detect_slurm_environment",
    "is_slurm_job",
    "get_slurm_task_id",
    "get_slurm_resources",
    "write_example_sbatch",
    "get_chunk_for_task",
    "EXAMPLE_SBATCH_HYPERSHELL",
    "EXAMPLE_SBATCH_PARALLEL",
]
