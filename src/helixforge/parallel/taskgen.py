"""Task-file generation from a user command template."""

from __future__ import annotations

import re
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from helixforge.parallel.plan import _region_bounds
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.parallel.plan import Chunk, Plan

_log = get_logger(__name__)

_PLACEHOLDER_RE = re.compile(r"\{(\w+)\}")


@dataclass
class TaskFile:
    """A generated task file: its path, the task count, and the template used."""

    path: Path
    n_tasks: int
    command_template: str
    wrapper_path: Path | None = None

    def preview(self, n: int = 3) -> list[str]:
        if not self.path.exists():
            return []
        out: list[str] = []
        with open(self.path) as fh:
            for i, line in enumerate(fh):
                if i >= n:
                    break
                out.append(line.rstrip("\n"))
        return out


def _chunk_variables(chunk: Chunk, output_dir: Path | None) -> dict[str, str]:
    """Resolve the placeholder values for one chunk (first region for coords)."""
    region0 = chunk.regions[0] if chunk.regions else chunk.chunk_id
    seqid, lo, hi = _region_bounds(region0)
    variables = {
        "chunk_id": chunk.chunk_id,
        "region": ",".join(chunk.regions),
        "seqid": seqid,
        "id_start": str(chunk.id_base),
        "novel_start": str(chunk.novel_base),
    }
    if lo is not None and hi is not None:
        variables.update(
            {
                "start": str(lo + 1),  # 1-based inclusive for readable CLI regions
                "end": str(hi),
                "start_0": str(lo),
                "end_0": str(hi),
                "size": str(hi - lo),
            }
        )
    if output_dir is not None:
        variables["output_dir"] = str(output_dir)
    return variables


def format_command(
    template: str,
    chunk: Chunk,
    *,
    output_dir: Path | None = None,
    extra_vars: dict[str, str] | None = None,
) -> str:
    """Expand ``template`` for one chunk; raise on an unknown placeholder."""
    variables = _chunk_variables(chunk, output_dir)
    if extra_vars:
        variables.update(extra_vars)
    unknown = set(_PLACEHOLDER_RE.findall(template)) - set(variables)
    if unknown:
        raise KeyError(
            f"unknown placeholder(s) in --command template: {sorted(unknown)}; "
            f"available: {sorted(variables)}"
        )
    return template.format(**variables)


def generate_task_file(
    plan: Plan,
    command_template: str,
    output_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    include_logging: bool = False,
    log_dir: str | Path | None = None,
    wrapper: str | Path | None = None,
    wrapper_setup: tuple[str, ...] | list[str] = (),
) -> TaskFile:
    """Expand ``command_template`` over every chunk → an executor-agnostic task file.

    One line per chunk. With ``include_logging`` each task redirects stdout/stderr
    to ``<log_dir>/<chunk_id>.log`` (``log_dir`` defaults to
    ``<output_dir>/logs`` or ``./logs``). With ``wrapper`` set, a reusable wrapper
    script is written (carrying ``wrapper_setup`` lines — module loads / conda
    activate) and each task line invokes ``bash <wrapper> '<expanded command>'``,
    so the setup runs once per task without bloating every line.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_dir = Path(output_dir) if output_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    log_dir_path: Path | None = None
    if include_logging:
        log_dir_path = (
            Path(log_dir)
            if log_dir
            else ((out_dir / "logs") if out_dir else Path("logs"))
        )
        log_dir_path.mkdir(parents=True, exist_ok=True)

    wrapper_path: Path | None = None
    if wrapper is not None:
        wrapper_path = _write_wrapper_script(wrapper, wrapper_setup)

    tasks: list[str] = []
    for chunk in plan.chunks:
        cmd = format_command(command_template, chunk, output_dir=out_dir)
        if wrapper_path is not None:
            cmd = f"bash {_sh_quote(str(wrapper_path))} {_sh_quote(cmd)}"
        if log_dir_path is not None:
            cmd = f"{cmd} > {log_dir_path / (chunk.chunk_id + '.log')} 2>&1"
        tasks.append(cmd)

    output_path.write_text("\n".join(tasks) + ("\n" if tasks else ""))
    _log.info(
        "wrote %d tasks → %s%s",
        len(tasks),
        output_path,
        f" (wrapper {wrapper_path})" if wrapper_path else "",
    )
    return TaskFile(
        path=output_path,
        n_tasks=len(tasks),
        command_template=command_template,
        wrapper_path=wrapper_path,
    )


def _sh_quote(s: str) -> str:
    import shlex

    return shlex.quote(s)


def _write_wrapper_script(
    wrapper: str | Path,
    setup: tuple[str, ...] | list[str],
) -> Path:
    """Write a reusable per-task wrapper that runs ``setup`` then ``"$@"``."""
    wrapper = Path(wrapper)
    wrapper.parent.mkdir(parents=True, exist_ok=True)
    setup_section = "\n".join(setup) if setup else "# (no setup commands)"
    script = f"""#!/bin/bash
set -euo pipefail

# HelixForge per-chunk wrapper — generated by `helixforge parallel tasks`.
# Runs the setup once (module loads / conda activate), then the chunk command
# passed as a single argument.

{setup_section}

eval "$1"
"""
    wrapper.write_text(script)
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    _log.info("wrote wrapper script → %s", wrapper)
    return wrapper


def estimate_parallelism(
    n_chunks: int,
    available_cores: int,
    overhead_factor: float = 1.5,
) -> int:
    """Recommend a worker count for a task file of ``n_chunks`` chunks.

    Never exceeds the chunk count or the core count; when chunks greatly
    outnumber cores it aims for a few chunks per worker so a slow chunk does not
    starve the pool. ``overhead_factor`` is accepted for call-site compatibility
    and does not change the result.
    """
    parallelism = min(n_chunks, available_cores)
    if n_chunks > available_cores * 4:
        parallelism = available_cores
    elif n_chunks > available_cores * 2:
        parallelism = max(1, available_cores // 2)
    return max(1, parallelism)


def generate_hypershell_command(
    task_file: str | Path,
    parallelism: int | None = None,
    timeout: int | None = None,
) -> str:
    """Render the HyperShell command that runs an executor-agnostic task file.

    The task file emitted by :func:`generate_task_file` is one shell command per
    line, which ``hs cluster`` consumes directly.
    """
    cmd = f"hs cluster {task_file}"
    if parallelism:
        cmd += f" --num-tasks {parallelism}"
    if timeout:
        cmd += f" --task-timeout {timeout}"
    return cmd


def get_optimal_workers(max_workers: int | None = None) -> int:
    """Pick a sane local worker count: available cores, capped by ``max_workers``."""
    import os

    cores = os.cpu_count() or 1
    return max(1, min(max_workers, cores) if max_workers else cores)


EXAMPLE_SBATCH_HYPERSHELL = """\
#!/bin/bash
#SBATCH --job-name=helixforge
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

# Example SLURM script for running a `helixforge parallel tasks` file with
# HyperShell. Modify partition, time, and resources for your cluster.

# Load required modules (site-specific)
# module load python/3.10
# module load conda
# conda activate helixforge

mkdir -p chunks logs
hs cluster tasks.txt --num-tasks ${SLURM_CPUS_PER_TASK}
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

# Example SLURM script for running a `helixforge parallel tasks` file with
# GNU Parallel. Modify partition, time, and resources for your cluster.

# Load required modules (site-specific)
# module load python/3.10
# module load parallel
# module load conda
# conda activate helixforge

mkdir -p chunks logs
parallel -j ${SLURM_CPUS_PER_TASK} < tasks.txt
"""


def write_example_sbatch(
    output_path: str | Path,
    executor: str = "hypershell",
) -> Path:
    """Write a copy-paste SBATCH wrapper that runs a task file under one node.

    This is the single-node "run the whole ``tasks.txt`` inside one allocation"
    pattern (HyperShell or GNU Parallel across ``$SLURM_CPUS_PER_TASK`` cores).
    For a true one-task-per-chunk Slurm *array* driven off a plan, use
    ``reconcile --hpc slurm`` / :func:`helixforge.parallel.tasks.write_slurm_array`.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if executor == "hypershell":
        content = EXAMPLE_SBATCH_HYPERSHELL
    elif executor == "parallel":
        content = EXAMPLE_SBATCH_PARALLEL
    else:
        raise ValueError(
            f"unknown executor {executor!r}; expected 'hypershell' or 'parallel'"
        )
    output_path.write_text(content)
    _log.info("wrote example sbatch (%s) → %s", executor, output_path)
    return output_path
