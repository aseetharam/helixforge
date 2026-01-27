"""Task file generation for embarrassingly parallel execution.

This module generates task files for parallel execution tools like:
- HyperShell (recommended): https://hypershell.readthedocs.io/
- GNU Parallel
- xargs
- Any tool that consumes one-command-per-line task files

HyperShell is the recommended approach as it provides scheduler-agnostic
parallelization that works with SLURM, PBS, SGE, or local execution.

Example:
    >>> from helixforge.parallel import ChunkPlan, TaskGenerator
    >>> plan = ChunkPlan.load("chunks.json")
    >>> gen = TaskGenerator(plan)
    >>> task_file = gen.generate(
    ...     command_template="helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o out/{chunk_id}.tsv",
    ...     output_path="tasks.txt",
    ... )
    >>> print(f"Generated {task_file.n_tasks} tasks")

    # Execute with HyperShell:
    # hs cluster tasks.txt --num-tasks 32
"""

from __future__ import annotations

import json
import logging
import re
import stat
from pathlib import Path
from typing import TYPE_CHECKING

import attrs

if TYPE_CHECKING:
    from helixforge.parallel.chunker import ChunkPlan, GenomicChunk

logger = logging.getLogger(__name__)


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class TaskFile:
    """Generated task file for parallel execution.

    Attributes:
        path: Path to the generated task file.
        n_tasks: Number of tasks in the file.
        chunk_plan_path: Path to the associated chunk plan JSON.
        command_template: Template used to generate commands.
    """

    path: Path
    n_tasks: int
    chunk_plan_path: Path | None
    command_template: str

    def preview(self, n: int = 5) -> list[str]:
        """Show first n tasks.

        Args:
            n: Number of tasks to show.

        Returns:
            List of first n task commands.
        """
        if not self.path.exists():
            return []

        lines = []
        with open(self.path) as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                lines.append(line.rstrip())
        return lines

    def read_all(self) -> list[str]:
        """Read all tasks.

        Returns:
            List of all task commands.
        """
        if not self.path.exists():
            return []
        return self.path.read_text().strip().split("\n")


# =============================================================================
# Command Formatting
# =============================================================================


def format_command(
    template: str,
    chunk: GenomicChunk,
    chunk_plan_path: Path | None = None,
    output_dir: Path | None = None,
    extra_vars: dict[str, str] | None = None,
) -> str:
    """Format command template with chunk variables.

    Available placeholders:
        {chunk_id}    - Unique chunk identifier
        {seqid}       - Scaffold/chromosome name
        {start}       - Start coordinate (1-based, for CLI commands)
        {end}         - End coordinate (1-based inclusive, for CLI commands)
        {start_0}     - Start coordinate (0-based, for internal use)
        {end_0}       - End coordinate (0-based exclusive, for internal use)
        {size}        - Chunk size in bases
        {chunk_plan}  - Path to chunk plan JSON
        {output_dir}  - Output directory path
        Any keys in extra_vars

    Note:
        {start} and {end} are converted from internal 0-based half-open
        coordinates to 1-based inclusive coordinates suitable for CLI
        commands. Use {start_0} and {end_0} for 0-based coordinates.

    Args:
        template: Command template with placeholders.
        chunk: GenomicChunk to format.
        chunk_plan_path: Optional path to chunk plan.
        output_dir: Optional output directory.
        extra_vars: Optional additional variables.

    Returns:
        Formatted command string.

    Raises:
        KeyError: If template contains unknown placeholder.

    Example:
        >>> format_command(
        ...     "cmd --region {seqid}:{start}-{end}",
        ...     GenomicChunk("chunk_0", "chr1", 999, 2000),
        ... )
        'cmd --region chr1:1000-2000'
    """
    # Convert 0-based half-open to 1-based inclusive for CLI
    # Internal: start=999, end=2000 (0-based half-open)
    # CLI:      start=1000, end=2000 (1-based inclusive)
    start_1based = chunk.start + 1
    end_1based = chunk.end  # half-open end equals 1-based inclusive end

    variables = {
        "chunk_id": chunk.chunk_id,
        "seqid": chunk.seqid,
        # 1-based coordinates for CLI commands
        "start": str(start_1based),
        "end": str(end_1based),
        # 0-based coordinates for internal use
        "start_0": str(chunk.start),
        "end_0": str(chunk.end),
        "size": str(chunk.size),
    }

    if chunk_plan_path is not None:
        variables["chunk_plan"] = str(chunk_plan_path)

    if output_dir is not None:
        variables["output_dir"] = str(output_dir)

    if extra_vars:
        variables.update(extra_vars)

    # Check for unknown placeholders
    placeholders = set(re.findall(r"\{(\w+)\}", template))
    unknown = placeholders - set(variables.keys())
    if unknown:
        raise KeyError(f"Unknown placeholders in template: {unknown}")

    return template.format(**variables)


# =============================================================================
# Task Generator
# =============================================================================


class TaskGenerator:
    """Generate task files for HyperShell or similar tools.

    Each chunk becomes one line/command in the task file.
    Users execute via their preferred parallel executor.

    Example:
        >>> gen = TaskGenerator(chunk_plan)
        >>> task_file = gen.generate(
        ...     command_template="echo {chunk_id}: {seqid}:{start}-{end}",
        ...     output_path="tasks.txt",
        ... )
        >>> # Execute with: hs cluster tasks.txt --num-tasks 32
    """

    def __init__(self, chunk_plan: ChunkPlan) -> None:
        """Initialize the generator.

        Args:
            chunk_plan: ChunkPlan containing chunks to process.
        """
        self.chunk_plan = chunk_plan

    def generate(
        self,
        command_template: str,
        output_path: Path | str,
        chunk_plan_output: Path | str | None = None,
        output_dir: Path | str | None = None,
        include_logging: bool = False,
        log_dir: Path | str | None = None,
        extra_vars: dict[str, str] | None = None,
    ) -> TaskFile:
        """Generate task file.

        Args:
            command_template: Command with placeholders:
                {chunk_id} - Chunk identifier
                {seqid} - Scaffold/chromosome name
                {start} - Start coordinate (1-based, for CLI commands)
                {end} - End coordinate (1-based inclusive, for CLI commands)
                {start_0} - Start coordinate (0-based, for internal use)
                {end_0} - End coordinate (0-based exclusive, for internal use)
                {size} - Chunk size in bases
                {chunk_plan} - Path to chunk plan JSON
                {output_dir} - Output directory
            output_path: Where to write task file.
            chunk_plan_output: Where to save chunk plan JSON.
                Required if template uses {chunk_plan}.
            output_dir: Output directory for results.
                Creates {output_dir} placeholder.
            include_logging: Add stdout/stderr redirection per task.
            log_dir: Directory for log files (defaults to output_dir/logs).
            extra_vars: Additional template variables.

        Returns:
            TaskFile with metadata.

        Note:
            The {start} and {end} placeholders output 1-based inclusive
            coordinates suitable for CLI commands like --region chr1:1000-2000.
            For 0-based half-open coordinates, use {start_0} and {end_0}.

        Example:
            >>> task_file = gen.generate(
            ...     "helixforge confidence --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o {output_dir}/{chunk_id}.tsv",
            ...     output_path="tasks.txt",
            ...     output_dir="outputs/",
            ... )
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        chunk_plan_path = None
        if chunk_plan_output is not None:
            chunk_plan_path = Path(chunk_plan_output)
            chunk_plan_path.parent.mkdir(parents=True, exist_ok=True)
            self.chunk_plan.save(chunk_plan_path)

        output_dir_path = Path(output_dir) if output_dir else None
        if output_dir_path:
            output_dir_path.mkdir(parents=True, exist_ok=True)

        if log_dir:
            log_dir_path = Path(log_dir)
        elif output_dir_path:
            log_dir_path = output_dir_path / "logs"
        else:
            log_dir_path = None

        if log_dir_path:
            log_dir_path.mkdir(parents=True, exist_ok=True)

        # Generate tasks
        tasks = []
        for chunk in self.chunk_plan.chunks:
            cmd = format_command(
                template=command_template,
                chunk=chunk,
                chunk_plan_path=chunk_plan_path,
                output_dir=output_dir_path,
                extra_vars=extra_vars,
            )

            if include_logging and log_dir_path:
                log_file = log_dir_path / f"{chunk.chunk_id}.log"
                cmd = f"{cmd} > {log_file} 2>&1"

            tasks.append(cmd)

        # Write task file
        with open(output_path, "w") as f:
            for task in tasks:
                f.write(task + "\n")

        logger.info(f"Generated {len(tasks)} tasks in {output_path}")

        return TaskFile(
            path=output_path,
            n_tasks=len(tasks),
            chunk_plan_path=chunk_plan_path,
            command_template=command_template,
        )

    def generate_with_wrapper(
        self,
        output_path: Path | str,
        wrapper_script: Path | str,
        chunk_plan_output: Path | str,
    ) -> TaskFile:
        """Generate task file that calls a wrapper script.

        Wrapper script receives chunk_id as argument and reads
        chunk details from the plan JSON. More flexible for
        complex workflows.

        Task file contains lines like:
            bash wrapper.sh chunk_0000
            bash wrapper.sh chunk_0001

        Args:
            output_path: Where to write task file.
            wrapper_script: Path to wrapper script.
            chunk_plan_output: Where to save chunk plan JSON.

        Returns:
            TaskFile with metadata.
        """
        output_path = Path(output_path)
        wrapper_script = Path(wrapper_script)
        chunk_plan_path = Path(chunk_plan_output)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        chunk_plan_path.parent.mkdir(parents=True, exist_ok=True)

        # Save chunk plan
        self.chunk_plan.save(chunk_plan_path)

        # Generate tasks
        tasks = []
        for chunk in self.chunk_plan.chunks:
            cmd = f"bash {wrapper_script} {chunk.chunk_id}"
            tasks.append(cmd)

        # Write task file
        with open(output_path, "w") as f:
            for task in tasks:
                f.write(task + "\n")

        logger.info(f"Generated {len(tasks)} wrapper tasks in {output_path}")

        return TaskFile(
            path=output_path,
            n_tasks=len(tasks),
            chunk_plan_path=chunk_plan_path,
            command_template=f"bash {wrapper_script} {{chunk_id}}",
        )

    @staticmethod
    def generate_wrapper_script(
        output_path: Path | str,
        chunk_plan_path: Path | str,
        commands: list[str],
        setup_commands: list[str] | None = None,
        output_dir: Path | str | None = None,
    ) -> Path:
        """Generate a wrapper script for complex per-chunk workflows.

        The script receives a chunk_id argument and sets up environment
        variables for the chunk ($CHUNK_ID, $SEQID, $START, $END).

        Args:
            output_path: Where to write wrapper script.
            chunk_plan_path: Path to chunk plan JSON.
            commands: Commands to run (can use $SEQID, $START, $END, $CHUNK_ID).
            setup_commands: Optional setup (module loads, conda activate, etc.).
            output_dir: Optional output directory.

        Returns:
            Path to generated wrapper script.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        setup_section = ""
        if setup_commands:
            setup_section = "\n".join(setup_commands)

        commands_section = "\n".join(commands)

        output_dir_line = ""
        if output_dir:
            output_dir_line = f'\nOUTPUT_DIR="{output_dir}"\nmkdir -p "$OUTPUT_DIR"'

        script = f'''#!/bin/bash
set -euo pipefail

# HelixForge chunk processing wrapper
# Generated by helixforge parallel

CHUNK_ID="${{1:?Chunk ID required}}"
CHUNK_PLAN="{chunk_plan_path}"
{output_dir_line}

# Parse chunk info from JSON
read SEQID START END < <(python3 -c "
import json, sys
with open('${{CHUNK_PLAN}}') as f:
    plan = json.load(f)
for chunk in plan['chunks']:
    if chunk['chunk_id'] == '${{CHUNK_ID}}':
        print(chunk['seqid'], chunk['start'], chunk['end'])
        sys.exit(0)
print('ERROR: Chunk not found', file=sys.stderr)
sys.exit(1)
")

export CHUNK_ID SEQID START END
{f'export OUTPUT_DIR' if output_dir else ''}

echo "[$(date)] Processing chunk ${{CHUNK_ID}}: ${{SEQID}}:${{START}}-${{END}}"

# Setup commands
{setup_section}

# Main commands
{commands_section}

echo "[$(date)] Chunk ${{CHUNK_ID}} complete"
'''

        output_path.write_text(script)
        output_path.chmod(output_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)

        logger.info(f"Generated wrapper script: {output_path}")

        return output_path


# =============================================================================
# Convenience Functions
# =============================================================================


def generate_simple_task_file(
    chunk_plan: ChunkPlan,
    helixforge_command: str,
    output_dir: Path | str,
    task_file: Path | str,
    chunk_plan_output: Path | str | None = None,
) -> TaskFile:
    """Quick generation for standard HelixForge commands.

    Assumes command supports --chunk-id and --region flags.
    Automatically adds --chunk-id {chunk_id} --region {seqid}:{start}-{end}
    and output redirection.

    Args:
        chunk_plan: ChunkPlan to generate tasks for.
        helixforge_command: Base command without chunk arguments.
            e.g., "helixforge confidence --helixer-h5 data.h5 --genome genome.fa"
        output_dir: Output directory for results.
        task_file: Where to write task file.
        chunk_plan_output: Optional path to save chunk plan.

    Returns:
        TaskFile with metadata.

    Example:
        >>> generate_simple_task_file(
        ...     plan,
        ...     "helixforge confidence --helixer-h5 data.h5 --helixer-gff data.gff3 --genome genome.fa",
        ...     Path("outputs/"),
        ...     Path("tasks.txt"),
        ... )
    """
    output_dir = Path(output_dir)

    # Build full command template
    template = (
        f"{helixforge_command} "
        f"--chunk-id {{chunk_id}} "
        f"--region {{seqid}}:{{start}}-{{end}} "
        f"-o {{output_dir}}/{{chunk_id}}.tsv"
    )

    gen = TaskGenerator(chunk_plan)
    return gen.generate(
        command_template=template,
        output_path=task_file,
        chunk_plan_output=chunk_plan_output,
        output_dir=output_dir,
    )


def estimate_parallelism(
    n_chunks: int,
    available_cores: int,
    overhead_factor: float = 1.5,
) -> int:
    """Estimate optimal parallelism for HyperShell.

    Args:
        n_chunks: Total number of chunks.
        available_cores: Available CPU cores.
        overhead_factor: Factor for overhead (default 1.5).

    Returns:
        Recommended parallelism value.
    """
    # Don't exceed chunk count or core count
    parallelism = min(n_chunks, available_cores)

    # For good load balancing, aim for 2-4 chunks per worker
    if n_chunks > available_cores * 4:
        parallelism = available_cores
    elif n_chunks > available_cores * 2:
        parallelism = max(1, available_cores // 2)

    return parallelism


def generate_hypershell_command(
    task_file: Path | str,
    parallelism: int | None = None,
    timeout: int | None = None,
) -> str:
    """Generate HyperShell cluster command.

    Args:
        task_file: Path to task file.
        parallelism: Number of parallel workers (default: auto).
        timeout: Task timeout in seconds.

    Returns:
        HyperShell command string.
    """
    cmd = f"hs cluster {task_file}"

    if parallelism:
        cmd += f" --num-tasks {parallelism}"

    if timeout:
        cmd += f" --task-timeout {timeout}"

    return cmd
