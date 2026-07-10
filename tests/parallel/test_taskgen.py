"""Tests for parallel/taskgen.py — v1 command-template model (realign).

Floor: 10. Covers placeholder expansion (including the v3-only ``{id_start}``),
the wrapper + setup-command + per-task logging emission, and the unknown-
placeholder guard. Concrete literal coordinates; both strands via the plan
fixture (chr1 plus-strand chunks + a chr2 minus-strand chunk).
"""

import pytest

from helixforge.parallel.plan import partition_by_strategy, reserve_id_ranges
from helixforge.parallel.taskgen import (
    format_command,
    generate_task_file,
)


@pytest.fixture
def plan(synthetic_genome):
    p = partition_by_strategy(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        strategy="genes", chunk_size=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(p)
    return p


# --------------------------------------------------------------------------
# format_command — placeholders
# --------------------------------------------------------------------------

def test_format_command_all_placeholders(plan):
    chunk = plan.chunks[0]  # chr1:0-3500, id_base 1
    cmd = format_command(
        "id={chunk_id} region={region} seqid={seqid} start={start} end={end} "
        "s0={start_0} e0={end_0} size={size} id_start={id_start} "
        "novel={novel_start}",
        chunk,
    )
    assert cmd == (
        "id=chunk_0000 region=chr1:0-3500 seqid=chr1 start=1 end=3500 "
        "s0=0 e0=3500 size=3500 id_start=1 novel=90000"
    )


def test_format_command_id_start_is_reserved_base(plan):
    # The second chunk's reserved id base is 3 (chunk 0 consumed [1,3)).
    assert format_command("{id_start}", plan.chunks[1]) == "3"
    assert format_command("{novel_start}", plan.chunks[1]) == "90002"


def test_format_command_minus_strand_chunk(plan):
    # chr2 carries the single minus-strand locus as a whole-scaffold chunk.
    chunk = plan.chunks[-1]
    cmd = format_command("{region} {seqid} {start}-{end} {id_start}", chunk)
    assert cmd == "chr2:0-5000 chr2 1-5000 6"


def test_format_command_output_dir(plan):
    from pathlib import Path
    cmd = format_command("-o {output_dir}/{chunk_id}.gff3", plan.chunks[0],
                         output_dir=Path("out"))
    assert cmd == "-o out/chunk_0000.gff3"


def test_format_command_unknown_placeholder_raises(plan):
    with pytest.raises(KeyError, match="unknown placeholder"):
        format_command("--bogus {not_a_field}", plan.chunks[0])


# --------------------------------------------------------------------------
# generate_task_file — one line per chunk
# --------------------------------------------------------------------------

def test_generate_task_file_one_line_per_chunk(plan, tmp_path):
    out = tmp_path / "tasks.txt"
    tf = generate_task_file(
        plan, "hf reconcile --region {region} --id-base {id_start}", out)
    lines = out.read_text().splitlines()
    assert tf.n_tasks == 4
    assert len(lines) == 4
    assert lines[0] == "hf reconcile --region chr1:0-3500 --id-base 1"
    assert lines[3] == "hf reconcile --region chr2:0-5000 --id-base 6"


def test_generate_task_file_output_dir_placeholder(plan, tmp_path):
    out = tmp_path / "tasks.txt"
    generate_task_file(plan, "run {output_dir}/{chunk_id}", out,
                       output_dir=tmp_path / "chunks")
    first = out.read_text().splitlines()[0]
    assert first == f"run {tmp_path / 'chunks'}/chunk_0000"
    assert (tmp_path / "chunks").is_dir()


def test_generate_task_file_include_logging(plan, tmp_path):
    out = tmp_path / "tasks.txt"
    generate_task_file(plan, "run {chunk_id}", out,
                       output_dir=tmp_path / "chunks", include_logging=True)
    first = out.read_text().splitlines()[0]
    assert first.endswith(f"> {tmp_path / 'chunks' / 'logs' / 'chunk_0000.log'} 2>&1")
    assert (tmp_path / "chunks" / "logs").is_dir()


# --------------------------------------------------------------------------
# Wrapper + setup commands
# --------------------------------------------------------------------------

def test_generate_task_file_wrapper_and_setup(plan, tmp_path):
    out = tmp_path / "tasks.txt"
    wrapper = tmp_path / "run_chunk.sh"
    tf = generate_task_file(
        plan, "hf reconcile --region {region}", out,
        wrapper=wrapper, wrapper_setup=("module load helixforge", "conda activate hf"),
    )
    assert tf.wrapper_path == wrapper
    # Each task invokes the wrapper with the expanded command as one argument.
    first = out.read_text().splitlines()[0]
    assert first.startswith(f"bash {wrapper} ")
    assert "hf reconcile --region chr1:0-3500" in first
    # The wrapper carries the setup lines and execs the passed command.
    wtext = wrapper.read_text()
    assert "module load helixforge" in wtext
    assert "conda activate hf" in wtext
    assert 'eval "$1"' in wtext


def test_generate_task_file_wrapper_is_executable(plan, tmp_path):
    import os
    out = tmp_path / "tasks.txt"
    wrapper = tmp_path / "w.sh"
    generate_task_file(plan, "echo {chunk_id}", out, wrapper=wrapper)
    assert os.access(wrapper, os.X_OK)


def test_generate_task_file_wrapper_with_logging(plan, tmp_path):
    out = tmp_path / "tasks.txt"
    wrapper = tmp_path / "w.sh"
    generate_task_file(plan, "echo {chunk_id}", out, output_dir=tmp_path / "o",
                       wrapper=wrapper, include_logging=True)
    first = out.read_text().splitlines()[0]
    assert first.startswith(f"bash {wrapper} ")
    assert first.endswith("chunk_0000.log 2>&1")


# --------------------------------------------------------------------------
# HyperShell / GNU-Parallel executor helpers (v1 surface parity)
# --------------------------------------------------------------------------


def test_generate_hypershell_command_bare():
    from helixforge.parallel.taskgen import generate_hypershell_command

    assert generate_hypershell_command("tasks.txt") == "hs cluster tasks.txt"


def test_generate_hypershell_command_full():
    from helixforge.parallel.taskgen import generate_hypershell_command

    cmd = generate_hypershell_command("tasks.txt", parallelism=8, timeout=3600)
    assert cmd == "hs cluster tasks.txt --num-tasks 8 --task-timeout 3600"


def test_estimate_parallelism_caps_at_cores():
    from helixforge.parallel.taskgen import estimate_parallelism

    # 100 chunks, 16 cores → capped at cores
    assert estimate_parallelism(100, 16) == 16
    # fewer chunks than cores → capped at chunk count
    assert estimate_parallelism(3, 16) == 3
    # never below 1
    assert estimate_parallelism(0, 16) == 1


def test_get_optimal_workers_respects_cap():
    from helixforge.parallel.taskgen import get_optimal_workers

    assert get_optimal_workers(1) == 1
    assert get_optimal_workers(1000) >= 1
    assert get_optimal_workers() >= 1


def test_write_example_sbatch_hypershell(tmp_path):
    from helixforge.parallel.taskgen import write_example_sbatch

    path = write_example_sbatch(tmp_path / "run.sbatch")
    text = path.read_text()
    assert text.startswith("#!/bin/bash")
    assert "#SBATCH --job-name=helixforge" in text
    assert "hs cluster tasks.txt --num-tasks ${SLURM_CPUS_PER_TASK}" in text


def test_write_example_sbatch_parallel(tmp_path):
    from helixforge.parallel.taskgen import write_example_sbatch

    path = write_example_sbatch(tmp_path / "run.sbatch", executor="parallel")
    text = path.read_text()
    assert "parallel -j ${SLURM_CPUS_PER_TASK} < tasks.txt" in text


def test_write_example_sbatch_rejects_unknown_executor(tmp_path):
    from helixforge.parallel.taskgen import write_example_sbatch

    with pytest.raises(ValueError):
        write_example_sbatch(tmp_path / "run.sbatch", executor="bogus")


# --------------------------------------------------------------------------
# Package re-export hub (v1-style convenience surface)
# --------------------------------------------------------------------------


def test_parallel_package_reexports_public_surface():
    import helixforge.parallel as p

    for name in (
        "partition_genome", "reserve_id_ranges", "Plan", "Chunk",
        "ChunkStrategy", "plan_chunks", "generate_task_file", "format_command",
        "generate_hypershell_command", "write_example_sbatch", "run_genome",
        "aggregate", "boundary_stitch", "suggest_plan",
    ):
        assert name in p.__all__, name
        assert hasattr(p, name), name
