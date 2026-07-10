"""Tests for parallel/tasks.py — per-chunk configs + dispatch (Phase 16 D2).

Floor: 7. ``run_pipeline`` is mocked (no real Mikado/tool chain). The local pool
test injects a ThreadPoolExecutor so the mock is observed in-process.
"""

from concurrent.futures import ThreadPoolExecutor
from unittest import mock

import pytest

from helixforge.parallel.plan import partition_genome, reserve_id_ranges
from helixforge.parallel.tasks import (
    build_chunk_configs,
    default_reconcile_template,
    pipeline_config_to_reconcile_argv,
    run_local,
    write_shell_driver,
    write_slurm_array,
)
from helixforge.reconcile.pipeline import PipelineConfig


@pytest.fixture
def plan(synthetic_genome):
    p = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(p)
    return p


@pytest.fixture
def base_config(synthetic_genome):
    return PipelineConfig(
        genome_fasta=synthetic_genome["fasta"],
        helixer_gff3=synthetic_genome["gff3"],
        output_prefix="run",
        procs=2,
    )


# --------------------------------------------------------------------------
# build_chunk_configs
# --------------------------------------------------------------------------

def test_build_chunk_configs_region_and_ids(plan, base_config):
    configs = build_chunk_configs(plan, base_config)
    assert len(configs) == 4
    assert configs[0].region == "chr1:0-3500"
    assert configs[0].id_base == 1
    assert configs[0].novel_base == 90000
    assert configs[1].region == "chr1:3500-10000"
    assert configs[1].id_base == 3
    assert configs[3].region == "chr2"
    assert configs[3].id_base == 6


def test_build_chunk_configs_per_chunk_paths(plan, base_config):
    configs = build_chunk_configs(plan, base_config)
    assert configs[0].output_prefix == "run.chunk_0000"
    assert configs[0].id_map_path == "run.chunk_0000.id_map.json"
    # report_path/work_dir re-derived from the per-chunk prefix (not shared).
    assert configs[0].report_path == "run.chunk_0000.report.tsv"
    assert configs[1].id_map_path == "run.chunk_0001.id_map.json"
    assert configs[0].chunk_id == "chunk_0000"


def test_build_chunk_configs_preserves_base_knobs(plan, base_config):
    configs = build_chunk_configs(plan, base_config)
    assert all(c.procs == 2 for c in configs)
    assert all(c.genome_fasta == base_config.genome_fasta for c in configs)


# --------------------------------------------------------------------------
# argv rendering
# --------------------------------------------------------------------------

def test_reconcile_argv_carries_region_and_id_base(plan, base_config):
    config = build_chunk_configs(plan, base_config)[1]
    argv = pipeline_config_to_reconcile_argv(config)
    assert argv[:2] == ["helixforge", "reconcile"]
    assert "--region" in argv and argv[argv.index("--region") + 1] == "chr1:3500-10000"
    assert "--id-base" in argv and argv[argv.index("--id-base") + 1] == "3"
    assert "--novel-base" in argv
    assert "--chunk-id" in argv and argv[argv.index("--chunk-id") + 1] == "chunk_0001"
    # bool flags render as their on/off form, never as "--flag value".
    assert "--pad" in argv
    assert "--no-pad" not in argv


def test_default_reconcile_template_has_placeholders(base_config):
    template = default_reconcile_template(base_config)
    # The per-chunk varying fields are placeholders, not baked values.
    assert "--region {region}" in template
    assert "--id-base {id_start}" in template
    assert "--novel-base {novel_start}" in template
    assert "--output-prefix {output_dir}/{chunk_id}" in template
    # The attached inputs ARE baked in (works out of the box).
    assert "--genome" in template and "--helixer" in template
    # No literal per-chunk values leaked through (e.g. the base output prefix).
    assert "--output-prefix run " not in template + " "
    assert "{id_start}" in template


def test_default_reconcile_template_expands_per_chunk(plan, base_config):
    from helixforge.parallel.taskgen import format_command
    template = default_reconcile_template(base_config)
    cmd = format_command(template, plan.chunks[1], output_dir="chunks")
    assert "--region chr1:3500-10000" in cmd
    assert "--id-base 3" in cmd
    assert "--output-prefix chunks/chunk_0001" in cmd
    assert "{" not in cmd  # every placeholder resolved


def test_reconcile_argv_skips_none_lists(base_config):
    argv = pipeline_config_to_reconcile_argv(base_config)
    # No StringTie/BAM evidence configured → those repeatable flags absent.
    assert "--stringtie" not in argv
    assert "--bam" not in argv
    # miniprot is None → flag omitted.
    assert "--miniprot" not in argv


# --------------------------------------------------------------------------
# Slurm array + shell driver
# --------------------------------------------------------------------------

def test_write_slurm_array_wellformed(plan, base_config, tmp_path):
    out = tmp_path / "array.sbatch"
    write_slurm_array(plan, base_config, out, cpus_per_task=4, mem="8G", time="04:00:00")
    text = out.read_text()
    assert text.startswith("#!/bin/bash")
    assert "#SBATCH --array=0-3" in text
    assert "#SBATCH --cpus-per-task=4" in text
    assert "#SBATCH --mem=8G" in text
    assert 'case "$SLURM_ARRAY_TASK_ID" in' in text
    assert "  0) helixforge reconcile" in text
    assert "  3) helixforge reconcile" in text


def test_write_slurm_array_max_concurrent(plan, base_config, tmp_path):
    out = tmp_path / "array.sbatch"
    write_slurm_array(plan, base_config, out, max_concurrent=2)
    assert "#SBATCH --array=0-3%2" in out.read_text()


def test_write_shell_driver_lists_all_chunks(plan, base_config, tmp_path):
    out = tmp_path / "driver.sh"
    write_shell_driver(plan, base_config, out)
    text = out.read_text()
    assert text.startswith("#!/bin/bash")
    assert text.count("helixforge reconcile") == 4
    assert "--region chr1:0-3500" in text


def test_write_shell_driver_parallel_gate(plan, base_config, tmp_path):
    out = tmp_path / "driver.sh"
    write_shell_driver(plan, base_config, out, workers=3)
    text = out.read_text()
    assert "MAX_JOBS=3" in text
    assert "wait -n" in text


# --------------------------------------------------------------------------
# Local pool (mocked run_pipeline)
# --------------------------------------------------------------------------

def test_run_local_dispatches_each_chunk(plan, base_config):
    with mock.patch("helixforge.parallel.tasks.run_pipeline") as run:
        prefixes = run_local(plan, base_config, workers=2,
                             _executor_cls=ThreadPoolExecutor)
    assert run.call_count == 4
    assert {str(p) for p in prefixes} == {
        "run.chunk_0000", "run.chunk_0001", "run.chunk_0002", "run.chunk_0003",
    }
    # Each call received a distinct per-chunk config.
    dispatched_regions = {c.args[0].region for c in run.call_args_list}
    assert dispatched_regions == {"chr1:0-3500", "chr1:3500-10000", "chr1:10000-20000", "chr2"}


def test_run_local_propagates_failure(plan, base_config):
    with mock.patch("helixforge.parallel.tasks.run_pipeline", side_effect=RuntimeError("boom")):
        with pytest.raises(RuntimeError):
            run_local(plan, base_config, workers=2, _executor_cls=ThreadPoolExecutor)
