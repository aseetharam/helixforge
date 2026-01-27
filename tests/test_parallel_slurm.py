"""Tests for helixforge.parallel.slurm module.

Tests cover minimal SLURM utilities:
- Environment detection (detect_slurm_environment, is_slurm_job, get_slurm_task_id)
- Resource access (get_slurm_resources)
- Chunk access (get_chunk_for_task)
- Example SBATCH generation (write_example_sbatch)
"""

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from helixforge.parallel.slurm import (
    detect_slurm_environment,
    is_slurm_job,
    get_slurm_task_id,
    get_slurm_resources,
    get_chunk_for_task,
    write_example_sbatch,
    EXAMPLE_SBATCH_HYPERSHELL,
    EXAMPLE_SBATCH_PARALLEL,
)
from helixforge.parallel.chunker import ChunkPlan, ChunkStrategy, GenomicChunk


# =============================================================================
# Environment Detection Tests
# =============================================================================


class TestDetectSlurmEnvironment:
    """Tests for detect_slurm_environment function."""

    def test_in_slurm_job(self):
        """Test detection when in SLURM environment."""
        slurm_env = {
            "SLURM_JOB_ID": "12345",
            "SLURM_ARRAY_TASK_ID": "5",
            "SLURM_CPUS_PER_TASK": "4",
            "SLURM_MEM_PER_NODE": "16G",
        }
        with patch.dict(os.environ, slurm_env, clear=False):
            env = detect_slurm_environment()

            assert env is not None
            assert env.get("SLURM_JOB_ID") == "12345"
            assert env.get("SLURM_ARRAY_TASK_ID") == "5"
            assert env.get("SLURM_CPUS_PER_TASK") == "4"

    def test_not_in_slurm(self):
        """Test detection when not in SLURM environment."""
        # Remove SLURM variables
        clean_env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
        with patch.dict(os.environ, clean_env, clear=True):
            env = detect_slurm_environment()
            assert env is None

    def test_partial_slurm_env(self):
        """Test with only some SLURM variables."""
        with patch.dict(
            os.environ,
            {"SLURM_JOB_ID": "12345"},
            clear=False,
        ):
            # First ensure we clear other SLURM vars for test isolation
            test_env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
            test_env["SLURM_JOB_ID"] = "12345"
            with patch.dict(os.environ, test_env, clear=True):
                env = detect_slurm_environment()
                assert env is not None
                assert "SLURM_JOB_ID" in env
                assert env.get("SLURM_ARRAY_TASK_ID") is None


class TestIsSlurmJob:
    """Tests for is_slurm_job function."""

    def test_is_slurm_job_true(self):
        """Test when in SLURM job."""
        with patch.dict(os.environ, {"SLURM_JOB_ID": "12345"}, clear=False):
            assert is_slurm_job() is True

    def test_is_slurm_job_false(self):
        """Test when not in SLURM job."""
        clean_env = {k: v for k, v in os.environ.items() if k != "SLURM_JOB_ID"}
        with patch.dict(os.environ, clean_env, clear=True):
            assert is_slurm_job() is False


class TestGetSlurmTaskId:
    """Tests for get_slurm_task_id function."""

    def test_get_task_id(self):
        """Test getting SLURM array task ID."""
        with patch.dict(os.environ, {"SLURM_ARRAY_TASK_ID": "7"}, clear=False):
            assert get_slurm_task_id() == 7

    def test_task_id_not_set(self):
        """Test when not in array job."""
        clean_env = {k: v for k, v in os.environ.items() if k != "SLURM_ARRAY_TASK_ID"}
        with patch.dict(os.environ, clean_env, clear=True):
            assert get_slurm_task_id() is None


# =============================================================================
# Resource Access Tests
# =============================================================================


class TestGetSlurmResources:
    """Tests for get_slurm_resources function."""

    def test_get_resources_in_slurm(self):
        """Test getting resources when in SLURM."""
        slurm_env = {
            "SLURM_JOB_ID": "12345",
            "SLURM_CPUS_PER_TASK": "8",
            "SLURM_MEM_PER_NODE": "32G",
            "SLURM_ARRAY_TASK_ID": "3",
            "SLURM_NODELIST": "node01",
        }
        with patch.dict(os.environ, slurm_env, clear=False):
            resources = get_slurm_resources()

            assert resources["cpus"] == 8
            assert resources["memory_mb"] == 32 * 1024  # 32G in MB
            assert resources["task_id"] == 3
            assert resources["node"] == "node01"

    def test_get_resources_not_in_slurm(self):
        """Test getting resources when not in SLURM."""
        clean_env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
        with patch.dict(os.environ, clean_env, clear=True):
            resources = get_slurm_resources()

            assert resources["cpus"] == 1
            assert resources["memory_mb"] is None
            assert resources["task_id"] is None
            assert resources["node"] is None

    def test_memory_parsing_mb(self):
        """Test parsing memory in MB."""
        slurm_env = {
            "SLURM_JOB_ID": "12345",
            "SLURM_MEM_PER_NODE": "16000M",
        }
        with patch.dict(os.environ, slurm_env, clear=False):
            resources = get_slurm_resources()
            assert resources["memory_mb"] == 16000

    def test_memory_parsing_kb(self):
        """Test parsing memory in KB."""
        slurm_env = {
            "SLURM_JOB_ID": "12345",
            "SLURM_MEM_PER_NODE": "16384000K",
        }
        with patch.dict(os.environ, slurm_env, clear=False):
            resources = get_slurm_resources()
            assert resources["memory_mb"] == 16000  # 16384000K / 1024

    def test_memory_parsing_no_unit(self):
        """Test parsing memory without unit (assumes MB)."""
        slurm_env = {
            "SLURM_JOB_ID": "12345",
            "SLURM_MEM_PER_NODE": "8000",
        }
        with patch.dict(os.environ, slurm_env, clear=False):
            resources = get_slurm_resources()
            assert resources["memory_mb"] == 8000


# =============================================================================
# Chunk Access Tests
# =============================================================================


class TestGetChunkForTask:
    """Tests for get_chunk_for_task function."""

    @pytest.fixture
    def chunk_plan_file(self, tmp_path) -> Path:
        """Create a chunk plan file."""
        chunks = [
            GenomicChunk("chunk_0000", "chr1", 0, 1000000),
            GenomicChunk("chunk_0001", "chr1", 1000000, 2000000),
            GenomicChunk("chunk_0002", "chr2", 0, 500000),
        ]
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SCAFFOLD,
            chunks=chunks,
            total_bases=2500000,
        )
        path = tmp_path / "plan.json"
        plan.save(path)
        return path

    def test_get_chunk_explicit_id(self, chunk_plan_file):
        """Test getting chunk with explicit task ID."""
        chunk = get_chunk_for_task(chunk_plan_file, task_id=1)

        assert chunk.chunk_id == "chunk_0001"
        assert chunk.seqid == "chr1"
        assert chunk.start == 1000000
        assert chunk.end == 2000000

    def test_get_chunk_first(self, chunk_plan_file):
        """Test getting first chunk."""
        chunk = get_chunk_for_task(chunk_plan_file, task_id=0)

        assert chunk.chunk_id == "chunk_0000"
        assert chunk.seqid == "chr1"
        assert chunk.start == 0

    def test_get_chunk_from_env(self, chunk_plan_file):
        """Test getting chunk from SLURM environment."""
        with patch.dict(os.environ, {"SLURM_ARRAY_TASK_ID": "2"}, clear=False):
            chunk = get_chunk_for_task(chunk_plan_file)

            assert chunk.chunk_id == "chunk_0002"
            assert chunk.seqid == "chr2"

    def test_get_chunk_no_task_id(self, chunk_plan_file):
        """Test error when no task ID available."""
        clean_env = {k: v for k, v in os.environ.items() if k != "SLURM_ARRAY_TASK_ID"}
        with patch.dict(os.environ, clean_env, clear=True):
            with pytest.raises(RuntimeError, match="Task ID not provided"):
                get_chunk_for_task(chunk_plan_file)

    def test_get_chunk_invalid_id(self, chunk_plan_file):
        """Test error for out-of-range task ID."""
        with pytest.raises(IndexError, match="exceeds chunk count"):
            get_chunk_for_task(chunk_plan_file, task_id=100)


# =============================================================================
# Example SBATCH Tests
# =============================================================================


class TestWriteExampleSbatch:
    """Tests for write_example_sbatch function."""

    def test_write_hypershell_script(self, tmp_path):
        """Test writing HyperShell SBATCH script."""
        output_path = tmp_path / "run.sbatch"
        result = write_example_sbatch(output_path, executor="hypershell")

        assert result == output_path
        assert output_path.exists()

        content = output_path.read_text()
        assert "#!/bin/bash" in content
        assert "#SBATCH" in content
        assert "hs cluster" in content
        assert "SLURM_CPUS_PER_TASK" in content

    def test_write_parallel_script(self, tmp_path):
        """Test writing GNU Parallel SBATCH script."""
        output_path = tmp_path / "run.sbatch"
        result = write_example_sbatch(output_path, executor="parallel")

        assert result == output_path
        assert output_path.exists()

        content = output_path.read_text()
        assert "#!/bin/bash" in content
        assert "#SBATCH" in content
        assert "parallel -j" in content
        assert "SLURM_CPUS_PER_TASK" in content

    def test_write_creates_parent_dirs(self, tmp_path):
        """Test that parent directories are created."""
        output_path = tmp_path / "deep" / "nested" / "run.sbatch"
        result = write_example_sbatch(output_path)

        assert output_path.exists()

    def test_write_unknown_executor(self, tmp_path):
        """Test error for unknown executor."""
        with pytest.raises(ValueError, match="Unknown executor"):
            write_example_sbatch(tmp_path / "run.sbatch", executor="unknown")


class TestExampleSbatchTemplates:
    """Tests for SBATCH template constants."""

    def test_hypershell_template_content(self):
        """Test HyperShell template has expected content."""
        assert "#!/bin/bash" in EXAMPLE_SBATCH_HYPERSHELL
        assert "#SBATCH --job-name=helixforge" in EXAMPLE_SBATCH_HYPERSHELL
        assert "hs cluster" in EXAMPLE_SBATCH_HYPERSHELL
        assert "tasks.txt" in EXAMPLE_SBATCH_HYPERSHELL

    def test_parallel_template_content(self):
        """Test GNU Parallel template has expected content."""
        assert "#!/bin/bash" in EXAMPLE_SBATCH_PARALLEL
        assert "#SBATCH --job-name=helixforge" in EXAMPLE_SBATCH_PARALLEL
        assert "parallel -j" in EXAMPLE_SBATCH_PARALLEL
        assert "tasks.txt" in EXAMPLE_SBATCH_PARALLEL

    def test_templates_have_resource_requests(self):
        """Test templates have proper SLURM resource requests."""
        for template in [EXAMPLE_SBATCH_HYPERSHELL, EXAMPLE_SBATCH_PARALLEL]:
            assert "#SBATCH --partition" in template
            assert "#SBATCH --time" in template
            assert "#SBATCH --mem" in template
            assert "#SBATCH --cpus-per-task" in template
