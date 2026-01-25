"""Tests for helixforge.parallel.taskgen module.

Tests cover:
- TaskFile data structure
- TaskGenerator class
- format_command function
- generate_simple_task_file function
- estimate_parallelism function
- generate_hypershell_command function
"""

import json
from pathlib import Path

import pytest

from helixforge.parallel.taskgen import (
    TaskFile,
    TaskGenerator,
    format_command,
    generate_simple_task_file,
    estimate_parallelism,
    generate_hypershell_command,
)
from helixforge.parallel.chunker import ChunkPlan, ChunkStrategy, GenomicChunk


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def sample_chunks() -> list[GenomicChunk]:
    """Create sample genomic chunks."""
    return [
        GenomicChunk("chunk_0000", "chr1", 0, 1000000),
        GenomicChunk("chunk_0001", "chr1", 1000000, 2000000),
        GenomicChunk("chunk_0002", "chr2", 0, 500000),
    ]


@pytest.fixture
def sample_chunk_plan(sample_chunks) -> ChunkPlan:
    """Create sample chunk plan."""
    return ChunkPlan(
        strategy=ChunkStrategy.BY_SCAFFOLD,
        chunks=sample_chunks,
        total_bases=2500000,
    )


@pytest.fixture
def task_generator(sample_chunk_plan) -> TaskGenerator:
    """Create task generator instance."""
    return TaskGenerator(sample_chunk_plan)


# =============================================================================
# TaskFile Tests
# =============================================================================


class TestTaskFile:
    """Tests for TaskFile data structure."""

    def test_task_file_creation(self, tmp_path):
        """Test TaskFile creation."""
        task_path = tmp_path / "tasks.txt"
        task_path.write_text("echo task1\necho task2\necho task3\n")

        task_file = TaskFile(
            path=task_path,
            n_tasks=3,
            chunk_plan_path=tmp_path / "chunks.json",
            command_template="echo {chunk_id}",
        )

        assert task_file.n_tasks == 3
        assert task_file.path == task_path
        assert task_file.command_template == "echo {chunk_id}"

    def test_task_file_preview(self, tmp_path):
        """Test TaskFile.preview method."""
        task_path = tmp_path / "tasks.txt"
        task_path.write_text("line1\nline2\nline3\nline4\nline5\n")

        task_file = TaskFile(
            path=task_path,
            n_tasks=5,
            chunk_plan_path=None,
            command_template="test",
        )

        preview = task_file.preview(3)
        assert len(preview) == 3
        assert preview == ["line1", "line2", "line3"]

    def test_task_file_preview_empty(self, tmp_path):
        """Test TaskFile.preview on nonexistent file."""
        task_file = TaskFile(
            path=tmp_path / "nonexistent.txt",
            n_tasks=0,
            chunk_plan_path=None,
            command_template="test",
        )

        assert task_file.preview() == []

    def test_task_file_read_all(self, tmp_path):
        """Test TaskFile.read_all method."""
        task_path = tmp_path / "tasks.txt"
        task_path.write_text("cmd1\ncmd2\ncmd3\n")

        task_file = TaskFile(
            path=task_path,
            n_tasks=3,
            chunk_plan_path=None,
            command_template="test",
        )

        all_tasks = task_file.read_all()
        assert len(all_tasks) == 3
        assert all_tasks == ["cmd1", "cmd2", "cmd3"]


# =============================================================================
# format_command Tests
# =============================================================================


class TestFormatCommand:
    """Tests for format_command function."""

    def test_basic_formatting(self, sample_chunks):
        """Test basic placeholder replacement."""
        chunk = sample_chunks[0]
        result = format_command(
            template="echo {chunk_id} {seqid}:{start}-{end}",
            chunk=chunk,
        )
        assert result == "echo chunk_0000 chr1:0-1000000"

    def test_size_placeholder(self, sample_chunks):
        """Test size placeholder."""
        chunk = sample_chunks[0]
        result = format_command(
            template="echo size={size}",
            chunk=chunk,
        )
        assert result == "echo size=1000000"

    def test_chunk_plan_placeholder(self, sample_chunks, tmp_path):
        """Test chunk_plan placeholder."""
        chunk = sample_chunks[0]
        plan_path = tmp_path / "chunks.json"

        result = format_command(
            template="--chunk-plan {chunk_plan}",
            chunk=chunk,
            chunk_plan_path=plan_path,
        )
        assert str(plan_path) in result

    def test_output_dir_placeholder(self, sample_chunks, tmp_path):
        """Test output_dir placeholder."""
        chunk = sample_chunks[0]
        output_dir = tmp_path / "outputs"

        result = format_command(
            template="-o {output_dir}/{chunk_id}.tsv",
            chunk=chunk,
            output_dir=output_dir,
        )
        assert str(output_dir) in result
        assert "chunk_0000.tsv" in result

    def test_extra_vars(self, sample_chunks):
        """Test extra variables."""
        chunk = sample_chunks[0]
        result = format_command(
            template="--genome {genome} --bam {bam}",
            chunk=chunk,
            extra_vars={"genome": "genome.fa", "bam": "rnaseq.bam"},
        )
        assert result == "--genome genome.fa --bam rnaseq.bam"

    def test_unknown_placeholder_error(self, sample_chunks):
        """Test error for unknown placeholder."""
        chunk = sample_chunks[0]
        with pytest.raises(KeyError, match="Unknown placeholders"):
            format_command(
                template="echo {unknown}",
                chunk=chunk,
            )

    def test_complex_template(self, sample_chunks, tmp_path):
        """Test complex command template."""
        chunk = sample_chunks[1]
        output_dir = tmp_path / "outputs"

        result = format_command(
            template="helixforge confidence --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o {output_dir}/{chunk_id}.tsv",
            chunk=chunk,
            output_dir=output_dir,
        )

        assert "helixforge confidence" in result
        assert "--chunk-id chunk_0001" in result
        assert "--region chr1:1000000-2000000" in result
        assert f"-o {output_dir}/chunk_0001.tsv" in result


# =============================================================================
# TaskGenerator Tests
# =============================================================================


class TestTaskGenerator:
    """Tests for TaskGenerator class."""

    def test_generator_creation(self, sample_chunk_plan):
        """Test generator initialization."""
        gen = TaskGenerator(sample_chunk_plan)
        assert gen.chunk_plan == sample_chunk_plan

    def test_generate_basic(self, task_generator, tmp_path):
        """Test basic task file generation."""
        task_file = task_generator.generate(
            command_template="echo {chunk_id}",
            output_path=tmp_path / "tasks.txt",
        )

        assert task_file.n_tasks == 3
        assert task_file.path.exists()

        content = task_file.path.read_text()
        assert "echo chunk_0000" in content
        assert "echo chunk_0001" in content
        assert "echo chunk_0002" in content

    def test_generate_with_output_dir(self, task_generator, tmp_path):
        """Test generation with output directory."""
        output_dir = tmp_path / "outputs"

        task_file = task_generator.generate(
            command_template="echo {chunk_id} -o {output_dir}/{chunk_id}.txt",
            output_path=tmp_path / "tasks.txt",
            output_dir=output_dir,
        )

        assert output_dir.exists()

        content = task_file.path.read_text()
        assert str(output_dir) in content

    def test_generate_with_chunk_plan_output(self, task_generator, tmp_path):
        """Test generation with chunk plan saving."""
        task_file = task_generator.generate(
            command_template="echo {chunk_id}",
            output_path=tmp_path / "tasks.txt",
            chunk_plan_output=tmp_path / "plan.json",
        )

        assert task_file.chunk_plan_path == tmp_path / "plan.json"
        assert task_file.chunk_plan_path.exists()

        # Verify JSON content
        with open(task_file.chunk_plan_path) as f:
            data = json.load(f)
        assert data["n_chunks"] == 3

    def test_generate_with_logging(self, task_generator, tmp_path):
        """Test generation with logging enabled."""
        output_dir = tmp_path / "outputs"

        task_file = task_generator.generate(
            command_template="echo {chunk_id}",
            output_path=tmp_path / "tasks.txt",
            output_dir=output_dir,
            include_logging=True,
        )

        content = task_file.path.read_text()
        # Each line should have output redirection
        lines = content.strip().split("\n")
        for line in lines:
            assert "> " in line
            assert "2>&1" in line

    def test_generate_with_extra_vars(self, task_generator, tmp_path):
        """Test generation with extra variables."""
        task_file = task_generator.generate(
            command_template="--genome {genome} --chunk {chunk_id}",
            output_path=tmp_path / "tasks.txt",
            extra_vars={"genome": "genome.fa"},
        )

        content = task_file.path.read_text()
        assert "--genome genome.fa" in content

    def test_generate_creates_parent_dirs(self, task_generator, tmp_path):
        """Test that generate creates parent directories."""
        task_file = task_generator.generate(
            command_template="echo {chunk_id}",
            output_path=tmp_path / "deep" / "nested" / "tasks.txt",
        )

        assert task_file.path.exists()

    def test_generate_with_wrapper(self, task_generator, tmp_path):
        """Test wrapper script generation."""
        wrapper_path = tmp_path / "wrapper.sh"
        plan_path = tmp_path / "chunks.json"

        # First save the chunk plan
        task_generator.chunk_plan.save(plan_path)

        task_file = task_generator.generate_with_wrapper(
            output_path=tmp_path / "tasks.txt",
            wrapper_script=wrapper_path,
            chunk_plan_output=plan_path,
        )

        assert task_file.n_tasks == 3

        content = task_file.path.read_text()
        lines = content.strip().split("\n")
        assert lines[0] == f"bash {wrapper_path} chunk_0000"
        assert lines[1] == f"bash {wrapper_path} chunk_0001"
        assert lines[2] == f"bash {wrapper_path} chunk_0002"

    def test_generate_wrapper_script(self, tmp_path):
        """Test static wrapper script generation."""
        wrapper_path = TaskGenerator.generate_wrapper_script(
            output_path=tmp_path / "wrapper.sh",
            chunk_plan_path=tmp_path / "chunks.json",
            commands=[
                "helixforge confidence --region ${SEQID}:${START}-${END}",
                "echo done",
            ],
            setup_commands=["module load python/3.10"],
            output_dir=tmp_path / "outputs",
        )

        assert wrapper_path.exists()

        content = wrapper_path.read_text()
        assert "#!/bin/bash" in content
        assert "set -euo pipefail" in content
        assert "CHUNK_ID=" in content
        assert "helixforge confidence" in content
        assert "module load python/3.10" in content
        assert "echo done" in content

        # Check executable permission
        import stat
        mode = wrapper_path.stat().st_mode
        assert mode & stat.S_IXUSR


# =============================================================================
# Convenience Function Tests
# =============================================================================


class TestGenerateSimpleTaskFile:
    """Tests for generate_simple_task_file function."""

    def test_simple_task_file(self, sample_chunk_plan, tmp_path):
        """Test simple task file generation."""
        task_file = generate_simple_task_file(
            chunk_plan=sample_chunk_plan,
            helixforge_command="helixforge confidence --genome genome.fa",
            output_dir=tmp_path / "outputs",
            task_file=tmp_path / "tasks.txt",
        )

        assert task_file.n_tasks == 3
        assert task_file.path.exists()

        content = task_file.path.read_text()
        # Check command structure
        assert "helixforge confidence" in content
        assert "--chunk-id chunk_0000" in content
        assert "--region chr1:0-1000000" in content
        assert "-o " in content

    def test_simple_task_file_with_plan_output(self, sample_chunk_plan, tmp_path):
        """Test simple task file with chunk plan output."""
        task_file = generate_simple_task_file(
            chunk_plan=sample_chunk_plan,
            helixforge_command="helixforge test",
            output_dir=tmp_path / "outputs",
            task_file=tmp_path / "tasks.txt",
            chunk_plan_output=tmp_path / "plan.json",
        )

        assert task_file.chunk_plan_path == tmp_path / "plan.json"
        assert task_file.chunk_plan_path.exists()


class TestEstimateParallelism:
    """Tests for estimate_parallelism function."""

    def test_basic_estimate(self):
        """Test basic parallelism estimate."""
        result = estimate_parallelism(n_chunks=100, available_cores=8)
        assert result == 8  # Limited by cores

    def test_fewer_chunks_than_cores(self):
        """Test when chunks < cores."""
        result = estimate_parallelism(n_chunks=4, available_cores=16)
        assert result == 4  # Limited by chunks

    def test_many_chunks(self):
        """Test with many chunks for load balancing."""
        result = estimate_parallelism(n_chunks=500, available_cores=32)
        assert result == 32  # Use all cores

    def test_moderate_chunks(self):
        """Test with moderate chunk count."""
        result = estimate_parallelism(n_chunks=100, available_cores=32)
        assert result == 16  # Half cores for better load balancing


class TestGenerateHypershellCommand:
    """Tests for generate_hypershell_command function."""

    def test_basic_command(self, tmp_path):
        """Test basic HyperShell command."""
        task_file = tmp_path / "tasks.txt"
        cmd = generate_hypershell_command(task_file)

        assert "hs launch" in cmd
        assert f"< {task_file}" in cmd

    def test_with_parallelism(self, tmp_path):
        """Test with parallelism option."""
        task_file = tmp_path / "tasks.txt"
        cmd = generate_hypershell_command(task_file, parallelism=16)

        assert "--parallelism 16" in cmd

    def test_with_timeout(self, tmp_path):
        """Test with timeout option."""
        task_file = tmp_path / "tasks.txt"
        cmd = generate_hypershell_command(task_file, timeout="2h")

        assert "--timeout 2h" in cmd

    def test_with_all_options(self, tmp_path):
        """Test with all options."""
        task_file = tmp_path / "tasks.txt"
        cmd = generate_hypershell_command(
            task_file,
            parallelism=32,
            timeout="1h",
        )

        assert "hs launch" in cmd
        assert "--parallelism 32" in cmd
        assert "--timeout 1h" in cmd
        assert f"< {task_file}" in cmd
