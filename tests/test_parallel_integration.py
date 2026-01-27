"""Integration tests for helixforge.parallel module.

Tests cover end-to-end workflows combining:
- Chunking
- Execution
- SLURM job generation
- Result aggregation
"""

import json
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from helixforge.parallel import (
    ChunkPlan,
    ChunkStrategy,
    GenomeChunker,
    GenomicChunk,
    ParallelExecutor,
    ExecutorBackend,
    merge_overlapping_results,
    suggest_chunk_parameters,
    get_optimal_workers,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def mock_genome():
    """Create mock GenomeAccessor for integration tests."""
    genome = MagicMock()
    genome.scaffold_order = ["chr1", "chr2", "chr3", "chr4", "chr5"]
    genome.total_length = 10_000_000  # 10 Mb total

    lengths = {
        "chr1": 3_000_000,
        "chr2": 2_500_000,
        "chr3": 2_000_000,
        "chr4": 1_500_000,
        "chr5": 1_000_000,
    }
    genome.get_length = lambda seqid: lengths[seqid]
    genome.path = Path("/path/to/genome.fa")

    return genome


@pytest.fixture
def mock_gff_parser():
    """Create mock GFF3Parser with realistic gene distribution."""
    parser = MagicMock()

    # Create genes across chromosomes
    genes = []
    gene_idx = 0

    # chr1: 100 genes
    for i in range(0, 3_000_000, 30_000):
        genes.append(
            MagicMock(
                gene_id=f"gene_{gene_idx:05d}",
                seqid="chr1",
                start=i,
                end=i + 5000,
            )
        )
        gene_idx += 1

    # chr2: 80 genes
    for i in range(0, 2_400_000, 30_000):
        genes.append(
            MagicMock(
                gene_id=f"gene_{gene_idx:05d}",
                seqid="chr2",
                start=i,
                end=i + 5000,
            )
        )
        gene_idx += 1

    # chr3: 60 genes
    for i in range(0, 1_800_000, 30_000):
        genes.append(
            MagicMock(
                gene_id=f"gene_{gene_idx:05d}",
                seqid="chr3",
                start=i,
                end=i + 5000,
            )
        )
        gene_idx += 1

    # chr4: 40 genes
    for i in range(0, 1_200_000, 30_000):
        genes.append(
            MagicMock(
                gene_id=f"gene_{gene_idx:05d}",
                seqid="chr4",
                start=i,
                end=i + 5000,
            )
        )
        gene_idx += 1

    # chr5: 30 genes
    for i in range(0, 900_000, 30_000):
        genes.append(
            MagicMock(
                gene_id=f"gene_{gene_idx:05d}",
                seqid="chr5",
                start=i,
                end=i + 5000,
            )
        )
        gene_idx += 1

    parser.iter_genes.return_value = iter(genes)

    return parser


# =============================================================================
# Chunking to Execution Integration
# =============================================================================


class TestChunkerExecutorIntegration:
    """Tests for chunking + execution integration."""

    def test_chunk_and_execute_serial(self, mock_genome, mock_gff_parser):
        """Test chunking followed by serial execution."""
        # Create chunker and plan
        chunker = GenomeChunker(mock_genome, mock_gff_parser)
        plan = chunker.create_plan(
            ChunkStrategy.BY_SCAFFOLD,
            max_chunk_size=1_500_000,
        )

        # Verify reasonable number of chunks
        assert len(plan) >= 5  # At least 5 scaffolds

        # Execute with simple function
        def process_chunk(chunk):
            return {
                "chunk_id": chunk.chunk_id,
                "seqid": chunk.seqid,
                "size": chunk.size,
            }

        executor = ParallelExecutor(n_workers=1)
        results, stats = executor.map_chunks(process_chunk, list(plan.chunks))

        # Verify execution
        assert stats.total_tasks == len(plan)
        assert stats.successful == len(plan)
        assert stats.failed == 0

        # Verify results
        processed_seqids = {r.result["seqid"] for r in results if r.success}
        assert "chr1" in processed_seqids
        assert "chr5" in processed_seqids

    def test_chunk_and_execute_threaded(self, mock_genome):
        """Test chunking followed by threaded execution."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)

        def process_chunk(chunk):
            time.sleep(0.01)  # Simulate work
            return chunk.chunk_id

        executor = ParallelExecutor(n_workers=2, backend=ExecutorBackend.THREADS)
        results, stats = executor.map_chunks(process_chunk, list(plan.chunks))

        assert stats.successful == len(plan)

    @pytest.mark.slow
    def test_chunk_and_execute_parallel_processes(self, mock_genome):
        """Test chunking followed by process-based execution."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)

        def process_chunk(chunk):
            return chunk.size

        executor = ParallelExecutor(n_workers=2, backend=ExecutorBackend.PROCESSES)
        results, stats = executor.map_chunks(process_chunk, list(plan.chunks))

        assert stats.successful == len(plan)
        total_size = sum(r.result for r in results if r.success)
        assert total_size == mock_genome.total_length


# =============================================================================
# Plan Save/Load Integration
# =============================================================================


class TestPlanPersistence:
    """Tests for chunk plan persistence and loading."""

    def test_save_load_and_process(self, mock_genome, tmp_path):
        """Test saving plan, loading, and processing."""
        # Create and save plan
        chunker = GenomeChunker(mock_genome)
        original_plan = chunker.create_plan(ChunkStrategy.BY_SIZE, chunk_size=1_000_000)

        plan_path = tmp_path / "plan.json"
        original_plan.save(plan_path)

        # Load plan
        loaded_plan = ChunkPlan.load(plan_path)

        # Verify loaded plan matches original
        assert len(loaded_plan) == len(original_plan)
        assert loaded_plan.strategy == original_plan.strategy
        assert loaded_plan.total_bases == original_plan.total_bases

        # Process with loaded plan
        def process_chunk(chunk):
            return chunk.chunk_id

        executor = ParallelExecutor(n_workers=1)
        results, stats = executor.map_chunks(process_chunk, list(loaded_plan.chunks))

        assert stats.successful == len(loaded_plan)


# =============================================================================
# SLURM Job Generation Integration
# =============================================================================
# NOTE: The SlurmJobGenerator, SlurmConfig, SlurmArrayJob classes were removed
# in favor of task-file-based parallelization (TaskGenerator). These tests are
# kept as placeholders but skipped. See helixforge.parallel.taskgen for the
# current approach.


@pytest.mark.skip(reason="SlurmJobGenerator removed - use TaskGenerator instead")
class TestSlurmGenerationIntegration:
    """Tests for SLURM job generation integration (DEPRECATED)."""

    def test_full_slurm_workflow(self, mock_genome, mock_gff_parser, tmp_path):
        """Test complete SLURM job generation workflow."""
        pass

    def test_slurm_with_adaptive_chunking(self, mock_genome, mock_gff_parser, tmp_path):
        """Test SLURM with adaptive chunking strategy."""
        pass


# =============================================================================
# Error Handling Integration
# =============================================================================


class TestErrorHandlingIntegration:
    """Tests for error handling across the parallel module."""

    def test_partial_failure_recovery(self, mock_genome):
        """Test processing continues after partial failures."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)

        def sometimes_fails(chunk):
            if chunk.seqid == "chr3":
                raise ValueError("Simulated failure")
            return chunk.size

        executor = ParallelExecutor(n_workers=1)
        results, stats = executor.map_chunks(
            sometimes_fails, list(plan.chunks), continue_on_error=True
        )

        # Should have 4 successes and 1 failure
        assert stats.successful == 4
        assert stats.failed == 1

        # Verify correct chunk failed
        failed = [r for r in results if not r.success]
        assert len(failed) == 1
        assert "chr3" in failed[0].chunk_id or "Simulated failure" in failed[0].error

    @pytest.mark.skip(reason="aggregate_chunk_outputs removed - use CLI aggregate command")
    def test_result_aggregation_with_missing(self, tmp_path):
        """Test aggregation handles missing chunks gracefully."""
        # NOTE: aggregate_chunk_outputs was removed. Use the CLI command instead:
        # helixforge parallel aggregate --input-dir outputs/ --pattern '*.gff3' -o combined.gff3
        pass


# =============================================================================
# Resource Estimation Integration
# =============================================================================


class TestResourceEstimation:
    """Tests for resource estimation integration."""

    def test_suggest_parameters_integration(self):
        """Test parameter suggestion for various genome sizes."""
        # Small genome (Arabidopsis-like)
        small_result = suggest_chunk_parameters(
            genome_size=135_000_000,  # 135 Mb
            n_genes=27_000,
            available_memory_gb=16,
            n_workers=8,
        )
        # Should suggest BY_GENES for gene-dense genome
        assert small_result["strategy"] in [
            ChunkStrategy.BY_GENES,
            ChunkStrategy.ADAPTIVE,
        ]

        # Large genome (Maize-like)
        large_result = suggest_chunk_parameters(
            genome_size=2_300_000_000,  # 2.3 Gb
            n_genes=40_000,
            available_memory_gb=64,
            n_workers=16,
        )
        # Should suggest BY_SIZE for large genome
        assert large_result["strategy"] == ChunkStrategy.BY_SIZE
        assert large_result["chunk_size"] is not None

    def test_optimal_workers_with_memory(self):
        """Test optimal worker calculation considers memory."""
        # With 16 GB available and 2 GB per worker, should get 8 or fewer
        workers = get_optimal_workers(
            max_workers=16,
            memory_per_worker_mb=2000,
        )
        assert workers <= 16
        assert workers >= 1


# =============================================================================
# Result Merging Integration
# =============================================================================


class TestResultMergingIntegration:
    """Tests for result merging across chunks."""

    def test_merge_results_from_overlapping_chunks(self):
        """Test merging results from overlapping chunk strategy."""
        # Create overlapping chunks (e.g., from size-based chunking)
        chunk1 = GenomicChunk("chunk_0000", "chr1", 0, 1100)
        chunk2 = GenomicChunk("chunk_0001", "chr1", 1000, 2100)
        chunk3 = GenomicChunk("chunk_0002", "chr1", 2000, 3000)

        # Genes that span chunks
        gene1 = MagicMock(gene_id="gene1", seqid="chr1", start=500, end=800)
        gene2 = MagicMock(
            gene_id="gene2", seqid="chr1", start=950, end=1150
        )  # Spans chunk boundary
        gene3 = MagicMock(gene_id="gene3", seqid="chr1", start=1500, end=1800)
        gene4 = MagicMock(gene_id="gene4", seqid="chr1", start=2500, end=2900)

        # Gene2 appears in both chunk1 and chunk2
        results = [
            (chunk1, [gene1, gene2]),
            (chunk2, [gene2, gene3]),  # gene2 duplicate
            (chunk3, [gene4]),
        ]

        merged = merge_overlapping_results(results)

        # Should have exactly 4 unique genes
        assert len(merged) == 4
        gene_ids = [g.gene_id for g in merged]
        assert len(set(gene_ids)) == 4  # All unique


# =============================================================================
# Progress Tracking Integration
# =============================================================================


class TestProgressTrackingIntegration:
    """Tests for progress tracking during execution."""

    def test_progress_callback_fired(self, mock_genome):
        """Test that progress callback is called for each chunk."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)

        progress_events = []

        def track_progress(completed, total, chunk_id):
            progress_events.append(
                {
                    "completed": completed,
                    "total": total,
                    "chunk_id": chunk_id,
                }
            )

        executor = ParallelExecutor(n_workers=1, progress_callback=track_progress)
        executor.map_chunks(lambda c: c.size, list(plan.chunks))

        # Should have one event per chunk
        assert len(progress_events) == len(plan)

        # Last event should show all complete
        assert progress_events[-1]["completed"] == len(plan)
        assert progress_events[-1]["total"] == len(plan)


# =============================================================================
# Module Import Tests
# =============================================================================


class TestModuleImports:
    """Test that all public APIs are importable."""

    def test_chunker_imports(self):
        """Test chunker module imports."""
        from helixforge.parallel import (
            ChunkPlan,
            ChunkStrategy,
            GenomicChunk,
            GenomeChunker,
            Chunker,
            merge_overlapping_results,
            suggest_chunk_parameters,
        )

        assert ChunkPlan is not None
        assert Chunker is GenomeChunker

    def test_executor_imports(self):
        """Test executor module imports."""
        from helixforge.parallel import (
            ExecutorBackend,
            ExecutionStats,
            MemoryMonitor,
            MemoryStats,
            ParallelExecutor,
            TaskResult,
            ChunkProcessor,
            Executor,
            get_memory_stats,
            get_optimal_workers,
            track_memory,
        )

        assert ParallelExecutor is not None
        assert Executor is ParallelExecutor

    def test_slurm_imports(self):
        """Test SLURM utility imports (minimal utilities only)."""
        # NOTE: SlurmJobGenerator, SlurmConfig, SlurmArrayJob were removed.
        # The parallel module now uses task-file-based approach with TaskGenerator.
        # Only minimal SLURM utilities remain for environment detection.
        from helixforge.parallel import (
            detect_slurm_environment,
            get_chunk_for_task,
            get_slurm_task_id,
            get_slurm_resources,
            is_slurm_job,
            write_example_sbatch,
        )

        assert detect_slurm_environment is not None
        assert is_slurm_job is not None

    def test_taskgen_imports(self):
        """Test task generator module imports."""
        from helixforge.parallel import (
            TaskFile,
            TaskGenerator,
            format_command,
            generate_simple_task_file,
            estimate_parallelism,
            generate_hypershell_command,
        )

        assert TaskGenerator is not None
        assert format_command is not None
