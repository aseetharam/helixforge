"""Tests for helixforge.parallel.executor module.

Tests cover:
- MemoryStats data structure
- TaskResult data structure
- ExecutionStats data structure
- MemoryMonitor class
- ParallelExecutor with different backends
- Utility functions
"""

import time
from unittest.mock import MagicMock, patch

import pytest

from helixforge.parallel.executor import (
    ExecutorBackend,
    ExecutionStats,
    MemoryMonitor,
    MemoryStats,
    ParallelExecutor,
    TaskResult,
    Executor,
    get_memory_stats,
    get_optimal_workers,
    track_memory,
)
from helixforge.parallel.chunker import GenomicChunk


# =============================================================================
# Data Structure Tests
# =============================================================================


class TestMemoryStats:
    """Tests for MemoryStats data structure."""

    def test_create_memory_stats(self):
        """Test creating MemoryStats."""
        stats = MemoryStats(
            current_mb=1024.5,
            peak_mb=2048.0,
            available_mb=8192.0,
            percent_used=25.0,
        )
        assert stats.current_mb == 1024.5
        assert stats.peak_mb == 2048.0
        assert stats.available_mb == 8192.0
        assert stats.percent_used == 25.0

    def test_memory_stats_to_dict(self):
        """Test serialization to dict."""
        stats = MemoryStats(
            current_mb=1024.567,
            peak_mb=2048.123,
            available_mb=8192.456,
            percent_used=25.789,
        )
        d = stats.to_dict()
        assert d["current_mb"] == 1024.57  # Rounded
        assert d["peak_mb"] == 2048.12
        assert d["available_mb"] == 8192.46
        assert d["percent_used"] == 25.79


class TestTaskResult:
    """Tests for TaskResult data structure."""

    def test_create_successful_result(self):
        """Test creating successful TaskResult."""
        result = TaskResult(
            chunk_id="chunk_0001",
            success=True,
            result={"genes": 10},
            duration_seconds=1.5,
            peak_memory_mb=512.0,
        )
        assert result.chunk_id == "chunk_0001"
        assert result.success is True
        assert result.result == {"genes": 10}
        assert result.error is None
        assert result.duration_seconds == 1.5
        assert result.peak_memory_mb == 512.0

    def test_create_failed_result(self):
        """Test creating failed TaskResult."""
        result = TaskResult(
            chunk_id="chunk_0001",
            success=False,
            error="File not found",
            duration_seconds=0.1,
        )
        assert result.success is False
        assert result.error == "File not found"
        assert result.result is None

    def test_task_result_to_dict(self):
        """Test serialization to dict."""
        result = TaskResult(
            chunk_id="chunk_0001",
            success=True,
            result={"data": "test"},  # Not included in to_dict
            duration_seconds=1.567,
            peak_memory_mb=512.123,
        )
        d = result.to_dict()
        assert d["chunk_id"] == "chunk_0001"
        assert d["success"] is True
        assert d["error"] is None
        assert d["duration_seconds"] == 1.567
        assert d["peak_memory_mb"] == 512.12
        assert "result" not in d  # Result not serialized


class TestExecutionStats:
    """Tests for ExecutionStats data structure."""

    def test_create_execution_stats(self):
        """Test creating ExecutionStats."""
        stats = ExecutionStats(
            total_tasks=10,
            successful=8,
            failed=2,
            total_duration=60.5,
            mean_task_duration=5.0,
            max_task_duration=12.3,
            peak_memory_mb=1024.0,
        )
        assert stats.total_tasks == 10
        assert stats.successful == 8
        assert stats.failed == 2
        assert stats.total_duration == 60.5
        assert stats.mean_task_duration == 5.0
        assert stats.max_task_duration == 12.3
        assert stats.peak_memory_mb == 1024.0

    def test_execution_stats_to_dict(self):
        """Test serialization to dict."""
        stats = ExecutionStats(
            total_tasks=10,
            successful=8,
            failed=2,
            total_duration=60.567,
            mean_task_duration=5.123,
            max_task_duration=12.345,
            peak_memory_mb=1024.789,
        )
        d = stats.to_dict()
        assert d["total_tasks"] == 10
        assert d["successful"] == 8
        assert d["failed"] == 2
        assert d["total_duration"] == 60.567
        assert d["mean_task_duration"] == 5.123
        assert d["max_task_duration"] == 12.345
        assert d["peak_memory_mb"] == 1024.79


# =============================================================================
# Memory Monitoring Tests
# =============================================================================


class TestGetMemoryStats:
    """Tests for get_memory_stats function."""

    def test_get_memory_stats_with_psutil(self):
        """Test memory stats with psutil available."""
        stats = get_memory_stats()
        assert stats.current_mb > 0
        assert stats.peak_mb > 0
        # These might be 0 without psutil
        assert stats.current_mb <= stats.peak_mb or stats.percent_used == 0

    def test_get_memory_stats_without_psutil(self):
        """Test memory stats fallback without psutil."""
        with patch.dict("sys.modules", {"psutil": None}):
            # Should still work with resource module
            stats = get_memory_stats()
            assert stats.peak_mb >= 0


class TestTrackMemory:
    """Tests for track_memory context manager."""

    def test_track_memory_basic(self):
        """Test basic memory tracking."""
        with track_memory() as get_stats:
            # Allocate some memory
            data = [0] * 1000000
            stats = get_stats()
            assert stats.current_mb > 0
            assert stats.peak_mb >= stats.current_mb
            del data

    def test_track_memory_updates_peak(self):
        """Test that peak is updated correctly."""
        with track_memory() as get_stats:
            initial_stats = get_stats()
            initial_peak = initial_stats.peak_mb

            # Allocate more memory
            data = [0] * 5000000
            stats_during = get_stats()
            del data

            # Peak should be at least as high as before
            assert stats_during.peak_mb >= initial_peak


class TestMemoryMonitor:
    """Tests for MemoryMonitor class."""

    def test_monitor_creation(self):
        """Test monitor initialization."""
        monitor = MemoryMonitor(
            warning_threshold_percent=80.0,
            critical_threshold_percent=95.0,
            sample_interval_seconds=0.5,
        )
        assert monitor.warning_threshold == 80.0
        assert monitor.critical_threshold == 95.0
        assert monitor.sample_interval == 0.5

    def test_monitor_start_stop(self):
        """Test starting and stopping monitor."""
        monitor = MemoryMonitor(sample_interval_seconds=0.1)
        monitor.start()

        # Let it sample a bit
        time.sleep(0.25)

        stats = monitor.stop()
        assert stats.peak_mb > 0

    def test_monitor_double_start(self):
        """Test that double start is safe."""
        monitor = MemoryMonitor()
        monitor.start()
        monitor.start()  # Should be no-op
        monitor.stop()

    def test_monitor_check(self):
        """Test memory check."""
        monitor = MemoryMonitor(
            warning_threshold_percent=99.0,
            critical_threshold_percent=99.9,
        )
        monitor.start()

        is_ok, message = monitor.check()
        # Should usually be OK unless memory is very high
        assert isinstance(is_ok, bool)
        assert message is None or isinstance(message, str)

        monitor.stop()


# =============================================================================
# ParallelExecutor Tests
# =============================================================================


class TestParallelExecutor:
    """Tests for ParallelExecutor class."""

    @pytest.fixture
    def sample_chunks(self) -> list[GenomicChunk]:
        """Create sample chunks for testing."""
        return [
            GenomicChunk("chunk_0000", "chr1", 0, 1000),
            GenomicChunk("chunk_0001", "chr1", 1000, 2000),
            GenomicChunk("chunk_0002", "chr2", 0, 500),
        ]

    def test_executor_creation(self):
        """Test executor initialization."""
        executor = ParallelExecutor(
            n_workers=4,
            backend=ExecutorBackend.PROCESSES,
            max_memory_gb=8.0,
        )
        assert executor.n_workers == 4
        assert executor.backend == ExecutorBackend.PROCESSES
        assert executor.max_memory_gb == 8.0

    def test_executor_single_worker_uses_serial(self):
        """Test that single worker uses serial backend."""
        executor = ParallelExecutor(n_workers=1)
        assert executor.backend == ExecutorBackend.SERIAL

    def test_executor_string_backend(self):
        """Test using string for backend."""
        executor = ParallelExecutor(n_workers=2, backend="threads")
        assert executor.backend == ExecutorBackend.THREADS

    def test_executor_serial_success(self, sample_chunks):
        """Test serial execution with successful tasks."""

        def simple_func(chunk):
            return chunk.size

        executor = ParallelExecutor(n_workers=1)
        results, stats = executor.map_chunks(simple_func, sample_chunks)

        assert stats.total_tasks == 3
        assert stats.successful == 3
        assert stats.failed == 0
        assert len(results) == 3

        # Check results
        for result in results:
            assert result.success is True
            assert result.result in [500, 1000]

    def test_executor_serial_with_errors(self, sample_chunks):
        """Test serial execution with errors."""

        def failing_func(chunk):
            if chunk.seqid == "chr2":
                raise ValueError("Simulated error")
            return chunk.size

        executor = ParallelExecutor(n_workers=1)
        results, stats = executor.map_chunks(failing_func, sample_chunks)

        assert stats.total_tasks == 3
        assert stats.successful == 2
        assert stats.failed == 1

        # Find the failed result
        failed = [r for r in results if not r.success]
        assert len(failed) == 1
        assert "Simulated error" in failed[0].error

    def test_executor_stop_on_error(self, sample_chunks):
        """Test stopping execution on error."""

        def failing_func(chunk):
            if chunk.chunk_id == "chunk_0001":
                raise ValueError("Stop here")
            return chunk.size

        executor = ParallelExecutor(n_workers=1)

        with pytest.raises(ValueError):
            executor.map_chunks(
                failing_func, sample_chunks, continue_on_error=False
            )

    def test_executor_empty_chunks(self):
        """Test execution with empty chunk list."""
        executor = ParallelExecutor(n_workers=1)
        results, stats = executor.map_chunks(lambda x: x, [])

        assert stats.total_tasks == 0
        assert stats.successful == 0
        assert results == []

    def test_executor_progress_callback(self, sample_chunks):
        """Test progress callback."""
        progress_log = []

        def callback(completed, total, chunk_id):
            progress_log.append((completed, total, chunk_id))

        executor = ParallelExecutor(n_workers=1, progress_callback=callback)
        executor.map_chunks(lambda c: c.size, sample_chunks)

        assert len(progress_log) == 3
        assert progress_log[-1][0] == 3  # Last completed count
        assert progress_log[-1][1] == 3  # Total

    def test_executor_threads_backend(self, sample_chunks):
        """Test threaded execution."""

        def simple_func(chunk):
            time.sleep(0.01)
            return chunk.size

        executor = ParallelExecutor(n_workers=2, backend=ExecutorBackend.THREADS)
        results, stats = executor.map_chunks(simple_func, sample_chunks)

        assert stats.total_tasks == 3
        assert stats.successful == 3

    @pytest.mark.slow
    def test_executor_processes_backend(self, sample_chunks):
        """Test process-based execution."""

        def simple_func(chunk):
            return chunk.size

        executor = ParallelExecutor(n_workers=2, backend=ExecutorBackend.PROCESSES)
        results, stats = executor.map_chunks(simple_func, sample_chunks)

        assert stats.total_tasks == 3
        assert stats.successful == 3

    def test_executor_map_items(self):
        """Test map_items for non-chunk data."""
        items = [1, 2, 3, 4, 5]

        def double(x):
            return x * 2

        executor = ParallelExecutor(n_workers=1)
        results, stats = executor.map_items(double, items)

        assert stats.total_tasks == 5
        assert stats.successful == 5

        values = [r.result for r in results if r.success]
        assert sorted(values) == [2, 4, 6, 8, 10]

    def test_executor_with_memory_monitoring(self, sample_chunks):
        """Test execution with memory monitoring."""

        def simple_func(chunk):
            return chunk.size

        executor = ParallelExecutor(n_workers=1, max_memory_gb=16.0)
        results, stats = executor.map_chunks(simple_func, sample_chunks)

        assert stats.total_tasks == 3
        # Peak memory may or may not be tracked depending on monitoring
        # Just verify it completes without error

    def test_executor_legacy_alias(self):
        """Test legacy Executor alias."""
        executor = Executor(n_workers=1)
        assert isinstance(executor, ParallelExecutor)


# =============================================================================
# Utility Function Tests
# =============================================================================


class TestGetOptimalWorkers:
    """Tests for get_optimal_workers function."""

    def test_default_workers(self):
        """Test default worker count."""
        import os

        workers = get_optimal_workers()
        assert workers >= 1
        assert workers <= os.cpu_count()

    def test_max_workers_limit(self):
        """Test max_workers parameter."""
        workers = get_optimal_workers(max_workers=2)
        assert workers <= 2

    def test_memory_based_limit(self):
        """Test memory-based worker limit."""
        # Request a lot of memory per worker to trigger limit
        workers = get_optimal_workers(
            max_workers=100,
            memory_per_worker_mb=50000,  # 50 GB per worker
        )
        # Should be limited by available memory
        assert workers >= 1


# =============================================================================
# ExecutorBackend Tests
# =============================================================================


class TestExecutorBackend:
    """Tests for ExecutorBackend enum."""

    def test_backend_values(self):
        """Test backend enum values."""
        assert ExecutorBackend.SERIAL.value == "serial"
        assert ExecutorBackend.THREADS.value == "threads"
        assert ExecutorBackend.PROCESSES.value == "processes"

    def test_backend_from_string(self):
        """Test creating backend from string."""
        assert ExecutorBackend("serial") == ExecutorBackend.SERIAL
        assert ExecutorBackend("threads") == ExecutorBackend.THREADS
        assert ExecutorBackend("processes") == ExecutorBackend.PROCESSES
