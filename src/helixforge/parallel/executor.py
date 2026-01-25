"""Local parallel execution using multiprocessing.

This module provides tools for executing HelixForge operations in
parallel using Python's multiprocessing or threading.

Features:
    - Multiple execution backends (serial, threads, processes)
    - Progress tracking with rich
    - Error handling and recovery
    - Memory monitoring

Example:
    >>> from helixforge.parallel.executor import ParallelExecutor, ExecutorBackend
    >>> executor = ParallelExecutor(n_workers=8)
    >>> results, stats = executor.map_chunks(process_chunk, chunks)
"""

from __future__ import annotations

import logging
import os
import resource
import time
import threading
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    Future,
    as_completed,
)
from contextlib import contextmanager
from enum import Enum
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Generator,
    Iterator,
    TypeVar,
)

import attrs

if TYPE_CHECKING:
    from helixforge.parallel.chunker import GenomicChunk, ChunkPlan
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser

logger = logging.getLogger(__name__)

T = TypeVar("T")
R = TypeVar("R")


# =============================================================================
# Enums
# =============================================================================


class ExecutorBackend(Enum):
    """Available execution backends."""

    SERIAL = "serial"
    THREADS = "threads"
    PROCESSES = "processes"


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class MemoryStats:
    """Memory usage statistics."""

    current_mb: float
    peak_mb: float
    available_mb: float
    percent_used: float

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "current_mb": round(self.current_mb, 2),
            "peak_mb": round(self.peak_mb, 2),
            "available_mb": round(self.available_mb, 2),
            "percent_used": round(self.percent_used, 2),
        }


@attrs.define(slots=True)
class TaskResult:
    """Result from a parallel task."""

    chunk_id: str
    success: bool
    result: Any | None = None
    error: str | None = None
    duration_seconds: float = 0.0
    peak_memory_mb: float | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "chunk_id": self.chunk_id,
            "success": self.success,
            "error": self.error,
            "duration_seconds": round(self.duration_seconds, 3),
            "peak_memory_mb": (
                round(self.peak_memory_mb, 2) if self.peak_memory_mb else None
            ),
        }


@attrs.define(slots=True)
class ExecutionStats:
    """Statistics from parallel execution."""

    total_tasks: int
    successful: int
    failed: int
    total_duration: float
    mean_task_duration: float
    max_task_duration: float
    peak_memory_mb: float | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "total_tasks": self.total_tasks,
            "successful": self.successful,
            "failed": self.failed,
            "total_duration": round(self.total_duration, 3),
            "mean_task_duration": round(self.mean_task_duration, 3),
            "max_task_duration": round(self.max_task_duration, 3),
            "peak_memory_mb": (
                round(self.peak_memory_mb, 2) if self.peak_memory_mb else None
            ),
        }


# =============================================================================
# Memory Monitoring
# =============================================================================


def get_memory_stats() -> MemoryStats:
    """Get current memory usage statistics.

    Returns:
        MemoryStats with current memory information.
    """
    try:
        import psutil

        process = psutil.Process()
        mem_info = process.memory_info()
        virtual = psutil.virtual_memory()

        # Get peak memory from resource module (more reliable on Linux)
        try:
            peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            # On Linux, ru_maxrss is in KB; on macOS it's in bytes
            import platform

            if platform.system() == "Darwin":
                peak_mb = peak_kb / 1024 / 1024
            else:
                peak_mb = peak_kb / 1024
        except Exception:
            peak_mb = mem_info.rss / 1024 / 1024

        return MemoryStats(
            current_mb=mem_info.rss / 1024 / 1024,
            peak_mb=peak_mb,
            available_mb=virtual.available / 1024 / 1024,
            percent_used=virtual.percent,
        )
    except ImportError:
        # psutil not available, use basic resource module
        usage = resource.getrusage(resource.RUSAGE_SELF)
        peak_kb = usage.ru_maxrss
        import platform

        if platform.system() == "Darwin":
            peak_mb = peak_kb / 1024 / 1024
        else:
            peak_mb = peak_kb / 1024

        return MemoryStats(
            current_mb=peak_mb,  # Approximation
            peak_mb=peak_mb,
            available_mb=0.0,
            percent_used=0.0,
        )


@contextmanager
def track_memory() -> Generator[Callable[[], MemoryStats], None, None]:
    """Context manager to track memory usage.

    Yields function to get current stats.
    Records peak usage automatically.

    Example:
        >>> with track_memory() as get_stats:
        ...     # do work
        ...     stats = get_stats()
        ...     print(f"Peak memory: {stats.peak_mb:.1f} MB")
    """
    start_stats = get_memory_stats()
    peak_mb = start_stats.current_mb

    def get_current_stats() -> MemoryStats:
        nonlocal peak_mb
        stats = get_memory_stats()
        peak_mb = max(peak_mb, stats.peak_mb)
        return MemoryStats(
            current_mb=stats.current_mb,
            peak_mb=peak_mb,
            available_mb=stats.available_mb,
            percent_used=stats.percent_used,
        )

    yield get_current_stats


class MemoryMonitor:
    """Monitor memory usage during parallel execution.

    Features:
    - Periodic sampling
    - Peak tracking
    - Warning thresholds
    - Automatic logging

    Example:
        >>> monitor = MemoryMonitor(warning_threshold_percent=80.0)
        >>> monitor.start()
        >>> # do work
        >>> stats = monitor.stop()
        >>> print(f"Peak: {stats.peak_mb:.1f} MB")
    """

    def __init__(
        self,
        warning_threshold_percent: float = 80.0,
        critical_threshold_percent: float = 95.0,
        sample_interval_seconds: float = 1.0,
    ) -> None:
        """Initialize the monitor.

        Args:
            warning_threshold_percent: Memory usage percentage to trigger warning.
            critical_threshold_percent: Memory usage percentage to trigger critical.
            sample_interval_seconds: Seconds between samples.
        """
        self.warning_threshold = warning_threshold_percent
        self.critical_threshold = critical_threshold_percent
        self.sample_interval = sample_interval_seconds

        self._running = False
        self._thread: threading.Thread | None = None
        self._peak_mb = 0.0
        self._samples: list[MemoryStats] = []
        self._lock = threading.Lock()

    def start(self) -> None:
        """Start background monitoring thread."""
        if self._running:
            return

        self._running = True
        self._peak_mb = 0.0
        self._samples = []
        self._thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self._thread.start()
        logger.debug("Memory monitor started")

    def stop(self) -> MemoryStats:
        """Stop monitoring and return final stats.

        Returns:
            Final MemoryStats with peak usage.
        """
        self._running = False
        if self._thread is not None:
            self._thread.join(timeout=2.0)
            self._thread = None

        final_stats = get_memory_stats()
        with self._lock:
            self._peak_mb = max(self._peak_mb, final_stats.peak_mb)
            peak = self._peak_mb

        logger.debug(f"Memory monitor stopped. Peak: {peak:.1f} MB")
        return MemoryStats(
            current_mb=final_stats.current_mb,
            peak_mb=peak,
            available_mb=final_stats.available_mb,
            percent_used=final_stats.percent_used,
        )

    def check(self) -> tuple[bool, str | None]:
        """Check current memory status.

        Returns:
            (is_ok, warning_message) tuple.
        """
        stats = get_memory_stats()
        with self._lock:
            self._peak_mb = max(self._peak_mb, stats.current_mb)

        if stats.percent_used >= self.critical_threshold:
            return False, f"CRITICAL: Memory at {stats.percent_used:.1f}%"
        elif stats.percent_used >= self.warning_threshold:
            return True, f"WARNING: Memory at {stats.percent_used:.1f}%"
        return True, None

    def _monitor_loop(self) -> None:
        """Background monitoring loop."""
        while self._running:
            try:
                stats = get_memory_stats()
                with self._lock:
                    self._peak_mb = max(self._peak_mb, stats.current_mb)
                    self._samples.append(stats)

                # Log warnings
                if stats.percent_used >= self.critical_threshold:
                    logger.warning(
                        f"CRITICAL: Memory usage at {stats.percent_used:.1f}%"
                    )
                elif stats.percent_used >= self.warning_threshold:
                    logger.warning(f"Memory usage at {stats.percent_used:.1f}%")

            except Exception as e:
                logger.debug(f"Memory monitoring error: {e}")

            time.sleep(self.sample_interval)


# =============================================================================
# Parallel Executor
# =============================================================================


class ParallelExecutor:
    """Execute tasks across genomic chunks in parallel.

    Features:
    - Multiple backends (serial, threads, processes)
    - Progress tracking with callbacks
    - Memory monitoring
    - Graceful error handling (continue on failure)
    - Result aggregation

    Example:
        >>> executor = ParallelExecutor(n_workers=4, backend="processes")
        >>> results, stats = executor.map_chunks(process_func, chunks)
        >>> print(f"Processed {stats.successful}/{stats.total_tasks} chunks")
    """

    def __init__(
        self,
        n_workers: int = 1,
        backend: ExecutorBackend | str = ExecutorBackend.PROCESSES,
        max_memory_gb: float | None = None,
        progress_callback: Callable[[int, int, str], None] | None = None,
    ) -> None:
        """Initialize executor.

        Args:
            n_workers: Number of parallel workers (1 = serial).
            backend: Execution backend.
            max_memory_gb: Memory limit (None = no limit).
            progress_callback: Called with (completed, total, chunk_id).
        """
        self.n_workers = max(1, n_workers)
        self.backend = (
            ExecutorBackend(backend) if isinstance(backend, str) else backend
        )
        self.max_memory_gb = max_memory_gb
        self.progress_callback = progress_callback

        # Auto-select serial if n_workers=1
        if self.n_workers == 1:
            self.backend = ExecutorBackend.SERIAL

        self._memory_monitor: MemoryMonitor | None = None

    def map_chunks(
        self,
        func: Callable[[GenomicChunk], R],
        chunks: list[GenomicChunk],
        continue_on_error: bool = True,
    ) -> tuple[list[TaskResult], ExecutionStats]:
        """Apply function to each chunk.

        Args:
            func: Function taking GenomicChunk, returning result.
            chunks: List of chunks to process.
            continue_on_error: If True, continue processing after failures.

        Returns:
            Tuple of (results_list, execution_stats).
        """
        if not chunks:
            return [], ExecutionStats(
                total_tasks=0,
                successful=0,
                failed=0,
                total_duration=0.0,
                mean_task_duration=0.0,
                max_task_duration=0.0,
            )

        logger.info(
            f"Processing {len(chunks)} chunks with {self.n_workers} workers "
            f"(backend={self.backend.value})"
        )

        start_time = time.time()

        # Start memory monitoring if limit set
        if self.max_memory_gb:
            self._memory_monitor = MemoryMonitor()
            self._memory_monitor.start()

        try:
            if self.backend == ExecutorBackend.SERIAL:
                results = self._execute_serial(func, chunks, continue_on_error)
            elif self.backend == ExecutorBackend.THREADS:
                results = self._execute_threaded(func, chunks, continue_on_error)
            else:
                results = self._execute_parallel(func, chunks, continue_on_error)
        finally:
            peak_memory = None
            if self._memory_monitor:
                final_stats = self._memory_monitor.stop()
                peak_memory = final_stats.peak_mb
                self._memory_monitor = None

        # Compute statistics
        total_duration = time.time() - start_time
        successful = sum(1 for r in results if r.success)
        failed = len(results) - successful
        durations = [r.duration_seconds for r in results]

        stats = ExecutionStats(
            total_tasks=len(results),
            successful=successful,
            failed=failed,
            total_duration=total_duration,
            mean_task_duration=sum(durations) / len(durations) if durations else 0,
            max_task_duration=max(durations) if durations else 0,
            peak_memory_mb=peak_memory,
        )

        logger.info(
            f"Completed: {successful}/{len(chunks)} chunks, "
            f"duration={total_duration:.1f}s"
        )

        return results, stats

    def map_items(
        self,
        func: Callable[[T], R],
        items: list[T],
        desc: str = "Processing",
        continue_on_error: bool = True,
    ) -> tuple[list[TaskResult], ExecutionStats]:
        """Apply function to items with automatic batching.

        For processing lists of genes, junctions, etc.

        Args:
            func: Function to apply to each item.
            items: Items to process.
            desc: Description for progress.
            continue_on_error: Continue after errors.

        Returns:
            Tuple of (results, stats).
        """
        # Wrap items to look like chunks
        wrapped_items = [
            _ItemWrapper(f"item_{i:06d}", item) for i, item in enumerate(items)
        ]

        def wrapped_func(item_wrapper: _ItemWrapper) -> R:
            return func(item_wrapper.item)

        return self.map_chunks(wrapped_func, wrapped_items, continue_on_error)

    def _execute_serial(
        self,
        func: Callable,
        chunks: list,
        continue_on_error: bool,
    ) -> list[TaskResult]:
        """Serial execution with progress tracking."""
        results = []
        total = len(chunks)

        for i, chunk in enumerate(chunks):
            chunk_id = getattr(chunk, "chunk_id", str(i))
            start_time = time.time()

            try:
                result = func(chunk)
                duration = time.time() - start_time
                results.append(
                    TaskResult(
                        chunk_id=chunk_id,
                        success=True,
                        result=result,
                        duration_seconds=duration,
                    )
                )
            except Exception as e:
                duration = time.time() - start_time
                results.append(
                    TaskResult(
                        chunk_id=chunk_id,
                        success=False,
                        error=str(e),
                        duration_seconds=duration,
                    )
                )
                if not continue_on_error:
                    logger.error(f"Task {chunk_id} failed: {e}")
                    raise

            if self.progress_callback:
                self.progress_callback(i + 1, total, chunk_id)

            # Check memory limit
            if self._memory_monitor and self.max_memory_gb:
                stats = get_memory_stats()
                if stats.current_mb > self.max_memory_gb * 1024:
                    logger.warning(
                        f"Memory limit exceeded: {stats.current_mb:.0f} MB > "
                        f"{self.max_memory_gb * 1024:.0f} MB"
                    )

        return results

    def _execute_parallel(
        self,
        func: Callable,
        chunks: list,
        continue_on_error: bool,
    ) -> list[TaskResult]:
        """Parallel execution using ProcessPoolExecutor."""
        results = []
        total = len(chunks)
        completed = 0

        # Create wrapper function that includes timing
        def timed_func(chunk):
            chunk_id = getattr(chunk, "chunk_id", "unknown")
            start_time = time.time()
            try:
                result = func(chunk)
                return TaskResult(
                    chunk_id=chunk_id,
                    success=True,
                    result=result,
                    duration_seconds=time.time() - start_time,
                )
            except Exception as e:
                return TaskResult(
                    chunk_id=chunk_id,
                    success=False,
                    error=str(e),
                    duration_seconds=time.time() - start_time,
                )

        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            futures: dict[Future, Any] = {}

            for chunk in chunks:
                future = executor.submit(timed_func, chunk)
                futures[future] = chunk

            for future in as_completed(futures):
                completed += 1
                try:
                    task_result = future.result()
                    results.append(task_result)

                    if self.progress_callback:
                        self.progress_callback(completed, total, task_result.chunk_id)

                    if not task_result.success and not continue_on_error:
                        logger.error(f"Task failed: {task_result.error}")
                        executor.shutdown(wait=False, cancel_futures=True)
                        raise RuntimeError(task_result.error)

                except Exception as e:
                    if not continue_on_error:
                        raise

        return results

    def _execute_threaded(
        self,
        func: Callable,
        chunks: list,
        continue_on_error: bool,
    ) -> list[TaskResult]:
        """Threaded execution (for I/O-bound tasks)."""
        results = []
        total = len(chunks)
        completed = 0
        results_lock = threading.Lock()

        def timed_func(chunk):
            chunk_id = getattr(chunk, "chunk_id", "unknown")
            start_time = time.time()
            try:
                result = func(chunk)
                return TaskResult(
                    chunk_id=chunk_id,
                    success=True,
                    result=result,
                    duration_seconds=time.time() - start_time,
                )
            except Exception as e:
                return TaskResult(
                    chunk_id=chunk_id,
                    success=False,
                    error=str(e),
                    duration_seconds=time.time() - start_time,
                )

        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures: dict[Future, Any] = {}

            for chunk in chunks:
                future = executor.submit(timed_func, chunk)
                futures[future] = chunk

            for future in as_completed(futures):
                completed += 1
                try:
                    task_result = future.result()
                    with results_lock:
                        results.append(task_result)

                    if self.progress_callback:
                        self.progress_callback(completed, total, task_result.chunk_id)

                    if not task_result.success and not continue_on_error:
                        raise RuntimeError(task_result.error)

                except Exception as e:
                    if not continue_on_error:
                        raise

        return results


@attrs.define
class _ItemWrapper:
    """Wrapper to give items a chunk_id for map_items."""

    chunk_id: str
    item: Any


# =============================================================================
# Chunk Processor
# =============================================================================


class ChunkProcessor:
    """High-level interface for processing genome in chunks.

    Combines chunking strategy with parallel execution.

    Example:
        >>> processor = ChunkProcessor(genome, gff_parser, n_workers=4)
        >>> results = processor.process(my_processing_func)
    """

    def __init__(
        self,
        genome: GenomeAccessor,
        gff_parser: GFF3Parser,
        n_workers: int = 1,
        chunk_strategy: str = "scaffold",
        chunk_size: int | None = None,
        memory_limit_gb: float | None = None,
    ) -> None:
        """Initialize the processor.

        Args:
            genome: GenomeAccessor for the genome.
            gff_parser: GFF3Parser for gene annotations.
            n_workers: Number of parallel workers.
            chunk_strategy: Chunking strategy name.
            chunk_size: Size parameter for chunking.
            memory_limit_gb: Memory limit in GB.
        """
        from helixforge.parallel.chunker import GenomeChunker, ChunkStrategy

        self.genome = genome
        self.gff_parser = gff_parser
        self.n_workers = n_workers
        self.memory_limit_gb = memory_limit_gb

        # Create chunker and plan
        self.chunker = GenomeChunker(genome, gff_parser)
        self.plan = self.chunker.create_plan(
            strategy=ChunkStrategy(chunk_strategy),
            chunk_size=chunk_size,
        )

        # Create executor
        self.executor = ParallelExecutor(
            n_workers=n_workers,
            backend=ExecutorBackend.PROCESSES if n_workers > 1 else ExecutorBackend.SERIAL,
            max_memory_gb=memory_limit_gb,
        )

    def process(
        self,
        processor_func: Callable[[GenomicChunk], list],
        aggregator_func: Callable[[list[list]], Any] | None = None,
        show_progress: bool = True,
    ) -> Any:
        """Process genome with given function.

        Args:
            processor_func: Function to apply to each chunk.
            aggregator_func: Optional function to combine results.
            show_progress: Show progress bar.

        Returns:
            Aggregated results.
        """
        # Set up progress callback
        progress_callback = None
        if show_progress:
            try:
                from rich.progress import Progress, SpinnerColumn, TextColumn

                progress = Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    transient=True,
                )
                task_id = progress.add_task("Processing...", total=len(self.plan))

                def callback(completed: int, total: int, chunk_id: str) -> None:
                    progress.update(
                        task_id,
                        completed=completed,
                        description=f"Processing {chunk_id}...",
                    )

                progress_callback = callback
                progress.start()
            except ImportError:
                pass

        try:
            results, stats = self.executor.map_chunks(
                processor_func,
                list(self.plan.chunks),
                continue_on_error=True,
            )
        finally:
            if show_progress and "progress" in locals():
                progress.stop()

        # Extract actual results
        chunk_results = [r.result for r in results if r.success and r.result is not None]

        # Aggregate if function provided
        if aggregator_func is not None:
            return aggregator_func(chunk_results)

        return chunk_results


# =============================================================================
# Utility Functions
# =============================================================================


def get_optimal_workers(
    max_workers: int | None = None,
    memory_per_worker_mb: int = 1000,
) -> int:
    """Determine optimal number of workers based on system resources.

    Args:
        max_workers: Maximum workers (defaults to CPU count).
        memory_per_worker_mb: Expected memory per worker in MB.

    Returns:
        Optimal number of workers.
    """
    cpu_count = os.cpu_count() or 1

    if max_workers is None:
        max_workers = cpu_count

    # Try to limit based on available memory
    try:
        stats = get_memory_stats()
        memory_limited = int(stats.available_mb / memory_per_worker_mb)
        max_workers = min(max_workers, max(1, memory_limited))
    except Exception:
        pass

    return min(max_workers, cpu_count)


def create_progress_bar(
    total: int,
    description: str = "Processing",
) -> Any:
    """Create rich progress bar for parallel execution.

    Args:
        total: Total number of items.
        description: Description text.

    Returns:
        Rich Progress object or None if rich not available.
    """
    try:
        from rich.progress import (
            Progress,
            SpinnerColumn,
            BarColumn,
            TextColumn,
            TimeRemainingColumn,
            MofNCompleteColumn,
        )

        return Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeRemainingColumn(),
        )
    except ImportError:
        return None


# =============================================================================
# Function Wrapper for Pickling
# =============================================================================


class _FunctionWrapper:
    """Wrapper to make instance methods picklable.

    Used when processing with ProcessPoolExecutor.
    """

    def __init__(self, obj: Any, method_name: str) -> None:
        """Initialize wrapper.

        Args:
            obj: Object instance.
            method_name: Name of method to wrap.
        """
        self.obj = obj
        self.method_name = method_name

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        """Call the wrapped method."""
        method = getattr(self.obj, self.method_name)
        return method(*args, **kwargs)


# =============================================================================
# Legacy Compatibility
# =============================================================================

# Alias for backward compatibility
Executor = ParallelExecutor
