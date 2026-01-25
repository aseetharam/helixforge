"""Logging configuration for HelixForge.

This module provides structured logging setup for HelixForge,
with support for console and file output.

Features:
    - Structured logging with rich formatting
    - File logging for debugging
    - Configurable verbosity levels
    - Progress logging for long operations

Example:
    >>> from helixforge.utils.logging import setup_logging, get_logger
    >>> setup_logging(verbosity=2)
    >>> logger = get_logger(__name__)
    >>> logger.info("Processing started")

TODO:
    - Implement structured logging
    - Add file handler configuration
    - Add progress logging
    - Support for JSON log format
    - Add timing utilities
"""

import logging
from pathlib import Path
from typing import Any

# =============================================================================
# Constants
# =============================================================================

# Default log format
DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Rich console format (when using rich handler)
RICH_FORMAT = "%(message)s"

# Log levels by verbosity
VERBOSITY_LEVELS = {
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


# =============================================================================
# Setup Functions
# =============================================================================


def setup_logging(
    verbosity: int = 1,
    log_file: Path | str | None = None,
    use_rich: bool = True,
) -> None:
    """Configure logging for HelixForge.

    Args:
        verbosity: Verbosity level (0=warning, 1=info, 2=debug).
        log_file: Optional file to log to.
        use_rich: Use rich for console output.
    """
    level = VERBOSITY_LEVELS.get(verbosity, logging.DEBUG)

    # Get root logger for helixforge
    logger = logging.getLogger("helixforge")
    logger.setLevel(level)

    # Remove existing handlers
    logger.handlers.clear()

    # Console handler
    if use_rich:
        try:
            from rich.logging import RichHandler

            console_handler = RichHandler(
                rich_tracebacks=True,
                show_time=False,
                show_path=False,
            )
            console_handler.setFormatter(logging.Formatter(RICH_FORMAT))
        except ImportError:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(logging.Formatter(DEFAULT_FORMAT))
    else:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter(DEFAULT_FORMAT))

    console_handler.setLevel(level)
    logger.addHandler(console_handler)

    # File handler
    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(DEFAULT_FORMAT))
        file_handler.setLevel(logging.DEBUG)  # Always log debug to file
        logger.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    """Get a logger for a module.

    Args:
        name: Module name (typically __name__).

    Returns:
        Logger instance.
    """
    return logging.getLogger(name)


# =============================================================================
# Progress Logging
# =============================================================================


class ProgressLogger:
    """Logger for long-running operations with progress tracking.

    Provides periodic updates for operations processing many items.

    Attributes:
        logger: The underlying logger.
        total: Total number of items.
        interval: Logging interval.

    Example:
        >>> progress = ProgressLogger(logger, total=1000)
        >>> for item in items:
        ...     process(item)
        ...     progress.update()
    """

    def __init__(
        self,
        logger: logging.Logger,
        total: int,
        interval: int = 100,
        description: str = "Processing",
    ) -> None:
        """Initialize progress logger.

        Args:
            logger: Logger to use.
            total: Total number of items.
            interval: Items between log messages.
            description: Description of the operation.
        """
        self.logger = logger
        self.total = total
        self.interval = interval
        self.description = description
        self.count = 0

    def update(self, n: int = 1) -> None:
        """Update progress counter.

        Args:
            n: Number of items completed.
        """
        self.count += n
        if self.count % self.interval == 0 or self.count == self.total:
            pct = 100 * self.count / self.total if self.total > 0 else 100
            self.logger.info(f"{self.description}: {self.count}/{self.total} ({pct:.1f}%)")

    def finish(self) -> None:
        """Mark progress as complete."""
        self.logger.info(f"{self.description}: Complete ({self.total} items)")


# =============================================================================
# Timing Utilities
# =============================================================================


class Timer:
    """Context manager for timing operations.

    Example:
        >>> with Timer("Processing", logger):
        ...     process_data()
        # Logs: "Processing completed in 1.23s"
    """

    def __init__(self, description: str, logger: logging.Logger | None = None) -> None:
        """Initialize timer.

        Args:
            description: Description of the operation.
            logger: Logger for output (uses print if None).
        """
        self.description = description
        self.logger = logger
        self.start_time: float = 0
        self.elapsed: float = 0

    def __enter__(self) -> "Timer":
        """Start timing."""
        import time

        self.start_time = time.perf_counter()
        return self

    def __exit__(self, *args: Any) -> None:
        """Stop timing and log result."""
        import time

        self.elapsed = time.perf_counter() - self.start_time

        message = f"{self.description} completed in {self.elapsed:.2f}s"
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
