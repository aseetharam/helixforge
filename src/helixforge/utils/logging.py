"""Logging setup (stdlib only)."""

from __future__ import annotations

import logging
import sys

_DEFAULT_NAME = "helixforge"


def setup_logging(
    level: str = "INFO",
    log_file: str | None = None,
    chunk_id: str | int | None = None,
) -> logging.Logger:
    """Configure the ``helixforge`` logger and return it.

    A ``chunk_id`` (e.g. a per-region worker id) is woven into the log format
    so parallel runs are distinguishable. Re-calling replaces existing handlers
    so logging stays idempotent across invocations.
    """
    logger = logging.getLogger(_DEFAULT_NAME)
    logger.setLevel(level)
    logger.handlers.clear()
    logger.propagate = False

    if chunk_id is not None:
        fmt = f"%(asctime)s [%(levelname)s] [{chunk_id}] %(name)s: %(message)s"
    else:
        fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    formatter = logging.Formatter(fmt)

    stream = logging.StreamHandler(sys.stderr)
    stream.setFormatter(formatter)
    logger.addHandler(stream)

    if log_file is not None:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_logger(name: str = _DEFAULT_NAME) -> logging.Logger:
    """Return a logger; non-root names become children of ``helixforge``."""
    if name == _DEFAULT_NAME or name.startswith(_DEFAULT_NAME + "."):
        return logging.getLogger(name)
    return logging.getLogger(f"{_DEFAULT_NAME}.{name}")
