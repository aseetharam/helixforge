"""Parallelization utilities for HelixForge.

This module provides tools for parallel execution of HelixForge
operations:

- Genome chunking for parallel processing
- Local multiprocessing execution
- SLURM cluster submission

Example:
    >>> from helixforge.parallel import Chunker, Executor
    >>> chunker = Chunker(chunk_size=1_000_000)
    >>> chunks = chunker.chunk_genome(genome)
    >>> executor = Executor(workers=8)
    >>> results = executor.map(process_chunk, chunks)
"""

# TODO: Import and expose main classes once implemented
# from helixforge.parallel.chunker import Chunker
# from helixforge.parallel.executor import Executor
# from helixforge.parallel.slurm import SlurmSubmitter

__all__: list[str] = []
