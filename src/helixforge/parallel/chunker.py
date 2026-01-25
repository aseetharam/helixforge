"""Genome chunking for parallel processing.

This module provides tools for dividing genomes into chunks suitable
for parallel processing on local systems or HPC clusters.

Features:
    - Multiple chunking strategies (scaffold, size, genes, adaptive)
    - Memory estimation for chunk processing
    - JSON serialization for SLURM array jobs
    - Gene-aware chunking to avoid splitting genes

Example:
    >>> from helixforge.parallel.chunker import GenomeChunker, ChunkStrategy
    >>> from helixforge.io.fasta import GenomeAccessor
    >>> genome = GenomeAccessor("genome.fa")
    >>> chunker = GenomeChunker(genome)
    >>> plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)
    >>> for chunk in plan:
    ...     print(f"{chunk.seqid}:{chunk.start}-{chunk.end}")
"""

from __future__ import annotations

import json
import logging
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Iterator

import attrs

if TYPE_CHECKING:
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser, GeneModel

logger = logging.getLogger(__name__)


# =============================================================================
# Enums
# =============================================================================


class ChunkStrategy(Enum):
    """Available chunking strategies."""

    BY_SCAFFOLD = "scaffold"  # One chunk per scaffold
    BY_SIZE = "size"  # Fixed base-pair windows
    BY_GENES = "genes"  # Fixed number of genes per chunk
    ADAPTIVE = "adaptive"  # Balance chunk sizes automatically


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class GenomicChunk:
    """Represents a processing unit.

    Attributes:
        chunk_id: Unique identifier for this chunk.
        seqid: Scaffold/chromosome name.
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
        gene_ids: Optional list of gene IDs in this chunk.
    """

    chunk_id: str
    seqid: str
    start: int  # 0-based
    end: int  # 0-based, exclusive
    gene_ids: list[str] | None = None

    @property
    def size(self) -> int:
        """Get chunk size in base pairs."""
        return self.end - self.start

    def overlaps(self, seqid: str, start: int, end: int) -> bool:
        """Check if region overlaps this chunk.

        Args:
            seqid: Scaffold name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            True if regions overlap.
        """
        if self.seqid != seqid:
            return False
        return self.start < end and start < self.end

    def contains(self, seqid: str, start: int, end: int) -> bool:
        """Check if this chunk fully contains a region.

        Args:
            seqid: Scaffold name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            True if this chunk contains the region.
        """
        if self.seqid != seqid:
            return False
        return self.start <= start and end <= self.end

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        result = {
            "chunk_id": self.chunk_id,
            "seqid": self.seqid,
            "start": self.start,
            "end": self.end,
        }
        if self.gene_ids is not None:
            result["gene_ids"] = self.gene_ids
        return result

    @classmethod
    def from_dict(cls, data: dict) -> GenomicChunk:
        """Create from dictionary."""
        return cls(
            chunk_id=data["chunk_id"],
            seqid=data["seqid"],
            start=data["start"],
            end=data["end"],
            gene_ids=data.get("gene_ids"),
        )

    def __str__(self) -> str:
        """String representation."""
        return f"{self.seqid}:{self.start}-{self.end}"


@attrs.define(slots=True)
class ChunkPlan:
    """Complete chunking plan for a genome.

    Attributes:
        strategy: Chunking strategy used.
        chunks: List of genomic chunks.
        total_bases: Total genome size in bases.
        total_genes: Total number of genes (if available).
        parameters: Strategy-specific parameters.
    """

    strategy: ChunkStrategy
    chunks: list[GenomicChunk]
    total_bases: int
    total_genes: int | None = None
    parameters: dict = attrs.Factory(dict)

    def __len__(self) -> int:
        """Return number of chunks."""
        return len(self.chunks)

    def __iter__(self) -> Iterator[GenomicChunk]:
        """Iterate over chunks."""
        return iter(self.chunks)

    def __getitem__(self, index: int) -> GenomicChunk:
        """Get chunk by index."""
        return self.chunks[index]

    def get_chunk(self, chunk_id: str) -> GenomicChunk | None:
        """Get chunk by ID.

        Args:
            chunk_id: Chunk identifier.

        Returns:
            GenomicChunk if found, None otherwise.
        """
        for chunk in self.chunks:
            if chunk.chunk_id == chunk_id:
                return chunk
        return None

    def get_chunk_by_index(self, index: int) -> GenomicChunk:
        """Get chunk by 0-based index.

        Args:
            index: Chunk index.

        Returns:
            GenomicChunk at that index.

        Raises:
            IndexError: If index out of range.
        """
        return self.chunks[index]

    def save(self, path: Path | str) -> None:
        """Save plan to JSON for SLURM array jobs.

        Args:
            path: Output file path.
        """
        path = Path(path)
        data = {
            "strategy": self.strategy.value,
            "total_bases": self.total_bases,
            "total_genes": self.total_genes,
            "parameters": self.parameters,
            "n_chunks": len(self.chunks),
            "chunks": [chunk.to_dict() for chunk in self.chunks],
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        logger.info(f"Saved chunk plan to {path}")

    @classmethod
    def load(cls, path: Path | str) -> ChunkPlan:
        """Load plan from JSON.

        Args:
            path: Input file path.

        Returns:
            ChunkPlan instance.
        """
        path = Path(path)
        with open(path) as f:
            data = json.load(f)

        chunks = [GenomicChunk.from_dict(c) for c in data["chunks"]]

        return cls(
            strategy=ChunkStrategy(data["strategy"]),
            chunks=chunks,
            total_bases=data["total_bases"],
            total_genes=data.get("total_genes"),
            parameters=data.get("parameters", {}),
        )

    def summary(self) -> dict:
        """Return summary statistics.

        Returns:
            Dictionary with chunk statistics.
        """
        if not self.chunks:
            return {
                "n_chunks": 0,
                "total_bases": self.total_bases,
                "mean_chunk_size": 0,
                "min_chunk_size": 0,
                "max_chunk_size": 0,
            }

        sizes = [c.size for c in self.chunks]
        return {
            "n_chunks": len(self.chunks),
            "total_bases": self.total_bases,
            "total_genes": self.total_genes,
            "mean_chunk_size": sum(sizes) / len(sizes),
            "median_chunk_size": sorted(sizes)[len(sizes) // 2],
            "min_chunk_size": min(sizes),
            "max_chunk_size": max(sizes),
            "strategy": self.strategy.value,
        }

    def find_chunk_for_region(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> GenomicChunk | None:
        """Find chunk that contains a region.

        Args:
            seqid: Scaffold name.
            start: Start position.
            end: End position.

        Returns:
            Chunk containing the region, or None.
        """
        for chunk in self.chunks:
            if chunk.contains(seqid, start, end):
                return chunk
        return None


# =============================================================================
# Main Chunker Class
# =============================================================================


class GenomeChunker:
    """Partition genome into processing chunks.

    Supports multiple strategies optimized for different scenarios:
    - BY_SCAFFOLD: Simple, good for genomes with many scaffolds
    - BY_SIZE: Fixed windows, good for few large chromosomes
    - BY_GENES: Balance gene count, good for annotation tasks
    - ADAPTIVE: Automatically balance based on content

    Attributes:
        genome: GenomeAccessor for the target genome.
        gff_parser: Optional GFF3Parser for gene-aware chunking.

    Example:
        >>> chunker = GenomeChunker(genome)
        >>> plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)
        >>> print(plan.summary())
    """

    def __init__(
        self,
        genome: GenomeAccessor,
        gff_parser: GFF3Parser | None = None,
    ) -> None:
        """Initialize the chunker.

        Args:
            genome: GenomeAccessor for the genome.
            gff_parser: Optional GFF3Parser for gene-aware chunking.
        """
        self.genome = genome
        self.gff_parser = gff_parser
        self._genes_by_scaffold: dict[str, list[GeneModel]] | None = None

    def _load_genes(self) -> dict[str, list[GeneModel]]:
        """Load and index genes by scaffold.

        Returns:
            Dictionary mapping seqid to list of genes.
        """
        if self._genes_by_scaffold is not None:
            return self._genes_by_scaffold

        if self.gff_parser is None:
            self._genes_by_scaffold = {}
            return self._genes_by_scaffold

        self._genes_by_scaffold = {}
        for gene in self.gff_parser.iter_genes():
            if gene.seqid not in self._genes_by_scaffold:
                self._genes_by_scaffold[gene.seqid] = []
            self._genes_by_scaffold[gene.seqid].append(gene)

        # Sort genes by position on each scaffold
        for seqid in self._genes_by_scaffold:
            self._genes_by_scaffold[seqid].sort(key=lambda g: (g.start, g.end))

        return self._genes_by_scaffold

    def create_plan(
        self,
        strategy: ChunkStrategy | str,
        chunk_size: int | None = None,
        min_chunk_size: int = 100_000,
        max_chunk_size: int | None = None,
        target_chunks: int | None = None,
    ) -> ChunkPlan:
        """Create chunking plan.

        Args:
            strategy: Chunking strategy to use.
            chunk_size: Size parameter (bases for BY_SIZE, genes for BY_GENES).
            min_chunk_size: Don't create chunks smaller than this.
            max_chunk_size: Split chunks larger than this.
            target_chunks: Target number of chunks for ADAPTIVE.

        Returns:
            ChunkPlan ready for parallel execution.

        Raises:
            ValueError: If strategy requires gff_parser but none provided.
        """
        if isinstance(strategy, str):
            strategy = ChunkStrategy(strategy)

        parameters = {
            "chunk_size": chunk_size,
            "min_chunk_size": min_chunk_size,
            "max_chunk_size": max_chunk_size,
            "target_chunks": target_chunks,
        }

        total_genes = None
        if strategy == ChunkStrategy.BY_SCAFFOLD:
            chunks = self._chunk_by_scaffold(max_chunk_size)
        elif strategy == ChunkStrategy.BY_SIZE:
            if chunk_size is None:
                chunk_size = 10_000_000  # 10 Mb default
            chunks = self._chunk_by_size(chunk_size, min_chunk_size)
        elif strategy == ChunkStrategy.BY_GENES:
            if self.gff_parser is None:
                raise ValueError("BY_GENES strategy requires gff_parser")
            if chunk_size is None:
                chunk_size = 100  # 100 genes default
            chunks = self._chunk_by_genes(chunk_size)
            genes_by_scaffold = self._load_genes()
            total_genes = sum(len(genes) for genes in genes_by_scaffold.values())
        elif strategy == ChunkStrategy.ADAPTIVE:
            if target_chunks is None:
                # Default: aim for chunks that fit in reasonable memory
                target_chunks = max(1, self.genome.total_length // 50_000_000)
            chunks = self._chunk_adaptive(target_chunks, min_chunk_size)
            if self.gff_parser is not None:
                genes_by_scaffold = self._load_genes()
                total_genes = sum(len(genes) for genes in genes_by_scaffold.values())
        else:
            raise ValueError(f"Unknown strategy: {strategy}")

        logger.info(
            f"Created chunk plan: {len(chunks)} chunks, "
            f"strategy={strategy.value}, total_bases={self.genome.total_length:,}"
        )

        return ChunkPlan(
            strategy=strategy,
            chunks=chunks,
            total_bases=self.genome.total_length,
            total_genes=total_genes,
            parameters=parameters,
        )

    def _chunk_by_scaffold(
        self,
        max_chunk_size: int | None,
    ) -> list[GenomicChunk]:
        """One chunk per scaffold, optionally split large scaffolds.

        Args:
            max_chunk_size: Maximum chunk size; split larger scaffolds.

        Returns:
            List of GenomicChunk objects.
        """
        chunks = []
        chunk_idx = 0

        for seqid in self.genome.scaffold_order:
            length = self.genome.get_length(seqid)

            if max_chunk_size is None or length <= max_chunk_size:
                # Single chunk for this scaffold
                chunks.append(
                    GenomicChunk(
                        chunk_id=f"chunk_{chunk_idx:04d}",
                        seqid=seqid,
                        start=0,
                        end=length,
                    )
                )
                chunk_idx += 1
            else:
                # Split scaffold into multiple chunks
                pos = 0
                while pos < length:
                    end = min(pos + max_chunk_size, length)
                    chunks.append(
                        GenomicChunk(
                            chunk_id=f"chunk_{chunk_idx:04d}",
                            seqid=seqid,
                            start=pos,
                            end=end,
                        )
                    )
                    chunk_idx += 1
                    pos = end

        return chunks

    def _chunk_by_size(
        self,
        chunk_size: int,
        min_chunk_size: int,
    ) -> list[GenomicChunk]:
        """Fixed-size windows across genome.

        Handles scaffold boundaries (don't span scaffolds).
        Last chunk per scaffold may be smaller.

        Args:
            chunk_size: Target chunk size in bases.
            min_chunk_size: Minimum chunk size.

        Returns:
            List of GenomicChunk objects.
        """
        chunks = []
        chunk_idx = 0

        for seqid in self.genome.scaffold_order:
            length = self.genome.get_length(seqid)

            if length < min_chunk_size:
                # Small scaffold becomes single chunk
                chunks.append(
                    GenomicChunk(
                        chunk_id=f"chunk_{chunk_idx:04d}",
                        seqid=seqid,
                        start=0,
                        end=length,
                    )
                )
                chunk_idx += 1
                continue

            # Split into chunks of chunk_size
            pos = 0
            while pos < length:
                end = min(pos + chunk_size, length)

                # Check if remaining is too small
                if length - end < min_chunk_size and length - end > 0:
                    # Extend this chunk to the end
                    end = length

                chunks.append(
                    GenomicChunk(
                        chunk_id=f"chunk_{chunk_idx:04d}",
                        seqid=seqid,
                        start=pos,
                        end=end,
                    )
                )
                chunk_idx += 1
                pos = end

        return chunks

    def _chunk_by_genes(
        self,
        genes_per_chunk: int,
    ) -> list[GenomicChunk]:
        """Group genes into chunks.

        Respects scaffold boundaries.
        Chunk coordinates span from first to last gene in group.

        Args:
            genes_per_chunk: Target number of genes per chunk.

        Returns:
            List of GenomicChunk objects.
        """
        genes_by_scaffold = self._load_genes()
        chunks = []
        chunk_idx = 0

        for seqid in self.genome.scaffold_order:
            genes = genes_by_scaffold.get(seqid, [])
            length = self.genome.get_length(seqid)

            if not genes:
                # Scaffold has no genes, create a single chunk
                chunks.append(
                    GenomicChunk(
                        chunk_id=f"chunk_{chunk_idx:04d}",
                        seqid=seqid,
                        start=0,
                        end=length,
                        gene_ids=[],
                    )
                )
                chunk_idx += 1
                continue

            # Group genes
            for i in range(0, len(genes), genes_per_chunk):
                batch = genes[i : i + genes_per_chunk]
                gene_ids = [g.gene_id for g in batch]

                # Chunk spans from first gene start to last gene end
                chunk_start = min(g.start for g in batch)
                chunk_end = max(g.end for g in batch)

                chunks.append(
                    GenomicChunk(
                        chunk_id=f"chunk_{chunk_idx:04d}",
                        seqid=seqid,
                        start=chunk_start,
                        end=chunk_end,
                        gene_ids=gene_ids,
                    )
                )
                chunk_idx += 1

        return chunks

    def _chunk_adaptive(
        self,
        target_chunks: int,
        min_chunk_size: int,
    ) -> list[GenomicChunk]:
        """Automatically balance chunks.

        Considers both size and gene density.
        Aims for roughly equal processing time per chunk.

        Args:
            target_chunks: Target number of chunks.
            min_chunk_size: Minimum chunk size in bases.

        Returns:
            List of GenomicChunk objects.
        """
        total_size = self.genome.total_length
        target_chunk_size = max(min_chunk_size, total_size // target_chunks)

        # Get gene counts per scaffold if available
        genes_by_scaffold = {}
        if self.gff_parser is not None:
            genes_by_scaffold = self._load_genes()

        chunks = []
        chunk_idx = 0

        for seqid in self.genome.scaffold_order:
            length = self.genome.get_length(seqid)
            n_genes = len(genes_by_scaffold.get(seqid, []))

            # Adjust chunk size based on gene density
            if n_genes > 0:
                # Gene-dense regions get smaller chunks
                gene_density = n_genes / length
                avg_gene_density = (
                    sum(len(g) for g in genes_by_scaffold.values()) / total_size
                    if genes_by_scaffold
                    else 0
                )

                if avg_gene_density > 0:
                    density_factor = gene_density / avg_gene_density
                    adjusted_size = int(target_chunk_size / max(0.5, density_factor))
                    adjusted_size = max(min_chunk_size, adjusted_size)
                else:
                    adjusted_size = target_chunk_size
            else:
                adjusted_size = target_chunk_size

            if length <= adjusted_size:
                # Single chunk for scaffold
                gene_ids = (
                    [g.gene_id for g in genes_by_scaffold.get(seqid, [])]
                    if genes_by_scaffold
                    else None
                )
                chunks.append(
                    GenomicChunk(
                        chunk_id=f"chunk_{chunk_idx:04d}",
                        seqid=seqid,
                        start=0,
                        end=length,
                        gene_ids=gene_ids,
                    )
                )
                chunk_idx += 1
            else:
                # Split scaffold
                genes = genes_by_scaffold.get(seqid, [])
                pos = 0
                gene_idx = 0

                while pos < length:
                    end = min(pos + adjusted_size, length)

                    # Try to avoid splitting genes
                    if genes and gene_idx < len(genes):
                        # Find genes in this range
                        chunk_gene_ids = []
                        while gene_idx < len(genes) and genes[gene_idx].start < end:
                            gene = genes[gene_idx]
                            if gene.end <= end or gene.start >= pos:
                                chunk_gene_ids.append(gene.gene_id)
                            # If gene spans boundary, extend chunk
                            if gene.start < end < gene.end:
                                end = gene.end
                            gene_idx += 1
                    else:
                        chunk_gene_ids = None

                    chunks.append(
                        GenomicChunk(
                            chunk_id=f"chunk_{chunk_idx:04d}",
                            seqid=seqid,
                            start=pos,
                            end=end,
                            gene_ids=chunk_gene_ids,
                        )
                    )
                    chunk_idx += 1
                    pos = end

        return chunks

    @staticmethod
    def estimate_memory(
        chunk: GenomicChunk,
        bytes_per_base: float = 20.0,
    ) -> int:
        """Estimate memory requirement for processing a chunk.

        Args:
            chunk: Chunk to estimate.
            bytes_per_base: Memory multiplier (HDF5 + working data).

        Returns:
            Estimated bytes required.
        """
        return int(chunk.size * bytes_per_base)


# =============================================================================
# Utility Functions
# =============================================================================


def suggest_chunk_parameters(
    genome_size: int,
    n_genes: int,
    available_memory_gb: float,
    n_workers: int,
) -> dict:
    """Suggest chunking parameters based on resources.

    Args:
        genome_size: Total genome size in bases.
        n_genes: Total number of genes.
        available_memory_gb: Available memory in GB.
        n_workers: Number of parallel workers.

    Returns:
        Dict with recommended strategy, chunk_size, and rationale.
    """
    # Memory per worker
    memory_per_worker_mb = (available_memory_gb * 1024) / n_workers

    # Conservative estimate: 20 bytes per base for processing
    bytes_per_base = 20.0
    max_chunk_bases = int((memory_per_worker_mb * 1024 * 1024) / bytes_per_base)

    # Aim for at least 2x workers worth of chunks for load balancing
    min_chunks = n_workers * 2

    result = {
        "memory_per_worker_mb": memory_per_worker_mb,
        "max_chunk_bases": max_chunk_bases,
    }

    if genome_size < 100_000_000:  # < 100 Mb
        # Small genome
        result["strategy"] = ChunkStrategy.BY_SCAFFOLD
        result["chunk_size"] = None
        result["rationale"] = "Small genome; one chunk per scaffold is efficient"
    elif n_genes > 0 and genome_size / n_genes < 50_000:
        # Gene-dense genome
        genes_per_chunk = max(50, n_genes // (min_chunks * 2))
        result["strategy"] = ChunkStrategy.BY_GENES
        result["chunk_size"] = genes_per_chunk
        result["rationale"] = "Gene-dense genome; balance by gene count"
    elif genome_size > 1_000_000_000:  # > 1 Gb
        # Large genome
        chunk_size = min(max_chunk_bases, 50_000_000)  # Max 50 Mb chunks
        result["strategy"] = ChunkStrategy.BY_SIZE
        result["chunk_size"] = chunk_size
        result["rationale"] = "Large genome; fixed-size chunks for memory control"
    else:
        # Medium genome, use adaptive
        target_chunks = max(min_chunks, genome_size // 10_000_000)
        result["strategy"] = ChunkStrategy.ADAPTIVE
        result["target_chunks"] = target_chunks
        result["chunk_size"] = None
        result["rationale"] = "Medium genome; adaptive chunking for balance"

    return result


def merge_overlapping_results(
    results: list[tuple[GenomicChunk, list[GeneModel]]],
) -> list[GeneModel]:
    """Merge results from overlapping chunks.

    Handles genes that span chunk boundaries by keeping the
    version from the chunk where the gene is fully contained.

    Args:
        results: List of (chunk, genes) tuples.

    Returns:
        Merged list of genes without duplicates.
    """
    seen_gene_ids: set[str] = set()
    merged_genes: list[GeneModel] = []

    # Sort by chunk position for consistent ordering
    sorted_results = sorted(results, key=lambda x: (x[0].seqid, x[0].start))

    for chunk, genes in sorted_results:
        for gene in genes:
            if gene.gene_id in seen_gene_ids:
                continue

            # Prefer genes fully contained in chunk
            if chunk.contains(gene.seqid, gene.start, gene.end):
                merged_genes.append(gene)
                seen_gene_ids.add(gene.gene_id)

    # Second pass: add genes that weren't fully contained anywhere
    for chunk, genes in sorted_results:
        for gene in genes:
            if gene.gene_id not in seen_gene_ids:
                merged_genes.append(gene)
                seen_gene_ids.add(gene.gene_id)

    # Sort by position
    merged_genes.sort(key=lambda g: (g.seqid, g.start, g.end))

    return merged_genes


# =============================================================================
# Legacy Compatibility
# =============================================================================

# Alias for backward compatibility
Chunker = GenomeChunker
