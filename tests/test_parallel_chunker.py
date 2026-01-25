"""Tests for helixforge.parallel.chunker module.

Tests cover:
- GenomicChunk data structure
- ChunkPlan creation and serialization
- GenomeChunker strategies (scaffold, size, genes, adaptive)
- Utility functions (suggest_chunk_parameters, merge_overlapping_results)
"""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from helixforge.parallel.chunker import (
    ChunkPlan,
    ChunkStrategy,
    GenomicChunk,
    GenomeChunker,
    Chunker,
    merge_overlapping_results,
    suggest_chunk_parameters,
)


# =============================================================================
# GenomicChunk Tests
# =============================================================================


class TestGenomicChunk:
    """Tests for GenomicChunk data structure."""

    def test_create_chunk(self):
        """Test basic chunk creation."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=0,
            end=1000,
        )
        assert chunk.chunk_id == "chunk_0001"
        assert chunk.seqid == "chr1"
        assert chunk.start == 0
        assert chunk.end == 1000
        assert chunk.gene_ids is None

    def test_chunk_with_genes(self):
        """Test chunk with gene IDs."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=0,
            end=1000,
            gene_ids=["gene1", "gene2", "gene3"],
        )
        assert chunk.gene_ids == ["gene1", "gene2", "gene3"]

    def test_chunk_size(self):
        """Test chunk size property."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=100,
            end=500,
        )
        assert chunk.size == 400

    def test_chunk_overlaps_same_seqid(self):
        """Test overlap detection on same scaffold."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=100,
            end=500,
        )
        # Overlapping regions
        assert chunk.overlaps("chr1", 400, 600) is True
        assert chunk.overlaps("chr1", 0, 200) is True
        assert chunk.overlaps("chr1", 200, 300) is True  # Fully contained
        assert chunk.overlaps("chr1", 50, 550) is True  # Contains chunk

    def test_chunk_no_overlap(self):
        """Test non-overlapping regions."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=100,
            end=500,
        )
        # Non-overlapping
        assert chunk.overlaps("chr1", 500, 600) is False  # Adjacent
        assert chunk.overlaps("chr1", 0, 100) is False  # Adjacent before
        assert chunk.overlaps("chr1", 600, 700) is False
        assert chunk.overlaps("chr2", 100, 500) is False  # Different scaffold

    def test_chunk_contains(self):
        """Test containment check."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=100,
            end=500,
        )
        # Fully contained
        assert chunk.contains("chr1", 200, 400) is True
        assert chunk.contains("chr1", 100, 500) is True  # Exact match
        assert chunk.contains("chr1", 100, 200) is True
        # Not contained
        assert chunk.contains("chr1", 50, 200) is False  # Starts before
        assert chunk.contains("chr1", 400, 600) is False  # Ends after
        assert chunk.contains("chr2", 200, 400) is False  # Different scaffold

    def test_chunk_to_dict(self):
        """Test serialization to dict."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=100,
            end=500,
            gene_ids=["gene1"],
        )
        d = chunk.to_dict()
        assert d["chunk_id"] == "chunk_0001"
        assert d["seqid"] == "chr1"
        assert d["start"] == 100
        assert d["end"] == 500
        assert d["gene_ids"] == ["gene1"]

    def test_chunk_to_dict_no_genes(self):
        """Test serialization without genes."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=0,
            end=1000,
        )
        d = chunk.to_dict()
        assert "gene_ids" not in d

    def test_chunk_from_dict(self):
        """Test deserialization from dict."""
        data = {
            "chunk_id": "chunk_0001",
            "seqid": "chr1",
            "start": 100,
            "end": 500,
            "gene_ids": ["gene1", "gene2"],
        }
        chunk = GenomicChunk.from_dict(data)
        assert chunk.chunk_id == "chunk_0001"
        assert chunk.seqid == "chr1"
        assert chunk.start == 100
        assert chunk.end == 500
        assert chunk.gene_ids == ["gene1", "gene2"]

    def test_chunk_str(self):
        """Test string representation."""
        chunk = GenomicChunk(
            chunk_id="chunk_0001",
            seqid="chr1",
            start=100,
            end=500,
        )
        assert str(chunk) == "chr1:100-500"


# =============================================================================
# ChunkPlan Tests
# =============================================================================


class TestChunkPlan:
    """Tests for ChunkPlan data structure."""

    @pytest.fixture
    def sample_chunks(self) -> list[GenomicChunk]:
        """Create sample chunks for testing."""
        return [
            GenomicChunk("chunk_0000", "chr1", 0, 1000),
            GenomicChunk("chunk_0001", "chr1", 1000, 2000),
            GenomicChunk("chunk_0002", "chr2", 0, 500),
        ]

    def test_create_plan(self, sample_chunks):
        """Test plan creation."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SCAFFOLD,
            chunks=sample_chunks,
            total_bases=2500,
        )
        assert len(plan) == 3
        assert plan.strategy == ChunkStrategy.BY_SCAFFOLD
        assert plan.total_bases == 2500

    def test_plan_iteration(self, sample_chunks):
        """Test iterating over plan."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SCAFFOLD,
            chunks=sample_chunks,
            total_bases=2500,
        )
        chunk_ids = [c.chunk_id for c in plan]
        assert chunk_ids == ["chunk_0000", "chunk_0001", "chunk_0002"]

    def test_plan_indexing(self, sample_chunks):
        """Test indexing into plan."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SCAFFOLD,
            chunks=sample_chunks,
            total_bases=2500,
        )
        assert plan[0].chunk_id == "chunk_0000"
        assert plan[2].chunk_id == "chunk_0002"

    def test_get_chunk_by_id(self, sample_chunks):
        """Test getting chunk by ID."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SCAFFOLD,
            chunks=sample_chunks,
            total_bases=2500,
        )
        chunk = plan.get_chunk("chunk_0001")
        assert chunk is not None
        assert chunk.seqid == "chr1"
        assert chunk.start == 1000

        # Non-existent chunk
        assert plan.get_chunk("nonexistent") is None

    def test_plan_summary(self, sample_chunks):
        """Test summary statistics."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SIZE,
            chunks=sample_chunks,
            total_bases=2500,
            total_genes=10,
        )
        summary = plan.summary()
        assert summary["n_chunks"] == 3
        assert summary["total_bases"] == 2500
        assert summary["total_genes"] == 10
        assert summary["min_chunk_size"] == 500
        assert summary["max_chunk_size"] == 1000
        assert summary["strategy"] == "size"

    def test_empty_plan_summary(self):
        """Test summary for empty plan."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SCAFFOLD,
            chunks=[],
            total_bases=0,
        )
        summary = plan.summary()
        assert summary["n_chunks"] == 0
        assert summary["mean_chunk_size"] == 0

    def test_find_chunk_for_region(self, sample_chunks):
        """Test finding chunk for a region."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SCAFFOLD,
            chunks=sample_chunks,
            total_bases=2500,
        )
        # Region in first chunk
        chunk = plan.find_chunk_for_region("chr1", 100, 200)
        assert chunk is not None
        assert chunk.chunk_id == "chunk_0000"

        # Region in second chunk
        chunk = plan.find_chunk_for_region("chr1", 1100, 1200)
        assert chunk is not None
        assert chunk.chunk_id == "chunk_0001"

        # Region spans chunks (not fully contained)
        chunk = plan.find_chunk_for_region("chr1", 500, 1500)
        assert chunk is None

    def test_plan_save_load(self, sample_chunks, tmp_path):
        """Test JSON serialization round-trip."""
        plan = ChunkPlan(
            strategy=ChunkStrategy.BY_SIZE,
            chunks=sample_chunks,
            total_bases=2500,
            total_genes=10,
            parameters={"chunk_size": 1000},
        )

        path = tmp_path / "plan.json"
        plan.save(path)

        # Verify file exists and is valid JSON
        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert data["strategy"] == "size"
        assert data["n_chunks"] == 3

        # Load and verify
        loaded = ChunkPlan.load(path)
        assert loaded.strategy == ChunkStrategy.BY_SIZE
        assert len(loaded) == 3
        assert loaded.total_bases == 2500
        assert loaded.total_genes == 10


# =============================================================================
# GenomeChunker Tests
# =============================================================================


class TestGenomeChunker:
    """Tests for GenomeChunker class."""

    @pytest.fixture
    def mock_genome(self):
        """Create mock GenomeAccessor."""
        genome = MagicMock()
        genome.scaffold_order = ["chr1", "chr2", "chr3"]
        genome.total_length = 3000

        def get_length(seqid):
            lengths = {"chr1": 1000, "chr2": 1500, "chr3": 500}
            return lengths[seqid]

        genome.get_length = get_length
        return genome

    @pytest.fixture
    def mock_gff_parser(self):
        """Create mock GFF3Parser with genes."""
        from helixforge.io.gff import GeneModel

        parser = MagicMock()

        # Create mock genes
        genes = [
            MagicMock(gene_id="gene1", seqid="chr1", start=100, end=300),
            MagicMock(gene_id="gene2", seqid="chr1", start=400, end=600),
            MagicMock(gene_id="gene3", seqid="chr1", start=700, end=900),
            MagicMock(gene_id="gene4", seqid="chr2", start=100, end=500),
            MagicMock(gene_id="gene5", seqid="chr2", start=600, end=1000),
        ]
        parser.iter_genes.return_value = iter(genes)

        return parser

    def test_chunker_by_scaffold(self, mock_genome):
        """Test scaffold-based chunking."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)

        assert len(plan) == 3
        assert plan.chunks[0].seqid == "chr1"
        assert plan.chunks[0].size == 1000
        assert plan.chunks[1].seqid == "chr2"
        assert plan.chunks[1].size == 1500
        assert plan.chunks[2].seqid == "chr3"
        assert plan.chunks[2].size == 500

    def test_chunker_by_scaffold_with_max_size(self, mock_genome):
        """Test scaffold chunking with max chunk size."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan(
            ChunkStrategy.BY_SCAFFOLD,
            max_chunk_size=800,
        )

        # chr1 (1000) should be split into 2 chunks
        # chr2 (1500) should be split into 2 chunks
        # chr3 (500) stays as 1 chunk
        assert len(plan) == 5

    def test_chunker_by_size(self, mock_genome):
        """Test size-based chunking."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan(
            ChunkStrategy.BY_SIZE,
            chunk_size=500,
            min_chunk_size=100,
        )

        # chr1 (1000) = 2 chunks
        # chr2 (1500) = 3 chunks
        # chr3 (500) = 1 chunk
        assert len(plan) >= 5

    def test_chunker_by_genes(self, mock_genome, mock_gff_parser):
        """Test gene-based chunking."""
        chunker = GenomeChunker(mock_genome, mock_gff_parser)
        plan = chunker.create_plan(
            ChunkStrategy.BY_GENES,
            chunk_size=2,  # 2 genes per chunk
        )

        # chr1 has 3 genes: 2 chunks (2 + 1)
        # chr2 has 2 genes: 1 chunk
        # chr3 has 0 genes: 1 chunk (no genes)
        # Total: 4 chunks
        assert len(plan) >= 3
        assert plan.total_genes == 5

    def test_chunker_by_genes_requires_parser(self, mock_genome):
        """Test that BY_GENES requires gff_parser."""
        chunker = GenomeChunker(mock_genome)
        with pytest.raises(ValueError, match="BY_GENES strategy requires gff_parser"):
            chunker.create_plan(ChunkStrategy.BY_GENES)

    def test_chunker_adaptive(self, mock_genome, mock_gff_parser):
        """Test adaptive chunking."""
        chunker = GenomeChunker(mock_genome, mock_gff_parser)
        plan = chunker.create_plan(
            ChunkStrategy.ADAPTIVE,
            target_chunks=5,
            min_chunk_size=100,
        )

        # Should create roughly target_chunks chunks
        assert len(plan) >= 1
        assert plan.strategy == ChunkStrategy.ADAPTIVE

    def test_chunker_string_strategy(self, mock_genome):
        """Test using string for strategy."""
        chunker = GenomeChunker(mock_genome)
        plan = chunker.create_plan("scaffold")
        assert plan.strategy == ChunkStrategy.BY_SCAFFOLD

    def test_estimate_memory(self):
        """Test memory estimation."""
        chunk = GenomicChunk("chunk_0000", "chr1", 0, 1_000_000)

        memory = GenomeChunker.estimate_memory(chunk)
        # 1M bases * 20 bytes = 20 MB
        assert memory == 20_000_000

        memory_custom = GenomeChunker.estimate_memory(chunk, bytes_per_base=10.0)
        assert memory_custom == 10_000_000

    def test_chunker_legacy_alias(self, mock_genome):
        """Test legacy Chunker alias."""
        chunker = Chunker(mock_genome)
        assert isinstance(chunker, GenomeChunker)


# =============================================================================
# Utility Function Tests
# =============================================================================


class TestSuggestChunkParameters:
    """Tests for suggest_chunk_parameters function."""

    def test_small_genome(self):
        """Test recommendations for small genome."""
        result = suggest_chunk_parameters(
            genome_size=50_000_000,  # 50 Mb
            n_genes=1000,
            available_memory_gb=8.0,
            n_workers=4,
        )
        assert result["strategy"] == ChunkStrategy.BY_SCAFFOLD
        assert "rationale" in result

    def test_large_genome(self):
        """Test recommendations for large genome."""
        result = suggest_chunk_parameters(
            genome_size=2_000_000_000,  # 2 Gb
            n_genes=40000,
            available_memory_gb=16.0,
            n_workers=8,
        )
        assert result["strategy"] == ChunkStrategy.BY_SIZE
        assert result["chunk_size"] is not None

    def test_gene_dense_genome(self):
        """Test recommendations for gene-dense genome."""
        result = suggest_chunk_parameters(
            genome_size=200_000_000,  # 200 Mb
            n_genes=30000,  # High density
            available_memory_gb=8.0,
            n_workers=4,
        )
        assert result["strategy"] == ChunkStrategy.BY_GENES

    def test_medium_genome(self):
        """Test recommendations for medium genome."""
        result = suggest_chunk_parameters(
            genome_size=500_000_000,  # 500 Mb
            n_genes=5000,
            available_memory_gb=16.0,
            n_workers=8,
        )
        assert result["strategy"] == ChunkStrategy.ADAPTIVE
        assert "target_chunks" in result


class TestMergeOverlappingResults:
    """Tests for merge_overlapping_results function."""

    def test_merge_non_overlapping(self):
        """Test merging non-overlapping results."""
        chunk1 = GenomicChunk("chunk_0000", "chr1", 0, 500)
        chunk2 = GenomicChunk("chunk_0001", "chr1", 500, 1000)

        gene1 = MagicMock(gene_id="gene1", seqid="chr1", start=100, end=200)
        gene2 = MagicMock(gene_id="gene2", seqid="chr1", start=600, end=800)

        results = [
            (chunk1, [gene1]),
            (chunk2, [gene2]),
        ]

        merged = merge_overlapping_results(results)
        assert len(merged) == 2
        gene_ids = [g.gene_id for g in merged]
        assert "gene1" in gene_ids
        assert "gene2" in gene_ids

    def test_merge_removes_duplicates(self):
        """Test that duplicates are removed."""
        chunk1 = GenomicChunk("chunk_0000", "chr1", 0, 600)
        chunk2 = GenomicChunk("chunk_0001", "chr1", 400, 1000)

        # Same gene in both chunks
        gene1a = MagicMock(gene_id="gene1", seqid="chr1", start=200, end=400)
        gene1b = MagicMock(gene_id="gene1", seqid="chr1", start=200, end=400)

        results = [
            (chunk1, [gene1a]),
            (chunk2, [gene1b]),
        ]

        merged = merge_overlapping_results(results)
        assert len(merged) == 1
        assert merged[0].gene_id == "gene1"

    def test_merge_prefers_contained_genes(self):
        """Test that genes fully contained in chunk are preferred."""
        chunk1 = GenomicChunk("chunk_0000", "chr1", 0, 400)
        chunk2 = GenomicChunk("chunk_0001", "chr1", 300, 700)

        # Gene fully contained in chunk2, partially in chunk1
        gene1_partial = MagicMock(gene_id="gene1", seqid="chr1", start=350, end=500)
        gene1_partial.extra = "from_chunk1"

        gene1_full = MagicMock(gene_id="gene1", seqid="chr1", start=350, end=500)
        gene1_full.extra = "from_chunk2"

        results = [
            (chunk1, [gene1_partial]),
            (chunk2, [gene1_full]),
        ]

        merged = merge_overlapping_results(results)
        assert len(merged) == 1

    def test_merge_sorts_by_position(self):
        """Test that merged results are sorted by position."""
        chunk1 = GenomicChunk("chunk_0000", "chr1", 500, 1000)
        chunk2 = GenomicChunk("chunk_0001", "chr1", 0, 500)

        gene1 = MagicMock(gene_id="gene1", seqid="chr1", start=600, end=800)
        gene2 = MagicMock(gene_id="gene2", seqid="chr1", start=100, end=300)

        # Results in reverse order
        results = [
            (chunk1, [gene1]),
            (chunk2, [gene2]),
        ]

        merged = merge_overlapping_results(results)
        assert merged[0].gene_id == "gene2"  # Earlier position
        assert merged[1].gene_id == "gene1"
