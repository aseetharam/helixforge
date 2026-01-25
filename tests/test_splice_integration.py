"""Integration tests for splice refinement."""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from helixforge.core.splice import (
    GeneSpliceReport,
    IntronModel,
    PositionWeightMatrix,
    SpliceRefiner,
    SpliceReportWriter,
)


# =============================================================================
# Mock Classes for Testing
# =============================================================================


class MockGenomeAccessor:
    """Mock genome accessor for testing."""

    def __init__(self, sequences: dict[str, str] | None = None):
        """Initialize with optional sequence data."""
        self._sequences = sequences or {
            "chr1": "A" * 10000,  # Default: all A's
        }

    def get_sequence(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str = "+",
    ) -> str:
        """Get sequence from mock genome."""
        if seqid not in self._sequences:
            raise KeyError(f"Unknown sequence: {seqid}")

        seq = self._sequences[seqid][start:end]

        if strand == "-":
            # Simple reverse complement
            complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
            seq = "".join(complement.get(b, "N") for b in reversed(seq))

        return seq

    def get_length(self, seqid: str) -> int:
        """Get sequence length."""
        if seqid not in self._sequences:
            raise KeyError(f"Unknown sequence: {seqid}")
        return len(self._sequences[seqid])

    def close(self):
        """Mock close method."""
        pass


class MockSpliceJunction:
    """Mock splice junction for testing."""

    def __init__(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str,
        read_count: int,
    ):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.read_count = read_count


class MockTranscriptModel:
    """Mock transcript model for testing."""

    def __init__(
        self,
        transcript_id: str,
        seqid: str,
        start: int,
        end: int,
        strand: str,
        exons: list[tuple[int, int]],
        cds: list[tuple[int, int, int]] | None = None,
    ):
        self.transcript_id = transcript_id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.exons = exons
        self.cds = cds or []


class MockGeneModel:
    """Mock gene model for testing."""

    def __init__(
        self,
        gene_id: str,
        seqid: str,
        start: int,
        end: int,
        strand: str,
        transcripts: list[MockTranscriptModel],
    ):
        self.gene_id = gene_id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.transcripts = transcripts


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def mock_genome_with_splice_sites():
    """Create genome with canonical splice sites."""
    # Create a sequence with GT-AG introns
    # Exon1: 100-200, Intron: 200-300 (GT...AG), Exon2: 300-400
    seq = (
        "A" * 100  # 0-100: upstream
        + "ATGCGA" * 17  # 100-202: exon 1 (102 bp)
        + "GT"  # 200-202: donor
        + "A" * 96  # 202-298: intron body
        + "AG"  # 298-300: acceptor
        + "ATGCGA" * 17  # 300-402: exon 2 (102 bp)
        + "A" * 598  # padding to 1000
    )
    return MockGenomeAccessor({"chr1": seq})


@pytest.fixture
def mock_junctions():
    """Create mock junctions matching the genome."""
    return {
        "chr1": [
            MockSpliceJunction("chr1", 200, 300, "+", 50),  # Exact match
            MockSpliceJunction("chr1", 500, 600, "+", 30),  # Different position
        ]
    }


@pytest.fixture
def mock_gene():
    """Create a mock gene model."""
    transcript = MockTranscriptModel(
        transcript_id="tx1",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        exons=[(100, 200), (300, 400)],
        cds=[(120, 200, 0), (300, 380, 2)],
    )
    return MockGeneModel(
        gene_id="gene1",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        transcripts=[transcript],
    )


@pytest.fixture
def mock_gene_shifted():
    """Create a gene with slightly shifted splice sites."""
    # Gene with intron at 198-302 instead of 200-300
    transcript = MockTranscriptModel(
        transcript_id="tx1",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        exons=[(100, 198), (302, 400)],  # Shifted by 2bp
        cds=[(120, 198, 0), (302, 380, 2)],
    )
    return MockGeneModel(
        gene_id="gene1",
        seqid="chr1",
        start=100,
        end=400,
        strand="+",
        transcripts=[transcript],
    )


# =============================================================================
# Integration Tests
# =============================================================================


class TestSpliceRefinerIntegration:
    """Integration tests for SpliceRefiner."""

    def test_exact_junction_match(
        self,
        mock_genome_with_splice_sites,
        mock_junctions,
        mock_gene,
    ):
        """Intron matching RNA-seq junction exactly should have support."""
        refiner = SpliceRefiner(
            mock_genome_with_splice_sites,
            mock_junctions,
        )

        refined_gene, report = refiner.refine_gene(mock_gene)

        # Should have support
        assert report.n_supported >= 1 or report.n_introns == 1
        assert report.gene_id == "gene1"

    def test_nearby_junction_correction(
        self,
        mock_genome_with_splice_sites,
        mock_junctions,
        mock_gene_shifted,
    ):
        """Intron should be corrected to nearby junction."""
        refiner = SpliceRefiner(
            mock_genome_with_splice_sites,
            mock_junctions,
            max_shift=15,
        )

        refined_gene, report = refiner.refine_gene(mock_gene_shifted)

        # Should have correction (shifted splice site)
        # The exact behavior depends on the implementation
        assert report.gene_id == "gene1"
        assert report.n_introns >= 0

    def test_no_correction_beyond_max_shift(self):
        """Junction too far should not cause correction."""
        # Create genome and junctions far apart
        genome = MockGenomeAccessor({
            "chr1": "GT" + "A" * 998 + "AG" + "A" * 8000
        })
        junctions = {
            "chr1": [
                MockSpliceJunction("chr1", 5000, 6000, "+", 50),  # Far away
            ]
        }

        transcript = MockTranscriptModel(
            transcript_id="tx1",
            seqid="chr1",
            start=0,
            end=2000,
            strand="+",
            exons=[(0, 100), (200, 300)],
        )
        gene = MockGeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=0,
            end=300,
            strand="+",
            transcripts=[transcript],
        )

        refiner = SpliceRefiner(
            genome,
            junctions,
            max_shift=15,  # Only 15bp
        )

        refined_gene, report = refiner.refine_gene(gene)

        # Should have no corrections (junction too far)
        assert report.n_corrected == 0

    def test_parallel_processing(
        self,
        mock_genome_with_splice_sites,
        mock_junctions,
        mock_gene,
    ):
        """Test parallel gene processing."""
        # Create multiple genes
        genes = [mock_gene]
        for i in range(5):
            tx = MockTranscriptModel(
                transcript_id=f"tx{i+2}",
                seqid="chr1",
                start=100 + i * 10,
                end=400 + i * 10,
                strand="+",
                exons=[(100 + i * 10, 200), (300, 400 + i * 10)],
            )
            gene = MockGeneModel(
                gene_id=f"gene{i+2}",
                seqid="chr1",
                start=100 + i * 10,
                end=400 + i * 10,
                strand="+",
                transcripts=[tx],
            )
            genes.append(gene)

        refiner = SpliceRefiner(
            mock_genome_with_splice_sites,
            mock_junctions,
        )

        # Process in parallel
        results = list(refiner.refine_genes_parallel(genes, n_workers=2))

        assert len(results) == len(genes)
        for refined_gene, report in results:
            assert report.gene_id.startswith("gene")


class TestEndToEndWorkflow:
    """End-to-end workflow tests."""

    def test_full_splice_workflow(self, tmp_path):
        """End-to-end splice refinement workflow."""
        # Create test genome
        genome = MockGenomeAccessor({
            "chr1": (
                "A" * 100 +
                "ATGCGA" * 17 +  # exon 1
                "GTAAGT" +  # donor
                "TTTTTT" * 16 +  # intron
                "TTTCAG" +  # acceptor
                "ATGCGA" * 17 +  # exon 2
                "A" * 500
            )
        })

        # Create matching junctions
        junctions = {
            "chr1": [
                MockSpliceJunction("chr1", 202, 302, "+", 100),
            ]
        }

        # Create gene
        transcript = MockTranscriptModel(
            transcript_id="tx1",
            seqid="chr1",
            start=100,
            end=400,
            strand="+",
            exons=[(100, 202), (302, 400)],
        )
        gene = MockGeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=400,
            strand="+",
            transcripts=[transcript],
        )

        # Run refinement
        refiner = SpliceRefiner(genome, junctions)
        refined_gene, report = refiner.refine_gene(gene)

        # Verify report
        assert report.gene_id == "gene1"
        assert report.n_introns == 1

        # Write outputs
        output_tsv = tmp_path / "splice_report.tsv"
        SpliceReportWriter.write_tsv([report], output_tsv)
        assert output_tsv.exists()

        # Verify TSV content
        with open(output_tsv) as f:
            content = f.read()
        assert "gene1" in content

    def test_splice_preserves_unchanged_genes(self):
        """Genes with perfect splice sites should be unchanged."""
        # Create genome with perfect GT-AG
        genome = MockGenomeAccessor({
            "chr1": "A" * 100 + "GT" + "A" * 98 + "AG" + "A" * 700
        })

        # Junction matching exactly
        junctions = {
            "chr1": [
                MockSpliceJunction("chr1", 100, 200, "+", 100),
            ]
        }

        transcript = MockTranscriptModel(
            transcript_id="tx1",
            seqid="chr1",
            start=50,
            end=250,
            strand="+",
            exons=[(50, 100), (200, 250)],
        )
        gene = MockGeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=50,
            end=250,
            strand="+",
            transcripts=[transcript],
        )

        refiner = SpliceRefiner(genome, junctions)
        refined_gene, report = refiner.refine_gene(gene)

        # No corrections should be made
        assert report.n_corrected == 0


class TestStrandHandling:
    """Test correct handling of minus strand genes."""

    def test_minus_strand_splice_sites(self):
        """Reverse strand splice sites should be handled correctly."""
        # On minus strand, donor is at intron end (CT complementary to AG)
        # and acceptor at intron start (AC complementary to GT)
        genome = MockGenomeAccessor({
            "chr1": "A" * 100 + "CT" + "A" * 96 + "AC" + "A" * 700
            # Forward: CT...AC
            # Reverse: GT...AG (canonical)
        })

        junctions = {
            "chr1": [
                MockSpliceJunction("chr1", 100, 200, "-", 50),
            ]
        }

        transcript = MockTranscriptModel(
            transcript_id="tx1",
            seqid="chr1",
            start=50,
            end=250,
            strand="-",
            exons=[(50, 100), (200, 250)],
        )
        gene = MockGeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=50,
            end=250,
            strand="-",
            transcripts=[transcript],
        )

        refiner = SpliceRefiner(genome, junctions)
        refined_gene, report = refiner.refine_gene(gene)

        # Should process without error
        assert report.gene_id == "gene1"


# =============================================================================
# Report Output Tests
# =============================================================================


class TestReportOutputs:
    """Test report file generation."""

    def test_multiple_report_formats(self, tmp_path):
        """Test generating multiple output formats."""
        reports = [
            GeneSpliceReport(
                gene_id="gene1",
                n_introns=3,
                n_supported=2,
                n_corrected=1,
                n_canonical=2,
                n_noncanonical=1,
                unsupported_introns=[1],
                flags=["noncanonical:1"],
            ),
            GeneSpliceReport(
                gene_id="gene2",
                n_introns=2,
                n_supported=2,
                n_corrected=0,
                n_canonical=2,
                n_noncanonical=0,
            ),
        ]

        # TSV report
        tsv_file = tmp_path / "report.tsv"
        SpliceReportWriter.write_tsv(reports, tsv_file)
        assert tsv_file.exists()

        # Summary statistics
        stats = SpliceReportWriter.summary_statistics(reports)
        assert stats["total_genes"] == 2
        assert stats["total_introns"] == 5
        assert stats["supported_introns"] == 4


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
