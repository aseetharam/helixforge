"""Tests for the refine pipeline module."""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from helixforge.core.refine import (
    RefineConfig,
    RefinedGene,
    RefineReportWriter,
    DEFAULT_MAX_SHIFT,
    DEFAULT_MIN_JUNCTION_READS,
    DEFAULT_MIN_TISSUES,
    DEFAULT_BOUNDARY_WINDOW,
    DEFAULT_CONFIDENCE_THRESHOLD,
)


# =============================================================================
# Test RefineConfig
# =============================================================================


class TestRefineConfig:
    """Test RefineConfig data structure."""

    def test_default_values(self):
        """Test default configuration values."""
        config = RefineConfig()

        assert config.max_shift == DEFAULT_MAX_SHIFT
        assert config.min_junction_reads == DEFAULT_MIN_JUNCTION_READS
        assert config.min_tissues == DEFAULT_MIN_TISSUES
        assert config.boundary_search_window == DEFAULT_BOUNDARY_WINDOW
        assert config.confidence_threshold == DEFAULT_CONFIDENCE_THRESHOLD
        assert config.adjust_boundaries is True
        assert config.stranded is False
        assert config.skip_coverage is False

    def test_custom_values(self):
        """Test custom configuration values."""
        config = RefineConfig(
            max_shift=20,
            min_junction_reads=5,
            min_tissues=2,
            adjust_boundaries=False,
            confidence_threshold=0.8,
            min_exon_coverage=10,
            boundary_tolerance=15,
            skip_coverage=True,
        )

        assert config.max_shift == 20
        assert config.min_junction_reads == 5
        assert config.min_tissues == 2
        assert config.adjust_boundaries is False
        assert config.confidence_threshold == 0.8
        assert config.min_exon_coverage == 10
        assert config.boundary_tolerance == 15
        assert config.skip_coverage is True


# =============================================================================
# Test RefinedGene
# =============================================================================


class TestRefinedGene:
    """Test RefinedGene data structure."""

    @pytest.fixture
    def mock_gene(self):
        """Create a mock gene model."""
        gene = MagicMock()
        gene.gene_id = "gene1"
        gene.seqid = "chr1"
        gene.start = 1000
        gene.end = 5000
        gene.strand = "+"
        gene.attributes = {}

        transcript = MagicMock()
        transcript.exons = [MagicMock(), MagicMock(), MagicMock()]
        gene.transcripts = [transcript]

        return gene

    @pytest.fixture
    def mock_confidence(self):
        """Create a mock confidence result."""
        conf = MagicMock()
        conf.overall_score = 0.85
        conf.mean_confidence = 0.85
        conf.worst_exon_score = 0.75
        conf.flags = []
        return conf

    @pytest.fixture
    def mock_evidence(self):
        """Create a mock evidence result."""
        ev = MagicMock()
        ev.aed = 0.1
        ev.n_introns = 2
        ev.n_introns_supported = 2
        ev.junction_support_ratio = 1.0
        ev.mean_exon_coverage = 50.0
        ev.flags = []
        ev.evidence_level = MagicMock()
        ev.evidence_level.value = "full"
        return ev

    def test_refined_gene_creation(self, mock_gene):
        """Test creating a RefinedGene instance."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
            splice_corrections=2,
            boundary_adjusted=True,
            flags=["SPLICE_CORRECTED"],
        )

        assert refined.gene.gene_id == "gene1"
        assert refined.splice_corrections == 2
        assert refined.boundary_adjusted is True
        assert "SPLICE_CORRECTED" in refined.flags

    def test_confidence_score_property(self, mock_gene, mock_confidence):
        """Test confidence_score property."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
            confidence=mock_confidence,
        )

        assert refined.confidence_score == 0.85

    def test_confidence_score_none(self, mock_gene):
        """Test confidence_score when no confidence is set."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
        )

        assert refined.confidence_score is None

    def test_evidence_score_property(self, mock_gene, mock_evidence):
        """Test evidence_score property."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
            evidence=mock_evidence,
        )

        # evidence_score = 1 - AED
        assert refined.evidence_score == 0.9

    def test_aed_property(self, mock_gene, mock_evidence):
        """Test aed property."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
            evidence=mock_evidence,
        )

        assert refined.aed == 0.1

    def test_junction_support_str(self, mock_gene, mock_evidence):
        """Test junction_support_str property."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
            evidence=mock_evidence,
        )

        assert refined.junction_support_str == "2/2"

    def test_junction_support_str_no_evidence(self, mock_gene):
        """Test junction_support_str when no evidence."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
        )

        assert refined.junction_support_str == "N/A"

    def test_to_report_dict(self, mock_gene, mock_confidence, mock_evidence):
        """Test to_report_dict method."""
        refined = RefinedGene(
            gene=mock_gene,
            original_gene=mock_gene,
            confidence=mock_confidence,
            evidence=mock_evidence,
            splice_corrections=1,
            boundary_adjusted=True,
            flags=["SPLICE_CORRECTED", "BOUNDARY_ADJUSTED"],
        )

        report = refined.to_report_dict()

        assert report["gene_id"] == "gene1"
        assert report["seqid"] == "chr1"
        assert report["start"] == 1000
        assert report["end"] == 5000
        assert report["strand"] == "+"
        assert report["n_exons"] == 3
        assert report["confidence_score"] == 0.85
        assert report["aed"] == 0.1
        assert report["splice_corrections"] == 1
        assert report["boundary_adjusted"] is True
        assert "SPLICE_CORRECTED" in report["flags"]


# =============================================================================
# Test RefineReportWriter
# =============================================================================


class TestRefineReportWriter:
    """Test RefineReportWriter class."""

    @pytest.fixture
    def mock_results(self):
        """Create mock refined gene results."""
        results = []

        for i in range(3):
            gene = MagicMock()
            gene.gene_id = f"gene{i+1}"
            gene.seqid = "chr1"
            gene.start = i * 10000
            gene.end = i * 10000 + 5000
            gene.strand = "+"
            gene.attributes = {}

            transcript = MagicMock()
            transcript.exons = [MagicMock()] * (i + 2)
            gene.transcripts = [transcript]

            confidence = MagicMock()
            confidence.overall_score = 0.8 + i * 0.05
            confidence.mean_confidence = 0.8 + i * 0.05
            confidence.worst_exon_score = 0.7
            confidence.flags = []

            evidence = MagicMock()
            evidence.aed = 0.1 + i * 0.1
            evidence.n_introns = i + 1
            evidence.n_introns_supported = i + 1
            evidence.junction_support_ratio = 1.0
            evidence.mean_exon_coverage = 50.0
            evidence.flags = []
            evidence.evidence_level = MagicMock()
            evidence.evidence_level.value = "full"

            refined = RefinedGene(
                gene=gene,
                original_gene=gene,
                confidence=confidence,
                evidence=evidence,
                splice_corrections=i,
                boundary_adjusted=i > 0,
                flags=["SPLICE_CORRECTED"] if i > 0 else [],
            )
            results.append(refined)

        return results

    def test_write_report(self, mock_results):
        """Test writing report TSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            RefineReportWriter.write_report(mock_results, output_path)

            assert output_path.exists()
            lines = output_path.read_text().strip().split("\n")

            # Header + 3 data lines
            assert len(lines) == 4

            # Check header
            header = lines[0].split("\t")
            assert "gene_id" in header
            assert "confidence_score" in header
            assert "aed" in header

            # Check first data line
            data = lines[1].split("\t")
            assert data[0] == "gene1"

        finally:
            output_path.unlink(missing_ok=True)

    def test_write_summary(self, mock_results):
        """Test generating summary statistics."""
        summary = RefineReportWriter.write_summary(mock_results)

        assert summary["total_genes"] == 3
        assert summary["mean_confidence"] is not None
        assert summary["mean_aed"] is not None
        assert summary["genes_splice_corrected"] == 2  # gene2 and gene3
        assert summary["boundary_adjusted"] == 2

    def test_write_summary_empty(self):
        """Test summary for empty results."""
        summary = RefineReportWriter.write_summary([])

        assert summary["total_genes"] == 0

    def test_format_summary(self, mock_results):
        """Test formatting summary as text."""
        summary = RefineReportWriter.write_summary(mock_results)
        text = RefineReportWriter.format_summary(summary)

        assert "REFINEMENT SUMMARY" in text
        assert "Total genes:" in text
        assert "3" in text


# =============================================================================
# Test Edge Cases
# =============================================================================


class TestEdgeCases:
    """Test edge cases for refine module."""

    def test_refined_gene_no_transcripts(self):
        """Test RefinedGene with gene having no transcripts."""
        gene = MagicMock()
        gene.gene_id = "gene1"
        gene.seqid = "chr1"
        gene.start = 1000
        gene.end = 2000
        gene.strand = "+"
        gene.attributes = {}
        gene.transcripts = []

        refined = RefinedGene(
            gene=gene,
            original_gene=gene,
        )

        report = refined.to_report_dict()
        assert report["n_exons"] == 0

    def test_refined_gene_single_exon(self):
        """Test RefinedGene with single-exon gene."""
        gene = MagicMock()
        gene.gene_id = "gene1"
        gene.seqid = "chr1"
        gene.start = 1000
        gene.end = 2000
        gene.strand = "+"
        gene.attributes = {}

        transcript = MagicMock()
        transcript.exons = [MagicMock()]
        gene.transcripts = [transcript]

        evidence = MagicMock()
        evidence.aed = 0.2
        evidence.n_introns = 0
        evidence.n_introns_supported = 0
        evidence.junction_support_ratio = 1.0  # Single-exon default
        evidence.mean_exon_coverage = 30.0
        evidence.flags = []
        evidence.evidence_level = MagicMock()
        evidence.evidence_level.value = "partial"

        refined = RefinedGene(
            gene=gene,
            original_gene=gene,
            evidence=evidence,
        )

        # Single-exon gene should show N/A for junction support
        assert refined.junction_support_str == "N/A"

    def test_config_boundary_edge_values(self):
        """Test config with edge case values."""
        config = RefineConfig(
            max_shift=0,  # No shift allowed
            min_junction_reads=1,  # Minimum reads
            min_tissues=1,
            boundary_search_window=0,  # No search
        )

        assert config.max_shift == 0
        assert config.min_junction_reads == 1
        assert config.boundary_search_window == 0
