"""Tests for evidence output writers."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from helixforge.core.evidence import (
    EvidenceLevel,
    EvidenceScore,
    ExonEvidence,
    JunctionEvidence,
)
from helixforge.core.evidence_output import (
    generate_evidence_summary_text,
    load_evidence_report_tsv,
    write_evidence_report_tsv,
    write_exon_details_tsv,
    write_gene_evidence_summary_tsv,
    write_junction_details_tsv,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def sample_scores():
    """Create sample EvidenceScore objects for testing."""
    return [
        EvidenceScore(
            gene_id="gene1",
            transcript_id="tx1",
            seqid="chr1",
            strand="+",
            n_introns=3,
            n_introns_supported=3,
            n_introns_exact=2,
            junction_support_ratio=1.0,
            n_exons=4,
            n_exons_expressed=4,
            mean_exon_coverage=50.0,
            exon_coverage_ratio=1.0,
            start_supported=True,
            stop_supported=True,
            strand_consistent=True,
            evidence_level=EvidenceLevel.FULL,
            aed=0.0,
            junction_evidence=[
                JunctionEvidence(1000, 2000, True, True, 50, 45, True),
                JunctionEvidence(2500, 3500, True, True, 30, 28, True),
                JunctionEvidence(4000, 5000, True, True, 20, 18, False, 2, -1),
            ],
            exon_evidence=[
                ExonEvidence(100, 1000, 60.0, 55.0, 10, 100, 0.95, 0.3),
                ExonEvidence(2000, 2500, 45.0, 42.0, 5, 80, 0.90, 0.35),
                ExonEvidence(3500, 4000, 50.0, 48.0, 8, 90, 0.92, 0.28),
                ExonEvidence(5000, 5500, 40.0, 38.0, 3, 70, 0.85, 0.40),
            ],
            flags=[],
        ),
        EvidenceScore(
            gene_id="gene2",
            transcript_id="tx2",
            seqid="chr1",
            strand="-",
            n_introns=2,
            n_introns_supported=1,
            n_introns_exact=1,
            junction_support_ratio=0.5,
            n_exons=3,
            n_exons_expressed=2,
            mean_exon_coverage=25.0,
            exon_coverage_ratio=0.67,
            start_supported=True,
            stop_supported=False,
            strand_consistent=True,
            evidence_level=EvidenceLevel.PARTIAL,
            aed=0.35,
            junction_evidence=[
                JunctionEvidence(10000, 11000, True, True, 40, 38, True),
                JunctionEvidence(12000, 13000, True, False, 0, 0, False),
            ],
            exon_evidence=[
                ExonEvidence(9000, 10000, 30.0, 28.0, 5, 60, 0.88, 0.32),
                ExonEvidence(11000, 12000, 25.0, 22.0, 2, 50, 0.75, 0.45),
                ExonEvidence(13000, 14000, 20.0, 18.0, 1, 40, 0.60, 0.50),
            ],
            flags=["PARTIAL_JUNCTION_SUPPORT", "STOP_UNSUPPORTED"],
        ),
        EvidenceScore(
            gene_id="gene3",
            transcript_id="tx3",
            seqid="chr2",
            strand="+",
            n_introns=0,
            n_introns_supported=0,
            n_introns_exact=0,
            junction_support_ratio=1.0,  # Single-exon
            n_exons=1,
            n_exons_expressed=0,
            mean_exon_coverage=2.0,
            exon_coverage_ratio=0.0,
            start_supported=False,
            stop_supported=False,
            strand_consistent=True,
            evidence_level=EvidenceLevel.NONE,
            aed=0.85,
            junction_evidence=[],
            exon_evidence=[
                ExonEvidence(500, 1500, 2.0, 1.0, 0, 5, 0.1, 0.8),
            ],
            flags=["NO_EXON_EXPRESSION", "START_UNSUPPORTED", "STOP_UNSUPPORTED"],
        ),
    ]


# =============================================================================
# Test Report Writers
# =============================================================================


class TestWriteEvidenceReportTsv:
    """Test write_evidence_report_tsv function."""

    def test_write_report(self, sample_scores):
        """Test writing evidence report TSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_evidence_report_tsv(sample_scores, output_path)

            # Verify file exists and has content
            assert output_path.exists()
            lines = output_path.read_text().strip().split("\n")

            # Should have header + 3 data lines
            assert len(lines) == 4

            # Check header
            header = lines[0].split("\t")
            assert "gene_id" in header
            assert "aed" in header
            assert "evidence_level" in header

            # Check first data line
            data = lines[1].split("\t")
            assert data[0] == "gene1"

        finally:
            output_path.unlink(missing_ok=True)

    def test_write_empty_report(self):
        """Test writing empty report."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_evidence_report_tsv([], output_path)

            lines = output_path.read_text().strip().split("\n")
            # Should have header only
            assert len(lines) == 1

        finally:
            output_path.unlink(missing_ok=True)


class TestWriteGeneSummaryTsv:
    """Test write_gene_evidence_summary_tsv function."""

    def test_write_summary(self, sample_scores):
        """Test writing gene summary TSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_gene_evidence_summary_tsv(sample_scores, output_path)

            assert output_path.exists()
            lines = output_path.read_text().strip().split("\n")

            # Header + 3 data lines
            assert len(lines) == 4

            # Verify columns
            header = lines[0].split("\t")
            assert "gene_id" in header
            assert "evidence_level" in header
            assert "aed" in header
            assert "n_flags" in header

        finally:
            output_path.unlink(missing_ok=True)


class TestWriteJunctionDetailsTsv:
    """Test write_junction_details_tsv function."""

    def test_write_junction_details(self, sample_scores):
        """Test writing junction details TSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_junction_details_tsv(sample_scores, output_path)

            assert output_path.exists()
            lines = output_path.read_text().strip().split("\n")

            # Header + 5 junctions (3 from gene1, 2 from gene2, 0 from gene3)
            assert len(lines) == 6

            header = lines[0].split("\t")
            assert "gene_id" in header
            assert "intron_start" in header
            assert "read_count" in header
            assert "exact_match" in header

        finally:
            output_path.unlink(missing_ok=True)


class TestWriteExonDetailsTsv:
    """Test write_exon_details_tsv function."""

    def test_write_exon_details(self, sample_scores):
        """Test writing exon details TSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_exon_details_tsv(sample_scores, output_path)

            assert output_path.exists()
            lines = output_path.read_text().strip().split("\n")

            # Header + 8 exons (4 from gene1, 3 from gene2, 1 from gene3)
            assert len(lines) == 9

            header = lines[0].split("\t")
            assert "gene_id" in header
            assert "exon_start" in header
            assert "mean_coverage" in header
            assert "is_expressed" in header

        finally:
            output_path.unlink(missing_ok=True)


# =============================================================================
# Test Loading
# =============================================================================


class TestLoadEvidenceReportTsv:
    """Test load_evidence_report_tsv function."""

    def test_roundtrip(self, sample_scores):
        """Test writing and loading back evidence report."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            # Write
            write_evidence_report_tsv(sample_scores, output_path)

            # Load back
            loaded_scores = load_evidence_report_tsv(output_path)

            # Verify
            assert len(loaded_scores) == len(sample_scores)

            for original, loaded in zip(sample_scores, loaded_scores):
                assert loaded.gene_id == original.gene_id
                assert loaded.transcript_id == original.transcript_id
                assert loaded.evidence_level == original.evidence_level
                assert loaded.aed == pytest.approx(original.aed, abs=0.0001)
                assert loaded.n_introns == original.n_introns
                assert loaded.n_exons == original.n_exons

        finally:
            output_path.unlink(missing_ok=True)


# =============================================================================
# Test Summary Text
# =============================================================================


class TestGenerateSummaryText:
    """Test generate_evidence_summary_text function."""

    def test_summary_text(self, sample_scores):
        """Test generating summary text."""
        text = generate_evidence_summary_text(sample_scores)

        assert "Evidence Score Summary" in text
        assert "Total genes:" in text
        assert "3" in text  # 3 genes
        assert "Full support:" in text
        assert "Mean AED:" in text

    def test_summary_empty(self):
        """Test summary for empty scores."""
        text = generate_evidence_summary_text([])
        assert "No genes to summarize" in text

    def test_summary_single_exon_stats(self, sample_scores):
        """Test single-exon gene stats in summary."""
        text = generate_evidence_summary_text(sample_scores)

        assert "Single-exon genes:" in text
        assert "Multi-exon genes:" in text


# =============================================================================
# Test Edge Cases
# =============================================================================


class TestEdgeCases:
    """Test edge cases for output writers."""

    def test_empty_flags(self):
        """Test score with empty flags."""
        scores = [
            EvidenceScore(
                gene_id="gene1",
                transcript_id="tx1",
                seqid="chr1",
                strand="+",
                evidence_level=EvidenceLevel.FULL,
                aed=0.0,
                flags=[],
            )
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_evidence_report_tsv(scores, output_path)
            loaded = load_evidence_report_tsv(output_path)

            assert loaded[0].flags == []

        finally:
            output_path.unlink(missing_ok=True)

    def test_multiple_flags(self):
        """Test score with multiple flags."""
        scores = [
            EvidenceScore(
                gene_id="gene1",
                transcript_id="tx1",
                seqid="chr1",
                strand="+",
                evidence_level=EvidenceLevel.MINIMAL,
                aed=0.8,
                flags=["FLAG1", "FLAG2", "FLAG3"],
            )
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_evidence_report_tsv(scores, output_path)
            loaded = load_evidence_report_tsv(output_path)

            assert len(loaded[0].flags) == 3
            assert "FLAG1" in loaded[0].flags
            assert "FLAG2" in loaded[0].flags
            assert "FLAG3" in loaded[0].flags

        finally:
            output_path.unlink(missing_ok=True)

    def test_special_characters_in_gene_id(self):
        """Test handling gene IDs with special characters."""
        scores = [
            EvidenceScore(
                gene_id="gene.1-a_special",
                transcript_id="tx.1-a_special",
                seqid="chr1",
                strand="+",
                evidence_level=EvidenceLevel.FULL,
                aed=0.1,
            )
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            output_path = Path(f.name)

        try:
            write_evidence_report_tsv(scores, output_path)
            loaded = load_evidence_report_tsv(output_path)

            assert loaded[0].gene_id == "gene.1-a_special"

        finally:
            output_path.unlink(missing_ok=True)
