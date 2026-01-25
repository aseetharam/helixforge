"""Unit tests for helixforge.io.bam module.

Tests cover:
- SpliceJunction and CoverageProfile data structures
- JunctionExtractor initialization
- CIGAR parsing for junction extraction
- Coverage calculation

Note: These tests use mocking since creating real BAM files requires
complex setup. Integration tests will use real BAM files.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from helixforge.io.bam import (
    CANONICAL_ACCEPTORS,
    CANONICAL_DONORS,
    CIGAR_D,
    CIGAR_I,
    CIGAR_M,
    CIGAR_N,
    CIGAR_S,
    CoverageProfile,
    JunctionExtractor,
    SpliceJunction,
    is_canonical_junction,
)


# =============================================================================
# Data Structure Tests
# =============================================================================


class TestSpliceJunction:
    """Tests for SpliceJunction named tuple."""

    def test_creation(self) -> None:
        """Test creating a SpliceJunction."""
        junction = SpliceJunction(
            seqid="chr1",
            start=100,
            end=200,
            strand="+",
            read_count=10,
            unique_count=8,
            is_canonical=True,
        )

        assert junction.seqid == "chr1"
        assert junction.start == 100
        assert junction.end == 200
        assert junction.strand == "+"
        assert junction.read_count == 10
        assert junction.unique_count == 8
        assert junction.is_canonical is True

    def test_default_canonical(self) -> None:
        """Test default is_canonical value."""
        junction = SpliceJunction(
            seqid="chr1",
            start=100,
            end=200,
            strand="+",
            read_count=10,
            unique_count=8,
        )

        assert junction.is_canonical is None

    def test_immutable(self) -> None:
        """Test that SpliceJunction is immutable."""
        junction = SpliceJunction(
            seqid="chr1",
            start=100,
            end=200,
            strand="+",
            read_count=10,
            unique_count=8,
        )

        with pytest.raises(AttributeError):
            junction.start = 150  # type: ignore


class TestCoverageProfile:
    """Tests for CoverageProfile named tuple."""

    def test_creation(self) -> None:
        """Test creating a CoverageProfile."""
        coverage = np.array([10, 15, 20, 15, 10], dtype=np.int32)
        profile = CoverageProfile(
            seqid="chr1",
            start=100,
            end=105,
            strand="+",
            values=coverage,
        )

        assert profile.seqid == "chr1"
        assert profile.start == 100
        assert profile.end == 105
        assert len(profile.values) == 5

    def test_length_property(self) -> None:
        """Test coverage array length matches region."""
        coverage = np.array([10, 15, 20, 15, 10], dtype=np.int32)
        profile = CoverageProfile(
            seqid="chr1",
            start=100,
            end=105,
            strand="+",
            values=coverage,
        )

        assert len(profile.values) == profile.end - profile.start


# =============================================================================
# CIGAR Parsing Tests
# =============================================================================


class TestCIGARConstants:
    """Tests for CIGAR operation constants."""

    def test_cigar_constants(self) -> None:
        """Test CIGAR operation constants match pysam values."""
        assert CIGAR_M == 0  # Match
        assert CIGAR_I == 1  # Insertion
        assert CIGAR_D == 2  # Deletion
        assert CIGAR_N == 3  # Skipped (intron)
        assert CIGAR_S == 4  # Soft clip


class TestIsCanonicalJunction:
    """Tests for is_canonical_junction function."""

    def test_gt_ag_canonical(self) -> None:
        """Test GT-AG is canonical."""
        assert is_canonical_junction("GT", "AG") is True

    def test_gc_ag_canonical(self) -> None:
        """Test GC-AG is canonical."""
        assert is_canonical_junction("GC", "AG") is True

    def test_at_ac_canonical(self) -> None:
        """Test AT-AC (minor spliceosome) is canonical."""
        assert is_canonical_junction("AT", "AC") is True

    def test_non_canonical(self) -> None:
        """Test non-canonical combinations."""
        assert is_canonical_junction("AA", "AA") is False
        assert is_canonical_junction("GT", "AC") is False
        assert is_canonical_junction("CT", "AG") is False

    def test_case_insensitive(self) -> None:
        """Test case insensitivity."""
        assert is_canonical_junction("gt", "ag") is True
        assert is_canonical_junction("Gt", "Ag") is True


# =============================================================================
# Constants Tests
# =============================================================================


class TestCanonicalSpliceSites:
    """Tests for canonical splice site constants."""

    def test_canonical_donors(self) -> None:
        """Test canonical donor sites."""
        assert "GT" in CANONICAL_DONORS
        assert "GC" in CANONICAL_DONORS

    def test_canonical_acceptors(self) -> None:
        """Test canonical acceptor sites."""
        assert "AG" in CANONICAL_ACCEPTORS


# =============================================================================
# JunctionExtractor Tests with Mocking
# =============================================================================


class TestJunctionExtractorInit:
    """Tests for JunctionExtractor initialization."""

    def test_init_missing_file(self, tmp_path: Path) -> None:
        """Test error for missing BAM file."""
        with pytest.raises(FileNotFoundError, match="BAM file not found"):
            JunctionExtractor(tmp_path / "missing.bam")

    def test_init_missing_index(self, tmp_path: Path) -> None:
        """Test error for missing index file."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()

        with pytest.raises(ValueError, match="BAM index not found"):
            JunctionExtractor(bam_path)

    @patch("helixforge.io.bam.pysam.AlignmentFile")
    def test_init_success(self, mock_alignment_file: MagicMock, tmp_path: Path) -> None:
        """Test successful initialization with mock."""
        # Create empty BAM file and index to pass existence check
        bam_path = tmp_path / "test.bam"
        bam_path.touch()
        (tmp_path / "test.bam.bai").touch()

        # Mock pysam.AlignmentFile
        mock_file = MagicMock()
        mock_file.references = ["chr1", "chr2"]
        mock_file.lengths = [1000, 500]
        mock_alignment_file.return_value = mock_file

        extractor = JunctionExtractor(bam_path)

        assert extractor.path == bam_path
        assert "chr1" in extractor.reference_lengths
        assert extractor.reference_lengths["chr1"] == 1000

        extractor.close()

    @patch("helixforge.io.bam.pysam.AlignmentFile")
    def test_context_manager(
        self, mock_alignment_file: MagicMock, tmp_path: Path
    ) -> None:
        """Test context manager usage."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()
        (tmp_path / "test.bam.bai").touch()

        mock_file = MagicMock()
        mock_file.references = ["chr1"]
        mock_file.lengths = [1000]
        mock_alignment_file.return_value = mock_file

        with JunctionExtractor(bam_path) as extractor:
            assert extractor.path == bam_path


class TestJunctionExtractorExtraction:
    """Tests for junction extraction functionality."""

    @pytest.fixture
    def mock_extractor(self, tmp_path: Path):
        """Create a JunctionExtractor with mocked BAM file."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()
        # Create fake index
        (tmp_path / "test.bam.bai").touch()

        with patch("helixforge.io.bam.pysam.AlignmentFile") as mock_alignment_file:
            mock_file = MagicMock()
            mock_file.references = ["chr1", "chr2"]
            mock_file.lengths = [1000, 500]

            # Create mock reads
            def create_mock_read(
                name: str,
                start: int,
                cigar: list,
                is_reverse: bool = False,
                mapq: int = 60,
            ):
                read = MagicMock()
                read.query_name = name
                read.reference_start = start
                read.cigartuples = cigar
                read.is_reverse = is_reverse
                read.mapping_quality = mapq
                read.is_unmapped = False
                read.is_secondary = False
                read.is_supplementary = False
                read.is_paired = False
                read.has_tag.return_value = False
                read.get_tag.side_effect = KeyError("XS")
                return read

            # Mock reads with splice junctions
            mock_reads = [
                # Read 1: junction at 150-250
                create_mock_read("read1", 100, [(0, 50), (3, 100), (0, 50)]),
                # Read 2: same junction
                create_mock_read("read2", 110, [(0, 40), (3, 100), (0, 60)]),
                # Read 3: different junction
                create_mock_read("read3", 400, [(0, 30), (3, 50), (0, 30)], is_reverse=True),
                # Read 4: no junction
                create_mock_read("read4", 100, [(0, 100)]),
            ]

            mock_file.fetch.return_value = iter(mock_reads)
            mock_alignment_file.return_value = mock_file

            extractor = JunctionExtractor(bam_path)
            yield extractor
            extractor.close()

    def test_extract_region(self, mock_extractor: JunctionExtractor) -> None:
        """Test extracting junctions from a region."""
        junctions = mock_extractor.extract_region("chr1", 0, 1000)

        # Should find junctions from the mock reads
        assert isinstance(junctions, list)


class TestJunctionExtractorCoverage:
    """Tests for coverage calculation."""

    @pytest.fixture
    def coverage_mock_extractor(self, tmp_path: Path):
        """Create extractor with mock for coverage testing."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()
        (tmp_path / "test.bam.bai").touch()

        with patch("helixforge.io.bam.pysam.AlignmentFile") as mock_alignment_file:
            mock_file = MagicMock()
            mock_file.references = ["chr1"]
            mock_file.lengths = [1000]

            # Mock fetch to return reads with coverage
            def create_mock_read(start: int, end: int, mapq: int = 60):
                read = MagicMock()
                read.reference_start = start
                read.is_unmapped = False
                read.is_secondary = False
                read.is_supplementary = False
                read.is_reverse = False
                read.mapping_quality = mapq
                read.get_blocks.return_value = [(start, end)]
                read.get_tag.side_effect = KeyError("XS")
                return read

            # Create reads covering positions 100-105
            mock_reads = [
                create_mock_read(100, 110),
                create_mock_read(98, 108),
            ]
            mock_file.fetch.return_value = iter(mock_reads)

            mock_alignment_file.return_value = mock_file

            extractor = JunctionExtractor(bam_path)
            yield extractor
            extractor.close()

    def test_get_coverage(
        self, coverage_mock_extractor: JunctionExtractor
    ) -> None:
        """Test coverage calculation."""
        profile = coverage_mock_extractor.get_coverage("chr1", 100, 105)

        assert profile.seqid == "chr1"
        assert profile.start == 100
        assert profile.end == 105
        assert len(profile.values) == 5


class TestJunctionExtractorBED:
    """Tests for BED output functionality."""

    @pytest.fixture
    def bed_mock_extractor(self, tmp_path: Path):
        """Create extractor with mock for BED testing."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()
        (tmp_path / "test.bam.bai").touch()

        with patch("helixforge.io.bam.pysam.AlignmentFile") as mock_alignment_file:
            mock_file = MagicMock()
            mock_file.references = ["chr1"]
            mock_file.lengths = [1000]
            mock_alignment_file.return_value = mock_file

            extractor = JunctionExtractor(bam_path)
            yield extractor
            extractor.close()

    def test_to_bed(self, bed_mock_extractor: JunctionExtractor, tmp_path: Path) -> None:
        """Test writing junctions to BED format."""
        output_path = tmp_path / "junctions.bed"

        junctions = [
            SpliceJunction("chr1", 150, 250, "+", 10, 8, True),
            SpliceJunction("chr1", 430, 480, "-", 5, 5, True),
        ]

        # Call instance method
        bed_mock_extractor.to_bed(junctions, output_path)

        # Verify output
        content = output_path.read_text()
        lines = [l for l in content.strip().split("\n") if not l.startswith("#")]

        assert len(lines) == 2

        # Check first line format
        parts = lines[0].split("\t")
        assert parts[0] == "chr1"
        assert parts[1] == "150"  # start
        assert parts[2] == "250"  # end
        # Score is min(1000, read_count)
        assert parts[4] == "10"
        assert parts[5] == "+"    # strand

    def test_to_bed_empty(self, bed_mock_extractor: JunctionExtractor, tmp_path: Path) -> None:
        """Test writing empty junction list."""
        output_path = tmp_path / "empty.bed"

        bed_mock_extractor.to_bed([], output_path)

        content = output_path.read_text()
        # Should just have the header line
        assert "#chrom" in content or content.strip() == ""


# =============================================================================
# Integration-like Tests (still using mocks)
# =============================================================================


class TestJunctionExtractorWorkflow:
    """Tests for typical usage workflow."""

    @patch("helixforge.io.bam.pysam.AlignmentFile")
    def test_full_workflow(
        self, mock_alignment_file: MagicMock, tmp_path: Path
    ) -> None:
        """Test complete workflow: open, extract, write BED."""
        bam_path = tmp_path / "test.bam"
        bam_path.touch()
        (tmp_path / "test.bam.bai").touch()
        bed_path = tmp_path / "junctions.bed"

        mock_file = MagicMock()
        mock_file.references = ["chr1"]
        mock_file.lengths = [1000]
        mock_file.fetch.return_value = iter([])  # No reads
        mock_alignment_file.return_value = mock_file

        with JunctionExtractor(bam_path) as extractor:
            junctions = extractor.extract_region("chr1", 0, 1000)
            extractor.to_bed(junctions, bed_path)

        assert bed_path.exists()
