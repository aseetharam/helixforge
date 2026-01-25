"""Tests for region parsing utilities.

Tests the genomic region parsing and validation functions used
for chunk-aware processing in parallel execution.
"""

import pytest

from helixforge.utils.regions import (
    GenomicRegion,
    format_region_for_cli,
    parse_region,
    region_from_scaffold,
    region_to_str,
    validate_region,
)


# =============================================================================
# Test GenomicRegion
# =============================================================================


class TestGenomicRegion:
    """Tests for GenomicRegion NamedTuple."""

    def test_creation(self):
        """Create a GenomicRegion."""
        region = GenomicRegion("chr1", 100, 200)
        assert region.seqid == "chr1"
        assert region.start == 100
        assert region.end == 200

    def test_length(self):
        """Test length property."""
        region = GenomicRegion("chr1", 100, 200)
        assert region.length == 100

    def test_length_zero(self):
        """Test zero-length region."""
        region = GenomicRegion("chr1", 100, 100)
        assert region.length == 0

    def test_str_representation(self):
        """Test string representation (1-based)."""
        region = GenomicRegion("chr1", 999, 2000)
        assert str(region) == "chr1:1000-2000"

    def test_contains_position_inside(self):
        """Test contains with position inside region."""
        region = GenomicRegion("chr1", 100, 200)
        assert region.contains("chr1", 150)
        assert region.contains("chr1", 100)  # Start is inclusive
        assert not region.contains("chr1", 200)  # End is exclusive

    def test_contains_position_outside(self):
        """Test contains with position outside region."""
        region = GenomicRegion("chr1", 100, 200)
        assert not region.contains("chr1", 50)
        assert not region.contains("chr1", 250)

    def test_contains_wrong_scaffold(self):
        """Test contains with different scaffold."""
        region = GenomicRegion("chr1", 100, 200)
        assert not region.contains("chr2", 150)

    def test_overlaps_yes(self):
        """Test overlapping regions."""
        r1 = GenomicRegion("chr1", 100, 200)
        r2 = GenomicRegion("chr1", 150, 250)
        assert r1.overlaps(r2)
        assert r2.overlaps(r1)

    def test_overlaps_adjacent(self):
        """Test adjacent regions (should not overlap)."""
        r1 = GenomicRegion("chr1", 100, 200)
        r2 = GenomicRegion("chr1", 200, 300)
        assert not r1.overlaps(r2)
        assert not r2.overlaps(r1)

    def test_overlaps_different_scaffold(self):
        """Test regions on different scaffolds."""
        r1 = GenomicRegion("chr1", 100, 200)
        r2 = GenomicRegion("chr2", 100, 200)
        assert not r1.overlaps(r2)

    def test_overlaps_contained(self):
        """Test contained region."""
        r1 = GenomicRegion("chr1", 100, 300)
        r2 = GenomicRegion("chr1", 150, 200)
        assert r1.overlaps(r2)
        assert r2.overlaps(r1)

    def test_contains_region_yes(self):
        """Test region containment."""
        r1 = GenomicRegion("chr1", 100, 300)
        r2 = GenomicRegion("chr1", 150, 200)
        assert r1.contains_region(r2)
        assert not r2.contains_region(r1)

    def test_contains_region_same(self):
        """Test region contains itself."""
        r1 = GenomicRegion("chr1", 100, 200)
        assert r1.contains_region(r1)

    def test_contains_region_partial(self):
        """Test partial overlap is not containment."""
        r1 = GenomicRegion("chr1", 100, 200)
        r2 = GenomicRegion("chr1", 150, 250)
        assert not r1.contains_region(r2)
        assert not r2.contains_region(r1)


# =============================================================================
# Test parse_region
# =============================================================================


class TestParseRegion:
    """Tests for parse_region function."""

    def test_standard_format(self):
        """Parse chr1:1000-2000 format."""
        region = parse_region("chr1:1000-2000")
        assert region.seqid == "chr1"
        assert region.start == 999  # 0-based
        assert region.end == 2000  # half-open

    def test_dotted_format(self):
        """Parse chr1:1000..2000 format."""
        region = parse_region("chr1:1000..2000")
        assert region.seqid == "chr1"
        assert region.start == 999
        assert region.end == 2000

    def test_underscore_seqid(self):
        """Handle scaffold_123:100-200."""
        region = parse_region("scaffold_123:100-200")
        assert region.seqid == "scaffold_123"
        assert region.start == 99
        assert region.end == 200

    def test_complex_seqid(self):
        """Handle complex scaffold names."""
        region = parse_region("chr1_random:100-200")
        assert region.seqid == "chr1_random"

    def test_single_base(self):
        """Parse single-base region."""
        region = parse_region("chr1:100-100")
        assert region.start == 99
        assert region.end == 100
        assert region.length == 1

    def test_large_coordinates(self):
        """Parse large genomic coordinates."""
        region = parse_region("chr1:1000000-2000000")
        assert region.start == 999999
        assert region.end == 2000000
        assert region.length == 1000001

    def test_invalid_format_no_colon(self):
        """Invalid format: missing colon."""
        with pytest.raises(ValueError, match="Invalid region format"):
            parse_region("chr1-1000-2000")

    def test_invalid_format_no_dash(self):
        """Invalid format: missing dash."""
        with pytest.raises(ValueError, match="Invalid region format"):
            parse_region("chr1:1000")

    def test_invalid_format_text_only(self):
        """Invalid format: text only."""
        with pytest.raises(ValueError, match="Invalid region format"):
            parse_region("invalid")

    def test_invalid_start_zero(self):
        """Invalid: 0 start position in 1-based input."""
        with pytest.raises(ValueError, match="Start position must be >= 1"):
            parse_region("chr1:0-100")

    def test_invalid_end_less_than_start(self):
        """Invalid: end < start."""
        with pytest.raises(ValueError, match="End must be >= start"):
            parse_region("chr1:200-100")


# =============================================================================
# Test validate_region
# =============================================================================


class TestValidateRegion:
    """Tests for validate_region function."""

    def test_valid_region(self, synthetic_fasta):
        """Valid region within scaffold bounds."""
        from helixforge.io.fasta import GenomeAccessor

        genome = GenomeAccessor(synthetic_fasta)
        region = GenomicRegion("chr1", 100, 500)

        # Should not raise
        validate_region(region, genome)
        genome.close()

    def test_invalid_scaffold(self, synthetic_fasta):
        """Invalid: scaffold doesn't exist."""
        from helixforge.io.fasta import GenomeAccessor

        genome = GenomeAccessor(synthetic_fasta)
        region = GenomicRegion("chrX", 100, 500)

        with pytest.raises(ValueError, match="not found in genome"):
            validate_region(region, genome)

        genome.close()

    def test_region_exceeds_length(self, synthetic_fasta):
        """Invalid: region exceeds scaffold length."""
        from helixforge.io.fasta import GenomeAccessor

        genome = GenomeAccessor(synthetic_fasta)
        region = GenomicRegion("chr1", 100, 2000)  # chr1 is only 1000bp

        with pytest.raises(ValueError, match="exceeds scaffold length"):
            validate_region(region, genome)

        genome.close()


# =============================================================================
# Test region_to_str
# =============================================================================


class TestRegionToStr:
    """Tests for region_to_str function."""

    def test_one_based(self):
        """Convert to 1-based string."""
        region = GenomicRegion("chr1", 999, 2000)
        assert region_to_str(region, one_based=True) == "chr1:1000-2000"

    def test_zero_based(self):
        """Convert to 0-based string."""
        region = GenomicRegion("chr1", 999, 2000)
        assert region_to_str(region, one_based=False) == "chr1:999-2000"


# =============================================================================
# Test region_from_scaffold
# =============================================================================


class TestRegionFromScaffold:
    """Tests for region_from_scaffold function."""

    def test_create_scaffold_region(self):
        """Create region covering entire scaffold."""
        region = region_from_scaffold("chr1", 100000)
        assert region.seqid == "chr1"
        assert region.start == 0
        assert region.end == 100000
        assert region.length == 100000


# =============================================================================
# Test format_region_for_cli
# =============================================================================


class TestFormatRegionForCli:
    """Tests for format_region_for_cli function."""

    def test_from_zero_based(self):
        """Format 0-based coords to 1-based CLI string."""
        result = format_region_for_cli("chr1", 999, 2000, from_zero_based=True)
        assert result == "chr1:1000-2000"

    def test_from_one_based(self):
        """Format already 1-based coords."""
        result = format_region_for_cli("chr1", 1000, 2000, from_zero_based=False)
        assert result == "chr1:1000-2000"


# =============================================================================
# Integration Tests
# =============================================================================


class TestRegionRoundTrip:
    """Test round-trip parsing and formatting."""

    def test_parse_and_format(self):
        """Parse region and format back to string."""
        original = "chr1:1000-2000"
        region = parse_region(original)
        result = str(region)
        assert result == original

    def test_format_and_parse(self):
        """Format region and parse back."""
        region = GenomicRegion("chr1", 999, 2000)
        formatted = region_to_str(region, one_based=True)
        parsed = parse_region(formatted)
        assert parsed == region
