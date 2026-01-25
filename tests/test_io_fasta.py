"""Unit tests for helixforge.io.fasta module.

Tests cover:
- GenomeAccessor initialization and indexing
- Sequence extraction with strand awareness
- Coordinate validation
- Utility functions (reverse_complement, GC content)
"""

from pathlib import Path

import numpy as np
import pytest

from helixforge.io.fasta import (
    COMPLEMENT_TABLE,
    GenomeAccessor,
    calculate_gc_content,
    extract_sequences_for_regions,
    load_genome,
    reverse_complement,
)


# =============================================================================
# Utility Function Tests
# =============================================================================


class TestReverseComplement:
    """Tests for reverse_complement function."""

    def test_simple_sequence(self) -> None:
        """Test reverse complement of a simple sequence."""
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GCGC") == "GCGC"

    def test_longer_sequence(self) -> None:
        """Test reverse complement of a longer sequence."""
        assert reverse_complement("ATCGATCG") == "CGATCGAT"

    def test_case_preservation(self) -> None:
        """Test that case is preserved."""
        assert reverse_complement("AcGt") == "aCgT"
        assert reverse_complement("acgt") == "acgt"

    def test_n_bases(self) -> None:
        """Test handling of N bases."""
        assert reverse_complement("ACNGT") == "ACNGT"
        assert reverse_complement("NNN") == "NNN"

    def test_iupac_ambiguity(self) -> None:
        """Test IUPAC ambiguity codes."""
        # R (A/G) complements to Y (C/T)
        assert reverse_complement("R") == "Y"
        # S (G/C) complements to S
        assert reverse_complement("S") == "S"

    def test_empty_sequence(self) -> None:
        """Test empty sequence."""
        assert reverse_complement("") == ""


class TestCalculateGCContent:
    """Tests for calculate_gc_content function."""

    def test_all_gc(self) -> None:
        """Test sequence that is all GC."""
        assert calculate_gc_content("GCGCGC") == 1.0

    def test_no_gc(self) -> None:
        """Test sequence with no GC."""
        assert calculate_gc_content("ATATAT") == 0.0

    def test_half_gc(self) -> None:
        """Test sequence that is half GC."""
        assert calculate_gc_content("ATGC") == 0.5

    def test_case_insensitive(self) -> None:
        """Test case insensitivity."""
        assert calculate_gc_content("gcgc") == 1.0
        assert calculate_gc_content("GcGc") == 1.0

    def test_n_excluded(self) -> None:
        """Test that N bases are excluded from calculation."""
        # 2 G, 2 C, 2 A, 2 T, 2 N = 4/8 = 0.5
        assert calculate_gc_content("GCATGCATNN") == 0.5

    def test_all_n(self) -> None:
        """Test sequence that is all N."""
        assert calculate_gc_content("NNNNNN") == 0.0

    def test_empty_sequence(self) -> None:
        """Test empty sequence."""
        assert calculate_gc_content("") == 0.0


# =============================================================================
# GenomeAccessor Tests
# =============================================================================


class TestGenomeAccessorInit:
    """Tests for GenomeAccessor initialization."""

    def test_init_success(self, synthetic_fasta: Path) -> None:
        """Test successful initialization."""
        genome = GenomeAccessor(synthetic_fasta)

        assert genome.path == synthetic_fasta
        assert genome.n_scaffolds == 2
        assert "chr1" in genome.scaffold_lengths
        assert "chr2" in genome.scaffold_lengths

        genome.close()

    def test_init_creates_index(self, synthetic_fasta: Path) -> None:
        """Test that .fai index is created if missing."""
        genome = GenomeAccessor(synthetic_fasta)

        fai_path = genome.get_fai_path()
        assert fai_path.exists()

        genome.close()

    def test_init_missing_file(self, tmp_path: Path) -> None:
        """Test error for missing FASTA file."""
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            GenomeAccessor(tmp_path / "missing.fa")

    def test_context_manager(self, synthetic_fasta: Path) -> None:
        """Test context manager usage."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert genome.n_scaffolds == 2

        # File should be closed after context exit
        assert genome._fasta is None


class TestGenomeAccessorProperties:
    """Tests for GenomeAccessor properties."""

    def test_scaffold_lengths(self, synthetic_fasta: Path) -> None:
        """Test scaffold_lengths property."""
        with GenomeAccessor(synthetic_fasta) as genome:
            lengths = genome.scaffold_lengths

            assert lengths["chr1"] == 1000
            assert lengths["chr2"] == 500

            # Should be a copy
            lengths["chr1"] = 0
            assert genome.scaffold_lengths["chr1"] == 1000

    def test_scaffold_order(self, synthetic_fasta: Path) -> None:
        """Test scaffold_order property."""
        with GenomeAccessor(synthetic_fasta) as genome:
            order = genome.scaffold_order

            assert order == ["chr1", "chr2"]

            # Should be a copy
            order.append("chr3")
            assert len(genome.scaffold_order) == 2

    def test_total_length(self, synthetic_fasta: Path) -> None:
        """Test total_length property."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert genome.total_length == 1500

    def test_n_scaffolds(self, synthetic_fasta: Path) -> None:
        """Test n_scaffolds property."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert genome.n_scaffolds == 2


class TestGenomeAccessorGetSequence:
    """Tests for GenomeAccessor.get_sequence method."""

    def test_simple_extraction(self, synthetic_fasta: Path) -> None:
        """Test simple sequence extraction."""
        with GenomeAccessor(synthetic_fasta) as genome:
            seq = genome.get_sequence("chr1", 0, 10)

            assert len(seq) == 10
            assert all(c in "ACGTacgt" for c in seq)

    def test_full_scaffold(self, synthetic_fasta: Path) -> None:
        """Test extracting entire scaffold."""
        with GenomeAccessor(synthetic_fasta) as genome:
            seq = genome.get_sequence("chr1", 0, 1000)
            assert len(seq) == 1000

    def test_reverse_strand(self, synthetic_fasta: Path) -> None:
        """Test extraction with reverse complement."""
        with GenomeAccessor(synthetic_fasta) as genome:
            forward = genome.get_sequence("chr1", 0, 10, strand="+")
            reverse = genome.get_sequence("chr1", 0, 10, strand="-")

            # Reverse complement should be different
            assert reverse == reverse_complement(forward)

    def test_chr2_extraction(self, synthetic_fasta: Path) -> None:
        """Test extraction from second scaffold."""
        with GenomeAccessor(synthetic_fasta) as genome:
            seq = genome.get_sequence("chr2", 100, 200)
            assert len(seq) == 100

    def test_invalid_seqid(self, synthetic_fasta: Path) -> None:
        """Test error for unknown scaffold."""
        with GenomeAccessor(synthetic_fasta) as genome:
            with pytest.raises(KeyError, match="Unknown scaffold"):
                genome.get_sequence("chrX", 0, 100)

    def test_negative_start(self, synthetic_fasta: Path) -> None:
        """Test error for negative start position."""
        with GenomeAccessor(synthetic_fasta) as genome:
            with pytest.raises(ValueError, match="cannot be negative"):
                genome.get_sequence("chr1", -1, 100)

    def test_end_beyond_scaffold(self, synthetic_fasta: Path) -> None:
        """Test error for end beyond scaffold length."""
        with GenomeAccessor(synthetic_fasta) as genome:
            with pytest.raises(ValueError, match="exceeds scaffold length"):
                genome.get_sequence("chr1", 0, 2000)

    def test_start_ge_end(self, synthetic_fasta: Path) -> None:
        """Test error when start >= end."""
        with GenomeAccessor(synthetic_fasta) as genome:
            with pytest.raises(ValueError, match="must be less than"):
                genome.get_sequence("chr1", 100, 100)

            with pytest.raises(ValueError, match="must be less than"):
                genome.get_sequence("chr1", 200, 100)


class TestGenomeAccessorGetLength:
    """Tests for GenomeAccessor.get_length method."""

    def test_get_length(self, synthetic_fasta: Path) -> None:
        """Test getting scaffold length."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert genome.get_length("chr1") == 1000
            assert genome.get_length("chr2") == 500

    def test_invalid_seqid(self, synthetic_fasta: Path) -> None:
        """Test error for unknown scaffold."""
        with GenomeAccessor(synthetic_fasta) as genome:
            with pytest.raises(KeyError, match="Unknown scaffold"):
                genome.get_length("chrX")


class TestGenomeAccessorValidateRegion:
    """Tests for GenomeAccessor.validate_region method."""

    def test_valid_region(self, synthetic_fasta: Path) -> None:
        """Test validation of valid region."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert genome.validate_region("chr1", 0, 100) is True
            assert genome.validate_region("chr1", 500, 1000) is True
            assert genome.validate_region("chr2", 0, 500) is True

    def test_invalid_seqid(self, synthetic_fasta: Path) -> None:
        """Test validation fails for unknown scaffold."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert genome.validate_region("chrX", 0, 100) is False

    def test_invalid_coords(self, synthetic_fasta: Path) -> None:
        """Test validation fails for invalid coordinates."""
        with GenomeAccessor(synthetic_fasta) as genome:
            # Negative start
            assert genome.validate_region("chr1", -1, 100) is False
            # End beyond scaffold
            assert genome.validate_region("chr1", 0, 2000) is False
            # Start >= end
            assert genome.validate_region("chr1", 100, 100) is False


class TestGenomeAccessorBuildOffsetIndex:
    """Tests for GenomeAccessor.build_offset_index method."""

    def test_build_offset_index(self, synthetic_fasta: Path) -> None:
        """Test building cumulative offset index."""
        with GenomeAccessor(synthetic_fasta) as genome:
            offsets = genome.build_offset_index()

            assert offsets["chr1"] == 0
            assert offsets["chr2"] == 1000


class TestGenomeAccessorGCContent:
    """Tests for GenomeAccessor.get_gc_content method."""

    def test_gc_content_full_scaffold(
        self, synthetic_fasta_with_n: Path
    ) -> None:
        """Test GC content for full scaffold."""
        with GenomeAccessor(synthetic_fasta_with_n) as genome:
            gc = genome.get_gc_content("chr1")
            # 10G + 10C out of 40 valid bases = 0.5
            assert gc == pytest.approx(0.5)

    def test_gc_content_region(self, synthetic_fasta_with_n: Path) -> None:
        """Test GC content for a region."""
        with GenomeAccessor(synthetic_fasta_with_n) as genome:
            # First 20 bases: 10G + 10C = 100% GC
            gc = genome.get_gc_content("chr1", 0, 20)
            assert gc == pytest.approx(1.0)


class TestGenomeAccessorDinucleotide:
    """Tests for GenomeAccessor.get_dinucleotide method."""

    def test_get_dinucleotide(self, synthetic_fasta: Path) -> None:
        """Test dinucleotide extraction."""
        with GenomeAccessor(synthetic_fasta) as genome:
            di = genome.get_dinucleotide("chr1", 100)
            assert len(di) == 2

    def test_get_dinucleotide_reverse(self, synthetic_fasta: Path) -> None:
        """Test dinucleotide extraction on reverse strand."""
        with GenomeAccessor(synthetic_fasta) as genome:
            forward = genome.get_dinucleotide("chr1", 100, strand="+")
            reverse = genome.get_dinucleotide("chr1", 100, strand="-")

            assert reverse == reverse_complement(forward)

    def test_out_of_bounds(self, synthetic_fasta: Path) -> None:
        """Test error for out of bounds dinucleotide."""
        with GenomeAccessor(synthetic_fasta) as genome:
            with pytest.raises(ValueError, match="out of bounds"):
                genome.get_dinucleotide("chr1", 999)


class TestGenomeAccessorIterScaffolds:
    """Tests for GenomeAccessor.iter_scaffolds method."""

    def test_iter_scaffolds(self, synthetic_fasta: Path) -> None:
        """Test iterating over scaffolds."""
        with GenomeAccessor(synthetic_fasta) as genome:
            scaffolds = list(genome.iter_scaffolds())

            assert len(scaffolds) == 2
            assert scaffolds[0][0] == "chr1"
            assert len(scaffolds[0][1]) == 1000
            assert scaffolds[1][0] == "chr2"
            assert len(scaffolds[1][1]) == 500


class TestGenomeAccessorDunderMethods:
    """Tests for GenomeAccessor dunder methods."""

    def test_contains(self, synthetic_fasta: Path) -> None:
        """Test __contains__ method."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert "chr1" in genome
            assert "chr2" in genome
            assert "chrX" not in genome

    def test_len(self, synthetic_fasta: Path) -> None:
        """Test __len__ method."""
        with GenomeAccessor(synthetic_fasta) as genome:
            assert len(genome) == 2


# =============================================================================
# Convenience Function Tests
# =============================================================================


class TestLoadGenome:
    """Tests for load_genome convenience function."""

    def test_load_genome(self, synthetic_fasta: Path) -> None:
        """Test load_genome function."""
        genome = load_genome(synthetic_fasta)

        assert isinstance(genome, GenomeAccessor)
        assert genome.n_scaffolds == 2

        genome.close()


class TestExtractSequencesForRegions:
    """Tests for extract_sequences_for_regions function."""

    def test_extract_multiple_regions(self, synthetic_fasta: Path) -> None:
        """Test extracting sequences for multiple regions."""
        with GenomeAccessor(synthetic_fasta) as genome:
            regions = [
                ("chr1", 0, 10, "+"),
                ("chr1", 100, 110, "-"),
                ("chr2", 50, 60, "+"),
            ]
            seqs = extract_sequences_for_regions(genome, regions)

            assert len(seqs) == 3
            assert all(len(s) == 10 for s in seqs)

    def test_empty_regions(self, synthetic_fasta: Path) -> None:
        """Test with empty regions list."""
        with GenomeAccessor(synthetic_fasta) as genome:
            seqs = extract_sequences_for_regions(genome, [])
            assert seqs == []
