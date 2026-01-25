"""Tests for splice refinement module."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from helixforge.core.splice import (
    GeneSpliceReport,
    IntronModel,
    PositionWeightMatrix,
    SpliceCorrection,
    SpliceSite,
    SpliceSiteType,
    SpliceReportWriter,
    classify_splice_type,
    get_splice_dinucleotides,
    is_canonical_splice_site,
)


# =============================================================================
# Test SpliceSiteClassification
# =============================================================================


class TestSpliceSiteClassification:
    """Test splice site type detection."""

    def test_canonical_gt_ag(self):
        """GT-AG should be classified as canonical GT_AG."""
        site_type = classify_splice_type("GT", "AG")
        assert site_type == SpliceSiteType.GT_AG

    def test_minor_gc_ag(self):
        """GC-AG should be classified as minor canonical."""
        site_type = classify_splice_type("GC", "AG")
        assert site_type == SpliceSiteType.GC_AG

    def test_u12_at_ac(self):
        """AT-AC should be classified as U12 minor spliceosome."""
        site_type = classify_splice_type("AT", "AC")
        assert site_type == SpliceSiteType.AT_AC

    def test_noncanonical_detection(self):
        """Non-canonical splice sites should be detected."""
        # Various non-canonical pairs
        assert classify_splice_type("AA", "AG") == SpliceSiteType.OTHER
        assert classify_splice_type("GT", "AA") == SpliceSiteType.OTHER
        assert classify_splice_type("CC", "TT") == SpliceSiteType.OTHER

    def test_is_canonical_function(self):
        """Test is_canonical_splice_site helper function."""
        assert is_canonical_splice_site("GT", "AG") is True
        assert is_canonical_splice_site("GC", "AG") is True
        assert is_canonical_splice_site("AT", "AC") is True  # U12 is also "canonical"
        assert is_canonical_splice_site("AA", "AG") is False

    def test_get_splice_dinucleotides(self):
        """Test extraction of dinucleotides from sequence."""
        sequence = "ATCGTAAGTTGCAG"
        #                ^^ donor  ^^ acceptor
        # Intron: positions 4-12 (GTAAGTTG)
        donor, acceptor = get_splice_dinucleotides(sequence, 4, 12)
        assert donor == "GT"
        assert acceptor == "AG"


class TestSpliceSiteModel:
    """Test SpliceSite data structure."""

    def test_splice_site_creation(self):
        """Test creating a SpliceSite instance."""
        site = SpliceSite(
            seqid="chr1",
            position=1000,
            site_type="donor",
            strand="+",
            dinucleotide="GT",
        )
        assert site.seqid == "chr1"
        assert site.position == 1000
        assert site.site_type == "donor"
        assert site.is_canonical is True

    def test_splice_site_type_property(self):
        """Test splice_site_type property."""
        donor_gt = SpliceSite("chr1", 100, "donor", "+", "GT")
        assert donor_gt.splice_site_type == SpliceSiteType.GT_AG

        donor_gc = SpliceSite("chr1", 100, "donor", "+", "GC")
        assert donor_gc.splice_site_type == SpliceSiteType.GC_AG

        acceptor_ag = SpliceSite("chr1", 200, "acceptor", "+", "AG")
        assert acceptor_ag.splice_site_type == SpliceSiteType.GT_AG


class TestIntronModel:
    """Test IntronModel data structure."""

    def test_intron_creation(self):
        """Test creating an IntronModel."""
        intron = IntronModel(
            seqid="chr1",
            start=1000,
            end=2000,
            strand="+",
            parent_transcript="tx1",
            donor_dinucleotide="GT",
            acceptor_dinucleotide="AG",
        )
        assert intron.length == 1000
        assert intron.is_canonical is True
        assert intron.splice_type == SpliceSiteType.GT_AG

    def test_intron_noncanonical(self):
        """Test non-canonical intron detection."""
        intron = IntronModel(
            seqid="chr1",
            start=1000,
            end=2000,
            strand="+",
            parent_transcript="tx1",
            donor_dinucleotide="AA",
            acceptor_dinucleotide="TT",
        )
        assert intron.is_canonical is False
        assert intron.splice_type == SpliceSiteType.OTHER

    def test_intron_with_evidence(self):
        """Test intron with RNA-seq evidence."""
        intron = IntronModel(
            seqid="chr1",
            start=1000,
            end=2000,
            strand="+",
            parent_transcript="tx1",
            donor_dinucleotide="GT",
            acceptor_dinucleotide="AG",
            rnaseq_support=50,
            pwm_score_donor=3.5,
            pwm_score_acceptor=2.8,
        )
        assert intron.rnaseq_support == 50
        assert intron.pwm_score_donor == pytest.approx(3.5)


# =============================================================================
# Test PWM Scoring
# =============================================================================


class TestPWMScoring:
    """Test position weight matrix operations."""

    def test_pwm_from_sequences(self):
        """Build PWM from aligned sequences."""
        # Simple donor sequences (9 positions)
        sequences = [
            "CAGGTAAGT",
            "AAGGTAAGT",
            "GAGGTAAGT",
            "CAGGTGAGT",
        ]
        pwm = PositionWeightMatrix.from_sequences(sequences, "donor")

        assert pwm.site_type == "donor"
        assert pwm.window_size == 9
        assert pwm.matrix.shape == (9, 4)

        # Position 3 (G of GT) should have high G frequency
        g_idx = 2  # G is index 2 in A, C, G, T
        assert pwm.matrix[3, g_idx] > 0.8  # Should be nearly 1.0

        # Position 4 (T of GT) should have high T frequency
        t_idx = 3  # T is index 3
        assert pwm.matrix[4, t_idx] > 0.8

    def test_pwm_scoring_perfect_match(self):
        """Consensus sequence should score highest."""
        # Create PWM from known sequences
        sequences = [
            "CAGGTAAGT",
            "CAGGTAAGT",
            "CAGGTAAGT",
        ]
        pwm = PositionWeightMatrix.from_sequences(sequences, "donor")

        # Score the consensus - should be high
        consensus_score = pwm.score("CAGGTAAGT")
        assert consensus_score > 0  # Positive log-likelihood ratio

    def test_pwm_scoring_mismatch(self):
        """Non-consensus should score lower."""
        sequences = [
            "CAGGTAAGT",
            "CAGGTAAGT",
            "CAGGTAAGT",
        ]
        pwm = PositionWeightMatrix.from_sequences(sequences, "donor")

        consensus_score = pwm.score("CAGGTAAGT")
        mismatch_score = pwm.score("TATTTCCCC")

        assert mismatch_score < consensus_score

    def test_pwm_scoring_wrong_length(self):
        """Wrong length sequence should raise error."""
        pwm = PositionWeightMatrix.from_sequences(
            ["CAGGTAAGT"] * 3, "donor"
        )

        with pytest.raises(ValueError, match="length"):
            pwm.score("CAGGT")  # Too short

        with pytest.raises(ValueError, match="length"):
            pwm.score("CAGGTAAGTTTT")  # Too long

    def test_score_all_positions(self):
        """Sliding window scoring."""
        sequences = ["CAGGTAAGT"] * 3
        pwm = PositionWeightMatrix.from_sequences(sequences, "donor")

        # Score a longer sequence
        long_seq = "NNNNCAGGTAAGTNNN"
        scores = pwm.score_all_positions(long_seq)

        # Should have positions for sliding window
        expected_positions = len(long_seq) - pwm.window_size + 1
        assert len(scores) == expected_positions

        # Best score should be at position 4 (where CAGGTAAGT starts)
        best_pos = np.argmax(scores)
        assert best_pos == 4

    def test_pwm_serialization(self):
        """Test PWM to_dict and from_dict."""
        original = PositionWeightMatrix.from_sequences(
            ["CAGGTAAGT"] * 5, "donor"
        )

        # Serialize
        data = original.to_dict()
        assert "matrix" in data
        assert "site_type" in data
        assert data["site_type"] == "donor"

        # Deserialize
        restored = PositionWeightMatrix.from_dict(data)
        assert restored.site_type == original.site_type
        assert restored.offset == original.offset
        np.testing.assert_array_almost_equal(restored.matrix, original.matrix)

    def test_load_plant_defaults(self):
        """Test loading default plant PWMs."""
        donor_pwm, acceptor_pwm = PositionWeightMatrix.load_plant_defaults()

        assert donor_pwm.site_type == "donor"
        assert acceptor_pwm.site_type == "acceptor"

        # Donor should be 9 positions
        assert donor_pwm.window_size == 9

        # Acceptor should be 16 positions
        assert acceptor_pwm.window_size == 16

        # G at position 3 (offset) should be invariant
        g_idx = 2
        assert donor_pwm.matrix[donor_pwm.offset, g_idx] > 0.9


# =============================================================================
# Test Splice Correction and Reports
# =============================================================================


class TestSpliceCorrection:
    """Test SpliceCorrection data structure."""

    def test_correction_creation(self):
        """Test creating a SpliceCorrection."""
        corr = SpliceCorrection(
            gene_id="gene1",
            transcript_id="tx1",
            intron_index=0,
            site="donor",
            original_position=1000,
            corrected_position=1003,
            shift=3,
            reason="rnaseq_junction",
            confidence=0.95,
            rnaseq_reads=50,
        )
        assert corr.shift == 3
        assert corr.reason == "rnaseq_junction"

    def test_correction_to_dict(self):
        """Test correction serialization."""
        corr = SpliceCorrection(
            gene_id="gene1",
            transcript_id="tx1",
            intron_index=0,
            site="acceptor",
            original_position=2000,
            corrected_position=1998,
            shift=-2,
            reason="pwm_optimization",
            confidence=0.8,
            rnaseq_reads=0,
        )

        data = corr.to_dict()
        assert data["gene_id"] == "gene1"
        assert data["shift"] == -2
        assert data["reason"] == "pwm_optimization"


class TestGeneSpliceReport:
    """Test GeneSpliceReport data structure."""

    def test_report_creation(self):
        """Test creating a GeneSpliceReport."""
        report = GeneSpliceReport(
            gene_id="gene1",
            n_introns=5,
            n_supported=4,
            n_corrected=1,
            n_canonical=4,
            n_noncanonical=1,
        )

        assert report.support_ratio == 0.8
        assert report.fully_supported is False

    def test_report_fully_supported(self):
        """Test fully supported gene detection."""
        report = GeneSpliceReport(
            gene_id="gene1",
            n_introns=3,
            n_supported=3,
            n_corrected=0,
            n_canonical=3,
            n_noncanonical=0,
        )

        assert report.fully_supported is True
        assert report.support_ratio == 1.0

    def test_report_no_introns(self):
        """Test report for single-exon gene."""
        report = GeneSpliceReport(
            gene_id="gene1",
            n_introns=0,
            n_supported=0,
            n_corrected=0,
            n_canonical=0,
            n_noncanonical=0,
        )

        assert report.support_ratio == 1.0  # Avoid division by zero

    def test_report_to_dict(self):
        """Test report serialization."""
        report = GeneSpliceReport(
            gene_id="gene1",
            n_introns=5,
            n_supported=4,
            n_corrected=1,
            n_canonical=4,
            n_noncanonical=1,
            unsupported_introns=[2],
            flags=["noncanonical_introns:1"],
        )

        data = report.to_dict()
        assert data["gene_id"] == "gene1"
        assert data["support_ratio"] == 0.8
        assert "noncanonical" in data["flags"]


# =============================================================================
# Test SpliceReportWriter
# =============================================================================


class TestSpliceReportWriter:
    """Test output formats."""

    @pytest.fixture
    def sample_reports(self):
        """Create sample reports for testing."""
        corrections = [
            SpliceCorrection(
                gene_id="gene1",
                transcript_id="tx1",
                intron_index=0,
                site="donor",
                original_position=1000,
                corrected_position=1003,
                shift=3,
                reason="rnaseq_junction",
                confidence=0.95,
                rnaseq_reads=50,
            ),
        ]

        reports = [
            GeneSpliceReport(
                gene_id="gene1",
                n_introns=5,
                n_supported=4,
                n_corrected=1,
                n_canonical=4,
                n_noncanonical=1,
                corrections=corrections,
                unsupported_introns=[2],
                flags=["noncanonical_introns:1"],
            ),
            GeneSpliceReport(
                gene_id="gene2",
                n_introns=3,
                n_supported=3,
                n_corrected=0,
                n_canonical=3,
                n_noncanonical=0,
            ),
        ]
        return reports

    def test_tsv_output(self, sample_reports, tmp_path):
        """Test TSV output generation."""
        output_file = tmp_path / "splice_report.tsv"

        SpliceReportWriter.write_tsv(sample_reports, output_file)

        assert output_file.exists()

        # Read and validate
        with open(output_file) as f:
            lines = f.readlines()

        assert len(lines) == 3  # Header + 2 reports
        assert "gene_id" in lines[0]
        assert "gene1" in lines[1]
        assert "gene2" in lines[2]

    def test_corrections_detail(self, sample_reports, tmp_path):
        """Test detailed corrections output."""
        output_file = tmp_path / "corrections.tsv"

        SpliceReportWriter.write_corrections_detail(sample_reports, output_file)

        assert output_file.exists()

        with open(output_file) as f:
            lines = f.readlines()

        assert len(lines) == 2  # Header + 1 correction
        assert "gene1" in lines[1]
        assert "donor" in lines[1]

    def test_summary_statistics(self, sample_reports):
        """Test summary statistics computation."""
        stats = SpliceReportWriter.summary_statistics(sample_reports)

        assert stats["total_genes"] == 2
        assert stats["total_introns"] == 8
        assert stats["supported_introns"] == 7
        assert stats["corrections_made"] == 1
        assert "canonical_rate" in stats


# =============================================================================
# Test Boundary Adjuster
# =============================================================================


class TestBoundaryAdjuster:
    """Test start/stop codon adjustment."""

    def test_find_start_codon(self):
        """Test start codon finder."""
        # This would require a mock GenomeAccessor
        # Placeholder for when integration tests are available
        pass

    def test_find_stop_codon(self):
        """Test stop codon finder."""
        pass

    def test_phase_consistency_check(self):
        """Test phase consistency validation."""
        pass


# =============================================================================
# Additional Edge Case Tests
# =============================================================================


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_sequences_pwm(self):
        """Building PWM from empty sequences should raise error."""
        with pytest.raises(ValueError, match="empty"):
            PositionWeightMatrix.from_sequences([], "donor")

    def test_unequal_length_sequences(self):
        """PWM from unequal length sequences should raise error."""
        with pytest.raises(ValueError, match="length"):
            PositionWeightMatrix.from_sequences(
                ["CAGGTAAGT", "CAGGT"],  # Different lengths
                "donor",
            )

    def test_lowercase_handling(self):
        """PWM should handle lowercase sequences."""
        pwm = PositionWeightMatrix.from_sequences(
            ["caggtaagt", "CAGGTAAGT", "CaGgTaAgT"],
            "donor",
        )
        # Should work without error
        score = pwm.score("caggtaagt")
        assert isinstance(score, float)

    def test_n_base_handling(self):
        """PWM scoring should handle N bases gracefully."""
        pwm = PositionWeightMatrix.from_sequences(
            ["CAGGTAAGT"] * 3, "donor"
        )
        # Score with N - should not raise
        score = pwm.score("NNNGTAAGT")
        assert isinstance(score, float)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
