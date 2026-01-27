"""Tests for evidence scoring module."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from helixforge.core.evidence import (
    EvidenceLevel,
    EvidenceScore,
    EvidenceScorer,
    EvidenceScorerConfig,
    ExonEvidence,
    JunctionEvidence,
    summarize_evidence_scores,
)


# =============================================================================
# Test EvidenceLevel Enum
# =============================================================================


class TestEvidenceLevel:
    """Test EvidenceLevel enum."""

    def test_level_values(self):
        """Test evidence level string values."""
        assert EvidenceLevel.FULL.value == "full"
        assert EvidenceLevel.PARTIAL.value == "partial"
        assert EvidenceLevel.MINIMAL.value == "minimal"
        assert EvidenceLevel.NONE.value == "none"

    def test_level_comparison(self):
        """Test evidence level ordering."""
        assert EvidenceLevel.NONE < EvidenceLevel.MINIMAL
        assert EvidenceLevel.MINIMAL < EvidenceLevel.PARTIAL
        assert EvidenceLevel.PARTIAL < EvidenceLevel.FULL

        assert EvidenceLevel.FULL > EvidenceLevel.PARTIAL
        assert EvidenceLevel.PARTIAL >= EvidenceLevel.PARTIAL
        assert EvidenceLevel.NONE <= EvidenceLevel.NONE


# =============================================================================
# Test JunctionEvidence
# =============================================================================


class TestJunctionEvidence:
    """Test JunctionEvidence data structure."""

    def test_junction_creation(self):
        """Test creating a JunctionEvidence instance."""
        je = JunctionEvidence(
            intron_start=1000,
            intron_end=2000,
            predicted=True,
            observed=True,
            read_count=50,
            unique_count=45,
            exact_match=True,
        )
        assert je.intron_start == 1000
        assert je.intron_end == 2000
        assert je.intron_length == 1000
        assert je.is_supported is True
        assert je.exact_match is True

    def test_junction_unsupported(self):
        """Test unsupported junction."""
        je = JunctionEvidence(
            intron_start=1000,
            intron_end=2000,
            predicted=True,
            observed=False,
            read_count=0,
        )
        assert je.is_supported is False
        assert je.observed is False

    def test_junction_low_support(self):
        """Test junction with low read support."""
        je = JunctionEvidence(
            intron_start=1000,
            intron_end=2000,
            predicted=True,
            observed=True,
            read_count=2,  # Below default threshold of 3
        )
        assert je.observed is True
        assert je.is_supported is False  # Below threshold

    def test_junction_with_shift(self):
        """Test junction with coordinate shift."""
        je = JunctionEvidence(
            intron_start=1000,
            intron_end=2000,
            predicted=True,
            observed=True,
            read_count=10,
            exact_match=False,
            shift_donor=3,
            shift_acceptor=-2,
        )
        assert je.exact_match is False
        assert je.shift_donor == 3
        assert je.shift_acceptor == -2


# =============================================================================
# Test ExonEvidence
# =============================================================================


class TestExonEvidence:
    """Test ExonEvidence data structure."""

    def test_exon_creation(self):
        """Test creating an ExonEvidence instance."""
        ee = ExonEvidence(
            exon_start=100,
            exon_end=300,
            mean_coverage=50.5,
            median_coverage=48.0,
            min_coverage=10,
            max_coverage=100,
            fraction_covered=0.95,
            coverage_uniformity=0.3,
        )
        assert ee.exon_start == 100
        assert ee.exon_end == 300
        assert ee.exon_length == 200
        assert ee.is_expressed is True  # median >= 5

    def test_exon_not_expressed(self):
        """Test non-expressed exon."""
        ee = ExonEvidence(
            exon_start=100,
            exon_end=300,
            mean_coverage=2.0,
            median_coverage=1.0,
            min_coverage=0,
            max_coverage=5,
            fraction_covered=0.1,
        )
        assert ee.is_expressed is False  # median < 5

    def test_exon_zero_coverage(self):
        """Test exon with no coverage."""
        ee = ExonEvidence(
            exon_start=100,
            exon_end=300,
        )
        assert ee.mean_coverage == 0.0
        assert ee.is_expressed is False


# =============================================================================
# Test EvidenceScore
# =============================================================================


class TestEvidenceScore:
    """Test EvidenceScore data structure."""

    def test_score_creation(self):
        """Test creating an EvidenceScore instance."""
        score = EvidenceScore(
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
        )
        assert score.has_full_junction_support is True
        assert score.is_single_exon is False

    def test_single_exon_gene(self):
        """Test single-exon gene detection."""
        score = EvidenceScore(
            gene_id="gene1",
            transcript_id="tx1",
            seqid="chr1",
            strand="+",
            n_introns=0,
            junction_support_ratio=1.0,  # No introns = 1.0 by convention
        )
        assert score.is_single_exon is True
        assert score.has_full_junction_support is True

    def test_partial_junction_support(self):
        """Test partial junction support."""
        score = EvidenceScore(
            gene_id="gene1",
            transcript_id="tx1",
            seqid="chr1",
            strand="+",
            n_introns=4,
            n_introns_supported=2,
            junction_support_ratio=0.5,
        )
        assert score.has_full_junction_support is False
        assert score.has_partial_junction_support is True

    def test_no_junction_support(self):
        """Test no junction support."""
        score = EvidenceScore(
            gene_id="gene1",
            transcript_id="tx1",
            seqid="chr1",
            strand="+",
            n_introns=4,
            n_introns_supported=0,
            junction_support_ratio=0.0,
        )
        assert score.has_full_junction_support is False
        assert score.has_partial_junction_support is False


# =============================================================================
# Test EvidenceScorerConfig
# =============================================================================


class TestEvidenceScorerConfig:
    """Test EvidenceScorerConfig."""

    def test_default_config(self):
        """Test default configuration values."""
        config = EvidenceScorerConfig()
        assert config.min_junction_reads == 3
        assert config.min_exon_coverage == 5
        assert config.boundary_tolerance == 10
        assert config.junction_weight == 0.5
        assert config.coverage_weight == 0.3
        assert config.boundary_weight == 0.2

    def test_custom_config(self):
        """Test custom configuration."""
        config = EvidenceScorerConfig(
            min_junction_reads=5,
            min_exon_coverage=10,
            boundary_tolerance=15,
        )
        assert config.min_junction_reads == 5
        assert config.min_exon_coverage == 10
        assert config.boundary_tolerance == 15

    def test_weight_normalization_warning(self, caplog):
        """Test warning for non-normalized weights."""
        import logging

        with caplog.at_level(logging.WARNING):
            config = EvidenceScorerConfig(
                junction_weight=0.5,
                coverage_weight=0.5,
                boundary_weight=0.5,
            )
        # Should warn that weights sum to 1.5, not 1.0
        assert "weights sum to" in caplog.text


# =============================================================================
# Test EvidenceScorer
# =============================================================================


class TestEvidenceScorer:
    """Test EvidenceScorer class."""

    def test_scorer_creation(self):
        """Test creating an EvidenceScorer."""
        scorer = EvidenceScorer()
        assert scorer.config.min_junction_reads == 3

    def test_scorer_with_config(self):
        """Test creating scorer with custom config."""
        config = EvidenceScorerConfig(min_junction_reads=10)
        scorer = EvidenceScorer(config)
        assert scorer.config.min_junction_reads == 10

    def test_build_junction_index(self):
        """Test building junction index."""
        from helixforge.io.bam import SpliceJunction

        junctions = {
            "chr1": [
                SpliceJunction("chr1", 1000, 2000, "+", 50, 45),
                SpliceJunction("chr1", 3000, 4000, "+", 30, 28),
            ],
            "chr2": [
                SpliceJunction("chr2", 500, 1500, "-", 20, 18),
            ],
        }

        scorer = EvidenceScorer()
        scorer.build_junction_index(junctions)

        assert "chr1" in scorer._junction_index
        assert "chr2" in scorer._junction_index
        assert (1000, 2000) in scorer._junction_index["chr1"]

    def test_find_junction_exact(self):
        """Test finding junction with exact match."""
        from helixforge.io.bam import SpliceJunction

        junctions = {
            "chr1": [
                SpliceJunction("chr1", 1000, 2000, "+", 50, 45),
            ],
        }

        scorer = EvidenceScorer()
        scorer.build_junction_index(junctions)

        junction, donor_shift, acceptor_shift = scorer._find_junction(
            "chr1", 1000, 2000
        )
        assert junction is not None
        assert junction.read_count == 50
        assert donor_shift == 0
        assert acceptor_shift == 0

    def test_find_junction_with_tolerance(self):
        """Test finding junction within tolerance."""
        from helixforge.io.bam import SpliceJunction

        junctions = {
            "chr1": [
                SpliceJunction("chr1", 1000, 2000, "+", 50, 45),
            ],
        }

        scorer = EvidenceScorer()
        scorer.build_junction_index(junctions)

        # Query slightly shifted coordinates
        junction, donor_shift, acceptor_shift = scorer._find_junction(
            "chr1", 1003, 1998, tolerance=5
        )
        assert junction is not None
        assert donor_shift == -3  # junction at 1000, query at 1003
        assert acceptor_shift == 2  # junction at 2000, query at 1998

    def test_find_junction_not_found(self):
        """Test junction not found."""
        scorer = EvidenceScorer()
        scorer.build_junction_index({})

        junction, donor_shift, acceptor_shift = scorer._find_junction(
            "chr1", 1000, 2000
        )
        assert junction is None

    def test_score_junction_supported(self):
        """Test scoring a supported junction."""
        from helixforge.io.bam import SpliceJunction

        junctions = {
            "chr1": [
                SpliceJunction("chr1", 1000, 2000, "+", 50, 45),
            ],
        }

        scorer = EvidenceScorer()
        scorer.build_junction_index(junctions)

        evidence = scorer.score_junction(1000, 2000, "chr1")
        assert evidence.observed is True
        assert evidence.read_count == 50
        assert evidence.exact_match is True
        assert evidence.is_supported is True

    def test_score_junction_unsupported(self):
        """Test scoring an unsupported junction."""
        scorer = EvidenceScorer()
        scorer.build_junction_index({})

        evidence = scorer.score_junction(1000, 2000, "chr1")
        assert evidence.observed is False
        assert evidence.read_count == 0
        assert evidence.is_supported is False


class TestEvidenceScorerAED:
    """Test AED calculation."""

    def test_aed_perfect_support(self):
        """Test AED of 0 for perfect support."""
        scorer = EvidenceScorer()
        aed = scorer._calculate_aed(
            junction_ratio=1.0,
            coverage_ratio=1.0,
            boundary_ratio=1.0,
        )
        assert aed == pytest.approx(0.0)

    def test_aed_no_support(self):
        """Test AED of 1 for no support."""
        scorer = EvidenceScorer()
        aed = scorer._calculate_aed(
            junction_ratio=0.0,
            coverage_ratio=0.0,
            boundary_ratio=0.0,
        )
        assert aed == pytest.approx(1.0)

    def test_aed_partial_support(self):
        """Test AED for partial support."""
        scorer = EvidenceScorer()
        # With default weights (0.5, 0.3, 0.2):
        # AED = 0.5*(1-0.5) + 0.3*(1-0.5) + 0.2*(1-0.5)
        #     = 0.5*0.5 + 0.3*0.5 + 0.2*0.5
        #     = 0.25 + 0.15 + 0.1 = 0.5
        aed = scorer._calculate_aed(
            junction_ratio=0.5,
            coverage_ratio=0.5,
            boundary_ratio=0.5,
        )
        assert aed == pytest.approx(0.5)


class TestEvidenceLevelClassification:
    """Test evidence level classification."""

    def test_classify_full_multi_exon(self):
        """Test full support classification for multi-exon gene."""
        scorer = EvidenceScorer()
        level = scorer._classify_evidence_level(
            junction_ratio=1.0,
            coverage_ratio=0.9,
            is_single_exon=False,
        )
        assert level == EvidenceLevel.FULL

    def test_classify_partial_multi_exon(self):
        """Test partial support classification."""
        scorer = EvidenceScorer()
        level = scorer._classify_evidence_level(
            junction_ratio=0.8,
            coverage_ratio=0.6,
            is_single_exon=False,
        )
        assert level == EvidenceLevel.PARTIAL

    def test_classify_minimal_multi_exon(self):
        """Test minimal support classification."""
        scorer = EvidenceScorer()
        level = scorer._classify_evidence_level(
            junction_ratio=0.2,
            coverage_ratio=0.1,
            is_single_exon=False,
        )
        assert level == EvidenceLevel.MINIMAL

    def test_classify_none_multi_exon(self):
        """Test no support classification."""
        scorer = EvidenceScorer()
        level = scorer._classify_evidence_level(
            junction_ratio=0.0,
            coverage_ratio=0.0,
            is_single_exon=False,
        )
        assert level == EvidenceLevel.NONE

    def test_classify_full_single_exon(self):
        """Test full support for single-exon gene."""
        scorer = EvidenceScorer()
        level = scorer._classify_evidence_level(
            junction_ratio=1.0,  # Ignored for single-exon
            coverage_ratio=0.9,
            is_single_exon=True,
        )
        assert level == EvidenceLevel.FULL

    def test_classify_none_single_exon(self):
        """Test no support for single-exon gene."""
        scorer = EvidenceScorer()
        level = scorer._classify_evidence_level(
            junction_ratio=1.0,
            coverage_ratio=0.0,
            is_single_exon=True,
        )
        assert level == EvidenceLevel.NONE


# =============================================================================
# Test Summary Functions
# =============================================================================


class TestSummarizeEvidenceScores:
    """Test summarize_evidence_scores function."""

    def test_empty_scores(self):
        """Test summary with no scores."""
        summary = summarize_evidence_scores([])
        assert summary["n_genes"] == 0
        assert summary["mean_aed"] == 0.0

    def test_single_score(self):
        """Test summary with one score."""
        scores = [
            EvidenceScore(
                gene_id="gene1",
                transcript_id="tx1",
                seqid="chr1",
                strand="+",
                evidence_level=EvidenceLevel.FULL,
                aed=0.1,
            )
        ]
        summary = summarize_evidence_scores(scores)
        assert summary["n_genes"] == 1
        assert summary["mean_aed"] == pytest.approx(0.1)
        assert summary["n_full_support"] == 1

    def test_multiple_scores(self):
        """Test summary with multiple scores."""
        scores = [
            EvidenceScore(
                gene_id="gene1",
                transcript_id="tx1",
                seqid="chr1",
                strand="+",
                evidence_level=EvidenceLevel.FULL,
                aed=0.0,
            ),
            EvidenceScore(
                gene_id="gene2",
                transcript_id="tx2",
                seqid="chr1",
                strand="+",
                evidence_level=EvidenceLevel.PARTIAL,
                aed=0.3,
            ),
            EvidenceScore(
                gene_id="gene3",
                transcript_id="tx3",
                seqid="chr1",
                strand="+",
                evidence_level=EvidenceLevel.NONE,
                aed=1.0,
            ),
        ]
        summary = summarize_evidence_scores(scores)
        assert summary["n_genes"] == 3
        assert summary["mean_aed"] == pytest.approx((0.0 + 0.3 + 1.0) / 3)
        assert summary["n_full_support"] == 1
        assert summary["n_partial_support"] == 1
        assert summary["n_no_support"] == 1
