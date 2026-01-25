"""Unit tests for helixforge.core.confidence module.

Tests cover:
- GeneConfidence and RegionConfidence data structures
- ConfidenceCalculator metrics calculation
- ConfidenceWriter output formats
- Classification and flagging logic
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from helixforge.core.confidence import (
    CLASS_CDS,
    CLASS_INTERGENIC,
    CLASS_INTRON,
    CLASS_UTR,
    ConfidenceCalculator,
    ConfidenceClass,
    ConfidenceMetric,
    ConfidenceWriter,
    DEFAULT_WEIGHTS,
    GeneConfidence,
    RegionConfidence,
)


# =============================================================================
# GeneConfidence Tests
# =============================================================================


class TestGeneConfidence:
    """Tests for GeneConfidence data structure."""

    def test_creation(self) -> None:
        """Test creating a GeneConfidence object."""
        conf = GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            mean_prob=0.85,
            min_prob=0.60,
            median_prob=0.87,
            entropy=0.5,
            boundary_sharpness=0.8,
            coding_consistency=0.9,
        )

        assert conf.gene_id == "gene1"
        assert conf.seqid == "chr1"
        assert conf.start == 100
        assert conf.end == 500
        assert conf.strand == "+"
        assert conf.mean_prob == 0.85

    def test_overall_score_calculation(self) -> None:
        """Test weighted overall score calculation."""
        conf = GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            mean_prob=0.8,
            min_prob=0.6,
            median_prob=0.85,
            entropy=0.5,
            boundary_sharpness=0.7,
            coding_consistency=0.9,
            worst_exon_score=0.75,
        )

        # Manual calculation based on DEFAULT_WEIGHTS
        expected = (
            DEFAULT_WEIGHTS["mean_prob"] * 0.8
            + DEFAULT_WEIGHTS["min_prob"] * 0.6
            + DEFAULT_WEIGHTS["boundary_sharpness"] * 0.7
            + DEFAULT_WEIGHTS["coding_consistency"] * 0.9
            + DEFAULT_WEIGHTS["worst_exon_score"] * 0.75
        )

        assert abs(conf.overall_score - expected) < 1e-6

    def test_n_low_conf_regions(self) -> None:
        """Test low confidence regions count property."""
        conf = GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            mean_prob=0.7,
            min_prob=0.4,
            median_prob=0.75,
            entropy=0.8,
            boundary_sharpness=0.6,
            coding_consistency=0.7,
            low_confidence_regions=[(150, 200, 0.5), (300, 350, 0.55)],
        )

        assert conf.n_low_conf_regions == 2

    def test_to_dict(self) -> None:
        """Test conversion to dictionary."""
        conf = GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            mean_prob=0.85,
            min_prob=0.60,
            median_prob=0.87,
            entropy=0.5,
            boundary_sharpness=0.8,
            coding_consistency=0.9,
            exon_scores=[0.85, 0.90],
            confidence_class="high",
            flags=["weak_exon"],
        )

        d = conf.to_dict()

        assert d["gene_id"] == "gene1"
        assert d["seqid"] == "chr1"
        assert d["mean_prob"] == 0.85
        assert d["confidence_class"] == "high"
        assert d["flags"] == "weak_exon"
        assert d["n_exons"] == 2

    def test_default_values(self) -> None:
        """Test default values for optional fields."""
        conf = GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            mean_prob=0.85,
            min_prob=0.60,
            median_prob=0.87,
            entropy=0.5,
            boundary_sharpness=0.8,
            coding_consistency=0.9,
        )

        assert conf.exon_scores == []
        assert conf.intron_scores == []
        assert conf.low_confidence_regions == []
        assert conf.confidence_class == "medium"
        assert conf.flags == []


class TestRegionConfidence:
    """Tests for RegionConfidence data structure."""

    def test_creation(self) -> None:
        """Test creating a RegionConfidence object."""
        entropy = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        max_prob = np.array([0.9, 0.85, 0.8, 0.75, 0.7])
        smoothed = np.array([0.88, 0.85, 0.82, 0.78, 0.74])

        region = RegionConfidence(
            seqid="chr1",
            start=100,
            end=105,
            per_base_entropy=entropy,
            per_base_max_prob=max_prob,
            smoothed_confidence=smoothed,
        )

        assert region.seqid == "chr1"
        assert region.length == 5
        assert len(region.per_base_entropy) == 5

    def test_mean_entropy(self) -> None:
        """Test mean entropy property."""
        entropy = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        max_prob = np.array([0.9, 0.85, 0.8, 0.75, 0.7])
        smoothed = np.array([0.88, 0.85, 0.82, 0.78, 0.74])

        region = RegionConfidence(
            seqid="chr1",
            start=100,
            end=105,
            per_base_entropy=entropy,
            per_base_max_prob=max_prob,
            smoothed_confidence=smoothed,
        )

        assert abs(region.mean_entropy - 0.7) < 1e-6

    def test_mean_confidence(self) -> None:
        """Test mean confidence property."""
        entropy = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        max_prob = np.array([0.9, 0.85, 0.8, 0.75, 0.7])
        smoothed = np.array([0.88, 0.85, 0.82, 0.78, 0.74])

        region = RegionConfidence(
            seqid="chr1",
            start=100,
            end=105,
            per_base_entropy=entropy,
            per_base_max_prob=max_prob,
            smoothed_confidence=smoothed,
        )

        assert abs(region.mean_confidence - 0.8) < 1e-6


# =============================================================================
# ConfidenceCalculator Tests
# =============================================================================


class TestConfidenceCalculator:
    """Tests for ConfidenceCalculator class."""

    @pytest.fixture
    def mock_reader(self):
        """Create a mock HDF5 reader."""
        reader = MagicMock()
        reader.scaffold_names = ["chr1", "chr2"]
        reader.scaffold_lengths = {"chr1": 1000, "chr2": 500}

        # Create mock predictions (high confidence)
        def mock_predictions(seqid, start, end):
            length = end - start
            # Create softmax-like probabilities
            np.random.seed(42)
            raw = np.random.random((length, 4)).astype(np.float32)
            # Make predictions confident (one class dominates)
            raw[:, 2] += 2.0  # Boost CDS class
            preds = raw / raw.sum(axis=1, keepdims=True)
            return preds

        reader.get_predictions_for_region.side_effect = mock_predictions
        return reader

    @pytest.fixture
    def mock_gene(self):
        """Create a mock gene model."""
        from helixforge.io.gff import GeneModel, TranscriptModel

        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 200), (300, 500)],
            cds=[(100, 200, 0), (300, 500, 0)],
        )

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[transcript],
            source="helixer",
        )

        return gene

    def test_calculator_init(self, mock_reader) -> None:
        """Test calculator initialization."""
        calc = ConfidenceCalculator(mock_reader)

        assert calc.hdf5_reader == mock_reader
        assert calc.low_conf_threshold == 0.7
        assert calc.window_size == 20

    def test_calculator_custom_threshold(self, mock_reader) -> None:
        """Test calculator with custom threshold."""
        calc = ConfidenceCalculator(
            mock_reader, low_conf_threshold=0.6, window_size=30
        )

        assert calc.low_conf_threshold == 0.6
        assert calc.window_size == 30

    def test_score_gene(self, mock_reader, mock_gene) -> None:
        """Test scoring a single gene."""
        calc = ConfidenceCalculator(mock_reader)
        score = calc.score_gene(mock_gene)

        assert score.gene_id == "gene1"
        assert score.seqid == "chr1"
        assert score.start == 100
        assert score.end == 500
        assert 0.0 <= score.mean_prob <= 1.0
        assert 0.0 <= score.min_prob <= 1.0
        assert score.entropy >= 0.0
        assert score.confidence_class in ["high", "medium", "low"]

    def test_score_gene_exon_scores(self, mock_reader, mock_gene) -> None:
        """Test that exon scores are calculated."""
        calc = ConfidenceCalculator(mock_reader)
        score = calc.score_gene(mock_gene)

        # Gene has 2 exons
        assert len(score.exon_scores) == 2
        assert all(0.0 <= s <= 1.0 for s in score.exon_scores)

    def test_get_region_confidence(self, mock_reader) -> None:
        """Test getting region confidence."""
        calc = ConfidenceCalculator(mock_reader)
        region = calc.get_region_confidence("chr1", 100, 200)

        assert region.seqid == "chr1"
        assert region.start == 100
        assert region.end == 200
        assert region.length == 100
        assert len(region.per_base_entropy) == 100
        assert len(region.per_base_max_prob) == 100

    def test_calculate_entropy(self, mock_reader) -> None:
        """Test entropy calculation."""
        calc = ConfidenceCalculator(mock_reader)

        # Uniform distribution has max entropy
        uniform = np.array([[0.25, 0.25, 0.25, 0.25]])
        entropy = calc._calculate_entropy(uniform)
        assert entropy[0] == pytest.approx(2.0, abs=0.01)  # log2(4) = 2

        # Certain prediction has 0 entropy
        certain = np.array([[1.0, 0.0, 0.0, 0.0]])
        entropy = calc._calculate_entropy(certain)
        assert entropy[0] == pytest.approx(0.0, abs=0.01)

    def test_identify_low_confidence_regions(self, mock_reader) -> None:
        """Test low confidence region identification."""
        calc = ConfidenceCalculator(mock_reader)

        # Create array with low-confidence region in middle
        max_probs = np.ones(100) * 0.9
        max_probs[30:50] = 0.5  # Low confidence region

        regions = calc._identify_low_confidence_regions(
            max_probs, gene_start=100, threshold=0.7, min_length=10
        )

        assert len(regions) == 1
        assert regions[0][0] == 130  # genomic start
        assert regions[0][1] == 150  # genomic end
        assert 0.4 <= regions[0][2] <= 0.6  # mean score

    def test_classify_confidence_high(self, mock_reader) -> None:
        """Test high confidence classification."""
        calc = ConfidenceCalculator(mock_reader)

        conf = GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            mean_prob=0.95,
            min_prob=0.85,
            median_prob=0.93,
            entropy=0.3,
            boundary_sharpness=0.9,
            coding_consistency=0.95,
            exon_scores=[0.9, 0.92],
            worst_exon_score=0.9,
        )

        conf_class, flags = calc._classify_confidence(conf)

        assert conf_class == "high"
        assert len(flags) == 0

    def test_classify_confidence_low(self, mock_reader) -> None:
        """Test low confidence classification."""
        calc = ConfidenceCalculator(mock_reader)

        conf = GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            mean_prob=0.5,
            min_prob=0.2,
            median_prob=0.55,
            entropy=1.8,
            boundary_sharpness=0.3,
            coding_consistency=0.4,
            exon_scores=[0.4, 0.5],
            worst_exon_score=0.4,
            low_confidence_regions=[(150, 200, 0.3), (300, 350, 0.35), (400, 450, 0.4), (460, 490, 0.45)],
        )

        conf_class, flags = calc._classify_confidence(conf)

        assert conf_class == "low"
        assert "weak_exon" in flags
        assert "high_entropy" in flags
        assert "uncertain_boundary" in flags

    def test_score_genes_parallel(self, mock_reader, mock_gene) -> None:
        """Test parallel gene scoring."""
        calc = ConfidenceCalculator(mock_reader)
        genes = [mock_gene, mock_gene]

        scores = list(calc.score_genes_parallel(genes, n_workers=1))

        assert len(scores) == 2
        assert all(isinstance(s, GeneConfidence) for s in scores)


# =============================================================================
# ConfidenceWriter Tests
# =============================================================================


class TestConfidenceWriter:
    """Tests for ConfidenceWriter class."""

    @pytest.fixture
    def sample_scores(self) -> list[GeneConfidence]:
        """Create sample confidence scores."""
        return [
            GeneConfidence(
                gene_id="gene1",
                seqid="chr1",
                start=100,
                end=500,
                strand="+",
                mean_prob=0.9,
                min_prob=0.8,
                median_prob=0.88,
                entropy=0.4,
                boundary_sharpness=0.85,
                coding_consistency=0.9,
                exon_scores=[0.88, 0.92],
                confidence_class="high",
                flags=[],
            ),
            GeneConfidence(
                gene_id="gene2",
                seqid="chr1",
                start=600,
                end=900,
                strand="-",
                mean_prob=0.6,
                min_prob=0.4,
                median_prob=0.62,
                entropy=1.2,
                boundary_sharpness=0.5,
                coding_consistency=0.55,
                exon_scores=[0.55, 0.65],
                confidence_class="low",
                flags=["weak_exon", "high_entropy"],
                low_confidence_regions=[(650, 700, 0.45)],
            ),
        ]

    def test_to_tsv(self, sample_scores, tmp_path: Path) -> None:
        """Test TSV output."""
        output_path = tmp_path / "scores.tsv"

        ConfidenceWriter.to_tsv(sample_scores, output_path)

        assert output_path.exists()
        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Header + 2 data rows
        assert len(lines) == 3

        # Check header
        header = lines[0].split("\t")
        assert "gene_id" in header
        assert "mean_prob" in header
        assert "confidence_class" in header

        # Check data
        row1 = lines[1].split("\t")
        assert row1[0] == "gene1"

    def test_to_bed(self, sample_scores, tmp_path: Path) -> None:
        """Test BED output."""
        output_path = tmp_path / "scores.bed"

        ConfidenceWriter.to_bed(sample_scores, output_path)

        assert output_path.exists()
        content = output_path.read_text()
        lines = [l for l in content.strip().split("\n") if not l.startswith("track")]

        assert len(lines) == 2

        # Check first line
        parts = lines[0].split("\t")
        assert parts[0] == "chr1"  # chrom
        assert parts[1] == "100"  # start
        assert parts[2] == "500"  # end
        assert parts[3] == "gene1"  # name

    def test_low_confidence_regions_bed(self, sample_scores, tmp_path: Path) -> None:
        """Test low-confidence regions BED output."""
        output_path = tmp_path / "low_conf.bed"

        ConfidenceWriter.low_confidence_regions_bed(sample_scores, output_path)

        assert output_path.exists()
        content = output_path.read_text()
        lines = [l for l in content.strip().split("\n") if not l.startswith("track")]

        # Only gene2 has low-conf regions
        assert len(lines) == 1
        parts = lines[0].split("\t")
        assert parts[0] == "chr1"
        assert "gene2_lcr" in parts[3]

    def test_to_dataframe(self, sample_scores) -> None:
        """Test DataFrame conversion."""
        df = ConfidenceWriter.to_dataframe(sample_scores)

        assert len(df) == 2
        assert "gene_id" in df.columns
        assert "mean_prob" in df.columns
        assert df.iloc[0]["gene_id"] == "gene1"
        assert df.iloc[1]["confidence_class"] == "low"


# =============================================================================
# Constants and Enums Tests
# =============================================================================


class TestConstants:
    """Tests for module constants."""

    def test_class_indices(self) -> None:
        """Test Helixer class index constants."""
        assert CLASS_INTERGENIC == 0
        assert CLASS_UTR == 1
        assert CLASS_CDS == 2
        assert CLASS_INTRON == 3

    def test_default_weights_sum(self) -> None:
        """Test that default weights sum to 1."""
        total = sum(DEFAULT_WEIGHTS.values())
        assert abs(total - 1.0) < 1e-6

    def test_confidence_metric_enum(self) -> None:
        """Test ConfidenceMetric enum values."""
        assert ConfidenceMetric.MEAN_PROB.value == "mean_prob"
        assert ConfidenceMetric.MIN_PROB.value == "min_prob"
        assert ConfidenceMetric.ENTROPY.value == "entropy"

    def test_confidence_class_enum(self) -> None:
        """Test ConfidenceClass enum values."""
        assert ConfidenceClass.HIGH.value == "high"
        assert ConfidenceClass.MEDIUM.value == "medium"
        assert ConfidenceClass.LOW.value == "low"


# =============================================================================
# Edge Cases Tests
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    @pytest.fixture
    def mock_reader(self):
        """Create a mock HDF5 reader."""
        reader = MagicMock()
        reader.scaffold_names = ["chr1"]
        reader.scaffold_lengths = {"chr1": 1000}

        def mock_predictions(seqid, start, end):
            length = end - start
            np.random.seed(42)
            raw = np.random.random((length, 4)).astype(np.float32)
            raw[:, 2] += 2.0
            preds = raw / raw.sum(axis=1, keepdims=True)
            return preds

        reader.get_predictions_for_region.side_effect = mock_predictions
        return reader

    def test_gene_without_transcripts(self, mock_reader) -> None:
        """Test scoring a gene without transcripts."""
        from helixforge.io.gff import GeneModel

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=200,
            strand="+",
            transcripts=[],
            source="helixer",
        )

        calc = ConfidenceCalculator(mock_reader)
        score = calc.score_gene(gene)

        assert score.gene_id == "gene1"
        assert len(score.exon_scores) == 0
        assert len(score.intron_scores) == 0

    def test_single_exon_gene(self, mock_reader) -> None:
        """Test scoring a single-exon gene."""
        from helixforge.io.gff import GeneModel, TranscriptModel

        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=300,
            strand="+",
            exons=[(100, 300)],
            cds=[(100, 300, 0)],
        )

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=300,
            strand="+",
            transcripts=[transcript],
            source="helixer",
        )

        calc = ConfidenceCalculator(mock_reader)
        score = calc.score_gene(gene)

        assert len(score.exon_scores) == 1
        assert len(score.intron_scores) == 0

    def test_very_short_gene(self, mock_reader) -> None:
        """Test scoring a very short gene."""
        from helixforge.io.gff import GeneModel, TranscriptModel

        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=110,
            strand="+",
            exons=[(100, 110)],
            cds=[(100, 110, 0)],
        )

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=110,
            strand="+",
            transcripts=[transcript],
            source="helixer",
        )

        calc = ConfidenceCalculator(mock_reader)
        score = calc.score_gene(gene)

        assert score.gene_id == "gene1"
        # Short gene with low confidence might get flagged
        assert score.confidence_class in ["high", "medium", "low"]

    def test_empty_low_confidence_regions(self, mock_reader) -> None:
        """Test when no low-confidence regions exist."""
        calc = ConfidenceCalculator(mock_reader)

        # All high confidence
        max_probs = np.ones(100) * 0.95

        regions = calc._identify_low_confidence_regions(
            max_probs, gene_start=100, threshold=0.7
        )

        assert len(regions) == 0

    def test_entire_region_low_confidence(self, mock_reader) -> None:
        """Test when entire region is low confidence."""
        calc = ConfidenceCalculator(mock_reader)

        max_probs = np.ones(100) * 0.4

        regions = calc._identify_low_confidence_regions(
            max_probs, gene_start=100, threshold=0.7
        )

        assert len(regions) == 1
        assert regions[0][0] == 100
        assert regions[0][1] == 200
