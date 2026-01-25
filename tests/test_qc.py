"""Tests for the QC module.

Tests for flags, aggregation, filtering, and report generation.
"""

import json
import tempfile
from pathlib import Path

import pytest


class TestQCFlag:
    """Tests for QCFlag class."""

    def test_flag_creation(self):
        """Test QCFlag creation."""
        from helixforge.qc import FlagCategory, FlagSeverity, QCFlag

        flag = QCFlag(
            code="TEST_FLAG",
            name="Test Flag",
            description="A test flag",
            category=FlagCategory.CONFIDENCE,
            severity=FlagSeverity.WARNING,
        )

        assert flag.code == "TEST_FLAG"
        assert flag.name == "Test Flag"
        assert flag.category == FlagCategory.CONFIDENCE
        assert flag.severity == FlagSeverity.WARNING

    def test_flag_hashable(self):
        """Test that QCFlag is hashable for use in sets."""
        from helixforge.qc import Flags

        # Should work in sets
        flag_set = {Flags.LOW_CONFIDENCE, Flags.HIGH_ENTROPY}
        assert len(flag_set) == 2

        # Add same flag twice
        flag_set.add(Flags.LOW_CONFIDENCE)
        assert len(flag_set) == 2

    def test_flag_equality(self):
        """Test QCFlag equality based on code."""
        from helixforge.qc import FlagCategory, FlagSeverity, QCFlag

        flag1 = QCFlag(
            code="SAME",
            name="Name1",
            description="Desc1",
            category=FlagCategory.CONFIDENCE,
            severity=FlagSeverity.WARNING,
        )
        flag2 = QCFlag(
            code="SAME",
            name="Name2",
            description="Desc2",
            category=FlagCategory.SPLICE,
            severity=FlagSeverity.ERROR,
        )

        assert flag1 == flag2


class TestFlags:
    """Tests for Flags registry."""

    def test_get_all_flags(self):
        """Test getting all registered flags."""
        from helixforge.qc import Flags

        all_flags = Flags.get_all()
        assert len(all_flags) > 0
        assert all(hasattr(f, "code") for f in all_flags)

    def test_get_by_category(self):
        """Test filtering flags by category."""
        from helixforge.qc import FlagCategory, Flags

        confidence_flags = Flags.get_by_category(FlagCategory.CONFIDENCE)
        assert len(confidence_flags) > 0
        assert all(f.category == FlagCategory.CONFIDENCE for f in confidence_flags)

    def test_get_by_severity(self):
        """Test filtering flags by severity."""
        from helixforge.qc import FlagSeverity, Flags

        error_flags = Flags.get_by_severity(FlagSeverity.ERROR)
        assert len(error_flags) > 0
        assert all(f.severity == FlagSeverity.ERROR for f in error_flags)

    def test_get_by_code(self):
        """Test getting flag by code."""
        from helixforge.qc import Flags

        flag = Flags.get_by_code("LOW_CONF")
        assert flag is not None
        assert flag.code == "LOW_CONF"

        # Unknown code
        unknown = Flags.get_by_code("NONEXISTENT")
        assert unknown is None


class TestFlagSeverity:
    """Tests for FlagSeverity comparison."""

    def test_severity_ordering(self):
        """Test severity comparison operators."""
        from helixforge.qc import FlagSeverity

        assert FlagSeverity.INFO < FlagSeverity.WARNING
        assert FlagSeverity.WARNING < FlagSeverity.ERROR
        assert FlagSeverity.ERROR < FlagSeverity.CRITICAL

        assert FlagSeverity.CRITICAL > FlagSeverity.ERROR
        assert FlagSeverity.WARNING >= FlagSeverity.INFO
        assert FlagSeverity.ERROR <= FlagSeverity.CRITICAL


class TestGeneQC:
    """Tests for GeneQC class."""

    def test_gene_qc_creation(self):
        """Test GeneQC creation."""
        from helixforge.qc import GeneQC

        qc = GeneQC(gene_id="gene001")
        assert qc.gene_id == "gene001"
        assert qc.tier == "unclassified"
        assert len(qc.flags) == 0

    def test_add_flags(self):
        """Test adding flags to GeneQC."""
        from helixforge.qc import Flags, GeneQC

        qc = GeneQC(gene_id="gene001")
        qc.add_flag(Flags.LOW_CONFIDENCE)
        qc.add_flag(Flags.WEAK_EXON)

        assert qc.flag_count == 2
        assert qc.has_flag(Flags.LOW_CONFIDENCE)
        assert qc.has_flag(Flags.WEAK_EXON)
        assert not qc.has_flag(Flags.INTERNAL_STOP)

    def test_max_severity(self):
        """Test max_severity property."""
        from helixforge.qc import FlagSeverity, Flags, GeneQC

        qc = GeneQC(gene_id="gene001")
        assert qc.max_severity is None

        qc.add_flag(Flags.LOW_CONFIDENCE)  # WARNING
        assert qc.max_severity == FlagSeverity.WARNING

        qc.add_flag(Flags.INTERNAL_STOP)  # CRITICAL
        assert qc.max_severity == FlagSeverity.CRITICAL

    def test_classify_tier_high(self):
        """Test tier classification - high."""
        from helixforge.qc import GeneQC

        qc = GeneQC(gene_id="gene001", confidence_score=0.90)
        qc.classify_tier()
        assert qc.tier == "high"

    def test_classify_tier_medium(self):
        """Test tier classification - medium."""
        from helixforge.qc import Flags, GeneQC

        # By confidence
        qc = GeneQC(gene_id="gene001", confidence_score=0.75)
        qc.classify_tier()
        assert qc.tier == "medium"

        # By warning flag
        qc2 = GeneQC(gene_id="gene002", confidence_score=0.90)
        qc2.add_flag(Flags.LOW_CONFIDENCE)  # WARNING severity
        qc2.classify_tier()
        assert qc2.tier == "medium"

    def test_classify_tier_low(self):
        """Test tier classification - low."""
        from helixforge.qc import Flags, GeneQC

        # By confidence
        qc = GeneQC(gene_id="gene001", confidence_score=0.55)
        qc.classify_tier()
        assert qc.tier == "low"

        # By error flag
        qc2 = GeneQC(gene_id="gene002", confidence_score=0.90)
        qc2.add_flag(Flags.CHIMERIC)  # ERROR severity
        qc2.classify_tier()
        assert qc2.tier == "low"

    def test_classify_tier_reject(self):
        """Test tier classification - reject."""
        from helixforge.qc import Flags, GeneQC

        # By confidence
        qc = GeneQC(gene_id="gene001", confidence_score=0.40)
        qc.classify_tier()
        assert qc.tier == "reject"

        # By critical flag
        qc2 = GeneQC(gene_id="gene002", confidence_score=0.90)
        qc2.add_flag(Flags.INTERNAL_STOP)  # CRITICAL severity
        qc2.classify_tier()
        assert qc2.tier == "reject"

    def test_to_dict(self):
        """Test GeneQC serialization."""
        from helixforge.qc import Flags, GeneQC

        qc = GeneQC(gene_id="gene001", confidence_score=0.85)
        qc.add_flag(Flags.LOW_CONFIDENCE)
        qc.classify_tier()

        data = qc.to_dict()
        assert data["gene_id"] == "gene001"
        assert data["confidence_score"] == 0.85
        assert "LOW_CONF" in data["flag_codes"]

    def test_from_dict(self):
        """Test GeneQC deserialization."""
        from helixforge.qc import GeneQC

        data = {
            "gene_id": "gene001",
            "tier": "medium",
            "confidence_score": 0.75,
            "flag_codes": ["LOW_CONF"],
        }

        qc = GeneQC.from_dict(data)
        assert qc.gene_id == "gene001"
        assert qc.tier == "medium"
        assert qc.confidence_score == 0.75
        assert qc.flag_count == 1


class TestQCAggregator:
    """Tests for QCAggregator."""

    def test_aggregate_confidence(self):
        """Test aggregating confidence results."""
        from helixforge.qc import QCAggregator

        aggregator = QCAggregator()

        confidence_results = {
            "gene001": {"overall_score": 0.90, "exon_min": 0.85, "entropy": 0.5},
            "gene002": {"overall_score": 0.45, "exon_min": 0.30, "entropy": 2.0},
        }

        gene_qcs = aggregator.aggregate(confidence_results=confidence_results)

        assert len(gene_qcs) == 2
        assert gene_qcs["gene001"].confidence_score == 0.90
        assert gene_qcs["gene001"].tier == "high"
        assert gene_qcs["gene002"].tier == "reject"  # Very low confidence

    def test_aggregate_from_files(self):
        """Test aggregating from TSV files."""
        from helixforge.qc import QCAggregator

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test TSV
            conf_path = Path(tmpdir) / "confidence.tsv"
            with open(conf_path, "w") as f:
                f.write("gene_id\toverall_score\texon_min\n")
                f.write("gene001\t0.90\t0.85\n")
                f.write("gene002\t0.65\t0.55\n")

            aggregator = QCAggregator()
            gene_qcs = aggregator.aggregate_from_files(confidence_tsv=conf_path)

            assert len(gene_qcs) == 2
            assert gene_qcs["gene001"].confidence_score == 0.90


class TestGeneFilter:
    """Tests for GeneFilter."""

    def test_filter_by_confidence(self):
        """Test filtering by confidence threshold."""
        from helixforge.qc import FilterCriteria, GeneFilter, GeneQC

        gene_qcs = {
            "gene001": GeneQC(gene_id="gene001", confidence_score=0.90),
            "gene002": GeneQC(gene_id="gene002", confidence_score=0.70),
            "gene003": GeneQC(gene_id="gene003", confidence_score=0.50),
        }

        criteria = FilterCriteria(min_confidence=0.75)
        gene_filter = GeneFilter(criteria)
        result = gene_filter.apply(gene_qcs)

        assert result.pass_count == 1
        assert result.fail_count == 2
        assert "gene001" in result.passed_ids()

    def test_filter_by_tier(self):
        """Test filtering by tier."""
        from helixforge.qc import FilterCriteria, GeneFilter, GeneQC

        gene_qcs = {
            "gene001": GeneQC(gene_id="gene001", tier="high"),
            "gene002": GeneQC(gene_id="gene002", tier="medium"),
            "gene003": GeneQC(gene_id="gene003", tier="low"),
        }

        criteria = FilterCriteria(allowed_tiers=["high", "medium"])
        gene_filter = GeneFilter(criteria)
        result = gene_filter.apply(gene_qcs)

        assert result.pass_count == 2
        assert "gene003" not in result.passed_ids()

    def test_preset_high_confidence(self):
        """Test high_confidence preset."""
        from helixforge.qc import GeneFilter, GeneQC

        gene_qcs = {
            "gene001": GeneQC(gene_id="gene001", tier="high", confidence_score=0.90),
            "gene002": GeneQC(gene_id="gene002", tier="medium", confidence_score=0.75),
        }

        gene_filter = GeneFilter.high_confidence()
        result = gene_filter.apply(gene_qcs)

        assert result.pass_count == 1
        assert "gene001" in result.passed_ids()

    def test_filter_by_excluded_flags(self):
        """Test filtering with excluded flags."""
        from helixforge.qc import FilterCriteria, Flags, GeneFilter, GeneQC

        qc1 = GeneQC(gene_id="gene001")
        qc2 = GeneQC(gene_id="gene002")
        qc2.add_flag(Flags.INTERNAL_STOP)

        gene_qcs = {"gene001": qc1, "gene002": qc2}

        criteria = FilterCriteria(exclude_flags=[Flags.INTERNAL_STOP])
        gene_filter = GeneFilter(criteria)
        result = gene_filter.apply(gene_qcs)

        assert result.pass_count == 1
        assert "gene001" in result.passed_ids()


class TestQCReportGenerator:
    """Tests for QCReportGenerator."""

    def test_generate_html_report(self):
        """Test HTML report generation."""
        from helixforge.qc import GeneQC, QCReportGenerator

        gene_qcs = {
            "gene001": GeneQC(gene_id="gene001", tier="high", confidence_score=0.90),
            "gene002": GeneQC(gene_id="gene002", tier="medium", confidence_score=0.75),
            "gene003": GeneQC(gene_id="gene003", tier="low", confidence_score=0.55),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "report.html"
            generator = QCReportGenerator()
            generator.generate(gene_qcs, output_path, title="Test Report")

            assert output_path.exists()
            content = output_path.read_text()
            assert "Test Report" in content
            assert "Chart.js" in content
            assert "gene001" in content

    def test_generate_json_report(self):
        """Test JSON report generation."""
        from helixforge.qc import GeneQC, generate_json_report

        gene_qcs = {
            "gene001": GeneQC(gene_id="gene001", tier="high"),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "report.json"
            generate_json_report(gene_qcs, output_path)

            assert output_path.exists()
            data = json.loads(output_path.read_text())
            assert "summary" in data
            assert "genes" in data
            assert "gene001" in data["genes"]


class TestUtilityFunctions:
    """Tests for utility functions."""

    def test_summarize_flags(self):
        """Test summarize_flags function."""
        from helixforge.qc import Flags, summarize_flags

        flags = [Flags.LOW_CONFIDENCE, Flags.WEAK_EXON, Flags.INTERNAL_STOP]
        summary = summarize_flags(flags)

        assert summary["confidence"] == 2  # LOW_CONFIDENCE and WEAK_EXON
        assert summary["structure"] == 1  # INTERNAL_STOP

    def test_max_severity_function(self):
        """Test max_severity function."""
        from helixforge.qc import FlagSeverity, Flags, max_severity

        flags = [Flags.LOW_CONFIDENCE, Flags.INTERNAL_STOP]  # WARNING and CRITICAL
        result = max_severity(flags)
        assert result == FlagSeverity.CRITICAL

        # Empty list
        assert max_severity([]) is None

    def test_split_by_tier(self):
        """Test split_by_tier function."""
        from helixforge.qc import GeneQC, split_by_tier

        gene_qcs = [
            GeneQC(gene_id="gene001", tier="high"),
            GeneQC(gene_id="gene002", tier="high"),
            GeneQC(gene_id="gene003", tier="medium"),
        ]

        tiered = split_by_tier(gene_qcs)
        assert len(tiered["high"]) == 2
        assert len(tiered["medium"]) == 1

    def test_summarize_qc_results(self):
        """Test summarize_qc_results function."""
        from helixforge.qc import Flags, GeneQC, summarize_qc_results

        qc1 = GeneQC(gene_id="gene001", tier="high", confidence_score=0.90)
        qc2 = GeneQC(gene_id="gene002", tier="medium", confidence_score=0.75)
        qc2.add_flag(Flags.LOW_CONFIDENCE)

        gene_qcs = {"gene001": qc1, "gene002": qc2}
        summary = summarize_qc_results(gene_qcs)

        assert summary["total_genes"] == 2
        assert summary["tier_counts"]["high"] == 1
        assert summary["tier_counts"]["medium"] == 1
        assert summary["avg_confidence"] == 0.825


class TestTieredOutput:
    """Tests for tiered output functions."""

    def test_write_tiered_gene_lists(self):
        """Test writing tiered gene lists."""
        from helixforge.qc import GeneQC, write_tiered_gene_lists

        gene_qcs = {
            "gene001": GeneQC(gene_id="gene001", tier="high"),
            "gene002": GeneQC(gene_id="gene002", tier="medium"),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            output_files = write_tiered_gene_lists(gene_qcs, tmpdir, "test")

            assert "high" in output_files
            assert "medium" in output_files

            high_content = output_files["high"].read_text()
            assert "gene001" in high_content

            medium_content = output_files["medium"].read_text()
            assert "gene002" in medium_content

    def test_export_qc_tsv(self):
        """Test exporting QC results to TSV."""
        from helixforge.qc import Flags, GeneQC, export_qc_tsv

        qc = GeneQC(gene_id="gene001", tier="high", confidence_score=0.90)
        qc.add_flag(Flags.UNCERTAIN_BOUNDARY)

        gene_qcs = {"gene001": qc}

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "qc_results.tsv"
            export_qc_tsv(gene_qcs, output_path)

            assert output_path.exists()
            content = output_path.read_text()
            assert "gene001" in content
            assert "high" in content
            assert "0.9" in content
