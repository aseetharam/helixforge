"""Tests for the homology module.

Tests for homology search, validation, and database management.
"""

import pytest
from pathlib import Path

from helixforge.homology.search import (
    HomologyHit,
    GeneHomology,
    SearchTool,
    translate_sequence,
    CODON_TABLE,
)
from helixforge.homology.validate import (
    HomologyStatus,
    HomologyValidator,
    ValidationThresholds,
    ValidationResult,
    ChimericEvidence,
    FragmentGroup,
    TEInterval,
    calculate_te_overlap,
    summarize_validation,
    get_genes_by_status,
    get_flagged_genes,
)
from helixforge.homology.databases import (
    DatabaseType,
    DatabaseInfo,
    list_databases,
)


# =============================================================================
# Test HomologyHit
# =============================================================================


class TestHomologyHit:
    """Tests for HomologyHit class."""

    def test_create_hit(self):
        """Create a basic HomologyHit."""
        hit = HomologyHit(
            query_id="gene1",
            subject_id="sp|P12345|PROT_ARATH",
            identity=85.0,
            alignment_length=300,
            mismatches=45,
            gap_opens=2,
            query_start=1,
            query_end=300,
            subject_start=1,
            subject_end=298,
            evalue=1e-100,
            bitscore=450.0,
            query_length=350,
            subject_length=320,
        )
        assert hit.query_id == "gene1"
        assert hit.identity == 85.0
        assert hit.evalue == 1e-100

    def test_query_coverage(self):
        """Test query coverage calculation."""
        hit = HomologyHit(
            query_id="gene1",
            subject_id="prot1",
            identity=90.0,
            alignment_length=200,
            mismatches=20,
            gap_opens=1,
            query_start=51,
            query_end=250,
            subject_start=1,
            subject_end=200,
            evalue=1e-50,
            bitscore=300.0,
            query_length=400,
            subject_length=200,
        )
        # Aligned: 250 - 51 + 1 = 200
        # Coverage: 200 / 400 = 0.5
        assert hit.query_coverage == pytest.approx(0.5, rel=0.01)

    def test_subject_coverage(self):
        """Test subject coverage calculation."""
        hit = HomologyHit(
            query_id="gene1",
            subject_id="prot1",
            identity=90.0,
            alignment_length=200,
            mismatches=20,
            gap_opens=1,
            query_start=1,
            query_end=200,
            subject_start=1,
            subject_end=150,
            evalue=1e-50,
            bitscore=300.0,
            query_length=200,
            subject_length=200,
        )
        # Aligned: 150 - 1 + 1 = 150
        # Coverage: 150 / 200 = 0.75
        assert hit.subject_coverage == pytest.approx(0.75, rel=0.01)

    def test_coverage_none_when_no_length(self):
        """Coverage returns None when length not available."""
        hit = HomologyHit(
            query_id="gene1",
            subject_id="prot1",
            identity=90.0,
            alignment_length=200,
            mismatches=20,
            gap_opens=1,
            query_start=1,
            query_end=200,
            subject_start=1,
            subject_end=200,
            evalue=1e-50,
            bitscore=300.0,
        )
        assert hit.query_coverage is None
        assert hit.subject_coverage is None

    def test_from_diamond_line(self):
        """Parse Diamond output line."""
        line = "gene1\tsp|P12345|PROT\t95.2\t250\t12\t1\t1\t250\t5\t254\t1e-120\t500\t260\t280"
        hit = HomologyHit.from_diamond_line(line)

        assert hit.query_id == "gene1"
        assert hit.subject_id == "sp|P12345|PROT"
        assert hit.identity == 95.2
        assert hit.alignment_length == 250
        assert hit.query_start == 1
        assert hit.query_end == 250
        assert hit.evalue == 1e-120
        assert hit.query_length == 260
        assert hit.subject_length == 280

    def test_from_diamond_line_minimal(self):
        """Parse Diamond line with 12 fields only."""
        line = "gene1\tprot1\t80.0\t100\t20\t2\t10\t109\t1\t100\t1e-30\t200"
        hit = HomologyHit.from_diamond_line(line)

        assert hit.query_id == "gene1"
        assert hit.identity == 80.0
        assert hit.query_length is None
        assert hit.subject_length is None

    def test_to_dict(self):
        """Convert hit to dictionary."""
        hit = HomologyHit(
            query_id="gene1",
            subject_id="prot1",
            identity=90.0,
            alignment_length=200,
            mismatches=20,
            gap_opens=1,
            query_start=1,
            query_end=200,
            subject_start=1,
            subject_end=200,
            evalue=1e-50,
            bitscore=300.0,
            query_length=200,
            subject_length=200,
        )
        d = hit.to_dict()
        assert d["query_id"] == "gene1"
        assert d["identity"] == 90.0
        assert "query_coverage" in d
        assert "subject_coverage" in d


# =============================================================================
# Test GeneHomology
# =============================================================================


class TestGeneHomology:
    """Tests for GeneHomology class."""

    def test_from_hits_empty(self):
        """Create GeneHomology with no hits."""
        homology = GeneHomology.from_hits("gene1", [])
        assert homology.gene_id == "gene1"
        assert homology.n_hits == 0
        assert homology.has_homology is False
        assert homology.status == "no_hit"

    def test_from_hits_with_hits(self):
        """Create GeneHomology with hits."""
        hits = [
            HomologyHit(
                query_id="gene1", subject_id="prot1",
                identity=90.0, alignment_length=200, mismatches=20, gap_opens=1,
                query_start=1, query_end=200, subject_start=1, subject_end=200,
                evalue=1e-50, bitscore=300.0,
            ),
            HomologyHit(
                query_id="gene1", subject_id="prot2",
                identity=80.0, alignment_length=150, mismatches=30, gap_opens=2,
                query_start=1, query_end=150, subject_start=1, subject_end=150,
                evalue=1e-30, bitscore=200.0,
            ),
        ]
        homology = GeneHomology.from_hits("gene1", hits)

        assert homology.gene_id == "gene1"
        assert homology.n_hits == 2
        assert homology.has_homology is True
        assert homology.best_hit_id == "prot1"  # Lowest e-value
        assert homology.best_hit_evalue == 1e-50

    def test_unique_subjects(self):
        """Get unique subject IDs."""
        hits = [
            HomologyHit(
                query_id="gene1", subject_id="prot1",
                identity=90.0, alignment_length=200, mismatches=20, gap_opens=1,
                query_start=1, query_end=100, subject_start=1, subject_end=100,
                evalue=1e-50, bitscore=300.0,
            ),
            HomologyHit(
                query_id="gene1", subject_id="prot1",
                identity=80.0, alignment_length=100, mismatches=20, gap_opens=1,
                query_start=150, query_end=250, subject_start=150, subject_end=250,
                evalue=1e-30, bitscore=200.0,
            ),
            HomologyHit(
                query_id="gene1", subject_id="prot2",
                identity=70.0, alignment_length=100, mismatches=30, gap_opens=2,
                query_start=1, query_end=100, subject_start=1, subject_end=100,
                evalue=1e-20, bitscore=150.0,
            ),
        ]
        homology = GeneHomology.from_hits("gene1", hits)
        assert homology.unique_subjects == {"prot1", "prot2"}


# =============================================================================
# Test ValidationThresholds
# =============================================================================


class TestValidationThresholds:
    """Tests for ValidationThresholds class."""

    def test_default_thresholds(self):
        """Test default thresholds."""
        t = ValidationThresholds.default()
        assert t.min_identity == 30.0
        assert t.max_evalue == 1e-5
        assert t.min_query_coverage == 0.5

    def test_strict_thresholds(self):
        """Test strict thresholds."""
        t = ValidationThresholds.strict()
        assert t.min_identity == 50.0
        assert t.max_evalue == 1e-10
        assert t.complete_query_coverage == 0.9

    def test_relaxed_thresholds(self):
        """Test relaxed thresholds."""
        t = ValidationThresholds.relaxed()
        assert t.min_identity == 20.0
        assert t.max_evalue == 1e-3
        assert t.complete_query_coverage == 0.6


# =============================================================================
# Test HomologyValidator
# =============================================================================


class TestHomologyValidator:
    """Tests for HomologyValidator class."""

    def test_validate_no_hits(self):
        """Validate gene with no hits."""
        validator = HomologyValidator()
        homology = GeneHomology.from_hits("gene1", [])
        results = validator.validate_genes({"gene1": homology})

        assert "gene1" in results
        assert results["gene1"].status == HomologyStatus.NO_HIT
        assert results["gene1"].has_homology is False

    def test_validate_complete(self):
        """Validate gene with good coverage."""
        validator = HomologyValidator()
        hits = [
            HomologyHit(
                query_id="gene1", subject_id="prot1",
                identity=90.0, alignment_length=280, mismatches=28, gap_opens=1,
                query_start=1, query_end=280, subject_start=1, subject_end=280,
                evalue=1e-100, bitscore=500.0,
                query_length=300, subject_length=300,
            ),
        ]
        homology = GeneHomology.from_hits("gene1", hits)
        results = validator.validate_genes({"gene1": homology})

        assert results["gene1"].status == HomologyStatus.COMPLETE
        assert results["gene1"].is_valid is True
        assert results["gene1"].identity == 90.0

    def test_validate_partial(self):
        """Validate gene with partial coverage."""
        validator = HomologyValidator()
        hits = [
            HomologyHit(
                query_id="gene1", subject_id="prot1",
                identity=90.0, alignment_length=100, mismatches=10, gap_opens=1,
                query_start=1, query_end=100, subject_start=1, subject_end=100,
                evalue=1e-50, bitscore=250.0,
                query_length=300, subject_length=300,  # Only 33% coverage
            ),
        ]
        homology = GeneHomology.from_hits("gene1", hits)
        results = validator.validate_genes({"gene1": homology})

        assert results["gene1"].status == HomologyStatus.PARTIAL
        assert "low_query_coverage" in results["gene1"].flags

    def test_detect_chimera(self):
        """Detect chimeric gene."""
        validator = HomologyValidator()
        hits = [
            # N-terminal hit to protein A
            HomologyHit(
                query_id="gene1", subject_id="protA",
                identity=80.0, alignment_length=150, mismatches=30, gap_opens=2,
                query_start=1, query_end=150, subject_start=1, subject_end=150,
                evalue=1e-50, bitscore=300.0,
                query_length=400, subject_length=200,
            ),
            # C-terminal hit to different protein B
            HomologyHit(
                query_id="gene1", subject_id="protB",
                identity=75.0, alignment_length=150, mismatches=37, gap_opens=3,
                query_start=200, query_end=350, subject_start=1, subject_end=150,
                evalue=1e-40, bitscore=250.0,
                query_length=400, subject_length=200,
            ),
        ]
        homology = GeneHomology.from_hits("gene1", hits)
        results = validator.validate_genes({"gene1": homology})

        assert results["gene1"].status == HomologyStatus.CHIMERIC
        assert results["gene1"].chimeric_evidence is not None
        assert results["gene1"].chimeric_evidence.subject_a == "protA"
        assert results["gene1"].chimeric_evidence.subject_b == "protB"

    def test_detect_fragments(self):
        """Detect fragmented genes."""
        validator = HomologyValidator()

        # Two genes both hitting same protein in non-overlapping regions
        gene1_hits = [
            HomologyHit(
                query_id="gene1", subject_id="full_prot",
                identity=90.0, alignment_length=150, mismatches=15, gap_opens=1,
                query_start=1, query_end=150, subject_start=1, subject_end=150,
                evalue=1e-80, bitscore=400.0,
                query_length=160, subject_length=350,
            ),
        ]
        gene2_hits = [
            HomologyHit(
                query_id="gene2", subject_id="full_prot",
                identity=85.0, alignment_length=180, mismatches=27, gap_opens=2,
                query_start=1, query_end=180, subject_start=160, subject_end=340,
                evalue=1e-70, bitscore=350.0,
                query_length=190, subject_length=350,
            ),
        ]

        gene_homology = {
            "gene1": GeneHomology.from_hits("gene1", gene1_hits),
            "gene2": GeneHomology.from_hits("gene2", gene2_hits),
        }

        results = validator.validate_genes(gene_homology)

        # Both should be marked as fragmented
        assert results["gene1"].status == HomologyStatus.FRAGMENTED
        assert results["gene2"].status == HomologyStatus.FRAGMENTED
        assert results["gene1"].fragment_group is not None
        assert results["gene1"].fragment_group == results["gene2"].fragment_group


# =============================================================================
# Test TE Overlap
# =============================================================================


class TestTEOverlap:
    """Tests for TE overlap detection."""

    def test_no_overlap(self):
        """No TE overlap."""
        te_intervals = {
            "chr1": [
                TEInterval("chr1", 1000, 2000, "LTR"),
            ]
        }
        overlap, tes = calculate_te_overlap("chr1", 100, 500, te_intervals)
        assert overlap == 0.0
        assert tes == []

    def test_partial_overlap(self):
        """Partial TE overlap."""
        te_intervals = {
            "chr1": [
                TEInterval("chr1", 400, 600, "LTR"),
            ]
        }
        # Gene 100-500 overlaps TE 400-600 by 100bp
        overlap, tes = calculate_te_overlap("chr1", 100, 500, te_intervals)
        assert overlap == pytest.approx(0.25, rel=0.01)  # 100/400
        assert len(tes) == 1

    def test_full_overlap(self):
        """Gene fully within TE."""
        te_intervals = {
            "chr1": [
                TEInterval("chr1", 100, 1000, "DNA"),
            ]
        }
        overlap, tes = calculate_te_overlap("chr1", 200, 400, te_intervals)
        assert overlap == 1.0  # Fully contained
        assert len(tes) == 1

    def test_different_scaffold(self):
        """TE on different scaffold."""
        te_intervals = {
            "chr2": [
                TEInterval("chr2", 100, 1000, "LTR"),
            ]
        }
        overlap, tes = calculate_te_overlap("chr1", 100, 500, te_intervals)
        assert overlap == 0.0
        assert tes == []


# =============================================================================
# Test Summary Functions
# =============================================================================


class TestSummaryFunctions:
    """Tests for summary and utility functions."""

    def test_summarize_validation_empty(self):
        """Summarize empty results."""
        summary = summarize_validation({})
        assert summary["total"] == 0
        assert summary["pct_with_homology"] == 0.0

    def test_summarize_validation(self):
        """Summarize validation results."""
        results = {
            "gene1": ValidationResult("gene1", HomologyStatus.COMPLETE),
            "gene2": ValidationResult("gene2", HomologyStatus.PARTIAL),
            "gene3": ValidationResult("gene3", HomologyStatus.NO_HIT),
            "gene4": ValidationResult("gene4", HomologyStatus.CHIMERIC),
        }
        summary = summarize_validation(results)

        assert summary["total"] == 4
        assert summary["complete"] == 1
        assert summary["partial"] == 1
        assert summary["no_hit"] == 1
        assert summary["chimeric"] == 1
        assert summary["pct_with_homology"] == 75.0

    def test_get_genes_by_status(self):
        """Get genes by status."""
        results = {
            "gene1": ValidationResult("gene1", HomologyStatus.COMPLETE),
            "gene2": ValidationResult("gene2", HomologyStatus.NO_HIT),
            "gene3": ValidationResult("gene3", HomologyStatus.COMPLETE),
        }
        complete_genes = get_genes_by_status(results, HomologyStatus.COMPLETE)
        assert set(complete_genes) == {"gene1", "gene3"}

    def test_get_flagged_genes(self):
        """Get genes with specific flag."""
        result1 = ValidationResult("gene1", HomologyStatus.PARTIAL)
        result1.flags.append("low_query_coverage")

        result2 = ValidationResult("gene2", HomologyStatus.PARTIAL)
        result2.flags.append("low_subject_coverage")

        result3 = ValidationResult("gene3", HomologyStatus.PARTIAL)
        result3.flags.append("low_query_coverage")

        results = {"gene1": result1, "gene2": result2, "gene3": result3}
        flagged = get_flagged_genes(results, "low_query_coverage")
        assert set(flagged) == {"gene1", "gene3"}


# =============================================================================
# Test Translation
# =============================================================================


class TestTranslation:
    """Tests for sequence translation."""

    def test_translate_simple(self):
        """Translate simple sequence."""
        # ATG (M) GGT (G) AAA (K) TGA (*)
        cds = "ATGGGTAAATGA"
        protein = translate_sequence(cds)
        assert protein == "MGK*"

    def test_translate_remove_stop(self):
        """Translate and remove terminal stop."""
        cds = "ATGGGTAAATGA"
        protein = translate_sequence(cds, stop_symbol="*")
        assert protein == "MGK*"
        # Caller typically removes terminal stop
        assert protein.rstrip("*") == "MGK"

    def test_translate_with_n(self):
        """Translate sequence with N (unknown base)."""
        # ATG (M) NNN (X) AAA (K)
        cds = "ATGNNNAAATGA"
        protein = translate_sequence(cds)
        assert protein == "MXK*"

    def test_translate_partial_codon(self):
        """Partial codon at end is ignored."""
        cds = "ATGGGTAA"  # Last AA is incomplete
        protein = translate_sequence(cds)
        assert protein == "MG"


# =============================================================================
# Test Database Info
# =============================================================================


class TestDatabaseInfo:
    """Tests for DatabaseInfo class."""

    def test_create_database_info(self):
        """Create DatabaseInfo object."""
        info = DatabaseInfo(
            name="Swiss-Prot",
            db_type=DatabaseType.SWISSPROT,
            path=Path("/data/swissprot.fasta"),
            n_sequences=500000,
        )
        assert info.name == "Swiss-Prot"
        assert info.db_type == DatabaseType.SWISSPROT
        assert info.is_downloaded is True  # Path exists not checked in test

    def test_list_databases(self):
        """List available databases."""
        dbs = list_databases()
        assert "swissprot" in dbs
        assert "swissprot_plants" in dbs
        assert "uniref90" in dbs


# =============================================================================
# Test SearchTool Enum
# =============================================================================


class TestSearchTool:
    """Tests for SearchTool enum."""

    def test_search_tool_values(self):
        """Test enum values."""
        assert SearchTool.DIAMOND.value == "diamond"
        assert SearchTool.MMSEQS2.value == "mmseqs2"

    def test_search_tool_from_string(self):
        """Create from string value."""
        assert SearchTool("diamond") == SearchTool.DIAMOND
        assert SearchTool("mmseqs2") == SearchTool.MMSEQS2
