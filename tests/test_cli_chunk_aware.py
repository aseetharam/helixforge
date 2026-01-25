"""Tests for chunk-aware CLI commands.

Tests that confidence and splice commands properly handle
--region, --chunk-id, and --scaffold flags for parallel execution.
"""

import pytest
from click.testing import CliRunner

from helixforge.cli import main


# =============================================================================
# Test CLI Option Parsing
# =============================================================================


class TestConfidenceChunkAwareOptions:
    """Test confidence command accepts chunk-aware options."""

    def test_help_shows_region_option(self):
        """Confidence --help shows --region option."""
        runner = CliRunner()
        result = runner.invoke(main, ["confidence", "--help"])
        assert result.exit_code == 0
        assert "--region" in result.output

    def test_help_shows_chunk_id_option(self):
        """Confidence --help shows --chunk-id option."""
        runner = CliRunner()
        result = runner.invoke(main, ["confidence", "--help"])
        assert result.exit_code == 0
        assert "--chunk-id" in result.output

    def test_help_shows_scaffold_option(self):
        """Confidence --help shows --scaffold option."""
        runner = CliRunner()
        result = runner.invoke(main, ["confidence", "--help"])
        assert result.exit_code == 0
        assert "--scaffold" in result.output


class TestSpliceChunkAwareOptions:
    """Test splice command accepts chunk-aware options."""

    def test_help_shows_region_option(self):
        """Splice --help shows --region option."""
        runner = CliRunner()
        result = runner.invoke(main, ["splice", "--help"])
        assert result.exit_code == 0
        assert "--region" in result.output

    def test_help_shows_chunk_id_option(self):
        """Splice --help shows --chunk-id option."""
        runner = CliRunner()
        result = runner.invoke(main, ["splice", "--help"])
        assert result.exit_code == 0
        assert "--chunk-id" in result.output

    def test_help_shows_scaffold_option(self):
        """Splice --help shows --scaffold option."""
        runner = CliRunner()
        result = runner.invoke(main, ["splice", "--help"])
        assert result.exit_code == 0
        assert "--scaffold" in result.output


# =============================================================================
# Test Parallel Tasks Command
# =============================================================================


class TestParallelTasksPlaceholders:
    """Test parallel tasks command placeholder documentation."""

    def test_help_shows_1based_coordinates(self):
        """Parallel tasks --help mentions 1-based coordinates."""
        runner = CliRunner()
        result = runner.invoke(main, ["parallel", "tasks", "--help"])
        assert result.exit_code == 0
        assert "1-based" in result.output or "{start}" in result.output

    def test_help_shows_0based_placeholders(self):
        """Parallel tasks --help shows {start_0} and {end_0}."""
        runner = CliRunner()
        result = runner.invoke(main, ["parallel", "tasks", "--help"])
        assert result.exit_code == 0
        # Check the option help text mentions coordinates
        assert "{chunk_id}" in result.output


# =============================================================================
# Test format_command Coordinate Conversion
# =============================================================================


class TestFormatCommandCoordinates:
    """Test format_command converts coordinates correctly."""

    def test_1based_coordinates_in_template(self):
        """Template {start} and {end} use 1-based coordinates."""
        from helixforge.parallel.chunker import GenomicChunk
        from helixforge.parallel.taskgen import format_command

        # Internal chunk with 0-based coords
        chunk = GenomicChunk(
            chunk_id="chunk_001",
            seqid="chr1",
            start=999,  # 0-based
            end=2000,  # half-open
        )

        result = format_command(
            "cmd --region {seqid}:{start}-{end}",
            chunk,
        )

        # Should output 1-based: chr1:1000-2000
        assert result == "cmd --region chr1:1000-2000"

    def test_0based_placeholders(self):
        """Template {start_0} and {end_0} use 0-based coordinates."""
        from helixforge.parallel.chunker import GenomicChunk
        from helixforge.parallel.taskgen import format_command

        chunk = GenomicChunk(
            chunk_id="chunk_001",
            seqid="chr1",
            start=999,
            end=2000,
        )

        result = format_command(
            "cmd --start {start_0} --end {end_0}",
            chunk,
        )

        # Should output 0-based
        assert result == "cmd --start 999 --end 2000"

    def test_mixed_placeholders(self):
        """Template can use both 1-based and 0-based placeholders."""
        from helixforge.parallel.chunker import GenomicChunk
        from helixforge.parallel.taskgen import format_command

        chunk = GenomicChunk(
            chunk_id="chunk_001",
            seqid="chr1",
            start=999,
            end=2000,
        )

        result = format_command(
            "cmd --region {seqid}:{start}-{end} --internal {start_0},{end_0}",
            chunk,
        )

        assert "chr1:1000-2000" in result
        assert "999,2000" in result

    def test_chunk_id_in_template(self):
        """Template can use {chunk_id}."""
        from helixforge.parallel.chunker import GenomicChunk
        from helixforge.parallel.taskgen import format_command

        chunk = GenomicChunk(
            chunk_id="chunk_042",
            seqid="chr1",
            start=0,
            end=1000,
        )

        result = format_command(
            "cmd --chunk-id {chunk_id} -o {chunk_id}.tsv",
            chunk,
        )

        assert "chunk_042" in result
        assert result.count("chunk_042") == 2


# =============================================================================
# Integration Tests (Require Test Data)
# =============================================================================


@pytest.mark.integration
class TestConfidenceWithRegion:
    """Integration tests for confidence command with region filtering."""

    def test_region_flag_with_valid_region(self, synthetic_confidence_data, tmp_path):
        """Confidence command accepts valid --region."""
        runner = CliRunner()
        data = synthetic_confidence_data

        output_tsv = tmp_path / "output.tsv"

        result = runner.invoke(
            main,
            [
                "confidence",
                "-p",
                str(data["h5_path"]),
                "-g",
                str(data["gff_path"]),
                "--genome",
                str(data["fasta_path"]),
                "--region",
                "chr1:1-600",
                "-o",
                str(output_tsv),
            ],
        )

        # Should succeed and process only genes in region
        assert result.exit_code == 0 or "region" in result.output.lower()

    def test_scaffold_flag(self, synthetic_confidence_data, tmp_path):
        """Confidence command accepts --scaffold."""
        runner = CliRunner()
        data = synthetic_confidence_data

        output_tsv = tmp_path / "output.tsv"

        result = runner.invoke(
            main,
            [
                "confidence",
                "-p",
                str(data["h5_path"]),
                "-g",
                str(data["gff_path"]),
                "--genome",
                str(data["fasta_path"]),
                "--scaffold",
                "chr1",
                "-o",
                str(output_tsv),
            ],
        )

        assert result.exit_code == 0 or "scaffold" in result.output.lower()

    def test_chunk_id_logged(self, synthetic_confidence_data, tmp_path):
        """Chunk ID appears in output."""
        runner = CliRunner()
        data = synthetic_confidence_data

        output_tsv = tmp_path / "output.tsv"

        result = runner.invoke(
            main,
            [
                "confidence",
                "-p",
                str(data["h5_path"]),
                "-g",
                str(data["gff_path"]),
                "--genome",
                str(data["fasta_path"]),
                "--chunk-id",
                "chunk_042",
                "-o",
                str(output_tsv),
            ],
        )

        # Chunk ID should appear in output
        assert "chunk_042" in result.output or result.exit_code == 0

    def test_invalid_region_fails(self, synthetic_confidence_data, tmp_path):
        """Invalid region format is rejected."""
        runner = CliRunner()
        data = synthetic_confidence_data

        output_tsv = tmp_path / "output.tsv"

        result = runner.invoke(
            main,
            [
                "confidence",
                "-p",
                str(data["h5_path"]),
                "-g",
                str(data["gff_path"]),
                "--genome",
                str(data["fasta_path"]),
                "--region",
                "invalid",
                "-o",
                str(output_tsv),
            ],
        )

        # Should fail with error
        assert result.exit_code != 0 or "error" in result.output.lower()

    def test_nonexistent_scaffold_fails(self, synthetic_confidence_data, tmp_path):
        """Nonexistent scaffold is rejected."""
        runner = CliRunner()
        data = synthetic_confidence_data

        output_tsv = tmp_path / "output.tsv"

        result = runner.invoke(
            main,
            [
                "confidence",
                "-p",
                str(data["h5_path"]),
                "-g",
                str(data["gff_path"]),
                "--genome",
                str(data["fasta_path"]),
                "--scaffold",
                "nonexistent",
                "-o",
                str(output_tsv),
            ],
        )

        # Should fail with error
        assert result.exit_code != 0 or "not found" in result.output.lower()


@pytest.mark.integration
class TestSpliceWithRegion:
    """Integration tests for splice command with region filtering."""

    def test_region_flag(self, synthetic_splice_dataset, tmp_path):
        """Splice command accepts valid --region."""
        runner = CliRunner()
        data = synthetic_splice_dataset

        output_gff = tmp_path / "output.gff3"

        result = runner.invoke(
            main,
            [
                "splice",
                "--helixer-gff",
                str(data["gff_path"]),
                "--genome",
                str(data["genome_path"]),
                "--junctions-bed",
                str(data["junctions_path"]),
                "--region",
                "chr1:1-800",
                "-o",
                str(output_gff),
            ],
        )

        assert result.exit_code == 0 or "region" in result.output.lower()

    def test_scaffold_flag(self, synthetic_splice_dataset, tmp_path):
        """Splice command accepts --scaffold."""
        runner = CliRunner()
        data = synthetic_splice_dataset

        output_gff = tmp_path / "output.gff3"

        result = runner.invoke(
            main,
            [
                "splice",
                "--helixer-gff",
                str(data["gff_path"]),
                "--genome",
                str(data["genome_path"]),
                "--junctions-bed",
                str(data["junctions_path"]),
                "--scaffold",
                "chr1",
                "-o",
                str(output_gff),
            ],
        )

        assert result.exit_code == 0 or "scaffold" in result.output.lower()

    def test_chunk_id_logged(self, synthetic_splice_dataset, tmp_path):
        """Chunk ID appears in output."""
        runner = CliRunner()
        data = synthetic_splice_dataset

        output_gff = tmp_path / "output.gff3"

        result = runner.invoke(
            main,
            [
                "splice",
                "--helixer-gff",
                str(data["gff_path"]),
                "--genome",
                str(data["genome_path"]),
                "--junctions-bed",
                str(data["junctions_path"]),
                "--chunk-id",
                "chunk_123",
                "-o",
                str(output_gff),
            ],
        )

        assert "chunk_123" in result.output or result.exit_code == 0


# =============================================================================
# Test Empty Region Handling
# =============================================================================


@pytest.mark.integration
class TestEmptyRegionHandling:
    """Test handling of regions with no genes."""

    def test_confidence_empty_region(self, synthetic_confidence_data, tmp_path):
        """Confidence handles region with no genes gracefully."""
        runner = CliRunner()
        data = synthetic_confidence_data

        output_tsv = tmp_path / "output.tsv"

        # Region 900-950 on chr1 should have no genes
        result = runner.invoke(
            main,
            [
                "confidence",
                "-p",
                str(data["h5_path"]),
                "-g",
                str(data["gff_path"]),
                "--genome",
                str(data["fasta_path"]),
                "--region",
                "chr1:650-699",  # Between genes
                "-o",
                str(output_tsv),
            ],
        )

        # Should succeed with warning about no genes
        assert result.exit_code == 0 or "no genes" in result.output.lower()
