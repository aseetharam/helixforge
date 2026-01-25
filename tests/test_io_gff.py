"""Unit tests for helixforge.io.gff module.

Tests cover:
- GFF3Parser initialization and parsing
- GeneModel and TranscriptModel data structures
- GFF3Writer output generation
- Coordinate conversion (1-based GFF3 to 0-based internal)
- Attribute parsing
"""

from pathlib import Path

import pytest

from helixforge.io.gff import (
    GFF3Parser,
    GFF3Writer,
    GeneModel,
    TranscriptModel,
    format_attributes,
    parse_attributes,
    read_gff,
)


# =============================================================================
# Utility Function Tests
# =============================================================================


class TestParseAttributes:
    """Tests for attribute string parsing."""

    def test_simple_attributes(self) -> None:
        """Test parsing simple attribute string."""
        attr_str = "ID=gene1;Name=TestGene"
        attrs = parse_attributes(attr_str)

        assert attrs["ID"] == "gene1"
        assert attrs["Name"] == "TestGene"

    def test_url_encoded_attributes(self) -> None:
        """Test parsing URL-encoded attributes."""
        attr_str = "ID=gene1;Note=A%3Bgene"
        attrs = parse_attributes(attr_str)

        assert attrs["Note"] == "A;gene"

    def test_empty_value(self) -> None:
        """Test attribute with empty value."""
        attr_str = "ID=gene1;Note="
        attrs = parse_attributes(attr_str)

        assert attrs["ID"] == "gene1"
        assert attrs.get("Note") == ""

    def test_multiple_values(self) -> None:
        """Test attribute with comma-separated values."""
        attr_str = "ID=gene1;Parent=mRNA1,mRNA2"
        attrs = parse_attributes(attr_str)

        # Parent should contain comma-separated value
        assert "mRNA1" in attrs["Parent"]

    def test_empty_string(self) -> None:
        """Test parsing empty attribute string."""
        attrs = parse_attributes("")
        assert attrs == {}

        attrs = parse_attributes(".")
        assert attrs == {}


class TestFormatAttributes:
    """Tests for attribute formatting."""

    def test_simple_format(self) -> None:
        """Test formatting simple attributes."""
        attrs = {"ID": "gene1", "Name": "TestGene"}
        result = format_attributes(attrs)

        assert "ID=gene1" in result
        assert "Name=TestGene" in result
        assert ";" in result

    def test_special_chars_encoded(self) -> None:
        """Test that special characters are encoded."""
        attrs = {"ID": "gene1", "Note": "test;value"}
        result = format_attributes(attrs)

        assert "%3B" in result

    def test_empty_attrs(self) -> None:
        """Test formatting empty attributes."""
        result = format_attributes({})
        assert result == "."


# =============================================================================
# TranscriptModel Tests
# =============================================================================


class TestTranscriptModel:
    """Tests for TranscriptModel data class."""

    def test_creation(self) -> None:
        """Test creating a TranscriptModel."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 200), (300, 500)],
            cds=[(100, 200, 0), (300, 500, 2)],
        )

        assert transcript.transcript_id == "mRNA1"
        assert transcript.parent_gene == "gene1"
        assert len(transcript.exons) == 2

    def test_exon_count(self) -> None:
        """Test n_exons property."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 200), (300, 500)],
        )

        assert transcript.n_exons == 2

    def test_introns_forward(self) -> None:
        """Test intron coordinate calculation for forward strand."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=600,
            strand="+",
            exons=[(100, 200), (300, 400), (500, 600)],
        )

        introns = transcript.introns
        assert introns == [(200, 300), (400, 500)]

    def test_introns_reverse(self) -> None:
        """Test intron coordinate calculation for reverse strand."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=600,
            strand="-",
            exons=[(500, 600), (300, 400), (100, 200)],
        )

        introns = transcript.introns
        # Introns are calculated from sorted exons
        assert len(introns) == 2

    def test_transcript_length(self) -> None:
        """Test transcript length calculation."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=600,
            strand="+",
            exons=[(100, 200), (300, 500)],  # 100 + 200 = 300
        )

        assert transcript.transcript_length == 300

    def test_cds_length(self) -> None:
        """Test CDS length calculation."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 300), (400, 500)],
            cds=[(150, 300, 0), (400, 450, 2)],  # 150 + 50 = 200
        )

        assert transcript.cds_length == 200

    def test_cds_length_no_cds(self) -> None:
        """Test CDS length for transcript without CDS."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 500)],
            cds=[],
        )

        assert transcript.cds_length == 0

    def test_single_exon(self) -> None:
        """Test single-exon transcript."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 500)],
        )

        assert transcript.n_exons == 1
        assert transcript.introns == []


# =============================================================================
# GeneModel Tests
# =============================================================================


class TestGeneModel:
    """Tests for GeneModel data class."""

    def test_creation(self) -> None:
        """Test creating a GeneModel."""
        transcript = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 200), (300, 500)],
        )

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[transcript],
        )

        assert gene.gene_id == "gene1"
        assert len(gene.transcripts) == 1

    def test_n_transcripts(self) -> None:
        """Test n_transcripts property."""
        t1 = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 500)],
        )
        t2 = TranscriptModel(
            transcript_id="mRNA2",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=400,
            strand="+",
            exons=[(100, 400)],
        )

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[t1, t2],
        )

        assert gene.n_transcripts == 2

    def test_gene_length(self) -> None:
        """Test gene_length property."""
        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[],
        )

        assert gene.gene_length == 400

    def test_qc_flags_default(self) -> None:
        """Test QC flags default to empty list."""
        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[],
        )

        assert gene.qc_flags == []

    def test_confidence_default(self) -> None:
        """Test confidence defaults to None."""
        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[],
        )

        assert gene.confidence is None

    def test_primary_transcript(self) -> None:
        """Test primary_transcript property."""
        t1 = TranscriptModel(
            transcript_id="mRNA1",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            exons=[(100, 500)],
            cds=[(100, 400, 0)],  # 300 bp CDS
        )
        t2 = TranscriptModel(
            transcript_id="mRNA2",
            parent_gene="gene1",
            seqid="chr1",
            start=100,
            end=400,
            strand="+",
            exons=[(100, 400)],
            cds=[(100, 200, 0)],  # 100 bp CDS
        )

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[t1, t2],
        )

        # Primary should be t1 (longer CDS)
        assert gene.primary_transcript == t1


# =============================================================================
# GFF3Parser Tests
# =============================================================================


class TestGFF3ParserInit:
    """Tests for GFF3Parser initialization."""

    def test_init_success(self, synthetic_gff3: Path) -> None:
        """Test successful initialization."""
        parser = GFF3Parser(synthetic_gff3)
        assert parser.path == synthetic_gff3

    def test_init_missing_file(self, tmp_path: Path) -> None:
        """Test error for missing file."""
        with pytest.raises(FileNotFoundError, match="GFF3 file not found"):
            GFF3Parser(tmp_path / "missing.gff3")


class TestGFF3ParserParsing:
    """Tests for GFF3Parser parsing functionality."""

    def test_iter_genes(self, synthetic_gff3: Path) -> None:
        """Test iterating over genes."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        assert len(genes) == 3
        assert genes[0].gene_id == "gene1"
        assert genes[1].gene_id == "gene2"
        assert genes[2].gene_id == "gene3"

    def test_gene_coordinates_converted(self, synthetic_gff3: Path) -> None:
        """Test that coordinates are converted to 0-based."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        # GFF3 has gene1 at 101-500 (1-based)
        # Should be 100-500 (0-based, half-open)
        gene1 = genes[0]
        assert gene1.start == 100
        assert gene1.end == 500

    def test_transcript_parsing(self, synthetic_gff3: Path) -> None:
        """Test transcript parsing."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        gene1 = genes[0]
        assert gene1.n_transcripts == 1

        transcript = gene1.transcripts[0]
        assert transcript.transcript_id == "mRNA1"
        assert len(transcript.exons) == 2

    def test_exon_coordinates(self, synthetic_gff3: Path) -> None:
        """Test exon coordinate parsing."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        transcript = genes[0].transcripts[0]
        # GFF3: exon1 101-200, exon2 301-500
        # 0-based: (100, 200), (300, 500)
        assert transcript.exons == [(100, 200), (300, 500)]

    def test_cds_coordinates(self, synthetic_gff3: Path) -> None:
        """Test CDS coordinate parsing."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        transcript = genes[0].transcripts[0]
        # CDS includes phase: (start, end, phase)
        assert len(transcript.cds) == 2
        assert transcript.cds[0][:2] == (100, 200)
        assert transcript.cds[1][:2] == (300, 500)

    def test_phase_parsing(self, synthetic_gff3: Path) -> None:
        """Test CDS phase parsing."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        transcript = genes[0].transcripts[0]
        # Phase is stored in the CDS tuples
        phases = [cds[2] for cds in transcript.cds]
        assert phases == [0, 2]

    def test_strand_parsing(self, synthetic_gff3: Path) -> None:
        """Test strand parsing."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        assert genes[0].strand == "+"
        assert genes[1].strand == "-"

    def test_source_parsing(self, synthetic_gff3: Path) -> None:
        """Test source field parsing."""
        parser = GFF3Parser(synthetic_gff3)
        genes = list(parser.iter_genes())

        assert genes[0].source == "helixer"

    def test_get_genes_in_region(self, synthetic_gff3: Path) -> None:
        """Test region-based gene retrieval."""
        parser = GFF3Parser(synthetic_gff3)
        # Region overlapping gene1 and gene2
        genes = parser.get_genes_in_region("chr1", 0, 1000)

        assert len(genes) == 2
        gene_ids = {g.gene_id for g in genes}
        assert gene_ids == {"gene1", "gene2"}

    def test_get_genes_in_region_partial_overlap(
        self, synthetic_gff3: Path
    ) -> None:
        """Test region with partial overlap."""
        parser = GFF3Parser(synthetic_gff3)
        # Region overlapping only gene1
        genes = parser.get_genes_in_region("chr1", 0, 300)

        assert len(genes) == 1
        assert genes[0].gene_id == "gene1"

    def test_get_genes_in_region_different_scaffold(
        self, synthetic_gff3: Path
    ) -> None:
        """Test region on different scaffold."""
        parser = GFF3Parser(synthetic_gff3)
        genes = parser.get_genes_in_region("chr2", 0, 500)

        assert len(genes) == 1
        assert genes[0].gene_id == "gene3"

    def test_url_encoded_attributes(self, synthetic_gff3_multiline: Path) -> None:
        """Test parsing of URL-encoded attributes."""
        parser = GFF3Parser(synthetic_gff3_multiline)
        genes = list(parser.iter_genes())

        # Name is stored in attributes dict
        assert genes[0].attributes.get("Name") == "Test Gene"


# =============================================================================
# GFF3Writer Tests
# =============================================================================


class TestGFF3Writer:
    """Tests for GFF3Writer class."""

    def test_write_single_gene(
        self, tmp_path: Path, sample_gene_model
    ) -> None:
        """Test writing a single gene."""
        output_path = tmp_path / "output.gff3"

        with GFF3Writer(output_path) as writer:
            writer.write_gene(sample_gene_model)

        # Verify output
        content = output_path.read_text()
        assert "##gff-version 3" in content
        assert "gene1" in content
        assert "mRNA1" in content

    def test_write_multiple_genes(
        self, tmp_path: Path, sample_gene_model
    ) -> None:
        """Test writing multiple genes."""
        output_path = tmp_path / "output.gff3"

        # Create second gene
        gene2 = GeneModel(
            gene_id="gene2",
            seqid="chr1",
            start=600,
            end=900,
            strand="-",
            transcripts=[],
        )

        with GFF3Writer(output_path) as writer:
            writer.write_gene(sample_gene_model)
            writer.write_gene(gene2)

        content = output_path.read_text()
        assert "gene1" in content
        assert "gene2" in content

    def test_coordinates_converted_to_1based(
        self, tmp_path: Path, sample_gene_model
    ) -> None:
        """Test that coordinates are converted to 1-based GFF3."""
        output_path = tmp_path / "output.gff3"

        with GFF3Writer(output_path) as writer:
            writer.write_gene(sample_gene_model)

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Find the gene line
        gene_line = [l for l in lines if "\tgene\t" in l][0]
        parts = gene_line.split("\t")

        # Start should be 101 (0-based 100 + 1)
        assert parts[3] == "101"
        # End should be 500 (same as 0-based end)
        assert parts[4] == "500"

    def test_write_with_confidence(self, tmp_path: Path) -> None:
        """Test writing gene with confidence score."""
        output_path = tmp_path / "output.gff3"

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[],
            confidence=0.95,
        )

        with GFF3Writer(output_path) as writer:
            writer.write_gene(gene)

        content = output_path.read_text()
        assert "confidence=0.95" in content

    def test_write_with_qc_flags(self, tmp_path: Path) -> None:
        """Test writing gene with QC flags."""
        output_path = tmp_path / "output.gff3"

        gene = GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=100,
            end=500,
            strand="+",
            transcripts=[],
            qc_flags=["low_confidence", "no_homology"],
        )

        with GFF3Writer(output_path) as writer:
            writer.write_gene(gene)

        content = output_path.read_text()
        assert "qc_flags=" in content

    def test_write_genes_iterable(
        self, tmp_path: Path, sample_gene_model
    ) -> None:
        """Test write_genes with iterable."""
        output_path = tmp_path / "output.gff3"

        genes = [sample_gene_model]

        with GFF3Writer(output_path) as writer:
            writer.write_genes(genes)

        content = output_path.read_text()
        assert "gene1" in content


# =============================================================================
# Convenience Function Tests
# =============================================================================


class TestReadGFF:
    """Tests for read_gff convenience function."""

    def test_read_gff(self, synthetic_gff3: Path) -> None:
        """Test read_gff function."""
        genes = read_gff(synthetic_gff3)

        assert len(genes) == 3
        assert all(isinstance(g, GeneModel) for g in genes)


# =============================================================================
# Roundtrip Tests
# =============================================================================


class TestGFF3Roundtrip:
    """Tests for GFF3 read/write roundtrip."""

    def test_roundtrip_preserves_structure(
        self, synthetic_gff3: Path, tmp_path: Path
    ) -> None:
        """Test that write/read preserves gene structure."""
        # Read original
        with GFF3Parser(synthetic_gff3) as parser:
            original_genes = list(parser.iter_genes())

        # Write to new file
        output_path = tmp_path / "roundtrip.gff3"
        with GFF3Writer(output_path) as writer:
            writer.write_genes(original_genes)

        # Read back
        with GFF3Parser(output_path) as parser:
            roundtrip_genes = list(parser.iter_genes())

        # Compare
        assert len(roundtrip_genes) == len(original_genes)

        for orig, roundtrip in zip(original_genes, roundtrip_genes):
            assert orig.gene_id == roundtrip.gene_id
            assert orig.seqid == roundtrip.seqid
            assert orig.start == roundtrip.start
            assert orig.end == roundtrip.end
            assert orig.strand == roundtrip.strand
            assert orig.n_transcripts == roundtrip.n_transcripts
