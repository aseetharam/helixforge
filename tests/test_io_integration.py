"""Integration tests for helixforge.io modules.

These tests verify that the I/O modules work together correctly
with synthetic test data that mimics real-world scenarios.

Tests cover:
- End-to-end workflows using multiple I/O modules
- Coordinate consistency across modules
- Data roundtrip preservation
"""

from pathlib import Path

import h5py
import numpy as np
import pytest


# =============================================================================
# Fixtures for Integration Tests
# =============================================================================


@pytest.fixture
def integration_test_data(tmp_path: Path) -> dict:
    """Create a complete synthetic dataset for integration testing.

    Creates:
    - FASTA file with two scaffolds
    - FAI index
    - HDF5 predictions file
    - GFF3 annotations file

    All with consistent coordinates.
    """
    data = {}

    # Create FASTA file
    fasta_path = tmp_path / "genome.fa"
    np.random.seed(42)

    chr1_seq = "".join(np.random.choice(list("ACGT"), 1000))
    chr2_seq = "".join(np.random.choice(list("ACGT"), 500))

    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(chr1_seq), 80):
            f.write(chr1_seq[i : i + 80] + "\n")
        f.write(">chr2\n")
        for i in range(0, len(chr2_seq), 80):
            f.write(chr2_seq[i : i + 80] + "\n")

    data["fasta_path"] = fasta_path

    # Create FAI index (format: name, length, offset, linebases, linewidth)
    fai_path = tmp_path / "genome.fa.fai"
    with open(fai_path, "w") as f:
        f.write("chr1\t1000\t6\t80\t81\n")
        # chr1 takes: header (6) + 13 lines * 81 - 1 (last line no newline at end)
        # Actually for proper offset calculation: 6 + 12*81 + 41 + 1 = 1020 for header
        # Simplified: just use reasonable offset
        f.write("chr2\t500\t1025\t80\t81\n")

    data["fai_path"] = fai_path

    # Create HDF5 predictions
    h5_path = tmp_path / "predictions.h5"
    total_length = 1500  # chr1 (1000) + chr2 (500)

    np.random.seed(42)
    raw = np.random.random((total_length, 4)).astype(np.float32)
    predictions = raw / raw.sum(axis=1, keepdims=True)

    with h5py.File(h5_path, "w") as f:
        f.create_dataset("predictions", data=predictions, dtype="float32")

    data["h5_path"] = h5_path

    # Create GFF3 file
    gff_path = tmp_path / "annotations.gff3"
    gff_content = """\
##gff-version 3
##sequence-region chr1 1 1000
##sequence-region chr2 1 500
chr1\thelixer\tgene\t101\t500\t.\t+\t.\tID=gene1;Name=TestGene1
chr1\thelixer\tmRNA\t101\t500\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\thelixer\texon\t101\t200\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\thelixer\texon\t301\t500\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\thelixer\tCDS\t101\t200\t.\t+\t0\tID=cds1;Parent=mRNA1
chr1\thelixer\tCDS\t301\t500\t.\t+\t0\tID=cds2;Parent=mRNA1
chr2\thelixer\tgene\t51\t250\t.\t+\t.\tID=gene2;Name=TestGene2
chr2\thelixer\tmRNA\t51\t250\t.\t+\t.\tID=mRNA2;Parent=gene2
chr2\thelixer\texon\t51\t250\t.\t+\t.\tID=exon3;Parent=mRNA2
chr2\thelixer\tCDS\t51\t250\t.\t+\t0\tID=cds3;Parent=mRNA2
"""
    gff_path.write_text(gff_content)
    data["gff_path"] = gff_path

    data["chr1_seq"] = chr1_seq
    data["chr2_seq"] = chr2_seq

    return data


# =============================================================================
# FASTA and HDF5 Integration
# =============================================================================


class TestFastaHDF5Integration:
    """Test FASTA and HDF5 modules working together."""

    def test_coordinate_consistency(self, integration_test_data: dict) -> None:
        """Verify coordinates are consistent between FASTA and HDF5."""
        from helixforge.io.fasta import GenomeAccessor
        from helixforge.io.hdf5 import HelixerHDF5Reader

        fasta_path = integration_test_data["fasta_path"]
        h5_path = integration_test_data["h5_path"]
        fai_path = integration_test_data["fai_path"]

        with GenomeAccessor(fasta_path) as genome:
            with HelixerHDF5Reader(h5_path, fai_path) as reader:
                # Check scaffolds match
                assert set(genome.scaffold_order) == set(reader.scaffold_names)

                # Check lengths match
                for seqid in genome.scaffold_order:
                    assert genome.get_length(seqid) == reader.scaffold_lengths[seqid]

    def test_extract_sequence_and_predictions(
        self, integration_test_data: dict
    ) -> None:
        """Test extracting sequence and predictions for same region."""
        from helixforge.io.fasta import GenomeAccessor
        from helixforge.io.hdf5 import HelixerHDF5Reader

        fasta_path = integration_test_data["fasta_path"]
        h5_path = integration_test_data["h5_path"]
        fai_path = integration_test_data["fai_path"]

        with GenomeAccessor(fasta_path) as genome:
            with HelixerHDF5Reader(h5_path, fai_path) as reader:
                # Extract same region from both
                region = ("chr1", 100, 200)
                seq = genome.get_sequence(*region)
                preds = reader.get_predictions_for_region(*region)

                # Verify consistent sizes
                assert len(seq) == preds.shape[0]
                assert preds.shape == (100, 4)


# =============================================================================
# GFF3 and FASTA Integration
# =============================================================================


class TestGFF3FastaIntegration:
    """Test GFF3 and FASTA modules working together."""

    def test_extract_gene_sequences(self, integration_test_data: dict) -> None:
        """Test extracting sequences for gene regions."""
        from helixforge.io.fasta import GenomeAccessor
        from helixforge.io.gff import GFF3Parser

        fasta_path = integration_test_data["fasta_path"]
        gff_path = integration_test_data["gff_path"]

        parser = GFF3Parser(gff_path)
        genes = list(parser.iter_genes())

        with GenomeAccessor(fasta_path) as genome:
            for gene in genes:
                # Extract gene sequence
                seq = genome.get_sequence(
                    gene.seqid, gene.start, gene.end, gene.strand
                )

                # Verify length matches gene model
                assert len(seq) == gene.gene_length

    def test_extract_exon_sequences(self, integration_test_data: dict) -> None:
        """Test extracting sequences for exon regions."""
        from helixforge.io.fasta import GenomeAccessor
        from helixforge.io.gff import GFF3Parser

        fasta_path = integration_test_data["fasta_path"]
        gff_path = integration_test_data["gff_path"]

        parser = GFF3Parser(gff_path)
        genes = list(parser.iter_genes())

        with GenomeAccessor(fasta_path) as genome:
            for gene in genes:
                for transcript in gene.transcripts:
                    spliced_seq = ""
                    for exon_start, exon_end in transcript.exons:
                        exon_seq = genome.get_sequence(
                            transcript.seqid, exon_start, exon_end
                        )
                        spliced_seq += exon_seq

                    # Verify spliced length matches transcript
                    assert len(spliced_seq) == transcript.transcript_length


# =============================================================================
# GFF3 and HDF5 Integration
# =============================================================================


class TestGFF3HDF5Integration:
    """Test GFF3 and HDF5 modules working together."""

    def test_get_predictions_for_genes(self, integration_test_data: dict) -> None:
        """Test extracting predictions for gene regions."""
        from helixforge.io.gff import GFF3Parser
        from helixforge.io.hdf5 import HelixerHDF5Reader

        gff_path = integration_test_data["gff_path"]
        h5_path = integration_test_data["h5_path"]
        fai_path = integration_test_data["fai_path"]

        parser = GFF3Parser(gff_path)
        genes = list(parser.iter_genes())

        with HelixerHDF5Reader(h5_path, fai_path) as reader:
            for gene in genes:
                # Get predictions for gene region
                preds = reader.get_predictions_for_region(
                    gene.seqid, gene.start, gene.end
                )

                # Verify shape
                assert preds.shape[0] == gene.gene_length
                assert preds.shape[1] == 4

                # Verify probabilities sum to ~1
                assert np.allclose(preds.sum(axis=1), 1.0, atol=1e-5)


# =============================================================================
# Full Pipeline Integration
# =============================================================================


class TestFullPipelineIntegration:
    """Test complete workflow across all I/O modules."""

    def test_gene_analysis_workflow(self, integration_test_data: dict) -> None:
        """Test a typical gene analysis workflow.

        1. Load gene annotations from GFF3
        2. Get predictions from HDF5
        3. Extract sequences from FASTA
        4. Write refined genes back to GFF3
        """
        from helixforge.io.fasta import GenomeAccessor
        from helixforge.io.gff import GFF3Parser, GFF3Writer, GeneModel
        from helixforge.io.hdf5 import HelixerHDF5Reader

        fasta_path = integration_test_data["fasta_path"]
        gff_path = integration_test_data["gff_path"]
        h5_path = integration_test_data["h5_path"]
        fai_path = integration_test_data["fai_path"]

        # Step 1: Load genes
        parser = GFF3Parser(gff_path)
        original_genes = list(parser.iter_genes())

        assert len(original_genes) == 2

        # Step 2 & 3: Analyze each gene
        refined_genes = []

        with GenomeAccessor(fasta_path) as genome:
            with HelixerHDF5Reader(h5_path, fai_path) as reader:
                for gene in original_genes:
                    # Get predictions
                    preds = reader.get_predictions_for_region(
                        gene.seqid, gene.start, gene.end
                    )

                    # Calculate mean confidence (max probability)
                    confidence = float(np.mean(np.max(preds, axis=1)))

                    # Get sequence
                    seq = genome.get_sequence(
                        gene.seqid, gene.start, gene.end, gene.strand
                    )

                    # Calculate GC content
                    gc = genome.get_gc_content(gene.seqid, gene.start, gene.end)

                    # Create refined gene with confidence
                    refined = GeneModel(
                        gene_id=gene.gene_id,
                        seqid=gene.seqid,
                        start=gene.start,
                        end=gene.end,
                        strand=gene.strand,
                        source="HelixForge",
                        transcripts=gene.transcripts,
                        confidence=confidence,
                        qc_flags=[],
                    )
                    refined_genes.append(refined)

        assert len(refined_genes) == 2
        assert all(g.confidence is not None for g in refined_genes)
        assert all(0.25 <= g.confidence <= 1.0 for g in refined_genes)

        # Step 4: Write refined genes
        output_path = integration_test_data["fasta_path"].parent / "refined.gff3"

        with GFF3Writer(output_path) as writer:
            writer.write_genes(refined_genes)

        # Verify output
        assert output_path.exists()
        content = output_path.read_text()
        assert "##gff-version 3" in content
        assert "gene1" in content
        assert "gene2" in content
        assert "confidence=" in content

    def test_roundtrip_gff3(self, integration_test_data: dict, tmp_path: Path) -> None:
        """Test GFF3 read-write roundtrip preserves structure."""
        from helixforge.io.gff import GFF3Parser, GFF3Writer

        gff_path = integration_test_data["gff_path"]
        output_path = tmp_path / "roundtrip.gff3"

        # Read original
        parser = GFF3Parser(gff_path)
        original_genes = list(parser.iter_genes())

        # Write to new file
        with GFF3Writer(output_path) as writer:
            writer.write_genes(original_genes)

        # Read back
        parser2 = GFF3Parser(output_path)
        roundtrip_genes = list(parser2.iter_genes())

        # Compare
        assert len(roundtrip_genes) == len(original_genes)

        for orig, roundtrip in zip(original_genes, roundtrip_genes):
            assert orig.gene_id == roundtrip.gene_id
            assert orig.seqid == roundtrip.seqid
            assert orig.start == roundtrip.start
            assert orig.end == roundtrip.end
            assert orig.strand == roundtrip.strand


# =============================================================================
# Edge Cases
# =============================================================================


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_base_region(self, integration_test_data: dict) -> None:
        """Test handling of single-base regions."""
        from helixforge.io.fasta import GenomeAccessor
        from helixforge.io.hdf5 import HelixerHDF5Reader

        fasta_path = integration_test_data["fasta_path"]
        h5_path = integration_test_data["h5_path"]
        fai_path = integration_test_data["fai_path"]

        with GenomeAccessor(fasta_path) as genome:
            # Single base extraction (but needs start < end)
            seq = genome.get_sequence("chr1", 100, 101)
            assert len(seq) == 1

        with HelixerHDF5Reader(h5_path, fai_path) as reader:
            preds = reader.get_predictions_for_region("chr1", 100, 101)
            assert preds.shape == (1, 4)

    def test_scaffold_boundary(self, integration_test_data: dict) -> None:
        """Test extraction at scaffold boundaries."""
        from helixforge.io.fasta import GenomeAccessor
        from helixforge.io.hdf5 import HelixerHDF5Reader

        fasta_path = integration_test_data["fasta_path"]
        h5_path = integration_test_data["h5_path"]
        fai_path = integration_test_data["fai_path"]

        with GenomeAccessor(fasta_path) as genome:
            # Extract from end of chr1
            seq = genome.get_sequence("chr1", 990, 1000)
            assert len(seq) == 10

            # Extract from start of chr2
            seq = genome.get_sequence("chr2", 0, 10)
            assert len(seq) == 10

        with HelixerHDF5Reader(h5_path, fai_path) as reader:
            # Predictions at scaffold boundaries
            preds1 = reader.get_predictions_for_region("chr1", 990, 1000)
            preds2 = reader.get_predictions_for_region("chr2", 0, 10)

            assert preds1.shape == (10, 4)
            assert preds2.shape == (10, 4)

    def test_empty_gene_list(self, integration_test_data: dict, tmp_path: Path) -> None:
        """Test writing empty gene list."""
        from helixforge.io.gff import GFF3Writer

        output_path = tmp_path / "empty.gff3"

        with GFF3Writer(output_path) as writer:
            writer.write_genes([])

        content = output_path.read_text()
        assert "##gff-version 3" in content
