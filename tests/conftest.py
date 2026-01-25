"""Pytest configuration and shared fixtures for HelixForge tests.

This module contains fixtures that are shared across multiple test modules.
Fixtures are organized by category:

- Path fixtures: Provide paths to test data files
- Data fixtures: Provide loaded test data objects
- Synthetic data fixtures: Generate test data programmatically
"""

from pathlib import Path
from typing import Generator

import h5py
import numpy as np
import pytest


# =============================================================================
# Path Fixtures
# =============================================================================


@pytest.fixture
def test_data_dir() -> Path:
    """Return the path to the test data directory."""
    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(exist_ok=True)
    return data_dir


# =============================================================================
# FASTA Fixtures
# =============================================================================


@pytest.fixture
def synthetic_fasta(tmp_path: Path) -> Path:
    """Create a synthetic FASTA file for testing.

    Creates a small genome with two scaffolds:
    - chr1: 1000 bp
    - chr2: 500 bp
    """
    fasta_path = tmp_path / "test_genome.fa"

    # Generate reproducible sequences
    np.random.seed(42)

    sequences = {
        "chr1": "".join(np.random.choice(list("ACGT"), 1000)),
        "chr2": "".join(np.random.choice(list("ACGT"), 500)),
    }

    with open(fasta_path, "w") as f:
        for seqid, seq in sequences.items():
            f.write(f">{seqid}\n")
            # Write in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")

    return fasta_path


@pytest.fixture
def synthetic_fasta_with_n(tmp_path: Path) -> Path:
    """Create a FASTA file with N bases for GC content testing."""
    fasta_path = tmp_path / "test_genome_with_n.fa"

    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        # 50 bases: 10 G, 10 C, 10 A, 10 T, 10 N = 50% GC of valid bases
        f.write("GGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTNNNNNNNNNN\n")

    return fasta_path


@pytest.fixture
def synthetic_fai(synthetic_fasta: Path) -> Path:
    """Create a FAI index file for the synthetic FASTA.

    Returns path to the .fai file.
    """
    fai_path = synthetic_fasta.with_suffix(".fa.fai")

    # FAI format: name, length, offset, linebases, linewidth
    # Calculate offsets based on FASTA format
    with open(fai_path, "w") as f:
        # chr1: 1000 bp, starts at offset 6 (">chr1\n" = 6 bytes)
        # Line bases = 80, line width = 81 (80 + newline)
        f.write("chr1\t1000\t6\t80\t81\n")
        # chr2: 500 bp, offset after chr1 section
        # chr1 section: 6 + (13 lines * 81) - 1 + 1 (last line has no full padding)
        # Actually: 6 + 12*81 + 41 = 6 + 972 + 41 = 1019
        # Then ">chr2\n" = 6 more bytes = 1025
        f.write("chr2\t500\t1025\t80\t81\n")

    return fai_path


# =============================================================================
# HDF5 Fixtures
# =============================================================================


@pytest.fixture
def synthetic_hdf5(tmp_path: Path, synthetic_fai: Path) -> Path:
    """Create a synthetic Helixer HDF5 file for testing.

    Creates predictions for 1500 positions (matching chr1: 1000 + chr2: 500)
    with 4 classes (intergenic, UTR, CDS, intron).
    """
    h5_path = tmp_path / "test_predictions.h5"

    # Total length from FAI
    total_length = 1500  # chr1 (1000) + chr2 (500)

    # Generate synthetic predictions (softmax-like probabilities)
    np.random.seed(42)
    raw = np.random.random((total_length, 4)).astype(np.float32)
    # Normalize to sum to 1 (softmax-like)
    predictions = raw / raw.sum(axis=1, keepdims=True)

    with h5py.File(h5_path, "w") as f:
        f.create_dataset("predictions", data=predictions, dtype="float32")

    return h5_path


@pytest.fixture
def synthetic_hdf5_3d(tmp_path: Path) -> Path:
    """Create a synthetic 3D HDF5 file (batched format)."""
    h5_path = tmp_path / "test_predictions_3d.h5"

    # 3 samples of 500 positions each = 1500 total
    np.random.seed(42)
    raw = np.random.random((3, 500, 4)).astype(np.float32)
    predictions = raw / raw.sum(axis=2, keepdims=True)

    with h5py.File(h5_path, "w") as f:
        f.create_dataset("predictions", data=predictions, dtype="float32")

    return h5_path


@pytest.fixture
def synthetic_hdf5_alt_name(tmp_path: Path) -> Path:
    """Create HDF5 with alternative dataset name (y_pred)."""
    h5_path = tmp_path / "test_predictions_alt.h5"

    np.random.seed(42)
    raw = np.random.random((1500, 4)).astype(np.float32)
    predictions = raw / raw.sum(axis=1, keepdims=True)

    with h5py.File(h5_path, "w") as f:
        f.create_dataset("y_pred", data=predictions, dtype="float32")

    return h5_path


# =============================================================================
# GFF3 Fixtures
# =============================================================================


@pytest.fixture
def synthetic_gff3(tmp_path: Path) -> Path:
    """Create a synthetic GFF3 file for testing."""
    gff_path = tmp_path / "test_annotations.gff3"

    content = """\
##gff-version 3
##sequence-region chr1 1 1000
##sequence-region chr2 1 500
chr1\thelixer\tgene\t101\t500\t.\t+\t.\tID=gene1;Name=TestGene1
chr1\thelixer\tmRNA\t101\t500\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\thelixer\texon\t101\t200\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\thelixer\texon\t301\t500\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\thelixer\tCDS\t101\t200\t.\t+\t0\tID=cds1;Parent=mRNA1
chr1\thelixer\tCDS\t301\t500\t.\t+\t2\tID=cds2;Parent=mRNA1
chr1\thelixer\tgene\t600\t900\t.\t-\t.\tID=gene2;Name=TestGene2
chr1\thelixer\tmRNA\t600\t900\t.\t-\t.\tID=mRNA2;Parent=gene2
chr1\thelixer\texon\t600\t700\t.\t-\t.\tID=exon3;Parent=mRNA2
chr1\thelixer\texon\t800\t900\t.\t-\t.\tID=exon4;Parent=mRNA2
chr1\thelixer\tCDS\t600\t700\t.\t-\t0\tID=cds3;Parent=mRNA2
chr1\thelixer\tCDS\t800\t900\t.\t-\t1\tID=cds4;Parent=mRNA2
chr2\thelixer\tgene\t51\t250\t.\t+\t.\tID=gene3;Name=TestGene3
chr2\thelixer\tmRNA\t51\t250\t.\t+\t.\tID=mRNA3;Parent=gene3
chr2\thelixer\texon\t51\t250\t.\t+\t.\tID=exon5;Parent=mRNA3
chr2\thelixer\tCDS\t51\t250\t.\t+\t0\tID=cds5;Parent=mRNA3
"""
    gff_path.write_text(content)
    return gff_path


@pytest.fixture
def synthetic_gff3_multiline(tmp_path: Path) -> Path:
    """Create a GFF3 with multi-line attributes and special characters."""
    gff_path = tmp_path / "test_multiline.gff3"

    content = """\
##gff-version 3
chr1\thelixer\tgene\t101\t500\t.\t+\t.\tID=gene1;Name=Test%20Gene;Note=A gene with spaces%3B and semicolons
chr1\thelixer\tmRNA\t101\t500\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\thelixer\texon\t101\t500\t.\t+\t.\tID=exon1;Parent=mRNA1
"""
    gff_path.write_text(content)
    return gff_path


# =============================================================================
# BAM Fixtures
# =============================================================================


@pytest.fixture
def synthetic_bam_data() -> dict:
    """Return data for creating synthetic BAM alignments.

    Returns dictionary with:
    - reference_lengths: Dict of reference lengths
    - alignments: List of alignment dictionaries
    """
    return {
        "reference_lengths": {"chr1": 1000, "chr2": 500},
        "alignments": [
            # Spliced read on chr1: exon 100-150, intron 150-250, exon 250-300
            {
                "name": "read1",
                "reference": "chr1",
                "start": 100,
                "cigar": "50M100N50M",
                "strand": "+",
                "mapq": 60,
            },
            # Another spliced read supporting same junction
            {
                "name": "read2",
                "reference": "chr1",
                "start": 110,
                "cigar": "40M100N60M",
                "strand": "+",
                "mapq": 60,
            },
            # Different junction
            {
                "name": "read3",
                "reference": "chr1",
                "start": 400,
                "cigar": "30M50N30M",
                "strand": "-",
                "mapq": 60,
            },
            # Unspliced read
            {
                "name": "read4",
                "reference": "chr1",
                "start": 100,
                "cigar": "100M",
                "strand": "+",
                "mapq": 60,
            },
            # Read on chr2
            {
                "name": "read5",
                "reference": "chr2",
                "start": 50,
                "cigar": "50M80N50M",
                "strand": "+",
                "mapq": 60,
            },
        ],
    }


# =============================================================================
# Gene Model Fixtures
# =============================================================================


@pytest.fixture
def sample_gene_model():
    """Create a sample GeneModel for testing."""
    from helixforge.io.gff import GeneModel, TranscriptModel

    transcript = TranscriptModel(
        transcript_id="mRNA1",
        parent_gene="gene1",
        seqid="chr1",
        start=100,  # 0-based
        end=500,
        strand="+",
        exons=[(100, 200), (300, 500)],
        cds=[(100, 200, 0), (300, 500, 2)],  # (start, end, phase)
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


# =============================================================================
# Cleanup Fixtures
# =============================================================================


@pytest.fixture(autouse=True)
def cleanup_temp_files(tmp_path: Path) -> Generator[None, None, None]:
    """Clean up temporary files after each test."""
    yield
    # Cleanup happens automatically with tmp_path fixture


# =============================================================================
# Marker Configuration
# =============================================================================


def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )


# =============================================================================
# Confidence Scoring Fixtures
# =============================================================================


@pytest.fixture
def synthetic_confidence_data(tmp_path: Path) -> dict:
    """Create synthetic data for confidence scoring tests.

    Creates a complete test dataset including:
    - FASTA genome with known sequences
    - FAI index
    - HDF5 predictions with controlled probabilities
    - GFF3 annotations with various gene structures

    The predictions are designed to create specific confidence
    scenarios for testing.
    """
    data = {}
    np.random.seed(42)

    # Create FASTA file
    fasta_path = tmp_path / "genome.fa"
    chr1_seq = "ATGC" * 250  # 1000 bp
    chr2_seq = "GCTA" * 125  # 500 bp

    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(chr1_seq), 80):
            f.write(chr1_seq[i : i + 80] + "\n")
        f.write(">chr2\n")
        for i in range(0, len(chr2_seq), 80):
            f.write(chr2_seq[i : i + 80] + "\n")

    data["fasta_path"] = fasta_path
    data["chr1_seq"] = chr1_seq
    data["chr2_seq"] = chr2_seq

    # Create FAI index
    fai_path = tmp_path / "genome.fa.fai"
    with open(fai_path, "w") as f:
        f.write("chr1\t1000\t6\t80\t81\n")
        f.write("chr2\t500\t1025\t80\t81\n")

    data["fai_path"] = fai_path

    # Create HDF5 predictions with controlled confidence patterns
    h5_path = tmp_path / "predictions.h5"
    total_length = 1500  # chr1 (1000) + chr2 (500)

    # Create predictions with known patterns:
    # - High confidence region (0-200): mostly class 2 (CDS) with prob > 0.9
    # - Medium confidence region (200-400): class 2 with prob 0.7-0.85
    # - Low confidence region (400-600): mixed classes with prob 0.4-0.6
    # - Intergenic (600-700): class 0 with high confidence
    # - Another gene (700-1000): varying confidence

    predictions = np.zeros((total_length, 4), dtype=np.float32)

    # Helper to set softmax-like probabilities
    def set_probs(start, end, dominant_class, confidence):
        for i in range(start, end):
            # Base uniform
            probs = np.array([0.1, 0.1, 0.1, 0.1])
            # Boost dominant class
            probs[dominant_class] = confidence
            # Normalize
            probs = probs / probs.sum()
            predictions[i] = probs

    # High-confidence CDS region
    set_probs(0, 200, 2, 0.9)  # CDS, high conf

    # Medium confidence CDS region (gene boundary area)
    set_probs(200, 400, 2, 0.7)  # CDS, medium conf

    # Low confidence region (uncertain)
    for i in range(400, 600):
        # Random low-confidence predictions
        probs = np.random.random(4) * 0.3 + 0.1
        predictions[i] = probs / probs.sum()

    # Intergenic region (high conf)
    set_probs(600, 700, 0, 0.85)  # Intergenic

    # Second gene with intron
    set_probs(700, 800, 2, 0.88)  # CDS exon 1
    set_probs(800, 900, 3, 0.82)  # Intron
    set_probs(900, 1000, 2, 0.87)  # CDS exon 2

    # chr2 - simple high confidence gene
    set_probs(1000, 1200, 2, 0.91)  # High conf CDS
    set_probs(1200, 1500, 0, 0.88)  # Intergenic

    with h5py.File(h5_path, "w") as f:
        f.create_dataset("predictions", data=predictions, dtype="float32")

    data["h5_path"] = h5_path
    data["predictions"] = predictions

    # Create GFF3 file with test genes
    gff_path = tmp_path / "annotations.gff3"
    gff_content = """\
##gff-version 3
##sequence-region chr1 1 1000
##sequence-region chr2 1 500
chr1\thelixer\tgene\t1\t600\t.\t+\t.\tID=gene1;Name=HighConfGene
chr1\thelixer\tmRNA\t1\t600\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\thelixer\texon\t1\t200\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\thelixer\texon\t201\t400\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\thelixer\texon\t401\t600\t.\t+\t.\tID=exon3;Parent=mRNA1
chr1\thelixer\tCDS\t1\t200\t.\t+\t0\tID=cds1a;Parent=mRNA1
chr1\thelixer\tCDS\t201\t400\t.\t+\t0\tID=cds1b;Parent=mRNA1
chr1\thelixer\tCDS\t401\t600\t.\t+\t0\tID=cds1c;Parent=mRNA1
chr1\thelixer\tgene\t701\t1000\t.\t+\t.\tID=gene2;Name=SplicedGene
chr1\thelixer\tmRNA\t701\t1000\t.\t+\t.\tID=mRNA2;Parent=gene2
chr1\thelixer\texon\t701\t800\t.\t+\t.\tID=exon2a;Parent=mRNA2
chr1\thelixer\texon\t901\t1000\t.\t+\t.\tID=exon2b;Parent=mRNA2
chr1\thelixer\tCDS\t701\t800\t.\t+\t0\tID=cds2a;Parent=mRNA2
chr1\thelixer\tCDS\t901\t1000\t.\t+\t0\tID=cds2b;Parent=mRNA2
chr2\thelixer\tgene\t1\t200\t.\t+\t.\tID=gene3;Name=SimpleGene
chr2\thelixer\tmRNA\t1\t200\t.\t+\t.\tID=mRNA3;Parent=gene3
chr2\thelixer\texon\t1\t200\t.\t+\t.\tID=exon3a;Parent=mRNA3
chr2\thelixer\tCDS\t1\t200\t.\t+\t0\tID=cds3;Parent=mRNA3
"""
    gff_path.write_text(gff_content)
    data["gff_path"] = gff_path

    return data


@pytest.fixture
def sample_gene_models():
    """Create multiple sample GeneModel objects for batch testing.

    Returns a list of genes with varying structures:
    - Single-exon gene
    - Two-exon gene with intron
    - Three-exon gene
    - Negative strand gene
    """
    from helixforge.io.gff import GeneModel, TranscriptModel

    genes = []

    # Single-exon gene
    t1 = TranscriptModel(
        transcript_id="mRNA1",
        parent_gene="gene1",
        seqid="chr1",
        start=0,
        end=200,
        strand="+",
        exons=[(0, 200)],
        cds=[(0, 200, 0)],
    )
    genes.append(
        GeneModel(
            gene_id="gene1",
            seqid="chr1",
            start=0,
            end=200,
            strand="+",
            transcripts=[t1],
            source="helixer",
        )
    )

    # Two-exon gene
    t2 = TranscriptModel(
        transcript_id="mRNA2",
        parent_gene="gene2",
        seqid="chr1",
        start=700,
        end=1000,
        strand="+",
        exons=[(700, 800), (900, 1000)],
        cds=[(700, 800, 0), (900, 1000, 0)],
    )
    genes.append(
        GeneModel(
            gene_id="gene2",
            seqid="chr1",
            start=700,
            end=1000,
            strand="+",
            transcripts=[t2],
            source="helixer",
        )
    )

    # Three-exon gene on chr2
    t3 = TranscriptModel(
        transcript_id="mRNA3",
        parent_gene="gene3",
        seqid="chr2",
        start=0,
        end=200,
        strand="+",
        exons=[(0, 60), (80, 140), (160, 200)],
        cds=[(0, 60, 0), (80, 140, 0), (160, 200, 0)],
    )
    genes.append(
        GeneModel(
            gene_id="gene3",
            seqid="chr2",
            start=0,
            end=200,
            strand="+",
            transcripts=[t3],
            source="helixer",
        )
    )

    return genes


@pytest.fixture
def sample_confidence_scores() -> list:
    """Create sample GeneConfidence objects for testing.

    Returns confidence scores representing high, medium, and low
    confidence genes.
    """
    from helixforge.core.confidence import GeneConfidence

    return [
        # High confidence gene
        GeneConfidence(
            gene_id="gene1",
            seqid="chr1",
            start=0,
            end=200,
            strand="+",
            mean_prob=0.92,
            min_prob=0.85,
            median_prob=0.91,
            entropy=0.35,
            boundary_sharpness=0.88,
            coding_consistency=0.94,
            exon_scores=[0.92],
            worst_exon_idx=0,
            worst_exon_score=0.92,
            confidence_class="high",
            flags=[],
        ),
        # Medium confidence gene
        GeneConfidence(
            gene_id="gene2",
            seqid="chr1",
            start=700,
            end=1000,
            strand="+",
            mean_prob=0.78,
            min_prob=0.65,
            median_prob=0.76,
            entropy=0.72,
            boundary_sharpness=0.75,
            coding_consistency=0.80,
            exon_scores=[0.82, 0.74],
            intron_scores=[0.68],
            worst_exon_idx=1,
            worst_exon_score=0.74,
            confidence_class="medium",
            flags=["uncertain_boundary"],
        ),
        # Low confidence gene
        GeneConfidence(
            gene_id="gene3",
            seqid="chr2",
            start=0,
            end=200,
            strand="+",
            mean_prob=0.55,
            min_prob=0.35,
            median_prob=0.58,
            entropy=1.4,
            boundary_sharpness=0.45,
            coding_consistency=0.52,
            exon_scores=[0.58, 0.42, 0.65],
            intron_scores=[0.40, 0.38],
            worst_exon_idx=1,
            worst_exon_score=0.42,
            low_confidence_regions=[(80, 140, 0.42)],
            confidence_class="low",
            flags=["weak_exon", "high_entropy", "uncertain_boundary"],
        ),
    ]


@pytest.fixture
def mock_hdf5_reader():
    """Create a mock HDF5 reader for confidence testing.

    The mock returns predictions with known patterns for testing.
    """
    from unittest.mock import MagicMock

    reader = MagicMock()
    reader.scaffold_names = ["chr1", "chr2"]
    reader.scaffold_lengths = {"chr1": 1000, "chr2": 500}

    def mock_predictions(seqid, start, end):
        """Return controlled predictions based on position."""
        length = end - start

        # Map positions to confidence levels based on test design
        if seqid == "chr1":
            offset = start
        else:
            offset = start + 1000

        predictions = np.zeros((length, 4), dtype=np.float32)

        for i in range(length):
            pos = offset + i

            if pos < 200:
                # High confidence CDS
                probs = [0.03, 0.02, 0.92, 0.03]
            elif pos < 400:
                # Medium confidence CDS
                probs = [0.08, 0.05, 0.77, 0.10]
            elif pos < 600:
                # Low confidence (uncertain)
                np.random.seed(pos)
                probs = np.random.random(4) * 0.3 + 0.1
                probs = probs / probs.sum()
            elif pos < 700:
                # Intergenic
                probs = [0.85, 0.05, 0.05, 0.05]
            elif pos < 800:
                # CDS exon
                probs = [0.04, 0.03, 0.88, 0.05]
            elif pos < 900:
                # Intron
                probs = [0.05, 0.05, 0.08, 0.82]
            elif pos < 1000:
                # CDS exon
                probs = [0.04, 0.03, 0.87, 0.06]
            elif pos < 1200:
                # chr2 high conf CDS
                probs = [0.03, 0.02, 0.91, 0.04]
            else:
                # chr2 intergenic
                probs = [0.88, 0.04, 0.04, 0.04]

            if isinstance(probs, list):
                predictions[i] = np.array(probs, dtype=np.float32)
            else:
                predictions[i] = probs.astype(np.float32)

        return predictions

    reader.get_predictions_for_region.side_effect = mock_predictions
    reader.__enter__ = MagicMock(return_value=reader)
    reader.__exit__ = MagicMock(return_value=False)

    return reader


# =============================================================================
# Splice Refinement Fixtures
# =============================================================================


@pytest.fixture
def synthetic_splice_dataset(tmp_path: Path) -> dict:
    """Create dataset with known splice site patterns.

    Includes:
    - Gene with perfect GT-AG introns matching junctions
    - Gene with shifted splice sites (correctable)
    - Gene with non-canonical splice sites
    - Gene with no RNA-seq support
    - Junctions at various positions

    Returns a dict with:
    - genome_path: Path to FASTA file
    - gff_path: Path to GFF3 file
    - junctions_path: Path to junctions BED file
    - expected_corrections: List of expected corrections
    """
    data = {}

    # Create genome FASTA with known splice sites
    # Gene 1: Perfect GT-AG at 200-300
    # Gene 2: Shifted GT-AG (should be corrected)
    # Gene 3: Non-canonical AT-TC
    # Gene 4: No junction support

    seq_chr1 = (
        "A" * 100  # 0-100: upstream
        + "ATGCGATCGA" * 10  # 100-200: exon 1
        + "GTAAGT"  # 200-206: canonical donor
        + "T" * 88  # 206-294: intron body
        + "TTTCAG"  # 294-300: canonical acceptor
        + "ATGCGATCGA" * 10  # 300-400: exon 2
        + "A" * 100  # 400-500: between genes

        # Gene 2 (shifted splice sites)
        + "ATGCGATCGA" * 10  # 500-600: exon 1
        + "GTAAGT"  # 600-606: donor (gene predicts 598)
        + "T" * 88  # 606-694: intron
        + "TTTCAG"  # 694-700: acceptor (gene predicts 702)
        + "ATGCGATCGA" * 10  # 700-800: exon 2
        + "A" * 200  # 800-1000: padding

        # Gene 3 (non-canonical)
        + "ATGCGATCGA" * 10  # 1000-1100: exon 1
        + "ATAAGT"  # 1100-1106: AT donor (non-canonical)
        + "T" * 88  # 1106-1194: intron
        + "TTTTCT"  # 1194-1200: TC acceptor (non-canonical)
        + "ATGCGATCGA" * 10  # 1200-1300: exon 2
        + "A" * 700  # padding to 2000
    )

    fasta_path = tmp_path / "genome.fa"
    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(seq_chr1), 80):
            f.write(seq_chr1[i : i + 80] + "\n")

    data["genome_path"] = fasta_path

    # Create FAI index
    fai_path = tmp_path / "genome.fa.fai"
    with open(fai_path, "w") as f:
        f.write(f"chr1\t{len(seq_chr1)}\t6\t80\t81\n")

    data["fai_path"] = fai_path

    # Create GFF3 with genes
    gff_path = tmp_path / "annotations.gff3"
    gff_content = """\
##gff-version 3
##sequence-region chr1 1 2000
# Gene 1: Perfect splice sites
chr1\thelixer\tgene\t101\t400\t.\t+\t.\tID=gene1;Name=PerfectSplice
chr1\thelixer\tmRNA\t101\t400\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\thelixer\texon\t101\t200\t.\t+\t.\tID=exon1a;Parent=mRNA1
chr1\thelixer\texon\t301\t400\t.\t+\t.\tID=exon1b;Parent=mRNA1
chr1\thelixer\tCDS\t101\t200\t.\t+\t0\tID=cds1a;Parent=mRNA1
chr1\thelixer\tCDS\t301\t400\t.\t+\t0\tID=cds1b;Parent=mRNA1
# Gene 2: Shifted splice sites (prediction is 598-702, junction is 600-700)
chr1\thelixer\tgene\t501\t800\t.\t+\t.\tID=gene2;Name=ShiftedSplice
chr1\thelixer\tmRNA\t501\t800\t.\t+\t.\tID=mRNA2;Parent=gene2
chr1\thelixer\texon\t501\t598\t.\t+\t.\tID=exon2a;Parent=mRNA2
chr1\thelixer\texon\t702\t800\t.\t+\t.\tID=exon2b;Parent=mRNA2
chr1\thelixer\tCDS\t501\t598\t.\t+\t0\tID=cds2a;Parent=mRNA2
chr1\thelixer\tCDS\t702\t800\t.\t+\t0\tID=cds2b;Parent=mRNA2
# Gene 3: Non-canonical splice sites
chr1\thelixer\tgene\t1001\t1300\t.\t+\t.\tID=gene3;Name=NonCanonical
chr1\thelixer\tmRNA\t1001\t1300\t.\t+\t.\tID=mRNA3;Parent=gene3
chr1\thelixer\texon\t1001\t1100\t.\t+\t.\tID=exon3a;Parent=mRNA3
chr1\thelixer\texon\t1201\t1300\t.\t+\t.\tID=exon3b;Parent=mRNA3
chr1\thelixer\tCDS\t1001\t1100\t.\t+\t0\tID=cds3a;Parent=mRNA3
chr1\thelixer\tCDS\t1201\t1300\t.\t+\t0\tID=cds3b;Parent=mRNA3
# Gene 4: No RNA-seq support
chr1\thelixer\tgene\t1501\t1800\t.\t+\t.\tID=gene4;Name=NoSupport
chr1\thelixer\tmRNA\t1501\t1800\t.\t+\t.\tID=mRNA4;Parent=gene4
chr1\thelixer\texon\t1501\t1600\t.\t+\t.\tID=exon4a;Parent=mRNA4
chr1\thelixer\texon\t1701\t1800\t.\t+\t.\tID=exon4b;Parent=mRNA4
chr1\thelixer\tCDS\t1501\t1600\t.\t+\t0\tID=cds4a;Parent=mRNA4
chr1\thelixer\tCDS\t1701\t1800\t.\t+\t0\tID=cds4b;Parent=mRNA4
"""
    gff_path.write_text(gff_content)
    data["gff_path"] = gff_path

    # Create junctions BED file
    # BED format: chrom, start, end, name, score (read count), strand
    junctions_path = tmp_path / "junctions.bed"
    junctions_content = """\
chr1\t200\t300\tjunc1\t100\t+
chr1\t600\t700\tjunc2\t50\t+
chr1\t1100\t1200\tjunc3\t30\t+
"""
    junctions_path.write_text(junctions_content)
    data["junctions_path"] = junctions_path

    # Define expected corrections
    data["expected_corrections"] = [
        # Gene 1: No correction needed (exact match)
        {"gene_id": "gene1", "corrected": False},
        # Gene 2: Should be corrected (598->600, 702->700)
        {
            "gene_id": "gene2",
            "corrected": True,
            "donor_shift": 2,
            "acceptor_shift": -2,
        },
        # Gene 3: Non-canonical, may or may not be corrected
        {"gene_id": "gene3", "noncanonical": True},
        # Gene 4: No support
        {"gene_id": "gene4", "unsupported": True},
    ]

    return data


@pytest.fixture
def mock_splice_junction():
    """Factory fixture for creating mock SpliceJunction objects."""

    def _create_junction(
        seqid: str = "chr1",
        start: int = 100,
        end: int = 200,
        strand: str = "+",
        read_count: int = 50,
    ):
        from unittest.mock import MagicMock

        junction = MagicMock()
        junction.seqid = seqid
        junction.start = start
        junction.end = end
        junction.strand = strand
        junction.read_count = read_count
        return junction

    return _create_junction


@pytest.fixture
def sample_pwm_data() -> dict:
    """Create sample PWM data for testing."""
    # Donor: 9 positions (-3 to +6)
    donor_matrix = [
        [0.35, 0.35, 0.15, 0.15],  # -3
        [0.60, 0.10, 0.20, 0.10],  # -2
        [0.10, 0.05, 0.80, 0.05],  # -1
        [0.00, 0.00, 1.00, 0.00],  # +1 (G)
        [0.00, 0.00, 0.00, 1.00],  # +2 (T)
        [0.55, 0.05, 0.25, 0.15],  # +3
        [0.65, 0.05, 0.15, 0.15],  # +4
        [0.15, 0.05, 0.60, 0.20],  # +5
        [0.20, 0.10, 0.20, 0.50],  # +6
    ]

    # Acceptor: 16 positions (-14 to +1)
    acceptor_matrix = [
        [0.15, 0.35, 0.15, 0.35],  # -14
        [0.15, 0.35, 0.15, 0.35],  # -13
        [0.15, 0.35, 0.15, 0.35],  # -12
        [0.15, 0.35, 0.15, 0.35],  # -11
        [0.15, 0.35, 0.15, 0.35],  # -10
        [0.15, 0.35, 0.15, 0.35],  # -9
        [0.15, 0.35, 0.15, 0.35],  # -8
        [0.12, 0.38, 0.12, 0.38],  # -7
        [0.12, 0.38, 0.12, 0.38],  # -6
        [0.12, 0.38, 0.12, 0.38],  # -5
        [0.15, 0.40, 0.10, 0.35],  # -4
        [0.25, 0.25, 0.25, 0.25],  # -3
        [0.15, 0.60, 0.10, 0.15],  # -2
        [0.70, 0.10, 0.10, 0.10],  # -1 (A)
        [0.00, 0.00, 1.00, 0.00],  # 0 (G)
        [0.15, 0.15, 0.50, 0.20],  # +1
    ]

    return {
        "donor": {
            "matrix": donor_matrix,
            "site_type": "donor",
            "offset": 3,
            "positions": [-3, -2, -1, 1, 2, 3, 4, 5, 6],
        },
        "acceptor": {
            "matrix": acceptor_matrix,
            "site_type": "acceptor",
            "offset": 14,
            "positions": list(range(-14, 2)),
        },
    }
