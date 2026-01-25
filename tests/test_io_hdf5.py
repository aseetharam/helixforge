"""Unit tests for helixforge.io.hdf5 module.

Tests cover:
- HelixerHDF5Reader initialization and schema detection
- CoordinateIndex mapping functions
- Chunk iteration and region extraction
- Error handling for invalid inputs
"""

from pathlib import Path

import h5py
import numpy as np
import pytest

from helixforge.io.hdf5 import (
    CLASS_CDS,
    CLASS_INTERGENIC,
    CLASS_INTRON,
    CLASS_NAMES,
    CLASS_UTR,
    CoordinateIndex,
    HelixerFormatError,
    HelixerHDF5Reader,
    detect_helixer_schema,
    validate_helixer_file,
)


# =============================================================================
# CoordinateIndex Tests
# =============================================================================


class TestCoordinateIndex:
    """Tests for the CoordinateIndex class."""

    def test_init_with_valid_fai(self, synthetic_fai: Path) -> None:
        """Test initialization with a valid FAI file."""
        index = CoordinateIndex(synthetic_fai)

        assert len(index.scaffold_order) == 2
        assert index.scaffold_order == ["chr1", "chr2"]
        assert index.scaffold_lengths["chr1"] == 1000
        assert index.scaffold_lengths["chr2"] == 500
        assert index.total_length == 1500

    def test_init_missing_file(self, tmp_path: Path) -> None:
        """Test initialization raises FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError, match="FAI index not found"):
            CoordinateIndex(tmp_path / "nonexistent.fai")

    def test_scaffold_offsets(self, synthetic_fai: Path) -> None:
        """Test scaffold offsets are calculated correctly."""
        index = CoordinateIndex(synthetic_fai)

        assert index.scaffold_offsets["chr1"] == 0
        assert index.scaffold_offsets["chr2"] == 1000

    def test_genomic_to_array(self, synthetic_fai: Path) -> None:
        """Test conversion from genomic to array coordinates."""
        index = CoordinateIndex(synthetic_fai)

        # chr1 region
        start, end = index.genomic_to_array("chr1", 100, 200)
        assert start == 100
        assert end == 200

        # chr2 region (offset by chr1 length)
        start, end = index.genomic_to_array("chr2", 50, 100)
        assert start == 1050  # 1000 + 50
        assert end == 1100  # 1000 + 100

    def test_genomic_to_array_invalid_seqid(self, synthetic_fai: Path) -> None:
        """Test error for unknown scaffold."""
        index = CoordinateIndex(synthetic_fai)

        with pytest.raises(KeyError, match="Unknown sequence"):
            index.genomic_to_array("chrX", 0, 100)

    def test_genomic_to_array_invalid_coords(self, synthetic_fai: Path) -> None:
        """Test error for invalid coordinates."""
        index = CoordinateIndex(synthetic_fai)

        # Negative start
        with pytest.raises(ValueError, match="Invalid coordinates"):
            index.genomic_to_array("chr1", -1, 100)

        # End beyond scaffold
        with pytest.raises(ValueError, match="Invalid coordinates"):
            index.genomic_to_array("chr1", 0, 2000)

        # Start >= end
        with pytest.raises(ValueError, match="Invalid coordinates"):
            index.genomic_to_array("chr1", 100, 100)

    def test_array_to_genomic(self, synthetic_fai: Path) -> None:
        """Test conversion from array to genomic coordinates."""
        index = CoordinateIndex(synthetic_fai)

        # Index in chr1
        seqid, pos = index.array_to_genomic(500)
        assert seqid == "chr1"
        assert pos == 500

        # Index in chr2
        seqid, pos = index.array_to_genomic(1200)
        assert seqid == "chr2"
        assert pos == 200

    def test_array_to_genomic_out_of_bounds(self, synthetic_fai: Path) -> None:
        """Test error for out of bounds array index."""
        index = CoordinateIndex(synthetic_fai)

        with pytest.raises(ValueError, match="out of bounds"):
            index.array_to_genomic(2000)

        with pytest.raises(ValueError, match="out of bounds"):
            index.array_to_genomic(-1)

    def test_get_scaffold_range(self, synthetic_fai: Path) -> None:
        """Test getting scaffold range."""
        index = CoordinateIndex(synthetic_fai)

        start, end = index.get_scaffold_range("chr1")
        assert start == 0
        assert end == 1000

        start, end = index.get_scaffold_range("chr2")
        assert start == 1000
        assert end == 1500


# =============================================================================
# Schema Detection Tests
# =============================================================================


class TestSchemaDetection:
    """Tests for Helixer schema detection."""

    def test_detect_standard_schema(self, synthetic_hdf5: Path) -> None:
        """Test detection of standard predictions dataset."""
        with h5py.File(synthetic_hdf5, "r") as f:
            schema = detect_helixer_schema(f)

        assert schema["predictions_path"] == "predictions"
        assert schema["shape"] == (1500, 4)
        assert schema["n_classes"] == 4

    def test_detect_alt_name_schema(self, synthetic_hdf5_alt_name: Path) -> None:
        """Test detection of alternative dataset name (y_pred)."""
        with h5py.File(synthetic_hdf5_alt_name, "r") as f:
            schema = detect_helixer_schema(f)

        assert schema["predictions_path"] == "y_pred"
        assert schema["n_classes"] == 4

    def test_detect_3d_schema(self, synthetic_hdf5_3d: Path) -> None:
        """Test detection of 3D batched format."""
        with h5py.File(synthetic_hdf5_3d, "r") as f:
            schema = detect_helixer_schema(f)

        assert schema["predictions_path"] == "predictions"
        assert schema["shape"] == (3, 500, 4)
        assert schema["n_classes"] == 4

    def test_detect_invalid_schema(self, tmp_path: Path) -> None:
        """Test error for HDF5 without predictions."""
        h5_path = tmp_path / "invalid.h5"
        with h5py.File(h5_path, "w") as f:
            f.create_dataset("other_data", data=[1, 2, 3])

        with h5py.File(h5_path, "r") as f:
            with pytest.raises(HelixerFormatError, match="No predictions dataset"):
                detect_helixer_schema(f)

    def test_validate_helixer_file(self, synthetic_hdf5: Path) -> None:
        """Test the validate_helixer_file convenience function."""
        schema = validate_helixer_file(synthetic_hdf5)
        assert "predictions_path" in schema
        assert schema["n_classes"] == 4


# =============================================================================
# HelixerHDF5Reader Tests
# =============================================================================


class TestHelixerHDF5Reader:
    """Tests for the main HelixerHDF5Reader class."""

    def test_init_success(self, synthetic_hdf5: Path, synthetic_fai: Path) -> None:
        """Test successful initialization."""
        reader = HelixerHDF5Reader(synthetic_hdf5, synthetic_fai)

        assert reader.path == synthetic_hdf5
        assert reader.n_positions == 1500
        assert len(reader.scaffold_names) == 2
        assert reader.scaffold_lengths["chr1"] == 1000

        reader.close()

    def test_init_missing_hdf5(self, tmp_path: Path, synthetic_fai: Path) -> None:
        """Test error for missing HDF5 file."""
        with pytest.raises(FileNotFoundError, match="HDF5 file not found"):
            HelixerHDF5Reader(tmp_path / "missing.h5", synthetic_fai)

    def test_context_manager(self, synthetic_hdf5: Path, synthetic_fai: Path) -> None:
        """Test context manager usage."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            assert reader.n_positions == 1500

        # File should be closed after context exit
        assert reader._file is None

    def test_get_predictions_for_region(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test extracting predictions for a region."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            preds = reader.get_predictions_for_region("chr1", 100, 200)

            assert preds.shape == (100, 4)
            assert preds.dtype == np.float32
            # Probabilities should sum to ~1
            assert np.allclose(preds.sum(axis=1), 1.0, atol=1e-5)

    def test_get_predictions_chr2(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test extracting predictions from second scaffold."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            preds = reader.get_predictions_for_region("chr2", 0, 100)

            assert preds.shape == (100, 4)

    def test_get_predictions_invalid_seqid(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test error for invalid scaffold."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            with pytest.raises(KeyError, match="Unknown sequence"):
                reader.get_predictions_for_region("chrX", 0, 100)

    def test_get_predictions_invalid_coords(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test error for invalid coordinates."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            with pytest.raises(ValueError):
                reader.get_predictions_for_region("chr1", -10, 100)

            with pytest.raises(ValueError):
                reader.get_predictions_for_region("chr1", 0, 2000)

    def test_iter_chunks(self, synthetic_hdf5: Path, synthetic_fai: Path) -> None:
        """Test iterating in chunks."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai, chunk_size=300) as reader:
            chunks = list(reader.iter_chunks())

            # chr1 (1000) = 4 chunks of 300 + partial
            # chr2 (500) = 2 chunks
            assert len(chunks) >= 5

            # Check chunk structure
            seqid, start, end, preds = chunks[0]
            assert seqid == "chr1"
            assert start == 0
            assert end == 300
            assert preds.shape == (300, 4)

    def test_iter_chunks_single_scaffold(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test iterating chunks for a single scaffold."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai, chunk_size=300) as reader:
            chunks = list(reader.iter_chunks(seqid="chr2"))

            # All chunks should be from chr2
            for seqid, _, _, _ in chunks:
                assert seqid == "chr2"

            # chr2 is 500 bp, so 2 chunks at 300
            assert len(chunks) == 2

    def test_iter_scaffolds(self, synthetic_hdf5: Path, synthetic_fai: Path) -> None:
        """Test iterating complete scaffolds."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            scaffolds = list(reader.iter_scaffolds())

            assert len(scaffolds) == 2

            seqid1, preds1 = scaffolds[0]
            assert seqid1 == "chr1"
            assert preds1.shape == (1000, 4)

            seqid2, preds2 = scaffolds[1]
            assert seqid2 == "chr2"
            assert preds2.shape == (500, 4)

    def test_get_base_probabilities(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test getting probabilities for specific positions."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            positions = np.array([100, 200, 300, 400])
            probs = reader.get_base_probabilities("chr1", positions)

            assert probs.shape == (4, 4)
            assert probs.dtype == np.float32

    def test_get_base_probabilities_empty(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test with empty positions array."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            positions = np.array([], dtype=np.int64)
            probs = reader.get_base_probabilities("chr1", positions)

            assert probs.shape == (0, 4)

    def test_get_class_predictions(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test getting argmax class predictions."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            classes = reader.get_class_predictions("chr1", 0, 100)

            assert classes.shape == (100,)
            assert classes.dtype == np.int8
            assert np.all((classes >= 0) & (classes <= 3))

    def test_get_max_probabilities(
        self, synthetic_hdf5: Path, synthetic_fai: Path
    ) -> None:
        """Test getting maximum probabilities."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            max_probs = reader.get_max_probabilities("chr1", 0, 100)

            assert max_probs.shape == (100,)
            assert max_probs.dtype == np.float32
            assert np.all((max_probs >= 0.25) & (max_probs <= 1.0))

    def test_get_schema_info(self, synthetic_hdf5: Path, synthetic_fai: Path) -> None:
        """Test getting schema information."""
        with HelixerHDF5Reader(synthetic_hdf5, synthetic_fai) as reader:
            info = reader.get_schema_info()

            assert "predictions_path" in info
            assert "shape" in info
            assert info["scaffolds"] == 2
            assert info["total_bases"] == 1500


# =============================================================================
# 3D Array Tests
# =============================================================================


class TestHelixerHDF5Reader3D:
    """Tests for 3D batched HDF5 format."""

    @pytest.fixture
    def reader_3d(
        self, synthetic_hdf5_3d: Path, tmp_path: Path
    ) -> HelixerHDF5Reader:
        """Create reader with 3D HDF5 and matching FAI."""
        # Create FAI matching 3D structure (3 * 500 = 1500 total)
        fai_path = tmp_path / "genome_3d.fa.fai"
        with open(fai_path, "w") as f:
            f.write("chr1\t1000\t6\t80\t81\n")
            f.write("chr2\t500\t1025\t80\t81\n")

        reader = HelixerHDF5Reader(synthetic_hdf5_3d, fai_path)
        yield reader
        reader.close()

    def test_init_3d(self, reader_3d: HelixerHDF5Reader) -> None:
        """Test initialization with 3D array."""
        assert reader_3d.schema["shape"] == (3, 500, 4)

    def test_get_predictions_within_sample(
        self, reader_3d: HelixerHDF5Reader
    ) -> None:
        """Test region extraction within a single sample."""
        preds = reader_3d.get_predictions_for_region("chr1", 100, 200)
        assert preds.shape == (100, 4)

    def test_get_predictions_spanning_samples(
        self, reader_3d: HelixerHDF5Reader
    ) -> None:
        """Test region extraction spanning multiple samples."""
        # This spans from sample 0 to sample 1 (500 boundary)
        preds = reader_3d.get_predictions_for_region("chr1", 400, 600)
        assert preds.shape == (200, 4)


# =============================================================================
# Constants Tests
# =============================================================================


class TestConstants:
    """Tests for module constants."""

    def test_class_indices(self) -> None:
        """Test class index constants."""
        assert CLASS_INTERGENIC == 0
        assert CLASS_UTR == 1
        assert CLASS_CDS == 2
        assert CLASS_INTRON == 3

    def test_class_names(self) -> None:
        """Test class name list."""
        assert CLASS_NAMES == ["intergenic", "UTR", "CDS", "intron"]
        assert CLASS_NAMES[CLASS_CDS] == "CDS"
