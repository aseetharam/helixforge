"""HDF5 file handling for Helixer predictions.

This module provides readers and writers for HDF5 files, specifically
the format used by Helixer for storing predictions.

Helixer HDF5 Format:
    The Helixer output contains softmax probabilities for each genomic
    position across four classes: intergenic, UTR, CDS, and intron.
    Sequences are concatenated; coordinate mapping is required via
    genome FASTA index.

Example:
    >>> from helixforge.io.hdf5 import HelixerHDF5Reader
    >>> reader = HelixerHDF5Reader("helixer_output.h5", "genome.fa.fai")
    >>> for seqid, start, end, preds in reader.iter_chunks():
    ...     process(seqid, start, end, preds)
"""

from __future__ import annotations

import logging
import threading
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator

import h5py
import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

# =============================================================================
# Constants
# =============================================================================

# Helixer class indices
CLASS_INTERGENIC = 0
CLASS_UTR = 1
CLASS_CDS = 2
CLASS_INTRON = 3

# Class names for reference
CLASS_NAMES = ["intergenic", "UTR", "CDS", "intron"]

# Expected dataset names in Helixer output (may vary by version)
DATASET_PREDICTIONS = "predictions"
DATASET_SEQUENCE_IDS = "seqids"
DATASET_START_ENDS = "start_ends"

# Alternative dataset names for different Helixer versions
ALT_DATASET_NAMES = {
    "predictions": ["predictions", "y_pred", "data/predictions"],
    "seqids": ["seqids", "sequence_ids", "data/seqids"],
}

# Default chunk size for iteration (50 Mb worth of predictions)
DEFAULT_CHUNK_SIZE = 50_000_000

logger = logging.getLogger(__name__)


# =============================================================================
# Exceptions
# =============================================================================


class HelixerFormatError(ValueError):
    """Raised when HDF5 file doesn't conform to Helixer format."""

    pass


# =============================================================================
# Schema Detection
# =============================================================================


def detect_helixer_schema(h5file: h5py.File) -> dict[str, Any]:
    """Detect the schema of a Helixer HDF5 file.

    Args:
        h5file: Open HDF5 file handle.

    Returns:
        Dictionary with schema information including:
        - predictions_path: Path to predictions dataset
        - seqids_path: Path to sequence IDs dataset (if present)
        - shape: Shape of predictions array
        - n_classes: Number of prediction classes
        - has_seqids: Whether sequence IDs are embedded

    Raises:
        HelixerFormatError: If no valid predictions dataset found.
    """
    schema: dict[str, Any] = {
        "predictions_path": None,
        "seqids_path": None,
        "shape": None,
        "n_classes": None,
        "has_seqids": False,
        "version": "unknown",
    }

    # Find predictions dataset
    for name in ALT_DATASET_NAMES["predictions"]:
        if name in h5file:
            schema["predictions_path"] = name
            dataset = h5file[name]
            schema["shape"] = dataset.shape
            schema["n_classes"] = dataset.shape[-1] if len(dataset.shape) > 1 else 4
            schema["dtype"] = str(dataset.dtype)
            break

    if schema["predictions_path"] is None:
        # Try to find any dataset that looks like predictions
        for key in h5file.keys():
            dataset = h5file[key]
            if isinstance(dataset, h5py.Dataset):
                shape = dataset.shape
                if len(shape) >= 2 and shape[-1] == 4:
                    schema["predictions_path"] = key
                    schema["shape"] = shape
                    schema["n_classes"] = 4
                    schema["dtype"] = str(dataset.dtype)
                    logger.warning(f"Using non-standard predictions dataset: {key}")
                    break

    if schema["predictions_path"] is None:
        raise HelixerFormatError(
            "No predictions dataset found in HDF5 file. "
            f"Available datasets: {list(h5file.keys())}"
        )

    # Find sequence IDs dataset
    for name in ALT_DATASET_NAMES["seqids"]:
        if name in h5file:
            schema["seqids_path"] = name
            schema["has_seqids"] = True
            break

    return schema


# =============================================================================
# Coordinate Index
# =============================================================================


class CoordinateIndex:
    """Maps genomic coordinates to HDF5 array indices.

    Built from FASTA index (.fai) file to enable random access to
    genomic regions within the concatenated HDF5 predictions.

    Attributes:
        scaffold_lengths: Dict mapping scaffold names to lengths.
        scaffold_offsets: Dict mapping scaffold names to start offsets.
        scaffold_order: List of scaffolds in order.
        total_length: Total genome length.
    """

    def __init__(self, fai_path: Path | str) -> None:
        """Initialize coordinate index from FAI file.

        Args:
            fai_path: Path to .fai index file.

        Raises:
            FileNotFoundError: If FAI file doesn't exist.
            ValueError: If FAI file is malformed.
        """
        self.fai_path = Path(fai_path)
        if not self.fai_path.exists():
            raise FileNotFoundError(f"FAI index not found: {self.fai_path}")

        self.scaffold_lengths: dict[str, int] = {}
        self.scaffold_offsets: dict[str, int] = {}
        self.scaffold_order: list[str] = []

        self._parse_fai()

    def _parse_fai(self) -> None:
        """Parse FAI file and build coordinate index."""
        offset = 0
        with open(self.fai_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                seqid = parts[0]
                length = int(parts[1])

                self.scaffold_lengths[seqid] = length
                self.scaffold_offsets[seqid] = offset
                self.scaffold_order.append(seqid)
                offset += length

        self.total_length = offset
        logger.debug(
            f"Loaded FAI index: {len(self.scaffold_order)} scaffolds, "
            f"{self.total_length:,} total bases"
        )

    def genomic_to_array(self, seqid: str, start: int, end: int) -> tuple[int, int]:
        """Convert genomic coordinates to array indices.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).

        Returns:
            Tuple of (array_start, array_end) indices.

        Raises:
            KeyError: If seqid not in index.
            ValueError: If coordinates are out of bounds.
        """
        if seqid not in self.scaffold_offsets:
            raise KeyError(f"Unknown sequence: {seqid}")

        scaffold_length = self.scaffold_lengths[seqid]
        if start < 0 or end > scaffold_length or start >= end:
            raise ValueError(
                f"Invalid coordinates for {seqid}: {start}-{end} "
                f"(scaffold length: {scaffold_length})"
            )

        offset = self.scaffold_offsets[seqid]
        return offset + start, offset + end

    def array_to_genomic(self, array_idx: int) -> tuple[str, int]:
        """Convert array index to genomic coordinate.

        Args:
            array_idx: Index into concatenated array.

        Returns:
            Tuple of (seqid, position).

        Raises:
            ValueError: If index is out of bounds.
        """
        if array_idx < 0 or array_idx >= self.total_length:
            raise ValueError(f"Array index out of bounds: {array_idx}")

        for seqid in self.scaffold_order:
            offset = self.scaffold_offsets[seqid]
            length = self.scaffold_lengths[seqid]
            if array_idx < offset + length:
                return seqid, array_idx - offset

        raise ValueError(f"Could not map array index: {array_idx}")

    def get_scaffold_range(self, seqid: str) -> tuple[int, int]:
        """Get array index range for a scaffold.

        Args:
            seqid: Scaffold name.

        Returns:
            Tuple of (start_idx, end_idx).
        """
        offset = self.scaffold_offsets[seqid]
        length = self.scaffold_lengths[seqid]
        return offset, offset + length


# =============================================================================
# Main Reader Class
# =============================================================================


class HelixerHDF5Reader:
    """Chunked, parallel-capable reader for Helixer HDF5 prediction files.

    Supports:
    - Configurable chunk sizes (by bases or by scaffold)
    - Memory-mapped reading for large files
    - Thread-safe parallel access
    - Iterator interface for streaming processing

    Attributes:
        path: Path to HDF5 file.
        coord_index: Coordinate index for genomic mapping.
        schema: Detected HDF5 schema information.

    Example:
        >>> reader = HelixerHDF5Reader("predictions.h5", "genome.fa.fai")
        >>> for seqid, start, end, preds in reader.iter_chunks():
        ...     process_chunk(seqid, start, end, preds)
    """

    def __init__(
        self,
        h5_path: Path | str,
        fasta_index: Path | str,
        chunk_size: int = DEFAULT_CHUNK_SIZE,
        use_mmap: bool = False,
        n_workers: int = 1,
    ) -> None:
        """Initialize the Helixer HDF5 reader.

        Args:
            h5_path: Path to Helixer HDF5 predictions file.
            fasta_index: Path to .fai file for coordinate mapping.
            chunk_size: Number of bases per chunk for iteration.
            use_mmap: Use memory-mapped I/O (for very large files).
            n_workers: Number of workers for parallel access.

        Raises:
            FileNotFoundError: If HDF5 or FAI file doesn't exist.
            HelixerFormatError: If HDF5 file is not valid Helixer format.
        """
        self.path = Path(h5_path)
        self.chunk_size = chunk_size
        self.use_mmap = use_mmap
        self.n_workers = n_workers

        if not self.path.exists():
            raise FileNotFoundError(f"HDF5 file not found: {self.path}")

        # Build coordinate index from FAI
        self.coord_index = CoordinateIndex(fasta_index)

        # Open HDF5 and detect schema
        self._lock = threading.Lock()
        self._file: h5py.File | None = None
        self._open_file()

        # Validate array size matches genome
        predictions_shape = self.schema["shape"]
        # Handle both 2D (positions, classes) and 3D (samples, positions, classes) shapes
        if len(predictions_shape) == 2:
            array_length = predictions_shape[0]
        elif len(predictions_shape) == 3:
            array_length = predictions_shape[0] * predictions_shape[1]
        else:
            array_length = predictions_shape[0]

        if array_length != self.coord_index.total_length:
            logger.warning(
                f"HDF5 array length ({array_length:,}) does not match "
                f"genome length ({self.coord_index.total_length:,}). "
                "Coordinate mapping may be incorrect."
            )

    def _open_file(self) -> None:
        """Open the HDF5 file and detect schema."""
        driver = "core" if self.use_mmap else None
        backing_store = not self.use_mmap

        self._file = h5py.File(
            self.path,
            "r",
            driver=driver,
            backing_store=backing_store if driver == "core" else True,
        )
        self.schema = detect_helixer_schema(self._file)
        logger.info(
            f"Opened HDF5 file: {self.path.name}, "
            f"shape={self.schema['shape']}, "
            f"dtype={self.schema.get('dtype', 'unknown')}"
        )

    @property
    def predictions_dataset(self) -> h5py.Dataset:
        """Get the predictions dataset."""
        if self._file is None:
            raise RuntimeError("HDF5 file is not open")
        return self._file[self.schema["predictions_path"]]

    @property
    def n_positions(self) -> int:
        """Total number of genomic positions."""
        return self.coord_index.total_length

    @property
    def scaffold_names(self) -> list[str]:
        """List of scaffold names in order."""
        return self.coord_index.scaffold_order.copy()

    @property
    def scaffold_lengths(self) -> dict[str, int]:
        """Dictionary of scaffold lengths."""
        return self.coord_index.scaffold_lengths.copy()

    def __enter__(self) -> HelixerHDF5Reader:
        """Context manager entry."""
        return self

    def __exit__(self, *args: Any) -> None:
        """Context manager exit."""
        self.close()

    def close(self) -> None:
        """Close the HDF5 file."""
        with self._lock:
            if self._file is not None:
                self._file.close()
                self._file = None

    def get_schema_info(self) -> dict[str, Any]:
        """Get schema information about the HDF5 file.

        Returns:
            Dictionary with schema details.
        """
        info = self.schema.copy()
        info["scaffolds"] = len(self.coord_index.scaffold_order)
        info["total_bases"] = self.coord_index.total_length
        return info

    def get_predictions_for_region(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> NDArray[np.float32]:
        """Extract predictions for a genomic region.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).

        Returns:
            Array of shape (end-start, 4) with class probabilities.

        Raises:
            KeyError: If seqid not in genome.
            ValueError: If coordinates are invalid.
        """
        array_start, array_end = self.coord_index.genomic_to_array(seqid, start, end)

        with self._lock:
            dataset = self.predictions_dataset
            shape = dataset.shape

            # Handle different array shapes
            if len(shape) == 2:
                # Shape: (positions, classes)
                predictions = dataset[array_start:array_end, :]
            elif len(shape) == 3:
                # Shape: (samples, positions, classes) - flatten first two dims
                # This handles batched predictions
                flat_idx_start = array_start
                flat_idx_end = array_end
                sample_size = shape[1]

                sample_start = flat_idx_start // sample_size
                pos_start = flat_idx_start % sample_size
                sample_end = (flat_idx_end - 1) // sample_size
                pos_end = (flat_idx_end - 1) % sample_size + 1

                if sample_start == sample_end:
                    predictions = dataset[sample_start, pos_start:pos_end, :]
                else:
                    # Spanning multiple samples - need to concatenate
                    chunks = []
                    for s in range(sample_start, sample_end + 1):
                        if s == sample_start:
                            chunks.append(dataset[s, pos_start:, :])
                        elif s == sample_end:
                            chunks.append(dataset[s, :pos_end, :])
                        else:
                            chunks.append(dataset[s, :, :])
                    predictions = np.concatenate(chunks, axis=0)
            else:
                raise HelixerFormatError(
                    f"Unexpected predictions shape: {shape}. "
                    "Expected 2D (positions, classes) or 3D (samples, positions, classes)."
                )

        return predictions.astype(np.float32)

    def iter_chunks(
        self,
        seqid: str | None = None,
        chunk_size: int | None = None,
    ) -> Iterator[tuple[str, int, int, NDArray[np.float32]]]:
        """Iterate over genome in chunks.

        Yields predictions in fixed-size chunks for memory-efficient
        processing of large genomes.

        Args:
            seqid: If provided, only yield chunks from this scaffold.
            chunk_size: Override default chunk size.

        Yields:
            Tuples of (seqid, start, end, predictions_array).
        """
        chunk_size = chunk_size or self.chunk_size

        scaffolds = [seqid] if seqid else self.coord_index.scaffold_order

        for scaffold in scaffolds:
            scaffold_length = self.coord_index.scaffold_lengths[scaffold]

            for start in range(0, scaffold_length, chunk_size):
                end = min(start + chunk_size, scaffold_length)
                predictions = self.get_predictions_for_region(scaffold, start, end)
                yield scaffold, start, end, predictions

    def iter_scaffolds(self) -> Iterator[tuple[str, NDArray[np.float32]]]:
        """Iterate over complete scaffolds.

        Yields:
            Tuples of (seqid, predictions_array).

        Warning:
            May use significant memory for large scaffolds.
        """
        for seqid in self.coord_index.scaffold_order:
            length = self.coord_index.scaffold_lengths[seqid]
            predictions = self.get_predictions_for_region(seqid, 0, length)
            yield seqid, predictions

    def get_base_probabilities(
        self,
        seqid: str,
        positions: NDArray[np.int64],
    ) -> NDArray[np.float32]:
        """Get probabilities for specific positions.

        Useful for extracting predictions at gene boundaries or
        specific features for scoring.

        Args:
            seqid: Scaffold name.
            positions: Array of 0-based positions.

        Returns:
            Array of shape (len(positions), 4) with class probabilities.

        Raises:
            ValueError: If any position is out of bounds.
        """
        if len(positions) == 0:
            return np.empty((0, 4), dtype=np.float32)

        # Validate positions
        scaffold_length = self.coord_index.scaffold_lengths.get(seqid)
        if scaffold_length is None:
            raise KeyError(f"Unknown sequence: {seqid}")

        if np.any(positions < 0) or np.any(positions >= scaffold_length):
            raise ValueError(
                f"Positions out of bounds for {seqid} (length={scaffold_length})"
            )

        # Convert to array indices
        offset = self.coord_index.scaffold_offsets[seqid]
        array_indices = positions + offset

        with self._lock:
            dataset = self.predictions_dataset
            shape = dataset.shape

            if len(shape) == 2:
                # Direct indexing for 2D array
                predictions = dataset[array_indices, :]
            else:
                # For 3D, need to compute sample/position indices
                sample_size = shape[1]
                results = []
                for idx in array_indices:
                    s = idx // sample_size
                    p = idx % sample_size
                    results.append(dataset[s, p, :])
                predictions = np.array(results)

        return predictions.astype(np.float32)

    def get_class_predictions(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> NDArray[np.int8]:
        """Get argmax class predictions for a region.

        Args:
            seqid: Scaffold name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            Array of shape (end-start,) with class indices (0-3).
        """
        predictions = self.get_predictions_for_region(seqid, start, end)
        return np.argmax(predictions, axis=1).astype(np.int8)

    def get_max_probabilities(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> NDArray[np.float32]:
        """Get maximum class probability for each position.

        Args:
            seqid: Scaffold name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            Array of shape (end-start,) with max probabilities.
        """
        predictions = self.get_predictions_for_region(seqid, start, end)
        return np.max(predictions, axis=1)


# =============================================================================
# Convenience Functions
# =============================================================================


def load_helixer_predictions(
    h5_path: Path | str,
    fasta_index: Path | str,
) -> HelixerHDF5Reader:
    """Convenience function to load Helixer predictions.

    Args:
        h5_path: Path to HDF5 file.
        fasta_index: Path to FAI file.

    Returns:
        HelixerHDF5Reader instance.
    """
    return HelixerHDF5Reader(h5_path, fasta_index)


def validate_helixer_file(h5_path: Path | str) -> dict[str, Any]:
    """Validate a Helixer HDF5 file without coordinate index.

    Args:
        h5_path: Path to HDF5 file.

    Returns:
        Schema information dictionary.

    Raises:
        HelixerFormatError: If file is not valid Helixer format.
    """
    with h5py.File(h5_path, "r") as f:
        return detect_helixer_schema(f)


# =============================================================================
# Legacy Compatibility
# =============================================================================

# Alias for backward compatibility with stub
HelixerReader = HelixerHDF5Reader
