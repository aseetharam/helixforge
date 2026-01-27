"""HDF5 file handling for Helixer predictions.

This module provides readers and writers for HDF5 files, specifically
the format used by Helixer for storing predictions.

Helixer HDF5 Format:
    The Helixer output contains softmax probabilities for each genomic
    position across four classes: intergenic, UTR, CDS, and intron.

    Helixer processes genomes in chunks (~21kb) and stores predictions
    for both strands separately:
    - Forward strand chunks: coordinates increase (start < end)
    - Reverse strand chunks: coordinates decrease (start > end)

    The companion input HDF5 file contains the chunk-to-coordinate mapping
    (seqids, start_ends) which is needed for proper strand-aware access.

Example:
    >>> from helixforge.io.hdf5 import HelixerHDF5Reader
    >>> # With strand-aware mapping (recommended)
    >>> reader = HelixerHDF5Reader(
    ...     "helixer_predictions.h5",
    ...     "genome.fa.fai",
    ...     input_h5_path="helixer_input.h5"
    ... )
    >>> # Get predictions for a gene on the minus strand
    >>> preds = reader.get_predictions_for_region("chr1", 1000, 2000, strand="-")
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
# Helixer Chunk Index (Strand-Aware)
# =============================================================================


class HelixerChunkIndex:
    """Strand-aware coordinate mapping for Helixer HDF5 predictions.

    Helixer stores predictions in chunks for both strands:
    - Forward strand: chunks where start < end (coordinates increase)
    - Reverse strand: chunks where start > end (coordinates decrease)

    This class loads the chunk mapping from the companion input HDF5 file
    and provides strand-aware coordinate lookups.

    Attributes:
        seqids: Array of sequence IDs for each chunk.
        start_ends: Array of (start, end) coordinates for each chunk.
        chunk_size: Size of each chunk (positions per chunk).
        is_forward: Boolean array indicating forward strand chunks.
    """

    def __init__(self, input_h5_path: Path | str) -> None:
        """Initialize chunk index from Helixer input HDF5.

        Args:
            input_h5_path: Path to the Helixer input HDF5 file containing
                the data/seqids and data/start_ends datasets.

        Raises:
            FileNotFoundError: If input HDF5 file doesn't exist.
            HelixerFormatError: If required datasets are missing.
        """
        self.input_path = Path(input_h5_path)
        if not self.input_path.exists():
            raise FileNotFoundError(f"Helixer input HDF5 not found: {self.input_path}")

        self._load_chunk_mapping()

    def _load_chunk_mapping(self) -> None:
        """Load chunk-to-coordinate mapping from input HDF5."""
        with h5py.File(self.input_path, "r") as f:
            # Check for required datasets
            if "data/seqids" not in f or "data/start_ends" not in f:
                raise HelixerFormatError(
                    f"Input HDF5 missing required datasets (data/seqids, data/start_ends). "
                    f"Available: {list(f.keys())}"
                )

            # Load mapping data
            self.seqids = f["data/seqids"][:].astype(str)
            self.start_ends = f["data/start_ends"][:]
            self.n_chunks = len(self.seqids)

            # Determine chunk size from first chunk
            self.chunk_size = abs(self.start_ends[0, 1] - self.start_ends[0, 0])

            # Determine strand for each chunk
            # Forward: start < end, Reverse: start > end
            self.is_forward = self.start_ends[:, 0] < self.start_ends[:, 1]

            # Build per-scaffold chunk indices for fast lookup
            self._build_scaffold_indices()

        logger.info(
            f"Loaded Helixer chunk index: {self.n_chunks} chunks, "
            f"chunk_size={self.chunk_size}, "
            f"{np.sum(self.is_forward)} forward, {np.sum(~self.is_forward)} reverse"
        )

    def _build_scaffold_indices(self) -> None:
        """Build indices for fast per-scaffold, per-strand lookup."""
        self.scaffold_chunks: dict[str, dict[str, list[int]]] = {}

        for chunk_idx in range(self.n_chunks):
            seqid = self.seqids[chunk_idx]
            strand = "+" if self.is_forward[chunk_idx] else "-"

            if seqid not in self.scaffold_chunks:
                self.scaffold_chunks[seqid] = {"+": [], "-": []}
            self.scaffold_chunks[seqid][strand].append(chunk_idx)

        # Sort chunks by start coordinate for binary search
        for seqid in self.scaffold_chunks:
            for strand in ["+", "-"]:
                chunks = self.scaffold_chunks[seqid][strand]
                if strand == "+":
                    # Forward: sort by start ascending
                    chunks.sort(key=lambda i: self.start_ends[i, 0])
                else:
                    # Reverse: sort by start descending (higher coords first)
                    chunks.sort(key=lambda i: self.start_ends[i, 0], reverse=True)

    def get_chunks_for_region(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str = "+",
    ) -> list[tuple[int, int, int]]:
        """Find chunks overlapping a genomic region on a specific strand.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).
            strand: Strand ("+" or "-").

        Returns:
            List of (chunk_idx, chunk_start_in_region, chunk_end_in_region) tuples.
            The chunk_start/end values are relative positions within the chunk
            that correspond to the requested region.

        Raises:
            KeyError: If seqid not found in index.
        """
        if seqid not in self.scaffold_chunks:
            raise KeyError(f"Unknown sequence: {seqid}")

        if strand not in ["+", "-"]:
            raise ValueError(f"Invalid strand: {strand}. Must be '+' or '-'.")

        chunk_indices = self.scaffold_chunks[seqid][strand]
        overlapping = []

        for chunk_idx in chunk_indices:
            chunk_start = self.start_ends[chunk_idx, 0]
            chunk_end = self.start_ends[chunk_idx, 1]

            if strand == "+":
                # Forward strand: normal coordinate comparison
                # Chunk covers [chunk_start, chunk_end)
                if chunk_end <= start:
                    continue  # Chunk entirely before region
                if chunk_start >= end:
                    break  # Chunk entirely after region (and sorted ascending)

                # Calculate overlap within chunk
                overlap_start = max(start, chunk_start)
                overlap_end = min(end, chunk_end)

                # Convert to positions within the chunk array
                pos_start = overlap_start - chunk_start
                pos_end = overlap_end - chunk_start

                overlapping.append((chunk_idx, pos_start, pos_end))
            else:
                # Reverse strand: chunk_start > chunk_end
                # Chunk covers [chunk_end, chunk_start) in genomic coords
                # But stored in reverse order in the array
                chunk_genomic_start = chunk_end  # Lower coordinate
                chunk_genomic_end = chunk_start  # Higher coordinate

                if chunk_genomic_end <= start:
                    break  # Chunk entirely before region (sorted descending)
                if chunk_genomic_start >= end:
                    continue  # Chunk entirely after region

                # Calculate overlap in genomic coordinates
                overlap_start = max(start, chunk_genomic_start)
                overlap_end = min(end, chunk_genomic_end)

                # Convert to positions within the chunk array
                # Reverse strand arrays go from high to low coords
                # Position 0 in array = chunk_start (high coord)
                # Position chunk_size-1 = chunk_end (low coord)
                chunk_length = chunk_genomic_end - chunk_genomic_start

                # overlap_start/end are in genomic coords (low to high)
                # Need to convert to array positions (high to low)
                pos_start = chunk_genomic_end - overlap_end
                pos_end = chunk_genomic_end - overlap_start

                overlapping.append((chunk_idx, pos_start, pos_end))

        return overlapping

    def get_scaffold_names(self) -> list[str]:
        """Get list of unique scaffold names."""
        return list(self.scaffold_chunks.keys())

    def has_dual_strand(self) -> bool:
        """Check if this index has both forward and reverse strand data."""
        n_forward = np.sum(self.is_forward)
        n_reverse = np.sum(~self.is_forward)
        return n_forward > 0 and n_reverse > 0

    def get_scaffold_lengths(self) -> dict[str, int]:
        """Get scaffold lengths derived from chunk coordinates.

        Returns:
            Dictionary mapping scaffold names to their lengths.
        """
        scaffold_lengths = {}
        for seqid in self.scaffold_chunks:
            # Get max coordinate from forward strand chunks
            fwd_chunks = self.scaffold_chunks[seqid]["+"]
            if fwd_chunks:
                max_coord = max(self.start_ends[i, 1] for i in fwd_chunks)
                scaffold_lengths[seqid] = max_coord
        return scaffold_lengths

    def to_coordinate_index(self) -> "CoordinateIndex":
        """Create a CoordinateIndex from this chunk index.

        This allows using the chunk index as a drop-in replacement
        for FAI-based coordinate mapping.

        Returns:
            CoordinateIndex with scaffold information from chunk data.
        """
        # Create a minimal CoordinateIndex without requiring FAI file
        coord_index = object.__new__(CoordinateIndex)
        coord_index.fai_path = None

        coord_index.scaffold_lengths = self.get_scaffold_lengths()
        coord_index.scaffold_order = list(self.scaffold_chunks.keys())

        # Build offsets
        offset = 0
        coord_index.scaffold_offsets = {}
        for seqid in coord_index.scaffold_order:
            coord_index.scaffold_offsets[seqid] = offset
            offset += coord_index.scaffold_lengths[seqid]

        coord_index.total_length = offset

        logger.debug(
            f"Created coordinate index from chunks: "
            f"{len(coord_index.scaffold_order)} scaffolds, "
            f"{coord_index.total_length:,} total bases"
        )

        return coord_index


def find_helixer_input_file(predictions_path: Path) -> Path | None:
    """Try to find the companion Helixer input HDF5 file.

    Helixer typically creates both an input and predictions file in the
    same directory with similar naming patterns.

    Args:
        predictions_path: Path to the predictions HDF5 file.

    Returns:
        Path to input HDF5 if found, None otherwise.
    """
    pred_dir = predictions_path.parent
    pred_name = predictions_path.stem

    # Common patterns for input files
    patterns = [
        pred_name.replace("_predictions", "_input") + ".h5",
        pred_name.replace("predictions", "input") + ".h5",
        pred_name + "_input.h5",
    ]

    for pattern in patterns:
        input_path = pred_dir / pattern
        if input_path.exists():
            logger.info(f"Found Helixer input file: {input_path}")
            return input_path

    # Also check if referenced in predictions file
    try:
        with h5py.File(predictions_path, "r") as f:
            if "test_data_path" in f.attrs:
                ref_path = f.attrs["test_data_path"]
                if isinstance(ref_path, bytes):
                    ref_path = ref_path.decode()
                # Try relative to predictions file
                input_path = pred_dir / Path(ref_path).name
                if input_path.exists():
                    logger.info(f"Found Helixer input file from attribute: {input_path}")
                    return input_path
    except Exception:
        pass

    return None


# =============================================================================
# Main Reader Class
# =============================================================================


class HelixerHDF5Reader:
    """Chunked, parallel-capable reader for Helixer HDF5 prediction files.

    Supports:
    - Strand-aware coordinate mapping (when input HDF5 is available)
    - Configurable chunk sizes (by bases or by scaffold)
    - Memory-mapped reading for large files
    - Thread-safe parallel access
    - Iterator interface for streaming processing

    Attributes:
        path: Path to HDF5 file.
        coord_index: Coordinate index for genomic mapping (FAI-based fallback).
        chunk_index: Strand-aware chunk index (when input HDF5 available).
        schema: Detected HDF5 schema information.

    Example:
        >>> # With strand-aware mapping (recommended)
        >>> reader = HelixerHDF5Reader(
        ...     "predictions.h5", "genome.fa.fai",
        ...     input_h5_path="input.h5"
        ... )
        >>> preds = reader.get_predictions_for_region("chr1", 1000, 2000, strand="-")
    """

    def __init__(
        self,
        h5_path: Path | str,
        fasta_index: Path | str | None = None,
        chunk_size: int = DEFAULT_CHUNK_SIZE,
        use_mmap: bool = False,
        n_workers: int = 1,
        input_h5_path: Path | str | None = None,
    ) -> None:
        """Initialize the Helixer HDF5 reader.

        Args:
            h5_path: Path to Helixer HDF5 predictions file.
            fasta_index: Path to .fai file for coordinate mapping. Optional if
                input_h5_path is provided or can be auto-detected.
            chunk_size: Number of bases per chunk for iteration.
            use_mmap: Use memory-mapped I/O (for very large files).
            n_workers: Number of workers for parallel access.
            input_h5_path: Path to Helixer input HDF5 file for strand-aware
                coordinate mapping. If not provided, will attempt to auto-detect.
                When provided, this is preferred over fasta_index.

        Raises:
            FileNotFoundError: If HDF5 file doesn't exist.
            HelixerFormatError: If HDF5 file is not valid Helixer format.
            ValueError: If neither fasta_index nor input_h5_path is available.
        """
        self.path = Path(h5_path)
        self.chunk_size = chunk_size
        self.use_mmap = use_mmap
        self.n_workers = n_workers

        if not self.path.exists():
            raise FileNotFoundError(f"HDF5 file not found: {self.path}")

        # Open HDF5 and detect schema
        self._lock = threading.Lock()
        self._file: h5py.File | None = None
        self._open_file()

        # Try to load strand-aware chunk index
        self.chunk_index: HelixerChunkIndex | None = None
        self._has_dual_strand = False

        # First try explicit input_h5_path
        if input_h5_path is not None:
            try:
                self.chunk_index = HelixerChunkIndex(input_h5_path)
                self._has_dual_strand = self.chunk_index.has_dual_strand()
            except (FileNotFoundError, HelixerFormatError) as e:
                logger.warning(f"Could not load chunk index from {input_h5_path}: {e}")

        # Then try auto-detect
        if self.chunk_index is None:
            input_path = find_helixer_input_file(self.path)
            if input_path is not None:
                try:
                    self.chunk_index = HelixerChunkIndex(input_path)
                    self._has_dual_strand = self.chunk_index.has_dual_strand()
                except (FileNotFoundError, HelixerFormatError) as e:
                    logger.warning(f"Could not load auto-detected chunk index: {e}")

        # Build coordinate index - prefer chunk index over FAI
        if self.chunk_index is not None:
            # Use chunk index to build coordinate index (no FAI needed)
            self.coord_index = self.chunk_index.to_coordinate_index()
        elif fasta_index is not None:
            # Fallback to FAI-based coordinate index
            self.coord_index = CoordinateIndex(fasta_index)
        else:
            raise ValueError(
                "No coordinate mapping available. Provide either fasta_index (.fai file) "
                "or input_h5_path (Helixer input HDF5), or ensure the input HDF5 file "
                "can be auto-detected alongside the predictions file."
            )

        # Validate array size vs genome
        predictions_shape = self.schema["shape"]
        if len(predictions_shape) == 2:
            array_length = predictions_shape[0]
        elif len(predictions_shape) == 3:
            array_length = predictions_shape[0] * predictions_shape[1]
        else:
            array_length = predictions_shape[0]

        genome_length = self.coord_index.total_length

        if array_length != genome_length:
            ratio = array_length / genome_length
            if self._has_dual_strand and 1.9 < ratio < 2.1:
                logger.info(
                    f"HDF5 contains dual-strand predictions "
                    f"({array_length:,} positions = ~2x genome length {genome_length:,}). "
                    f"Using strand-aware coordinate mapping."
                )
            else:
                logger.warning(
                    f"HDF5 array length ({array_length:,}) does not match "
                    f"genome length ({genome_length:,}). "
                    "Coordinate mapping may be incorrect."
                )

    def _open_file(self) -> None:
        """Open the HDF5 file and detect schema."""
        # Use core driver for memory-mapped access
        if self.use_mmap:
            self._file = h5py.File(
                self.path,
                "r",
                driver="core",
                backing_store=False,  # Don't write changes back
            )
        else:
            self._file = h5py.File(self.path, "r")
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
        strand: str = "+",
    ) -> NDArray[np.float32]:
        """Extract predictions for a genomic region.

        When strand-aware chunk index is available (dual-strand HDF5),
        retrieves predictions from the appropriate strand-specific chunks.
        For reverse strand, the predictions are returned in genomic order
        (low to high coordinates), matching the forward strand convention.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).
            strand: Strand ("+" or "-"). Used when dual-strand data available.

        Returns:
            Array of shape (end-start, 4) with class probabilities.
            For reverse strand, predictions are flipped to genomic order.

        Raises:
            KeyError: If seqid not in genome.
            ValueError: If coordinates are invalid.
        """
        # Use strand-aware chunk index if available
        if self.chunk_index is not None and self._has_dual_strand:
            return self._get_predictions_stranded(seqid, start, end, strand)

        # Fallback to original FAI-based method (no strand awareness)
        return self._get_predictions_legacy(seqid, start, end)

    def _get_predictions_stranded(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str,
    ) -> NDArray[np.float32]:
        """Get predictions using strand-aware chunk mapping.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).
            strand: Strand ("+" or "-").

        Returns:
            Array of shape (end-start, 4) with class probabilities.
        """
        assert self.chunk_index is not None

        # Get overlapping chunks for this region and strand
        chunk_info = self.chunk_index.get_chunks_for_region(seqid, start, end, strand)

        if not chunk_info:
            raise ValueError(
                f"No chunks found for region {seqid}:{start}-{end} on strand {strand}"
            )

        with self._lock:
            dataset = self.predictions_dataset
            chunks_data = []

            for chunk_idx, pos_start, pos_end in chunk_info:
                chunk_preds = dataset[chunk_idx, pos_start:pos_end, :]
                chunks_data.append(chunk_preds)

            predictions = np.concatenate(chunks_data, axis=0)

        # For reverse strand, the predictions are stored in reverse order
        # (high coords to low coords). Flip to match genomic order.
        if strand == "-":
            predictions = np.flip(predictions, axis=0)

        return predictions.astype(np.float32)

    def _get_predictions_legacy(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> NDArray[np.float32]:
        """Get predictions using legacy FAI-based coordinate mapping.

        This method doesn't account for strand and assumes predictions
        are concatenated linearly. Use when chunk index is not available.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (0-based, exclusive).

        Returns:
            Array of shape (end-start, 4) with class probabilities.
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

    @property
    def has_dual_strand(self) -> bool:
        """Check if this reader has dual-strand prediction data."""
        return self._has_dual_strand

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
