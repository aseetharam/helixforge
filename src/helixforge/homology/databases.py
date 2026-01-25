"""Database management for homology searches.

This module provides utilities for downloading, managing, and formatting
protein reference databases for homology validation.

Supported databases:
- Swiss-Prot (curated): High-quality, manually annotated proteins
- TrEMBL (automatic): Computationally annotated proteins
- UniRef90/50/100: Clustered protein sequences
- OrthoDB: Orthologous groups

Example:
    >>> from helixforge.homology.databases import DatabaseManager
    >>> db_manager = DatabaseManager(cache_dir="~/.helixforge/databases")
    >>> swissprot_path = db_manager.get_database("swissprot_plants")
    >>> # Database will be downloaded if not present
"""

from __future__ import annotations

import gzip
import hashlib
import logging
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Iterator
from urllib.request import urlopen, urlretrieve
from urllib.error import URLError, HTTPError

logger = logging.getLogger(__name__)


# =============================================================================
# Constants
# =============================================================================


# UniProt download URLs
UNIPROT_BASE = "https://ftp.uniprot.org/pub/databases/uniprot"
SWISSPROT_URL = f"{UNIPROT_BASE}/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
TREMBL_URL = f"{UNIPROT_BASE}/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"

# UniProt taxonomy-specific URLs
SWISSPROT_TAXONOMY_URL = (
    f"{UNIPROT_BASE}/current_release/knowledgebase/taxonomic_divisions/"
    "uniprot_sprot_{division}.dat.gz"
)

# UniRef URLs
UNIREF_URL = f"{UNIPROT_BASE}/uniref/uniref{{level}}/uniref{{level}}.fasta.gz"

# Taxonomy division codes
TAXONOMY_DIVISIONS = {
    "archaea": "archaea",
    "bacteria": "bacteria",
    "fungi": "fungi",
    "human": "human",
    "invertebrates": "invertebrates",
    "mammals": "mammals",
    "plants": "plants",
    "rodents": "rodents",
    "vertebrates": "vertebrates",
    "viruses": "viruses",
}

# Plant taxonomy IDs for filtering
PLANT_TAXONOMY_IDS = [
    3193,   # Embryophyta (land plants)
    33090,  # Viridiplantae (green plants)
    35493,  # Streptophyta
]

# Default cache directory
DEFAULT_CACHE_DIR = Path.home() / ".helixforge" / "databases"


# =============================================================================
# Enums
# =============================================================================


class DatabaseType(Enum):
    """Types of protein databases."""

    SWISSPROT = "swissprot"
    """Swiss-Prot: curated, high-quality annotations."""

    TREMBL = "trembl"
    """TrEMBL: automated annotations."""

    UNIREF100 = "uniref100"
    """UniRef100: clustered at 100% identity."""

    UNIREF90 = "uniref90"
    """UniRef90: clustered at 90% identity."""

    UNIREF50 = "uniref50"
    """UniRef50: clustered at 50% identity."""

    CUSTOM = "custom"
    """Custom user-provided database."""


# =============================================================================
# Data Structures
# =============================================================================


@dataclass
class DatabaseInfo:
    """Information about a protein database.

    Attributes:
        name: Human-readable database name.
        db_type: Database type enum.
        path: Local path to database file.
        formatted_path: Path to formatted (Diamond/MMseqs2) database.
        download_url: URL for downloading.
        download_date: When the database was downloaded.
        n_sequences: Number of sequences (if known).
        taxonomy_filter: Taxonomy filter applied (if any).
        checksum: MD5 checksum of downloaded file.
    """

    name: str
    db_type: DatabaseType
    path: Path | None = None
    formatted_path: Path | None = None
    download_url: str | None = None
    download_date: datetime | None = None
    n_sequences: int | None = None
    taxonomy_filter: str | None = None
    checksum: str | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "name": self.name,
            "db_type": self.db_type.value,
            "path": str(self.path) if self.path else None,
            "formatted_path": str(self.formatted_path) if self.formatted_path else None,
            "download_url": self.download_url,
            "download_date": self.download_date.isoformat() if self.download_date else None,
            "n_sequences": self.n_sequences,
            "taxonomy_filter": self.taxonomy_filter,
            "checksum": self.checksum,
        }

    @property
    def is_downloaded(self) -> bool:
        """Check if database is downloaded."""
        return self.path is not None and self.path.exists()

    @property
    def is_formatted(self) -> bool:
        """Check if database is formatted for searching."""
        if self.formatted_path is None:
            return False
        # Check for Diamond .dmnd extension
        dmnd_path = self.formatted_path.with_suffix(".dmnd")
        if dmnd_path.exists():
            return True
        # Check for raw formatted path
        return self.formatted_path.exists()


# =============================================================================
# Download Functions
# =============================================================================


def download_file(
    url: str,
    output_path: Path,
    show_progress: bool = True,
    chunk_size: int = 8192,
) -> Path:
    """Download a file from URL.

    Args:
        url: URL to download from.
        output_path: Local path to save file.
        show_progress: Whether to log progress.
        chunk_size: Download chunk size in bytes.

    Returns:
        Path to downloaded file.

    Raises:
        URLError: If download fails.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading {url}")
    logger.info(f"  -> {output_path}")

    try:
        with urlopen(url, timeout=30) as response:
            # Get file size if available
            total_size = response.headers.get("Content-Length")
            if total_size:
                total_size = int(total_size)
                logger.info(f"  Size: {total_size / 1024 / 1024:.1f} MB")

            downloaded = 0
            last_progress = 0

            with open(output_path, "wb") as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)

                    # Log progress every 10%
                    if show_progress and total_size:
                        progress = int(100 * downloaded / total_size)
                        if progress >= last_progress + 10:
                            logger.info(f"  Progress: {progress}%")
                            last_progress = progress

        logger.info(f"  Downloaded: {downloaded / 1024 / 1024:.1f} MB")
        return output_path

    except (URLError, HTTPError) as e:
        if output_path.exists():
            output_path.unlink()
        raise URLError(f"Failed to download {url}: {e}") from e


def calculate_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """Calculate MD5 checksum of a file.

    Args:
        file_path: Path to file.
        chunk_size: Read chunk size.

    Returns:
        MD5 hex digest.
    """
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def decompress_gzip(
    gz_path: Path,
    output_path: Path | None = None,
    keep_original: bool = False,
) -> Path:
    """Decompress a gzip file.

    Args:
        gz_path: Path to gzip file.
        output_path: Output path. Auto-determined if None.
        keep_original: Keep the original gzip file.

    Returns:
        Path to decompressed file.
    """
    gz_path = Path(gz_path)

    if output_path is None:
        # Remove .gz extension
        output_path = gz_path.with_suffix("")
        if output_path.suffix == "":
            output_path = gz_path.parent / (gz_path.stem + ".fasta")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Decompressing {gz_path.name}")

    with gzip.open(gz_path, "rb") as f_in:
        with open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    if not keep_original:
        gz_path.unlink()
        logger.debug(f"Removed {gz_path}")

    logger.info(f"Decompressed to {output_path}")
    return output_path


# =============================================================================
# Swiss-Prot Functions
# =============================================================================


def download_swissprot(
    output_dir: Path | str,
    taxonomy: str | None = None,
    decompress: bool = True,
) -> DatabaseInfo:
    """Download Swiss-Prot database.

    Args:
        output_dir: Directory to save database.
        taxonomy: Taxonomy filter (e.g., "plants", "fungi").
        decompress: Decompress the downloaded file.

    Returns:
        DatabaseInfo with paths and metadata.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine URL and filename
    if taxonomy and taxonomy.lower() in TAXONOMY_DIVISIONS:
        # Use taxonomy-specific file
        division = TAXONOMY_DIVISIONS[taxonomy.lower()]
        url = SWISSPROT_TAXONOMY_URL.format(division=division)
        filename = f"uniprot_sprot_{division}.fasta.gz"
        name = f"Swiss-Prot ({taxonomy})"
    else:
        # Use full Swiss-Prot
        url = SWISSPROT_URL
        filename = "uniprot_sprot.fasta.gz"
        name = "Swiss-Prot (complete)"

    gz_path = output_dir / filename

    # Download
    download_file(url, gz_path)

    # Calculate checksum
    checksum = calculate_md5(gz_path)

    # Decompress if requested
    if decompress:
        fasta_path = decompress_gzip(gz_path)
    else:
        fasta_path = gz_path

    # Count sequences
    n_sequences = count_fasta_sequences(fasta_path)

    return DatabaseInfo(
        name=name,
        db_type=DatabaseType.SWISSPROT,
        path=fasta_path,
        download_url=url,
        download_date=datetime.now(),
        n_sequences=n_sequences,
        taxonomy_filter=taxonomy,
        checksum=checksum,
    )


def download_swissprot_plants(
    output_dir: Path | str,
    decompress: bool = True,
) -> DatabaseInfo:
    """Download Swiss-Prot plant proteins.

    Convenience function for plant genomics workflows.

    Args:
        output_dir: Directory to save database.
        decompress: Decompress the downloaded file.

    Returns:
        DatabaseInfo with paths and metadata.
    """
    return download_swissprot(output_dir, taxonomy="plants", decompress=decompress)


# =============================================================================
# UniRef Functions
# =============================================================================


def download_uniref(
    output_dir: Path | str,
    level: int = 90,
    decompress: bool = True,
) -> DatabaseInfo:
    """Download UniRef database.

    Args:
        output_dir: Directory to save database.
        level: Clustering level (100, 90, or 50).
        decompress: Decompress the downloaded file.

    Returns:
        DatabaseInfo with paths and metadata.

    Raises:
        ValueError: If level is not 100, 90, or 50.
    """
    if level not in (100, 90, 50):
        raise ValueError(f"UniRef level must be 100, 90, or 50, got {level}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    url = UNIREF_URL.format(level=level)
    filename = f"uniref{level}.fasta.gz"
    gz_path = output_dir / filename

    # Download
    download_file(url, gz_path)

    # Calculate checksum
    checksum = calculate_md5(gz_path)

    # Decompress if requested
    if decompress:
        fasta_path = decompress_gzip(gz_path)
    else:
        fasta_path = gz_path

    # Determine database type
    db_type = {
        100: DatabaseType.UNIREF100,
        90: DatabaseType.UNIREF90,
        50: DatabaseType.UNIREF50,
    }[level]

    return DatabaseInfo(
        name=f"UniRef{level}",
        db_type=db_type,
        path=fasta_path,
        download_url=url,
        download_date=datetime.now(),
        taxonomy_filter=None,
        checksum=checksum,
    )


# =============================================================================
# FASTA Utilities
# =============================================================================


def count_fasta_sequences(fasta_path: Path | str) -> int:
    """Count sequences in a FASTA file.

    Args:
        fasta_path: Path to FASTA file.

    Returns:
        Number of sequences.
    """
    fasta_path = Path(fasta_path)
    count = 0

    # Handle gzipped files
    if fasta_path.suffix == ".gz":
        opener = gzip.open
        mode = "rt"
    else:
        opener = open
        mode = "r"

    with opener(fasta_path, mode) as f:
        for line in f:
            if line.startswith(">"):
                count += 1

    return count


def iter_fasta(fasta_path: Path | str) -> Iterator[tuple[str, str, str]]:
    """Iterate over FASTA sequences.

    Args:
        fasta_path: Path to FASTA file.

    Yields:
        (sequence_id, description, sequence) tuples.
    """
    fasta_path = Path(fasta_path)

    # Handle gzipped files
    if fasta_path.suffix == ".gz":
        opener = gzip.open
        mode = "rt"
    else:
        opener = open
        mode = "r"

    with opener(fasta_path, mode) as f:
        seq_id = None
        description = ""
        sequence_parts = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Yield previous sequence
                if seq_id is not None:
                    yield seq_id, description, "".join(sequence_parts)

                # Parse new header
                header = line[1:]
                parts = header.split(None, 1)
                seq_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                sequence_parts = []
            else:
                sequence_parts.append(line)

        # Yield last sequence
        if seq_id is not None:
            yield seq_id, description, "".join(sequence_parts)


def filter_fasta_by_taxonomy(
    input_fasta: Path | str,
    output_fasta: Path | str,
    taxonomy_ids: list[int],
) -> int:
    """Filter FASTA by taxonomy IDs in UniProt headers.

    UniProt FASTA headers include taxonomy info like:
    >sp|P12345|NAME_ARATH ... OX=3702 ...

    Args:
        input_fasta: Input FASTA path.
        output_fasta: Output FASTA path.
        taxonomy_ids: List of NCBI taxonomy IDs to keep.

    Returns:
        Number of sequences written.
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    taxonomy_set = set(taxonomy_ids)
    n_written = 0

    with open(output_fasta, "w") as out:
        for seq_id, description, sequence in iter_fasta(input_fasta):
            # Parse taxonomy ID from description
            # Format: ... OX=12345 ...
            ox_match = None
            for part in description.split():
                if part.startswith("OX="):
                    try:
                        ox_match = int(part[3:])
                    except ValueError:
                        pass
                    break

            if ox_match and ox_match in taxonomy_set:
                out.write(f">{seq_id} {description}\n")
                # Write sequence in 60-char lines
                for i in range(0, len(sequence), 60):
                    out.write(sequence[i:i + 60] + "\n")
                n_written += 1

    logger.info(f"Filtered {n_written} sequences by taxonomy")
    return n_written


def subset_fasta(
    input_fasta: Path | str,
    output_fasta: Path | str,
    sequence_ids: set[str],
) -> int:
    """Subset FASTA file by sequence IDs.

    Args:
        input_fasta: Input FASTA path.
        output_fasta: Output FASTA path.
        sequence_ids: Set of sequence IDs to keep.

    Returns:
        Number of sequences written.
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0

    with open(output_fasta, "w") as out:
        for seq_id, description, sequence in iter_fasta(input_fasta):
            if seq_id in sequence_ids:
                out.write(f">{seq_id} {description}\n")
                for i in range(0, len(sequence), 60):
                    out.write(sequence[i:i + 60] + "\n")
                n_written += 1

    logger.info(f"Subset {n_written} sequences")
    return n_written


# =============================================================================
# Database Manager
# =============================================================================


class DatabaseManager:
    """Manages protein reference databases.

    Handles downloading, caching, and formatting databases
    for homology searches.

    Example:
        >>> manager = DatabaseManager()
        >>> db_info = manager.get_database("swissprot_plants")
        >>> print(db_info.formatted_path)  # Path to Diamond database
    """

    # Predefined database configurations
    PREDEFINED_DATABASES = {
        "swissprot": {
            "name": "Swiss-Prot (complete)",
            "db_type": DatabaseType.SWISSPROT,
            "download_func": download_swissprot,
            "download_args": {},
        },
        "swissprot_plants": {
            "name": "Swiss-Prot (plants)",
            "db_type": DatabaseType.SWISSPROT,
            "download_func": download_swissprot,
            "download_args": {"taxonomy": "plants"},
        },
        "swissprot_fungi": {
            "name": "Swiss-Prot (fungi)",
            "db_type": DatabaseType.SWISSPROT,
            "download_func": download_swissprot,
            "download_args": {"taxonomy": "fungi"},
        },
        "uniref90": {
            "name": "UniRef90",
            "db_type": DatabaseType.UNIREF90,
            "download_func": download_uniref,
            "download_args": {"level": 90},
        },
        "uniref50": {
            "name": "UniRef50",
            "db_type": DatabaseType.UNIREF50,
            "download_func": download_uniref,
            "download_args": {"level": 50},
        },
    }

    def __init__(
        self,
        cache_dir: Path | str | None = None,
        auto_format: bool = True,
    ) -> None:
        """Initialize database manager.

        Args:
            cache_dir: Directory for caching databases.
            auto_format: Automatically format databases after download.
        """
        self.cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.auto_format = auto_format
        self._databases: dict[str, DatabaseInfo] = {}

        logger.info(f"Database cache directory: {self.cache_dir}")

    def list_available(self) -> list[str]:
        """List available predefined databases.

        Returns:
            List of database names.
        """
        return list(self.PREDEFINED_DATABASES.keys())

    def list_cached(self) -> list[DatabaseInfo]:
        """List databases in the cache.

        Returns:
            List of DatabaseInfo objects for cached databases.
        """
        cached = []
        for db_name, db_info in self._databases.items():
            if db_info.is_downloaded:
                cached.append(db_info)
        return cached

    def get_database(
        self,
        name: str,
        force_download: bool = False,
    ) -> DatabaseInfo:
        """Get a database, downloading if necessary.

        Args:
            name: Database name (predefined or custom path).
            force_download: Force re-download even if cached.

        Returns:
            DatabaseInfo with database paths.

        Raises:
            ValueError: If database name is unknown.
        """
        # Check if it's a predefined database
        if name in self.PREDEFINED_DATABASES:
            return self._get_predefined(name, force_download)

        # Check if it's a path to a custom database
        custom_path = Path(name)
        if custom_path.exists():
            return self._register_custom(name, custom_path)

        raise ValueError(
            f"Unknown database '{name}'. "
            f"Available: {', '.join(self.list_available())}"
        )

    def _get_predefined(
        self,
        name: str,
        force_download: bool = False,
    ) -> DatabaseInfo:
        """Get a predefined database.

        Args:
            name: Predefined database name.
            force_download: Force re-download.

        Returns:
            DatabaseInfo with database paths.
        """
        config = self.PREDEFINED_DATABASES[name]

        # Check cache
        if name in self._databases and not force_download:
            db_info = self._databases[name]
            if db_info.is_downloaded:
                logger.info(f"Using cached database: {name}")
                return db_info

        # Download
        logger.info(f"Downloading database: {name}")
        db_dir = self.cache_dir / name
        download_func = config["download_func"]
        download_args = config["download_args"].copy()
        download_args["output_dir"] = db_dir

        db_info = download_func(**download_args)

        # Format if requested
        if self.auto_format:
            self._format_database(db_info)

        self._databases[name] = db_info
        return db_info

    def _register_custom(
        self,
        name: str,
        path: Path,
    ) -> DatabaseInfo:
        """Register a custom database.

        Args:
            name: Name to register database under.
            path: Path to FASTA file.

        Returns:
            DatabaseInfo with database paths.
        """
        path = Path(path)

        db_info = DatabaseInfo(
            name=f"Custom: {path.name}",
            db_type=DatabaseType.CUSTOM,
            path=path,
        )

        if self.auto_format:
            self._format_database(db_info)

        self._databases[name] = db_info
        return db_info

    def _format_database(self, db_info: DatabaseInfo) -> None:
        """Format database for Diamond searching.

        Args:
            db_info: Database info to format.
        """
        if db_info.path is None:
            logger.warning("Cannot format database: no path set")
            return

        from helixforge.homology.search import HomologySearch, SearchTool

        # Use Diamond to format
        output_path = db_info.path.with_suffix("")
        searcher = HomologySearch(tool=SearchTool.DIAMOND)

        try:
            formatted_path = searcher.format_database(db_info.path, output_path)
            db_info.formatted_path = formatted_path
            logger.info(f"Formatted database: {formatted_path}")
        except Exception as e:
            logger.warning(f"Failed to format database: {e}")

    def format_database(
        self,
        name: str,
        tool: str = "diamond",
    ) -> Path:
        """Manually format a cached database.

        Args:
            name: Database name.
            tool: Formatting tool ("diamond" or "mmseqs").

        Returns:
            Path to formatted database.
        """
        if name not in self._databases:
            raise ValueError(f"Database not cached: {name}")

        db_info = self._databases[name]
        if db_info.path is None:
            raise ValueError(f"Database has no path: {name}")

        from helixforge.homology.search import HomologySearch, SearchTool

        search_tool = SearchTool(tool)
        output_path = db_info.path.with_suffix("")

        searcher = HomologySearch(tool=search_tool)
        formatted_path = searcher.format_database(db_info.path, output_path)

        db_info.formatted_path = formatted_path
        return formatted_path


# =============================================================================
# CLI Helpers
# =============================================================================


def setup_database(
    database_name: str,
    cache_dir: Path | str | None = None,
    force: bool = False,
) -> DatabaseInfo:
    """High-level function to setup a database for searching.

    Convenience function for CLI usage.

    Args:
        database_name: Name of database to setup.
        cache_dir: Cache directory. Uses default if None.
        force: Force re-download.

    Returns:
        DatabaseInfo with formatted database paths.
    """
    manager = DatabaseManager(cache_dir=cache_dir)
    return manager.get_database(database_name, force_download=force)


def list_databases() -> dict[str, str]:
    """List available databases with descriptions.

    Returns:
        {name: description} mapping.
    """
    descriptions = {
        "swissprot": "Swiss-Prot complete (manually curated, all organisms)",
        "swissprot_plants": "Swiss-Prot plants only (Viridiplantae)",
        "swissprot_fungi": "Swiss-Prot fungi only",
        "uniref90": "UniRef90 clustered proteins (large)",
        "uniref50": "UniRef50 clustered proteins (smaller)",
    }
    return descriptions
