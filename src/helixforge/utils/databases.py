"""Reference protein database download, caching, and Diamond formatting."""

from __future__ import annotations

import gzip
import hashlib
import shutil
import subprocess
from enum import Enum
from pathlib import Path

from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import attrs

from helixforge.utils.logging import get_logger

_log = get_logger(__name__)

UNIPROT_BASE = "https://ftp.uniprot.org/pub/databases/uniprot"
SWISSPROT_URL = (
    f"{UNIPROT_BASE}/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
)
TREMBL_URL = (
    f"{UNIPROT_BASE}/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
)
UNIREF_URL = f"{UNIPROT_BASE}/uniref/uniref{{level}}/uniref{{level}}.fasta.gz"

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

DEFAULT_CACHE_DIR = Path.home() / ".helixforge" / "databases"

TAX_DIVISION_BASE = (
    f"{UNIPROT_BASE}/current_release/knowledgebase/taxonomic_divisions"
)

# UniProt publishes full Swiss-Prot/TrEMBL and the UniRef clusters as FASTA
# (.fasta.gz), but the taxonomic-division subsets only as flat files
# (.dat.gz). The "format" field records which: "dat" entries are converted to
# FASTA via convert_dat_to_fasta before Diamond formatting; "fasta" entries
# are passed to Diamond directly.
PREDEFINED_DATABASES: dict[str, dict[str, str]] = {
    "swissprot": {
        "name": "Swiss-Prot (complete)",
        "description": "Curated, high-quality annotations — all organisms",
        "url": SWISSPROT_URL,
        "filename": "uniprot_sprot.fasta.gz",
        "format": "fasta",
    },
    "swissprot_plants": {
        "name": "Swiss-Prot (plants)",
        "description": "Curated plant proteins (Viridiplantae)",
        "url": f"{TAX_DIVISION_BASE}/uniprot_sprot_plants.dat.gz",
        "filename": "uniprot_sprot_plants.dat.gz",
        "format": "dat",
    },
    "swissprot_fungi": {
        "name": "Swiss-Prot (fungi)",
        "description": "Curated fungal proteins",
        "url": f"{TAX_DIVISION_BASE}/uniprot_sprot_fungi.dat.gz",
        "filename": "uniprot_sprot_fungi.dat.gz",
        "format": "dat",
    },
    "uniref90": {
        "name": "UniRef90",
        "description": "Clustered at 90% identity — medium size",
        "url": UNIREF_URL.format(level=90),
        "filename": "uniref90.fasta.gz",
        "format": "fasta",
    },
    "uniref50": {
        "name": "UniRef50",
        "description": "Clustered at 50% identity — smaller",
        "url": UNIREF_URL.format(level=50),
        "filename": "uniref50.fasta.gz",
        "format": "fasta",
    },
    "trembl": {
        "name": "TrEMBL (complete)",
        "description": "Automated annotations — all organisms (very large)",
        "url": TREMBL_URL,
        "filename": "uniprot_trembl.fasta.gz",
        "format": "fasta",
    },
}


class DatabaseType(Enum):
    """Types of protein databases."""

    SWISSPROT = "swissprot"
    TREMBL = "trembl"
    UNIREF100 = "uniref100"
    UNIREF90 = "uniref90"
    UNIREF50 = "uniref50"
    CUSTOM = "custom"


def _db_type_from_name(name: str) -> DatabaseType:
    """Infer ``DatabaseType`` from a predefined database name."""
    if name.startswith("swissprot"):
        return DatabaseType.SWISSPROT
    if name == "trembl":
        return DatabaseType.TREMBL
    _uniref_map: dict[str, DatabaseType] = {
        "uniref100": DatabaseType.UNIREF100,
        "uniref90": DatabaseType.UNIREF90,
        "uniref50": DatabaseType.UNIREF50,
    }
    return _uniref_map.get(name, DatabaseType.CUSTOM)


@attrs.define
class DatabaseInfo:
    """Metadata for a protein database."""

    name: str
    db_type: DatabaseType
    path: Path | None = None
    formatted_path: Path | None = None
    download_url: str | None = None
    n_sequences: int | None = None
    checksum: str | None = None

    @property
    def is_downloaded(self) -> bool:
        """True if the FASTA exists on disk."""
        return self.path is not None and self.path.exists()

    @property
    def is_formatted(self) -> bool:
        """True if a Diamond ``.dmnd`` file exists."""
        if self.formatted_path is None:
            return False
        dmnd = self.formatted_path.with_suffix(".dmnd")
        return dmnd.exists() or self.formatted_path.exists()


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------


def download_file(url: str, output_path: Path, chunk_size: int = 8192) -> Path:
    """Download ``url`` to ``output_path`` with progress logging."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _log.info("downloading %s", url)
    _log.info("  -> %s", output_path)

    try:
        with urlopen(url, timeout=30) as resp:
            total = resp.headers.get("Content-Length")
            total = int(total) if total else None
            if total:
                _log.info("  size: %.1f MB", total / 1024 / 1024)
            downloaded = 0
            last_pct = 0
            with open(output_path, "wb") as fh:
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    fh.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = int(100 * downloaded / total)
                        if pct >= last_pct + 10:
                            _log.info("  progress: %d%%", pct)
                            last_pct = pct
            _log.info("  downloaded: %.1f MB", downloaded / 1024 / 1024)
            return output_path
    except (URLError, HTTPError):
        if output_path.exists():
            output_path.unlink()
        raise


def calculate_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """Return the MD5 hex digest of ``file_path``."""
    md5 = hashlib.md5()
    with open(file_path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def decompress_gzip(
    gz_path: Path,
    output_path: Path | None = None,
    keep_original: bool = False,
) -> Path:
    """Decompress a ``.gz`` file. Removes the original unless ``keep_original``."""
    gz_path = Path(gz_path)
    if output_path is None:
        output_path = gz_path.with_suffix("")
        if output_path.suffix == "":
            output_path = gz_path.parent / (gz_path.stem + ".fasta")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _log.info("decompressing %s", gz_path.name)
    with gzip.open(gz_path, "rb") as fin, open(output_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    if not keep_original:
        gz_path.unlink()
    _log.info("decompressed to %s", output_path)
    return output_path


def count_fasta_sequences(fasta_path: Path) -> int:
    """Count ``>`` header lines in a (possibly gzipped) FASTA."""
    fasta_path = Path(fasta_path)
    count = 0
    if fasta_path.suffix == ".gz":
        fh = gzip.open(fasta_path, "rt")
    else:
        fh = open(fasta_path, "r")
    try:
        for line in fh:
            if line.startswith(">"):
                count += 1
    finally:
        fh.close()
    return count


def convert_dat_to_fasta(dat_path: Path, fasta_path: Path) -> int:
    """Convert UniProt flat-file (.dat) to FASTA.

    Args:
        dat_path: Input .dat file.
        fasta_path: Output .fasta file.

    Returns:
        Number of sequences written.
    """
    dat_path = Path(dat_path)
    fasta_path = Path(fasta_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    fin = gzip.open(dat_path, "rt") if dat_path.suffix == ".gz" else open(dat_path, "r")

    entry_name = ""
    accession = ""
    description = ""
    organism = ""
    seq_parts: list[str] = []
    in_sequence = False
    n_written = 0

    def _reset() -> None:
        nonlocal entry_name, accession, description, organism, seq_parts, in_sequence
        entry_name = ""
        accession = ""
        description = ""
        organism = ""
        seq_parts = []
        in_sequence = False

    try:
        with open(fasta_path, "w") as out:
            for line in fin:
                if line.startswith("ID   "):
                    _reset()
                    tokens = line[5:].split()
                    entry_name = tokens[0] if tokens else ""
                elif line.startswith("AC   ") and not accession:
                    accession = line[5:].split(";")[0].strip()
                elif line.startswith("DE   ") and "RecName: Full=" in line and not description:
                    description = line.split("RecName: Full=", 1)[1].split(";")[0].strip()
                elif line.startswith("OS   ") and not organism:
                    organism = line[5:].strip().rstrip(".")
                elif line.startswith("SQ   "):
                    in_sequence = True
                elif line.startswith("     ") and in_sequence:
                    seq_parts.append("".join(c for c in line if c.isalpha()))
                elif line.startswith("//"):
                    seq = "".join(seq_parts)
                    if accession and seq:
                        header = f">sp|{accession}|{entry_name}"
                        if description:
                            header += f" {description}"
                        if organism:
                            header += f" OS={organism}"
                        out.write(header + "\n")
                        for i in range(0, len(seq), 60):
                            out.write(seq[i : i + 60] + "\n")
                        n_written += 1
                    _reset()
    finally:
        fin.close()

    _log.info("converted %d entries from %s", n_written, dat_path.name)
    return n_written


# ---------------------------------------------------------------------------
# DatabaseManager
# ---------------------------------------------------------------------------


class DatabaseManager:
    """Download, cache, and Diamond-format reference protein databases."""

    def __init__(
        self,
        cache_dir: Path | str | None = None,
        diamond_bin: str = "diamond",
    ) -> None:
        self.cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.diamond_bin = diamond_bin
        self._databases: dict[str, DatabaseInfo] = {}
        _log.info("database cache: %s", self.cache_dir)

    def list_available(self) -> list[str]:
        """Return the names of all predefined databases."""
        return list(PREDEFINED_DATABASES.keys())

    def get_database(
        self,
        name: str,
        force_download: bool = False,
    ) -> DatabaseInfo:
        """Get a database by name or path, downloading/formatting as needed.

        Args:
            name: Predefined name (e.g. ``swissprot``) or a local FASTA path.
            force_download: Re-download even when the file is already cached.

        Returns:
            ``DatabaseInfo`` with ``path`` and ``formatted_path`` populated.

        Raises:
            ValueError: Unknown name and not an existing file path.
        """
        if name in PREDEFINED_DATABASES:
            return self._get_predefined(name, force_download)
        custom_path = Path(name)
        if custom_path.exists():
            return self._register_custom(name, custom_path)
        raise ValueError(
            f"Unknown database {name!r}. Available: {', '.join(self.list_available())}"
        )

    # -- internals ----------------------------------------------------------

    def _get_predefined(self, name: str, force: bool) -> DatabaseInfo:
        if name in self._databases and not force:
            info = self._databases[name]
            if info.is_downloaded:
                _log.info("using cached database: %s", name)
                return info

        cfg = PREDEFINED_DATABASES[name]
        db_dir = self.cache_dir / name
        db_dir.mkdir(parents=True, exist_ok=True)
        gz_path = db_dir / cfg["filename"]

        download_file(cfg["url"], gz_path)
        checksum = calculate_md5(gz_path)
        decompressed = decompress_gzip(gz_path)
        if cfg.get("format") == "dat":
            fasta_path = decompressed.with_suffix(".fasta")
            convert_dat_to_fasta(decompressed, fasta_path)
        else:
            fasta_path = decompressed
        n_seq = count_fasta_sequences(fasta_path)

        info = DatabaseInfo(
            name=cfg["name"],
            db_type=_db_type_from_name(name),
            path=fasta_path,
            download_url=cfg["url"],
            n_sequences=n_seq,
            checksum=checksum,
        )
        self.format_database(info)
        self._databases[name] = info
        return info

    def _register_custom(self, name: str, path: Path) -> DatabaseInfo:
        info = DatabaseInfo(
            name=f"Custom: {path.name}",
            db_type=DatabaseType.CUSTOM,
            path=path,
        )
        self.format_database(info)
        self._databases[name] = info
        return info

    def format_database(self, db_info: DatabaseInfo) -> None:
        """Run ``diamond makedb`` on the FASTA to produce a ``.dmnd`` file."""
        if db_info.path is None:
            _log.warning("cannot format database: no path set")
            return
        output = db_info.path.with_suffix("")
        argv = [
            self.diamond_bin,
            "makedb",
            "--in",
            str(db_info.path),
            "--db",
            str(output),
        ]
        _log.info("formatting: %s", " ".join(argv))
        try:
            subprocess.run(argv, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"diamond not found: {self.diamond_bin!r}. "
                "Install Diamond and put it on PATH, or pass --format-tool."
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"diamond makedb failed (exit {exc.returncode}).\n"
                f"argv: {argv}\nstderr: {(exc.stderr or '')[-2000:]}"
            ) from exc
        db_info.formatted_path = output.with_suffix(".dmnd")
        _log.info("formatted database: %s", db_info.formatted_path)
