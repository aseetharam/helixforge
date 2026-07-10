"""Per-format input sniffers with contextual errors."""

from __future__ import annotations

import os
from pathlib import Path

__all__ = [
    "FormatError",
    "sniff_fasta",
    "sniff_gff3",
    "sniff_gtf",
    "sniff_star_sj",
    "sniff_bam",
    "sniff_input",
]


class FormatError(ValueError):
    """A malformed/non-standard input file, located by file + line + value.

    Subclasses ``ValueError`` so callers that already catch ``ValueError`` (e.g.
    around ``parse_star_sj_tab``) are unaffected, while gaining the contextual
    ``path``/``line``/``value`` fields and a message that names them.
    """

    def __init__(
        self,
        message: str,
        *,
        path: str | os.PathLike[str] | None = None,
        line: int | None = None,
        value: object | None = None,
    ) -> None:
        self.path = str(path) if path is not None else None
        self.line = line
        self.value = value
        loc = []
        if self.path is not None:
            loc.append(self.path)
        if line is not None:
            loc.append(f"line {line}")
        prefix = f"{':'.join(loc)}: " if loc else ""
        suffix = f" (value: {value!r})" if value is not None else ""
        super().__init__(f"{prefix}{message}{suffix}")


# Strand / phase vocabularies (GFF3 spec; internally only +/- are used, but the
# spec allows '.' for strandless and '?' for unknown, accept both in the
# *sniffer* so a valid record is never falsely rejected).
_GFF_STRANDS = frozenset({"+", "-", ".", "?"})
_GFF_PHASES = frozenset({"0", "1", "2", "."})
_STAR_STRAND_CODES = frozenset({"0", "1", "2"})


def _require_exists(path: str | os.PathLike[str], kind: str) -> Path:
    p = Path(path)
    if not p.exists():
        raise FormatError(f"{kind} file not found", path=path)
    return p


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------


def sniff_fasta(path: str | os.PathLike[str]) -> None:
    """Validate a FASTA: first char ``>``, no bare CR, no empty records.

    * The first non-blank byte of the file must be ``>`` (a header).
    * No bare carriage return (``\\r`` not part of ``\\r\\n``), a classic
      Mac/Windows line-ending corruption that splits records mid-stream.
    * No empty record: every header must be followed by at least one
      non-blank sequence line before the next header / EOF.

    Raises :class:`FormatError` (with the offending line) on the first problem.
    """
    p = _require_exists(path, "FASTA")
    raw = p.read_bytes()
    if raw == b"":
        raise FormatError("empty FASTA file", path=path)

    # Bare-CR scan over the raw bytes (a \r not immediately followed by \n).
    for i, b in enumerate(raw):
        if b == 0x0D and (i + 1 >= len(raw) or raw[i + 1] != 0x0A):
            # Locate the 1-based line (count \n before this byte) + 1.
            line_no = raw.count(b"\n", 0, i) + 1
            raise FormatError(
                "bare carriage return (\\r), convert line endings to LF",
                path=path,
                line=line_no,
            )

    text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()
    seen_header = False
    have_seq_since_header = False
    header_line_no = 0
    for idx, line in enumerate(lines, start=1):
        if line.strip() == "":
            continue
        if not seen_header:
            # First content line must be a header.
            if not line.startswith(">"):
                raise FormatError(
                    "FASTA must start with a '>' header",
                    path=path,
                    line=idx,
                    value=line[:40],
                )
            seen_header = True
            have_seq_since_header = False
            header_line_no = idx
            continue
        if line.startswith(">"):
            if not have_seq_since_header:
                raise FormatError(
                    "empty FASTA record (header with no sequence)",
                    path=path,
                    line=header_line_no,
                )
            have_seq_since_header = False
            header_line_no = idx
        else:
            have_seq_since_header = True
    if seen_header and not have_seq_since_header:
        raise FormatError(
            "empty FASTA record (header with no sequence)",
            path=path,
            line=header_line_no,
        )


# ---------------------------------------------------------------------------
# GFF3 / GTF
# ---------------------------------------------------------------------------


def _sniff_gxf(
    path: str | os.PathLike[str], *, kind: str, require_version: bool
) -> None:
    """Shared GFF3/GTF feature-line validator.

    Every feature (non-comment, non-blank) line must have ≥ 9 tab fields, a
    strand in :data:`_GFF_STRANDS`, and a phase in :data:`_GFF_PHASES`. GFF3
    additionally must carry a ``##gff-version`` directive.
    """
    p = _require_exists(path, kind)
    saw_version = False
    saw_feature = False
    with p.open() as fh:
        for idx, raw_line in enumerate(fh, start=1):
            line = raw_line.rstrip("\n")
            if line == "":
                continue
            if line.startswith("#"):
                if line.startswith("##gff-version"):
                    saw_version = True
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                raise FormatError(
                    f"{kind} feature line has {len(cols)} tab fields, need 9",
                    path=path,
                    line=idx,
                    value=line[:60],
                )
            strand = cols[6]
            if strand not in _GFF_STRANDS:
                raise FormatError(
                    f"invalid strand (expected one of {sorted(_GFF_STRANDS)})",
                    path=path,
                    line=idx,
                    value=strand,
                )
            phase = cols[7]
            if phase not in _GFF_PHASES:
                raise FormatError(
                    f"invalid phase (expected one of {sorted(_GFF_PHASES)})",
                    path=path,
                    line=idx,
                    value=phase,
                )
            saw_feature = True
    if require_version and not saw_version:
        raise FormatError(
            "GFF3 missing the '##gff-version' directive",
            path=path,
        )
    if not saw_feature:
        raise FormatError(f"{kind} has no feature lines", path=path)


def sniff_gff3(path: str | os.PathLike[str]) -> None:
    """Validate a GFF3 (9 tab fields, strand/phase vocab, ``##gff-version``)."""
    _sniff_gxf(path, kind="GFF3", require_version=True)


def sniff_gtf(path: str | os.PathLike[str]) -> None:
    """Validate a GTF (9 tab fields, strand/phase vocab; no version directive)."""
    _sniff_gxf(path, kind="GTF", require_version=False)


# ---------------------------------------------------------------------------
# STAR SJ.out.tab
# ---------------------------------------------------------------------------


def _star_int(value: str, *, path: str | os.PathLike[str], line: int, col: int) -> int:
    """Parse a STAR SJ integer column or raise a located :class:`FormatError`."""
    try:
        return int(value)
    except ValueError:
        raise FormatError(
            f"non-numeric value in STAR SJ column {col}",
            path=path,
            line=line,
            value=value,
        ) from None


def sniff_star_sj(path: str | os.PathLike[str]) -> None:
    """Validate a STAR ``SJ.out.tab``: numeric columns, strand code ∈ {0,1,2}.

    Columns (1-based): chrom, intron_start, intron_end, strand(0/1/2), motif,
    annotated, unique-reads, multimap-reads, max-overhang. Columns 2,3,4,5,6,7,
    8,9 must be integers; column 4 (strand) must be 0, 1, or 2.
    """
    p = _require_exists(path, "STAR SJ.out.tab")
    saw_row = False
    with p.open() as fh:
        for idx, raw_line in enumerate(fh, start=1):
            line = raw_line.rstrip("\n")
            if line.strip() == "" or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                raise FormatError(
                    f"STAR SJ line has {len(cols)} columns, need 9 (truncated?)",
                    path=path,
                    line=idx,
                    value=line[:60],
                )
            # Numeric columns: 2,3 (coords), 5,6 (motif/annotated codes),
            # 7,8 (read counts), 9 (overhang). Column 1 is the chrom name.
            for col in (2, 3, 5, 6, 7, 8, 9):
                _star_int(cols[col - 1], path=path, line=idx, col=col)
            strand_code = cols[3]
            if strand_code not in _STAR_STRAND_CODES:
                raise FormatError(
                    "invalid STAR strand code (expected 0, 1, or 2)",
                    path=path,
                    line=idx,
                    value=strand_code,
                )
            saw_row = True
    if not saw_row:
        raise FormatError("STAR SJ.out.tab has no junction rows", path=path)


# ---------------------------------------------------------------------------
# BAM
# ---------------------------------------------------------------------------


def sniff_bam(path: str | os.PathLike[str]) -> None:
    """Validate a BAM header: ``@HD`` + ``SO:coordinate`` and ≥ 1 ``@SQ``.

    A name-sorted or unsorted BAM (or one missing its ``@SQ`` reference
    dictionary) cannot be coordinate-fetched; this catches it before the
    junction/coverage stages silently return nothing.
    """
    _require_exists(path, "BAM")
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - pysam is a core dep
        raise FormatError("pysam required to validate a BAM", path=path) from exc
    try:
        af = pysam.AlignmentFile(str(path), "rb", check_sq=False)
    except (OSError, ValueError) as exc:
        raise FormatError(f"cannot open BAM ({exc})", path=path) from exc
    try:
        header = af.header.to_dict()
    finally:
        af.close()
    sq = header.get("SQ") or []
    if not sq:
        raise FormatError(
            "BAM header has no @SQ reference sequences",
            path=path,
        )
    hd = header.get("HD") or {}
    so = hd.get("SO")
    if so != "coordinate":
        raise FormatError(
            "BAM is not coordinate-sorted (@HD SO:coordinate required)",
            path=path,
            value=so,
        )


# ---------------------------------------------------------------------------
# dispatcher
# ---------------------------------------------------------------------------

_SNIFFERS = {
    "fasta": sniff_fasta,
    "gff3": sniff_gff3,
    "gtf": sniff_gtf,
    "star_sj": sniff_star_sj,
    "bam": sniff_bam,
}


def sniff_input(path: str | os.PathLike[str], kind: str) -> None:
    """Dispatch to the sniffer named by ``kind`` (see :data:`_SNIFFERS`)."""
    try:
        sniffer = _SNIFFERS[kind]
    except KeyError:
        raise ValueError(f"unknown input kind {kind!r}") from None
    sniffer(path)
