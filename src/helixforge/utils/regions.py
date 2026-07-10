"""Region string parsing and coordinate conversion."""

from __future__ import annotations


def parse_region(region: str) -> tuple[str, int | None, int | None]:
    """Parse a region string ``"seqid:start-end"`` → ``(seqid, start, end)``.

    Coordinates are returned verbatim as integers (no semantic conversion —
    that is the job of ``gff3_to_internal`` / ``internal_to_gff3``). A bare
    ``"seqid"`` returns ``(seqid, None, None)``. Requires ``start <= end`` when
    both are present.
    """
    if not isinstance(region, str) or not region:
        raise ValueError(f"region must be a non-empty string, got {region!r}")

    if ":" not in region:
        return region, None, None

    seqid, _, span = region.rpartition(":")
    if not seqid:
        raise ValueError(f"missing seqid in region {region!r}")
    if "-" not in span:
        raise ValueError(f"malformed region span {span!r} in {region!r}")

    start_s, _, end_s = span.partition("-")
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        raise ValueError(f"non-integer coordinates in region {region!r}") from None

    if start > end:
        raise ValueError(f"region start {start} > end {end} in {region!r}")
    return seqid, start, end


def format_region(seqid: str, start: int | None = None, end: int | None = None) -> str:
    """Inverse of :func:`parse_region`."""
    if start is None and end is None:
        return str(seqid)
    if start is None or end is None:
        raise ValueError("both start and end must be provided (or neither)")
    if start > end:
        raise ValueError(f"start {start} > end {end}")
    return f"{seqid}:{start}-{end}"


def gff3_to_internal(start: int, end: int) -> tuple[int, int]:
    """1-based inclusive ``[start, end]`` → 0-based half-open ``[start-1, end)``."""
    if start < 1:
        raise ValueError(f"GFF3 start must be >= 1, got {start}")
    if end < start:
        raise ValueError(f"GFF3 end {end} < start {start}")
    return start - 1, end


def internal_to_gff3(start: int, end: int) -> tuple[int, int]:
    """0-based half-open ``[start, end)`` → 1-based inclusive ``[start+1, end]``."""
    if start < 0:
        raise ValueError(f"internal start must be >= 0, got {start}")
    if end <= start:
        raise ValueError(f"internal end {end} must be > start {start}")
    return start + 1, end
