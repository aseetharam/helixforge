"""samtools wrappers: sort / index / merge / faidx."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from helixforge.prep._subprocess import output_is_fresh, run_tool
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


def sort(
    in_bam: str | Path,
    out_bam: str | Path,
    threads: int = 4,
    samtools_bin: str = "samtools",
    force: bool = False,
) -> Path:
    """``samtools sort`` → coordinate-sorted ``out_bam`` (returned).

    Skips the sort when ``out_bam`` is already fresh (Phase 22 §3.1) unless
    ``force``.
    """
    out_bam = Path(out_bam)
    if not force and output_is_fresh(out_bam, [in_bam]):
        _log.info("skip samtools sort: %s up to date", out_bam.name)
        return out_bam
    run_tool(
        [samtools_bin, "sort", "-@", threads, "-o", out_bam, in_bam],
        log=_log,
    )
    return out_bam


def index(
    bam: str | Path,
    samtools_bin: str = "samtools",
    force: bool = False,
) -> Path:
    """``samtools index`` → the ``.bai`` index path (returned). Skips if fresh."""
    bam = Path(bam)
    bai = Path(f"{bam}.bai")
    if not force and output_is_fresh(bai, [bam]):
        _log.info("skip samtools index: %s up to date", bai.name)
        return bai
    run_tool([samtools_bin, "index", bam], log=_log)
    return bai


def merge(
    in_bams: Iterable[str | Path],
    out_bam: str | Path,
    threads: int = 4,
    samtools_bin: str = "samtools",
    force: bool = False,
) -> Path:
    """``samtools merge`` several BAMs → ``out_bam`` (returned; ``-f`` overwrite).

    Skips when ``out_bam`` is fresher than every input BAM unless ``force``.
    """
    out_bam = Path(out_bam)
    if not force and output_is_fresh(out_bam, in_bams):
        _log.info("skip samtools merge: %s up to date", out_bam.name)
        return out_bam
    run_tool(
        [samtools_bin, "merge", "-@", threads, "-f", out_bam, *in_bams],
        log=_log,
    )
    return out_bam


def faidx(
    fasta: str | Path,
    samtools_bin: str = "samtools",
    force: bool = False,
) -> Path:
    """``samtools faidx`` → the ``.fai`` index path (returned). Skips if fresh."""
    fasta = Path(fasta)
    fai = Path(f"{fasta}.fai")
    if not force and output_is_fresh(fai, [fasta]):
        _log.info("skip samtools faidx: %s up to date", fai.name)
        return fai
    run_tool([samtools_bin, "faidx", fasta], log=_log)
    return fai
