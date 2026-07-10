"""Read alignment wrappers: STAR and HISAT2."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from helixforge.prep._subprocess import (
    output_is_fresh,
    run_pipe,
    run_tool,
    with_retries,
)
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)

_FASTA_SUFFIXES = (
    ".fasta.gz",
    ".fna.gz",
    ".fa.gz",
    ".fasta",
    ".fna",
    ".fa",
)


@dataclass
class AlignResult:
    """Products of an alignment: the BAM, optional SJ tab, optional log."""

    bam: Path
    sj_tab: Path | None = None
    log: Path | None = None


def _is_fasta(path: str | Path) -> bool:
    """True if ``path`` looks like a FASTA (so we should build an index)."""
    return str(path).lower().endswith(_FASTA_SUFFIXES)


def _strip_fasta_suffix(path: str | Path) -> str:
    """Drop a trailing FASTA suffix to derive an index prefix from a FASTA."""
    s = str(path)
    low = s.lower()
    for suf in _FASTA_SUFFIXES:
        if low.endswith(suf):
            return s[: -len(suf)]
    return s


def _reads_have_gz(reads: Iterable[str | Path]) -> bool:
    return any(str(r).lower().endswith(".gz") for r in reads)


def build_star_index(
    genome_fasta: str | Path,
    index_dir: str | Path,
    threads: int = 8,
    star_bin: str = "STAR",
    sjdb_gtf: str | Path | None = None,
    force: bool = False,
) -> Path:
    """``STAR --runMode genomeGenerate`` → the genome index dir (returned).

    Skips the (expensive) build when the index already exists — STAR's
    ``SAindex`` marker is fresher than the genome FASTA — unless ``force``.
    """
    index_dir = Path(index_dir)
    if not force and output_is_fresh(index_dir / "SAindex", [genome_fasta]):
        _log.info("skip STAR genomeGenerate: index %s up to date", index_dir)
        return index_dir
    index_dir.mkdir(parents=True, exist_ok=True)
    argv: list[str | Path | int] = [
        star_bin,
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        index_dir,
        "--genomeFastaFiles",
        genome_fasta,
        "--runThreadN",
        threads,
    ]
    if sjdb_gtf is not None:
        argv += ["--sjdbGTFfile", sjdb_gtf]
    run_tool(argv, log=_log)
    return index_dir


def run_star(
    genome_dir_or_fasta: str | Path,
    reads: Iterable[str | Path],
    out_prefix: str | Path,
    threads: int = 8,
    star_bin: str = "STAR",
    sjdb_gtf: str | Path | None = None,
    two_pass: bool = False,
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
    retries: int = 0,
) -> AlignResult:
    """Align ``reads`` with STAR → :class:`AlignResult` (BAM + SJ.out.tab + log).

    If ``genome_dir_or_fasta`` is a FASTA, a genome index is built first into
    ``{out_prefix}_star_index`` and used for alignment; if it is an existing
    index directory it is used directly (no rebuild). ``reads`` is a list of
    FASTQ paths (two = paired, one = single); gzipped reads add
    ``--readFilesCommand zcat``. ``two_pass`` enables ``--twopassMode Basic``.

    Skips the whole sample (index build + alignment) when the sorted BAM is
    already fresh (Phase 22 §3.1) unless ``force``.
    """
    reads_list = list(reads)
    out_bam = Path(f"{out_prefix}Aligned.sortedByCoord.out.bam")
    if not force and output_is_fresh(out_bam, reads_list):
        _log.info("skip STAR align: %s up to date", out_bam.name)
        return AlignResult(
            bam=out_bam,
            sj_tab=Path(f"{out_prefix}SJ.out.tab"),
            log=Path(f"{out_prefix}Log.final.out"),
        )
    if _is_fasta(genome_dir_or_fasta):
        genome_dir: str | Path = build_star_index(
            genome_dir_or_fasta,
            f"{out_prefix}_star_index",
            threads=threads,
            star_bin=star_bin,
            sjdb_gtf=sjdb_gtf,
            force=force,
        )
    else:
        genome_dir = genome_dir_or_fasta

    argv: list[str | Path | int] = [
        star_bin,
        "--runMode",
        "alignReads",
        "--genomeDir",
        genome_dir,
        "--readFilesIn",
        *reads_list,
        "--runThreadN",
        threads,
        "--outFileNamePrefix",
        out_prefix,
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
    ]
    if _reads_have_gz(reads_list):
        argv += ["--readFilesCommand", "zcat"]
    if sjdb_gtf is not None:
        argv += ["--sjdbGTFfile", sjdb_gtf]
    if two_pass:
        argv += ["--twopassMode", "Basic"]
    if extra_args:
        argv += list(extra_args)
    # Alignment fully regenerates its BAM → idempotent → safe to retry (D5).
    run_tool(argv, log=_log, retries=retries)

    return AlignResult(
        bam=Path(f"{out_prefix}Aligned.sortedByCoord.out.bam"),
        sj_tab=Path(f"{out_prefix}SJ.out.tab"),
        log=Path(f"{out_prefix}Log.final.out"),
    )


def build_hisat2_index(
    genome_fasta: str | Path,
    index_prefix: str | Path,
    threads: int = 8,
    hisat2_bin: str = "hisat2",
    force: bool = False,
) -> str | Path:
    """``hisat2-build`` → the index prefix (returned). Skips if already built."""
    if not force and output_is_fresh(Path(f"{index_prefix}.1.ht2"), [genome_fasta]):
        _log.info("skip hisat2-build: index %s up to date", index_prefix)
        return index_prefix
    argv: list[str | Path | int] = [
        f"{hisat2_bin}-build",
        "-p",
        threads,
        genome_fasta,
        index_prefix,
    ]
    run_tool(argv, log=_log)
    return index_prefix


def run_hisat2(
    index_prefix: str | Path,
    reads: Iterable[str | Path],
    out_bam: str | Path,
    threads: int = 8,
    hisat2_bin: str = "hisat2",
    samtools_bin: str = "samtools",
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
    retries: int = 0,
) -> AlignResult:
    """Align ``reads`` with HISAT2, piped into ``samtools sort`` → BAM (indexed).

    If ``index_prefix`` is a FASTA, the index is built first (``hisat2-build``)
    and its prefix used. ``reads`` is a list of FASTQ paths (two = paired, one =
    single). HISAT2 does not emit a STAR-style ``SJ.out.tab``, so
    ``AlignResult.sj_tab`` is ``None``; ``log`` is the HISAT2 summary file.

    Skips alignment when ``out_bam`` is already fresh (Phase 22 §3.1) unless
    ``force``.
    """
    reads_list = list(reads)
    out_bam = Path(out_bam)
    summary = Path(f"{out_bam}.summary")
    if not force and output_is_fresh(out_bam, reads_list):
        _log.info("skip HISAT2 align: %s up to date", out_bam.name)
        return AlignResult(bam=out_bam, sj_tab=None, log=summary)
    idx: str | Path
    if _is_fasta(index_prefix):
        idx = _strip_fasta_suffix(index_prefix)
        build_hisat2_index(
            index_prefix, idx, threads=threads, hisat2_bin=hisat2_bin, force=force
        )
    else:
        idx = index_prefix

    read_args: list[str | Path]
    if len(reads_list) == 2:
        read_args = ["-1", reads_list[0], "-2", reads_list[1]]
    else:
        read_args = ["-U", reads_list[0]]
    hisat2_argv: list[str | Path | int] = [
        hisat2_bin,
        "-p",
        threads,
        "-x",
        idx,
        *read_args,
        "--summary-file",
        summary,
    ]
    if extra_args:
        hisat2_argv += list(extra_args)
    sort_argv: list[str | Path | int] = [
        samtools_bin,
        "sort",
        "-@",
        threads,
        "-o",
        out_bam,
        "-",
    ]
    # The align|sort pipe regenerates the BAM wholesale → idempotent → retryable.
    with_retries(lambda: run_pipe([hisat2_argv, sort_argv], log=_log), retries=retries)
    run_tool([samtools_bin, "index", out_bam], log=_log)

    return AlignResult(bam=out_bam, sj_tab=None, log=summary)
