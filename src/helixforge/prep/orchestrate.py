"""Evidence-prep orchestrator."""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from helixforge.prep.align import (
    AlignResult,
    build_hisat2_index,
    build_star_index,
    run_hisat2,
    run_star,
)
from helixforge.prep.assemble import assemble_samples
from helixforge.prep.protein_align import run_miniprot
from helixforge.prep.samtools import faidx
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


@dataclass
class PreppedInputs:
    """Staged evidence ready to feed ``PipelineConfig`` (Phase 8)."""

    genome_fasta: str
    bam_paths: list[str] = field(default_factory=list)
    stringtie_list: list[str] = field(default_factory=list)
    star_sj_paths: list[str] = field(default_factory=list)
    miniprot_gff: str | None = None

    def as_config_kwargs(self) -> dict[str, Any]:
        """Return a dict of the fields ``PipelineConfig`` accepts verbatim."""
        return {
            "genome_fasta": self.genome_fasta,
            "bam_paths": list(self.bam_paths),
            "stringtie_list": list(self.stringtie_list),
            "star_sj_paths": list(self.star_sj_paths),
            "miniprot_gff": self.miniprot_gff,
        }


def _normalize_reads(
    rnaseq_reads: Iterable[str | Path | Iterable[str | Path]] | None,
) -> list[list[str]]:
    """Normalise to a list of per-sample read lists.

    Each sample may be a single path (single-end) or a list/tuple of one or two
    paths (single/paired-end).
    """
    if not rnaseq_reads:
        return []
    samples: list[list[str]] = []
    for sample in rnaseq_reads:
        if isinstance(sample, (str, Path)):
            samples.append([str(sample)])
        else:
            samples.append([str(r) for r in sample])
    return samples


def _map_samples(
    fns: list[Callable[[], Any]],
    workers: int,
) -> list[Any]:
    """Run per-sample thunks serially or over a bounded thread pool (order kept).

    Each thunk wraps a subprocess wrapper (STAR/HISAT2/StringTie), which releases
    the GIL inside ``subprocess.run``, so a ``ThreadPoolExecutor`` gives real
    concurrency without pickling the wrappers. ``ThreadPoolExecutor.map`` yields
    results in **input order**, so the returned BAM/SJ lists are deterministic
    regardless of which sample finishes first.
    """
    if workers and workers > 1 and len(fns) > 1:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=workers) as ex:
            return list(ex.map(lambda f: f(), fns))
    return [f() for f in fns]


def prep_evidence(
    genome_fasta: str | Path,
    *,
    rnaseq_reads: Iterable[str | Path | Iterable[str | Path]] | None = None,
    aligner: str = "star",
    proteome: str | Path | None = None,
    out_dir: str | Path,
    threads: int = 8,
    workers: int = 1,
    force: bool = False,
    **bins: str,
) -> PreppedInputs:
    """Produce staged evidence from raw inputs → :class:`PreppedInputs`.

    Stages (each skipped when its inputs are absent), logged with timing:
    ``faidx`` (always) → per-sample align (STAR or HISAT2) → StringTie assemble
    → miniprot (if ``proteome``). For STAR the genome index is built once and
    reused across samples; STAR ``SJ.out.tab`` files are collected (HISAT2 emits
    none). Tool binaries override via keyword: ``star_bin``, ``hisat2_bin``,
    ``samtools_bin``, ``stringtie_bin``, ``miniprot_bin``.

    ``workers``: samples are independent jobs, so
    with ``workers > 1`` the per-sample align + assemble run over a bounded pool,
    wall-time drops from ``Σ(sample_times)`` toward ``max(sample_times)``. It is a
    **concurrency budget**: each aligner is itself ``threads``-parallel, so size
    ``workers × threads`` to the node's cores. Default 1 keeps the serial order
    byte-identical. The shared genome index is always built once, before the pool.

    ``force``: when False (default), each stage whose output is
    already present and up to date is skipped, a re-run of a long prep resumes
    instead of re-aligning every sample. ``force=True`` regenerates everything.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    star_bin: str = bins.get("star_bin", "STAR")
    hisat2_bin: str = bins.get("hisat2_bin", "hisat2")
    samtools_bin: str = bins.get("samtools_bin", "samtools")
    stringtie_bin: str = bins.get("stringtie_bin", "stringtie")
    miniprot_bin: str = bins.get("miniprot_bin", "miniprot")

    def _timed(label: str, fn: Callable[[], Any]) -> Any:
        t0 = time.perf_counter()
        result = fn()
        _log.info("prep stage %s done in %.2fs", label, time.perf_counter() - t0)
        return result

    _timed("faidx", lambda: faidx(genome_fasta, samtools_bin=samtools_bin, force=force))

    bam_paths: list[Path] = []
    sj_paths: list[Path] = []
    stringtie_list: list[Path] = []
    samples = _normalize_reads(rnaseq_reads)
    if samples:
        if aligner == "star":
            index_dir: Path = _timed(
                "star-index",
                lambda: build_star_index(
                    genome_fasta,
                    out_dir / "star_index",
                    threads=threads,
                    star_bin=star_bin,
                    force=force,
                ),
            )
            align_fns: list[Callable[[], Any]] = []
            for i, reads in enumerate(samples):
                prefix = str(out_dir / f"sample{i}_")

                def _star_fn(
                    _reads: list[str] = reads,
                    _prefix: str = prefix,
                ) -> AlignResult:
                    return run_star(
                        index_dir,
                        _reads,
                        _prefix,
                        threads=threads,
                        star_bin=star_bin,
                        force=force,
                    )

                align_fns.append(_star_fn)
            results: list[Any] = _timed(
                "star-align", lambda: _map_samples(align_fns, workers)
            )
            for res in results:
                bam_paths.append(res.bam)
                if res.sj_tab is not None:
                    sj_paths.append(res.sj_tab)
        elif aligner == "hisat2":
            index_prefix = str(out_dir / "hisat2_index")
            _timed(
                "hisat2-index",
                lambda: build_hisat2_index(
                    genome_fasta,
                    index_prefix,
                    threads=threads,
                    hisat2_bin=hisat2_bin,
                    force=force,
                ),
            )
            align_fns = []
            for i, reads in enumerate(samples):
                out_bam = out_dir / f"sample{i}.bam"

                def _hisat2_fn(
                    _reads: list[str] = reads,
                    _out_bam: Path = out_bam,
                ) -> AlignResult:
                    return run_hisat2(
                        index_prefix,
                        _reads,
                        _out_bam,
                        threads=threads,
                        hisat2_bin=hisat2_bin,
                        samtools_bin=samtools_bin,
                        force=force,
                    )

                align_fns.append(_hisat2_fn)
            results = _timed("hisat2-align", lambda: _map_samples(align_fns, workers))
            for res in results:
                bam_paths.append(res.bam)
        else:
            raise ValueError(f"unknown aligner: {aligner!r} (use 'star' or 'hisat2')")

        stringtie_list = _timed(
            "stringtie",
            lambda: assemble_samples(
                bam_paths,
                out_dir / "stringtie",
                threads=threads,
                workers=workers,
                stringtie_bin=stringtie_bin,
                force=force,
            ),
        )
    else:
        _log.info("prep: no RNA-seq reads, skipping align + assemble")

    miniprot_gff: Path | None = None
    if proteome:
        miniprot_gff = _timed(
            "miniprot",
            lambda: run_miniprot(
                genome_fasta,
                proteome,
                out_dir / "miniprot.gff",
                threads=threads,
                miniprot_bin=miniprot_bin,
                force=force,
            ),
        )
    else:
        _log.info("prep: no proteome, skipping miniprot")

    return PreppedInputs(
        genome_fasta=str(genome_fasta),
        bam_paths=[str(b) for b in bam_paths],
        stringtie_list=[str(g) for g in stringtie_list],
        star_sj_paths=[str(s) for s in sj_paths],
        miniprot_gff=(str(miniprot_gff) if miniprot_gff is not None else None),
    )
