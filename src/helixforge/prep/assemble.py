"""StringTie assembly wrappers."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from helixforge.prep._subprocess import output_is_fresh, run_tool
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


def run_stringtie(
    bam: str | Path,
    out_gtf: str | Path,
    sample_id: str,
    threads: int = 4,
    stringtie_bin: str = "stringtie",
    guide_gtf: str | Path | None = None,
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
    retries: int = 0,
) -> Path:
    """Assemble one BAM → ``out_gtf`` (returned). ``-l sample_id`` labels ids.

    ``guide_gtf`` (``-G``) supplies a reference annotation to guide assembly.
    Skips when ``out_gtf`` is already fresh unless ``force``.
    ``retries``: StringTie regenerates its GTF wholesale, so an
    allowlisted transient failure is retried; default 0 keeps one attempt.
    """
    out_gtf = Path(out_gtf)
    if not force and output_is_fresh(out_gtf, [bam]):
        _log.info("skip stringtie: %s up to date", out_gtf.name)
        return out_gtf
    argv: list[str | Path | int] = [
        stringtie_bin,
        bam,
        "-o",
        out_gtf,
        "-p",
        threads,
        "-l",
        sample_id,
    ]
    if guide_gtf is not None:
        argv += ["-G", guide_gtf]
    if extra_args:
        argv += list(extra_args)
    run_tool(argv, log=_log, retries=retries)
    return out_gtf


def assemble_samples(
    bams: Iterable[str | Path],
    out_dir: str | Path,
    threads: int = 4,
    workers: int = 1,
    **kw: object,
) -> list[Path]:
    """Assemble each BAM → one GTF in ``out_dir``; return the list of GTF paths.

    The sample id (and the GTF stem) is derived from each BAM's stem, so
    transcript labels are unique per sample. Extra keywords
    (``stringtie_bin``, ``guide_gtf``, ``extra_args``) pass through to
    :func:`run_stringtie`. The returned list feeds ``PipelineConfig.stringtie_list``.

    ``workers``: samples are independent, so with
    ``workers > 1`` StringTie runs over a bounded thread pool (each subprocess
    releases the GIL). Returned paths stay in input (BAM) order regardless of
    completion order. Default 1 keeps the serial behavior byte-identical.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bams_list = list(bams)
    out_gtfs = [out_dir / f"{Path(bam).stem}.gtf" for bam in bams_list]

    def _one(bam: str | Path, out_gtf: Path) -> Path:
        return run_stringtie(bam, out_gtf, Path(bam).stem, threads=threads, **kw)  # type: ignore[arg-type]  # **kw passes through stringtie kwargs

    if workers and workers > 1 and len(bams_list) > 1:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=workers) as ex:
            return list(ex.map(_one, bams_list, out_gtfs))
    return [_one(bam, out_gtf) for bam, out_gtf in zip(bams_list, out_gtfs)]
