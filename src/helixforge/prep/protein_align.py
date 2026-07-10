"""miniprot protein-to-genome alignment wrapper."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from helixforge.prep._subprocess import output_is_fresh, run_tool
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


def run_miniprot(
    genome_fasta: str | Path,
    proteome_fasta: str | Path,
    out_gff: str | Path,
    threads: int = 8,
    miniprot_bin: str = "miniprot",
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
) -> Path:
    """``miniprot -t <threads> --gff <genome> <proteome>`` → ``out_gff`` (returned).

    miniprot writes the GFF3 to stdout; it is redirected to ``out_gff`` via the
    subprocess (never a shell ``>``). ``extra_args`` (e.g. index/preset flags)
    are appended verbatim. Skips when ``out_gff`` is already fresh
    (Phase 22 §3.1) unless ``force``.
    """
    out_gff = Path(out_gff)
    if not force and output_is_fresh(out_gff, [genome_fasta, proteome_fasta]):
        _log.info("skip miniprot: %s up to date", out_gff.name)
        return out_gff
    argv: list[str | Path | int] = [
        miniprot_bin,
        "-t",
        threads,
        "--gff",
        genome_fasta,
        proteome_fasta,
    ]
    if extra_args:
        argv += list(extra_args)
    run_tool(argv, stdout_path=out_gff, log=_log)
    return out_gff
