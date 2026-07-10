"""OPTIONAL Helixer wrapper — GPU, off the default path."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from helixforge.prep._subprocess import run_tool
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


@dataclass
class HelixerResult:
    """The GFF3 and HDF5 produced by a Helixer run."""

    gff3: Path
    hdf5: Path


def run_helixer(
    genome_fasta: str | Path,
    out_dir: str | Path,
    model: str = "land_plant_v0.3",
    lineage: str = "land_plant",
    subseq_len: int = 64152,
    helixer_bin: str = "Helixer.py",
    extra_args: Iterable[str | Path] | None = None,
) -> HelixerResult:
    """Run Helixer on ``genome_fasta`` → :class:`HelixerResult` (gff3, hdf5).

    REQUIRES A GPU and is OFF the default prep path — most users already have
    Helixer output and pass it straight into ``PipelineConfig.helixer_gff3`` /
    ``helixer_h5``. ``lineage`` selects the bundled model weights (``model`` is
    recorded for provenance); ``subseq_len`` is Helixer's subsequence length.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gff3 = out_dir / "helixer.gff3"
    hdf5 = out_dir / "helixer.h5"
    _log.info("run_helixer (GPU): lineage=%s model=%s", lineage, model)
    argv: list[str | Path | int] = [
        helixer_bin,
        "--fasta-path",
        genome_fasta,
        "--lineage",
        lineage,
        "--subsequence-length",
        subseq_len,
        "--gff-output-path",
        gff3,
    ]
    if extra_args:
        argv += list(extra_args)
    run_tool(argv, log=_log)
    return HelixerResult(gff3=gff3, hdf5=hdf5)
