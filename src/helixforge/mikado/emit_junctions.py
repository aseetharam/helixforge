"""Junction emitters for Mikado serialise."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from helixforge.prep._subprocess import run_tool

if TYPE_CHECKING:
    from helixforge.reconcile.models import SpliceJunction


def junctions_to_portcullis_tab(
    junctions: list[SpliceJunction], out_path: str | Path
) -> Path:
    """Write ``SpliceJunction`` list as a Portcullis-style BED12 junctions file.

    Columns: chrom, chromStart, chromEnd, name, score (read_count, capped 1000),
    strand, thickStart, thickEnd, itemRgb, blockCount(2), blockSizes,
    blockStarts. ``samples`` is preserved in the name (``juncN_sM``). Returns the
    output ``Path``.
    """
    out_path = Path(out_path)
    with out_path.open("w") as fh:
        for i, j in enumerate(junctions):
            left_anchor = 1 if j.donor >= 1 else 0
            chrom_start = j.donor - left_anchor
            chrom_end = j.acceptor + 1
            block2_start = j.acceptor - chrom_start  # = acceptor - donor + left_anchor
            block_sizes = f"{left_anchor},1"
            block_starts = f"0,{block2_start}"
            name = f"junc{i}_s{j.samples}"
            score = min(1000, j.read_count)
            fh.write(
                f"{j.seqid}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{j.strand}\t"
                f"{chrom_start}\t{chrom_end}\t0\t2\t{block_sizes}\t{block_starts}\n"
            )
    return out_path


def run_portcullis(
    genome_fa: str | Path,
    bam_paths: list[str | Path],
    out_dir: str | Path,
    threads: int = 4,
    portcullis_bin: str = "portcullis",
) -> Path:
    """Run ``portcullis full`` and return the filtered pass-junctions ``.tab`` path.

    Portcullis output is preferred over :func:`junctions_to_portcullis_tab`
    (higher precision). Raises ``FileNotFoundError`` if the binary is absent.
    The subprocess goes through the shared :func:`run_tool` so
    a Portcullis failure surfaces with the argv + a stderr tail, and the full
    stderr is streamed to ``<out_dir>/portcullis.stderr.log`` — matching the
    instrumentation every other external tool already has. Mocked in tests.
    """
    if shutil.which(portcullis_bin) is None:
        raise FileNotFoundError(
            f"portcullis binary not found: {portcullis_bin!r}. "
            "Install Portcullis (bioconda) or use junctions_to_portcullis_tab()."
        )
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    argv: list[str] = [
        portcullis_bin,
        "full",
        "-t",
        str(threads),
        "-o",
        str(out_dir),
        str(Path(genome_fa).resolve()),
        *[str(Path(b).resolve()) for b in bam_paths],
    ]
    run_tool(argv, cwd=out_dir, stderr_log=out_dir / "portcullis.stderr.log")
    return out_dir / "3-filt" / "portcullis_filtered.pass.junctions.tab"
