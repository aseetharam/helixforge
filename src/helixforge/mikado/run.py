"""Mikado pipeline runner."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from helixforge.prep._subprocess import ToolError, tool_version, with_retries
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


def _parse_version_tuple(raw: str | None) -> tuple[int, ...] | None:
    """Pull a ``(major, minor[, patch])`` int tuple from a version string."""
    if not raw:
        return None
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", raw)
    if not m:
        return None
    return tuple(int(g) for g in m.groups() if g is not None)


def mikado_version(mikado_bin: str = "mikado") -> tuple[int, ...] | None:
    """Return Mikado's version as an int tuple (e.g. ``(2, 3, 4)``) or ``None``.

    Best-effort: a missing binary or unparseable output yields ``None`` (callers
    then fall back to the pinned-2.3.x flag spelling). Used to branch the known
    2.x CLI flag variants in :func:`run_serialise`.
    """
    return _parse_version_tuple(tool_version(mikado_bin, "--version"))


def _serialise_genome_flag(version: tuple[int, ...] | None) -> str:
    """Pick the serialise genome flag for a detected Mikado version.

    Mikado 2.3 renamed ``serialise --genome_fasta`` to ``--genome``. Older 2.x
    used ``--genome_fasta``; unknown/None defaults to the pinned-2.3.x ``--genome``.
    """
    if version is not None and version < (2, 3):
        return "--genome_fasta"
    return "--genome"


def _require_tool(bin_name: str | Path, label: str) -> None:
    """Raise a precise ``FileNotFoundError`` if ``bin_name`` does not resolve.

    A path-like value (contains a separator) must exist on disk; a bare command
    must be found on ``PATH``. Called immediately before a step actually runs, so
    a cold run fails fast naming the missing tool rather than mid-pipeline.
    """
    name = str(bin_name)
    looks_like_path = os.sep in name or (os.altsep and os.altsep in name)
    if looks_like_path:
        if not os.path.exists(name):
            raise FileNotFoundError(f"{label} not found at configured path: {name!r}")
    elif shutil.which(name) is None:
        raise FileNotFoundError(
            f"{label} not found on PATH: {name!r}. Install it or pass an "
            "absolute path via the corresponding *_bin argument."
        )


# Config/results containers use dataclasses (not attrs, which is reserved for
# biological data models).
@dataclass
class MikadoInputs:
    """The Phase-4-emitted inputs that drive a Mikado run."""

    configuration_yaml: str
    scoring_file: str
    genome_fa: str
    protein_db: str
    junctions_tab: str
    external_tsv: str


@dataclass
class MikadoRunResult:
    """All paths produced by a full Mikado run."""

    prepared_gtf: Path
    prepared_fasta: Path
    orfs_bed: Path
    blast_out: Path
    mikado_db: Path
    loci_gff3: Path
    loci_metrics_tsv: Path
    loci_scores_tsv: Path


def _run_cmd(
    argv: list[Any],
    step: str,
    cwd: str | Path | None = None,
    *,
    retries: int = 0,
    retry_exit_codes: Iterable[int] | None = None,
    backoff: float = 0.5,
    _sleep: Callable[[float], object] = time.sleep,
) -> object:
    """Run ``argv`` (list) with logging + timing; raise an informative error.

    ``cwd`` sets the subprocess working directory (needed for TransDecoder,
    which writes its final ``*.transdecoder.*`` outputs to the CWD rather than
    to ``--output_dir``).

    ``retries``: a bounded retry-with-backoff for the **idempotent** Mikado-chain
    externals (only DIAMOND opts in via :func:`run_diamond`). Default 0 ⇒ one
    attempt. **Never** pass ``retries`` for the non-idempotent ``serialise`` DB
    insertion — rely on checkpointing.
    """
    argv_str = [str(a) for a in argv]

    def _attempt() -> object:
        _log.info("[%s] start: %s", step, " ".join(argv_str))
        t0 = time.perf_counter()
        try:
            proc = subprocess.run(
                argv_str,
                check=True,
                capture_output=True,
                text=True,
                cwd=(str(cwd) if cwd is not None else None),
            )
        except subprocess.CalledProcessError as exc:
            stderr_tail = (exc.stderr or "")[-2000:]
            raise ToolError(
                f"Mikado step {step!r} failed (exit {exc.returncode}).\n"
                f"argv: {argv_str}\n"
                f"stderr tail:\n{stderr_tail}",
                returncode=exc.returncode,
                argv=argv_str,
            ) from exc
        _log.info("[%s] done in %.2fs", step, time.perf_counter() - t0)
        return proc

    return with_retries(
        _attempt,
        retries=retries,
        retry_exit_codes=retry_exit_codes,
        backoff=backoff,
        _sleep=_sleep,
        log=_log,
    )


def _bin(name: str, bin_dir: str | Path | None = None) -> str:
    return str(Path(bin_dir) / name) if bin_dir else name


def run_prepare(
    configuration_yaml: str | Path,
    out_dir: str | Path,
    procs: int = 1,
    mikado_bin: str = "mikado",
) -> tuple[Path, Path]:
    """``mikado prepare`` → ``(mikado_prepared.gtf, mikado_prepared.fasta)``."""
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    argv: list[Any] = [
        mikado_bin,
        "prepare",
        "--configuration",
        Path(configuration_yaml).resolve(),
        "--procs",
        procs,
        "--output-dir",
        out_dir,
    ]
    _run_cmd(argv, "prepare", cwd=out_dir)
    return out_dir / "mikado_prepared.gtf", out_dir / "mikado_prepared.fasta"


def run_transdecoder(
    prepared_fasta: str | Path,
    out_dir: str | Path,
    single_best_only: bool = True,
    transdecoder_bin_dir: str | Path | None = None,
) -> Path:
    """TransDecoder ``LongOrfs`` + ``Predict`` → the ``.transdecoder.bed`` (bed12).

    TransDecoder writes its final ``*.transdecoder.{bed,pep,gff3,cds}`` outputs
    to the working directory (NOT ``--output_dir``, which only holds the
    intermediate checkpoint dir), so both steps run with ``cwd=out_dir`` to keep
    everything self-contained and make the returned path correct.

    All paths are resolved to absolute before being passed to the subprocess —
    a relative ``--output_dir`` combined with ``cwd=out_dir`` caused
    TransDecoder to nest the path and fail its mkdir.
    """
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    prepared_fasta = Path(prepared_fasta).resolve()
    long_orfs = _bin("TransDecoder.LongOrfs", transdecoder_bin_dir)
    predict = _bin("TransDecoder.Predict", transdecoder_bin_dir)

    _run_cmd(
        [long_orfs, "-t", prepared_fasta, "--output_dir", out_dir],
        "transdecoder.LongOrfs",
        cwd=out_dir,
    )
    predict_argv: list[Any] = [predict, "-t", prepared_fasta, "--output_dir", out_dir]
    if single_best_only:
        predict_argv.append("--single_best_only")
    _run_cmd(predict_argv, "transdecoder.Predict", cwd=out_dir)

    return out_dir / f"{prepared_fasta.name}.transdecoder.bed"


def run_diamond(
    prepared_fasta: str | Path,
    protein_db_fasta: str | Path,
    out_dir: str | Path,
    threads: int = 4,
    evalue: float = 1e-5,
    max_target_seqs: int = 5,
    diamond_bin: str = "diamond",
    out_format: str = "xml",
    retries: int = 0,
    retry_exit_codes: Iterable[int] | None = None,
) -> Path:
    """DIAMOND ``makedb`` + ``blastx`` → homology output for ``mikado serialise``.

    ``out_format='xml'`` (BLAST-XML, ``--outfmt 5``) is the default Mikado reads;
    ``out_format='tabular'`` emits ``--outfmt 6`` with Mikado's expected columns.
    Confirm the accepted format against the pinned Mikado serialise.

    ``retries``: DIAMOND ``makedb``/``blastx`` both fully regenerate
    their output, so a transient cluster failure is retried with backoff. Default
    0 keeps one attempt (byte-identical). Serialise is **not** retried.
    """
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    prepared_fasta = Path(prepared_fasta).resolve()
    protein_db_fasta = Path(protein_db_fasta).resolve()
    db_path = out_dir / "diamond_db"
    _run_cmd(
        [diamond_bin, "makedb", "--in", protein_db_fasta, "--db", db_path],
        "diamond.makedb",
        cwd=out_dir,
        retries=retries,
        retry_exit_codes=retry_exit_codes,
    )

    ext: str
    fmt_args: list[str]
    if out_format == "xml":
        ext, fmt_args = "xml", ["--outfmt", "5"]
    else:
        ext = "tsv"
        fmt_args = [
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ]
    blast_out = out_dir / f"mikado_diamond.{ext}"
    _run_cmd(
        [
            diamond_bin,
            "blastx",
            "--query",
            prepared_fasta,
            "--db",
            db_path,
            "--out",
            blast_out,
            "--threads",
            threads,
            "--evalue",
            evalue,
            "--max-target-seqs",
            max_target_seqs,
            *fmt_args,
        ],
        "diamond.blastx",
        cwd=out_dir,
        retries=retries,
        retry_exit_codes=retry_exit_codes,
    )
    return blast_out


def run_serialise(
    configuration_yaml: str | Path,
    prepared_fasta: str | Path,
    orfs_bed: str | Path,
    blast_out: str | Path,
    blast_db: str | Path,
    junctions_tab: str | Path,
    external_tsv: str | Path,
    genome_fa: str | Path,
    out_dir: str | Path,
    mikado_bin: str = "mikado",
    version: tuple[int, ...] | None = None,
) -> Path:
    """``mikado serialise`` → ``mikado.db``.

    Put the SQLite DB on fast local scratch (single-threaded insertion); the
    caller controls this via ``out_dir`` (the pipeline passes a local ``work_dir``,
    never a network FS). The genome flag spelling is branched on the detected
    Mikado version (``--genome`` for 2.3+, ``--genome_fasta`` for older 2.x);
    ``version=None`` auto-detects via :func:`mikado_version`.
    """
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    db_path = out_dir / "mikado.db"
    if version is None:
        version = mikado_version(mikado_bin)
    genome_flag = _serialise_genome_flag(version)
    argv: list[Any] = [
        mikado_bin,
        "serialise",
        "--configuration",
        Path(configuration_yaml).resolve(),
        "--transcripts",
        Path(prepared_fasta).resolve(),
        "--orfs",
        Path(orfs_bed).resolve(),
        "--xml",
        Path(blast_out).resolve(),
        "--blast_targets",
        Path(blast_db).resolve(),
        "--junctions",
        Path(junctions_tab).resolve(),
        "--external-scores",
        Path(external_tsv).resolve(),
        genome_flag,
        Path(genome_fa).resolve(),
        "--output-dir",
        out_dir,
        db_path,
    ]
    _run_cmd(argv, "serialise", cwd=out_dir)
    return db_path


def run_pick(
    configuration_yaml: str | Path,
    scoring_file: str | Path,
    prepared_gtf: str | Path,
    out_dir: str | Path,
    procs: int = 1,
    mikado_bin: str = "mikado",
) -> Path:
    """``mikado pick`` → ``mikado.loci.gff3``."""
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    argv: list[Any] = [
        mikado_bin,
        "pick",
        "--configuration",
        Path(configuration_yaml).resolve(),
        "--scoring-file",
        Path(scoring_file).resolve(),
        "--procs",
        procs,
        "--output-dir",
        out_dir,
        Path(prepared_gtf).resolve(),
    ]
    _run_cmd(argv, "pick", cwd=out_dir)
    return out_dir / "mikado.loci.gff3"


def run_mikado(
    inputs: MikadoInputs,
    out_dir: str | Path,
    procs: int = 1,
    threads: int = 4,
    mikado_bin: str = "mikado",
    diamond_bin: str = "diamond",
    transdecoder_bin_dir: str | Path | None = None,
    diamond_format: str = "xml",
    resume: bool = False,
) -> MikadoRunResult:
    """Run the full chain and return a :class:`MikadoRunResult`.

    Steps: prepare → TransDecoder → DIAMOND → serialise → pick. Each step is
    logged with timing; nonzero exits raise with the argv and stderr tail.
    ``out_dir`` must be caller-provided **local scratch** (the SQLite DB lands
    here; never a network FS).

    ``resume`` makes the chain restartable as a first-class, tested argument
    (M3's monkeypatching of the ``run_*`` funcs is no longer required for a cold
    run): when ``True``, any step whose expected output already exists under
    ``out_dir`` is reused as-is and its external tool is never invoked; when
    ``False`` (default) the full chain runs. Each external tool's binary is
    resolved (and a precise error raised if missing) immediately before — and
    only when — that step actually runs.
    """
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    _log.info("run_mikado: out_dir=%s resume=%s", out_dir, resume)

    prepared_gtf = out_dir / "mikado_prepared.gtf"
    prepared_fasta = out_dir / "mikado_prepared.fasta"
    if resume and prepared_gtf.exists() and prepared_fasta.exists():
        _log.info("[prepare] resume: reusing %s", prepared_gtf.name)
    else:
        _require_tool(mikado_bin, "mikado")
        prepared_gtf, prepared_fasta = run_prepare(
            inputs.configuration_yaml, out_dir, procs=procs, mikado_bin=mikado_bin
        )

    orfs_bed = out_dir / f"{prepared_fasta.name}.transdecoder.bed"
    if resume and orfs_bed.exists():
        _log.info("[transdecoder] resume: reusing %s", orfs_bed.name)
    else:
        _require_tool(
            _bin("TransDecoder.LongOrfs", transdecoder_bin_dir), "TransDecoder"
        )
        orfs_bed = run_transdecoder(
            prepared_fasta, out_dir, transdecoder_bin_dir=transdecoder_bin_dir
        )

    blast_ext = "xml" if diamond_format == "xml" else "tsv"
    blast_out = out_dir / f"mikado_diamond.{blast_ext}"
    if resume and blast_out.exists():
        _log.info("[diamond] resume: reusing %s", blast_out.name)
    else:
        _require_tool(diamond_bin, "diamond")
        blast_out = run_diamond(
            prepared_fasta,
            inputs.protein_db,
            out_dir,
            threads=threads,
            diamond_bin=diamond_bin,
            out_format=diamond_format,
        )

    mikado_db = out_dir / "mikado.db"
    if resume and mikado_db.exists():
        _log.info("[serialise] resume: reusing %s", mikado_db.name)
    else:
        _require_tool(mikado_bin, "mikado")
        mikado_db = run_serialise(
            inputs.configuration_yaml,
            prepared_fasta,
            orfs_bed,
            blast_out,
            inputs.protein_db,
            inputs.junctions_tab,
            inputs.external_tsv,
            inputs.genome_fa,
            out_dir,
            mikado_bin=mikado_bin,
        )

    loci_gff3 = out_dir / "mikado.loci.gff3"
    if resume and loci_gff3.exists():
        _log.info("[pick] resume: reusing %s", loci_gff3.name)
    else:
        _require_tool(mikado_bin, "mikado")
        loci_gff3 = run_pick(
            inputs.configuration_yaml,
            inputs.scoring_file,
            prepared_gtf,
            out_dir,
            procs=procs,
            mikado_bin=mikado_bin,
        )

    return MikadoRunResult(
        prepared_gtf=prepared_gtf,
        prepared_fasta=prepared_fasta,
        orfs_bed=orfs_bed,
        blast_out=blast_out,
        mikado_db=mikado_db,
        loci_gff3=loci_gff3,
        loci_metrics_tsv=out_dir / "mikado.loci.metrics.tsv",
        loci_scores_tsv=out_dir / "mikado.loci.scores.tsv",
    )
