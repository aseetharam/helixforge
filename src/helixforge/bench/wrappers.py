"""Benchmarking tool wrappers."""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    import pandas as pd

_log = get_logger(__name__)


class BenchmarkError(RuntimeError):
    """A benchmarking tool was missing or exited non-zero.

    ``kind`` distinguishes the two failure modes so :func:`benchmark_all` can
    report ``tool-absent`` vs ``failed`` in its status column
    instead of silently dropping a tool's rows with no trace of *why*.
    """

    def __init__(self, *args: Any, kind: str = "failed") -> None:
        super().__init__(*args)
        self.kind = kind


# ---------------------------------------------------------------------------
# Subprocess runner
# ---------------------------------------------------------------------------


def _run_tool(
    argv: list[Any], step: str, cwd: Path | str | None = None
) -> subprocess.CompletedProcess[str]:
    """Run ``argv`` (list) with logging + timing; raise :class:`BenchmarkError`.

    A missing binary (``FileNotFoundError``) and a non-zero exit both surface as
    a ``BenchmarkError`` with an actionable message — ``benchmark_all`` catches
    these so one absent tool never sinks the whole table.
    """
    argv = [str(a) for a in argv]
    _log.info("[%s] start: %s", step, " ".join(argv))
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            argv,
            check=True,
            capture_output=True,
            text=True,
            cwd=(str(cwd) if cwd is not None else None),
        )
    except FileNotFoundError as exc:
        raise BenchmarkError(
            f"{step}: tool not found ({argv[0]!r}). Is it installed and on PATH?",
            kind="tool-absent",
        ) from exc
    except subprocess.CalledProcessError as exc:
        stderr_tail = (exc.stderr or "")[-2000:]
        raise BenchmarkError(
            f"{step} failed (exit {exc.returncode}).\nargv: {argv}\n"
            f"stderr tail:\n{stderr_tail}",
            kind="failed",
        ) from exc
    _log.info("[%s] done in %.2fs", step, time.perf_counter() - t0)
    return proc


def _slug(text: str) -> str:
    """Lower-case, collapse non-alphanumerics to ``_`` for stable dict keys."""
    return re.sub(r"[^a-z0-9]+", "_", text.strip().lower()).strip("_")


def _glob_one(root: Path | str, pattern: str, what: str) -> Path:
    """First file matching ``pattern`` under ``root`` (recursive); error if none."""
    hits = sorted(Path(root).rglob(pattern))
    if not hits:
        raise BenchmarkError(f"{what}: no {pattern!r} found under {root}")
    return hits[0]


# ---------------------------------------------------------------------------
# Parsers — tolerant, line-oriented; take a path, read synthetic stats
# ---------------------------------------------------------------------------

# "<label> level...: Sn  Pr  F1"  (mikado compare .stats)
_MIKADO_LEVEL_RE = re.compile(
    r"^(?P<label>.*?level[^:]*):\s*(?P<sn>[\d.]+)\s+(?P<pr>[\d.]+)\s+(?P<f1>[\d.]+)",
    re.IGNORECASE,
)
# "<label> level: Sn | Pr |"  (gffcompare .stats)
_GFFCOMPARE_LEVEL_RE = re.compile(
    r"^\s*(?P<label>.*?level):\s*(?P<sn>[\d.]+)\s*\|\s*(?P<pr>[\d.]+)",
    re.IGNORECASE,
)
# compleasm summary.txt: "S:99.22%, 253" / "N:255"
_COMPLEASM_RE = re.compile(r"^\s*(?P<k>[SDFIMN]):(?P<v>[\d.]+)%?(?:,\s*(?P<n>\d+))?")
# BUSCO short_summary: "C:97.8%[S:96.5%,D:1.3%],F:0.7%,M:1.5%,n:255"
_BUSCO_RE = re.compile(
    r"C:(?P<C>[\d.]+)%\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],"
    r"F:(?P<F>[\d.]+)%,M:(?P<M>[\d.]+)%,n:(?P<n>\d+)"
)
# OMArk .sum: tolerant of both the flat "Single:90.00% (13042)" form and the
# real bracketed form "Single:[S:90.40%, 13102]" / "Consistent:[C:85.00%, ...]"
# (the percentage may sit after a "[X:" prefix). Captures the label before the
# first colon and the first percentage anywhere on the line.
_OMARK_RE = re.compile(r"^\s*(?P<label>[A-Za-z][A-Za-z ,]*?):[^\d%]*(?P<v>[\d.]+)%")
# AGAT agat_sp_statistics: "Number of gene    30000" / "mean cds length (bp)  1234.5"
_AGAT_RE = re.compile(r"^(?P<label>[A-Za-z][\w()/ '%.-]*?\S)\s{2,}(?P<v>-?[\d.]+)\s*$")

_COMPLEASM_KEYS = {
    "S": "single",
    "D": "duplicated",
    "F": "fragmented",
    "I": "incomplete",
    "M": "missing",
}


def parse_mikado_compare_stats(path: Path | str) -> dict[str, dict[str, float]]:
    """Parse a ``mikado compare`` ``.stats`` → ``{level: {sn, pr, f1}}``.

    Keys are slugged levels (``base``, ``exon_stringent``, ``intron``,
    ``intron_chain``, ``transcript``, ``gene`` …).
    """
    out: dict[str, dict[str, float]] = {}
    for line in Path(path).read_text().splitlines():
        m = _MIKADO_LEVEL_RE.search(line)
        if m:
            key = _slug(m.group("label").lower().replace("level", ""))
            out[key] = {
                "sn": float(m.group("sn")),
                "pr": float(m.group("pr")),
                "f1": float(m.group("f1")),
            }
    return out


def parse_gffcompare_stats(path: Path | str) -> dict[str, dict[str, float]]:
    """Parse a gffcompare ``.stats`` → ``{level: {sn, pr}}`` (transcript-level Sn/Pr)."""
    out: dict[str, dict[str, float]] = {}
    for line in Path(path).read_text().splitlines():
        m = _GFFCOMPARE_LEVEL_RE.search(line)
        if m:
            key = _slug(m.group("label").lower().replace("level", ""))
            out[key] = {"sn": float(m.group("sn")), "pr": float(m.group("pr"))}
    return out


def parse_compleasm_summary(path: Path | str) -> dict[str, float | int]:
    """Parse a compleasm ``summary.txt`` → completeness dict (S/D/F/M %, n, complete)."""
    out: dict[str, float | int] = {}
    for line in Path(path).read_text().splitlines():
        m = _COMPLEASM_RE.match(line)
        if not m:
            continue
        if m.group("k") == "N":
            out["n"] = int(float(m.group("v")))
        else:
            out[_COMPLEASM_KEYS[m.group("k")]] = float(m.group("v"))
    if "single" in out and "duplicated" in out:
        out["complete"] = float(out["single"]) + float(out["duplicated"])
    return out


def parse_busco_summary(path: Path | str) -> dict[str, float | int]:
    """Parse a BUSCO ``short_summary*.txt`` → completeness dict (complete=C, S/D/F/M, n)."""
    m = _BUSCO_RE.search(Path(path).read_text())
    if not m:
        return {}
    return {
        "complete": float(m.group("C")),
        "single": float(m.group("S")),
        "duplicated": float(m.group("D")),
        "fragmented": float(m.group("F")),
        "missing": float(m.group("M")),
        "n": int(m.group("n")),
    }


def parse_omark_summary(path: Path | str) -> dict[str, float]:
    """Parse an OMArk ``.sum`` → percentages keyed by slug (+ ``complete``).

    Captures both the completeness axis (single/duplicated/missing) and the
    consistency axis (consistent/inconsistent/contaminants/unknown).
    """
    out: dict[str, float] = {}
    for line in Path(path).read_text().splitlines():
        m = _OMARK_RE.match(line)
        if m:
            out[_slug(m.group("label"))] = float(m.group("v"))
    if "single" in out and "duplicated" in out:
        out["complete"] = out["single"] + out["duplicated"]
    return out


def parse_agat_stats(path: Path | str) -> dict[str, float]:
    """Parse ``agat_sp_statistics`` output → ``{slug: float}`` of structural counts."""
    out: dict[str, float] = {}
    for line in Path(path).read_text().splitlines():
        m = _AGAT_RE.match(line)
        if m:
            out[_slug(m.group("label"))] = float(m.group("v"))
    return out


# ---------------------------------------------------------------------------
# Run wrappers — build argv, invoke, parse the output file
# ---------------------------------------------------------------------------


def run_mikado_compare(
    reference_gff3: str | Path,
    prediction_gff3: str | Path,
    out_prefix: str | Path,
    mikado_bin: str = "mikado",
) -> dict[str, dict[str, float]]:
    """``mikado compare -r ref -p pred -o out_prefix``; return the parsed ``.stats`` dict.

    (The ``.stats`` file is written next to ``out_prefix`` as ``<out_prefix>.stats``;
    the parsed Sn/Pr/F1-per-level dict is the useful return value.)
    """
    _run_tool(
        [
            mikado_bin,
            "compare",
            "-r",
            reference_gff3,
            "-p",
            prediction_gff3,
            "-o",
            out_prefix,
        ],
        "mikado.compare",
    )
    return parse_mikado_compare_stats(Path(f"{out_prefix}.stats"))


def run_gffcompare(
    reference_gff3: str | Path,
    prediction_gtf: str | Path,
    out_prefix: str | Path,
    gffcompare_bin: str = "gffcompare",
) -> dict[str, dict[str, float]]:
    """``gffcompare -r ref -o out_prefix pred``; return transcript-level Sn/Pr dict."""
    _run_tool(
        [gffcompare_bin, "-r", reference_gff3, "-o", out_prefix, prediction_gtf],
        "gffcompare",
    )
    return parse_gffcompare_stats(Path(f"{out_prefix}.stats"))


def run_compleasm(
    protein_or_genome: str | Path,
    lineage: str,
    out_dir: str | Path,
    threads: int = 4,
    compleasm_bin: str = "compleasm",
) -> dict[str, float | int]:
    """``compleasm protein`` on a protein set; return parsed ``summary.txt`` dict."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _run_tool(
        [
            compleasm_bin,
            "protein",
            "-p",
            protein_or_genome,
            "-l",
            lineage,
            "-o",
            out_dir,
            "-t",
            threads,
        ],
        "compleasm",
    )
    return parse_compleasm_summary(out_dir / "summary.txt")


def run_busco(
    protein_fa: str | Path,
    lineage: str,
    out_dir: str | Path,
    threads: int = 4,
    mode: str = "proteins",
    busco_bin: str = "busco",
) -> dict[str, float | int]:
    """``busco -m proteins`` on a protein set; return parsed short-summary dict."""
    out_dir = Path(out_dir)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    _run_tool(
        [
            busco_bin,
            "-i",
            protein_fa,
            "-l",
            lineage,
            "-m",
            mode,
            "-o",
            out_dir.name,
            "--out_path",
            out_dir.parent,
            "-c",
            threads,
            "-f",
        ],
        "busco",
    )
    return parse_busco_summary(_glob_one(out_dir, "short_summary*.txt", "busco"))


def run_omark(
    proteins_fa: str | Path,
    omadb: str | Path,
    out_dir: str | Path,
    threads: int = 4,
    omark_bin: str = "omark",
) -> dict[str, float]:
    """``omark -f proteins -d omadb -o out_dir``; return completeness/consistency dict."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _run_tool(
        [omark_bin, "-f", proteins_fa, "-d", omadb, "-o", out_dir],
        "omark",
    )
    return parse_omark_summary(_glob_one(out_dir, "*.sum", "omark"))


def run_agat_stats(
    gff3: str | Path,
    out_path: str | Path,
    agat_bin: str = "agat_sp_statistics.pl",
) -> dict[str, float]:
    """``agat_sp_statistics.pl --gff gff3 -o out_path``; return the parsed counts dict."""
    _run_tool([agat_bin, "--gff", gff3, "-o", out_path], "agat.statistics")
    return parse_agat_stats(out_path)


# ---------------------------------------------------------------------------
# benchmark_all — run the applicable tools → one comparable table
# ---------------------------------------------------------------------------


def _flatten(prefix: str, value: Any, out: dict[str, Any]) -> None:
    """Flatten nested ``{level: {sn: ..}}`` into ``out[f'{prefix}_{k}'] = v``."""
    if isinstance(value, dict):
        for k, v in value.items():
            _flatten(f"{prefix}_{k}" if prefix else str(k), v, out)
    else:
        out[prefix] = value


def benchmark_all(
    reconciled_gff3: str | Path,
    proteins_fa: str | Path | None,
    out_dir: str | Path,
    reference_gff3: str | Path | None = None,
    lineage: str | None = None,
    omadb: str | Path | None = None,
    threads: int = 4,
    tools: list[str] | None = None,
    bins: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Run the applicable benchmarking tools → one tidy ``pandas.DataFrame``.

    Columns: ``tool, metric, value`` (long form — easy to pivot for the ablation
    figure or merge across runs). Which tools run depends on the inputs:

    - ``agat`` always (structural sanity on ``reconciled_gff3``);
    - ``mikado_compare`` + ``gffcompare`` when ``reference_gff3`` is given;
    - ``compleasm`` + ``busco`` when ``lineage`` and ``proteins_fa`` are given;
    - ``omark`` when ``omadb`` and ``proteins_fa`` are given.

    ``tools`` optionally restricts to a subset of tool names. ``bins`` overrides
    individual binary paths (e.g. ``{"mikado": "/opt/mikado"}``). Every row carries
    a ``status`` column: ``ok`` for a metric a tool produced,
    and a single synthetic ``status``-metric row per attempted tool that failed —
    ``tool-absent`` when the binary was not found, ``failed`` on a non-zero exit —
    so the caller can tell "tool absent" from "metric genuinely zero" instead of
    seeing the tool silently vanish. This fills the BUSCO/compleasm/OMArk
    completeness stubs left in ``stats/before_after`` (Phase 9).
    """
    import pandas as pd

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bins = bins or {}
    rows: list[dict[str, Any]] = []

    def _want(name: str) -> bool:
        return tools is None or name in tools

    def _collect(tool: str, result: Any) -> None:
        flat: dict[str, Any] = {}
        _flatten("", result, flat)
        for metric, value in flat.items():
            rows.append(
                {"tool": tool, "metric": metric, "value": value, "status": "ok"}
            )

    def _try(name: str, thunk: Any) -> None:
        if not _want(name):
            return
        try:
            _collect(name, thunk())
        except BenchmarkError as exc:
            _log.warning("benchmark tool %r skipped (%s): %s", name, exc.kind, exc)
            rows.append(
                {
                    "tool": name,
                    "metric": "status",
                    "value": float("nan"),
                    "status": exc.kind,
                }
            )

    _try(
        "agat",
        lambda: run_agat_stats(
            reconciled_gff3,
            out_dir / "agat_stats.txt",
            agat_bin=bins.get("agat", "agat_sp_statistics.pl"),
        ),
    )

    if reference_gff3:
        _try(
            "mikado_compare",
            lambda: run_mikado_compare(
                reference_gff3,
                reconciled_gff3,
                str(out_dir / "mikado_compare"),
                mikado_bin=bins.get("mikado", "mikado"),
            ),
        )
        _try(
            "gffcompare",
            lambda: run_gffcompare(
                reference_gff3,
                reconciled_gff3,
                str(out_dir / "gffcompare"),
                gffcompare_bin=bins.get("gffcompare", "gffcompare"),
            ),
        )

    if lineage and proteins_fa:
        _try(
            "compleasm",
            lambda: run_compleasm(
                proteins_fa,
                lineage,
                out_dir / "compleasm",
                threads=threads,
                compleasm_bin=bins.get("compleasm", "compleasm"),
            ),
        )
        _try(
            "busco",
            lambda: run_busco(
                proteins_fa,
                lineage,
                out_dir / "busco",
                threads=threads,
                busco_bin=bins.get("busco", "busco"),
            ),
        )

    if omadb and proteins_fa:
        _try(
            "omark",
            lambda: run_omark(
                proteins_fa,
                omadb,
                out_dir / "omark",
                threads=threads,
                omark_bin=bins.get("omark", "omark"),
            ),
        )

    return pd.DataFrame(rows, columns=["tool", "metric", "value", "status"])
