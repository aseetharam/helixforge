"""End-to-end benchmark driver — the M7 validation harness."""

from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, Any

from helixforge.bench.isoform import (
    isoform_accuracy,
    isoform_result_to_rows,
    match_genes_1to1,
    matched_locus_isoform_accuracy,
)
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    import pandas as pd
    from helixforge.reconcile.pipeline import PipelineConfig
    from helixforge.reconcile.models import ReconciledGene

_log = get_logger(__name__)

__all__ = [
    "isoform_accuracy",
    "isoform_result_to_rows",
    "match_genes_1to1",
    "matched_locus_isoform_accuracy",
    "run_araport11_benchmark",
    "write_benchmark_md",
]


# ---------------------------------------------------------------------------
# Markdown rendering (pure)
# ---------------------------------------------------------------------------


def _fmt(v: Any) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        if math.isnan(v):
            return "—"
        return f"{v:.2f}" if v != int(v) else str(int(v))
    return str(v)


def _df_to_md(df: pd.DataFrame, columns: list[str] | None = None) -> str:
    """Render a ``pandas.DataFrame`` as a GitHub markdown table."""
    cols = list(columns or df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(_fmt(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_benchmark_md(
    out_path: str | Path,
    *,
    before_after: pd.DataFrame | None = None,
    benchmark: pd.DataFrame | None = None,
    ablation: pd.DataFrame | None = None,
    title: str = "HelixForge v3 — Araport11 benchmark",
    notes: str | None = None,
) -> str:
    """Render the benchmark doc from the (optional) result tables.

    Each section is emitted only when its table is supplied; absent sections show
    a "not run" note so the doc is honest about what was produced. Returns the
    written path.
    """
    parts: list[str] = [f"# {title}", ""]
    if notes:
        parts += [notes, ""]

    parts += ["## Before / after (Helixer → HelixForge)", ""]
    if before_after is not None and len(before_after):
        parts.append(
            _df_to_md(before_after, ["metric", "helixer", "helixforge", "delta"])
        )
    else:
        parts.append("_not run_")
    parts.append("")

    parts += ["## External benchmarks", ""]
    if benchmark is not None and len(benchmark):
        parts.append(_df_to_md(benchmark, ["tool", "metric", "value"]))
    else:
        parts.append("_not run (needs reference GFF3 / lineage / OMA db)_")
    parts.append("")

    parts += [
        "## Ablations (§11 — Helixer-coupling levers)",
        "",
        "The `no_helixer_support` row vs `full` isolates the project's novel "
        "contribution (the Helixer↔evidence external metric).",
        "",
    ]
    if ablation is not None and len(ablation):
        parts.append(_df_to_md(ablation))
    else:
        parts.append("_not run_")
    parts.append("")

    Path(out_path).write_text("\n".join(parts))
    _log.info("wrote benchmark doc → %s", out_path)
    return str(out_path)


# ---------------------------------------------------------------------------
# Orchestration (runs real tools; executed out of band for M7)
# ---------------------------------------------------------------------------


def run_araport11_benchmark(
    config: PipelineConfig,
    reference_gff3: str,
    proteins_fa: str,
    out_dir: str | Path,
    *,
    reconciled_gff3: str | None = None,
    genes: list[ReconciledGene] | None = None,
    junctions: Any = None,
    lineage: str | None = None,
    omadb: str | None = None,
    variants: tuple[str, ...] = (
        "full",
        "no_helixer_support",
        "no_reference_flag",
        "no_pad",
        "strict_vs_permissive",
    ),
    threads: int = 4,
    doc_path: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run the full validation harness → result tables + ``BENCHMARK.md``.

    ``config`` is a ``PipelineConfig`` (drives reconcile + the ablations).
    ``reconciled_gff3``/``genes`` may be supplied to skip re-running the pipeline
    for the before/after table; otherwise the pipeline is run once. The before/
    after Helixer-support column is driven by ``config.helixer_h5`` and the
    junction set by ``junctions``. Returns ``(before_after_df, benchmark_df,
    ablation_df)`` and writes the doc.
    """
    import pandas as pd

    from helixforge.bench.ablation import run_ablation
    from helixforge.bench.wrappers import BenchmarkError, benchmark_all
    from helixforge.reconcile.pipeline import run_pipeline
    from helixforge.stats.before_after import before_after_table

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    after: Any
    if genes is None and reconciled_gff3 is None:
        _log.info("driver: running reconcile for the before/after table")
        genes = run_pipeline(config)
        reconciled_gff3 = f"{config.output_prefix}.gff3"
    after = genes if genes is not None else reconciled_gff3

    before_after = before_after_table(
        config.helixer_gff3,
        after,
        h5_path=config.helixer_h5,
        junctions=junctions,
    )
    pred_gff3 = reconciled_gff3 or f"{config.output_prefix}.gff3"
    benchmark = benchmark_all(
        pred_gff3,
        proteins_fa,
        str(out_dir / "bench"),
        reference_gff3=reference_gff3,
        lineage=lineage,
        omadb=omadb,
        threads=threads,
    )

    # Isoform-level accuracy (the poster's "better isoforms" number). Needs a
    # reference; fault-tolerant like benchmark_all so an absent ``mikado`` never
    # sinks the table (a single ``status`` row records why instead).
    if reference_gff3:
        iso_dir = out_dir / "isoform"
        iso_dir.mkdir(parents=True, exist_ok=True)
        try:
            iso = isoform_accuracy(
                pred_gff3,
                reference_gff3,
                str(iso_dir / "iso"),
                helixer_gff3=config.helixer_gff3,
                proteins_fa=proteins_fa,
                omadb=omadb,
            )
            iso_rows = isoform_result_to_rows(iso)
        except BenchmarkError as exc:
            _log.warning("isoform metric skipped (%s): %s", exc.kind, exc)
            iso_rows = [
                {
                    "tool": "isoform",
                    "metric": "status",
                    "value": float("nan"),
                    "status": exc.kind,
                }
            ]
        benchmark = pd.concat(
            [
                benchmark,
                pd.DataFrame(iso_rows, columns=["tool", "metric", "value", "status"]),
            ],
            ignore_index=True,
        )

    ablation = run_ablation(
        config,
        list(variants),
        str(out_dir / "ablation"),
        reference_gff3=reference_gff3,
        lineage=lineage,
        omadb=omadb,
        threads=threads,
    )

    doc_path = doc_path or str(out_dir / "BENCHMARK.md")
    write_benchmark_md(
        doc_path,
        before_after=before_after,
        benchmark=benchmark,
        ablation=ablation,
    )
    return before_after, benchmark, ablation
