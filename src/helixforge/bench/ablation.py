"""Ablation engine — the manuscript figure."""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import TYPE_CHECKING, Any

from helixforge.bench.wrappers import benchmark_all
from helixforge.export.writers import write_protein_fasta
from helixforge.io.fasta import GenomeAccessor
from helixforge.reconcile.pipeline import run_pipeline
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    import pandas as pd
    from helixforge.reconcile.pipeline import PipelineConfig

_log = get_logger(__name__)


# variant name → PipelineConfig field overrides (DESIGN §11). ``full`` is the
# baseline (no override). Output paths are re-derived per variant in
# ``_variant_config`` so reruns don't clobber each other.
ABLATION_VARIANTS: dict[str, dict[str, Any]] = {
    "full": {},
    "no_helixer_support": {"helixer_support_weight": 0.0},
    "no_reference_flag": {"helixer_is_reference": False},
    "no_pad": {"pad": False},
    "strict_vs_permissive": {"scoring_profile": "permissive"},
}

# Knobs that change ONLY the pick / alternative-splicing stage — the prepared
# transcripts and the serialise DB are identical to the baseline, so a variant
# toggling only these can reuse the baseline's Mikado artifacts and re-run pick
# alone (``PipelineConfig.reuse_mikado_dir``). ``helixer_support_weight`` (external
# metric → serialise) and ``helixer_is_reference`` (prepare) are NOT here: they
# change upstream stages and must re-run Mikado in full.
PICK_ONLY_KNOBS: frozenset[str] = frozenset(
    {
        "pad",
        "scoring_profile",
        "only_confirmed_introns",
        "keep_retained_introns",
        "max_isoforms",
    }
)


def _is_pick_only(variant: str) -> bool:
    """True if ``variant`` only toggles pick-stage knobs (reuse-eligible)."""
    overrides = ABLATION_VARIANTS.get(variant, {})
    return bool(overrides) and set(overrides) <= PICK_ONLY_KNOBS


def _order_variants(variants: list[str]) -> list[str]:
    """Run the full-Mikado baseline first so pick-only variants can reuse it."""
    variants = list(variants)
    if "full" in variants:
        return ["full"] + [v for v in variants if v != "full"]
    return variants


def _variant_config(
    config: PipelineConfig,
    variant: str,
    out_dir: Path | str,
    reuse_mikado_dir: str | None = None,
) -> PipelineConfig:
    """Clone ``config`` with the variant's lever toggled and per-variant outputs.

    ``reuse_mikado_dir`` (set by :func:`run_ablation` for pick-only variants)
    points the variant at a baseline's ``mikado_run`` dir so only pick reruns.
    """
    if variant not in ABLATION_VARIANTS:
        raise ValueError(
            f"unknown ablation variant {variant!r}; "
            f"known variants: {sorted(ABLATION_VARIANTS)}"
        )
    vdir = Path(out_dir) / variant
    overrides: dict[str, Any] = dict(ABLATION_VARIANTS[variant])
    # Re-derive output/report/id_map/work paths under the variant dir: setting
    # them to None lets PipelineConfig.__post_init__ rebuild from the new prefix.
    overrides["output_prefix"] = str(
        vdir / f"{Path(config.output_prefix).name}_{variant}"
    )
    overrides["report_path"] = None
    overrides["id_map_path"] = None
    overrides["work_dir"] = None
    if reuse_mikado_dir is not None:
        overrides["reuse_mikado_dir"] = str(reuse_mikado_dir)
    return dataclasses.replace(config, **overrides)


def _export_proteins(
    genes: list[Any],
    config: PipelineConfig,
    out_dir: Path | str,
    variant: str,
) -> str | None:
    """Write the variant's protein FASTA (for compleasm/BUSCO/OMArk); ``None`` if no genome."""
    if not Path(config.genome_fasta).exists():
        return None
    vdir = Path(out_dir) / variant
    vdir.mkdir(parents=True, exist_ok=True)
    proteins = vdir / f"{variant}.proteins.fa"
    genome = GenomeAccessor(config.genome_fasta)
    try:
        write_protein_fasta(genes, genome, str(proteins))
    finally:
        genome.close()
    return str(proteins)


def run_ablation(
    config: PipelineConfig,
    variants: list[str],
    out_dir: Path | str,
    reference_gff3: str | None = None,
    lineage: str | None = None,
    omadb: str | None = None,
    threads: int = 4,
) -> pd.DataFrame:
    """Rerun the pipeline + benchmark under each ``variant``; return a DataFrame.

    One row per variant; columns are ``variant``, ``num_genes`` and the flattened
    benchmark metrics (``<tool>.<metric>``) from :func:`benchmark_all`. The
    ``reference_gff3`` / ``lineage`` / ``omadb`` benchmarking inputs are shared
    across variants so the rows are directly comparable.

    To keep the ablation tractable, the full-Mikado baseline (``full``) runs
    first; pick-only variants (``no_pad``, ``strict_vs_permissive`` — see
    :data:`PICK_ONLY_KNOBS`) then **reuse** its prepared transcripts + serialise
    DB and re-run only ``mikado pick``. Variants that change an upstream stage
    (``no_helixer_support``, ``no_reference_flag``) re-run Mikado in full. Output
    rows preserve the caller's ``variants`` order regardless of run order.
    """
    import pandas as pd

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    shared_mikado_dir: str | None = None
    results: dict[str, dict[str, Any]] = {}
    for variant in _order_variants(variants):
        reuse = shared_mikado_dir if _is_pick_only(variant) else None
        vcfg = _variant_config(config, variant, out_dir, reuse_mikado_dir=reuse)
        _log.info(
            "ablation[%s]: running pipeline%s",
            variant,
            " (reusing baseline Mikado)" if reuse else "",
        )
        genes = run_pipeline(vcfg)

        # The baseline's Mikado artifacts seed reuse for later pick-only variants.
        if variant == "full":
            assert vcfg.work_dir is not None
            shared_mikado_dir = str(Path(vcfg.work_dir) / "mikado_run")

        proteins = _export_proteins(genes, vcfg, out_dir, variant)
        bench_df = benchmark_all(
            f"{vcfg.output_prefix}.gff3",
            proteins,
            str(Path(out_dir) / variant / "bench"),
            reference_gff3=reference_gff3,
            lineage=lineage,
            omadb=omadb,
            threads=threads,
        )

        row: dict[str, Any] = {"variant": variant, "num_genes": len(genes)}
        for _, r in bench_df.iterrows():
            row[f"{r['tool']}.{r['metric']}"] = r["value"]
        results[variant] = row

    return pd.DataFrame([results[v] for v in variants])
