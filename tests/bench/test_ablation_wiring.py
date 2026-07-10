"""Phase 17 D3 — ablation runner wiring + shared-Mikado reuse.

Pipeline + benchmark are mocked. Asserts each variant toggles exactly its field,
the table has one row per variant, and pick-only variants reuse the baseline's
Mikado artifacts while stage-changing variants re-run Mikado in full. Floor: 5.
"""

from pathlib import Path

import pandas as pd
import pytest

from helixforge.bench import ablation
from helixforge.bench.ablation import (
    PICK_ONLY_KNOBS,
    _is_pick_only,
    _order_variants,
    run_ablation,
)
from helixforge.reconcile.pipeline import PipelineConfig


@pytest.fixture
def base_config():
    return PipelineConfig(genome_fasta="/nope/genome.fa", helixer_gff3="/nope/helixer.gff3",
                          output_prefix="run")


@pytest.fixture
def recorder(monkeypatch):
    seen = []

    def fake_run_pipeline(config):
        seen.append(config)
        return ["gene_a", "gene_b"]

    def fake_benchmark_all(reconciled_gff3, proteins_fa, out_dir, **kw):
        return pd.DataFrame([{"tool": "agat", "metric": "number_of_gene", "value": 2.0}])

    monkeypatch.setattr(ablation, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(ablation, "benchmark_all", fake_benchmark_all)
    return seen


def _by_variant(recorder):
    """Map a recorded config back to its variant via the output-prefix suffix."""
    out = {}
    for cfg in recorder:
        name = Path(cfg.output_prefix).name  # e.g. "run_no_pad"
        out[name.split("run_", 1)[-1]] = cfg
    return out


# ---------------------------------------------------------------------------
# Stage classification
# ---------------------------------------------------------------------------

def test_pick_only_classification():
    assert _is_pick_only("no_pad") is True
    assert _is_pick_only("strict_vs_permissive") is True
    assert _is_pick_only("no_helixer_support") is False  # external metric → serialise
    assert _is_pick_only("no_reference_flag") is False   # prepare stage
    assert _is_pick_only("full") is False                # no override


def test_pick_only_knobs_membership():
    assert "pad" in PICK_ONLY_KNOBS
    assert "scoring_profile" in PICK_ONLY_KNOBS
    assert "helixer_support_weight" not in PICK_ONLY_KNOBS


def test_order_variants_runs_full_first():
    assert _order_variants(["no_pad", "full", "no_helixer_support"])[0] == "full"
    # full absent → order preserved.
    assert _order_variants(["no_pad", "no_helixer_support"]) == ["no_pad", "no_helixer_support"]


# ---------------------------------------------------------------------------
# run_ablation wiring
# ---------------------------------------------------------------------------

def test_one_row_per_variant_preserves_order(base_config, tmp_path, recorder):
    variants = ["no_pad", "full", "no_helixer_support"]
    df = run_ablation(base_config, variants, tmp_path)
    # Output preserves the caller's order even though full runs first internally.
    assert list(df["variant"]) == variants
    assert len(df) == 3
    assert len(recorder) == 3


def test_pick_only_variant_reuses_baseline_mikado(base_config, tmp_path, recorder):
    run_ablation(base_config, ["full", "no_pad", "no_helixer_support"], tmp_path)
    cfgs = _by_variant(recorder)
    # no_pad is pick-only → reuse the full baseline's mikado_run dir.
    full_mikado = str(Path(cfgs["full"].work_dir) / "mikado_run")
    assert cfgs["no_pad"].reuse_mikado_dir == full_mikado
    # full itself never reuses.
    assert cfgs["full"].reuse_mikado_dir is None


def test_stage_changing_variant_does_not_reuse(base_config, tmp_path, recorder):
    run_ablation(base_config, ["full", "no_helixer_support", "no_reference_flag"], tmp_path)
    cfgs = _by_variant(recorder)
    assert cfgs["no_helixer_support"].reuse_mikado_dir is None
    assert cfgs["no_reference_flag"].reuse_mikado_dir is None


def test_each_variant_toggles_exactly_its_field(base_config, tmp_path, recorder):
    run_ablation(base_config, ["full", "no_pad", "strict_vs_permissive"], tmp_path)
    cfgs = _by_variant(recorder)
    assert cfgs["no_pad"].pad is False
    assert cfgs["no_pad"].scoring_profile == "strict"      # only pad moved
    assert cfgs["strict_vs_permissive"].scoring_profile == "permissive"
    assert cfgs["strict_vs_permissive"].pad is True        # only profile moved


def test_no_reuse_when_full_absent(base_config, tmp_path, recorder):
    # Without a baseline, pick-only variants cannot reuse anything.
    run_ablation(base_config, ["no_pad"], tmp_path)
    cfgs = _by_variant(recorder)
    assert cfgs["no_pad"].reuse_mikado_dir is None
