"""Phase 11 D2 — ablation engine.

The pipeline and benchmark are mocked; we assert each variant toggles exactly
the right ``PipelineConfig`` lever and that the table has one row per variant.
"""

import pandas as pd
import pytest

from helixforge.bench import ablation
from helixforge.bench.ablation import ABLATION_VARIANTS, _variant_config, run_ablation
from helixforge.reconcile.pipeline import PipelineConfig


@pytest.fixture
def base_config():
    # Fake paths: genome does not exist → _export_proteins returns None (no I/O).
    return PipelineConfig(genome_fasta="/nope/genome.fa", helixer_gff3="/nope/helixer.gff3",
                          output_prefix="run")


@pytest.fixture
def recorder(monkeypatch):
    """Patch run_pipeline (records configs) + benchmark_all (tiny table)."""
    seen = []

    def fake_run_pipeline(config):
        seen.append(config)
        return ["gene_a", "gene_b"]

    def fake_benchmark_all(reconciled_gff3, proteins_fa, out_dir, **kw):
        return pd.DataFrame(
            [{"tool": "agat", "metric": "number_of_gene", "value": 2.0}]
        )

    monkeypatch.setattr(ablation, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(ablation, "benchmark_all", fake_benchmark_all)
    return seen


def _config_for(seen, variant_index):
    return seen[variant_index]


def test_variant_full_changes_no_lever(base_config, tmp_path):
    cfg = _variant_config(base_config, "full", tmp_path)
    assert cfg.pad is True
    assert cfg.helixer_is_reference is True
    assert cfg.helixer_support_weight == 1.0
    assert cfg.scoring_profile == "strict"


def test_variant_no_helixer_support(base_config, tmp_path):
    cfg = _variant_config(base_config, "no_helixer_support", tmp_path)
    assert cfg.helixer_support_weight == 0.0
    assert cfg.helixer_is_reference is True  # only this lever moved


def test_variant_no_reference_flag(base_config, tmp_path):
    cfg = _variant_config(base_config, "no_reference_flag", tmp_path)
    assert cfg.helixer_is_reference is False
    assert cfg.pad is True


def test_variant_no_pad(base_config, tmp_path):
    cfg = _variant_config(base_config, "no_pad", tmp_path)
    assert cfg.pad is False


def test_variant_strict_vs_permissive(base_config, tmp_path):
    cfg = _variant_config(base_config, "strict_vs_permissive", tmp_path)
    assert cfg.scoring_profile == "permissive"


def test_variant_outputs_are_per_variant(base_config, tmp_path):
    a = _variant_config(base_config, "full", tmp_path)
    b = _variant_config(base_config, "no_pad", tmp_path)
    assert a.output_prefix != b.output_prefix
    assert "full" in a.output_prefix and "no_pad" in b.output_prefix
    # report/id_map/work paths re-derived under the new prefix
    assert a.report_path.startswith(a.output_prefix)
    assert a.work_dir != b.work_dir


def test_unknown_variant_raises(base_config, tmp_path):
    with pytest.raises(ValueError, match="unknown ablation variant"):
        _variant_config(base_config, "bogus", tmp_path)


def test_run_ablation_one_row_per_variant(base_config, tmp_path, recorder):
    variants = ["full", "no_helixer_support", "no_pad"]
    df = run_ablation(base_config, variants, tmp_path)
    assert list(df["variant"]) == variants
    assert len(df) == 3
    assert len(recorder) == 3  # pipeline run once per variant


def test_run_ablation_toggles_levers_per_row(base_config, tmp_path, recorder):
    run_ablation(base_config, ["full", "no_helixer_support", "no_reference_flag"], tmp_path)
    assert recorder[0].helixer_support_weight == 1.0
    assert recorder[1].helixer_support_weight == 0.0
    assert recorder[2].helixer_is_reference is False


def test_run_ablation_table_has_benchmark_metrics(base_config, tmp_path, recorder):
    df = run_ablation(base_config, ["full"], tmp_path)
    assert "agat.number_of_gene" in df.columns
    assert df["num_genes"].iloc[0] == 2
    assert df["agat.number_of_gene"].iloc[0] == 2.0


def test_full_variant_set_runs(base_config, tmp_path, recorder):
    df = run_ablation(base_config, sorted(ABLATION_VARIANTS), tmp_path)
    assert len(df) == len(ABLATION_VARIANTS)
