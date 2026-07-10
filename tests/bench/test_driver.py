"""Phase 17 D4 — benchmark doc renderer (pure)."""

import math

import pandas as pd

from helixforge.bench.driver import write_benchmark_md


def test_write_benchmark_md_renders_all_sections(tmp_path):
    before_after = pd.DataFrame([
        {"metric": "Genes", "helixer": 100, "helixforge": 110, "delta": 10},
        {"metric": "Mean Helixer support", "helixer": 0.80, "helixforge": 0.85, "delta": 0.05},
    ])
    benchmark = pd.DataFrame([{"tool": "agat", "metric": "number_of_gene", "value": 110.0}])
    ablation = pd.DataFrame([
        {"variant": "full", "num_genes": 110},
        {"variant": "no_helixer_support", "num_genes": 108},
    ])
    out = tmp_path / "BENCHMARK.md"
    write_benchmark_md(out, before_after=before_after, benchmark=benchmark, ablation=ablation)
    text = out.read_text()
    assert "# HelixForge v3, Araport11 benchmark" in text
    assert "Before / after" in text and "Genes" in text
    assert "External benchmarks" in text and "number_of_gene" in text
    assert "Ablations" in text and "no_helixer_support" in text
    assert "| 0.85 |" in text  # helixforge support rendered


def test_write_benchmark_md_marks_absent_sections(tmp_path):
    out = tmp_path / "BENCHMARK.md"
    write_benchmark_md(out)
    text = out.read_text()
    # No fabricated numbers — honest "not run" markers (phase: no fabrication).
    assert text.count("_not run") >= 2


def test_write_benchmark_md_nan_renders_dash(tmp_path):
    before_after = pd.DataFrame([
        {"metric": "Mean Helixer support", "helixer": float("nan"),
         "helixforge": float("nan"), "delta": float("nan")},
    ])
    out = tmp_path / "BENCHMARK.md"
    write_benchmark_md(out, before_after=before_after)
    line = [ln for ln in out.read_text().splitlines() if "Mean Helixer support" in ln][0]
    assert "—" in line
    assert not any(c.isdigit() for c in line)
