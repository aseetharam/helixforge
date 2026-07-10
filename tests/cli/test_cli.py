"""Phase 11 D3 — click CLI.

The Python API is mocked; we assert each command parses its options and calls
the right API function with the expected arguments. ``--help`` works for every
command and subcommand.
"""

import pandas as pd
import pytest
from click.testing import CliRunner

from helixforge import cli


@pytest.fixture
def runner():
    return CliRunner()


def _touch_inputs():
    """Create the genome + helixer files the pipeline options require to exist."""
    from pathlib import Path

    Path("genome.fa").write_text(">chr1\nACGTACGTAC\n")
    Path("helixer.gff3").write_text("##gff-version 3\n")


# ---------------------------------------------------------------------------
# --help
# ---------------------------------------------------------------------------


def test_main_help(runner):
    result = runner.invoke(cli.main, ["--help"])
    assert result.exit_code == 0
    # Intent-grouped survivors are listed…
    for cmd in ("reconcile", "parallel", "confidence", "evidence", "stats",
                "doctor", "viz", "utils"):
        assert cmd in result.output
    # …and the hidden command is not listed (`run` removal is checked by
    # test_run_command_removed; the word "run" appears in prose here).
    assert "benchmark" not in result.output  # hidden (dev/paper tooling)


def test_main_help_is_intent_grouped(runner):
    result = runner.invoke(cli.main, ["--help"])
    assert result.exit_code == 0
    for section in ("Annotate:", "Score & inspect:", "Preflight & utilities:"):
        assert section in result.output


def test_command_depth_at_most_two(runner):
    """Cap nesting at two levels: top-level command + at most one subcommand."""
    import click

    def depth(cmd, d=0):
        if isinstance(cmd, click.Group):
            return max([depth(c, d + 1) for c in cmd.commands.values()] or [d])
        return d

    assert depth(cli.main) <= 2


def test_run_command_removed(runner):
    result = runner.invoke(cli.main, ["run", "--help"])
    assert result.exit_code != 0  # no such command


def test_benchmark_hidden_but_invocable(runner):
    # Hidden from the listing, yet --help (and execution) still work.
    result = runner.invoke(cli.main, ["benchmark", "--help"])
    assert result.exit_code == 0
    assert "all" in result.output and "ablation" in result.output


def test_reconcile_help(runner):
    result = runner.invoke(cli.main, ["reconcile", "--help"])
    assert result.exit_code == 0
    assert "--genome" in result.output
    assert "--scoring-profile" in result.output


def test_protein_db_help_mentions_mikado_gating(runner):
    result = runner.invoke(cli.main, ["reconcile", "--help"])
    assert result.exit_code == 0
    text = result.output.lower()
    idx = text.find("--protein-db file")
    assert idx >= 0
    snippet = text[idx : idx + 250]
    assert "isoform" in snippet
    assert "mikado" in snippet or "reconcil" in snippet


def test_miniprot_help_not_substitute(runner):
    result = runner.invoke(cli.main, ["reconcile", "--help"])
    assert result.exit_code == 0
    text = result.output.lower()
    idx = text.find("--miniprot")
    assert idx >= 0
    snippet = text[idx : idx + 200]
    assert "not" in snippet and "protein-db" in snippet


def test_stats_help(runner):
    result = runner.invoke(cli.main, ["stats", "--help"])
    assert result.exit_code == 0
    assert "--helixer" in result.output


def test_benchmark_help_lists_subcommands(runner):
    result = runner.invoke(cli.main, ["benchmark", "--help"])
    assert result.exit_code == 0
    assert "all" in result.output and "ablation" in result.output


# ---------------------------------------------------------------------------
# reconcile
# ---------------------------------------------------------------------------


def test_reconcile_calls_pipeline(runner, monkeypatch):
    captured = {}

    def fake_run_pipeline(config):
        captured["config"] = config
        return []

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--output-prefix", "myrun",
        ])
    assert result.exit_code == 0, result.output
    cfg = captured["config"]
    assert cfg.genome_fasta == "genome.fa"
    assert cfg.helixer_gff3 == "helixer.gff3"
    assert cfg.output_prefix == "myrun"


def test_reconcile_toggles_pad_and_profile(runner, monkeypatch):
    captured = {}

    def fake_run_pipeline(config):
        captured["config"] = config
        return []

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--no-pad", "--scoring-profile", "permissive",
            "--helixer-support-weight", "0",
        ])
    assert result.exit_code == 0, result.output
    cfg = captured["config"]
    assert cfg.pad is False
    assert cfg.scoring_profile == "permissive"
    assert cfg.helixer_support_weight == 0.0


def test_reconcile_functional_flags_map_to_config(runner, monkeypatch):
    # §3.9: --functional-* flags flow into PipelineConfig (default path off).
    captured = {}

    def fake_run_pipeline(config):
        captured["config"] = config
        return []

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--functional-annotation", "--functional-tool", "both",
            "--functional-db", "/data/eggnog",
            "--interproscan-bin", "/abs/interproscan.sh",
            "--eggnog-bin", "/abs/emapper.py",
        ])
    assert result.exit_code == 0, result.output
    cfg = captured["config"]
    assert cfg.functional_annotation is True
    assert cfg.functional_tool == "both"
    assert cfg.functional_db == "/data/eggnog"
    assert cfg.interproscan_bin == "/abs/interproscan.sh"
    assert cfg.eggnog_bin == "/abs/emapper.py"


def test_reconcile_helixer_h5_flags_map_to_config(runner, monkeypatch):
    # --helixer-h5 + --helixer-input-h5 flow into PipelineConfig so reconcile can
    # resolve Helixer's split predictions/input halves (parity with `confidence`).
    captured = {}

    def fake_run_pipeline(config):
        captured["config"] = config
        return []

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    with runner.isolated_filesystem():
        _touch_inputs()
        from pathlib import Path

        Path("preds_predictions.h5").write_text("x")
        Path("preds_input.h5").write_text("x")
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--helixer-h5", "preds_predictions.h5",
            "--helixer-input-h5", "preds_input.h5",
        ])
    assert result.exit_code == 0, result.output
    cfg = captured["config"]
    assert cfg.helixer_h5 == "preds_predictions.h5"
    assert cfg.helixer_input_h5 == "preds_input.h5"


def test_reconcile_functional_annotation_off_by_default(runner, monkeypatch):
    captured = {}

    def fake_run_pipeline(config):
        captured["config"] = config
        return []

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
        ])
    assert result.exit_code == 0, result.output
    # Default reproduces the golden path: hook off, default tool/binaries.
    assert captured["config"].functional_annotation is False
    assert captured["config"].functional_tool == "interproscan"
    assert captured["config"].interproscan_bin == "interproscan.sh"


def test_reconcile_functional_tool_rejects_unknown(runner, monkeypatch):
    monkeypatch.setattr(cli, "run_pipeline", lambda config: [])
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--functional-tool", "emapper",   # not a valid selector (use 'eggnog')
        ])
    assert result.exit_code != 0
    assert "emapper" in result.output.lower() or "invalid" in result.output.lower()


def test_reconcile_collects_repeatable_inputs(runner, monkeypatch):
    captured = {}

    def fake_run_pipeline(config):
        captured["config"] = config
        return []

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    with runner.isolated_filesystem():
        from pathlib import Path

        _touch_inputs()
        Path("a.gtf").write_text("")
        Path("b.gtf").write_text("")
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--stringtie", "a.gtf", "--stringtie", "b.gtf",
        ])
    assert result.exit_code == 0, result.output
    assert captured["config"].stringtie_list == ["a.gtf", "b.gtf"]


def test_reconcile_missing_genome_errors(runner, monkeypatch):
    monkeypatch.setattr(cli, "run_pipeline", lambda config: [])
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, ["reconcile", "--helixer", "helixer.gff3"])
    assert result.exit_code != 0
    assert "genome" in result.output.lower()


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------


def test_stats_calls_before_after(runner, monkeypatch):
    calls = {}

    def fake_table(helixer, reconciled, genome=None, h5_path=None, **kw):
        calls["table"] = (helixer, reconciled)
        calls["h5_path"] = h5_path
        return pd.DataFrame([{"metric": "Genes", "helixer": 1, "helixforge": 1, "delta": 0}])

    def fake_write(df, path):
        calls["write"] = path
        return path

    monkeypatch.setattr(cli, "before_after_table", fake_table)
    monkeypatch.setattr(cli, "write_summary", fake_write)
    with runner.isolated_filesystem():
        from pathlib import Path

        Path("helixer.gff3").write_text("##gff-version 3\n")
        Path("recon.gff3").write_text("##gff-version 3\n")
        result = runner.invoke(cli.main, [
            "stats", "--helixer", "helixer.gff3", "--reconciled", "recon.gff3",
            "--out", "summary.md",
        ])
    assert result.exit_code == 0, result.output
    assert calls["table"] == ("helixer.gff3", "recon.gff3")
    assert calls["write"] == "summary.md"


# ---------------------------------------------------------------------------
# viz
# ---------------------------------------------------------------------------


class _FakeGene:
    def __init__(self, gid):
        self.gene_id = gid


def _patch_viz_pipeline(monkeypatch, genes=None):
    monkeypatch.setattr(cli, "run_pipeline", lambda config: genes or [_FakeGene("HFG_00001")])


def test_viz_static_calls_plot_loci(runner, monkeypatch):
    calls = {}
    _patch_viz_pipeline(monkeypatch)
    monkeypatch.setattr(cli, "plot_loci",
                        lambda genes, out_dir, **kw: calls.setdefault("plot_loci", (genes, out_dir)))
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "viz", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--out-dir", "figs", "--mode", "static",
        ])
    assert result.exit_code == 0, result.output
    assert "plot_loci" in calls


def test_viz_interactive_calls_index(runner, monkeypatch):
    calls = {}
    _patch_viz_pipeline(monkeypatch)
    monkeypatch.setattr(cli, "interactive_index",
                        lambda genes, out_dir: calls.setdefault("idx", out_dir))
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "viz", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--out-dir", "figs", "--mode", "interactive",
        ])
    assert result.exit_code == 0, result.output
    assert "idx" in calls


def test_viz_tracks_calls_write_bed12(runner, monkeypatch):
    calls = {}
    _patch_viz_pipeline(monkeypatch)
    monkeypatch.setattr(cli, "write_bed12",
                        lambda genes, path: calls.setdefault("bed", path))
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "viz", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--out-dir", "tracks", "--mode", "tracks",
        ])
    assert result.exit_code == 0, result.output
    assert "bed" in calls


def test_viz_unknown_gene_errors(runner, monkeypatch):
    _patch_viz_pipeline(monkeypatch, genes=[_FakeGene("HFG_00001")])
    monkeypatch.setattr(cli, "plot_loci", lambda *a, **k: None)
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "viz", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--out-dir", "figs", "--gene", "HFG_99999",
        ])
    assert result.exit_code != 0
    assert "not found" in result.output


# ---------------------------------------------------------------------------
# benchmark
# ---------------------------------------------------------------------------


def test_benchmark_all_calls_api(runner, monkeypatch):
    calls = {}

    def fake_benchmark_all(reconciled, proteins, out_dir, **kw):
        calls["args"] = (reconciled, proteins, out_dir, kw)
        return pd.DataFrame([{"tool": "agat", "metric": "n", "value": 1.0}])

    monkeypatch.setattr(cli, "benchmark_all", fake_benchmark_all)
    with runner.isolated_filesystem():
        from pathlib import Path

        Path("recon.gff3").write_text("##gff-version 3\n")
        Path("proteins.fa").write_text(">p\nMA\n")
        result = runner.invoke(cli.main, [
            "benchmark", "all", "--reconciled", "recon.gff3",
            "--proteins", "proteins.fa", "--out-dir", "bench",
            "--lineage", "brassicales_odb10",
        ])
    assert result.exit_code == 0, result.output
    recon, prot, out_dir, kw = calls["args"]
    assert recon == "recon.gff3" and prot == "proteins.fa"
    assert kw["lineage"] == "brassicales_odb10"


def test_benchmark_ablation_calls_api(runner, monkeypatch):
    calls = {}

    def fake_run_ablation(config, variants, out_dir, **kw):
        calls["args"] = (config, variants, out_dir)
        return pd.DataFrame([{"variant": v} for v in variants])

    monkeypatch.setattr(cli, "run_ablation", fake_run_ablation)
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "benchmark", "ablation", "--genome", "genome.fa",
            "--helixer", "helixer.gff3", "--out-dir", "abl",
            "--variants", "full,no_pad",
        ])
    assert result.exit_code == 0, result.output
    config, variants, out_dir = calls["args"]
    assert variants == ["full", "no_pad"]
    assert config.genome_fasta == "genome.fa"


# ---------------------------------------------------------------------------
# Chunked reconcile (`--scatter`) — folded-in `run`: preflight + driver wiring
# ---------------------------------------------------------------------------


def _write_genome(name="genome.fa", seqids=("chr1", "chr2")):
    from pathlib import Path
    Path(name).write_text("".join(f">{s}\n{'ACGT' * 16}\n" for s in seqids))


def _write_helixer(name="helixer.gff3", seqids=("chr1",)):
    from pathlib import Path
    lines = ["##gff-version 3"]
    for i, s in enumerate(seqids):
        lines.append(f"{s}\tHelixer\tgene\t100\t200\t.\t+\t.\tID=g{i}")
    Path(name).write_text("\n".join(lines) + "\n")


def test_reconcile_scatter_blocks_on_disjoint_seqids(runner):
    with runner.isolated_filesystem():
        _write_genome(seqids=("chr1", "chr2"))
        _write_helixer(seqids=("1",))  # disjoint naming → preflight error
        result = runner.invoke(
            cli.main,
            ["reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
             "--scatter", "auto"])
        assert result.exit_code != 0
        assert "preflight" in result.output.lower()


def test_reconcile_single_pass_skips_preflight(runner, monkeypatch):
    # The default single in-process pass does not run the preflight; it goes
    # straight to run_pipeline even with inconsistent seqids (doctor owns the gate).
    captured = {}

    def fake_run_pipeline(config):
        captured["ran"] = True
        return []

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    with runner.isolated_filesystem():
        _write_genome(seqids=("chr1",))
        _write_helixer(seqids=("1",))  # disjoint, but no preflight on single pass
        result = runner.invoke(
            cli.main, ["reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3"])
        assert result.exit_code == 0, result.output
        assert captured.get("ran")


def test_reconcile_scatter_skip_preflight_reaches_driver(runner, monkeypatch):
    from helixforge.parallel import run as run_mod

    called = {}

    def fake_run_genome(base_config, **kw):
        called["yes"] = True
        return run_mod.RunResult(
            mode="single", genes=[], manifest_path="m.json", out_prefix="hf")

    monkeypatch.setattr(run_mod, "run_genome", fake_run_genome)
    with runner.isolated_filesystem():
        _write_genome(seqids=("chr1",))
        _write_helixer(seqids=("1",))  # would fail preflight, but we skip it
        result = runner.invoke(
            cli.main,
            ["reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
             "--scatter", "auto", "--skip-preflight"])
        assert result.exit_code == 0, result.output
        assert called.get("yes")


def test_reconcile_scatter_passes_concordant_inputs(runner, monkeypatch):
    from helixforge.parallel import run as run_mod

    monkeypatch.setattr(
        run_mod, "run_genome",
        lambda base_config, **kw: run_mod.RunResult(
            mode="single", genes=[], manifest_path="m.json", out_prefix="hf"))
    with runner.isolated_filesystem():
        _write_genome(seqids=("chr1", "chr2"))
        _write_helixer(seqids=("chr1",))  # concordant subset
        result = runner.invoke(
            cli.main,
            ["reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
             "--scatter", "auto"])
        assert result.exit_code == 0, result.output


def test_doctor_runs_integrity_gate(runner, monkeypatch):
    from helixforge.prep import doctor as doc
    monkeypatch.setattr(doc, "tool_version", lambda b, a="--version": "1.0")
    with runner.isolated_filesystem():
        _write_genome(seqids=("chr1", "chr2"))
        _write_helixer(seqids=("1",))  # disjoint → preflight fails
        result = runner.invoke(
            cli.main, ["doctor", "--genome", "genome.fa", "--helixer", "helixer.gff3"])
        assert result.exit_code != 0
        assert "concordance" in result.output.lower() or "preflight" in result.output.lower()
