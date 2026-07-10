"""TRaCE CLI flags (Phase 33b D4). Floor: 4.

``--trace-primary`` + the ``--trace-*`` params parse into PipelineConfig, default
off, and appear in ``reconcile --help``.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from helixforge import cli


@pytest.fixture
def runner():
    return CliRunner()


def _touch_inputs():
    from pathlib import Path

    Path("genome.fa").write_text(">chr1\nACGTACGTAC\n")
    Path("helixer.gff3").write_text("##gff-version 3\n")


def test_trace_flags_map_to_config(runner, monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda config: captured.update(config=config) or [])
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--trace-primary",
            "--trace-max-aed", "0.3",
            "--trace-min-tpm", "1.0",
            "--trace-min-overlap", "0.4",
            "--trace-weight-domain", "12",
            "--trace-weight-protein", "4",
            "--trace-weight-cdna", "2",
            "--no-trace-use-domain",
        ])
    assert result.exit_code == 0, result.output
    cfg = captured["config"]
    assert cfg.trace_primary is True
    assert cfg.trace_max_aed == 0.3
    assert cfg.trace_min_tpm == 1.0
    assert cfg.trace_min_overlap == 0.4
    assert cfg.trace_weight_domain == 12.0
    assert cfg.trace_weight_protein == 4.0
    assert cfg.trace_weight_cdna == 2.0
    assert cfg.trace_use_domain is False


def test_trace_primary_off_by_default(runner, monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda config: captured.update(config=config) or [])
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
        ])
    assert result.exit_code == 0, result.output
    cfg = captured["config"]
    # Default reproduces the golden path: TRaCE off, TRaCE-paper defaults.
    assert cfg.trace_primary is False
    assert cfg.trace_max_aed == 0.5
    assert cfg.trace_weight_domain == 9.0
    assert cfg.trace_use_domain is True


def test_trace_help_lists_flags(runner):
    result = runner.invoke(cli.main, ["reconcile", "--help"])
    assert result.exit_code == 0
    for flag in ("--trace-primary", "--trace-max-aed", "--trace-weight-domain",
                 "--trace-use-domain"):
        assert flag in result.output


def test_trace_params_recorded_in_resolved_config(runner, monkeypatch):
    # D4: the params must be reproducible — they appear in the nested run config
    # (which provenance hashes).
    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda config: captured.update(config=config) or [])
    with runner.isolated_filesystem():
        _touch_inputs()
        runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--trace-primary", "--trace-weight-protein", "7",
        ])
    nested = captured["config"].to_nested_dict()
    assert nested["trace_primary"] is True
    assert nested["trace_weight_protein"] == 7.0
