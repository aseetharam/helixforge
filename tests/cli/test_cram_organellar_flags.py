"""CLI exposure of CRAM `--reference` + organellar genetic-code flags (Phase 33 D1/D2).

Floor: 9. Asserts the new flags parse, reach ``PipelineConfig`` / the evidence
opener, validate bad input with a clear error, and appear in ``--help``. The
Python API is patched so nothing real is opened.
"""

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from helixforge import cli


@pytest.fixture
def runner():
    return CliRunner()


def _touch_inputs():
    Path("genome.fa").write_text(">chrMt\nACGTACGTAC\n>chrPt\nACGTACGTAC\n")
    Path("helixer.gff3").write_text("##gff-version 3\n")


# --- D2: --transl-table / --transl-table-map → PipelineConfig ---------------

def test_transl_table_default_is_one(runner, monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda c: (captured.__setitem__("c", c), [])[1])
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
        ])
    assert result.exit_code == 0, result.output
    assert captured["c"].transl_table == 1
    assert captured["c"].transl_table_map is None


def test_transl_table_flows_to_config(runner, monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda c: (captured.__setitem__("c", c), [])[1])
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--transl-table", "11",
        ])
    assert result.exit_code == 0, result.output
    assert captured["c"].transl_table == 11


def test_transl_table_map_equals_form_parses(runner, monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda c: (captured.__setitem__("c", c), [])[1])
    with runner.isolated_filesystem():
        _touch_inputs()
        Path("codes.txt").write_text("# organellar codes\nchrMt=1\nchrPt=11\n")
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--transl-table-map", "codes.txt",
        ])
    assert result.exit_code == 0, result.output
    assert captured["c"].transl_table_map == {"chrMt": 1, "chrPt": 11}


def test_transl_table_map_whitespace_form_parses(runner, monkeypatch):
    captured = {}
    monkeypatch.setattr(cli, "run_pipeline", lambda c: (captured.__setitem__("c", c), [])[1])
    with runner.isolated_filesystem():
        _touch_inputs()
        Path("codes.tsv").write_text("chrMt\t1\nchrPt 11\n")
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--transl-table-map", "codes.tsv",
        ])
    assert result.exit_code == 0, result.output
    assert captured["c"].transl_table_map == {"chrMt": 1, "chrPt": 11}


def test_bad_transl_table_errors_clearly(runner, monkeypatch):
    monkeypatch.setattr(cli, "run_pipeline", lambda c: [])
    with runner.isolated_filesystem():
        _touch_inputs()
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--transl-table", "999",
        ])
    assert result.exit_code != 0
    assert "unknown NCBI transl_table id 999" in result.output


def test_bad_table_in_map_errors_clearly(runner, monkeypatch):
    monkeypatch.setattr(cli, "run_pipeline", lambda c: [])
    with runner.isolated_filesystem():
        _touch_inputs()
        Path("codes.txt").write_text("chrMt=999\n")
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--transl-table-map", "codes.txt",
        ])
    assert result.exit_code != 0
    assert "unknown NCBI transl_table id 999" in result.output


def test_bad_seqid_in_map_errors_against_genome(runner, monkeypatch):
    monkeypatch.setattr(cli, "run_pipeline", lambda c: [])
    with runner.isolated_filesystem():
        _touch_inputs()  # genome has chrMt / chrPt only
        Path("codes.txt").write_text("ChrMt=1\n")  # wrong case -> not in genome
        result = runner.invoke(cli.main, [
            "reconcile", "--genome", "genome.fa", "--helixer", "helixer.gff3",
            "--transl-table-map", "codes.txt",
        ])
    assert result.exit_code != 0
    assert "not in the genome" in result.output


# --- D1: evidence --reference reaches score_annotation ----------------------

def test_evidence_reference_reaches_scorer(runner, monkeypatch):
    import helixforge.score.evidence as ev

    captured = {}

    def fake_score(gff3_path, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame()

    monkeypatch.setattr(ev, "score_annotation", fake_score)
    monkeypatch.setattr(ev, "summarize_evidence", lambda df: {"rows": 0})
    monkeypatch.setattr(ev, "write_evidence_tsv", lambda df, p: Path(p).write_text(""))

    with runner.isolated_filesystem():
        Path("ann.gff3").write_text("##gff-version 3\n")
        Path("reads.cram").write_text("CRAMdummy")
        Path("genome.fa").write_text(">chr1\nACGT\n")
        result = runner.invoke(cli.main, [
            "evidence", "--gff3", "ann.gff3", "--bam", "reads.cram",
            "--reference", "genome.fa",
        ])
    assert result.exit_code == 0, result.output
    assert captured["reference_filename"] == "genome.fa"


def test_help_shows_new_flags(runner):
    ev_help = runner.invoke(cli.main, ["evidence", "--help"]).output
    assert "--reference" in ev_help
    rec_help = runner.invoke(cli.main, ["reconcile", "--help"]).output
    assert "--transl-table" in rec_help
    assert "--transl-table-map" in rec_help
    doc_help = runner.invoke(cli.main, ["doctor", "--help"]).output
    assert "--reference" in doc_help
