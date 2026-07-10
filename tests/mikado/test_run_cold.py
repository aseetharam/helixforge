"""Cold-run hardening of mikado/run.py (Phase 13 D2). Floor: 10.

Version parsing, the serialise genome-flag branch, tool-resolution preflight,
and the ``resume`` argument (short-circuit existing outputs). Subprocess and tool
resolution fully mocked — no real Mikado/TransDecoder/DIAMOND.
"""

import subprocess

import pytest

from helixforge.mikado import run as run_mod
from helixforge.mikado.run import (
    MikadoInputs,
    _parse_version_tuple,
    _require_tool,
    _serialise_genome_flag,
    mikado_version,
    run_mikado,
)


def _inputs():
    return MikadoInputs(
        configuration_yaml="config.yaml",
        scoring_file="scoring.yaml",
        genome_fa="genome.fa",
        protein_db="proteins.fa",
        junctions_tab="junctions.tab",
        external_tsv="external.tsv",
    )


# --- version parsing ---


def test_parse_version_three_components():
    assert _parse_version_tuple("Mikado v2.3.4") == (2, 3, 4)


def test_parse_version_two_components():
    assert _parse_version_tuple("diamond version 2.1") == (2, 1)


def test_parse_version_embedded():
    assert _parse_version_tuple("STAR 2.7.11b") == (2, 7, 11)


def test_parse_version_none_and_garbage():
    assert _parse_version_tuple(None) is None
    assert _parse_version_tuple("no digits here") is None


def test_mikado_version_uses_tool_version(monkeypatch):
    monkeypatch.setattr(run_mod, "tool_version", lambda b, a="--version": "Mikado 2.3.4")
    assert mikado_version("mikado") == (2, 3, 4)


def test_mikado_version_missing_returns_none(monkeypatch):
    monkeypatch.setattr(run_mod, "tool_version", lambda b, a="--version": None)
    assert mikado_version("mikado") is None


# --- serialise genome-flag branch ---


def test_genome_flag_modern():
    assert _serialise_genome_flag((2, 3, 4)) == "--genome"


def test_genome_flag_legacy():
    assert _serialise_genome_flag((2, 2, 0)) == "--genome_fasta"


def test_genome_flag_unknown_defaults_modern():
    assert _serialise_genome_flag(None) == "--genome"


def test_serialise_uses_legacy_flag_for_old_version(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(
        run_mod.subprocess, "run",
        lambda argv, check, capture_output, text, cwd=None: calls.append(argv)
        or subprocess.CompletedProcess(argv, 0, stdout="", stderr=""),
    )
    run_mod.run_serialise(
        "config.yaml", "prep.fasta", "orfs.bed", "blast.xml", "db.fa",
        "junctions.tab", "ext.tsv", "genome.fa", tmp_path, version=(2, 2, 0),
    )
    assert "--genome_fasta" in calls[0]
    assert "--genome" not in calls[0]


# --- tool resolution preflight ---


def test_require_tool_missing_raises(monkeypatch):
    monkeypatch.setattr(run_mod.shutil, "which", lambda name: None)
    with pytest.raises(FileNotFoundError, match="mikado"):
        _require_tool("mikado", "mikado")


def test_require_tool_found_ok(monkeypatch):
    monkeypatch.setattr(run_mod.shutil, "which", lambda name: "/usr/bin/mikado")
    _require_tool("mikado", "mikado")  # no raise


def test_require_tool_path_must_exist(tmp_path):
    missing = tmp_path / "nope" / "mikado"
    with pytest.raises(FileNotFoundError):
        _require_tool(str(missing), "mikado")


def test_run_mikado_missing_tool_precise_error(monkeypatch, tmp_path):
    monkeypatch.setattr(run_mod.shutil, "which", lambda name: None)
    monkeypatch.setattr(run_mod, "tool_version", lambda b, a="--version": None)
    with pytest.raises(FileNotFoundError, match="mikado"):
        run_mikado(_inputs(), tmp_path)


# --- resume short-circuits existing outputs ---


@pytest.fixture
def cold_env(monkeypatch):
    """Mock subprocess + resolution; record argvs run."""
    calls = []

    def fake_run(argv, check, capture_output, text, cwd=None):
        calls.append([str(a) for a in argv])
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(run_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(run_mod.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(run_mod, "tool_version", lambda b, a="--version": "mikado 2.3.4")
    return calls


def test_run_mikado_no_resume_runs_full_chain(cold_env, tmp_path):
    run_mikado(_inputs(), tmp_path, resume=False)
    steps = [c[1] if c[0].endswith(("mikado", "diamond")) else c[0] for c in cold_env]
    assert steps[0] == "prepare"
    assert "LongOrfs" in steps[1]
    assert "Predict" in steps[2]
    assert steps[3] == "makedb"
    assert steps[4] == "blastx"
    assert steps[5] == "serialise"
    assert steps[6] == "pick"


def test_run_mikado_resume_skips_existing(cold_env, tmp_path):
    # Pre-create every expected step output → resume must run NO subprocess.
    (tmp_path / "mikado_prepared.gtf").write_text("x")
    (tmp_path / "mikado_prepared.fasta").write_text("x")
    (tmp_path / "mikado_prepared.fasta.transdecoder.bed").write_text("x")
    (tmp_path / "mikado_diamond.xml").write_text("x")
    (tmp_path / "mikado.db").write_text("x")
    (tmp_path / "mikado.loci.gff3").write_text("x")
    result = run_mikado(_inputs(), tmp_path, resume=True)
    assert cold_env == []
    assert result.loci_gff3.name == "mikado.loci.gff3"


def test_run_mikado_resume_partial_runs_only_missing(cold_env, tmp_path):
    # prepare + transdecoder + diamond outputs exist; serialise + pick missing.
    (tmp_path / "mikado_prepared.gtf").write_text("x")
    (tmp_path / "mikado_prepared.fasta").write_text("x")
    (tmp_path / "mikado_prepared.fasta.transdecoder.bed").write_text("x")
    (tmp_path / "mikado_diamond.xml").write_text("x")
    run_mikado(_inputs(), tmp_path, resume=True)
    steps = [c[1] if c[0].endswith(("mikado", "diamond")) else c[0] for c in cold_env]
    assert steps == ["serialise", "pick"]
