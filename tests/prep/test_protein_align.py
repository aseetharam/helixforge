"""Tests for prep/protein_align.py (Phase 12). Floor: 5. Tools fully mocked."""

import subprocess

import pytest

from helixforge.prep import protein_align as pa


@pytest.fixture
def calls(monkeypatch):
    recorded = []

    def fake_run_tool(argv, **kw):
        recorded.append(([str(a) for a in argv], kw))
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(pa, "run_tool", fake_run_tool)
    return recorded


def test_miniprot_argv_and_path(calls, tmp_path):
    out = pa.run_miniprot("genome.fa", "prot.fa", tmp_path / "mp.gff", threads=8)
    argv, kw = calls[0]
    assert argv[0] == "miniprot"
    assert "--gff" in argv
    assert argv[argv.index("-t") + 1] == "8"
    assert argv[-2:] == ["genome.fa", "prot.fa"]
    assert str(out) == str(tmp_path / "mp.gff")


def test_miniprot_redirects_stdout_to_out_gff(calls, tmp_path):
    out_gff = tmp_path / "mp.gff"
    pa.run_miniprot("genome.fa", "prot.fa", out_gff)
    _, kw = calls[0]
    assert str(kw["stdout_path"]) == str(out_gff)


def test_miniprot_custom_threads(calls, tmp_path):
    pa.run_miniprot("g.fa", "p.fa", tmp_path / "mp.gff", threads=24)
    argv, _ = calls[0]
    assert argv[argv.index("-t") + 1] == "24"


def test_miniprot_custom_bin_and_extra_args(calls, tmp_path):
    pa.run_miniprot(
        "g.fa", "p.fa", tmp_path / "mp.gff",
        miniprot_bin="/opt/miniprot", extra_args=["-p", "0.5"],
    )
    argv, _ = calls[0]
    assert argv[0] == "/opt/miniprot"
    assert "-p" in argv and argv[argv.index("-p") + 1] == "0.5"


def test_miniprot_error_propagates(monkeypatch, tmp_path):
    def boom(argv, **kw):
        raise RuntimeError("miniprot failed (exit 1)")

    monkeypatch.setattr(pa, "run_tool", boom)
    with pytest.raises(RuntimeError, match="miniprot failed"):
        pa.run_miniprot("g.fa", "p.fa", tmp_path / "mp.gff")
