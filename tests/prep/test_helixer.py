"""Tests for prep/helixer.py (Phase 12). Floor: 3. OPTIONAL GPU tool, mocked."""

import subprocess

import pytest

from helixforge.prep import helixer as hx
from helixforge.prep.helixer import HelixerResult, run_helixer


@pytest.fixture
def calls(monkeypatch):
    recorded = []

    def fake_run_tool(argv, **kw):
        recorded.append([str(a) for a in argv])
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(hx, "run_tool", fake_run_tool)
    return recorded


def test_helixer_argv_and_result_paths(calls, tmp_path):
    res = run_helixer("genome.fa", tmp_path / "hx", lineage="land_plant", subseq_len=64152)
    argv = calls[0]
    assert argv[0] == "Helixer.py"
    assert argv[argv.index("--fasta-path") + 1] == "genome.fa"
    assert argv[argv.index("--lineage") + 1] == "land_plant"
    assert argv[argv.index("--subsequence-length") + 1] == "64152"
    assert isinstance(res, HelixerResult)
    assert str(res.gff3) == str(tmp_path / "hx" / "helixer.gff3")
    assert str(res.hdf5) == str(tmp_path / "hx" / "helixer.h5")
    assert argv[argv.index("--gff-output-path") + 1] == str(res.gff3)


def test_helixer_custom_bin(calls, tmp_path):
    run_helixer("g.fa", tmp_path / "hx", helixer_bin="/opt/Helixer.py")
    assert calls[0][0] == "/opt/Helixer.py"


def test_helixer_docstring_marks_gpu_and_off_default_path():
    doc = run_helixer.__doc__
    assert "GPU" in doc
    # documented as off the default prep path (most users supply Helixer output)
    assert "OFF the default" in doc or "off the default" in doc.lower()
