"""Tests for prep/samtools.py (Phase 12). Floor: 8. Tools fully mocked."""

import subprocess

import pytest

from helixforge.prep import samtools as st


@pytest.fixture
def calls(monkeypatch):
    recorded = []

    def fake_run_tool(argv, **kw):
        recorded.append([str(a) for a in argv])
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(st, "run_tool", fake_run_tool)
    return recorded


def test_sort_argv_and_path(calls, tmp_path):
    out = st.sort("in.bam", tmp_path / "sorted.bam", threads=4)
    argv = calls[0]
    assert argv[:2] == ["samtools", "sort"]
    assert argv[argv.index("-@") + 1] == "4"
    assert argv[argv.index("-o") + 1] == str(tmp_path / "sorted.bam")
    assert argv[-1] == "in.bam"
    assert str(out) == str(tmp_path / "sorted.bam")


def test_sort_custom_threads(calls, tmp_path):
    st.sort("in.bam", tmp_path / "s.bam", threads=12)
    assert calls[0][calls[0].index("-@") + 1] == "12"


def test_index_argv_and_bai_path(calls):
    bai = st.index("aln.bam")
    assert calls[0] == ["samtools", "index", "aln.bam"]
    assert str(bai) == "aln.bam.bai"


def test_merge_argv_and_path(calls, tmp_path):
    out = st.merge(["a.bam", "b.bam", "c.bam"], tmp_path / "merged.bam", threads=8)
    argv = calls[0]
    assert argv[:2] == ["samtools", "merge"]
    assert argv[argv.index("-@") + 1] == "8"
    assert "-f" in argv
    assert argv[-3:] == ["a.bam", "b.bam", "c.bam"]
    assert str(out) == str(tmp_path / "merged.bam")


def test_faidx_argv_and_fai_path(calls):
    fai = st.faidx("genome.fa")
    assert calls[0] == ["samtools", "faidx", "genome.fa"]
    assert str(fai) == "genome.fa.fai"


def test_sort_custom_bin(calls, tmp_path):
    st.sort("in.bam", tmp_path / "s.bam", samtools_bin="/opt/samtools")
    assert calls[0][0] == "/opt/samtools"


def test_faidx_custom_bin(calls):
    st.faidx("g.fa", samtools_bin="/opt/samtools")
    assert calls[0][0] == "/opt/samtools"


def test_sort_error_propagates(monkeypatch, tmp_path):
    def boom(argv, **kw):
        raise RuntimeError("samtools sort failed (exit 1)")

    monkeypatch.setattr(st, "run_tool", boom)
    with pytest.raises(RuntimeError, match="samtools sort failed"):
        st.sort("in.bam", tmp_path / "s.bam")
