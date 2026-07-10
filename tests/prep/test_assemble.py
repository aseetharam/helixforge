"""Tests for prep/assemble.py (Phase 12). Floor: 6. Tools fully mocked."""

import subprocess

import pytest

from helixforge.prep import assemble as asm


@pytest.fixture
def calls(monkeypatch):
    recorded = []

    def fake_run_tool(argv, **kw):
        recorded.append([str(a) for a in argv])
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(asm, "run_tool", fake_run_tool)
    return recorded


def test_stringtie_argv_label_and_path(calls, tmp_path):
    out = asm.run_stringtie("aln.bam", tmp_path / "s.gtf", "sampA", threads=4)
    argv = calls[0]
    assert argv[0] == "stringtie"
    assert argv[1] == "aln.bam"
    assert argv[argv.index("-o") + 1] == str(tmp_path / "s.gtf")
    assert argv[argv.index("-p") + 1] == "4"
    assert argv[argv.index("-l") + 1] == "sampA"  # unique label per sample
    assert str(out) == str(tmp_path / "s.gtf")


def test_stringtie_guide_gtf(calls, tmp_path):
    asm.run_stringtie("aln.bam", tmp_path / "s.gtf", "s1", guide_gtf="ref.gtf")
    argv = calls[0]
    assert argv[argv.index("-G") + 1] == "ref.gtf"


def test_stringtie_extra_args_and_custom_bin(calls, tmp_path):
    asm.run_stringtie(
        "aln.bam", tmp_path / "s.gtf", "s1",
        stringtie_bin="/opt/stringtie", extra_args=["-m", "200"],
    )
    argv = calls[0]
    assert argv[0] == "/opt/stringtie"
    assert "-m" in argv and argv[argv.index("-m") + 1] == "200"


def test_assemble_samples_one_gtf_per_bam(calls, tmp_path):
    gtfs = asm.assemble_samples(["x/sampA.bam", "y/sampB.bam"], tmp_path / "gtfs")
    assert len(gtfs) == 2
    assert len(calls) == 2
    assert str(gtfs[0]) == str(tmp_path / "gtfs" / "sampA.gtf")
    assert str(gtfs[1]) == str(tmp_path / "gtfs" / "sampB.gtf")


def test_assemble_samples_uses_bam_stem_as_label(calls, tmp_path):
    asm.assemble_samples(["dir/repB.bam"], tmp_path / "g")
    argv = calls[0]
    assert argv[argv.index("-l") + 1] == "repB"


def test_assemble_samples_passes_through_kwargs(calls, tmp_path):
    asm.assemble_samples(
        ["s1.bam"], tmp_path / "g", threads=7, stringtie_bin="/opt/st", guide_gtf="ref.gtf"
    )
    argv = calls[0]
    assert argv[0] == "/opt/st"
    assert argv[argv.index("-p") + 1] == "7"
    assert argv[argv.index("-G") + 1] == "ref.gtf"


def test_stringtie_error_propagates(monkeypatch, tmp_path):
    def boom(argv, **kw):
        raise RuntimeError("stringtie failed (exit 1)")

    monkeypatch.setattr(asm, "run_tool", boom)
    with pytest.raises(RuntimeError, match="stringtie failed"):
        asm.run_stringtie("aln.bam", tmp_path / "s.gtf", "s1")
