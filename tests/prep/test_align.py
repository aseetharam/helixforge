"""Tests for prep/align.py (Phase 12). Floor: 12. Tools fully mocked."""

import subprocess

import pytest

from helixforge.prep import align as align_mod
from helixforge.prep.align import AlignResult, run_hisat2, run_star


@pytest.fixture
def calls(monkeypatch):
    """Capture stringified argv from run_tool and run_pipe at the align module."""
    tool_calls = []
    pipe_calls = []

    def fake_run_tool(argv, **kw):
        tool_calls.append([str(a) for a in argv])
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    def fake_run_pipe(commands, **kw):
        pipe_calls.append([[str(a) for a in argv] for argv in commands])
        return []

    monkeypatch.setattr(align_mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(align_mod, "run_pipe", fake_run_pipe)
    return tool_calls, pipe_calls


# --- STAR ---

def test_star_builds_index_when_given_fasta(calls, tmp_path):
    tool_calls, _ = calls
    run_star(str(tmp_path / "genome.fa"), [str(tmp_path / "r1.fq")], str(tmp_path / "out_"))
    gen = tool_calls[0]
    assert gen[0] == "STAR"
    assert "--runMode" in gen and gen[gen.index("--runMode") + 1] == "genomeGenerate"
    assert "--genomeFastaFiles" in gen


def test_star_uses_dir_without_rebuild(calls, tmp_path):
    tool_calls, _ = calls
    run_star(str(tmp_path / "star_index"), [str(tmp_path / "r1.fq")], str(tmp_path / "out_"))
    # only the align command runs (no genomeGenerate)
    assert len(tool_calls) == 1
    assert "alignReads" in tool_calls[0]


def test_star_paired_reads(calls, tmp_path):
    tool_calls, _ = calls
    run_star(str(tmp_path / "idx"), ["r1.fq", "r2.fq"], str(tmp_path / "o_"))
    argv = tool_calls[0]
    i = argv.index("--readFilesIn")
    assert argv[i + 1] == "r1.fq" and argv[i + 2] == "r2.fq"


def test_star_single_read(calls, tmp_path):
    tool_calls, _ = calls
    run_star(str(tmp_path / "idx"), ["r1.fq"], str(tmp_path / "o_"))
    argv = tool_calls[0]
    i = argv.index("--readFilesIn")
    assert argv[i + 1] == "r1.fq"
    assert argv[i + 2] == "--runThreadN"  # only one read file


def test_star_gzip_reads_add_zcat(calls, tmp_path):
    tool_calls, _ = calls
    run_star(str(tmp_path / "idx"), ["r1.fq.gz", "r2.fq.gz"], str(tmp_path / "o_"))
    argv = tool_calls[0]
    assert "--readFilesCommand" in argv
    assert argv[argv.index("--readFilesCommand") + 1] == "zcat"


def test_star_two_pass(calls, tmp_path):
    tool_calls, _ = calls
    run_star(str(tmp_path / "idx"), ["r1.fq"], str(tmp_path / "o_"), two_pass=True)
    argv = tool_calls[0]
    assert "--twopassMode" in argv and argv[argv.index("--twopassMode") + 1] == "Basic"


def test_star_threads_sjdb_and_extra(calls, tmp_path):
    tool_calls, _ = calls
    run_star(
        str(tmp_path / "idx"), ["r1.fq"], str(tmp_path / "o_"),
        threads=16, sjdb_gtf="ann.gtf", extra_args=["--outFilterMultimapNmax", "1"],
    )
    argv = tool_calls[0]
    assert argv[argv.index("--runThreadN") + 1] == "16"
    assert argv[argv.index("--sjdbGTFfile") + 1] == "ann.gtf"
    assert "--outFilterMultimapNmax" in argv


def test_star_custom_bin(calls, tmp_path):
    tool_calls, _ = calls
    run_star(str(tmp_path / "idx"), ["r1.fq"], str(tmp_path / "o_"), star_bin="/opt/STAR")
    assert tool_calls[0][0] == "/opt/STAR"


def test_star_result_paths(calls, tmp_path):
    _, _ = calls
    res = run_star(str(tmp_path / "idx"), ["r1.fq"], "/tmp/run_")
    assert isinstance(res, AlignResult)
    assert str(res.bam) == "/tmp/run_Aligned.sortedByCoord.out.bam"
    assert str(res.sj_tab) == "/tmp/run_SJ.out.tab"
    assert str(res.log) == "/tmp/run_Log.final.out"


def test_star_error_propagates(monkeypatch, tmp_path):
    def boom(argv, **kw):
        raise RuntimeError("STAR failed (exit 1)")

    monkeypatch.setattr(align_mod, "run_tool", boom)
    with pytest.raises(RuntimeError, match="STAR failed"):
        run_star(str(tmp_path / "idx"), ["r1.fq"], str(tmp_path / "o_"))


# --- HISAT2 ---

def test_hisat2_paired_pipe_and_index(calls, tmp_path):
    tool_calls, pipe_calls = calls
    res = run_hisat2("idx", ["r1.fq", "r2.fq"], tmp_path / "out.bam")
    # one pipeline (hisat2 | samtools sort) + one index call
    assert len(pipe_calls) == 1
    hisat2_argv, sort_argv = pipe_calls[0]
    assert hisat2_argv[0] == "hisat2"
    assert "-1" in hisat2_argv and "-2" in hisat2_argv
    assert hisat2_argv[hisat2_argv.index("-x") + 1] == "idx"
    assert sort_argv[0] == "samtools" and sort_argv[-1] == "-"
    assert tool_calls[-1][:2] == ["samtools", "index"]
    assert res.sj_tab is None
    assert str(res.bam).endswith("out.bam")


def test_hisat2_single_read(calls, tmp_path):
    _, pipe_calls = calls
    run_hisat2("idx", ["r1.fq"], tmp_path / "out.bam")
    hisat2_argv, _ = pipe_calls[0]
    assert "-U" in hisat2_argv
    assert "-1" not in hisat2_argv


def test_hisat2_builds_index_from_fasta(calls, tmp_path):
    tool_calls, _ = calls
    run_hisat2(str(tmp_path / "genome.fa"), ["r1.fq"], tmp_path / "out.bam")
    build = tool_calls[0]
    assert build[0] == "hisat2-build"
    # index prefix is the FASTA path with its suffix stripped
    assert build[-1] == str(tmp_path / "genome")


def test_hisat2_custom_bins(calls, tmp_path):
    tool_calls, pipe_calls = calls
    run_hisat2(
        "idx", ["r1.fq"], tmp_path / "out.bam",
        hisat2_bin="/opt/hisat2", samtools_bin="/opt/samtools",
    )
    hisat2_argv, sort_argv = pipe_calls[0]
    assert hisat2_argv[0] == "/opt/hisat2"
    assert sort_argv[0] == "/opt/samtools"
    assert tool_calls[-1][0] == "/opt/samtools"
