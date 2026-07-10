"""Tests for prep/_subprocess.py (Phase 12). Subprocess fully mocked."""

import subprocess

import pytest

from helixforge.prep import _subprocess as sp


class _FakeProc:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def test_run_tool_success_stringifies_argv(monkeypatch):
    seen = {}

    def fake_run(argv, **kw):
        seen["argv"] = argv
        seen["kw"] = kw
        return subprocess.CompletedProcess(argv, 0, stdout="ok", stderr="")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    proc = sp.run_tool(["STAR", "--runThreadN", 8])
    assert seen["argv"] == ["STAR", "--runThreadN", "8"]  # ints stringified
    assert seen["kw"]["check"] is True
    assert proc.returncode == 0


def test_run_tool_nonzero_raises_with_argv_and_stderr_tail(monkeypatch):
    def fake_run(argv, **kw):
        raise subprocess.CalledProcessError(2, argv, output="", stderr="boom-detail")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError) as exc:
        sp.run_tool(["samtools", "sort"])
    msg = str(exc.value)
    assert "samtools" in msg
    assert "exit 2" in msg
    assert "boom-detail" in msg


def test_run_tool_missing_binary_raises_tool_not_found(monkeypatch):
    def fake_run(argv, **kw):
        raise FileNotFoundError(argv[0])

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    with pytest.raises(FileNotFoundError) as exc:
        sp.run_tool(["nope_tool", "--help"])
    assert "tool not found" in str(exc.value)
    assert "nope_tool" in str(exc.value)


def test_run_tool_stdout_path_redirects(monkeypatch, tmp_path):
    out = tmp_path / "result.gff"
    seen = {}

    def fake_run(argv, **kw):
        seen["kw"] = kw
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    sp.run_tool(["miniprot", "--gff"], stdout_path=out)
    # capture_output must NOT be set when redirecting; stdout fh is passed.
    assert "capture_output" not in seen["kw"]
    assert seen["kw"]["stderr"] is subprocess.PIPE
    assert out.exists()


def test_run_tool_cwd_passed(monkeypatch, tmp_path):
    seen = {}

    def fake_run(argv, **kw):
        seen["kw"] = kw
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    sp.run_tool(["stringtie"], cwd=tmp_path)
    assert seen["kw"]["cwd"] == str(tmp_path)


def test_tool_version_returns_first_nonempty_line(monkeypatch):
    def fake_run(argv, **kw):
        return _FakeProc(0, stdout="\n2.7.10a\nextra\n", stderr="")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    assert sp.tool_version("STAR") == "2.7.10a"


def test_tool_version_reads_stderr_when_stdout_empty(monkeypatch):
    def fake_run(argv, **kw):
        return _FakeProc(0, stdout="", stderr="stringtie 2.2.1")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    assert sp.tool_version("stringtie") == "stringtie 2.2.1"


def test_tool_version_missing_binary_returns_none(monkeypatch):
    def fake_run(argv, **kw):
        raise FileNotFoundError(argv[0])

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    assert sp.tool_version("ghost") is None


def test_tool_version_uses_errors_replace(monkeypatch):
    # A tool printing a non-UTF-8 version banner must not crash tool_version
    # (subprocess.run is asked for errors="replace"); regression for the
    # UnicodeDecodeError that crashed `helixforge doctor` live.
    seen = {}

    def fake_run(argv, **kw):
        seen.update(kw)
        return _FakeProc(0, stdout="miniprot 0.13", stderr="")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    assert sp.tool_version("miniprot") == "miniprot 0.13"
    assert seen.get("errors") == "replace"


def test_tool_version_decode_error_returns_none(monkeypatch):
    def fake_run(argv, **kw):
        raise UnicodeDecodeError("utf-8", b"\xab", 0, 1, "invalid start byte")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    assert sp.tool_version("weird") is None


def test_tool_version_custom_arg(monkeypatch):
    seen = {}

    def fake_run(argv, **kw):
        seen["argv"] = argv
        return _FakeProc(0, stdout="v1", stderr="")

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    assert sp.tool_version("samtools", version_arg="version") == "v1"
    assert seen["argv"] == ["samtools", "version"]


class _FakePipe:
    def close(self):
        pass


def test_run_pipe_builds_pipeline_and_checks_returncodes(monkeypatch):
    created = []

    class FakePopen:
        def __init__(self, argv, **kw):
            self.argv = argv
            self.kw = kw
            self.returncode = 0
            self.stdout = _FakePipe()
            created.append(self)

        def communicate(self):
            return ("", "")

    monkeypatch.setattr(sp.subprocess, "Popen", FakePopen)
    sp.run_pipe([["hisat2", "-p", 4], ["samtools", "sort", "-"]])
    assert len(created) == 2
    assert created[0].argv == ["hisat2", "-p", "4"]
    assert created[1].argv == ["samtools", "sort", "-"]
    # second command receives the first's stdout as stdin
    assert created[1].kw["stdin"] is created[0].stdout


def test_run_pipe_nonzero_raises(monkeypatch):
    class FakePopen:
        def __init__(self, argv, **kw):
            self.argv = argv
            self.returncode = 0 if argv[0] == "hisat2" else 1
            self.stdout = _FakePipe()

        def communicate(self):
            return ("", "sort failed")

    monkeypatch.setattr(sp.subprocess, "Popen", FakePopen)
    with pytest.raises(RuntimeError) as exc:
        sp.run_pipe([["hisat2"], ["samtools", "sort"]])
    assert "samtools" in str(exc.value)
    assert "sort failed" in str(exc.value)
