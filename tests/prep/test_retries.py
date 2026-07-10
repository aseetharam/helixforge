"""Phase 24 D5 — bounded subprocess retries. Floor: 2.

An idempotent step (alignment / StringTie / DIAMOND) is retried with backoff on
an **allowlisted** transient exit; a non-idempotent step (Mikado ``serialise``)
is **never** retried (rely on Phase 22 checkpointing). Subprocess is mocked; the
``_sleep`` seam keeps the backoff instant.
"""

import subprocess

import pytest

from helixforge.mikado import run as mrun
from helixforge.prep import _subprocess as sp


def _fail(code, argv):
    return subprocess.CalledProcessError(code, argv, stderr="transient FS hiccup")


class _OK:
    returncode = 0
    stderr = ""


# ---------------------------------------------------------------------------
# Idempotent step: retried
# ---------------------------------------------------------------------------

def test_transient_failure_retried_then_succeeds(monkeypatch):
    calls = []

    def fake_run(argv, **kw):
        calls.append(1)
        if len(calls) < 3:          # fail the first two attempts
            raise _fail(1, argv)
        return _OK()

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    proc = sp.run_tool(
        ["aligner", "in"], retries=3, retry_exit_codes={1}, _sleep=lambda s: None
    )
    assert proc.returncode == 0
    assert len(calls) == 3           # 2 transient failures + 1 success


def test_exhausts_retries_then_raises(monkeypatch):
    calls = []

    def fake_run(argv, **kw):
        calls.append(1)
        raise _fail(1, argv)

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    with pytest.raises(sp.ToolError):
        sp.run_tool(["aligner"], retries=2, retry_exit_codes={1}, _sleep=lambda s: None)
    assert len(calls) == 3           # initial attempt + 2 retries, then give up


def test_non_allowlisted_exit_not_retried(monkeypatch):
    calls = []

    def fake_run(argv, **kw):
        calls.append(1)
        raise _fail(42, argv)        # 42 is NOT in the allowlist

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    with pytest.raises(sp.ToolError):
        sp.run_tool(["aligner"], retries=3, retry_exit_codes={1}, _sleep=lambda s: None)
    assert len(calls) == 1           # a non-transient exit is fatal immediately


def test_diamond_blastx_retried(monkeypatch, tmp_path):
    # makedb succeeds; blastx fails once then succeeds → 3 subprocess calls total.
    seq = []

    def fake_run(argv, **kw):
        step = "blastx" if "blastx" in argv else "makedb"
        seq.append(step)
        if step == "blastx" and seq.count("blastx") == 1:
            raise _fail(1, argv)
        return _OK()

    monkeypatch.setattr(mrun.subprocess, "run", fake_run)
    out = mrun.run_diamond(
        "prepared.fa", "db.fa", tmp_path, retries=2, retry_exit_codes={1},
    )
    assert seq == ["makedb", "blastx", "blastx"]
    assert str(out).endswith("mikado_diamond.xml")


# ---------------------------------------------------------------------------
# Non-idempotent step: never retried
# ---------------------------------------------------------------------------

def test_default_run_tool_does_not_retry(monkeypatch):
    calls = []

    def fake_run(argv, **kw):
        calls.append(1)
        raise _fail(1, argv)

    monkeypatch.setattr(sp.subprocess, "run", fake_run)
    with pytest.raises(sp.ToolError):
        sp.run_tool(["tool"])        # retries defaults to 0
    assert len(calls) == 1


def test_serialise_is_never_retried(monkeypatch, tmp_path):
    calls = []

    def fake_run(argv, **kw):
        calls.append(1)
        raise _fail(1, argv)         # serialise DB insertion is non-idempotent

    monkeypatch.setattr(mrun.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError):
        mrun.run_serialise(
            "cfg.yaml", "prepared.fa", "orfs.bed", "blast.xml", "db.fa",
            "junc.tab", "ext.tsv", "genome.fa", tmp_path, version=(2, 3, 4),
        )
    assert len(calls) == 1           # exactly one attempt — no retry
