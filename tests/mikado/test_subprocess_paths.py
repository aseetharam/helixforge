"""Tests for absolute-path enforcement in external-tool subprocess invocations.

Every path argument passed to an external tool must be absolute.  Every
subprocess must run with an explicit cwd.  These tests inspect the argv and
cwd that the command builders would use — no real tools are invoked.
"""

import os
import subprocess
from pathlib import Path

import pytest

from helixforge.mikado import run as run_mod
from helixforge.mikado.run import (
    run_diamond,
    run_pick,
    run_prepare,
    run_serialise,
    run_transdecoder,
)


@pytest.fixture
def capture(monkeypatch):
    """Record every subprocess.run call (argv + cwd); return a list of dicts."""
    calls: list[dict] = []

    def fake_run(argv, check, capture_output=True, text=True, cwd=None, **kw):
        calls.append({"argv": argv, "cwd": cwd})
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(run_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(run_mod.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(
        run_mod, "tool_version", lambda b, a="--version": "mikado 2.3.4"
    )
    return calls


def _all_path_args(argv: list[str]) -> list[str]:
    """Extract every argv element that looks like a file/dir path.

    Skips the binary (argv[0]), subcommands, bare flags, and pure numbers.
    """
    skip_next = False
    paths = []
    for i, a in enumerate(argv):
        if skip_next:
            skip_next = False
            continue
        if i == 0:
            continue
        # skip flags and their values that are not paths
        if a.startswith("-"):
            # flags whose next arg is a non-path value
            if a in ("--procs", "--threads", "--evalue", "--max-target-seqs",
                      "--outfmt", "-p", "-t"):
                skip_next = True
            continue
        # skip bare numbers (procs, threads, evalue etc)
        try:
            float(a)
            continue
        except ValueError:
            pass
        # skip mikado subcommands and format column names
        if a in ("prepare", "serialise", "pick", "configure", "makedb",
                  "blastx", "qseqid", "sseqid", "pident", "length",
                  "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                  "evalue", "bitscore", "BAM", "SortedByCoordinate",
                  "--single_best_only"):
            continue
        paths.append(a)
    return paths


# --- run_transdecoder: the original bug ---


def test_transdecoder_output_dir_is_absolute(capture, tmp_path):
    """The bug: a relative out_dir + cwd=out_dir caused TransDecoder to nest."""
    rel_dir = "some_relative/work/mikado_run"
    # Use a relative path — the function must resolve it
    bed = run_transdecoder(tmp_path / "mikado_prepared.fasta", rel_dir)
    for call in capture:
        assert call["cwd"] is not None, "TransDecoder must have explicit cwd"
        assert os.path.isabs(call["cwd"]), (
            f"TransDecoder cwd must be absolute, got {call['cwd']!r}"
        )
        for p in _all_path_args(call["argv"]):
            assert os.path.isabs(p), (
                f"TransDecoder argv path must be absolute, got {p!r}"
            )


def test_transdecoder_creates_output_dir(capture, tmp_path):
    """The pipeline creates the output dir before invoking TransDecoder."""
    out_dir = tmp_path / "td_out"
    run_transdecoder(tmp_path / "p.fasta", out_dir)
    assert out_dir.exists(), "TransDecoder output dir must be pre-created"


# --- run_prepare ---


def test_prepare_absolute_paths_and_cwd(capture, tmp_path):
    run_prepare("config.yaml", tmp_path / "out")
    for call in capture:
        assert call["cwd"] is not None, "prepare must have explicit cwd"
        assert os.path.isabs(call["cwd"])
        for p in _all_path_args(call["argv"]):
            assert os.path.isabs(p), f"prepare argv path not absolute: {p!r}"


# --- run_diamond ---


def test_diamond_absolute_paths_and_cwd(capture, tmp_path):
    run_diamond(tmp_path / "p.fasta", "protein.fa", tmp_path / "out")
    for call in capture:
        assert call["cwd"] is not None, "diamond must have explicit cwd"
        assert os.path.isabs(call["cwd"])
        for p in _all_path_args(call["argv"]):
            assert os.path.isabs(p), f"diamond argv path not absolute: {p!r}"


# --- run_serialise ---


def test_serialise_absolute_paths_and_cwd(capture, tmp_path):
    run_serialise(
        "config.yaml", "p.fasta", "orfs.bed", "blast.xml", "prot.fa",
        "j.bed", "ext.tsv", "genome.fa", tmp_path / "out",
    )
    for call in capture:
        assert call["cwd"] is not None, "serialise must have explicit cwd"
        assert os.path.isabs(call["cwd"])
        for p in _all_path_args(call["argv"]):
            assert os.path.isabs(p), f"serialise argv path not absolute: {p!r}"


# --- run_pick ---


def test_pick_absolute_paths_and_cwd(capture, tmp_path):
    run_pick("config.yaml", "scoring.yaml", "prepared.gtf", tmp_path / "out")
    for call in capture:
        assert call["cwd"] is not None, "pick must have explicit cwd"
        assert os.path.isabs(call["cwd"])
        for p in _all_path_args(call["argv"]):
            assert os.path.isabs(p), f"pick argv path not absolute: {p!r}"


# --- relative out_dir regression: all functions resolve it ---


@pytest.mark.parametrize("func,args", [
    ("prepare", ("config.yaml", "rel_work/out")),
    ("transdecoder", ("rel_work/p.fasta", "rel_work/out")),
    ("diamond", ("rel_work/p.fasta", "prot.fa", "rel_work/out")),
    ("pick", ("config.yaml", "scoring.yaml", "prepared.gtf", "rel_work/out")),
])
def test_relative_out_dir_resolved(capture, func, args):
    """A relative out_dir is silently resolved to absolute — never passed raw."""
    fn = getattr(run_mod, f"run_{func}")
    fn(*args)
    for call in capture:
        assert call["cwd"] is not None
        assert os.path.isabs(call["cwd"]), (
            f"{func}: cwd must be absolute, got {call['cwd']!r}"
        )
