"""Tests for the broadened preflight: work-dir writability + tool-chain checks.

An unwritable or uncreatable work dir, or a missing tool in the configured
chain, must fail **before** the expensive Mikado stage (seconds, not hours).
"""

import os
import stat
from types import SimpleNamespace

import pytest

from helixforge.prep.preflight import (
    PreflightReport,
    check_tool_chain,
    check_work_dir,
    run_preflight,
)


# ---------------------------------------------------------------------------
# check_work_dir
# ---------------------------------------------------------------------------


def test_workdir_created_and_writable(tmp_path):
    """A valid work_dir is created with stage subdirs; no errors."""
    cfg = SimpleNamespace(work_dir=str(tmp_path / "new_work"))
    errors = check_work_dir(cfg)
    assert errors == []
    assert (tmp_path / "new_work").is_dir()
    assert (tmp_path / "new_work" / "mikado_inputs").is_dir()
    assert (tmp_path / "new_work" / "mikado_run").is_dir()


def test_workdir_no_attr():
    """A config without work_dir (e.g. doctor's SimpleNamespace) is fine."""
    cfg = SimpleNamespace()
    assert check_work_dir(cfg) == []


def test_workdir_none():
    """work_dir=None is fine (not configured yet)."""
    cfg = SimpleNamespace(work_dir=None)
    assert check_work_dir(cfg) == []


@pytest.mark.skipif(os.getuid() == 0, reason="root bypasses permissions")
def test_workdir_unwritable(tmp_path):
    """An unwritable parent dir fails preflight with a clear message."""
    blocked = tmp_path / "blocked"
    blocked.mkdir()
    blocked.chmod(stat.S_IRUSR | stat.S_IXUSR)  # r-x only
    try:
        cfg = SimpleNamespace(work_dir=str(blocked / "work"))
        errors = check_work_dir(cfg)
        assert len(errors) >= 1
        assert "cannot create" in errors[0].lower() or "not writable" in errors[0].lower()
    finally:
        blocked.chmod(stat.S_IRWXU)


@pytest.mark.skipif(os.getuid() == 0, reason="root bypasses permissions")
def test_workdir_unwritable_existing(tmp_path):
    """An existing but read-only work dir fails the write probe."""
    work = tmp_path / "readonly_work"
    work.mkdir()
    work.chmod(stat.S_IRUSR | stat.S_IXUSR)  # r-x only
    try:
        cfg = SimpleNamespace(work_dir=str(work))
        errors = check_work_dir(cfg)
        assert len(errors) >= 1
        assert any("not writable" in e.lower() or "cannot create" in e.lower()
                    for e in errors)
    finally:
        work.chmod(stat.S_IRWXU)


# ---------------------------------------------------------------------------
# check_tool_chain
# ---------------------------------------------------------------------------


def test_tool_chain_basic_mode():
    """No stringtie + no protein_db → basic mode → no tool requirements."""
    cfg = SimpleNamespace(
        stringtie_list=[],
        protein_db=None,
        mikado_bin="mikado",
        diamond_bin="diamond",
        transdecoder_bin_dir=None,
        portcullis_bin=None,
        bam_paths=[],
    )
    assert check_tool_chain(cfg) == []


def test_tool_chain_missing_mikado(monkeypatch):
    """Missing mikado in full mode fails preflight."""
    monkeypatch.setattr("shutil.which", lambda name: None)
    cfg = SimpleNamespace(
        stringtie_list=["sample.gtf"],
        protein_db="proteins.fa",
        mikado_bin="mikado",
        diamond_bin="diamond",
        transdecoder_bin_dir=None,
        portcullis_bin=None,
        bam_paths=[],
    )
    errors = check_tool_chain(cfg)
    assert len(errors) >= 1
    tool_names = " ".join(errors)
    assert "mikado" in tool_names


def test_tool_chain_all_present(monkeypatch):
    """All tools on PATH → no errors."""
    monkeypatch.setattr("shutil.which", lambda name: f"/usr/bin/{name}")
    cfg = SimpleNamespace(
        stringtie_list=["sample.gtf"],
        protein_db="proteins.fa",
        mikado_bin="mikado",
        diamond_bin="diamond",
        transdecoder_bin_dir=None,
        portcullis_bin=None,
        bam_paths=[],
    )
    assert check_tool_chain(cfg) == []


def test_tool_chain_custom_transdecoder_bin_dir(monkeypatch):
    """A custom transdecoder_bin_dir resolves the full path."""
    found = set()
    def fake_which(name):
        found.add(name)
        return f"/usr/bin/{name}"
    monkeypatch.setattr("shutil.which", fake_which)
    cfg = SimpleNamespace(
        stringtie_list=["sample.gtf"],
        protein_db="proteins.fa",
        mikado_bin="mikado",
        diamond_bin="diamond",
        transdecoder_bin_dir="/opt/td/bin",
        portcullis_bin=None,
        bam_paths=[],
    )
    check_tool_chain(cfg)
    assert "/opt/td/bin/TransDecoder.LongOrfs" in found


def test_tool_chain_portcullis_checked_when_configured(monkeypatch):
    """Portcullis is checked only when portcullis_bin + BAMs are configured."""
    missing = set()
    def fake_which(name):
        if name == "missing_portcullis":
            return None
        return f"/usr/bin/{name}"
    monkeypatch.setattr("shutil.which", fake_which)

    # With portcullis + BAMs → error
    cfg = SimpleNamespace(
        stringtie_list=[],
        protein_db=None,
        mikado_bin="mikado",
        diamond_bin="diamond",
        transdecoder_bin_dir=None,
        portcullis_bin="missing_portcullis",
        bam_paths=["sample.bam"],
    )
    errors = check_tool_chain(cfg)
    assert any("missing_portcullis" in e for e in errors)

    # Without BAMs → no portcullis check
    cfg2 = SimpleNamespace(
        stringtie_list=[],
        protein_db=None,
        mikado_bin="mikado",
        diamond_bin="diamond",
        transdecoder_bin_dir=None,
        portcullis_bin="missing_portcullis",
        bam_paths=[],
    )
    assert check_tool_chain(cfg2) == []


# ---------------------------------------------------------------------------
# run_preflight integration
# ---------------------------------------------------------------------------


def test_preflight_report_includes_workdir_and_tools(tmp_path):
    """The combined preflight report surfaces work-dir and tool errors."""
    cfg = SimpleNamespace(
        genome_fasta=None,
        helixer_gff3=None,
        helixer_h5=None,
        stringtie_list=[],
        bam_paths=[],
        star_sj_paths=[],
        miniprot_gff=None,
        work_dir=str(tmp_path / "ok_work"),
    )
    report = run_preflight(cfg)
    assert isinstance(report, PreflightReport)
    assert report.workdir_errors == []
    assert report.tool_errors == []


def test_preflight_report_renders_workdir_error(tmp_path, monkeypatch):
    """An unwritable work_dir appears in the rendered report."""
    # Force check_work_dir to return an error
    monkeypatch.setattr(
        "helixforge.prep.preflight.check_work_dir",
        lambda cfg: ["cannot create stage dir /bad/path: Permission denied"],
    )
    cfg = SimpleNamespace(
        genome_fasta=None,
        helixer_gff3=None,
        helixer_h5=None,
        stringtie_list=[],
        bam_paths=[],
        star_sj_paths=[],
        miniprot_gff=None,
        work_dir="/bad/path",
    )
    report = run_preflight(cfg)
    assert not report.ok()
    rendered = report.render()
    assert "work-dir" in rendered
    assert "Permission denied" in rendered
