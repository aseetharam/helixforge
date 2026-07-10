"""Tests for mikado/config.py (Phase 4). Floor: 12. Subprocess mocked."""

import pytest

from helixforge.mikado import config as config_mod
from helixforge.mikado.config import (
    install_scoring_profile,
    list_scoring_profiles,
    validate_scoring_profile,
    write_configuration,
    write_input_list,
)


def _parse_list(path):
    return [line.rstrip("\n").split("\t") for line in open(path) if line.strip()]


# --- write_input_list ---

def test_input_list_returns_path(tmp_path):
    entries = [{"file": "helixer.gtf", "label": "helixer", "is_reference": True}]
    out = write_input_list(entries, tmp_path / "list.txt")
    assert out.exists()


def test_input_list_helixer_is_reference(tmp_path):
    entries = [
        {"file": "helixer.gtf", "label": "helixer", "is_reference": True},
        {"file": "stringtie.gtf", "label": "st"},
    ]
    out = write_input_list(entries, tmp_path / "list.txt")
    rows = _parse_list(out)
    helixer_row = next(r for r in rows if r[1] == "helixer")
    st_row = next(r for r in rows if r[1] == "st")
    # column order: file,label,strand_specific,score,is_reference,exclude_redundant
    assert helixer_row[4] == "True"
    assert st_row[4] == "False"


def test_input_list_defaults(tmp_path):
    entries = [{"file": "st.gtf", "label": "st"}]
    out = write_input_list(entries, tmp_path / "list.txt")
    row = _parse_list(out)[0]
    assert row[2] == "True"   # strand_specific default
    assert row[3] == "0"      # score default
    assert row[5] == "False"  # exclude_redundant default


def test_input_list_exclude_redundant(tmp_path):
    entries = [{"file": "st.gtf", "label": "st", "exclude_redundant": True}]
    out = write_input_list(entries, tmp_path / "list.txt")
    assert _parse_list(out)[0][5] == "True"


# --- write_configuration (template mode) ---

def test_config_sets_as_knobs(tmp_path):
    out = write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                              tmp_path / "config.yaml")
    text = out.read_text()
    assert "report: true" in text
    assert "only_confirmed_introns: true" in text
    assert "pad: true" in text
    assert "keep_retained_introns: false" in text
    assert "max_isoforms: 5" in text


def test_config_sets_flank_and_chimera(tmp_path):
    out = write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                              tmp_path / "config.yaml")
    text = out.read_text()
    assert "flank: 200" in text
    assert "execute: true" in text


def test_config_references_inputs(tmp_path):
    out = write_configuration("genome.fa", "mylist.txt", "myscore.yaml", "myjunc.bed",
                              tmp_path / "config.yaml")
    text = out.read_text()
    assert "genome.fa" in text
    assert "myscore.yaml" in text
    assert "myjunc.bed" in text
    assert "mylist.txt" in text


def test_config_valid_ccodes_listed(tmp_path):
    out = write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                              tmp_path / "config.yaml")
    text = out.read_text()
    assert "valid_ccodes:" in text
    assert "- j" in text


def test_config_pick_override_max_isoforms(tmp_path):
    out = write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                              tmp_path / "config.yaml", max_isoforms=10)
    assert "max_isoforms: 10" in out.read_text()


def test_config_override_flank(tmp_path):
    out = write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                              tmp_path / "config.yaml", flank=500)
    assert "flank: 500" in out.read_text()


def test_config_override_keep_retained_introns(tmp_path):
    out = write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                              tmp_path / "config.yaml", keep_retained_introns=True)
    assert "keep_retained_introns: true" in out.read_text()


def test_config_subprocess_mode_argv(monkeypatch, tmp_path):
    monkeypatch.setattr(config_mod.shutil, "which", lambda b: "/usr/bin/mikado")
    captured = {}
    monkeypatch.setattr(config_mod.subprocess, "run",
                        lambda argv, check, **kw: captured.update(argv=argv, check=check, **kw))
    out = write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                              tmp_path / "config.yaml", use_subprocess=True)
    argv = captured["argv"]
    assert argv[0] == "mikado"
    assert argv[1] == "configure"
    assert "--list" in argv and "--reference" in argv
    assert "--scoring" in argv and "--junctions" in argv


def test_config_subprocess_missing_binary_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(config_mod.shutil, "which", lambda b: None)
    with pytest.raises(FileNotFoundError):
        write_configuration("g.fa", "list.txt", "score.yaml", "j.bed",
                            tmp_path / "config.yaml", use_subprocess=True)


# --- install_scoring_profile ---

def test_install_scoring_profile_copies(tmp_path):
    out = install_scoring_profile("strict", dest_dir=tmp_path)
    assert out.exists()
    text = out.read_text()
    assert "requirements:" in text
    assert "scoring:" in text


def test_install_scoring_profile_external_metrics_present(tmp_path):
    out = install_scoring_profile("strict", dest_dir=tmp_path)
    assert "external.helixer_support" in out.read_text()


def test_install_scoring_profile_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="unknown scoring profile"):
        install_scoring_profile("nonexistent", dest_dir=tmp_path)


# --- list_scoring_profiles ---

def test_list_scoring_profiles_returns_known():
    profiles = list_scoring_profiles()
    assert "strict" in profiles
    assert "permissive" in profiles


def test_list_scoring_profiles_sorted():
    profiles = list_scoring_profiles()
    assert profiles == sorted(profiles)


# --- validate_scoring_profile ---

def test_validate_scoring_profile_strict():
    validate_scoring_profile("strict")


def test_validate_scoring_profile_permissive():
    validate_scoring_profile("permissive")


def test_validate_scoring_profile_unknown_raises():
    with pytest.raises(FileNotFoundError, match="unknown scoring profile"):
        validate_scoring_profile("nonexistent_xyz")


def test_validate_scoring_profile_unknown_lists_valid():
    with pytest.raises(FileNotFoundError, match="strict") as exc_info:
        validate_scoring_profile("nonexistent_xyz")
    assert "permissive" in str(exc_info.value)


# --- install_scoring_profile (permissive) ---

def test_install_scoring_profile_permissive(tmp_path):
    out = install_scoring_profile("permissive", dest_dir=tmp_path)
    assert out.exists()
    text = out.read_text()
    assert "requirements:" in text
    assert "scoring:" in text


def test_install_scoring_profile_permissive_external_metrics(tmp_path):
    out = install_scoring_profile("permissive", dest_dir=tmp_path)
    assert "external.helixer_support" in out.read_text()


# --- _read_scoring_template does not use a malformed docs/scoring path ---

def test_no_malformed_docs_scoring_path():
    """The legacy ``_SCORING_DIR`` fallback to ``docs/scoring`` via ``parents[3]``
    is removed. Only importlib.resources resolution should exist."""
    import inspect
    source = inspect.getsource(config_mod._read_scoring_template)
    assert "parents[" not in source
    assert "docs" not in source
