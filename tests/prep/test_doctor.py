"""Tests for the doctor preflight (Phase 13 D3). Floor: 8. tool_version mocked.

Found / missing / mismatch rendering; required-missing ⇒ not OK (CLI nonzero);
bin overrides via dict and via a config-like object.
"""

import pytest
from click.testing import CliRunner

from helixforge.cli import main
from helixforge.prep import doctor as doc
from helixforge.prep.doctor import FOUND, MISMATCH, MISSING, check_tools

# bin → a version string that satisfies its pin (pin-less tools take any).
PIN_BY_BIN = {s.default_bin: (s.pinned or "1.0") for s in doc.TOOL_MATRIX}


def _all_found(bin_name, arg="--version"):
    return PIN_BY_BIN.get(bin_name, "1.0")


def _status(report, key):
    return next(s for s in report.statuses if s.key == key)


def test_all_found(monkeypatch):
    monkeypatch.setattr(doc, "tool_version", _all_found)
    report = check_tools()
    assert report.ok()
    assert all(s.status == FOUND for s in report.statuses)


def test_required_missing_not_ok(monkeypatch):
    def mock(bin_name, arg="--version"):
        return None if bin_name == "mikado" else _all_found(bin_name)

    monkeypatch.setattr(doc, "tool_version", mock)
    report = check_tools()
    assert not report.ok()
    assert _status(report, "mikado").status == MISSING


def test_required_missing_listed(monkeypatch):
    monkeypatch.setattr(
        doc, "tool_version",
        lambda b, a="--version": None if b == "diamond" else _all_found(b),
    )
    report = check_tools()
    missing_keys = {s.key for s in report.missing_required()}
    assert "diamond" in missing_keys


def test_version_mismatch(monkeypatch):
    monkeypatch.setattr(
        doc, "tool_version",
        lambda b, a="--version": "2.0.0" if b == "mikado" else _all_found(b),
    )
    report = check_tools()
    # 2.0.0 != pinned 2.3.4 → mismatch, but mismatch is not "missing" → still OK.
    assert _status(report, "mikado").status == MISMATCH
    assert report.ok()


def test_optional_missing_still_ok(monkeypatch):
    monkeypatch.setattr(
        doc, "tool_version",
        lambda b, a="--version": None if b == "STAR" else _all_found(b),
    )
    report = check_tools()
    assert _status(report, "star").status == MISSING
    assert report.ok()  # STAR is optional


def test_pinless_tool_any_version_found(monkeypatch):
    monkeypatch.setattr(
        doc, "tool_version",
        lambda b, a="--version": "0.12.6" if b == "gffcompare" else _all_found(b),
    )
    report = check_tools()
    assert _status(report, "gffcompare").status == FOUND


def test_render_contains_table(monkeypatch):
    monkeypatch.setattr(
        doc, "tool_version",
        lambda b, a="--version": None if b == "mikado" else _all_found(b),
    )
    text = check_tools().render()
    assert "TOOL" in text and "STATUS" in text
    assert "mikado" in text
    assert MISSING in text


def test_dict_override_bin(monkeypatch):
    seen = {}

    def mock(bin_name, arg="--version"):
        seen[bin_name] = True
        return _all_found("mikado")

    monkeypatch.setattr(doc, "tool_version", mock)
    report = check_tools({"mikado": "/opt/mikado/mikado"})
    assert _status(report, "mikado").bin == "/opt/mikado/mikado"
    assert "/opt/mikado/mikado" in seen


def test_config_object_override_bin(monkeypatch):
    monkeypatch.setattr(doc, "tool_version", lambda b, a="--version": _all_found("mikado"))

    class Cfg:
        mikado_bin = "/custom/mikado"
        diamond_bin = None

    report = check_tools(Cfg())
    assert _status(report, "mikado").bin == "/custom/mikado"


# --- CLI wiring ---


def test_cli_doctor_ok_exit_zero(monkeypatch):
    monkeypatch.setattr(doc, "tool_version", _all_found)
    result = CliRunner().invoke(main, ["doctor"])
    assert result.exit_code == 0
    assert "STATUS" in result.output


def test_cli_doctor_missing_required_nonzero(monkeypatch):
    monkeypatch.setattr(
        doc, "tool_version",
        lambda b, a="--version": None if b == "mikado" else _all_found(b),
    )
    result = CliRunner().invoke(main, ["doctor"])
    assert result.exit_code != 0


# --- emitted-config schema verification (Phase 23 D5) ---


def _write_config(path, sections):
    path.write_text("".join(f"{s}:\n  k: v\n" for s in sections))
    return str(path)


def test_verify_emitted_config_passes_good_schema(tmp_path, monkeypatch):
    # mocked mikado: pinned version detected.
    monkeypatch.setattr(doc, "mikado_version", lambda b="mikado": (2, 3, 4))
    cfg = _write_config(tmp_path / "configuration.yaml",
                        ["reference", "prepare", "serialise", "pick"])
    scoring = _write_config(tmp_path / "scoring.yaml", ["requirements", "scoring"])
    v = doc.verify_emitted_config(cfg, scoring_path=scoring)
    assert v.ok is True
    assert v.config_missing == [] and v.scoring_missing == []
    assert v.warnings == []                       # pinned version → no drift warning


def test_verify_emitted_config_flags_mismatched_schema(tmp_path, monkeypatch):
    monkeypatch.setattr(doc, "mikado_version", lambda b="mikado": (2, 3, 4))
    # deliberately missing the `pick` + `serialise` sections.
    cfg = _write_config(tmp_path / "configuration.yaml", ["reference", "prepare"])
    v = doc.verify_emitted_config(cfg)
    assert v.ok is False
    assert "pick" in v.config_missing and "serialise" in v.config_missing


def test_verify_emitted_config_warns_on_version_drift(tmp_path):
    cfg = _write_config(tmp_path / "configuration.yaml",
                        ["reference", "prepare", "serialise", "pick"])
    # an undetectable version is a loud warning but not a hard failure.
    v = doc.verify_emitted_config(cfg, version=None, mikado_bin="/nonexistent/mikado")
    assert v.ok is True
    assert any("version" in w for w in v.warnings)


def test_cli_doctor_check_config_bad_exits_nonzero(tmp_path, monkeypatch):
    monkeypatch.setattr(doc, "tool_version", _all_found)
    monkeypatch.setattr(doc, "mikado_version", lambda b="mikado": (2, 3, 4))
    cfg = _write_config(tmp_path / "configuration.yaml", ["reference"])  # missing 3
    result = CliRunner().invoke(main, ["doctor", "--check-config", cfg])
    assert result.exit_code != 0
    assert "schema" in result.output.lower()


# --- biology.v4 Phase 26 D4: tool-version drift warnings ---


def test_report_drift_warnings_for_mikado(monkeypatch):
    # Mikado resolves to a non-pinned version → MISMATCH → a drift warning that
    # calls out the scoring/config schema risk (CLAUDE.md §7).
    monkeypatch.setattr(
        doc, "tool_version",
        lambda b, a="--version": "2.0.0" if b == "mikado" else _all_found(b),
    )
    report = check_tools()
    warnings = report.drift_warnings()
    assert any(w.startswith("mikado:") and "schema" in w for w in warnings)
    assert "WARNING: version drift" in report.render()
    # drift is a warning, not a hard failure.
    assert report.ok()


def test_no_drift_warnings_when_all_pinned(monkeypatch):
    monkeypatch.setattr(doc, "tool_version", _all_found)
    report = check_tools()
    assert report.drift_warnings() == []


# --- Phase 33 D4: environment hygiene (shim + CRAM-reference readiness) ------

from helixforge.prep.doctor import check_cram_reference, check_shim


def _make_checkout(tmp_path):
    """A directory that looks like a helixforge checkout (has src/helixforge)."""
    pkg = tmp_path / "src" / "helixforge"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    return tmp_path, pkg


def test_check_shim_ok_when_module_in_checkout(tmp_path):
    root, pkg = _make_checkout(tmp_path)
    status = check_shim(root, module_file=str(pkg / "__init__.py"),
                        which=lambda n: "/some/bin/helixforge")
    assert status.applicable and not status.foreign
    assert "OK" in status.render()


def test_check_shim_foreign_when_module_outside_checkout(tmp_path):
    root, _pkg = _make_checkout(tmp_path)
    foreign = tmp_path / "other" / "helixforge"
    foreign.mkdir(parents=True)
    status = check_shim(root, module_file=str(foreign / "__init__.py"),
                        which=lambda n: "/conda/mikado-env/bin/helixforge")
    assert status.applicable and status.foreign
    assert "FOREIGN" in status.render()
    assert status.console_script == "/conda/mikado-env/bin/helixforge"


def test_check_shim_not_applicable_outside_a_checkout(tmp_path):
    # tmp_path has no src/helixforge → the check cannot apply (no false positive).
    status = check_shim(tmp_path, module_file=str(tmp_path / "helixforge" / "__init__.py"),
                        which=lambda n: None)
    assert not status.applicable
    assert not status.foreign
    assert "n/a" in status.render()


def test_check_cram_reference_ok_with_reference():
    status = check_cram_reference(["reads.cram"], reference="genome.fa", env={})
    assert status.has_cram and status.ok
    assert "OK" in status.render()


def test_check_cram_reference_missing_without_ref_or_cache():
    status = check_cram_reference(["reads.cram"], reference=None, env={})
    assert status.has_cram and not status.ok
    assert "MISSING" in status.render()


def test_check_cram_reference_ok_via_ref_cache_env():
    status = check_cram_reference(["reads.cram"], reference=None,
                                  env={"REF_CACHE": "/scratch/ref/%2s/%2s/%s"})
    assert status.ok and status.ref_cache_set
    assert "REF_CACHE" in status.render()


def test_check_cram_reference_na_for_plain_bam():
    status = check_cram_reference(["reads.bam"], reference=None, env={})
    assert not status.has_cram and status.ok  # no CRAM → nothing to warn about
    assert "n/a" in status.render()


# --- check_scoring_profiles ---

def test_check_scoring_profiles_ok():
    from helixforge.prep.doctor import check_scoring_profiles

    status = check_scoring_profiles()
    assert status.ok
    assert "strict" in status.profiles
    assert "permissive" in status.profiles
    assert all(status.profiles.values())


def test_check_scoring_profiles_render():
    from helixforge.prep.doctor import check_scoring_profiles

    status = check_scoring_profiles()
    rendered = status.render()
    assert "strict: OK" in rendered
    assert "permissive: OK" in rendered
