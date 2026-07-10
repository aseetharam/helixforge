"""Tests for stage checkpoint/resume + --mikado-loci in run_pipeline. Floor: 4.

Two layers: the :class:`Checkpoint` validation contract (unit), and the
``run_pipeline`` resume wiring (the expensive MIKADO stage is skipped + re-parsed
on a valid checkpoint, re-run on a missing/invalid one, and always re-run when
``resume=False``). Plus ``--mikado-loci`` tests: externally-supplied Mikado
output skips the PREP+MIKADO chain entirely.

The pipeline stages are mocked so the tests are fast and assert *which* stage
functions were called.
"""

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from helixforge.reconcile import pipeline as pl
from helixforge.reconcile.pipeline import PipelineConfig
from helixforge.utils.checkpoint import Checkpoint


# ---------------------------------------------------------------------------
# Checkpoint validation contract
# ---------------------------------------------------------------------------

def test_checkpoint_marks_and_validates_existing_outputs(tmp_path):
    out = tmp_path / "loci.gff3"
    out.write_text("##gff-version 3\n")
    cp = Checkpoint(tmp_path / "checkpoint.json", enabled=True)
    assert cp.is_complete("MIKADO", [out]) is False  # not yet marked
    cp.mark("MIKADO", [out])
    # A fresh reader sees the persisted manifest and validates.
    assert Checkpoint(tmp_path / "checkpoint.json", enabled=True).is_complete(
        "MIKADO", [out]
    ) is True


def test_checkpoint_missing_output_triggers_rerun(tmp_path):
    out = tmp_path / "loci.gff3"
    out.write_text("data\n")
    Checkpoint(tmp_path / "checkpoint.json", enabled=True).mark("MIKADO", [out])
    out.unlink()  # output vanished after the checkpoint was written
    cp = Checkpoint(tmp_path / "checkpoint.json", enabled=True)
    assert cp.is_complete("MIKADO", [out]) is False


def test_checkpoint_changed_output_triggers_rerun(tmp_path):
    out = tmp_path / "loci.gff3"
    out.write_text("original-content\n")
    Checkpoint(tmp_path / "checkpoint.json", enabled=True).mark("MIKADO", [out])
    out.write_text("tampered-different-length-content\n")  # content (hash) changed
    cp = Checkpoint(tmp_path / "checkpoint.json", enabled=True)
    assert cp.is_complete("MIKADO", [out]) is False


def test_checkpoint_disabled_is_inert(tmp_path):
    out = tmp_path / "loci.gff3"
    out.write_text("data\n")
    cp = Checkpoint(tmp_path / "checkpoint.json", enabled=False)
    cp.mark("MIKADO", [out])
    assert not (tmp_path / "checkpoint.json").exists()  # no manifest on a cold run
    assert cp.is_complete("MIKADO", [out]) is False


def test_checkpoint_corrupt_manifest_triggers_rerun(tmp_path):
    out = tmp_path / "loci.gff3"
    out.write_text("data\n")
    (tmp_path / "checkpoint.json").write_text("{ this is not valid json")
    cp = Checkpoint(tmp_path / "checkpoint.json", enabled=True)
    assert cp.is_complete("MIKADO", [out]) is False


# ---------------------------------------------------------------------------
# run_pipeline resume wiring
# ---------------------------------------------------------------------------

class _FakeGenome:
    def close(self):
        pass


@pytest.fixture
def mocked_pipeline(monkeypatch, tmp_path):
    """Patch every stage fn so run_pipeline runs fast; return the MIKADO mocks."""
    monkeypatch.setattr(pl, "_load_inputs", lambda config: {
        "genome": _FakeGenome(), "h5": None, "helixer_loci": [],
        "stringtie_agg": {}, "stringtie_all": [],
    })
    monkeypatch.setattr(pl, "classify_loci", lambda *a, **k: {})
    monkeypatch.setattr(pl, "reconcile", lambda *a, **k: ([], {}, []))
    monkeypatch.setattr(pl, "_save_id_map", lambda *a, **k: None)
    monkeypatch.setattr(pl, "_finalize_genes",
                        lambda genes, config, inputs, stats=None: [])

    def fake_write_outputs(genes, config, functional=None):
        # Phase 31 D1: _write_outputs gained an optional `functional` kwarg.
        from pathlib import Path
        Path(f"{config.output_prefix}.gff3").write_text("##gff-version 3\n")

    def fake_write_report(genes, report_path):
        from pathlib import Path
        Path(report_path).write_text("gene_id\n")

    monkeypatch.setattr(pl, "_write_outputs", fake_write_outputs)
    monkeypatch.setattr(pl, "_write_report", fake_write_report)

    run_mikado = MagicMock(return_value=[])
    reparse = MagicMock(return_value=[])
    monkeypatch.setattr(pl, "_run_mikado_stage", run_mikado)
    monkeypatch.setattr(pl, "_reparse_mikado_loci", reparse)
    # The startup preflight checks shutil.which for tools; stub it so the
    # mocked pipeline proceeds without real binaries.
    from helixforge.prep import preflight as _pf
    monkeypatch.setattr(_pf, "check_tool_chain", lambda cfg: [])
    return run_mikado, reparse


def _config(tmp_path, *, resume):
    work = tmp_path / "work"
    return PipelineConfig(
        genome_fasta=str(tmp_path / "g.fa"),
        helixer_gff3=str(tmp_path / "h.gff3"),
        stringtie_list=[str(tmp_path / "s.gtf")],   # + protein_db ⇒ Mikado runs
        protein_db=str(tmp_path / "db.fa"),
        output_prefix=str(tmp_path / "out"),
        report_path=str(tmp_path / "out.report.tsv"),
        id_map_path=str(tmp_path / "id_map.json"),
        work_dir=str(work),
        resume=resume,
    )


def _seed_mikado_checkpoint(config):
    """Create the 3 Mikado loci files + a valid MIKADO checkpoint manifest."""
    from pathlib import Path
    for p in pl._mikado_output_paths(config):
        p = Path(p)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("loci\n")
    Checkpoint(Path(config.work_dir) / "checkpoint.json", enabled=True).mark(
        "MIKADO", pl._mikado_output_paths(config)
    )


def test_resume_valid_checkpoint_skips_mikado(mocked_pipeline, tmp_path):
    run_mikado, reparse = mocked_pipeline
    config = _config(tmp_path, resume=True)
    _seed_mikado_checkpoint(config)
    pl.run_pipeline(config)
    assert run_mikado.call_count == 0          # external chain skipped
    assert reparse.call_count == 1             # loci re-parsed from disk


def test_resume_missing_checkpoint_reruns_mikado(mocked_pipeline, tmp_path):
    run_mikado, reparse = mocked_pipeline
    config = _config(tmp_path, resume=True)   # no checkpoint seeded
    pl.run_pipeline(config)
    assert run_mikado.call_count == 1
    assert reparse.call_count == 0


def test_resume_false_always_reruns_mikado(mocked_pipeline, tmp_path):
    run_mikado, reparse = mocked_pipeline
    config = _config(tmp_path, resume=False)
    _seed_mikado_checkpoint(config)            # valid checkpoint present...
    pl.run_pipeline(config)
    assert run_mikado.call_count == 1          # ...but resume=False ignores it
    assert reparse.call_count == 0


def test_cold_run_writes_no_checkpoint(mocked_pipeline, tmp_path):
    config = _config(tmp_path, resume=False)
    pl.run_pipeline(config)
    assert not (Path(config.work_dir) / "checkpoint.json").exists()


# ---------------------------------------------------------------------------
# --mikado-loci: external Mikado output
# ---------------------------------------------------------------------------

# Minimal valid Mikado loci GFF3 that parse_loci_gff3 can consume.
_MIKADO_GFF3 = """\
##gff-version 3
chr1\tmikado\tgene\t1001\t2000\t.\t+\t.\tID=mikado.chr1G1
chr1\tmikado\tmRNA\t1001\t2000\t.\t+\t.\tID=mikado.chr1G1.1;Parent=mikado.chr1G1
chr1\tmikado\texon\t1001\t1500\t.\t+\t.\tID=mikado.chr1G1.1.exon1;Parent=mikado.chr1G1.1
chr1\tmikado\texon\t1700\t2000\t.\t+\t.\tID=mikado.chr1G1.1.exon2;Parent=mikado.chr1G1.1
"""

_MIKADO_METRICS = "tid\tblast_score\thas_start_codon\thas_stop_codon\n"
_MIKADO_SCORES = "tid\tscore\n"


def _write_external_loci(tmp_path, *, with_metrics=True, with_scores=True):
    """Write a minimal external Mikado loci GFF3 (+companions) into tmp_path."""
    gff3 = tmp_path / "external_mikado.loci.gff3"
    gff3.write_text(_MIKADO_GFF3)
    if with_metrics:
        (tmp_path / "external_mikado.loci.metrics.tsv").write_text(_MIKADO_METRICS)
    if with_scores:
        (tmp_path / "external_mikado.loci.scores.tsv").write_text(_MIKADO_SCORES)
    return str(gff3)


def _config_external(tmp_path, *, mikado_loci_gff3, resume=False):
    return PipelineConfig(
        genome_fasta=str(tmp_path / "g.fa"),
        helixer_gff3=str(tmp_path / "h.gff3"),
        output_prefix=str(tmp_path / "out"),
        report_path=str(tmp_path / "out.report.tsv"),
        id_map_path=str(tmp_path / "id_map.json"),
        work_dir=str(tmp_path / "work"),
        mikado_loci_gff3=mikado_loci_gff3,
        resume=resume,
    )


def test_mikado_loci_skips_mikado_stage(mocked_pipeline, tmp_path):
    """--mikado-loci supplied → _run_mikado_stage never called."""
    run_mikado, reparse = mocked_pipeline
    gff3 = _write_external_loci(tmp_path)
    config = _config_external(tmp_path, mikado_loci_gff3=gff3)
    pl.run_pipeline(config)
    assert run_mikado.call_count == 0
    assert reparse.call_count == 0


def test_mikado_loci_parses_external_gff3(tmp_path):
    """_parse_external_mikado_loci returns MikadoLocus objects."""
    gff3 = _write_external_loci(tmp_path)
    config = _config_external(tmp_path, mikado_loci_gff3=gff3)
    loci = pl._parse_external_mikado_loci(config)
    assert len(loci) == 1
    assert loci[0].locus_id == "mikado.chr1G1"
    assert len(loci[0].transcripts) == 1


def test_mikado_loci_missing_gff3_raises(tmp_path):
    config = _config_external(
        tmp_path, mikado_loci_gff3=str(tmp_path / "nonexistent.gff3")
    )
    with pytest.raises(FileNotFoundError, match="not found"):
        pl._parse_external_mikado_loci(config)


def test_mikado_loci_empty_gff3_raises(tmp_path):
    gff3 = tmp_path / "empty.loci.gff3"
    gff3.write_text("")
    config = _config_external(tmp_path, mikado_loci_gff3=str(gff3))
    with pytest.raises(ValueError, match="empty"):
        pl._parse_external_mikado_loci(config)


def test_mikado_loci_without_companions_warns(tmp_path, caplog):
    """Missing metrics/scores TSVs log warnings but do not crash."""
    gff3 = _write_external_loci(tmp_path, with_metrics=False, with_scores=False)
    config = _config_external(tmp_path, mikado_loci_gff3=gff3)
    import logging
    with caplog.at_level(logging.WARNING):
        loci = pl._parse_external_mikado_loci(config)
    assert len(loci) == 1
    assert "metrics TSV not found" in caplog.text
    assert "scores TSV not found" in caplog.text


def test_mikado_loci_bypasses_should_run_mikado_gate(mocked_pipeline, tmp_path):
    """--mikado-loci works even without StringTie + protein_db (basic mode)."""
    run_mikado, reparse = mocked_pipeline
    gff3 = _write_external_loci(tmp_path)
    config = PipelineConfig(
        genome_fasta=str(tmp_path / "g.fa"),
        helixer_gff3=str(tmp_path / "h.gff3"),
        # No stringtie_list, no protein_db → _should_run_mikado returns False
        output_prefix=str(tmp_path / "out"),
        report_path=str(tmp_path / "out.report.tsv"),
        id_map_path=str(tmp_path / "id_map.json"),
        work_dir=str(tmp_path / "work"),
        mikado_loci_gff3=gff3,
    )
    pl.run_pipeline(config)
    assert run_mikado.call_count == 0
    assert reparse.call_count == 0


# ---------------------------------------------------------------------------
# PREP-stage resume
# ---------------------------------------------------------------------------

def test_prep_resume_skips_when_outputs_exist(monkeypatch, tmp_path):
    """When config.resume=True and PREP outputs exist, PREP is not re-run.

    We check that helixer_gff3_to_gtf (PREP's first action) is not called when
    the configuration.yaml + list.txt + a scoring YAML already exist in
    mikado_inputs/.
    """
    from helixforge.reconcile.pipeline import _run_mikado_stage

    work = tmp_path / "work"
    mik_in = work / "mikado_inputs"
    mik_in.mkdir(parents=True)
    (mik_in / "configuration.yaml").write_text("pick:\n  scoring_file: s.yaml\n")
    (mik_in / "list.txt").write_text("helixer.gtf\n")
    (mik_in / "scoring_strict.yaml").write_text("requirements:\nscoring:\n")
    (mik_in / "junctions.bed").write_text("")
    mik_run = work / "mikado_run"
    mik_run.mkdir(parents=True)

    config = PipelineConfig(
        genome_fasta=str(tmp_path / "g.fa"),
        helixer_gff3=str(tmp_path / "h.gff3"),
        stringtie_list=[str(tmp_path / "s.gtf")],
        protein_db=str(tmp_path / "db.fa"),
        work_dir=str(work),
        resume=True,
    )

    gtf_calls = []
    monkeypatch.setattr(pl, "helixer_gff3_to_gtf", lambda *a, **k: gtf_calls.append(1))
    monkeypatch.setattr(
        pl, "run_prepare", lambda *a, **k: (str(mik_run / "p.gtf"), str(mik_run / "p.fa"))
    )
    monkeypatch.setattr(pl, "run_transdecoder", lambda *a, **k: str(mik_run / "orfs.bed"))
    monkeypatch.setattr(pl, "run_diamond", lambda *a, **k: str(mik_run / "blast.out"))
    monkeypatch.setattr(pl, "run_serialise", lambda *a, **k: None)
    monkeypatch.setattr(pl, "run_pick", lambda *a, **k: str(mik_run / "mikado.loci.gff3"))
    monkeypatch.setattr(pl, "parse_loci_gff3", lambda *a, **k: [])
    monkeypatch.setattr(pl, "_build_external_scores", lambda *a, **k: {})
    monkeypatch.setattr(pl, "write_external_scores_tsv", lambda *a, **k: mik_in / "ext.tsv")

    inputs: dict = {
        "h5": None, "per_sample": {}, "struct_tpm": {},
        "global_max_tpm": 0.0, "junctions": [],
    }
    _run_mikado_stage(config, inputs, [])
    assert len(gtf_calls) == 0, "PREP should be skipped when outputs exist"


def test_prep_no_resume_always_runs(monkeypatch, tmp_path):
    """When config.resume=False, PREP always re-runs even if outputs exist."""
    from helixforge.reconcile.pipeline import _run_mikado_stage

    work = tmp_path / "work"
    mik_in = work / "mikado_inputs"
    mik_in.mkdir(parents=True)
    (mik_in / "configuration.yaml").write_text("pick:\n  scoring_file: s.yaml\n")
    (mik_in / "list.txt").write_text("helixer.gtf\n")
    (mik_in / "scoring_strict.yaml").write_text("requirements:\nscoring:\n")
    mik_run = work / "mikado_run"
    mik_run.mkdir(parents=True)

    config = PipelineConfig(
        genome_fasta=str(tmp_path / "g.fa"),
        helixer_gff3=str(tmp_path / "h.gff3"),
        stringtie_list=[str(tmp_path / "s.gtf")],
        protein_db=str(tmp_path / "db.fa"),
        work_dir=str(work),
        resume=False,
    )

    gtf_calls = []
    monkeypatch.setattr(
        pl, "helixer_gff3_to_gtf",
        lambda *a, **k: (gtf_calls.append(1), str(mik_in / "helixer.gtf"))[1],
    )
    monkeypatch.setattr(
        pl, "stringtie_to_labelled_gtf", lambda *a, **k: None,
    )
    monkeypatch.setattr(
        pl, "write_input_list", lambda *a, **k: mik_in / "list.txt",
    )
    monkeypatch.setattr(
        pl, "_emit_junctions_bed", lambda *a, **k: mik_in / "junctions.bed",
    )
    monkeypatch.setattr(
        pl, "install_scoring_profile", lambda *a, **k: mik_in / "scoring_strict.yaml",
    )
    monkeypatch.setattr(
        pl, "write_configuration", lambda *a, **k: str(mik_in / "configuration.yaml"),
    )
    monkeypatch.setattr(
        pl, "run_prepare", lambda *a, **k: (str(mik_run / "p.gtf"), str(mik_run / "p.fa")),
    )
    monkeypatch.setattr(pl, "run_transdecoder", lambda *a, **k: str(mik_run / "orfs.bed"))
    monkeypatch.setattr(pl, "run_diamond", lambda *a, **k: str(mik_run / "blast.out"))
    monkeypatch.setattr(pl, "run_serialise", lambda *a, **k: None)
    monkeypatch.setattr(pl, "run_pick", lambda *a, **k: str(mik_run / "mikado.loci.gff3"))
    monkeypatch.setattr(pl, "parse_loci_gff3", lambda *a, **k: [])
    monkeypatch.setattr(pl, "_build_external_scores", lambda *a, **k: {})
    monkeypatch.setattr(pl, "write_external_scores_tsv", lambda *a, **k: mik_in / "ext.tsv")

    inputs: dict = {
        "h5": None, "per_sample": {}, "struct_tpm": {},
        "global_max_tpm": 0.0, "junctions": [],
    }
    _run_mikado_stage(config, inputs, [])
    assert len(gtf_calls) == 1, "PREP should re-run when resume=False"
