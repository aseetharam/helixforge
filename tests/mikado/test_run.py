"""Tests for mikado/run.py (Phase 5). Floor: 14. Subprocess fully mocked."""

import subprocess

import pytest

from helixforge.mikado import run as run_mod
from helixforge.mikado.run import (
    MikadoInputs,
    MikadoRunResult,
    run_diamond,
    run_mikado,
    run_pick,
    run_prepare,
    run_serialise,
    run_transdecoder,
)


@pytest.fixture
def capture_argv(monkeypatch):
    """Record every subprocess.run argv; return a successful CompletedProcess."""
    class _Calls(list):
        cwds = None  # parallel list of cwd per call (None unless set)

    calls = _Calls()
    calls.cwds = []

    def fake_run(argv, check, capture_output, text, cwd=None):
        calls.append(argv)
        calls.cwds.append(cwd)
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(run_mod.subprocess, "run", fake_run)
    # Phase 13: run_mikado now resolves each tool's bin (shutil.which) and probes
    # the mikado version before running a step. Stub both so the chain tests stay
    # hermetic (no real binaries) and fast.
    monkeypatch.setattr(run_mod.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(run_mod, "tool_version", lambda b, a="--version": "mikado 2.3.4")
    return calls


# --- individual steps ---

def test_prepare_argv_and_paths(capture_argv, tmp_path):
    gtf, fasta = run_prepare("config.yaml", tmp_path, procs=4)
    argv = capture_argv[0]
    assert argv[0] == "mikado" and argv[1] == "prepare"
    assert "--configuration" in argv
    assert argv[argv.index("--procs") + 1] == "4"
    assert gtf.name == "mikado_prepared.gtf"
    assert fasta.name == "mikado_prepared.fasta"


def test_prepare_custom_bin(capture_argv, tmp_path):
    run_prepare("config.yaml", tmp_path, mikado_bin="/opt/mikado")
    assert capture_argv[0][0] == "/opt/mikado"


def test_transdecoder_two_steps_and_path(capture_argv, tmp_path):
    bed = run_transdecoder(tmp_path / "mikado_prepared.fasta", tmp_path)
    assert "TransDecoder.LongOrfs" in capture_argv[0][0]
    assert "TransDecoder.Predict" in capture_argv[1][0]
    assert "--single_best_only" in capture_argv[1]
    assert bed.name == "mikado_prepared.fasta.transdecoder.bed"
    # final *.transdecoder.* outputs land in out_dir -> both steps run there
    assert capture_argv.cwds[0] == str(tmp_path)
    assert capture_argv.cwds[1] == str(tmp_path)


def test_transdecoder_single_best_off(capture_argv, tmp_path):
    run_transdecoder(tmp_path / "p.fasta", tmp_path, single_best_only=False)
    assert "--single_best_only" not in capture_argv[1]


def test_transdecoder_bin_dir(capture_argv, tmp_path):
    run_transdecoder(tmp_path / "p.fasta", tmp_path, transdecoder_bin_dir="/td/bin")
    assert capture_argv[0][0] == "/td/bin/TransDecoder.LongOrfs"


def test_diamond_makedb_then_blastx_xml(capture_argv, tmp_path):
    out = run_diamond(tmp_path / "p.fasta", "prot.fa", tmp_path)
    assert capture_argv[0][1] == "makedb"
    assert capture_argv[1][1] == "blastx"
    assert "5" in capture_argv[1]  # --outfmt 5 (XML)
    assert out.name == "mikado_diamond.xml"


def test_diamond_tabular_format(capture_argv, tmp_path):
    out = run_diamond(tmp_path / "p.fasta", "prot.fa", tmp_path, out_format="tabular")
    blastx = capture_argv[1]
    assert "6" in blastx
    assert "qseqid" in blastx
    assert out.name == "mikado_diamond.tsv"


def test_diamond_params_passed(capture_argv, tmp_path):
    run_diamond(tmp_path / "p.fasta", "prot.fa", tmp_path, threads=8, max_target_seqs=10)
    blastx = capture_argv[1]
    assert blastx[blastx.index("--threads") + 1] == "8"
    assert blastx[blastx.index("--max-target-seqs") + 1] == "10"


def test_serialise_argv_and_path(capture_argv, tmp_path):
    db = run_serialise(
        "config.yaml", "p.fasta", "orfs.bed", "blast.xml", "prot.fa",
        "j.bed", "ext.tsv", "genome.fa", tmp_path,
    )
    argv = capture_argv[0]
    assert argv[1] == "serialise"
    assert "--external-scores" in argv
    assert "--junctions" in argv
    # Mikado >=2.3 serialise uses --genome (not the old --genome_fasta)
    assert "--genome" in argv and "--genome_fasta" not in argv
    assert db.name == "mikado.db"


def test_pick_argv_and_path(capture_argv, tmp_path):
    loci = run_pick("config.yaml", "scoring.yaml", "prepared.gtf", tmp_path, procs=2)
    argv = capture_argv[0]
    assert argv[1] == "pick"
    assert "--scoring-file" in argv
    assert argv[argv.index("--procs") + 1] == "2"
    assert loci.name == "mikado.loci.gff3"


# --- full chain ---

def _inputs():
    return MikadoInputs(
        configuration_yaml="config.yaml",
        scoring_file="scoring.yaml",
        genome_fa="genome.fa",
        protein_db="prot.fa",
        junctions_tab="j.bed",
        external_tsv="ext.tsv",
    )


def test_run_mikado_chains_all_steps(capture_argv, tmp_path):
    run_mikado(_inputs(), tmp_path)
    subcommands = [c[1] if c[0] in ("mikado", "diamond") else c[0] for c in capture_argv]
    # prepare, LongOrfs, Predict, makedb, blastx, serialise, pick
    assert subcommands[0] == "prepare"
    assert "LongOrfs" in subcommands[1]
    assert "Predict" in subcommands[2]
    assert subcommands[3] == "makedb"
    assert subcommands[4] == "blastx"
    assert subcommands[5] == "serialise"
    assert subcommands[6] == "pick"


def test_run_mikado_result_paths(capture_argv, tmp_path):
    result = run_mikado(_inputs(), tmp_path)
    assert isinstance(result, MikadoRunResult)
    assert result.loci_gff3.name == "mikado.loci.gff3"
    assert result.prepared_fasta.name == "mikado_prepared.fasta"
    assert result.mikado_db.name == "mikado.db"
    assert result.loci_metrics_tsv.name == "mikado.loci.metrics.tsv"
    assert result.loci_scores_tsv.name == "mikado.loci.scores.tsv"


# --- error handling ---

def test_nonzero_exit_raises_informative(monkeypatch, tmp_path):
    def fake_run(argv, check, capture_output, text, cwd=None):
        raise subprocess.CalledProcessError(
            returncode=1, cmd=argv, stderr="boom: traceback details here"
        )

    monkeypatch.setattr(run_mod.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError) as exc:
        run_prepare("config.yaml", tmp_path)
    msg = str(exc.value)
    assert "prepare" in msg
    assert "boom: traceback details here" in msg
    assert "exit 1" in msg


def test_argv_is_list_not_string(monkeypatch, tmp_path):
    seen = {}

    def fake_run(argv, check, capture_output, text, cwd=None):
        seen["argv"] = argv
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(run_mod.subprocess, "run", fake_run)
    run_prepare("config.yaml", tmp_path)
    assert isinstance(seen["argv"], list)
    assert all(isinstance(a, str) for a in seen["argv"])
