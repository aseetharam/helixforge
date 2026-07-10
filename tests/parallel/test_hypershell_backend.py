"""Tests for parallel/hypershell_backend.py — subprocess-first HyperShell path.

Floor: 8. `hs` is never really invoked — a fake runner records argv. Both
strands via the ``synthetic_genome`` fixture (chr1 plus + chr2 minus, whole
scaffold). Concrete literal coordinates. No ``import hypershell`` occurs: the
backend drives the CLI by subprocess (optional ``helixforge[hpc]`` extra).
"""

import json
import sys
from types import SimpleNamespace

import pytest

from helixforge.parallel.hypershell_backend import (
    run_hypershell,
    write_hypershell_plan,
)
from helixforge.parallel.plan import partition_genome, reserve_id_ranges
from helixforge.reconcile.pipeline import PipelineConfig


@pytest.fixture
def plan(synthetic_genome):
    p = partition_genome(
        synthetic_genome["fasta"], synthetic_genome["gff3"],
        target_loci_per_chunk=2, min_boundary_gap=1000,
    )
    reserve_id_ranges(p)
    return p


@pytest.fixture
def base_config(synthetic_genome, tmp_path):
    return PipelineConfig(
        genome_fasta=synthetic_genome["fasta"],
        helixer_gff3=synthetic_genome["gff3"],
        output_prefix=str(tmp_path / "run"),
        procs=2,
    )


class _FakeRunner:
    """Records argv; answers `hs --version` and `hs cluster` without a real tool."""

    def __init__(self, returncode: int = 0, stderr: str = ""):
        self.calls: list[list[str]] = []
        self.returncode = returncode
        self.stderr = stderr

    def __call__(self, argv, **kwargs):
        self.calls.append(list(argv))
        if len(argv) >= 2 and argv[1] == "--version":
            return SimpleNamespace(returncode=0, stdout="HyperShell 2.8.1", stderr="")
        return SimpleNamespace(returncode=self.returncode, stdout="", stderr=self.stderr)

    def cluster_call(self) -> list[str]:
        return [c for c in self.calls if len(c) >= 2 and c[1] == "cluster"][0]


# --------------------------------------------------------------------------
# run_hypershell — subprocess dispatch
# --------------------------------------------------------------------------


def test_run_hypershell_invokes_hs_cluster(plan, base_config, tmp_path):
    runner = _FakeRunner()
    prefixes = run_hypershell(
        plan, base_config, str(tmp_path / "run"), num_tasks=8, _runner=runner,
    )
    call = runner.cluster_call()
    assert call[0] == "hs"
    assert call[1] == "cluster"
    assert str(tmp_path / "run.tasks.txt") in call
    assert call[call.index("--num-tasks") + 1] == "8"
    assert len(prefixes) == len(plan.chunks)
    assert prefixes[0].name == plan.chunks[0].chunk_id


def test_run_hypershell_tasks_file_one_line_per_chunk(plan, base_config, tmp_path):
    runner = _FakeRunner()
    run_hypershell(plan, base_config, str(tmp_path / "run"), num_tasks=4, _runner=runner)
    lines = (tmp_path / "run.tasks.txt").read_text().splitlines()
    assert len(lines) == len(plan.chunks)
    assert all("reconcile" in ln for ln in lines)
    assert "--region chr2" in lines[-1]  # minus-strand whole scaffold


def test_run_hypershell_custom_bin_and_num_tasks(plan, base_config, tmp_path):
    runner = _FakeRunner()
    run_hypershell(
        plan, base_config, str(tmp_path / "r"),
        num_tasks=16, hs_bin="/opt/hs", _runner=runner,
    )
    call = runner.cluster_call()
    assert call[0] == "/opt/hs"
    assert call[call.index("--num-tasks") + 1] == "16"


def test_run_hypershell_nonzero_exit_raises(plan, base_config, tmp_path):
    runner = _FakeRunner(returncode=2, stderr="boom\nsomething failed")
    with pytest.raises(RuntimeError, match="exited 2"):
        run_hypershell(plan, base_config, str(tmp_path / "r"), num_tasks=4, _runner=runner)


def test_run_hypershell_missing_binary_raises(plan, base_config, tmp_path):
    def _missing(argv, **kwargs):
        raise FileNotFoundError(argv[0])

    with pytest.raises(RuntimeError, match="not found"):
        run_hypershell(plan, base_config, str(tmp_path / "r"), num_tasks=4, _runner=_missing)


# --------------------------------------------------------------------------
# write_hypershell_plan — #37 native-plan emitter
# --------------------------------------------------------------------------


def test_write_hypershell_plan_schema(plan, base_config, tmp_path):
    path, template = write_hypershell_plan(plan, base_config, tmp_path / "plan.hs.json")
    payload = json.loads(path.read_text())
    assert set(payload) == {"metadata", "chunks"}
    assert len(payload["chunks"]) == len(plan.chunks)
    c0 = payload["chunks"][0]
    assert set(c0) == {"chunk_id", "region", "id_start", "novel_start", "output_dir"}
    assert c0["chunk_id"] == plan.chunks[0].chunk_id
    assert c0["region"] == "chr1:0-3500"
    assert c0["id_start"] == "1"
    assert c0["novel_start"] == "90000"
    # minus-strand whole scaffold pre-formatted as a bare seqid region
    assert payload["chunks"][-1]["region"] == "chr2"
    # template only references pre-resolved chunk fields
    assert "{region}" in template
    assert "{id_start}" in template
    assert "{novel_start}" in template


def test_write_hypershell_plan_rejects_multiregion(base_config, tmp_path):
    from helixforge.parallel.plan import Chunk, Plan

    p = Plan(
        chunks=[Chunk(chunk_id="chunk_0000",
                      regions=["chr1:0-100", "chr2:0-100"], num_loci=2)],
        total_loci=2,
    )
    reserve_id_ranges(p)
    with pytest.raises(NotImplementedError, match="single-region"):
        write_hypershell_plan(p, base_config, tmp_path / "x.json")


# --------------------------------------------------------------------------
# no hard dependency on hypershell
# --------------------------------------------------------------------------


def test_backend_does_not_import_hypershell():
    import helixforge.parallel.hypershell_backend  # noqa: F401

    assert "hypershell" not in sys.modules
