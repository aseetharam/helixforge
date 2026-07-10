"""Phase 24 D1 — the ``helixforge run`` whole-genome driver. Floor: 5.

``scatter="off"`` must equal the single-process annotation byte-for-byte;
``scatter=N`` must plan → execute → aggregate into globally-unique, stable HFG
IDs with no Helixer locus lost. Chunk execution is run **in-process** via a
ThreadPoolExecutor seam (no real Mikado — the synthetic genome runs in basic
mode). Concrete literal coordinates from the shared ``synthetic_genome`` fixture.
"""

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from helixforge.io.fasta import GenomeAccessor
from helixforge.parallel import run as runmod
from helixforge.parallel.run import RunResult, run_genome
from helixforge.parallel.tasks import run_local as _real_run_local
from helixforge.reconcile.pipeline import PipelineConfig, run_pipeline


@pytest.fixture
def base_config(synthetic_genome, tmp_path):
    # Pre-build the .fai so concurrent chunk threads only read it (no rebuild race).
    GenomeAccessor(synthetic_genome["fasta"]).close()
    return PipelineConfig(
        genome_fasta=synthetic_genome["fasta"],
        helixer_gff3=synthetic_genome["gff3"],
        output_prefix=str(tmp_path / "genome"),
    )


@pytest.fixture
def in_process_pool(monkeypatch):
    """Run chunks in-process (ThreadPoolExecutor) instead of spawning processes."""
    def patched(plan, base_config, workers=4):
        return _real_run_local(plan, base_config, workers=workers,
                               _executor_cls=ThreadPoolExecutor)

    monkeypatch.setattr(runmod, "run_local", patched)


# ---------------------------------------------------------------------------
# scatter="off" == single-process
# ---------------------------------------------------------------------------

def test_scatter_off_equals_single_process(base_config, tmp_path):
    result = run_genome(base_config, scatter="off")
    assert isinstance(result, RunResult)
    assert result.mode == "single"
    driver_ids = [g.gene_id for g in result.genes]

    # A direct run_pipeline with a separate prefix must yield identical genes.
    direct_cfg = PipelineConfig(
        genome_fasta=base_config.genome_fasta,
        helixer_gff3=base_config.helixer_gff3,
        output_prefix=str(tmp_path / "direct"),
    )
    direct_ids = [g.gene_id for g in run_pipeline(direct_cfg)]
    assert driver_ids == direct_ids
    assert len(driver_ids) == 6           # 5 chr1 + 1 chr2 loci, none lost


def test_scatter_off_writes_manifest(base_config):
    result = run_genome(base_config, scatter="off")
    payload = json.loads(Path(result.manifest_path).read_text())
    assert payload["mode"] == "single"
    assert payload["num_genes"] == 6


# ---------------------------------------------------------------------------
# scatter=N local: plan → exec → aggregate
# ---------------------------------------------------------------------------

def test_scatter_local_unique_ids_no_gene_lost(base_config, in_process_pool, tmp_path):
    master = tmp_path / "master.id_map.json"
    result = run_genome(
        base_config, scatter=3, hpc="local", workers=2,
        master_id_map_path=str(master),
    )
    assert result.mode == "scatter-local"
    agg = result.aggregate
    # Every Helixer locus covered exactly once across all chunks.
    assert agg.num_loci == 6
    assert agg.num_genes == 6
    # Globally-unique HFG IDs in the merged map.
    hfgs = list(json.loads(master.read_text()).values())
    assert len(hfgs) == len(set(hfgs)) == 6


def test_scatter_local_manifest_counts(base_config, in_process_pool):
    result = run_genome(base_config, scatter=3, hpc="local", workers=2)
    payload = json.loads(Path(result.manifest_path).read_text())
    assert payload["mode"] == "scatter-local"
    assert payload["num_genes"] == 6
    assert payload["num_chunks"] == len(result.plan)
    assert len(payload["chunk_prefixes"]) == len(result.plan)


# ---------------------------------------------------------------------------
# scatter=N hypershell: plan → `hs cluster` → aggregate (hs run mocked)
# ---------------------------------------------------------------------------

@pytest.fixture
def in_process_hs(monkeypatch):
    """Stand in for `hs cluster` by running the chunks in-process, so the driver's
    hypershell branch (plan → execute → aggregate) is covered without a real `hs`."""
    import helixforge.parallel.hypershell_backend as hsmod

    def patched(plan, base_config, out_prefix, *, num_tasks, hs_bin="hs", **kw):
        return _real_run_local(plan, base_config, workers=num_tasks,
                               _executor_cls=ThreadPoolExecutor)

    monkeypatch.setattr(hsmod, "run_hypershell", patched)


def test_scatter_hypershell_unique_ids_no_gene_lost(base_config, in_process_hs, tmp_path):
    master = tmp_path / "master.id_map.json"
    result = run_genome(
        base_config, scatter=3, hpc="hypershell", num_tasks=2,
        master_id_map_path=str(master),
    )
    assert result.mode == "scatter-hypershell"
    agg = result.aggregate
    assert agg.num_loci == 6
    assert agg.num_genes == 6
    hfgs = list(json.loads(master.read_text()).values())
    assert len(hfgs) == len(set(hfgs)) == 6


def test_scatter_hypershell_manifest(base_config, in_process_hs):
    result = run_genome(base_config, scatter=3, hpc="hypershell", num_tasks=2)
    payload = json.loads(Path(result.manifest_path).read_text())
    assert payload["mode"] == "scatter-hypershell"
    assert payload["hpc"] == "hypershell"
    assert payload["num_genes"] == 6
    assert payload["num_chunks"] == len(result.plan)


def test_scatter_local_stable_ids_across_reruns(base_config, in_process_pool, tmp_path):
    master = str(tmp_path / "master.id_map.json")
    first = run_genome(base_config, scatter=3, hpc="local", workers=2,
                       master_id_map_path=master)
    first_map = json.loads(Path(master).read_text())
    # A second run seeds chunk ranges from the persisted master → identical HFGs.
    second = run_genome(base_config, scatter=3, hpc="local", workers=2,
                        master_id_map_path=master)
    second_map = json.loads(Path(master).read_text())
    assert first_map == second_map
    assert first.aggregate.num_genes == second.aggregate.num_genes == 6


# ---------------------------------------------------------------------------
# scatter slurm: emit array, do not run
# ---------------------------------------------------------------------------

def test_scatter_slurm_emits_array_without_running(base_config, tmp_path):
    result = run_genome(base_config, scatter=2, hpc="slurm")
    assert result.mode == "scatter-slurm"
    assert result.script_path.exists()
    text = result.script_path.read_text()
    assert "#SBATCH --array=" in text
    assert "helixforge reconcile" in text
    # No aggregated genome GFF3 was produced (chunks run under the scheduler).
    assert not Path(f"{base_config.output_prefix}.gff3").exists()
    payload = json.loads(Path(result.manifest_path).read_text())
    assert payload["mode"] == "scatter-slurm"
    assert "slurm_array" in payload


# ---------------------------------------------------------------------------
# CLI: `reconcile --scatter` (the folded-in `run`) == single-pass `reconcile`
# ---------------------------------------------------------------------------

def _gene_ids_from_gff3(path):
    ids = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9 or cols[2] != "gene":
            continue
        for part in cols[8].split(";"):
            if part.startswith("ID="):
                ids.append(part[3:])
                break
    return sorted(ids)


def test_cli_reconcile_chunked_equals_single_pass(
    synthetic_genome, in_process_pool, tmp_path
):
    """`reconcile --scatter N` and single-pass `reconcile` agree gene-for-gene.

    Chunking changes memory/speed, not the annotation: the merged chunked GFF3
    must carry exactly the same gene set as the single in-process pass.
    """
    from click.testing import CliRunner

    from helixforge import cli

    GenomeAccessor(synthetic_genome["fasta"]).close()  # prebuild .fai (no race)
    runner = CliRunner()

    single_prefix = str(tmp_path / "single")
    r1 = runner.invoke(cli.main, [
        "reconcile",
        "--genome", synthetic_genome["fasta"],
        "--helixer", synthetic_genome["gff3"],
        "--output-prefix", single_prefix,
    ])
    assert r1.exit_code == 0, r1.output

    chunk_prefix = str(tmp_path / "chunked")
    r2 = runner.invoke(cli.main, [
        "reconcile",
        "--genome", synthetic_genome["fasta"],
        "--helixer", synthetic_genome["gff3"],
        "--output-prefix", chunk_prefix,
        "--scatter", "3", "--workers", "2", "--skip-preflight",
    ])
    assert r2.exit_code == 0, r2.output

    single_ids = _gene_ids_from_gff3(f"{single_prefix}.gff3")
    chunked_ids = _gene_ids_from_gff3(f"{chunk_prefix}.gff3")
    assert single_ids == chunked_ids
    assert len(single_ids) == 6  # 5 chr1 + 1 chr2 loci, none lost
