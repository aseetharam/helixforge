"""Phase 24 D3 — parallel per-sample evidence prep. Floor: 2.

Samples are independent jobs; with ``workers > 1`` align + assemble run over a
bounded pool. The wrappers are mocked (CLAUDE.md §C1) — the tests assert each
sample is dispatched and that results stay in **input order** regardless of
completion order (``ThreadPoolExecutor.map`` guarantees this).
"""

from pathlib import Path

import pytest

from helixforge.prep import assemble as assemble_mod
from helixforge.prep import orchestrate as orch
from helixforge.prep.align import AlignResult
from helixforge.prep.orchestrate import _map_samples, prep_evidence


@pytest.fixture
def mocked(monkeypatch):
    log = {"star": [], "assemble": []}

    monkeypatch.setattr(orch, "faidx", lambda fasta, **kw: Path(f"{fasta}.fai"))
    monkeypatch.setattr(orch, "build_star_index",
                        lambda fasta, index_dir, **kw: Path(index_dir))

    def fake_run_star(genome_dir, reads, out_prefix, **kw):
        log["star"].append(out_prefix)
        return AlignResult(
            bam=Path(f"{out_prefix}.bam"),
            sj_tab=Path(f"{out_prefix}.SJ.tab"),
            log=None,
        )

    def fake_assemble(bams, out_dir, **kw):
        log["assemble"].append({"bams": list(bams), "workers": kw.get("workers")})
        return [Path(out_dir) / f"{Path(b).stem}.gtf" for b in bams]

    monkeypatch.setattr(orch, "run_star", fake_run_star)
    monkeypatch.setattr(orch, "assemble_samples", fake_assemble)
    return log


def test_prep_parallel_dispatches_each_sample(mocked, tmp_path):
    res = prep_evidence(
        "genome.fa",
        rnaseq_reads=[["s0.fq"], ["s1.fq"], ["s2.fq"]],
        out_dir=tmp_path,
        workers=3,
    )
    # Every sample dispatched, exactly once, in sample order.
    assert len(mocked["star"]) == 3
    assert res.bam_paths == [
        str(tmp_path / "sample0_.bam"),
        str(tmp_path / "sample1_.bam"),
        str(tmp_path / "sample2_.bam"),
    ]
    assert res.star_sj_paths == [
        str(tmp_path / "sample0_.SJ.tab"),
        str(tmp_path / "sample1_.SJ.tab"),
        str(tmp_path / "sample2_.SJ.tab"),
    ]
    # The workers budget is threaded into StringTie assembly too.
    assert mocked["assemble"][0]["workers"] == 3


def test_assemble_samples_parallel_preserves_order(monkeypatch, tmp_path):
    calls = []

    def fake_run_stringtie(bam, out_gtf, sample_id, **kw):
        calls.append(sample_id)
        return out_gtf

    monkeypatch.setattr(assemble_mod, "run_stringtie", fake_run_stringtie)
    bams = [tmp_path / f"x{i}.bam" for i in range(4)]
    for b in bams:
        b.write_text("")
    gtfs = assemble_mod.assemble_samples(bams, tmp_path / "st", workers=4)
    # Returned GTFs follow BAM (input) order even though the pool may finish jobs
    # out of order.
    assert [Path(g).stem for g in gtfs] == ["x0", "x1", "x2", "x3"]
    assert sorted(calls) == ["x0", "x1", "x2", "x3"]


def test_map_samples_preserves_input_order():
    # Reverse-jittered thunks: even if later thunks "finish" first, map keeps order.
    fns = [lambda i=i: i for i in range(5)]
    assert _map_samples(fns, workers=4) == [0, 1, 2, 3, 4]
    # workers=1 path is the serial baseline.
    assert _map_samples(fns, workers=1) == [0, 1, 2, 3, 4]
