"""Tests for prep/orchestrate.py (Phase 12). Floor: 7. Wrappers mocked."""

from pathlib import Path

import pytest

from helixforge.prep import orchestrate as orch
from helixforge.prep.align import AlignResult
from helixforge.prep.orchestrate import PreppedInputs, prep_evidence


@pytest.fixture
def mocked(monkeypatch):
    """Mock every prep wrapper at the orchestrate module; record the call log."""
    log = {"faidx": 0, "star_index": 0, "hisat2_index": 0,
           "star": [], "hisat2": [], "assemble": [], "miniprot": 0}

    def fake_faidx(fasta, **kw):
        log["faidx"] += 1
        return Path(f"{fasta}.fai")

    def fake_build_star_index(fasta, index_dir, **kw):
        log["star_index"] += 1
        return Path(index_dir)

    def fake_build_hisat2_index(fasta, index_prefix, **kw):
        log["hisat2_index"] += 1
        return index_prefix

    def fake_run_star(genome_dir, reads, out_prefix, **kw):
        log["star"].append((genome_dir, list(reads), out_prefix))
        return AlignResult(
            bam=Path(f"{out_prefix}Aligned.sortedByCoord.out.bam"),
            sj_tab=Path(f"{out_prefix}SJ.out.tab"),
            log=Path(f"{out_prefix}Log.final.out"),
        )

    def fake_run_hisat2(index_prefix, reads, out_bam, **kw):
        log["hisat2"].append((index_prefix, list(reads), Path(out_bam)))
        return AlignResult(bam=Path(out_bam), sj_tab=None, log=None)

    def fake_assemble_samples(bams, out_dir, **kw):
        log["assemble"].append(list(bams))
        return [Path(out_dir) / f"{Path(b).stem}.gtf" for b in bams]

    def fake_run_miniprot(genome, proteome, out_gff, **kw):
        log["miniprot"] += 1
        return Path(out_gff)

    monkeypatch.setattr(orch, "faidx", fake_faidx)
    monkeypatch.setattr(orch, "build_star_index", fake_build_star_index)
    monkeypatch.setattr(orch, "build_hisat2_index", fake_build_hisat2_index)
    monkeypatch.setattr(orch, "run_star", fake_run_star)
    monkeypatch.setattr(orch, "run_hisat2", fake_run_hisat2)
    monkeypatch.setattr(orch, "assemble_samples", fake_assemble_samples)
    monkeypatch.setattr(orch, "run_miniprot", fake_run_miniprot)
    return log


def test_full_star_chain(mocked, tmp_path):
    res = prep_evidence(
        "genome.fa",
        rnaseq_reads=[["a_1.fq", "a_2.fq"], ["b.fq"]],
        proteome="prot.fa",
        out_dir=tmp_path,
    )
    assert isinstance(res, PreppedInputs)
    assert mocked["faidx"] == 1
    assert mocked["star_index"] == 1            # index built once, reused
    assert len(mocked["star"]) == 2             # one align per sample
    assert len(res.bam_paths) == 2
    assert len(res.star_sj_paths) == 2          # STAR emits SJ.out.tab
    assert len(res.stringtie_list) == 2
    assert res.miniprot_gff == str(tmp_path / "miniprot.gff")
    assert res.genome_fasta == "genome.fa"


def test_hisat2_chain_has_no_sj(mocked, tmp_path):
    res = prep_evidence(
        "genome.fa", rnaseq_reads=[["a.fq"]], aligner="hisat2", out_dir=tmp_path,
    )
    assert mocked["hisat2_index"] == 1
    assert len(mocked["hisat2"]) == 1
    assert len(res.bam_paths) == 1
    assert res.star_sj_paths == []              # HISAT2 emits no STAR SJ tab
    assert len(res.stringtie_list) == 1


def test_skip_align_when_no_reads(mocked, tmp_path):
    res = prep_evidence("genome.fa", proteome="prot.fa", out_dir=tmp_path)
    assert mocked["faidx"] == 1
    assert mocked["star_index"] == 0
    assert mocked["star"] == []
    assert mocked["assemble"] == []
    assert res.bam_paths == []
    assert res.stringtie_list == []
    assert res.star_sj_paths == []
    assert res.miniprot_gff == str(tmp_path / "miniprot.gff")


def test_skip_miniprot_when_no_proteome(mocked, tmp_path):
    res = prep_evidence("genome.fa", rnaseq_reads=[["a.fq"]], out_dir=tmp_path)
    assert mocked["miniprot"] == 0
    assert res.miniprot_gff is None
    assert len(res.bam_paths) == 1


def test_faidx_always_runs_even_with_nothing(mocked, tmp_path):
    res = prep_evidence("genome.fa", out_dir=tmp_path)
    assert mocked["faidx"] == 1
    assert res.bam_paths == []
    assert res.stringtie_list == []
    assert res.miniprot_gff is None


def test_single_path_sample_normalized(mocked, tmp_path):
    # a bare string sample (not a list) is treated as single-end
    prep_evidence("genome.fa", rnaseq_reads=["a.fq"], out_dir=tmp_path)
    assert mocked["star"][0][1] == ["a.fq"]


def test_unknown_aligner_raises(mocked, tmp_path):
    with pytest.raises(ValueError, match="unknown aligner"):
        prep_evidence(
            "genome.fa", rnaseq_reads=[["a.fq"]], aligner="bowtie", out_dir=tmp_path,
        )


def test_as_config_kwargs_matches_pipeline_fields(mocked, tmp_path):
    res = prep_evidence(
        "genome.fa", rnaseq_reads=[["a.fq"]], proteome="p.fa", out_dir=tmp_path,
    )
    kwargs = res.as_config_kwargs()
    assert set(kwargs) == {
        "genome_fasta", "bam_paths", "stringtie_list",
        "star_sj_paths", "miniprot_gff",
    }
    # every key is a real PipelineConfig field
    from helixforge.reconcile.pipeline import PipelineConfig
    field_names = set(PipelineConfig.__dataclass_fields__)
    assert set(kwargs).issubset(field_names)
