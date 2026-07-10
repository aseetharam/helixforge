"""Tests for reconcile/pipeline.py (Phase 8). Floor: 20.

All heavy I/O and every Mikado/TransDecoder/DIAMOND/Portcullis subprocess is
mocked; the pipeline control flow is exercised for real. Synthetic fixtures
only (CLAUDE.md §12) — both strands represented. Concrete literal coordinates.
"""

from pathlib import Path

import pytest

from helixforge.reconcile import pipeline as pl
from helixforge.reconcile.models import Exon, SpliceJunction
from helixforge.reconcile.pipeline import PipelineConfig, run_pipeline


# ---------------------------------------------------------------------------
# Synthetic inputs
#   g1  chr1 +  101-200   single exon   (overlaps the mocked Mikado locus)
#   g2  chr1 -  301-600   two exons (intron 400-500)   multi-exon minus strand
#   g3  chr2 +  101-250   single exon   second scaffold
# ---------------------------------------------------------------------------

_HELIXER_GFF = """\
##gff-version 3
chr1\tHelixer\tgene\t101\t200\t.\t+\t.\tID=g1
chr1\tHelixer\tmRNA\t101\t200\t.\t+\t.\tID=g1.m;Parent=g1
chr1\tHelixer\texon\t101\t200\t.\t+\t.\tID=g1.e1;Parent=g1.m
chr1\tHelixer\tgene\t301\t600\t.\t-\t.\tID=g2
chr1\tHelixer\tmRNA\t301\t600\t.\t-\t.\tID=g2.m;Parent=g2
chr1\tHelixer\texon\t301\t400\t.\t-\t.\tID=g2.e1;Parent=g2.m
chr1\tHelixer\texon\t501\t600\t.\t-\t.\tID=g2.e2;Parent=g2.m
chr2\tHelixer\tgene\t101\t250\t.\t+\t.\tID=g3
chr2\tHelixer\tmRNA\t101\t250\t.\t+\t.\tID=g3.m;Parent=g3
chr2\tHelixer\texon\t101\t250\t.\t+\t.\tID=g3.e1;Parent=g3.m
"""

_STRINGTIE_GTF = """\
chr1\tStringTie\ttranscript\t101\t200\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; TPM "5.0";
chr1\tStringTie\texon\t101\t200\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; TPM "5.0";
"""

# Mocked `mikado pick` output — one coding locus 1:1 with Helixer g1.
_LOCI_GFF3 = """\
##gff-version 3
chr1\tMikado\tgene\t101\t200\t.\t+\t.\tID=mikado.1G0;Name=mikado.1G0
chr1\tMikado\tmRNA\t101\t200\t.\t+\t.\tID=mikado.1G0.1;Parent=mikado.1G0
chr1\tMikado\texon\t101\t200\t.\t+\t.\tID=mikado.1G0.1.e1;Parent=mikado.1G0.1
chr1\tMikado\tCDS\t101\t199\t.\t+\t0\tID=mikado.1G0.1.c1;Parent=mikado.1G0.1
"""
_METRICS = "tid\tblast_score\thas_start_codon\thas_stop_codon\nmikado.1G0.1\t150.0\tTrue\tTrue\n"
_SCORES = "tid\tscore\nmikado.1G0.1\t15.0\n"


@pytest.fixture
def genome_fasta(tmp_path):
    fa = tmp_path / "genome.fasta"
    fa.write_text(">chr1\n" + "ACGT" * 200 + "\n>chr2\n" + "ACGT" * 100 + "\n")
    return str(fa)


@pytest.fixture
def helixer_gff3(tmp_path):
    p = tmp_path / "helixer.gff3"
    p.write_text(_HELIXER_GFF)
    return str(p)


@pytest.fixture
def stringtie_gtf(tmp_path):
    p = tmp_path / "sampleA.gtf"
    p.write_text(_STRINGTIE_GTF)
    return str(p)


@pytest.fixture
def basic_config(tmp_path, genome_fasta, helixer_gff3):
    return PipelineConfig(
        genome_fasta=genome_fasta,
        helixer_gff3=helixer_gff3,
        output_prefix=str(tmp_path / "out"),
        work_dir=str(tmp_path / "work"),
    )


@pytest.fixture
def mock_mikado(monkeypatch):
    """Mock every Mikado subprocess step; `run_pick` drops a synthetic loci.gff3."""
    def fake_prepare(cfg, out_dir, procs=1, mikado_bin="mikado"):
        out_dir = Path(out_dir)
        gtf = out_dir / "mikado_prepared.gtf"
        gtf.write_text(
            'chr1\tMikado\texon\t101\t200\t.\t+\t.\ttranscript_id "mikado.1G0.1";\n'
        )
        fasta = out_dir / "mikado_prepared.fasta"
        fasta.write_text(">mikado.1G0.1\nACGTACGT\n")
        return gtf, fasta

    def fake_transdecoder(prepared_fasta, out_dir, **kw):
        bed = Path(out_dir) / "orfs.bed"
        bed.write_text("")
        return bed

    def fake_diamond(prepared_fasta, db, out_dir, **kw):
        xml = Path(out_dir) / "diamond.xml"
        xml.write_text("")
        return xml

    def fake_serialise(*a, **k):
        return Path("mikado.db")

    def fake_pick(cfg, scoring, prepared_gtf, out_dir, procs=1, mikado_bin="mikado"):
        out_dir = Path(out_dir)
        (out_dir / "mikado.loci.gff3").write_text(_LOCI_GFF3)
        (out_dir / "mikado.loci.metrics.tsv").write_text(_METRICS)
        (out_dir / "mikado.loci.scores.tsv").write_text(_SCORES)
        return out_dir / "mikado.loci.gff3"

    monkeypatch.setattr(pl, "run_prepare", fake_prepare)
    monkeypatch.setattr(pl, "run_transdecoder", fake_transdecoder)
    monkeypatch.setattr(pl, "run_diamond", fake_diamond)
    monkeypatch.setattr(pl, "run_serialise", fake_serialise)
    monkeypatch.setattr(pl, "run_pick", fake_pick)
    # The startup preflight checks shutil.which for tools before the mocked
    # run_* functions are reached; stub it so the test chain proceeds.
    from helixforge.prep import preflight as _pf
    monkeypatch.setattr(_pf, "check_tool_chain", lambda cfg: [])


@pytest.fixture
def full_config(tmp_path, genome_fasta, helixer_gff3, stringtie_gtf):
    return PipelineConfig(
        genome_fasta=genome_fasta,
        helixer_gff3=helixer_gff3,
        stringtie_list=[stringtie_gtf],
        protein_db=str(tmp_path / "proteins.fasta"),  # path only; diamond mocked
        output_prefix=str(tmp_path / "out"),
        work_dir=str(tmp_path / "work"),
        use_mikado_configure=False,  # templated YAML — no mikado binary needed
    )


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def test_config_required_fields(genome_fasta, helixer_gff3):
    cfg = PipelineConfig(genome_fasta=genome_fasta, helixer_gff3=helixer_gff3)
    assert cfg.genome_fasta == genome_fasta
    assert cfg.helixer_gff3 == helixer_gff3


def test_config_derives_output_paths():
    cfg = PipelineConfig(genome_fasta="g.fa", helixer_gff3="h.gff3",
                         output_prefix="/tmp/run")
    assert cfg.report_path == "/tmp/run.report.tsv"
    assert cfg.id_map_path == "/tmp/run.id_map.json"
    assert cfg.work_dir == "/tmp/run_work"


def test_config_explicit_paths_respected():
    cfg = PipelineConfig(genome_fasta="g.fa", helixer_gff3="h.gff3",
                         report_path="/r.tsv", id_map_path="/m.json",
                         work_dir="/w")
    assert cfg.report_path == "/r.tsv"
    assert cfg.id_map_path == "/m.json"
    assert cfg.work_dir == "/w"


def test_config_chunk_id_namespaces_work_dir():
    cfg = PipelineConfig(genome_fasta="g.fa", helixer_gff3="h.gff3",
                         output_prefix="run", chunk_id="chunk7")
    assert cfg.work_dir == "run_work/chunk7"


def test_config_default_knobs():
    cfg = PipelineConfig(genome_fasta="g.fa", helixer_gff3="h.gff3")
    assert cfg.scoring_profile == "strict"
    assert cfg.pad is True
    assert cfg.max_isoforms == 5
    assert cfg.reciprocal_overlap == 0.5
    assert cfg.short_cds_threshold == 300


def test_should_run_mikado_gate():
    base = dict(genome_fasta="g.fa", helixer_gff3="h.gff3")
    assert pl._should_run_mikado(PipelineConfig(**base)) is False
    assert pl._should_run_mikado(
        PipelineConfig(stringtie_list=["a.gtf"], **base)) is False
    assert pl._should_run_mikado(
        PipelineConfig(stringtie_list=["a.gtf"], protein_db="p.fa", **base)) is True


def test_should_run_mikado_gate_still_false_without_protein_db():
    """Gate logic unchanged: StringTie alone does not enable Mikado."""
    base = dict(genome_fasta="g.fa", helixer_gff3="h.gff3")
    cfg = PipelineConfig(stringtie_list=["a.gtf"], **base)
    assert pl._should_run_mikado(cfg) is False


def test_mikado_skip_warns_when_stringtie_without_protein_db(caplog):
    """When StringTie is present but --protein-db is missing, the log message
    names the missing input and says 'no isoform discovery'."""
    import logging

    base = dict(genome_fasta="g.fa", helixer_gff3="h.gff3")
    cfg = PipelineConfig(stringtie_list=["a.gtf"], **base)
    with caplog.at_level(logging.WARNING, logger="helixforge.reconcile.pipeline"):
        result = pl._run_mikado_stage(cfg, {}, [])
    assert result == []
    assert any("--protein-db" in m for m in caplog.messages)
    assert any("no isoform discovery" in m for m in caplog.messages)


def test_mikado_skip_info_when_no_stringtie_no_protein(caplog):
    """When neither StringTie nor protein_db is present, the log is INFO
    and does NOT claim StringTie is absent when it isn't."""
    import logging

    base = dict(genome_fasta="g.fa", helixer_gff3="h.gff3")
    cfg = PipelineConfig(**base)
    with caplog.at_level(logging.INFO, logger="helixforge.reconcile.pipeline"):
        result = pl._run_mikado_stage(cfg, {}, [])
    assert result == []
    assert any("basic mode" in m for m in caplog.messages)
    assert not any("--protein-db" in m for m in caplog.messages)


def test_run_pipeline_warns_stringtie_without_protein_db(tmp_path):
    """run_pipeline emits a UserWarning at startup when StringTie is present
    but --protein-db is missing, before the heavy stages."""
    import warnings

    cfg = PipelineConfig(
        genome_fasta=str(tmp_path / "nonexistent.fa"),
        helixer_gff3=str(tmp_path / "nonexistent.gff3"),
        stringtie_list=["a.gtf"],
        work_dir=str(tmp_path / "work"),
    )
    with pytest.warns(UserWarning, match="--protein-db not supplied"):
        with pytest.raises(FileNotFoundError):
            run_pipeline(cfg)


def test_run_pipeline_fails_fast_on_bad_scoring_profile(tmp_path):
    """Scoring profile is validated before any I/O, so even a missing genome
    raises FileNotFoundError for the profile, not for the genome file."""
    cfg = PipelineConfig(
        genome_fasta=str(tmp_path / "nonexistent.fa"),
        helixer_gff3=str(tmp_path / "nonexistent.gff3"),
        scoring_profile="totally_bogus_profile",
        work_dir=str(tmp_path / "work"),
    )
    with pytest.raises(FileNotFoundError, match="unknown scoring profile"):
        run_pipeline(cfg)


# ---------------------------------------------------------------------------
# Helper units
# ---------------------------------------------------------------------------

def test_structure_key_format():
    key = pl._structure_key("chr1", "+", [Exon(100, 200), Exon(300, 400)])
    assert key == "chr1:+:((100, 200), (300, 400))"


def test_dedup_flags():
    # Phase 19 D1: the four duplicated ``_dedup_flags`` copies were unified into
    # ``helixforge.qc.flags.dedup_flags`` (re-imported by pipeline.py).
    from helixforge.qc.flags import HELIXER_ONLY, NO_EXPRESSION, dedup_flags
    out = dedup_flags([HELIXER_ONLY, NO_EXPRESSION, HELIXER_ONLY])
    assert [f.name for f in out] == ["HELIXER_ONLY", "NO_EXPRESSION"]


def test_parse_prepared_gtf_both_strands(tmp_path):
    gtf = tmp_path / "prep.gtf"
    gtf.write_text(
        'chr1\tM\texon\t101\t200\t.\t+\t.\ttranscript_id "t1";\n'
        'chr1\tM\texon\t301\t400\t.\t+\t.\ttranscript_id "t1";\n'
        'chr2\tM\texon\t501\t600\t.\t-\t.\ttranscript_id "t2";\n'
    )
    out = pl._parse_prepared_gtf(str(gtf))
    assert out["t1"] == ("chr1", "+", [Exon(100, 200), Exon(300, 400)])
    assert out["t2"] == ("chr2", "-", [Exon(500, 600)])


def test_build_external_scores_no_h5_defaults_zero(tmp_path):
    gtf = tmp_path / "prep.gtf"
    gtf.write_text('chr1\tM\texon\t101\t200\t.\t+\t.\ttranscript_id "t1";\n')
    rows = pl._build_external_scores(str(gtf), None, {}, 0.0)
    assert rows["t1"] == {"helixer_support": 0.0, "helixer_locus_conf": 0.0, "tpm": 0.0}


def test_build_external_scores_tpm_from_struct(tmp_path):
    gtf = tmp_path / "prep.gtf"
    gtf.write_text('chr1\tM\texon\t101\t200\t.\t+\t.\ttranscript_id "t1";\n')
    struct_tpm = {"chr1:+:((100, 200),)": 8.0}
    rows = pl._build_external_scores(str(gtf), None, struct_tpm, 16.0)
    assert rows["t1"]["tpm"] == pytest.approx(0.5)


def test_merge_star_junctions_unions_samples(tmp_path):
    line = "chr1\t101\t200\t1\t0\t0\t5\t0\t30\n"
    a = tmp_path / "a_SJ.out.tab"
    b = tmp_path / "b_SJ.out.tab"
    a.write_text(line)
    b.write_text(line)
    merged = pl._merge_star_junctions([str(a), str(b)])
    assert len(merged) == 1
    j = merged[0]
    assert (j.seqid, j.donor, j.acceptor, j.strand) == ("chr1", 100, 200, "+")
    assert j.read_count == 10 and j.samples == 2


# ---------------------------------------------------------------------------
# Basic mode (Helixer + genome only → all backstop)
# ---------------------------------------------------------------------------

def test_basic_mode_returns_backstop_genes(basic_config):
    genes = run_pipeline(basic_config)
    assert len(genes) == 3
    assert all(g.origin == "helixer_backstop" for g in genes)


def test_basic_mode_gene_ids_anchored(basic_config):
    genes = run_pipeline(basic_config)
    assert all(g.gene_id.startswith("HFG_") for g in genes)


def test_basic_mode_writes_full_gff3(basic_config):
    run_pipeline(basic_config)
    gff = Path(f"{basic_config.output_prefix}.gff3")
    assert gff.exists()
    assert gff.read_text().startswith("##gff-version 3")


def test_basic_mode_writes_tier_files(basic_config):
    run_pipeline(basic_config)
    for n in (1, 2, 3):
        assert Path(f"{basic_config.output_prefix}.tier{n}.gff3").exists()


def test_basic_mode_multiexon_minus_junction_flag(basic_config):
    genes = run_pipeline(basic_config)
    g2 = next(g for g in genes if g.strand == "-")
    # multi-exon backstop with no junction evidence → flagged, fraction 0
    assert any(f.name == "NO_JUNCTION_SUPPORT" for f in g2.flags)
    assert g2.transcripts[0].junction_support_fraction == 0.0


def test_basic_mode_report_columns_and_rows(basic_config):
    genes = run_pipeline(basic_config)
    lines = Path(basic_config.report_path).read_text().splitlines()
    header = lines[0].split("\t")
    assert tuple(header) == pl._REPORT_COLUMNS
    assert len(lines) == 1 + len(genes)  # header + one row per gene


def test_basic_mode_report_one_row_per_gene(basic_config):
    genes = run_pipeline(basic_config)
    lines = Path(basic_config.report_path).read_text().splitlines()[1:]
    reported = {ln.split("\t")[0] for ln in lines}
    assert reported == {g.gene_id for g in genes}


def test_basic_mode_id_map_persisted(basic_config):
    run_pipeline(basic_config)
    p = Path(basic_config.id_map_path)
    assert p.exists()
    import json
    mapping = json.loads(p.read_text())
    assert all(v.startswith("HFG_") for v in mapping.values())


def test_id_map_stable_across_runs(basic_config):
    genes1 = run_pipeline(basic_config)
    ids1 = {(g.seqid, g.start): g.gene_id for g in genes1}
    genes2 = run_pipeline(basic_config)  # reuses persisted id_map
    ids2 = {(g.seqid, g.start): g.gene_id for g in genes2}
    assert ids1 == ids2


def test_region_filter_restricts_scaffold(tmp_path, genome_fasta, helixer_gff3):
    cfg = PipelineConfig(
        genome_fasta=genome_fasta, helixer_gff3=helixer_gff3,
        output_prefix=str(tmp_path / "out"), work_dir=str(tmp_path / "work"),
        region="chr2",
    )
    genes = run_pipeline(cfg)
    assert len(genes) == 1
    assert genes[0].seqid == "chr2"


def test_missing_required_input_raises(tmp_path, helixer_gff3):
    cfg = PipelineConfig(genome_fasta=str(tmp_path / "nope.fa"),
                         helixer_gff3=helixer_gff3,
                         output_prefix=str(tmp_path / "out"))
    with pytest.raises(FileNotFoundError):
        run_pipeline(cfg)


# ---------------------------------------------------------------------------
# Full mode (mocked Mikado end-to-end)
# ---------------------------------------------------------------------------

def test_full_mode_produces_mikado_gene(full_config, mock_mikado):
    genes = run_pipeline(full_config)
    origins = {g.origin for g in genes}
    assert "mikado_1to1" in origins
    # g1 reconciled with Mikado; g2/g3 carried through as backstop
    assert len(genes) == 3


def test_full_mode_mikado_gene_tier1(full_config, mock_mikado):
    genes = run_pipeline(full_config)
    mik = next(g for g in genes if g.origin == "mikado_1to1")
    assert mik.tier == 1                       # coding + blast homology
    assert mik.transcripts[0].cds is not None


def test_full_mode_tier1_file_has_mikado_gene(full_config, mock_mikado):
    genes = run_pipeline(full_config)
    mik = next(g for g in genes if g.origin == "mikado_1to1")
    tier1 = Path(f"{full_config.output_prefix}.tier1.gff3").read_text()
    assert mik.gene_id in tier1


def test_full_mode_writes_external_scores(full_config, mock_mikado):
    run_pipeline(full_config)
    ext = Path(full_config.work_dir) / "mikado_inputs" / "external_scores.tsv"
    assert ext.exists()
    header = ext.read_text().splitlines()[0].split("\t")
    assert header == ["tid", "helixer_support", "helixer_locus_conf", "tpm"]


def test_full_mode_writes_report_and_gff3(full_config, mock_mikado):
    genes = run_pipeline(full_config)
    assert Path(f"{full_config.output_prefix}.gff3").exists()
    rows = Path(full_config.report_path).read_text().splitlines()
    assert len(rows) == 1 + len(genes)


def test_full_mode_report_has_cds_flag(full_config, mock_mikado):
    genes = run_pipeline(full_config)
    all_lines = Path(full_config.report_path).read_text().splitlines()
    # Phase 30 added a "biotype" report column (shifts has_cds), so resolve the
    # column by header name rather than a hardcoded index.
    header = all_lines[0].split("\t")
    cds_idx = header.index("has_cds")
    has_cds = dict(
        (ln.split("\t")[0], ln.split("\t")[cds_idx]) for ln in all_lines[1:]
    )
    mik = next(g for g in genes if g.origin == "mikado_1to1")
    assert has_cds[mik.gene_id] == "true"


def test_tier_files_are_cumulative(full_config, mock_mikado):
    run_pipeline(full_config)
    t1 = Path(f"{full_config.output_prefix}.tier1.gff3").read_text()
    t2 = Path(f"{full_config.output_prefix}.tier2.gff3").read_text()
    t3 = Path(f"{full_config.output_prefix}.tier3.gff3").read_text()
    # every tier-1 gene id also appears in the cumulative tier-2 and tier-3 files
    for gid in [ln.split("ID=")[1].split(";")[0]
                for ln in t1.splitlines() if "\tgene\t" in ln]:
        assert gid in t2 and gid in t3
