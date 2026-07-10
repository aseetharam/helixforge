"""CLI plumbing for the v1-style evidence inputs (repeatable + comma + ``-list``).

Exercises the ``helixforge evidence`` flags through ``CliRunner`` without
decoding real BAM/SJ: a fake ``score_annotation`` records the merged path lists
the CLI passes, so we assert the expansion/merge wiring (not the scoring). The
``--stringtie`` deprecation shim is covered too.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

import helixforge.cli as cli
from helixforge.cli import main


@pytest.fixture
def gff3(tmp_path):
    p = tmp_path / "ann.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "chr1\t.\texon\t1\t100\t.\t+\t.\tID=g1.t1.e1;Parent=g1.t1\n"
    )
    return str(p)


@pytest.fixture
def captured(monkeypatch):
    """Patch evidence.score_annotation to capture the merged inputs (no decode)."""
    calls: dict[str, object] = {}

    def _fake_score_annotation(gff3_path, **kwargs):
        import pandas as pd

        calls["gff3_path"] = gff3_path
        calls.update(kwargs)
        return pd.DataFrame(columns=["transcript_id"])

    # The command imports these names locally from helixforge.score.evidence.
    import helixforge.score.evidence as ev

    monkeypatch.setattr(ev, "score_annotation", _fake_score_annotation)
    monkeypatch.setattr(ev, "write_evidence_tsv", lambda df, path: None)
    monkeypatch.setattr(ev, "summarize_evidence", lambda df: {})
    return calls


def _bam(tmp_path, name):
    p = tmp_path / name
    p.write_text("x")
    return str(p)


# ---------------------------------------------------------------------------
# repeated + comma forms
# ---------------------------------------------------------------------------


def test_bam_repeated_and_sj_comma(tmp_path, gff3, captured):
    a, b = _bam(tmp_path, "a.bam"), _bam(tmp_path, "b.bam")
    s1, s2 = _bam(tmp_path, "s1.SJ.out.tab"), _bam(tmp_path, "s2.SJ.out.tab")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3,
        "--bam", a, "--bam", b,          # repeated
        "--sj", f"{s1},{s2}",            # comma
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["bam_paths"] == [a, b]
    assert captured["star_sj_paths"] == [s1, s2]


def test_stringtie_repeated_and_comma(tmp_path, gff3, captured):
    g1, g2, g3 = (_bam(tmp_path, f"{n}.gtf") for n in ("g1", "g2", "g3"))
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3,
        "--stringtie", g1, "--stringtie", f"{g2},{g3}",
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["stringtie_gtfs"] == [g1, g2, g3]


# ---------------------------------------------------------------------------
# -list companions and merge
# ---------------------------------------------------------------------------


def test_bam_list_merges_with_direct_flag(tmp_path, gff3, captured):
    a, b, c = (_bam(tmp_path, f"{n}.bam") for n in ("a", "b", "c"))
    fofn = tmp_path / "bams.list"
    fofn.write_text(f"{b}\n{c}\n")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3,
        "--bam", a, "--bam-list", str(fofn),
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["bam_paths"] == [a, b, c]


def test_list_only_invocation(tmp_path, gff3, captured):
    s1, s2 = _bam(tmp_path, "s1.SJ.out.tab"), _bam(tmp_path, "s2.SJ.out.tab")
    fofn = tmp_path / "sj.list"
    fofn.write_text(f"{s1}\n{s2}\n")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3,
        "--sj-list", str(fofn),
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["star_sj_paths"] == [s1, s2]


def test_stringtie_list_companion(tmp_path, gff3, captured):
    g1, g2 = _bam(tmp_path, "g1.gtf"), _bam(tmp_path, "g2.gtf")
    fofn = tmp_path / "st.list"
    fofn.write_text(f"{g1}\n{g2}\n")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3,
        "--stringtie-list", str(fofn),
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["stringtie_gtfs"] == [g1, g2]


# ---------------------------------------------------------------------------
# --help and the deprecation shim
# ---------------------------------------------------------------------------


def test_new_list_flags_in_help():
    res = CliRunner().invoke(main, ["evidence", "--help"])
    assert res.exit_code == 0
    for flag in ("--bam-list", "--sj-list", "--stringtie-list"):
        assert flag in res.output


def test_stringtie_deprecation_shim_warns_and_works(tmp_path, gff3, captured):
    g1, g2 = _bam(tmp_path, "g1.gtf"), _bam(tmp_path, "g2.gtf")
    legacy = tmp_path / "samples.txt"  # non-GTF → looks like a FOFN
    legacy.write_text(f"{g1}\n{g2}\n")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3,
        "--stringtie", str(legacy),
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert "deprecated" in res.output
    assert "--stringtie-list" in res.output
    assert captured["stringtie_gtfs"] == [g1, g2]  # still resolves the GTFs


def test_no_evidence_source_errors(tmp_path, gff3, captured):
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3, "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code != 0
    assert "at least one evidence source" in res.output


# ---------------------------------------------------------------------------
# protein axis + threads
# ---------------------------------------------------------------------------


def test_proteins_repeated_and_list_merge(tmp_path, gff3, captured):
    p1, p2, p3 = (_bam(tmp_path, f"{n}.faa") for n in ("p1", "p2", "p3"))
    mp = _bam(tmp_path, "aln.gff")
    fofn = tmp_path / "prot.list"
    fofn.write_text(f"{p2}\n{p3}\n")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3,
        "--proteins", p1, "--proteins-list", str(fofn),
        "--miniprot-gff", mp,  # precomputed -> no --genome needed
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["proteins"] == [p1, p2, p3]
    assert captured["miniprot_gff"] == mp


def test_proteins_without_genome_or_gff_errors(tmp_path, gff3, captured):
    p1 = _bam(tmp_path, "p1.faa")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3, "--proteins", p1,
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code != 0
    assert "--proteins needs" in res.output


def test_proteins_with_genome_ok(tmp_path, gff3, captured):
    p1 = _bam(tmp_path, "p1.faa")
    genome = _bam(tmp_path, "genome.fa")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3, "--proteins", p1, "--genome", genome,
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["proteins"] == [p1]
    assert captured["genome"] == genome


def test_threads_flag_passed_through(tmp_path, gff3, captured):
    a = _bam(tmp_path, "a.bam")
    res = CliRunner().invoke(main, [
        "evidence", "--gff3", gff3, "--bam", a, "-j", "8",
        "--out", str(tmp_path / "out.tsv"),
    ])
    assert res.exit_code == 0, res.output
    assert captured["threads"] == 8


def test_protein_and_threads_in_help():
    res = CliRunner().invoke(main, ["evidence", "--help"])
    assert res.exit_code == 0
    for flag in ("--proteins", "--proteins-list", "--miniprot-gff", "--threads"):
        assert flag in res.output


# ---------------------------------------------------------------------------
# pipeline_options share the same convention (reconcile builds a PipelineConfig)
# ---------------------------------------------------------------------------


def test_pipeline_options_expand_via_build_config(tmp_path):
    a, b = _bam(tmp_path, "a.bam"), _bam(tmp_path, "b.bam")
    g1, g2 = _bam(tmp_path, "g1.gtf"), _bam(tmp_path, "g2.gtf")
    fofn = tmp_path / "bams.list"
    fofn.write_text(f"{b}\n")
    cfg = cli._build_config({
        "genome_fasta": _bam(tmp_path, "genome.fa"),
        "helixer_gff3": _bam(tmp_path, "helixer.gff3"),
        "bam": (a,), "bam_list": str(fofn),
        "stringtie": (f"{g1},{g2}",), "stringtie_list": None,
        "star_sj": (), "star_sj_list": None,
        "bigwig": (),
    })
    assert cfg.bam_paths == [a, b]
    assert cfg.stringtie_list == [g1, g2]


# ---------------------------------------------------------------------------
# --bigwig multi-file convention (comma + repeated + --bigwig-list)
# ---------------------------------------------------------------------------


def test_bigwig_comma_separated(tmp_path):
    a, b = _bam(tmp_path, "a.bw"), _bam(tmp_path, "b.bw")
    cfg = cli._build_config({
        "genome_fasta": _bam(tmp_path, "genome.fa"),
        "helixer_gff3": _bam(tmp_path, "helixer.gff3"),
        "bigwig": (f"{a},{b}",), "bigwig_list": None,
        "bam": (), "bam_list": None,
        "stringtie": (), "stringtie_list": None,
        "star_sj": (), "star_sj_list": None,
    })
    assert cfg.bigwig_paths == [a, b]


def test_bigwig_repeated(tmp_path):
    a, b = _bam(tmp_path, "a.bw"), _bam(tmp_path, "b.bw")
    cfg = cli._build_config({
        "genome_fasta": _bam(tmp_path, "genome.fa"),
        "helixer_gff3": _bam(tmp_path, "helixer.gff3"),
        "bigwig": (a, b), "bigwig_list": None,
        "bam": (), "bam_list": None,
        "stringtie": (), "stringtie_list": None,
        "star_sj": (), "star_sj_list": None,
    })
    assert cfg.bigwig_paths == [a, b]


def test_bigwig_list_fofn(tmp_path):
    a, b = _bam(tmp_path, "a.bw"), _bam(tmp_path, "b.bw")
    fofn = tmp_path / "bw.list"
    fofn.write_text(f"{a}\n{b}\n")
    cfg = cli._build_config({
        "genome_fasta": _bam(tmp_path, "genome.fa"),
        "helixer_gff3": _bam(tmp_path, "helixer.gff3"),
        "bigwig": (), "bigwig_list": str(fofn),
        "bam": (), "bam_list": None,
        "stringtie": (), "stringtie_list": None,
        "star_sj": (), "star_sj_list": None,
    })
    assert cfg.bigwig_paths == [a, b]


def test_bigwig_merge_direct_and_list(tmp_path):
    a, b, c = (_bam(tmp_path, f"{n}.bw") for n in ("a", "b", "c"))
    fofn = tmp_path / "bw.list"
    fofn.write_text(f"{b}\n{c}\n")
    cfg = cli._build_config({
        "genome_fasta": _bam(tmp_path, "genome.fa"),
        "helixer_gff3": _bam(tmp_path, "helixer.gff3"),
        "bigwig": (a,), "bigwig_list": str(fofn),
        "bam": (), "bam_list": None,
        "stringtie": (), "stringtie_list": None,
        "star_sj": (), "star_sj_list": None,
    })
    assert cfg.bigwig_paths == [a, b, c]


def test_bigwig_dedup_across_forms(tmp_path):
    a, b = _bam(tmp_path, "a.bw"), _bam(tmp_path, "b.bw")
    fofn = tmp_path / "bw.list"
    fofn.write_text(f"{a}\n")
    cfg = cli._build_config({
        "genome_fasta": _bam(tmp_path, "genome.fa"),
        "helixer_gff3": _bam(tmp_path, "helixer.gff3"),
        "bigwig": (a, b), "bigwig_list": str(fofn),
        "bam": (), "bam_list": None,
        "stringtie": (), "stringtie_list": None,
        "star_sj": (), "star_sj_list": None,
    })
    assert cfg.bigwig_paths == [a, b]


def test_bigwig_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="bigWig file"):
        cli._build_config({
            "genome_fasta": _bam(tmp_path, "genome.fa"),
            "helixer_gff3": _bam(tmp_path, "helixer.gff3"),
            "bigwig": (str(tmp_path / "nope.bw"),), "bigwig_list": None,
            "bam": (), "bam_list": None,
            "stringtie": (), "stringtie_list": None,
            "star_sj": (), "star_sj_list": None,
        })


def test_bigwig_flags_in_reconcile_help():
    res = CliRunner().invoke(main, ["reconcile", "--help"])
    assert res.exit_code == 0
    assert "--bigwig" in res.output
    assert "--bigwig-list" in res.output
