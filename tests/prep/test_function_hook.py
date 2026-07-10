"""Phase 31 D1/D2 — functional-annotation hook (InterProScan / eggNOG).

Concrete literal coordinates / values only (CLAUDE.md §12). The external tools
are NEVER run in CI — every run wrapper is exercised with ``run_tool`` monkey-
patched, asserting the argv is a list with the expected flags (CLAUDE.md §C1).
"""

from __future__ import annotations

import attrs
import pytest

from helixforge.prep import function as fn
from helixforge.prep.function import (
    FunctionalRecord,
    annotate_function,
    apply_domain_credibility,
    parse_eggnog_annotations,
    parse_interproscan_tsv,
)
from helixforge.qc.flags import DOMAIN_COMPLETE
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)

_IPR_COLS = 15


def _ipr_row(tid, analysis, sig, start, stop, seq_len=200, interpro="-", go="-"):
    cols = ["-"] * _IPR_COLS
    cols[0] = tid
    cols[2] = str(seq_len)
    cols[3] = analysis
    cols[4] = sig
    cols[6] = str(start)
    cols[7] = str(stop)
    cols[11] = interpro
    cols[13] = go
    return "\t".join(cols)


def _gene(gid, strand="+", cds=True, biotype="protein_coding"):
    exons = [Exon(1000, 1200), Exon(1300, 1500), Exon(1600, 1800)]
    cdsseg = (
        [CDSSegment(1050, 1200, 0), CDSSegment(1300, 1500, 0), CDSSegment(1600, 1649, 1)]
        if cds else None
    )
    tx = TranscriptCandidate(
        transcript_id=f"{gid}.1", locus_id=gid, source="mikado", seqid="chr1",
        start=1000, end=1800, strand=strand, exons=exons, cds=cdsseg,
        tpm=5.0, junction_support_fraction=1.0, confidence=0.9,
        combined_score=10.0, is_primary=True, biotype=biotype,
    )
    return ReconciledGene(
        gene_id=gid, seqid="chr1", start=1000, end=1800, strand=strand, tier=1,
        transcripts=[tx], primary_transcript_id=f"{gid}.1",
        classification=LocusClassification(locus_id=gid, status="EXPRESSED", max_tpm=5.0),
        origin="mikado_1to1", biotype=biotype,
    )


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def test_parse_interproscan_go_and_dbxref(tmp_path):
    p = tmp_path / "ipr.tsv"
    p.write_text(
        _ipr_row("HFG_00001.1", "Pfam", "PF00069", 30, 150,
                 interpro="IPR000719", go="GO:0004672|GO:0005524") + "\n"
    )
    rec = parse_interproscan_tsv(p)["HFG_00001.1"]
    assert rec.go_terms == ("GO:0004672", "GO:0005524")
    assert rec.dbxrefs == ("InterPro:IPR000719", "Pfam:PF00069")


def test_parse_interproscan_domain_complete_contained(tmp_path):
    # A Pfam match [30,150] fully inside a 200-aa protein → domain_complete True.
    p = tmp_path / "ipr.tsv"
    p.write_text(_ipr_row("HFG_00001.1", "Pfam", "PF00069", 30, 150, seq_len=200) + "\n")
    assert parse_interproscan_tsv(p)["HFG_00001.1"].domain_complete is True


def test_parse_interproscan_domain_not_complete_when_overrun(tmp_path):
    # Match stop (280) exceeds the protein length (200) → not fully contained.
    p = tmp_path / "ipr.tsv"
    p.write_text(_ipr_row("HFG_00001.1", "Pfam", "PF00069", 30, 280, seq_len=200) + "\n")
    assert parse_interproscan_tsv(p)["HFG_00001.1"].domain_complete is False


def test_parse_eggnog_annotations(tmp_path):
    p = tmp_path / "egg.annotations"
    p.write_text(
        "## emapper run\n"
        "#query\tseed_ortholog\tevalue\tscore\tGOs\tPFAMs\n"
        "HFG_00002.1\tsp|Q12345\t1e-10\t100\tGO:0008150,GO:0003674\tPF00001,PF00002\n"
    )
    rec = parse_eggnog_annotations(p)["HFG_00002.1"]
    assert rec.go_terms == ("GO:0003674", "GO:0008150")
    assert "Pfam:PF00001" in rec.dbxrefs and "UniProt:sp|Q12345" in rec.dbxrefs
    assert rec.domain_complete is True  # PFAMs present


# ---------------------------------------------------------------------------
# Run wrappers — argv asserted, tool mocked
# ---------------------------------------------------------------------------


def test_interproscan_argv_is_list_with_flags(tmp_path, monkeypatch):
    captured = {}

    def fake_run_tool(argv, **kw):
        captured["argv"] = list(argv)
        # write a minimal TSV so the parser has something
        out = argv[argv.index("-o") + 1]
        open(out, "w").write(_ipr_row("p.1", "Pfam", "PF00069", 1, 50) + "\n")

    monkeypatch.setattr(fn, "run_tool", fake_run_tool)
    proteins = tmp_path / "proteins.fa"
    proteins.write_text(">p.1\nMKPGF\n")
    out = fn.run_interproscan(proteins, tmp_path / "ipr.tsv", applications=["Pfam"])
    argv = captured["argv"]
    assert isinstance(argv, list)
    assert argv[0] == "interproscan.sh"
    assert "-f" in argv and argv[argv.index("-f") + 1] == "TSV"
    assert "-appl" in argv and argv[argv.index("-appl") + 1] == "Pfam"
    assert out.exists()


def test_eggnog_argv_is_list_with_flags(tmp_path, monkeypatch):
    captured = {}

    def fake_run_tool(argv, **kw):
        captured["argv"] = list(argv)

    monkeypatch.setattr(fn, "run_tool", fake_run_tool)
    proteins = tmp_path / "proteins.fa"
    proteins.write_text(">p.1\nMKPGF\n")
    fn.run_eggnog_mapper(proteins, tmp_path / "egg", db="/data/eggnog")
    argv = captured["argv"]
    assert argv[0] == "emapper.py"
    assert "--itype" in argv and argv[argv.index("--itype") + 1] == "proteins"
    assert "--data_dir" in argv and argv[argv.index("--data_dir") + 1] == "/data/eggnog"


def test_annotate_function_off_by_default_runs_nothing(tmp_path, monkeypatch):
    # enabled=False short-circuits: no subprocess, empty records.
    def boom(*a, **k):
        raise AssertionError("run_tool must not be called when disabled")

    monkeypatch.setattr(fn, "run_tool", boom)
    proteins = tmp_path / "proteins.fa"
    proteins.write_text(">p.1\nMKPGF\n")
    out = annotate_function([_gene("HFG_00001")], proteins,
                            out_dir=tmp_path / "fn", enabled=False)
    assert out == {}


def test_annotate_function_interproscan_end_to_end(tmp_path, monkeypatch):
    def fake_run_tool(argv, **kw):
        out = argv[argv.index("-o") + 1]
        open(out, "w").write(
            _ipr_row("HFG_00001.1", "Pfam", "PF00069", 10, 120, seq_len=130,
                     interpro="IPR000719", go="GO:0004672") + "\n"
        )

    monkeypatch.setattr(fn, "run_tool", fake_run_tool)
    proteins = tmp_path / "proteins.fa"
    proteins.write_text(">HFG_00001.1\nMKPGF\n")
    records = annotate_function([_gene("HFG_00001")], proteins,
                                out_dir=tmp_path / "fn", tool="interproscan", enabled=True)
    assert "HFG_00001.1" in records
    assert records["HFG_00001.1"].domain_complete is True


def test_annotate_function_rejects_unknown_tool(tmp_path):
    with pytest.raises(ValueError):
        annotate_function([], tmp_path / "p.fa", out_dir=tmp_path, tool="blastp")


# ---------------------------------------------------------------------------
# D2 — domain credibility flag (count-neutral: INFO flag only)
# ---------------------------------------------------------------------------


def test_apply_domain_credibility_flags_complete_both_strands():
    for strand in ("+", "-"):
        g = _gene("HFG_00007", strand=strand)
        func = {"HFG_00007.1": FunctionalRecord("HFG_00007.1", domain_complete=True)}
        out = apply_domain_credibility([g], func)
        assert DOMAIN_COMPLETE in out[0].flags


def test_apply_domain_credibility_empty_map_is_identity():
    g = _gene("HFG_00008")
    assert apply_domain_credibility([g], {}) == [g]


def test_apply_domain_credibility_no_flag_when_incomplete():
    g = _gene("HFG_00009")
    func = {"HFG_00009.1": FunctionalRecord("HFG_00009.1", domain_complete=False)}
    out = apply_domain_credibility([g], func)
    assert DOMAIN_COMPLETE not in out[0].flags
    # tier / biotype untouched (count-neutral)
    assert out[0].tier == 1 and out[0].biotype == "protein_coding"
