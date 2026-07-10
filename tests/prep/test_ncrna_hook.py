"""Tests for prep/ncrna.py — Phase 30 D4 structured-ncRNA hook (§2.6). Floor: 1.

The hook is opt-in (off by default) and every external tool is a subprocess CLI
(argv-as-list, mocked here — no real tRNAscan-SE / Infernal in CI). Verifies the
disabled default, the exact argv, and that parsed loci become biotyped genes
(both strands).
"""

from unittest.mock import patch

import pytest

from helixforge.prep import ncrna


def test_disabled_by_default_runs_no_subprocess():
    # enabled=False -> returns [] and never invokes the tool.
    with patch("helixforge.prep.ncrna.run_tool") as run:
        out = ncrna.scan_structured_ncrna("genome.fa", "/tmp/wd", enabled=False)
    assert out == []
    run.assert_not_called()


def test_trnascan_argv_is_a_list_with_expected_flags():
    calls = []
    with patch("helixforge.prep.ncrna.run_tool",
               side_effect=lambda argv, **k: calls.append([str(a) for a in argv])):
        with patch("helixforge.prep.ncrna.output_is_fresh", return_value=False):
            ncrna.run_trnascan_se("genome.fa", "out.txt", threads=4)
    assert calls == [["tRNAscan-SE", "-E", "--thread", "4", "-o", "out.txt", "genome.fa"]]


def test_infernal_cmscan_argv():
    calls = []
    with patch("helixforge.prep.ncrna.run_tool",
               side_effect=lambda argv, **k: calls.append([str(a) for a in argv])):
        with patch("helixforge.prep.ncrna.output_is_fresh", return_value=False):
            ncrna.run_infernal_cmscan("genome.fa", "Rfam.cm", "hits.tbl",
                                      clanin="Rfam.clanin", threads=8)
    argv = calls[0]
    assert argv[0] == "cmscan"
    assert "--cut_ga" in argv and "--rfam" in argv
    assert argv[argv.index("--clanin") + 1] == "Rfam.clanin"
    assert argv[argv.index("--cpu") + 1] == "8"
    assert argv[argv.index("--tblout") + 1] == "hits.tbl"
    assert argv[-2:] == ["Rfam.cm", "genome.fa"]


def test_infernal_requires_rfam_cm():
    with pytest.raises(ValueError, match="rfam_cm"):
        ncrna.scan_structured_ncrna("genome.fa", "/tmp/wd", enabled=True, tool="infernal")


def test_parse_trnascan_both_strands(tmp_path):
    tbl = (
        "Sequence\ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote\n"
        "Name\t\t\t\t\t\t\t\t\t\n"
        "--------\t----\t-----\t----\t----\t-----\t----\t----\t------\t----\n"
        "chr1\t1\t100\t172\tLys\tCTT\t0\t0\t75.3\t\n"
        "chr2\t1\t500\t428\tPro\tAGG\t0\t0\t60.0\tpseudo\n"
    )
    p = tmp_path / "trnascan.out"
    p.write_text(tbl)
    recs = ncrna.parse_trnascan(p)
    # plus strand: 1-based 100..172 -> internal (99, 172); minus: 500>428 -> (427, 500)
    assert recs == [
        ("chr1", 99, 172, "+", "tRNA"),
        ("chr2", 427, 500, "-", "pseudogene"),
    ]


def test_build_ncrna_genes_assigns_biotype_and_ids(tmp_path):
    recs = [
        ("chr1", 99, 172, "+", "tRNA"),
        ("chr1", 427, 500, "-", "rRNA"),
    ]
    genes = ncrna.build_ncrna_genes(recs)
    assert [g.gene_id for g in genes] == ["HFG_95000", "HFG_95001"]
    assert [g.biotype for g in genes] == ["tRNA", "rRNA"]
    assert [g.transcripts[0].biotype for g in genes] == ["tRNA", "rRNA"]
    assert {g.strand for g in genes} == {"+", "-"}
    assert all(g.origin == "novel" and g.tier == 4 for g in genes)
    assert all(g.transcripts[0].source == "ncrna" for g in genes)


def test_scan_trnascan_end_to_end_mocked(tmp_path):
    tbl = (
        "Sequence\ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote\n"
        "Name\t\t\t\t\t\t\t\t\t\n"
        "--------\t----\t-----\t----\t----\t-----\t----\t----\t------\t----\n"
        "chr1\t1\t100\t172\tLys\tCTT\t0\t0\t75.3\t\n"
    )

    def fake_run(argv, **kwargs):
        # tRNAscan writes to the -o path; emulate it.
        out = argv[argv.index("-o") + 1]
        from pathlib import Path
        Path(out).write_text(tbl)

    with patch("helixforge.prep.ncrna.run_tool", side_effect=fake_run):
        genes = ncrna.scan_structured_ncrna(
            "genome.fa", tmp_path / "wd", enabled=True, tool="trnascan",
        )
    assert len(genes) == 1
    assert genes[0].biotype == "tRNA"
