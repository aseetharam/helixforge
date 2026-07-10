"""Phase 41 — ``helixforge utils`` CLI subcommand group.

Tests that the group and its subcommands parse options correctly and fail on
missing required options. (``validate-inputs`` was folded into ``helixforge
doctor`` — input validation is a single preflight command.)
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from helixforge.cli import main


@pytest.fixture
def runner():
    return CliRunner()


# ---------------------------------------------------------------------------
# Group help
# ---------------------------------------------------------------------------

_SUBCOMMANDS = [
    "fetch-db",
    "extract-proteins",
    "align",
    "qc",
    "convert",
    "filter",
    "summarize",
]


def test_utils_help(runner):
    result = runner.invoke(main, ["utils", "--help"])
    assert result.exit_code == 0
    for name in _SUBCOMMANDS:
        assert name in result.output, f"{name} missing from utils --help"


# ---------------------------------------------------------------------------
# Per-subcommand --help
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("subcmd", _SUBCOMMANDS)
def test_subcommand_help(runner, subcmd):
    result = runner.invoke(main, ["utils", subcmd, "--help"])
    assert result.exit_code == 0, result.output


# ---------------------------------------------------------------------------
# fetch-db --list
# ---------------------------------------------------------------------------

def test_fetch_db_list(runner):
    result = runner.invoke(main, ["utils", "fetch-db", "--list"])
    assert result.exit_code == 0
    assert "swissprot" in result.output
    assert "Description" in result.output


def test_fetch_db_missing_db(runner):
    result = runner.invoke(main, ["utils", "fetch-db"])
    assert result.exit_code != 0


# ---------------------------------------------------------------------------
# extract-proteins missing required
# ---------------------------------------------------------------------------

def test_extract_proteins_missing_gff3(runner):
    result = runner.invoke(main, ["utils", "extract-proteins",
                                  "--genome", "g.fa", "--out", "p.fa"])
    assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Stub invocations — each prints "Not yet implemented"
# ---------------------------------------------------------------------------

def _touch(p: Path, content: str = "") -> str:
    p.write_text(content)
    return str(p)


def test_fetch_db_custom(runner, tmp_path):
    db_fa = _touch(tmp_path / "proteins.fa", ">p\nMKL\n")
    from unittest import mock
    with mock.patch("subprocess.run") as mock_run:
        mock_run.return_value = mock.MagicMock(returncode=0)
        result = runner.invoke(main, ["utils", "fetch-db", "--db", db_fa])
    assert result.exit_code == 0
    assert "Custom: proteins.fa" in result.output


def test_extract_proteins_cli(runner, tmp_path):
    gff3_text = (
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t12\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\t.\texon\t1\t12\t.\t+\t.\tID=t1.exon1;Parent=t1\n"
        "chr1\t.\tCDS\t1\t12\t.\t+\t0\tID=t1.CDS1;Parent=t1\n"
    )
    gff3 = _touch(tmp_path / "genes.gff3", gff3_text)
    genome = _touch(tmp_path / "genome.fa", ">chr1\nATGTTTCCCGGGAAAAAAAAAAAA\n")
    out = str(tmp_path / "proteins.fa")
    result = runner.invoke(main, ["utils", "extract-proteins",
                                  "--gff3", gff3, "--genome", genome,
                                  "--out", out, "--min-length", "1"])
    assert result.exit_code == 0
    assert "Extracted" in result.output


def test_align_implemented(runner, tmp_path):
    genome = _touch(tmp_path / "genome.fa", ">chr1\nACGT\n")
    proteins = _touch(tmp_path / "proteins.fa", ">p\nMKL\n")
    out = str(tmp_path / "out.gff")
    from unittest import mock
    with mock.patch("helixforge.prep.protein_align.run_tool") as mock_run, \
         mock.patch("helixforge.prep.protein_align.output_is_fresh", return_value=False):
        mock_run.return_value = mock.MagicMock(returncode=0)
        result = runner.invoke(main, ["utils", "align",
                                      "--genome", genome, "--proteins", proteins,
                                      "--out", out])
    assert result.exit_code == 0
    assert "miniprot" in result.output


def test_qc_implemented(runner, tmp_path):
    gff3 = _touch(tmp_path / "genes.gff3",
                   "##gff-version 3\n"
                   "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1;tier=1;origin=test\n"
                   "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
                   "chr1\t.\texon\t1\t100\t.\t+\t.\tID=t1.exon1;Parent=t1\n")
    out = str(tmp_path / "report.json")
    result = runner.invoke(main, ["utils", "qc",
                                  "--gff3", gff3, "--out", out,
                                  "--format", "json"])
    assert result.exit_code == 0
    assert "1 genes" in result.output


def test_validate_inputs_command_removed(runner, tmp_path):
    # Folded into `helixforge doctor`; no longer a utils subcommand.
    genome = _touch(tmp_path / "genome.fa", ">chr1\nACGT\n")
    result = runner.invoke(main, ["utils", "validate-inputs",
                                  "--genome", genome])
    assert result.exit_code != 0


def test_convert_implemented(runner, tmp_path):
    inp = _touch(tmp_path / "genes.gff3",
                 "##gff-version 3\nchr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                 "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
                 "chr1\t.\texon\t1\t100\t.\t+\t.\tID=t1.exon1;Parent=t1\n")
    out = str(tmp_path / "genes.gtf")
    result = runner.invoke(main, ["utils", "convert",
                                  "--input", inp, "--out", out])
    assert result.exit_code == 0
    assert "1 genes" in result.output


def test_filter_implemented(runner, tmp_path):
    gff3 = _touch(tmp_path / "genes.gff3",
                   "##gff-version 3\n"
                   "chr1\tHF\tgene\t1\t300\t.\t+\t.\tID=g1;tier=1;flags=\n"
                   "chr1\tHF\tmRNA\t1\t300\t.\t+\t.\tID=g1.1;Parent=g1\n"
                   "chr1\tHF\texon\t1\t300\t.\t+\t.\tID=g1.1.exon1;Parent=g1.1\n")
    out = str(tmp_path / "filtered.gff3")
    result = runner.invoke(main, ["utils", "filter",
                                  "--gff3", gff3, "--out", out])
    assert result.exit_code == 0
    assert "Kept 1" in result.output


def test_summarize_implemented(runner, tmp_path):
    gff3 = _touch(tmp_path / "genes.gff3",
                   "##gff-version 3\n"
                   "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                   "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
                   "chr1\t.\texon\t1\t100\t.\t+\t.\tID=t1.exon1;Parent=t1\n")
    result = runner.invoke(main, ["utils", "summarize",
                                  "--gff3", gff3])
    assert result.exit_code == 0
    assert "gene_count" in result.output


# ---------------------------------------------------------------------------
# main --help includes utils
# ---------------------------------------------------------------------------

def test_main_help_includes_utils(runner):
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "utils" in result.output
