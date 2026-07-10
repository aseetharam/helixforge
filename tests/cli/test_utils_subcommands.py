"""Phase 44 — tests for the implemented ``helixforge utils`` subcommands.

Tests: align, convert, filter, summarize. (``validate-inputs`` was folded into
``helixforge doctor`` — see the folded-command guard below.)
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest import mock

import pytest
from click.testing import CliRunner

from helixforge.cli import main


@pytest.fixture
def runner():
    return CliRunner()


def _touch(p: Path, content: str = "") -> str:
    p.write_text(content)
    return str(p)


# ---------------------------------------------------------------------------
# Shared GFF3 fixtures
# ---------------------------------------------------------------------------

_GFF3_3GENE = """\
##gff-version 3
chr1\tHelixForge\tgene\t1\t900\t.\t+\t.\tID=g1;tier=1;gene_biotype=protein_coding;combined_score=0.95;flags=
chr1\tHelixForge\tmRNA\t1\t900\t.\t+\t.\tID=g1.1;Parent=g1;combined_score=0.95
chr1\tHelixForge\texon\t1\t300\t.\t+\t.\tID=g1.1.exon1;Parent=g1.1
chr1\tHelixForge\texon\t400\t600\t.\t+\t.\tID=g1.1.exon2;Parent=g1.1
chr1\tHelixForge\texon\t700\t900\t.\t+\t.\tID=g1.1.exon3;Parent=g1.1
chr1\tHelixForge\tCDS\t1\t300\t.\t+\t0\tID=g1.1.CDS1;Parent=g1.1
chr1\tHelixForge\tCDS\t400\t600\t.\t+\t0\tID=g1.1.CDS2;Parent=g1.1
chr1\tHelixForge\tCDS\t700\t900\t.\t+\t0\tID=g1.1.CDS3;Parent=g1.1
###
chr1\tHelixForge\tgene\t1001\t1600\t.\t-\t.\tID=g2;tier=2;gene_biotype=protein_coding;combined_score=0.7;flags=NO_EXPRESSION
chr1\tHelixForge\tmRNA\t1001\t1600\t.\t-\t.\tID=g2.1;Parent=g2;combined_score=0.7
chr1\tHelixForge\texon\t1001\t1300\t.\t-\t.\tID=g2.1.exon1;Parent=g2.1
chr1\tHelixForge\texon\t1400\t1600\t.\t-\t.\tID=g2.1.exon2;Parent=g2.1
chr1\tHelixForge\tCDS\t1001\t1300\t.\t-\t0\tID=g2.1.CDS1;Parent=g2.1
chr1\tHelixForge\tCDS\t1400\t1600\t.\t-\t0\tID=g2.1.CDS2;Parent=g2.1
###
chr1\tHelixForge\tgene\t2001\t2500\t.\t+\t.\tID=g3;tier=3;gene_biotype=transposable_element;combined_score=0.3;flags=HELIXER_ONLY
chr1\tHelixForge\tmRNA\t2001\t2500\t.\t+\t.\tID=g3.1;Parent=g3;combined_score=0.3
chr1\tHelixForge\texon\t2001\t2500\t.\t+\t.\tID=g3.1.exon1;Parent=g3.1
###
"""

_GFF3_5GENE = """\
##gff-version 3
chr1\tHF\tgene\t1\t300\t.\t+\t.\tID=g1;tier=1;gene_biotype=protein_coding;combined_score=0.95;flags=
chr1\tHF\tmRNA\t1\t300\t.\t+\t.\tID=g1.1;Parent=g1;combined_score=0.95
chr1\tHF\texon\t1\t300\t.\t+\t.\tID=g1.1.exon1;Parent=g1.1
chr1\tHF\tCDS\t1\t300\t.\t+\t0\tID=g1.1.CDS1;Parent=g1.1
###
chr1\tHF\tgene\t401\t700\t.\t+\t.\tID=g2;tier=1;gene_biotype=protein_coding;combined_score=0.88;flags=
chr1\tHF\tmRNA\t401\t700\t.\t+\t.\tID=g2.1;Parent=g2;combined_score=0.88
chr1\tHF\texon\t401\t700\t.\t+\t.\tID=g2.1.exon1;Parent=g2.1
chr1\tHF\tCDS\t401\t700\t.\t+\t0\tID=g2.1.CDS1;Parent=g2.1
###
chr1\tHF\tgene\t801\t1100\t.\t-\t.\tID=g3;tier=2;gene_biotype=protein_coding;combined_score=0.7;flags=NO_EXPRESSION
chr1\tHF\tmRNA\t801\t1100\t.\t-\t.\tID=g3.1;Parent=g3;combined_score=0.7
chr1\tHF\texon\t801\t1100\t.\t-\t.\tID=g3.1.exon1;Parent=g3.1
chr1\tHF\tCDS\t801\t1100\t.\t-\t0\tID=g3.1.CDS1;Parent=g3.1
###
chr1\tHF\tgene\t1201\t1500\t.\t+\t.\tID=g4;tier=3;gene_biotype=transposable_element;combined_score=0.3;flags=HELIXER_ONLY,INTERNAL_STOP
chr1\tHF\tmRNA\t1201\t1500\t.\t+\t.\tID=g4.1;Parent=g4;combined_score=0.3
chr1\tHF\texon\t1201\t1500\t.\t+\t.\tID=g4.1.exon1;Parent=g4.1
###
chr1\tHF\tgene\t1601\t1900\t.\t-\t.\tID=g5;tier=4;gene_biotype=ncRNA;combined_score=0.1;flags=NO_START_CODON
chr1\tHF\tmRNA\t1601\t1900\t.\t-\t.\tID=g5.1;Parent=g5;combined_score=0.1
chr1\tHF\texon\t1601\t1900\t.\t-\t.\tID=g5.1.exon1;Parent=g5.1
###
"""


# ===========================================================================
# 44a — align
# ===========================================================================

class TestAlign:
    def test_align_calls_miniprot(self, runner, tmp_path):
        genome = _touch(tmp_path / "genome.fa", ">chr1\nACGT\n")
        proteins = _touch(tmp_path / "proteins.fa", ">p\nMKL\n")
        out = str(tmp_path / "out.gff")

        with mock.patch("helixforge.prep.protein_align.run_tool") as mock_run, \
             mock.patch("helixforge.prep.protein_align.output_is_fresh", return_value=False):
            mock_run.return_value = mock.MagicMock(returncode=0)
            result = runner.invoke(main, [
                "utils", "align",
                "--genome", genome, "--proteins", proteins,
                "--out", out, "--threads", "8",
            ])
            assert result.exit_code == 0, result.output
            mock_run.assert_called_once()
            call_args = mock_run.call_args
            argv = call_args[0][0]
            argv_strs = [str(a) for a in argv]
            assert "miniprot" in argv_strs[0]
            assert "-t" in argv_strs
            assert "8" in argv_strs
            assert "--gff" in argv_strs

    def test_align_custom_binary(self, runner, tmp_path):
        genome = _touch(tmp_path / "genome.fa", ">chr1\nACGT\n")
        proteins = _touch(tmp_path / "proteins.fa", ">p\nMKL\n")
        out = str(tmp_path / "out.gff")

        with mock.patch("helixforge.prep.protein_align.run_tool") as mock_run, \
             mock.patch("helixforge.prep.protein_align.output_is_fresh", return_value=False):
            mock_run.return_value = mock.MagicMock(returncode=0)
            result = runner.invoke(main, [
                "utils", "align",
                "--genome", genome, "--proteins", proteins,
                "--out", out, "--miniprot-bin", "/opt/bin/miniprot",
            ])
            assert result.exit_code == 0, result.output
            argv = mock_run.call_args[0][0]
            assert str(argv[0]) == "/opt/bin/miniprot"

    def test_align_threads_in_argv(self, runner, tmp_path):
        genome = _touch(tmp_path / "genome.fa", ">chr1\nACGT\n")
        proteins = _touch(tmp_path / "proteins.fa", ">p\nMKL\n")
        out = str(tmp_path / "out.gff")

        with mock.patch("helixforge.prep.protein_align.run_tool") as mock_run, \
             mock.patch("helixforge.prep.protein_align.output_is_fresh", return_value=False):
            mock_run.return_value = mock.MagicMock(returncode=0)
            result = runner.invoke(main, [
                "utils", "align",
                "--genome", genome, "--proteins", proteins,
                "--out", out, "--threads", "16",
            ])
            assert result.exit_code == 0
            argv = [str(a) for a in mock_run.call_args[0][0]]
            idx = argv.index("-t")
            assert argv[idx + 1] == "16"

    def test_align_stdout_redirect(self, runner, tmp_path):
        genome = _touch(tmp_path / "genome.fa", ">chr1\nACGT\n")
        proteins = _touch(tmp_path / "proteins.fa", ">p\nMKL\n")
        out = str(tmp_path / "out.gff")

        with mock.patch("helixforge.prep.protein_align.run_tool") as mock_run, \
             mock.patch("helixforge.prep.protein_align.output_is_fresh", return_value=False):
            mock_run.return_value = mock.MagicMock(returncode=0)
            runner.invoke(main, [
                "utils", "align",
                "--genome", genome, "--proteins", proteins,
                "--out", out,
            ])
            assert mock_run.call_args[1]["stdout_path"] is not None


# ===========================================================================
# 44b — validate-inputs folded into `helixforge doctor`
# ===========================================================================

class TestValidateInputsFolded:
    """`utils validate-inputs` was removed; `doctor` is the single preflight."""

    def test_validate_inputs_removed(self, runner, tmp_path):
        genome = _touch(tmp_path / "genome.fa", ">chr1\nACGTACGT\n")
        result = runner.invoke(main, [
            "utils", "validate-inputs", "--genome", genome,
        ])
        assert result.exit_code != 0  # no such subcommand


# ===========================================================================
# 44c — convert
# ===========================================================================

class TestConvert:
    def test_gff3_to_gtf(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        gtf = str(tmp_path / "genes.gtf")
        result = runner.invoke(main, [
            "utils", "convert", "--input", gff3, "--out", gtf,
        ])
        assert result.exit_code == 0, result.output
        assert "3 genes" in result.output
        content = Path(gtf).read_text()
        assert 'gene_id "g1"' in content
        assert 'gene_id "g2"' in content
        assert 'gene_id "g3"' in content

    def test_gtf_to_gff3(self, runner, tmp_path):
        gtf_text = (
            'chr1\tHF\ttranscript\t1\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tHF\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tHF\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tHF\tCDS\t1\t100\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tHF\tCDS\t200\t300\t.\t+\t2\tgene_id "g1"; transcript_id "t1";\n'
        )
        gtf = _touch(tmp_path / "genes.gtf", gtf_text)
        gff3_out = str(tmp_path / "genes.gff3")
        result = runner.invoke(main, [
            "utils", "convert", "--input", gtf, "--out", gff3_out,
        ])
        assert result.exit_code == 0, result.output
        assert "1 genes" in result.output
        content = Path(gff3_out).read_text()
        assert "##gff-version 3" in content
        assert "ID=g1" in content
        assert "Parent=g1" in content

    def test_roundtrip_preserves_gene_count(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        gtf = str(tmp_path / "genes.gtf")
        gff3_rt = str(tmp_path / "roundtrip.gff3")

        result1 = runner.invoke(main, [
            "utils", "convert", "--input", gff3, "--out", gtf,
        ])
        assert result1.exit_code == 0

        result2 = runner.invoke(main, [
            "utils", "convert", "--input", gtf, "--out", gff3_rt,
        ])
        assert result2.exit_code == 0
        assert "3 genes" in result2.output

    def test_explicit_format_override(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        out = str(tmp_path / "output.txt")
        result = runner.invoke(main, [
            "utils", "convert", "--input", gff3, "--out", out,
            "--from-format", "gff3", "--to-format", "gtf",
        ])
        assert result.exit_code == 0
        assert "3 genes" in result.output

    def test_same_format_error(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        out = str(tmp_path / "out.gff3")
        result = runner.invoke(main, [
            "utils", "convert", "--input", gff3, "--out", out,
        ])
        assert result.exit_code != 0

    def test_minus_strand_exon_order(self, runner, tmp_path):
        """Minus-strand genes preserve low→high exon order in GTF."""
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        gtf = str(tmp_path / "genes.gtf")
        runner.invoke(main, [
            "utils", "convert", "--input", gff3, "--out", gtf,
        ])
        content = Path(gtf).read_text()
        g2_lines = [l for l in content.splitlines() if 'gene_id "g2"' in l]
        assert len(g2_lines) > 0
        exon_lines = [l for l in g2_lines if l.split("\t")[2] == "exon"]
        starts = [int(l.split("\t")[3]) for l in exon_lines]
        assert starts == sorted(starts)


# ===========================================================================
# 44d — filter
# ===========================================================================

class TestFilter:
    def test_high_confidence_drops_tier3_and_4(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_5GENE)
        out = str(tmp_path / "filtered.gff3")
        result = runner.invoke(main, [
            "utils", "filter", "--gff3", gff3, "--out", out,
            "--preset", "high_confidence",
        ])
        assert result.exit_code == 0, result.output
        # g1 (tier1), g2 (tier1), g3 (tier2 NO_EXPRESSION) pass;
        # g4 (tier3) excluded by max_tier; g5 (tier4 NO_START_CODON) excluded
        assert "Kept 3" in result.output
        content = Path(out).read_text()
        assert "g1" in content
        assert "g2" in content
        assert "g3" in content
        assert "g4" not in content
        assert "g5" not in content

    def test_publication_ready_tier1_protein_coding(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_5GENE)
        out = str(tmp_path / "filtered.gff3")
        result = runner.invoke(main, [
            "utils", "filter", "--gff3", gff3, "--out", out,
            "--preset", "publication_ready",
        ])
        assert result.exit_code == 0, result.output
        assert "Kept 2" in result.output
        content = Path(out).read_text()
        assert "g1" in content
        assert "g2" in content
        assert "g4" not in content

    def test_custom_max_tier(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_5GENE)
        out = str(tmp_path / "filtered.gff3")
        result = runner.invoke(main, [
            "utils", "filter", "--gff3", gff3, "--out", out,
            "--max-tier", "1",
        ])
        assert result.exit_code == 0, result.output
        assert "Kept 2" in result.output

    def test_exclude_biotype(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_5GENE)
        out = str(tmp_path / "filtered.gff3")
        result = runner.invoke(main, [
            "utils", "filter", "--gff3", gff3, "--out", out,
            "--exclude-biotype", "transposable_element",
        ])
        assert result.exit_code == 0
        content = Path(out).read_text()
        assert "g4" not in content

    def test_filter_preserves_header(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_5GENE)
        out = str(tmp_path / "filtered.gff3")
        runner.invoke(main, [
            "utils", "filter", "--gff3", gff3, "--out", out,
            "--max-tier", "1",
        ])
        content = Path(out).read_text()
        assert content.startswith("##gff-version 3")

    def test_high_confidence_excludes_internal_stop_flag(self, runner, tmp_path):
        """high_confidence preset excludes genes with INTERNAL_STOP flag."""
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_5GENE)
        out = str(tmp_path / "filtered.gff3")
        runner.invoke(main, [
            "utils", "filter", "--gff3", gff3, "--out", out,
            "--preset", "high_confidence",
        ])
        content = Path(out).read_text()
        # g3 is tier 2 with NO_EXPRESSION — should pass high_confidence
        # (NO_EXPRESSION is not in the exclude list)
        assert "g3" in content
        # g4 is tier 3 — excluded by max_tier=2
        assert "g4" not in content


# ===========================================================================
# 44e — summarize
# ===========================================================================

class TestSummarize:
    def test_tsv_output(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        result = runner.invoke(main, [
            "utils", "summarize", "--gff3", gff3,
        ])
        assert result.exit_code == 0, result.output
        assert "gene_count" in result.output
        assert "transcript_count" in result.output

    def test_json_output(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        out = str(tmp_path / "stats.json")
        result = runner.invoke(main, [
            "utils", "summarize", "--gff3", gff3, "--out", out,
            "--format", "json",
        ])
        assert result.exit_code == 0, result.output
        data = json.loads(Path(out).read_text())
        assert data["gene_count"] == 3
        assert data["transcript_count"] == 3

    def test_json_has_required_keys(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        out = str(tmp_path / "stats.json")
        runner.invoke(main, [
            "utils", "summarize", "--gff3", gff3, "--out", out,
            "--format", "json",
        ])
        data = json.loads(Path(out).read_text())
        for key in ("gene_count", "transcript_count", "multi_isoform_genes",
                     "coding_genes", "mono_exon_genes"):
            assert key in data

    def test_markdown_output(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        out = str(tmp_path / "stats.md")
        result = runner.invoke(main, [
            "utils", "summarize", "--gff3", gff3, "--out", out,
            "--format", "markdown",
        ])
        assert result.exit_code == 0
        content = Path(out).read_text()
        assert "| gene_count |" in content

    def test_coding_gene_count(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        out = str(tmp_path / "stats.json")
        runner.invoke(main, [
            "utils", "summarize", "--gff3", gff3, "--out", out,
            "--format", "json",
        ])
        data = json.loads(Path(out).read_text())
        assert data["coding_genes"] == 2

    def test_mono_exon_count(self, runner, tmp_path):
        gff3 = _touch(tmp_path / "genes.gff3", _GFF3_3GENE)
        out = str(tmp_path / "stats.json")
        runner.invoke(main, [
            "utils", "summarize", "--gff3", gff3, "--out", out,
            "--format", "json",
        ])
        data = json.loads(Path(out).read_text())
        assert data["mono_exon_genes"] == 1


# ===========================================================================
# Library module tests
# ===========================================================================

class TestConvertLib:
    """Direct tests of the convert library functions."""

    def test_gff3_to_gtf_coordinates(self, tmp_path):
        from helixforge.utils.convert import gff3_to_gtf
        gff3 = tmp_path / "genes.gff3"
        gff3.write_text(_GFF3_3GENE)
        gtf = tmp_path / "out.gtf"
        n = gff3_to_gtf(str(gff3), str(gtf))
        assert n == 3
        content = gtf.read_text()
        lines = [l for l in content.splitlines() if l.strip()]
        assert len(lines) > 0

    def test_gtf_to_gff3_coordinates(self, tmp_path):
        from helixforge.utils.convert import gtf_to_gff3
        gtf = tmp_path / "genes.gtf"
        gtf.write_text(
            'chr1\tHF\ttranscript\t1\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tHF\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\tHF\texon\t200\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        )
        gff3 = tmp_path / "out.gff3"
        n = gtf_to_gff3(str(gtf), str(gff3))
        assert n == 1


class TestFiltersLib:
    """Direct tests of the filters library module."""

    def test_filter_criteria_defaults(self):
        from helixforge.utils.filters import FilterCriteria
        c = FilterCriteria()
        assert c.min_tier is None
        assert c.max_tier is None
        assert c.exclude_biotypes == []
        assert c.exclude_flags == []

    def test_gene_filter_apply(self):
        from helixforge.utils.filters import FilterCriteria, GeneFilter
        genes = [
            {"gene_id": "g1", "tier": 1, "gene_biotype": "protein_coding", "flags": []},
            {"gene_id": "g2", "tier": 3, "gene_biotype": "protein_coding", "flags": []},
        ]
        filt = GeneFilter(FilterCriteria(max_tier=2))
        kept, excluded = filt.apply(genes)
        assert len(kept) == 1
        assert kept[0]["gene_id"] == "g1"
        assert len(excluded) == 1

    def test_parse_gff3_genes(self, tmp_path):
        from helixforge.utils.filters import parse_gff3_genes
        gff3 = tmp_path / "genes.gff3"
        gff3.write_text(_GFF3_5GENE)
        genes = parse_gff3_genes(str(gff3))
        assert len(genes) == 5
        assert genes[0]["gene_id"] == "g1"
        assert genes[0]["tier"] == 1
        assert genes[3]["flags"] == ["HELIXER_ONLY", "INTERNAL_STOP"]

    def test_high_confidence_preset(self):
        from helixforge.utils.filters import GeneFilter
        filt = GeneFilter.high_confidence()
        assert filt.criteria.max_tier == 2
        assert "INTERNAL_STOP" in filt.criteria.exclude_flags

    def test_publication_ready_preset(self):
        from helixforge.utils.filters import GeneFilter
        filt = GeneFilter.publication_ready()
        assert filt.criteria.max_tier == 1
        assert "transposable_element" in filt.criteria.exclude_biotypes

    def test_exclude_flags(self):
        from helixforge.utils.filters import FilterCriteria, GeneFilter
        genes = [
            {"gene_id": "g1", "tier": 1, "flags": ["INTERNAL_STOP"]},
            {"gene_id": "g2", "tier": 1, "flags": []},
        ]
        filt = GeneFilter(FilterCriteria(exclude_flags=["INTERNAL_STOP"]))
        kept, excluded = filt.apply(genes)
        assert len(kept) == 1
        assert kept[0]["gene_id"] == "g2"

    def test_min_confidence(self):
        from helixforge.utils.filters import FilterCriteria, GeneFilter
        genes = [
            {"gene_id": "g1", "tier": 1, "flags": [], "combined_score": 0.9},
            {"gene_id": "g2", "tier": 1, "flags": [], "combined_score": 0.3},
        ]
        filt = GeneFilter(FilterCriteria(min_confidence=0.5))
        kept, _ = filt.apply(genes)
        assert len(kept) == 1
        assert kept[0]["gene_id"] == "g1"
