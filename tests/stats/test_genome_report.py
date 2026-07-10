"""Tests for stats.genome_report — standalone GFF3 QC report generation.

Concrete literal values only; synthetic fixtures (no real data files).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from helixforge.stats.genome_report import generate_qc_report


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_GFF3_CONTENT = """\
##gff-version 3
chr1\tHelixForge\tgene\t1001\t5000\t.\t+\t.\tID=HFG_00001;gene_biotype=protein_coding;tier=1;origin=mikado_1to1
chr1\tHelixForge\tmRNA\t1001\t5000\t.\t+\t.\tID=HFG_00001.1;Parent=HFG_00001;flags=ALL_JUNCTIONS_SUPPORTED;as_events=ES@1501-2000
chr1\tHelixForge\texon\t1001\t1500\t.\t+\t.\tID=HFG_00001.1.exon1;Parent=HFG_00001.1
chr1\tHelixForge\texon\t2001\t3000\t.\t+\t.\tID=HFG_00001.1.exon2;Parent=HFG_00001.1
chr1\tHelixForge\texon\t3501\t5000\t.\t+\t.\tID=HFG_00001.1.exon3;Parent=HFG_00001.1
chr1\tHelixForge\tCDS\t1101\t1500\t.\t+\t0\tID=HFG_00001.1.CDS1;Parent=HFG_00001.1
chr1\tHelixForge\tCDS\t2001\t3000\t.\t+\t0\tID=HFG_00001.1.CDS2;Parent=HFG_00001.1
chr1\tHelixForge\tCDS\t3501\t4200\t.\t+\t0\tID=HFG_00001.1.CDS3;Parent=HFG_00001.1
###
chr1\tHelixForge\tgene\t6001\t8000\t.\t-\t.\tID=HFG_00002;gene_biotype=protein_coding;tier=1;origin=mikado_1to1
chr1\tHelixForge\tmRNA\t6001\t8000\t.\t-\t.\tID=HFG_00002.1;Parent=HFG_00002;flags=NO_STOP
chr1\tHelixForge\texon\t6001\t7000\t.\t-\t.\tID=HFG_00002.1.exon1;Parent=HFG_00002.1
chr1\tHelixForge\texon\t7501\t8000\t.\t-\t.\tID=HFG_00002.1.exon2;Parent=HFG_00002.1
chr1\tHelixForge\tCDS\t6200\t7000\t.\t-\t0\tID=HFG_00002.1.CDS1;Parent=HFG_00002.1
chr1\tHelixForge\tCDS\t7501\t7800\t.\t-\t2\tID=HFG_00002.1.CDS2;Parent=HFG_00002.1
###
chr1\tHelixForge\tgene\t9001\t12000\t.\t+\t.\tID=HFG_00003;gene_biotype=protein_coding;tier=2;origin=split
chr1\tHelixForge\tmRNA\t9001\t12000\t.\t+\t.\tID=HFG_00003.1;Parent=HFG_00003;flags=PARTIAL_JUNCTION_SUPPORT
chr1\tHelixForge\texon\t9001\t10000\t.\t+\t.\tID=HFG_00003.1.exon1;Parent=HFG_00003.1
chr1\tHelixForge\texon\t10501\t12000\t.\t+\t.\tID=HFG_00003.1.exon2;Parent=HFG_00003.1
chr1\tHelixForge\tCDS\t9201\t10000\t.\t+\t0\tID=HFG_00003.1.CDS1;Parent=HFG_00003.1
chr1\tHelixForge\tCDS\t10501\t11500\t.\t+\t1\tID=HFG_00003.1.CDS2;Parent=HFG_00003.1
chr1\tHelixForge\tmRNA\t9001\t11500\t.\t+\t.\tID=HFG_00003.2;Parent=HFG_00003;flags=PARTIAL_JUNCTION_SUPPORT;as_events=A3@10001-10500
chr1\tHelixForge\texon\t9001\t9800\t.\t+\t.\tID=HFG_00003.2.exon1;Parent=HFG_00003.2
chr1\tHelixForge\texon\t10501\t11500\t.\t+\t.\tID=HFG_00003.2.exon2;Parent=HFG_00003.2
chr1\tHelixForge\tCDS\t9201\t9800\t.\t+\t0\tID=HFG_00003.2.CDS1;Parent=HFG_00003.2
chr1\tHelixForge\tCDS\t10501\t11200\t.\t+\t1\tID=HFG_00003.2.CDS2;Parent=HFG_00003.2
###
chr1\tHelixForge\tgene\t13001\t14000\t.\t+\t.\tID=HFG_00004;gene_biotype=lncRNA;tier=2;origin=mikado_1to1
chr1\tHelixForge\tmRNA\t13001\t14000\t.\t+\t.\tID=HFG_00004.1;Parent=HFG_00004;flags=NO_HOMOL
chr1\tHelixForge\texon\t13001\t14000\t.\t+\t.\tID=HFG_00004.1.exon1;Parent=HFG_00004.1
###
chr1\tHelixForge\tgene\t15001\t16500\t.\t-\t.\tID=HFG_00005;gene_biotype=protein_coding;tier=1;origin=merge
chr1\tHelixForge\tmRNA\t15001\t16500\t.\t-\t.\tID=HFG_00005.1;Parent=HFG_00005;flags=LOCUS_MERGE
chr1\tHelixForge\texon\t15001\t16000\t.\t-\t.\tID=HFG_00005.1.exon1;Parent=HFG_00005.1
chr1\tHelixForge\texon\t16201\t16500\t.\t-\t.\tID=HFG_00005.1.exon2;Parent=HFG_00005.1
chr1\tHelixForge\tCDS\t15001\t16000\t.\t-\t0\tID=HFG_00005.1.CDS1;Parent=HFG_00005.1
chr1\tHelixForge\tCDS\t16201\t16500\t.\t-\t0\tID=HFG_00005.1.CDS2;Parent=HFG_00005.1
###
chr2\tHelixForge\tgene\t1001\t3000\t.\t+\t.\tID=HFG_00006;gene_biotype=protein_coding;tier=3;origin=helixer_backstop
chr2\tHelixForge\tmRNA\t1001\t3000\t.\t+\t.\tID=HFG_00006.1;Parent=HFG_00006;flags=HELIXER_ONLY,BACKSTOP_RESCUED
chr2\tHelixForge\texon\t1001\t2000\t.\t+\t.\tID=HFG_00006.1.exon1;Parent=HFG_00006.1
chr2\tHelixForge\texon\t2501\t3000\t.\t+\t.\tID=HFG_00006.1.exon2;Parent=HFG_00006.1
chr2\tHelixForge\tCDS\t1201\t2000\t.\t+\t0\tID=HFG_00006.1.CDS1;Parent=HFG_00006.1
chr2\tHelixForge\tCDS\t2501\t2900\t.\t+\t2\tID=HFG_00006.1.CDS2;Parent=HFG_00006.1
###
chr2\tHelixForge\tgene\t4001\t5000\t.\t+\t.\tID=HFG_00007;gene_biotype=pseudogene;tier=4;origin=novel
chr2\tHelixForge\tmRNA\t4001\t5000\t.\t+\t.\tID=HFG_00007.1;Parent=HFG_00007;flags=PSEUDOGENE_CANDIDATE,NO_HOMOL
chr2\tHelixForge\texon\t4001\t5000\t.\t+\t.\tID=HFG_00007.1.exon1;Parent=HFG_00007.1
###
chr2\tHelixForge\tgene\t6001\t9000\t.\t-\t.\tID=HFG_00008;gene_biotype=protein_coding;tier=1;origin=mikado_1to1
chr2\tHelixForge\tmRNA\t6001\t9000\t.\t-\t.\tID=HFG_00008.1;Parent=HFG_00008;flags=ALL_JUNCTIONS_SUPPORTED;as_events=IR@7001-7500
chr2\tHelixForge\texon\t6001\t7000\t.\t-\t.\tID=HFG_00008.1.exon1;Parent=HFG_00008.1
chr2\tHelixForge\texon\t7501\t8000\t.\t-\t.\tID=HFG_00008.1.exon2;Parent=HFG_00008.1
chr2\tHelixForge\texon\t8501\t9000\t.\t-\t.\tID=HFG_00008.1.exon3;Parent=HFG_00008.1
chr2\tHelixForge\tCDS\t6001\t7000\t.\t-\t0\tID=HFG_00008.1.CDS1;Parent=HFG_00008.1
chr2\tHelixForge\tCDS\t7501\t8000\t.\t-\t2\tID=HFG_00008.1.CDS2;Parent=HFG_00008.1
chr2\tHelixForge\tCDS\t8501\t9000\t.\t-\t1\tID=HFG_00008.1.CDS3;Parent=HFG_00008.1
###
chr2\tHelixForge\tgene\t10001\t11000\t.\t+\t.\tID=HFG_00009;gene_biotype=lncRNA;tier=3;origin=mikado_1to1
chr2\tHelixForge\tmRNA\t10001\t11000\t.\t+\t.\tID=HFG_00009.1;Parent=HFG_00009;flags=NO_EXPRESSION
chr2\tHelixForge\texon\t10001\t10500\t.\t+\t.\tID=HFG_00009.1.exon1;Parent=HFG_00009.1
chr2\tHelixForge\texon\t10701\t11000\t.\t+\t.\tID=HFG_00009.1.exon2;Parent=HFG_00009.1
###
chr2\tHelixForge\tgene\t12001\t13500\t.\t-\t.\tID=HFG_00010;gene_biotype=protein_coding;tier=2;origin=mikado_1to1
chr2\tHelixForge\tmRNA\t12001\t13500\t.\t-\t.\tID=HFG_00010.1;Parent=HFG_00010;flags=NO_START
chr2\tHelixForge\texon\t12001\t12500\t.\t-\t.\tID=HFG_00010.1.exon1;Parent=HFG_00010.1
chr2\tHelixForge\texon\t13001\t13500\t.\t-\t.\tID=HFG_00010.1.exon2;Parent=HFG_00010.1
chr2\tHelixForge\tCDS\t12001\t12500\t.\t-\t0\tID=HFG_00010.1.CDS1;Parent=HFG_00010.1
chr2\tHelixForge\tCDS\t13001\t13200\t.\t-\t2\tID=HFG_00010.1.CDS2;Parent=HFG_00010.1
###
"""


@pytest.fixture
def gff3_path(tmp_path: Path) -> Path:
    p = tmp_path / "test_annotation.gff3"
    p.write_text(_GFF3_CONTENT)
    return p


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestGenerateQcReport:
    """Test the three output formats and edge cases."""

    def test_json_output(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        assert out.exists()
        parsed = json.loads(out.read_text())

        assert parsed["total_genes"] == 10
        assert parsed["total_transcripts"] == 11
        assert parsed["multi_isoform_genes"] == 1
        assert parsed["coding_genes"] == 7

        assert "tier" in parsed
        assert "biotype" in parsed
        assert "flags" in parsed
        assert "cds_length" in parsed
        assert "exon_count" in parsed
        assert "intron_length" in parsed

        assert stats["total_genes"] == 10

    def test_html_output(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.html"
        stats = generate_qc_report(gff3_path, out, fmt="html")

        assert out.exists()
        html = out.read_text()

        assert "<html" in html
        assert "Tier Distribution" in html
        assert "QC Flags" in html
        assert "Summary Statistics" in html
        assert "Biotype Distribution" in html
        assert "test_annotation" in html
        assert stats["total_genes"] == 10

    def test_tsv_output(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.tsv"
        stats = generate_qc_report(gff3_path, out, fmt="tsv")

        assert out.exists()
        lines = out.read_text().strip().split("\n")
        kv = {}
        for line in lines:
            parts = line.split("\t", 1)
            if len(parts) == 2:
                kv[parts[0]] = parts[1]

        assert kv["total_genes"] == "10"
        assert kv["total_transcripts"] == "11"
        assert kv["coding_genes"] == "7"
        assert stats["total_genes"] == 10

    def test_no_genome_no_h5(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(
            gff3_path, out, genome_path=None, helixer_h5_path=None, fmt="json"
        )
        assert "genome" not in stats
        assert "confidence" not in stats
        assert stats["total_genes"] == 10

    def test_tier_counts(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        assert stats["tier"]["1"] == 4
        assert stats["tier"]["2"] == 3
        assert stats["tier"]["3"] == 2
        assert stats["tier"]["4"] == 1

    def test_biotype_counts(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        assert stats["biotype"]["protein_coding"] == 7
        assert stats["biotype"]["lncRNA"] == 2
        assert stats["biotype"]["pseudogene"] == 1

    def test_origin_counts(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        assert stats["origin"]["mikado_1to1"] == 6
        assert stats["origin"]["split"] == 1
        assert stats["origin"]["merge"] == 1
        assert stats["origin"]["helixer_backstop"] == 1
        assert stats["origin"]["novel"] == 1

    def test_flag_counts(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        flags = stats["flags"]
        assert flags["ALL_JUNCTIONS_SUPPORTED"] == 2
        assert flags["PARTIAL_JUNCTION_SUPPORT"] == 2
        assert flags["NO_HOMOL"] == 2
        assert flags["HELIXER_ONLY"] == 1
        assert flags["BACKSTOP_RESCUED"] == 1
        assert flags["NO_STOP"] == 1
        assert flags["LOCUS_MERGE"] == 1
        assert flags["PSEUDOGENE_CANDIDATE"] == 1
        assert flags["NO_EXPRESSION"] == 1
        assert flags["NO_START"] == 1

    def test_as_event_counts(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        assert stats["as_events"]["ES"] == 1
        assert stats["as_events"]["A3"] == 1
        assert stats["as_events"]["IR"] == 1

    def test_cds_length_distribution(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        cds = stats["cds_length"]
        assert cds["count"] == 8
        assert cds["min"] > 0
        assert cds["max"] > cds["min"]

    def test_multi_isoform_gene(self, gff3_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "report.json"
        stats = generate_qc_report(gff3_path, out, fmt="json")

        assert stats["multi_isoform_genes"] == 1

    def test_returns_same_stats_regardless_of_format(
        self, gff3_path: Path, tmp_path: Path
    ) -> None:
        stats_html = generate_qc_report(gff3_path, tmp_path / "r.html", fmt="html")
        stats_json = generate_qc_report(gff3_path, tmp_path / "r.json", fmt="json")
        stats_tsv = generate_qc_report(gff3_path, tmp_path / "r.tsv", fmt="tsv")

        assert stats_html["total_genes"] == stats_json["total_genes"] == stats_tsv["total_genes"]
        assert stats_html["tier"] == stats_json["tier"] == stats_tsv["tier"]
