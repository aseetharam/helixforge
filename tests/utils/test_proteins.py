"""Tests for ``helixforge.utils.proteins``."""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

from helixforge.utils.proteins import extract_proteins


# ---------------------------------------------------------------------------
# Synthetic genome + GFF3 fixtures
# ---------------------------------------------------------------------------
#
# Scaffold layout (0-based):
#
# scaffold_1 (500 bp):
#   gene_a (+strand, 3-exon, valid CDS → 36 nt = 12 aa protein "MFLIVGPRSTKD")
#     exon1: 100–112  (12 bp)  CDS: 100–112
#     exon2: 150–162  (12 bp)  CDS: 150–162
#     exon3: 200–212  (12 bp)  CDS: 200–212
#     transcript_a2: same gene, shorter CDS (exons 1+2 only → 24 nt = 8 aa)
#   gene_d (+strand, 1-exon, CDS = 9 nt → 3 aa protein, below min_length=30)
#     exon: 300–309 (9 bp)  CDS: 300–309
#
# scaffold_2 (500 bp):
#   gene_b (-strand, 2-exon, valid CDS → 36 nt = 12 aa)
#     exon1: 100–118  (18 bp)  CDS: 100–118
#     exon2: 200–218  (18 bp)  CDS: 200–218
#   gene_c (+strand, 1-exon, CDS with internal stop → still translated faithfully)
#     exon: 300–318  (18 bp)  CDS: 300–318
#
# For gene_a (+strand): CDS nt = seq[100:112] + seq[150:162] + seq[200:212]
#   We engineer the sequence so this gives ATG + 11 coding codons.
#
# For gene_b (-strand): CDS nt = RC(seq[200:218] + seq[100:118])
#   reading order is high-to-low, reversed then RC'd.
#   We engineer the sequence so this gives ATG + 11 coding codons.

# We need to carefully design the genome sequences.
# Standard codon table:
#   ATG=M, TTT=F, CTG=L, ATC=I, GTG=V, GGC=G, CCC=P, CGT=R, TCC=S, ACC=T, AAG=K, GAT=D

# gene_a (+strand): 3 exons of 12bp each, concatenated = 36 nt
_GENE_A_CDS_SEQ = "ATGTTTCTGATCGTGGGCCCCCGTTCCACCAAGGAT"  # 36 nt, 12 aa: MFLIVGPRSTKD
# protein: M F L I V G P R S T K D

# gene_b (-strand): 2 exons of 18bp each. Reading order: exon2(200-218) then exon1(100-118),
# reversed and RC'd. So the plus-strand sequence for exon1(100-118) and exon2(200-218) on
# scaffold_2 must be chosen so that RC(exon2_seq + exon1_seq) = desired CDS.
# Desired CDS: ATGCTTAACGCCGGTTACTCGATCAATCCGGCCTAG  (36 nt)
# protein: M L N A G Y S I N P A *  -- wait, last 3 is stop codon but CDS excludes stop.
# Let's do 36 nt CDS (no stop included) = 12 aa protein.
_GENE_B_CDS_DESIRED = "ATGCTTAACGCCGGTTACTCGATCAATCCGGCCGAT"  # 36 nt, 12 aa: MLNAGYSINFAD
# protein: MLNAGYSINPAD

# For minus strand: reading order is exon2 (higher coords) → exon1 (lower coords)
# Plus-strand stored: exon1_seq || exon2_seq
# Reading order for -strand: RC(exon2_seq + exon1_seq) = _GENE_B_CDS_DESIRED
# So exon2_seq + exon1_seq = RC(_GENE_B_CDS_DESIRED)
# RC("ATGCTTAACGCCGGTTACTCGATCAATCCGGCCGAT") = "ATCGGCCGGATTGATCGAGTAACCGGCGTTAAGCAT"
# exon2 is 18bp (positions 200-218): first 18 chars of the RC
# exon1 is 18bp (positions 100-118): last 18 chars of the RC
_GENE_B_RC = "ATCGGCCGGATTGATCGAGTAACCGGCGTTAAGCAT"  # 36 chars
_GENE_B_EXON2_PLUS = _GENE_B_RC[:18]  # goes at scaffold_2[200:218]
_GENE_B_EXON1_PLUS = _GENE_B_RC[18:]  # goes at scaffold_2[100:118]

# gene_c (+strand): 1 exon of 18 bp CDS containing an internal stop (TAA at codon 3)
# codons: ATG TTT TAA ATG GGC CCC = MFLMGP  -- wait, TAA = *, so protein = MF*MGP (6 aa)
_GENE_C_CDS_SEQ = "ATGTTTTAAATGGGCCCC"  # 18 nt, 6 aa: MF*MGP

# gene_d (+strand): 1 exon of 9 bp CDS → 3 aa protein (below min_length=30)
_GENE_D_CDS_SEQ = "ATGTTTCCC"  # 9 nt, 3 aa: MFP


def _build_scaffold_1() -> str:
    """500 bp scaffold with gene_a (+) and gene_d (+)."""
    seq = list("A" * 500)
    # gene_a exon1 CDS: positions 100-112
    for i, nt in enumerate(_GENE_A_CDS_SEQ[:12]):
        seq[100 + i] = nt
    # gene_a exon2 CDS: positions 150-162
    for i, nt in enumerate(_GENE_A_CDS_SEQ[12:24]):
        seq[150 + i] = nt
    # gene_a exon3 CDS: positions 200-212
    for i, nt in enumerate(_GENE_A_CDS_SEQ[24:36]):
        seq[200 + i] = nt
    # gene_d exon/CDS: positions 300-309
    for i, nt in enumerate(_GENE_D_CDS_SEQ):
        seq[300 + i] = nt
    return "".join(seq)


def _build_scaffold_2() -> str:
    """500 bp scaffold with gene_b (-) and gene_c (+)."""
    seq = list("A" * 500)
    # gene_b exon1 (plus-strand stored): positions 100-118
    for i, nt in enumerate(_GENE_B_EXON1_PLUS):
        seq[100 + i] = nt
    # gene_b exon2 (plus-strand stored): positions 200-218
    for i, nt in enumerate(_GENE_B_EXON2_PLUS):
        seq[200 + i] = nt
    # gene_c exon/CDS: positions 300-318
    for i, nt in enumerate(_GENE_C_CDS_SEQ):
        seq[300 + i] = nt
    return "".join(seq)


_GFF3_CONTENT = """\
##gff-version 3
scaffold_1\t.\tgene\t101\t212\t.\t+\t.\tID=gene_a
scaffold_1\t.\tmRNA\t101\t212\t.\t+\t.\tID=tx_a1;Parent=gene_a
scaffold_1\t.\texon\t101\t112\t.\t+\t.\tID=tx_a1.exon1;Parent=tx_a1
scaffold_1\t.\texon\t151\t162\t.\t+\t.\tID=tx_a1.exon2;Parent=tx_a1
scaffold_1\t.\texon\t201\t212\t.\t+\t.\tID=tx_a1.exon3;Parent=tx_a1
scaffold_1\t.\tCDS\t101\t112\t.\t+\t0\tID=tx_a1.CDS1;Parent=tx_a1
scaffold_1\t.\tCDS\t151\t162\t.\t+\t0\tID=tx_a1.CDS2;Parent=tx_a1
scaffold_1\t.\tCDS\t201\t212\t.\t+\t0\tID=tx_a1.CDS3;Parent=tx_a1
scaffold_1\t.\tmRNA\t101\t162\t.\t+\t.\tID=tx_a2;Parent=gene_a
scaffold_1\t.\texon\t101\t112\t.\t+\t.\tID=tx_a2.exon1;Parent=tx_a2
scaffold_1\t.\texon\t151\t162\t.\t+\t.\tID=tx_a2.exon2;Parent=tx_a2
scaffold_1\t.\tCDS\t101\t112\t.\t+\t0\tID=tx_a2.CDS1;Parent=tx_a2
scaffold_1\t.\tCDS\t151\t162\t.\t+\t0\tID=tx_a2.CDS2;Parent=tx_a2
scaffold_1\t.\tgene\t301\t309\t.\t+\t.\tID=gene_d
scaffold_1\t.\tmRNA\t301\t309\t.\t+\t.\tID=tx_d;Parent=gene_d
scaffold_1\t.\texon\t301\t309\t.\t+\t.\tID=tx_d.exon1;Parent=tx_d
scaffold_1\t.\tCDS\t301\t309\t.\t+\t0\tID=tx_d.CDS1;Parent=tx_d
scaffold_2\t.\tgene\t101\t218\t.\t-\t.\tID=gene_b
scaffold_2\t.\tmRNA\t101\t218\t.\t-\t.\tID=tx_b;Parent=gene_b
scaffold_2\t.\texon\t101\t118\t.\t-\t.\tID=tx_b.exon1;Parent=tx_b
scaffold_2\t.\texon\t201\t218\t.\t-\t.\tID=tx_b.exon2;Parent=tx_b
scaffold_2\t.\tCDS\t101\t118\t.\t-\t0\tID=tx_b.CDS1;Parent=tx_b
scaffold_2\t.\tCDS\t201\t218\t.\t-\t0\tID=tx_b.CDS2;Parent=tx_b
scaffold_2\t.\tgene\t301\t318\t.\t+\t.\tID=gene_c
scaffold_2\t.\tmRNA\t301\t318\t.\t+\t.\tID=tx_c;Parent=gene_c
scaffold_2\t.\texon\t301\t318\t.\t+\t.\tID=tx_c.exon1;Parent=tx_c
scaffold_2\t.\tCDS\t301\t318\t.\t+\t0\tID=tx_c.CDS1;Parent=tx_c
"""


@pytest.fixture()
def test_data(tmp_path: Path) -> dict[str, Path]:
    """Write synthetic genome FASTA + GFF3; return paths."""
    genome_path = tmp_path / "genome.fa"
    s1 = _build_scaffold_1()
    s2 = _build_scaffold_2()
    genome_path.write_text(f">scaffold_1\n{s1}\n>scaffold_2\n{s2}\n")

    # Build .fai (pyfaidx does it on open, but writing manually is more reliable)
    import pyfaidx
    pyfaidx.Fasta(str(genome_path))  # creates .fai as side effect

    gff3_path = tmp_path / "genes.gff3"
    gff3_path.write_text(_GFF3_CONTENT)

    out_path = tmp_path / "proteins.fa"

    return {"genome": genome_path, "gff3": gff3_path, "out": out_path}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestExtractProteinsBasic:
    def test_plus_strand_protein(self, test_data: dict[str, Path]) -> None:
        results = extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=True, min_length=1,
        )
        # gene_a longest isoform (tx_a1) = 36 nt → 12 aa
        gene_a_key = [k for k in results if k.startswith("gene_a.")][0]
        assert results[gene_a_key] == 12

        text = test_data["out"].read_text()
        # The protein for the 36-nt CDS: MFLIVGPRSTKD
        from helixforge.utils.sequences import translate
        expected = translate(_GENE_A_CDS_SEQ)
        assert expected == "MFLIVGPRSTKD"
        assert expected in text

    def test_minus_strand_protein(self, test_data: dict[str, Path]) -> None:
        results = extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=True, min_length=1,
        )
        gene_b_key = [k for k in results if k.startswith("gene_b.")][0]
        assert results[gene_b_key] == 12

        text = test_data["out"].read_text()
        from helixforge.utils.sequences import translate
        expected = translate(_GENE_B_CDS_DESIRED)
        assert expected == "MLNAGYSINPAD"
        # Verify the minus-strand protein is correct: the actual protein
        # must match what we'd get from translating the designed CDS.
        assert expected in text

    def test_internal_stop_included(self, test_data: dict[str, Path]) -> None:
        results = extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=True, min_length=1,
        )
        gene_c_key = [k for k in results if k.startswith("gene_c.")][0]
        assert results[gene_c_key] == 6

        text = test_data["out"].read_text()
        from helixforge.utils.sequences import translate
        expected = translate(_GENE_C_CDS_SEQ)
        assert "*" in expected
        assert expected in text

    def test_short_protein_skipped(self, test_data: dict[str, Path]) -> None:
        results = extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=True, min_length=30,
        )
        gene_d_keys = [k for k in results if k.startswith("gene_d.")]
        assert gene_d_keys == []


class TestLongestOnly:
    def test_longest_only_keeps_one(self, test_data: dict[str, Path]) -> None:
        results = extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=True, min_length=1,
        )
        gene_a_keys = [k for k in results if k.startswith("gene_a.")]
        assert len(gene_a_keys) == 1
        assert results[gene_a_keys[0]] == 12

    def test_all_isoforms_keeps_both(self, test_data: dict[str, Path]) -> None:
        results = extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=False, min_length=1,
        )
        gene_a_keys = sorted(k for k in results if k.startswith("gene_a."))
        assert len(gene_a_keys) == 2
        lengths = sorted(results[k] for k in gene_a_keys)
        assert lengths == [8, 12]


class TestTranslTable:
    def test_alternative_table_reaches_translate(self, test_data: dict[str, Path]) -> None:
        from helixforge.utils.sequences import translate as _translate
        with mock.patch("helixforge.utils.proteins.translate", wraps=_translate) as m:
            extract_proteins(
                test_data["gff3"], test_data["genome"], test_data["out"],
                longest_only=True, min_length=1, transl_table=4,
            )
            assert m.call_count > 0
            for call in m.call_args_list:
                assert call.kwargs["transl_table"] == 4


class TestOutputFormat:
    def test_fasta_header_format(self, test_data: dict[str, Path]) -> None:
        extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=True, min_length=1,
        )
        lines = test_data["out"].read_text().splitlines()
        headers = [l for l in lines if l.startswith(">")]
        for h in headers:
            assert "length=" in h

    def test_fasta_line_wrapping(self, test_data: dict[str, Path]) -> None:
        extract_proteins(
            test_data["gff3"], test_data["genome"], test_data["out"],
            longest_only=True, min_length=1,
        )
        lines = test_data["out"].read_text().splitlines()
        seq_lines = [l for l in lines if not l.startswith(">")]
        for sl in seq_lines:
            assert len(sl) <= 60
