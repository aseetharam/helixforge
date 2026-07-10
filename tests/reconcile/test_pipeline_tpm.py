"""StringTie TPM reaches the reconcile pipeline (parity with `evidence`).

Mikado consumes StringTie TPM only as an external metric and never emits it back,
so reconciled transcripts used to carry tpm=None and the report's
frac_genes_tpm_pass was always 0 despite valid StringTie input. The pipeline now
assigns each transcript the TPM of the best exonic-overlap StringTie structure —
the same rule `evidence` uses. Both strands, synthetic fixtures, literal coords.
"""

import json
from pathlib import Path

import pytest

from helixforge.io.stringtie import best_overlapping_tpm, build_tpm_overlap_index
from helixforge.reconcile.models import Exon
from helixforge.reconcile.pipeline import PipelineConfig, run_pipeline

# g1 chr1 +  101-200 single exon   (overlaps a + StringTie tx)
# g2 chr1 -  301-600 two exons     (overlaps a - StringTie tx)
# g3 chr2 +  101-250 single exon   (no StringTie overlap)
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

# Both strands; TPM 5.0 (+) and 8.0 (-) both pass the 0.5 default.
_STRINGTIE_GTF = """\
chr1\tStringTie\ttranscript\t101\t200\t1000\t+\t.\tgene_id "S1"; transcript_id "S1.1"; TPM "5.0";
chr1\tStringTie\texon\t101\t200\t1000\t+\t.\tgene_id "S1"; transcript_id "S1.1"; TPM "5.0";
chr1\tStringTie\ttranscript\t301\t600\t1000\t-\t.\tgene_id "S2"; transcript_id "S2.1"; TPM "8.0";
chr1\tStringTie\texon\t301\t400\t1000\t-\t.\tgene_id "S2"; transcript_id "S2.1"; TPM "8.0";
chr1\tStringTie\texon\t501\t600\t1000\t-\t.\tgene_id "S2"; transcript_id "S2.1"; TPM "8.0";
"""


@pytest.fixture
def inputs(tmp_path):
    fa = tmp_path / "genome.fasta"
    fa.write_text(">chr1\n" + "ACGT" * 200 + "\n>chr2\n" + "ACGT" * 100 + "\n")
    gff = tmp_path / "helixer.gff3"
    gff.write_text(_HELIXER_GFF)
    gtf = tmp_path / "sampleA.gtf"
    gtf.write_text(_STRINGTIE_GTF)
    return str(fa), str(gff), str(gtf)


# --- shared overlap helpers (both strands) ----------------------------------

def test_overlap_index_and_lookup_both_strands(inputs):
    _, _, gtf = inputs
    index = build_tpm_overlap_index([gtf])
    # + strand model overlapping S1.1
    assert best_overlapping_tpm(index, "chr1", "+", [Exon(101, 200)]) == 5.0
    # - strand model overlapping S2.1
    assert best_overlapping_tpm(
        index, "chr1", "-", [Exon(301, 400), Exon(501, 600)]
    ) == 8.0
    # wrong strand -> no match
    assert best_overlapping_tpm(index, "chr1", "-", [Exon(101, 200)]) is None
    # no overlap -> None
    assert best_overlapping_tpm(index, "chr2", "+", [Exon(101, 250)]) is None


# --- end-to-end basic-mode run ----------------------------------------------

def test_reconcile_assigns_tpm_and_report_passes(tmp_path, inputs):
    fa, gff, gtf = inputs
    cfg = PipelineConfig(
        genome_fasta=fa,
        helixer_gff3=gff,
        stringtie_list=[gtf],  # no protein_db -> basic (backstop) mode
        output_prefix=str(tmp_path / "out"),
        work_dir=str(tmp_path / "work"),
    )
    genes = run_pipeline(cfg)
    by_id = {g.seqid + g.strand: g for g in genes}
    tpm = {(g.seqid, g.strand): g.transcripts[0].tpm for g in genes}
    # The + and - overlapping genes get their StringTie TPM; the chr2 gene none.
    assert tpm[("chr1", "+")] == 5.0
    assert tpm[("chr1", "-")] == 8.0
    assert tpm[("chr2", "+")] is None

    report = json.loads(Path(f"{cfg.output_prefix}.report.json").read_text())
    sup = report["support"]
    assert sup["genes_tpm_pass"] == 2
    assert sup["frac_genes_tpm_pass"] > 0.0
