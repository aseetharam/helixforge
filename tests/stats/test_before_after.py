"""Phase 9 D2 — before/after summary + delta table.

Concrete literal coordinates only (CLAUDE.md §12). ``annotation_summary`` is
checked on both a ``ReconciledGene`` list and a parsed Helixer GFF3 so the
before/after comparison is symmetric.
"""

import math

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.stats.before_after import (
    annotation_summary,
    before_after_table,
    busco_completeness,
    write_summary,
)
from helixforge.utils.sequences import reverse_complement


# chrP: plus ORF, CDS [3,18) complete (ATG...TAA).
_CHRP = "CCC" + "ATGAAACCCGGGTTT" + "TAA" + "G" * 176


class MockGenome:
    def __init__(self, seqs):
        self.seqs = dict(seqs)

    def get_sequence(self, seqid, start, end, strand="+"):
        seq = self.seqs[seqid][start:end]
        return reverse_complement(seq) if strand == "-" else seq


@pytest.fixture
def genome():
    return MockGenome({"chrP": _CHRP})


def make_tx(tid, seqid="chr1", strand="+", exon_bounds=((1000, 1200), (1300, 1500)),
            cds=None, cds_partial=False, tpm=10.0, junction_support=1.0):
    return TranscriptCandidate(
        transcript_id=tid,
        locus_id=tid.rsplit(".", 1)[0],
        source="mikado",
        seqid=seqid,
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        exons=[Exon(s, e) for s, e in exon_bounds],
        cds=[CDSSegment(*c) for c in cds] if cds else None,
        cds_partial=cds_partial,
        tpm=tpm,
        junction_support_fraction=junction_support,
    )


def make_gene(gid, txs, seqid="chr1", strand="+", origin="mikado_1to1", tier=1):
    return ReconciledGene(
        gene_id=gid,
        seqid=seqid,
        start=min(t.start for t in txs),
        end=max(t.end for t in txs),
        strand=strand,
        tier=tier,
        transcripts=txs,
        primary_transcript_id=txs[0].transcript_id,
        classification=LocusClassification(locus_id=gid, status="EXPRESSED", max_tpm=10.0),
        origin=origin,
    )


@pytest.fixture
def gene_set():
    # gene 1: two isoforms, multi-exon, with a complete CDS on the primary.
    g1 = make_gene(
        "HFG_00001",
        [
            make_tx("HFG_00001.1", seqid="chrP", exon_bounds=((0, 60),), cds=((3, 18, 0),)),
            make_tx("HFG_00001.2", seqid="chrP", exon_bounds=((0, 30), (40, 60)),
                    cds=((3, 18, 0),)),
        ],
        seqid="chrP",
    )
    # gene 2: single isoform, mono-exon, no CDS.
    g2 = make_gene("HFG_00002", [make_tx("HFG_00002.1", exon_bounds=((2000, 2300),))])
    # gene 3: single isoform, multi-exon, CDS.
    g3 = make_gene(
        "HFG_00003",
        [make_tx("HFG_00003.1", seqid="chrP", exon_bounds=((0, 60),), cds=((3, 18, 0),))],
        seqid="chrP",
    )
    return [g1, g2, g3]


_HELIXER_GFF3 = """##gff-version 3
chrP\tHelixer\tgene\t1\t60\t.\t+\t.\tID=g1
chrP\tHelixer\tmRNA\t1\t60\t.\t+\t.\tID=g1.t1;Parent=g1
chrP\tHelixer\texon\t1\t60\t.\t+\t.\tID=g1.t1.e1;Parent=g1.t1
chrP\tHelixer\tCDS\t4\t18\t.\t+\t0\tID=g1.t1.c1;Parent=g1.t1
chr1\tHelixer\tgene\t2001\t2300\t.\t+\t.\tID=g2
chr1\tHelixer\tmRNA\t2001\t2300\t.\t+\t.\tID=g2.t1;Parent=g2
chr1\tHelixer\texon\t2001\t2300\t.\t+\t.\tID=g2.t1.e1;Parent=g2.t1
"""


@pytest.fixture
def helixer_gff3(tmp_path):
    p = tmp_path / "helixer.gff3"
    p.write_text(_HELIXER_GFF3)
    return str(p)


# ---------------------------------------------------------------------------
# annotation_summary
# ---------------------------------------------------------------------------


def test_summary_counts(gene_set):
    s = annotation_summary(gene_set)
    assert s["gene_count"] == 3
    assert s["transcript_count"] == 4
    assert s["coding_gene_count"] == 2


def test_summary_isoform_distribution(gene_set):
    s = annotation_summary(gene_set)
    assert s["isoforms_per_gene_max"] == 2
    assert s["isoforms_per_gene_distribution"] == {2: 1, 1: 2}
    assert s["isoforms_per_gene_mean"] == pytest.approx(4 / 3)
    assert s["isoforms_per_gene_median"] == 1


def test_summary_mono_vs_multi(gene_set):
    s = annotation_summary(gene_set)
    assert s["mono_exon_genes"] == 2   # genes 2 and 3 (single-exon)
    assert s["multi_exon_genes"] == 1  # gene 1 (iso2 has 2 exons)


def test_summary_cds_lengths(gene_set):
    s = annotation_summary(gene_set)
    # Two coding genes, each primary CDS length 15.
    assert s["mean_cds_length"] == 15
    assert s["median_cds_length"] == 15


def test_summary_pct_complete_with_genome(gene_set, genome):
    s = annotation_summary(gene_set, genome=genome)
    # Both coding genes have a complete ORF (ATG..TAA, no internal stop).
    assert s["pct_complete_orfs"] == 100.0


def test_summary_pct_complete_structural_proxy(gene_set):
    # Without a genome, completeness falls back to mod-3 + no broken-ORF flag.
    s = annotation_summary(gene_set)
    assert s["pct_complete_orfs"] == 100.0


def test_summary_mean_aed_present_for_reconciled(gene_set):
    s = annotation_summary(gene_set)
    # junction_support=1.0 and tpm>0 → per-gene AED 0.0 → mean 0.0.
    assert s["mean_aed"] == 0.0


def test_summary_accepts_gff3(helixer_gff3):
    s = annotation_summary(helixer_gff3)
    assert s["gene_count"] == 2
    assert s["transcript_count"] == 2
    assert s["coding_gene_count"] == 1
    # GFF3 input has no evidence fields → mean AED is NaN.
    assert math.isnan(s["mean_aed"])


def test_summary_gff3_mono_multi(helixer_gff3):
    s = annotation_summary(helixer_gff3)
    # both GFF3 genes are single-exon.
    assert s["mono_exon_genes"] == 2
    assert s["multi_exon_genes"] == 0


# ---------------------------------------------------------------------------
# before_after_table + write_summary
# ---------------------------------------------------------------------------


def test_before_after_table_shape(helixer_gff3, gene_set):
    df = before_after_table(helixer_gff3, gene_set)
    assert {"metric", "helixer", "helixforge", "delta"} <= set(df.columns)
    row = df[df["metric"] == "Genes"].iloc[0]
    assert row["helixer"] == 2
    assert row["helixforge"] == 3
    assert row["delta"] == 1


def test_before_after_table_has_stub_rows(helixer_gff3, gene_set):
    df = before_after_table(helixer_gff3, gene_set)
    busco = df[df["metric"] == "BUSCO completeness %"].iloc[0]
    assert math.isnan(busco["helixer"])
    assert math.isnan(busco["helixforge"])


def test_busco_stub_is_nan():
    assert math.isnan(busco_completeness([]))


def test_write_summary_markdown(tmp_path, helixer_gff3, gene_set):
    df = before_after_table(helixer_gff3, gene_set)
    out = tmp_path / "summary.md"
    write_summary(df, out)
    text = out.read_text()
    assert "| Metric | Helixer | HelixForge | Δ |" in text
    assert "Genes" in text
