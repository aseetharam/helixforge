"""Phase 9 D1 — evidence concordance stats.

Concrete literal coordinates only (CLAUDE.md §12); both strands are exercised
for every coordinate/codon path. Synthetic fixtures + a tiny ``MockGenome``.
"""

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    SpliceJunction,
    TranscriptCandidate,
)
from helixforge.stats.evidence_concordance import (
    cds_completeness,
    compute_aed,
    concordance_table,
    evidence_agreement_matrix,
    intron_concordance,
    write_concordance_tsv,
)
from helixforge.utils.sequences import reverse_complement


# ---------------------------------------------------------------------------
# Genome + builders
# ---------------------------------------------------------------------------

# chrP: plus ORF.  CDS [3,18) = ATG AAA CCC GGG TTT ; stop TAA at [18,21).
_CHRP = "CCC" + "ATGAAACCCGGGTTT" + "TAA" + "G" * 176
# chrP_bad: internal stop.  CDS [3,21) = ATG TAA CCC GGG TTT (premature stop).
_CHRP_BAD = "CCC" + "ATGTAACCCGGGTTT" + "TGA" + "G" * 176
# chrP_nostop: codon after CDS [3,18) is GGG, not a stop.
_CHRP_NOSTOP = "CCC" + "ATGAAACCCGGGTTT" + "GGG" + "G" * 176
# chrM: minus ORF.  Reading '-' over CDS [3,18) yields ATG AAA CCC GGG TTT.
_CHRM = "TTA" + reverse_complement("ATGAAACCCGGGTTT") + "G" * 182


class MockGenome:
    def __init__(self, seqs):
        self.seqs = dict(seqs)

    def get_sequence(self, seqid, start, end, strand="+"):
        seq = self.seqs[seqid][start:end]
        return reverse_complement(seq) if strand == "-" else seq


@pytest.fixture
def genome():
    return MockGenome(
        {"chrP": _CHRP, "chrPbad": _CHRP_BAD, "chrPns": _CHRP_NOSTOP, "chrM": _CHRM}
    )


def _classification(locus_id="HFG_00001"):
    return LocusClassification(locus_id=locus_id, status="EXPRESSED", max_tpm=10.0)


def make_tx(
    tid="HFG_00001.1",
    seqid="chr1",
    strand="+",
    exon_bounds=((1000, 1200), (1300, 1500), (1600, 1800)),
    cds=None,
    cds_partial=False,
    tpm=10.0,
    junction_support=1.0,
    confidence=0.8,
    protein_id=None,
    blast_score=None,
    is_primary=True,
):
    return TranscriptCandidate(
        transcript_id=tid,
        locus_id="HFG_00001",
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
        confidence=confidence,
        protein_id=protein_id,
        blast_score=blast_score,
        is_primary=is_primary,
    )


def make_gene(txs, seqid="chr1", strand="+", origin="mikado_1to1", tier=1, flags=None):
    if not isinstance(txs, list):
        txs = [txs]
    return ReconciledGene(
        gene_id="HFG_00001",
        seqid=seqid,
        start=min(t.start for t in txs),
        end=max(t.end for t in txs),
        strand=strand,
        tier=tier,
        transcripts=txs,
        primary_transcript_id=txs[0].transcript_id,
        classification=_classification(),
        origin=origin,
        flags=flags or [],
    )


# ---------------------------------------------------------------------------
# intron_concordance
# ---------------------------------------------------------------------------


def test_intron_all_supported_plus():
    gene = make_gene(make_tx())
    junctions = [
        SpliceJunction("chr1", 1200, 1300, "+", 10),
        SpliceJunction("chr1", 1500, 1600, "+", 10),
    ]
    iso = intron_concordance(gene, junctions)["isoforms"]["HFG_00001.1"]
    assert iso["supported"] == 2
    assert iso["contradicted"] == 0
    assert iso["novel"] == 0
    assert iso["precision"] == 1.0
    assert iso["recall"] == 1.0
    assert iso["f1"] == 1.0


def test_intron_all_supported_minus():
    gene = make_gene(
        make_tx(strand="-", exon_bounds=((2000, 2200), (2300, 2500), (2600, 2800))),
        strand="-",
    )
    junctions = [
        SpliceJunction("chr1", 2200, 2300, "-", 8),
        SpliceJunction("chr1", 2500, 2600, "-", 8),
    ]
    iso = intron_concordance(gene, junctions)["isoforms"]["HFG_00001.1"]
    assert iso["supported"] == 2
    assert iso["precision"] == 1.0
    assert iso["recall"] == 1.0


def test_intron_contradicted_shares_one_boundary():
    gene = make_gene(make_tx())
    junctions = [
        SpliceJunction("chr1", 1200, 1300, "+", 10),   # supports intron 1
        SpliceJunction("chr1", 1500, 1650, "+", 10),   # shares donor of intron 2
    ]
    iso = intron_concordance(gene, junctions)["isoforms"]["HFG_00001.1"]
    assert iso["supported"] == 1
    assert iso["contradicted"] == 1
    assert iso["novel"] == 0
    assert iso["precision"] == 0.5
    assert iso["recall"] == 0.5


def test_intron_novel_when_no_junction_touches():
    gene = make_gene(make_tx())
    iso = intron_concordance(gene, [])["isoforms"]["HFG_00001.1"]
    assert iso["supported"] == 0
    assert iso["contradicted"] == 0
    assert iso["novel"] == 2
    assert iso["precision"] == 0.0
    assert iso["recall"] is None
    assert iso["f1"] is None


def test_intron_min_reads_filter():
    gene = make_gene(make_tx())
    junctions = [
        SpliceJunction("chr1", 1200, 1300, "+", 2),   # below min_reads=3
        SpliceJunction("chr1", 1500, 1600, "+", 5),
    ]
    iso = intron_concordance(gene, junctions, min_reads=3)["isoforms"]["HFG_00001.1"]
    assert iso["supported"] == 1
    assert iso["novel"] == 1


def test_intron_strand_mismatch_not_counted():
    gene = make_gene(make_tx())  # '+' gene
    junctions = [SpliceJunction("chr1", 1200, 1300, "-", 10)]
    iso = intron_concordance(gene, junctions)["isoforms"]["HFG_00001.1"]
    assert iso["supported"] == 0
    assert iso["novel"] == 2


def test_intron_mono_exon_has_no_introns():
    gene = make_gene(make_tx(exon_bounds=((1000, 1800),)))
    iso = intron_concordance(gene, [])["isoforms"]["HFG_00001.1"]
    assert iso["num_introns"] == 0
    assert iso["precision"] is None
    assert iso["recall"] is None
    assert iso["f1"] is None


# ---------------------------------------------------------------------------
# cds_completeness
# ---------------------------------------------------------------------------


def test_cds_completeness_complete_plus(genome):
    gene = make_gene(make_tx(seqid="chrP", exon_bounds=((0, 60),), cds=((3, 18, 0),)), seqid="chrP")
    cc = cds_completeness(gene, genome)["isoforms"]["HFG_00001.1"]
    assert cc["has_cds"] is True
    assert cc["has_start"] is True
    assert cc["has_stop"] is True
    assert cc["internal_stop_free"] is True
    assert cc["cds_length"] == 15


def test_cds_completeness_complete_minus(genome):
    gene = make_gene(
        make_tx(seqid="chrM", strand="-", exon_bounds=((0, 60),), cds=((3, 18, 0),)),
        seqid="chrM",
        strand="-",
    )
    cc = cds_completeness(gene, genome)["isoforms"]["HFG_00001.1"]
    assert cc["has_start"] is True
    assert cc["has_stop"] is True
    assert cc["internal_stop_free"] is True


def test_cds_completeness_internal_stop(genome):
    gene = make_gene(make_tx(seqid="chrPbad", exon_bounds=((0, 60),), cds=((3, 21, 0),)), seqid="chrPbad")
    cc = cds_completeness(gene, genome)["isoforms"]["HFG_00001.1"]
    assert cc["has_start"] is True
    assert cc["internal_stop_free"] is False


def test_cds_completeness_no_stop(genome):
    gene = make_gene(make_tx(seqid="chrPns", exon_bounds=((0, 60),), cds=((3, 18, 0),)), seqid="chrPns")
    cc = cds_completeness(gene, genome)["isoforms"]["HFG_00001.1"]
    assert cc["has_start"] is True
    assert cc["has_stop"] is False


def test_cds_completeness_no_cds(genome):
    gene = make_gene(make_tx(exon_bounds=((1000, 1800),)))
    cc = cds_completeness(gene, genome)["isoforms"]["HFG_00001.1"]
    assert cc["has_cds"] is False
    assert cc["has_start"] is None
    assert cc["cds_length"] == 0


def test_cds_completeness_protein_coverage(genome):
    gene = make_gene(
        make_tx(seqid="chrP", exon_bounds=((0, 60),), cds=((3, 18, 0),),
                protein_id="P1", blast_score=200.0),
        seqid="chrP",
    )
    # Translated protein is 5 aa (MKPGF); reference is 10 aa → coverage 0.5.
    cc = cds_completeness(gene, genome, protein_lengths={"P1": 10})["isoforms"]["HFG_00001.1"]
    assert cc["protein_coverage"] == 0.5


def test_cds_completeness_no_genome_reports_none():
    gene = make_gene(make_tx(seqid="chrP", exon_bounds=((0, 60),), cds=((3, 18, 0),)), seqid="chrP")
    cc = cds_completeness(gene, None)["isoforms"]["HFG_00001.1"]
    assert cc["has_start"] is None
    assert cc["has_stop"] is None
    assert cc["cds_length"] == 15


# ---------------------------------------------------------------------------
# compute_aed
# ---------------------------------------------------------------------------


def test_compute_aed_perfect():
    assert compute_aed(1.0, 1.0, 1.0) == 0.0


def test_compute_aed_worst():
    assert compute_aed(0.0, 0.0, 0.0) == 1.0


def test_compute_aed_partial_signals():
    assert compute_aed(0.5, None, None) == 0.5


def test_compute_aed_no_signals_is_none():
    assert compute_aed(None, None, None) is None


# ---------------------------------------------------------------------------
# evidence_agreement_matrix
# ---------------------------------------------------------------------------


def test_agreement_matrix_identical_source():
    tx = make_tx()
    gene = make_gene(tx)
    sources = {"helixer": [make_tx(tid="HX.1")]}
    m = evidence_agreement_matrix(gene, sources)["HFG_00001.1"]["helixer"]
    assert m["exon_overlap"] == 1.0
    assert m["intron_match"] == 1.0


def test_agreement_matrix_disjoint_source():
    gene = make_gene(make_tx())
    disjoint = make_tx(tid="HX.2", exon_bounds=((5000, 5200), (5300, 5500)))
    m = evidence_agreement_matrix(gene, {"st": [disjoint]})["HFG_00001.1"]["st"]
    assert m["exon_overlap"] == 0.0
    assert m["intron_match"] == 0.0


# ---------------------------------------------------------------------------
# concordance_table
# ---------------------------------------------------------------------------


def test_concordance_table_one_row_per_isoform(genome):
    g1 = make_gene(make_tx(seqid="chrP", exon_bounds=((0, 60),), cds=((3, 18, 0),)), seqid="chrP")
    iso2 = make_tx(tid="HFG_00001.2", exon_bounds=((1000, 1200), (1600, 1800)),
                   is_primary=False)
    iso1 = make_tx(tid="HFG_00001.1")
    g2 = make_gene([iso1, iso2])
    df = concordance_table([g1, g2], [], genome)
    assert len(df) == 3
    assert "intron_f1" in df.columns
    assert "aed" in df.columns


def test_concordance_table_with_sources_adds_columns(genome):
    gene = make_gene(make_tx())
    df = concordance_table([gene], [], genome, sources={"helixer": [make_tx(tid="HX.1")]})
    assert "agree_helixer_exon" in df.columns
    assert "agree_helixer_intron" in df.columns


def test_write_concordance_tsv(tmp_path, genome):
    gene = make_gene(make_tx(seqid="chrP", exon_bounds=((0, 60),), cds=((3, 18, 0),)), seqid="chrP")
    df = concordance_table([gene], [], genome)
    out = tmp_path / "concordance.tsv"
    write_concordance_tsv(df, out)
    text = out.read_text()
    assert "transcript_id" in text.splitlines()[0]
    assert "HFG_00001.1" in text
