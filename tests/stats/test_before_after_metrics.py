"""Phase 17 D1 — populated before/after AED + Helixer-support metrics.

Concrete literal coordinates; both strands. The Helixer support comes from a
synthetic HDF5 (same channel layout as the I/O fixtures); AED from concrete
junction / expression / protein-coverage inputs. Floor: 8.
"""

import math

import h5py
import numpy as np
import pytest

from helixforge.io.hdf5 import HDF5ConfidenceReader
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    SpliceJunction,
    TranscriptCandidate,
)
from helixforge.stats.before_after import annotation_summary, before_after_table
from helixforge.utils.sequences import reverse_complement

# Channels: [intergenic, UTR, CDS, intron].
_CDS = [0.05, 0.05, 0.90, 0.00]
_INTERGENIC = [0.85, 0.05, 0.05, 0.05]
_INTRON = [0.05, 0.05, 0.05, 0.85]
_CHUNK = 20
_LEN = 80


def _vec(pos):
    # chr1: CDS over [0,40), intron [40,60), intergenic [60,80).
    if pos < 40:
        return _CDS
    if pos < 60:
        return _INTRON
    return _INTERGENIC


@pytest.fixture
def h5_path(tmp_path):
    rows, seqids, start_ends = [], [], []
    for low in range(0, _LEN, _CHUNK):
        high = min(low + _CHUNK, _LEN)
        row = np.zeros((_CHUNK, 4), dtype=np.float32)
        for j in range(high - low):
            row[j] = _vec(low + j)
        rows.append(row)
        seqids.append(b"chr1")
        start_ends.append([low, high])
    p = tmp_path / "helixer.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("predictions", data=np.stack(rows))
        f.create_dataset("seqids", data=np.array(seqids))
        f.create_dataset("start_ends", data=np.array(start_ends, dtype=np.int64))
    return str(p)


# chrP: plus ORF, CDS [3,18) complete (ATG...TAA).
_CHRP = "CCC" + "ATGAAACCCGGGTTT" + "TAA" + "G" * 176


class MockGenome:
    def __init__(self, seqs):
        self.seqs = dict(seqs)

    def get_sequence(self, seqid, start, end, strand="+"):
        seq = self.seqs[seqid][start:end]
        return reverse_complement(seq) if strand == "-" else seq


def _tx(tid, seqid, strand, exons, cds=None, tpm=10.0, junction=1.0, protein=None):
    return TranscriptCandidate(
        transcript_id=tid, locus_id=tid.rsplit(".", 1)[0], source="mikado",
        seqid=seqid, start=exons[0][0], end=exons[-1][1], strand=strand,
        exons=[Exon(s, e) for s, e in exons],
        cds=[CDSSegment(*c) for c in cds] if cds else None,
        tpm=tpm, junction_support_fraction=junction, protein_id=protein,
    )


def _gene(gid, tx, seqid="chr1", strand="+"):
    return ReconciledGene(
        gene_id=gid, seqid=seqid, start=tx.start, end=tx.end, strand=strand,
        tier=1, transcripts=[tx], primary_transcript_id=tx.transcript_id,
        classification=LocusClassification(gid, "EXPRESSED"), origin="mikado_1to1",
    )


# Single-exon Helixer gene over the CDS-prob region (for support before-set).
_HELIXER_GFF3 = """##gff-version 3
chr1\tHelixer\tgene\t1\t40\t.\t+\t.\tID=g1
chr1\tHelixer\tmRNA\t1\t40\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tHelixer\texon\t1\t40\t.\t+\t.\tID=g1.t1.e1;Parent=g1.t1
chr1\tHelixer\tCDS\t1\t39\t.\t+\t0\tID=g1.t1.c1;Parent=g1.t1
"""

# Two-exon Helixer gene with an intron at internal donor=100, acceptor=200.
_HELIXER_GFF3_SPLICED = """##gff-version 3
chr1\tHelixer\tgene\t1\t300\t.\t+\t.\tID=g1
chr1\tHelixer\tmRNA\t1\t300\t.\t+\t.\tID=g1.t;Parent=g1
chr1\tHelixer\texon\t1\t100\t.\t+\t.\tID=g1.e1;Parent=g1.t
chr1\tHelixer\texon\t201\t300\t.\t+\t.\tID=g1.e2;Parent=g1.t
"""


@pytest.fixture
def helixer_gff3(tmp_path):
    p = tmp_path / "helixer.gff3"
    p.write_text(_HELIXER_GFF3)
    return str(p)


@pytest.fixture
def helixer_gff3_spliced(tmp_path):
    p = tmp_path / "helixer_spliced.gff3"
    p.write_text(_HELIXER_GFF3_SPLICED)
    return str(p)


# ---------------------------------------------------------------------------
# helixer_support
# ---------------------------------------------------------------------------

def test_helixer_support_populated_for_reconciled(h5_path):
    gene = _gene("HFG_00001", _tx("HFG_00001.1", "chr1", "+", [(0, 40)], cds=[(0, 39, 0)]))
    with HDF5ConfidenceReader(h5_path) as r:
        s = annotation_summary([gene], h5_reader=r)
    # exon sits entirely over the CDS-prob=0.90 region.
    assert s["mean_helixer_support"] == pytest.approx(0.90, abs=1e-5)


def test_helixer_support_populated_for_gff3_before_set(h5_path, helixer_gff3):
    with HDF5ConfidenceReader(h5_path) as r:
        s = annotation_summary(helixer_gff3, h5_reader=r)
    assert s["mean_helixer_support"] == pytest.approx(0.90, abs=1e-5)


def test_helixer_support_nan_without_h5():
    gene = _gene("HFG_00001", _tx("HFG_00001.1", "chr1", "+", [(0, 40)], cds=[(0, 39, 0)]))
    s = annotation_summary([gene])
    assert math.isnan(s["mean_helixer_support"])


def test_helixer_support_minus_strand(h5_path):
    gene = _gene("HFG_M", _tx("HFG_M.1", "chr1", "-", [(0, 40)], cds=[(0, 39, 0)]),
                 strand="-")
    with HDF5ConfidenceReader(h5_path) as r:
        s = annotation_summary([gene], h5_reader=r)
    # Helixer channels are strand-agnostic per position → same 0.90 support.
    assert s["mean_helixer_support"] == pytest.approx(0.90, abs=1e-5)


# ---------------------------------------------------------------------------
# AED
# ---------------------------------------------------------------------------

def test_aed_from_junction_and_expression_reconciled():
    # junction 0.5, expression present → aed = 1 - mean(0.5, 1.0) = 0.25.
    gene = _gene("HFG_00001",
                 _tx("HFG_00001.1", "chr1", "+", [(0, 100), (200, 300)], junction=0.5))
    s = annotation_summary([gene])
    assert s["mean_aed"] == pytest.approx(0.25)


def test_aed_for_gff3_from_junction_set(helixer_gff3_spliced):
    jset = {("chr1", 100, 200, "+")}  # internal donor=100, acceptor=200
    s = annotation_summary(helixer_gff3_spliced, junction_set=jset)
    # one intron, fully supported (junction=1.0), no other signal → aed = 0.0.
    assert s["mean_aed"] == pytest.approx(0.0)


def test_aed_gff3_nan_without_junctions(helixer_gff3):
    # Single-exon Helixer gene, no junctions → no signal → NaN.
    s = annotation_summary(helixer_gff3)
    assert math.isnan(s["mean_aed"])


def test_aed_protein_coverage_lowers_distance():
    genome = MockGenome({"chrP": _CHRP})
    tx = _tx("HFG_00001.1", "chrP", "+", [(0, 60)], cds=[(3, 18, 0)],
             junction=0.0, tpm=0.0, protein="P1")
    gene = _gene("HFG_00001", tx, seqid="chrP")
    # translated ORF length 5 aa; ref length 5 → coverage 1.0.
    s_cov = annotation_summary([gene], genome=genome, protein_lengths={"P1": 5})
    s_nocov = annotation_summary([gene], genome=genome)
    # junction 0 + expression 0 → aed 1.0 without coverage; coverage pulls it down.
    assert s_cov["mean_aed"] < s_nocov["mean_aed"]
    assert s_cov["mean_aed"] == pytest.approx(2 / 3)


# ---------------------------------------------------------------------------
# before_after_table wiring
# ---------------------------------------------------------------------------

def test_before_after_table_helixer_support_row(h5_path, helixer_gff3):
    gene = _gene("HFG_00001", _tx("HFG_00001.1", "chr1", "+", [(0, 40)], cds=[(0, 39, 0)]))
    df = before_after_table(helixer_gff3, [gene], h5_path=h5_path)
    row = df[df["metric"] == "Mean Helixer support"].iloc[0]
    assert row["helixer"] == pytest.approx(0.90, abs=1e-5)
    assert row["helixforge"] == pytest.approx(0.90, abs=1e-5)
    assert row["delta"] == pytest.approx(0.0, abs=1e-5)


def test_before_after_table_helixer_support_nan_without_h5(helixer_gff3):
    gene = _gene("HFG_00001", _tx("HFG_00001.1", "chr1", "+", [(0, 40)], cds=[(0, 39, 0)]))
    df = before_after_table(helixer_gff3, [gene])
    row = df[df["metric"] == "Mean Helixer support"].iloc[0]
    assert math.isnan(row["helixer"])
    assert math.isnan(row["helixforge"])


def test_before_after_table_aed_from_junctions(helixer_gff3):
    gene = _gene("HFG_00001",
                 _tx("HFG_00001.1", "chr1", "+", [(0, 100), (200, 300)], junction=1.0))
    junctions = [SpliceJunction("chr1", 100, 200, "+", read_count=10, samples=1)]
    df = before_after_table(helixer_gff3, [gene], junctions=junctions)
    row = df[df["metric"] == "Mean AED"].iloc[0]
    assert row["helixforge"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Defensive guard — a malformed GFF3 with a duplicate gene ID parses into a
# transcript-less gene; annotation_summary must skip it (warn), not crash with
# IndexError (the scorer:stats failure). Both strands.
# ---------------------------------------------------------------------------

_DUP_ID_GFF3 = """\
##gff-version 3
chr1\thf\tgene\t1000\t2000\t.\t{s}\t.\tID=HFG_00010
chr1\thf\tmRNA\t1000\t2000\t.\t{s}\t.\tID=HFG_00010.1;Parent=HFG_00010
chr1\thf\texon\t1000\t2000\t.\t{s}\t.\tID=HFG_00010.1.exon1;Parent=HFG_00010.1
chr1\thf\tgene\t3000\t4000\t.\t{s}\t.\tID=HFG_00010
chr1\thf\tmRNA\t3000\t4000\t.\t{s}\t.\tID=HFG_00010.1;Parent=HFG_00010
chr1\thf\texon\t3000\t4000\t.\t{s}\t.\tID=HFG_00010.1.exon2;Parent=HFG_00010.1
"""


@pytest.mark.parametrize("strand", ["+", "-"])
def test_annotation_summary_skips_transcriptless_gene(tmp_path, caplog, strand):
    p = tmp_path / "dup.gff3"
    p.write_text(_DUP_ID_GFF3.format(s=strand))
    import logging

    with caplog.at_level(logging.WARNING):
        s = annotation_summary(str(p))
    # Did not crash; at least one well-formed gene survived.
    assert s["gene_count"] >= 1
    assert any("transcript-less" in r.message for r in caplog.records)
