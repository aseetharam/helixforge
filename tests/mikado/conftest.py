"""Synthetic fixtures for Phase 4 Mikado-emitter tests. No real Mikado/Portcullis."""

import h5py
import numpy as np
import pytest

# Helixer GFF3 for emit_gtf:
#   gene1 chr1 + 101-300  exons 101-150,201-300  CDS 101-150,201-300
#   gene2 chr2 - 101-200  exon 101-200  (no CDS)
_HELIXER_GFF = """\
##gff-version 3
chr1\tHelixer\tgene\t101\t300\t.\t+\t.\tID=gene1
chr1\tHelixer\tmRNA\t101\t300\t.\t+\t.\tID=gene1.m;Parent=gene1
chr1\tHelixer\texon\t101\t150\t.\t+\t.\tID=gene1.e1;Parent=gene1.m
chr1\tHelixer\texon\t201\t300\t.\t+\t.\tID=gene1.e2;Parent=gene1.m
chr1\tHelixer\tCDS\t101\t150\t.\t+\t0\tID=gene1.c1;Parent=gene1.m
chr1\tHelixer\tCDS\t201\t300\t.\t+\t0\tID=gene1.c2;Parent=gene1.m
chr2\tHelixer\tgene\t101\t200\t.\t-\t.\tID=gene2
chr2\tHelixer\tmRNA\t101\t200\t.\t-\t.\tID=gene2.m;Parent=gene2
chr2\tHelixer\texon\t101\t200\t.\t-\t.\tID=gene2.e1;Parent=gene2.m
"""


@pytest.fixture
def helixer_gff_path(tmp_path):
    p = tmp_path / "helixer.gff3"
    p.write_text(_HELIXER_GFF)
    return str(p)


# HDF5 confidence for emit_external. chr1 length 512:
#   [0,100) intergenic, [100,150) CDS, [150,200) intron(high), [200,300) CDS,
#   [300,450) intergenic, [450,500) CDS, [500,512) intergenic
_INTERGENIC = [0.85, 0.05, 0.05, 0.05]
_CDS = [0.05, 0.05, 0.90, 0.00]
_INTRON = [0.05, 0.05, 0.05, 0.85]
_CHUNK = 128


def _vec(pos):
    if pos < 100:
        return _INTERGENIC
    if pos < 150:
        return _CDS
    if pos < 200:
        return _INTRON
    if pos < 300:
        return _CDS
    if pos < 450:
        return _INTERGENIC
    if pos < 500:
        return _CDS
    return _INTERGENIC


@pytest.fixture
def ext_h5_path(tmp_path):
    p = tmp_path / "helixer.h5"
    length = 512
    seqids, start_ends, rows = [], [], []
    for low in range(0, length, _CHUNK):
        high = min(low + _CHUNK, length)
        row = np.zeros((_CHUNK, 4), dtype=np.float32)
        for j in range(high - low):
            row[j] = _vec(low + j)
        rows.append(row)
        seqids.append(b"chr1")
        start_ends.append([low, high])
    with h5py.File(p, "w") as f:
        f.create_dataset("predictions", data=np.stack(rows))
        f.create_dataset("seqids", data=np.array(seqids))
        f.create_dataset("start_ends", data=np.array(start_ends, dtype=np.int64))
    return str(p)


# ==========================================================================
# Phase 5 — Mikado output fixtures (parse.py)
# ==========================================================================

# mikado.loci.gff3:
#   chr1 + locus mikado.1G  with two transcripts:
#     .1  exons 101-150,201-300  CDS 112-150,201-290 (39+90=129, mod-3)  protein_id
#     .2  single exon 101-300, no CDS
#   chr2 - locus mikado.2G  one transcript: exon 101-250, CDS 101-250 (150, mod-3)
_LOCI_GFF = """\
##gff-version 3
chr1\tMikado\tgene\t101\t300\t.\t+\t.\tID=mikado.1G
chr1\tMikado\tmRNA\t101\t300\t.\t+\t.\tID=mikado.1G.1;Parent=mikado.1G;protein_id=sp|P12345
chr1\tMikado\texon\t101\t150\t.\t+\t.\tParent=mikado.1G.1
chr1\tMikado\texon\t201\t300\t.\t+\t.\tParent=mikado.1G.1
chr1\tMikado\tCDS\t112\t150\t.\t+\t0\tParent=mikado.1G.1
chr1\tMikado\tCDS\t201\t290\t.\t+\t0\tParent=mikado.1G.1
chr1\tMikado\tmRNA\t101\t300\t.\t+\t.\tID=mikado.1G.2;Parent=mikado.1G
chr1\tMikado\texon\t101\t300\t.\t+\t.\tParent=mikado.1G.2
chr2\tMikado\tgene\t101\t250\t.\t-\t.\tID=mikado.2G
chr2\tMikado\tmRNA\t101\t250\t.\t-\t.\tID=mikado.2G.1;Parent=mikado.2G
chr2\tMikado\texon\t101\t250\t.\t-\t.\tParent=mikado.2G.1
chr2\tMikado\tCDS\t101\t250\t.\t-\t0\tParent=mikado.2G.1
"""

_METRICS_TSV = """\
tid\tparent\tcdna_length\tcombined_cds_length\tis_complete\tblast_score
mikado.1G.1\tmikado.1G\t140\t129\tTrue\t816.0
mikado.1G.2\tmikado.1G\t200\t0\tFalse\t0.0
mikado.2G.1\tmikado.2G\t150\t150\tTrue\t120.5
"""

_SCORES_TSV = """\
tid\tparent\tscore
mikado.1G.1\tmikado.1G\t18.5
mikado.1G.2\tmikado.1G\t12.0
mikado.2G.1\tmikado.2G\t20.0
"""


@pytest.fixture
def loci_gff_path(tmp_path):
    p = tmp_path / "mikado.loci.gff3"
    p.write_text(_LOCI_GFF)
    return str(p)


@pytest.fixture
def metrics_tsv_path(tmp_path):
    p = tmp_path / "mikado.loci.metrics.tsv"
    p.write_text(_METRICS_TSV)
    return str(p)


@pytest.fixture
def scores_tsv_path(tmp_path):
    p = tmp_path / "mikado.loci.scores.tsv"
    p.write_text(_SCORES_TSV)
    return str(p)
