"""Synthetic I/O fixtures (Phases 1-2).

Tiny FASTA / Helixer GFF3 / HDF5 / StringTie GTF / miniprot GFF3 / BAM / STAR SJ
with **known content at known positions**. No real biological data (CLAUDE.md
§12). Coordinates in comments are stated in both GFF3 (1-based inclusive) and
internal (0-based half-open) where relevant.
"""

import h5py
import numpy as np
import pysam
import pytest

# --- genome sequences (64 bp each) ---
# chr1 repeats AAAACCCCGGGGTTTT: pos 0-3 AAAA, 4-7 CCCC, 8-11 GGGG, 12-15 TTTT, ...
CHR1_SEQ = "AAAACCCCGGGGTTTT" * 4
# chr2 lowercase to exercise upper-casing.
CHR2_SEQ = "acgt" * 16


def _write_fasta(path, sequences):
    with open(path, "w") as fh:
        for sid, seq in sequences.items():
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


@pytest.fixture
def fasta_path(tmp_path):
    p = tmp_path / "genome.fa"
    _write_fasta(p, {"chr1": CHR1_SEQ, "chr2": CHR2_SEQ})
    return str(p)


# --- Helixer GFF3 ---
# gene1  chr1 +  : gene 11-40  exons 11-20,31-40  CDS 12-20,31-39
# gene2  chr2 -  : gene 11-40  exons 11-25,31-40  CDS 11-25,31-39
# gene3  chr2 +  : gene 50-58  exon 50-58  (no CDS — optional)
_HELIXER_GFF = """\
##gff-version 3
chr1\tHelixer\tgene\t11\t40\t.\t+\t.\tID=gene1
chr1\tHelixer\tmRNA\t11\t40\t.\t+\t.\tID=gene1.mRNA;Parent=gene1
chr1\tHelixer\texon\t11\t20\t.\t+\t.\tID=gene1.exon1;Parent=gene1.mRNA
chr1\tHelixer\texon\t31\t40\t.\t+\t.\tID=gene1.exon2;Parent=gene1.mRNA
chr1\tHelixer\tCDS\t12\t20\t.\t+\t0\tID=gene1.cds1;Parent=gene1.mRNA
chr1\tHelixer\tCDS\t31\t39\t.\t+\t0\tID=gene1.cds2;Parent=gene1.mRNA
chr2\tHelixer\tgene\t11\t40\t.\t-\t.\tID=gene2
chr2\tHelixer\tmRNA\t11\t40\t.\t-\t.\tID=gene2.mRNA;Parent=gene2
chr2\tHelixer\texon\t11\t25\t.\t-\t.\tID=gene2.exon1;Parent=gene2.mRNA
chr2\tHelixer\texon\t31\t40\t.\t-\t.\tID=gene2.exon2;Parent=gene2.mRNA
chr2\tHelixer\tCDS\t11\t25\t.\t-\t0\tID=gene2.cds1;Parent=gene2.mRNA
chr2\tHelixer\tCDS\t31\t39\t.\t-\t0\tID=gene2.cds2;Parent=gene2.mRNA
chr2\tHelixer\tgene\t50\t58\t.\t+\t.\tID=gene3
chr2\tHelixer\tmRNA\t50\t58\t.\t+\t.\tID=gene3.mRNA;Parent=gene3
chr2\tHelixer\texon\t50\t58\t.\t+\t.\tID=gene3.exon1;Parent=gene3.mRNA
"""


@pytest.fixture
def helixer_gff_path(tmp_path):
    p = tmp_path / "helixer.gff3"
    p.write_text(_HELIXER_GFF)
    return str(p)


# --- Helixer HDF5 ---
# Two scaffolds (chr1, chr2), 64 bp each, chunk length 20 → 4 chunks each.
# Channels: [intergenic, UTR, CDS, intron].
_INTERGENIC = [0.85, 0.05, 0.05, 0.05]
_UTR = [0.05, 0.85, 0.05, 0.05]
_CDS = [0.05, 0.05, 0.90, 0.00]
_INTRON = [0.05, 0.05, 0.05, 0.85]

# chr1: 0-9 intergenic, 10-19 CDS, 20-39 intron, 40-49 UTR, 50-63 intergenic
# chr2: 0-19 intergenic, 20-39 CDS, 40-59 intron, 60-63 UTR
def _vec(sid, pos):
    if sid == "chr1":
        if pos < 10:
            return _INTERGENIC
        if pos < 20:
            return _CDS
        if pos < 40:
            return _INTRON
        if pos < 50:
            return _UTR
        return _INTERGENIC
    if pos < 20:
        return _INTERGENIC
    if pos < 40:
        return _CDS
    if pos < 60:
        return _INTRON
    return _UTR


_CHUNK_LEN = 20
_SCAFFOLDS = {"chr1": 64, "chr2": 64}


def _build_predictions():
    seqids, start_ends, rows = [], [], []
    for sid, length in _SCAFFOLDS.items():
        for low in range(0, length, _CHUNK_LEN):
            high = min(low + _CHUNK_LEN, length)
            row = np.zeros((_CHUNK_LEN, 4), dtype=np.float32)
            for j in range(high - low):
                row[j] = _vec(sid, low + j)
            rows.append(row)
            seqids.append(sid.encode())
            start_ends.append([low, high])
    return (
        np.stack(rows),
        np.array(seqids),
        np.array(start_ends, dtype=np.int64),
    )


@pytest.fixture
def hdf5_path(tmp_path):
    p = tmp_path / "helixer.h5"
    predictions, seqids, start_ends = _build_predictions()
    with h5py.File(p, "w") as f:
        f.create_dataset("predictions", data=predictions)
        f.create_dataset("seqids", data=seqids)
        f.create_dataset("start_ends", data=start_ends)
    return str(p)


@pytest.fixture
def hdf5_reverse_path(tmp_path):
    """A single reverse-orientation chunk: chrR covers genomic [0,10), reversed.

    CDS is placed at genomic positions 0-4 (chunk-local indices 9..5).
    """
    p = tmp_path / "helixer_rev.h5"
    row = _reverse_row()
    with h5py.File(p, "w") as f:
        f.create_dataset("predictions", data=row[np.newaxis, :, :])
        f.create_dataset("seqids", data=np.array([b"chrR"]))
        f.create_dataset("start_ends", data=np.array([[10, 0]], dtype=np.int64))
    return str(p)


def _reverse_row():
    row = np.zeros((10, 4), dtype=np.float32)
    for j in range(10):
        # local j -> genomic 9 - j; CDS at genomic 0-4 -> local 9..5
        genomic = 9 - j
        row[j] = _CDS if genomic < 5 else _INTERGENIC
    return row


# --- Phase 13: native split-Helixer layout fixtures ---
# Same underlying values as the combined fixtures, but metadata in a
# *_input.h5 (data/seqids + data/start_ends) and softmax in a *_predictions.h5.


@pytest.fixture
def hdf5_split_paths(tmp_path):
    """(input_h5, predictions_h5) holding the same data as ``hdf5_path``."""
    predictions, seqids, start_ends = _build_predictions()
    inp = tmp_path / "Sample_input.h5"
    pred = tmp_path / "Sample_predictions.h5"
    with h5py.File(inp, "w") as f:
        g = f.create_group("data")
        g.create_dataset("seqids", data=seqids)
        g.create_dataset("start_ends", data=start_ends)
    with h5py.File(pred, "w") as f:
        f.create_dataset("predictions", data=predictions)
    return str(inp), str(pred)


@pytest.fixture
def hdf5_reverse_split_paths(tmp_path):
    """(input_h5, predictions_h5) holding the same data as ``hdf5_reverse_path``."""
    inp = tmp_path / "Rev_input.h5"
    pred = tmp_path / "Rev_predictions.h5"
    with h5py.File(inp, "w") as f:
        g = f.create_group("data")
        g.create_dataset("seqids", data=np.array([b"chrR"]))
        g.create_dataset("start_ends", data=np.array([[10, 0]], dtype=np.int64))
    with h5py.File(pred, "w") as f:
        f.create_dataset("predictions", data=_reverse_row()[np.newaxis, :, :])
    return str(inp), str(pred)


# ==========================================================================
# Phase 2 — evidence fixtures
# ==========================================================================

# --- StringTie GTF (sample A) ---
# STRG.1.1 chr1 +  TPM 5.5  exons 101-150, 201-300  -> internal (100,150),(200,300)
# STRG.2.1 chr2 -  TPM 0.0  (filtered)              exon 101-200
# STRG.3.1 chr1 +  TPM 0.3  single exon 501-600     -> internal (500,600)
_STRINGTIE_A = """\
chr1\tStringTie\ttranscript\t101\t300\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; cov "10.0"; TPM "5.5";
chr1\tStringTie\texon\t101\t150\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1";
chr1\tStringTie\texon\t201\t300\t1000\t+\t.\tgene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "2";
chr2\tStringTie\ttranscript\t101\t200\t0\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1"; cov "0.0"; TPM "0.0";
chr2\tStringTie\texon\t101\t200\t0\t-\t.\tgene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1";
chr1\tStringTie\ttranscript\t501\t600\t50\t+\t.\tgene_id "STRG.3"; transcript_id "STRG.3.1"; cov "0.5"; TPM "0.3";
chr1\tStringTie\texon\t501\t600\t50\t+\t.\tgene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1";
"""

# StringTie sample B: STRG.9.1 has identical structure to STRG.1.1 (TPM 7.0).
_STRINGTIE_B = """\
chr1\tStringTie\ttranscript\t101\t300\t1200\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.1"; cov "12.0"; TPM "7.0";
chr1\tStringTie\texon\t101\t150\t1200\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.1"; exon_number "1";
chr1\tStringTie\texon\t201\t300\t1200\t+\t.\tgene_id "STRG.9"; transcript_id "STRG.9.1"; exon_number "2";
"""


@pytest.fixture
def stringtie_gtf_a(tmp_path):
    p = tmp_path / "sampleA.gtf"
    p.write_text(_STRINGTIE_A)
    return str(p)


@pytest.fixture
def stringtie_gtf_b(tmp_path):
    p = tmp_path / "sampleB.gtf"
    p.write_text(_STRINGTIE_B)
    return str(p)


@pytest.fixture
def stringtie_sample_list(tmp_path, stringtie_gtf_a, stringtie_gtf_b):
    p = tmp_path / "samples.txt"
    p.write_text(f"# stringtie samples\n{stringtie_gtf_a}\n\n{stringtie_gtf_b}\n")
    return str(p)


# --- miniprot GFF3 ---
# MP1 chr1 +  Rank 1 Id 0.95  mRNA 101-300  CDS 101-150,201-300  Target protA 1 100
# MP4 chr1 +  Rank 2 Id 0.90  mRNA 120-280  CDS 120-280          Target protC 1 60
# MP3 chr1 +  Rank 2 Id 0.60  mRNA 401-500  CDS 401-500          Target protA 40 100
# MP2 chr2 -  Rank 1 Id 88.0% mRNA 101-250  CDS 101-250          Target protB 1 50
# protA max target_end = 100 (denominator); protB = 50; protC = 60
_MINIPROT = """\
##gff-version 3
chr1\tminiprot\tmRNA\t101\t300\t200\t+\t.\tID=MP1;Rank=1;Identity=0.95;Target=protA 1 100
chr1\tminiprot\tCDS\t101\t150\t100\t+\t0\tParent=MP1;Target=protA 1 17
chr1\tminiprot\tCDS\t201\t300\t100\t+\t0\tParent=MP1;Target=protA 18 50
chr1\tminiprot\tmRNA\t120\t280\t150\t+\t.\tID=MP4;Rank=2;Identity=0.90;Target=protC 1 60
chr1\tminiprot\tCDS\t120\t280\t150\t+\t0\tParent=MP4;Target=protC 1 53
chr1\tminiprot\tmRNA\t401\t500\t80\t+\t.\tID=MP3;Rank=2;Identity=0.60;Target=protA 40 100
chr1\tminiprot\tCDS\t401\t500\t80\t+\t0\tParent=MP3;Target=protA 40 70
chr2\tminiprot\tmRNA\t101\t250\t150\t-\t.\tID=MP2;Rank=1;Identity=88.0%;Target=protB 1 50
chr2\tminiprot\tCDS\t101\t250\t150\t-\t0\tParent=MP2;Target=protB 1 50
"""


@pytest.fixture
def miniprot_gff_path(tmp_path):
    p = tmp_path / "miniprot.gff3"
    p.write_text(_MINIPROT)
    return str(p)


# --- BAM builder ---
_BAM_HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}],
}

# Query-consuming CIGAR ops: M(0) I(1) S(4) =(7) X(8)
_QUERY_OPS = {0, 1, 4, 7, 8}


def _build_bam(path, reads):
    reads = sorted(reads, key=lambda r: (r["ref_id"], r["start"]))
    with pysam.AlignmentFile(path, "wb", header=_BAM_HEADER) as af:
        for spec in reads:
            a = pysam.AlignedSegment(af.header)
            a.query_name = spec["name"]
            a.reference_id = spec["ref_id"]
            a.reference_start = spec["start"]
            a.mapping_quality = spec["mapq"]
            a.cigartuples = spec["cigar"]
            qlen = sum(l for op, l in spec["cigar"] if op in _QUERY_OPS)
            a.query_sequence = "A" * qlen
            a.flag = spec.get("flag", 0)
            if spec.get("xs"):
                a.set_tag("XS", spec["xs"], "A")
            af.write(a)
    pysam.index(path)


_SPLICED = [(0, 50), (3, 100), (0, 50)]  # 50M 100N 50M


def _sample1_reads():
    reads = []
    # 5 plus-strand spliced reads: junction donor 150, acceptor 250
    for i in range(5):
        reads.append(dict(name=f"p{i}", ref_id=0, start=100, cigar=_SPLICED, mapq=60, xs="+"))
    # 3 minus-strand spliced reads: junction donor 450, acceptor 550
    for i in range(3):
        reads.append(dict(name=f"m{i}", ref_id=0, start=400, cigar=_SPLICED, mapq=60, xs="-"))
    # low-MAPQ spliced read (mapq 5): junction donor 650, acceptor 750
    reads.append(dict(name="lowmapq", ref_id=0, start=600, cigar=_SPLICED, mapq=5, xs="+"))
    # short-overhang read (left 5M): junction donor 805, acceptor 905
    reads.append(dict(name="shortoh", ref_id=0, start=800, cigar=[(0, 5), (3, 100), (0, 50)], mapq=60, xs="+"))
    # secondary spliced read (flag 0x100): must be skipped
    reads.append(dict(name="sec", ref_id=0, start=110, cigar=_SPLICED, mapq=60, xs="+", flag=256))
    # spliced read with no XS tag on chr2: unknown strand -> excluded
    reads.append(dict(name="noxs", ref_id=1, start=100, cigar=_SPLICED, mapq=60))
    return reads


def _sample2_reads():
    reads = []
    # plus junction donor 150 acceptor 250 (2 reads) -> overlaps sample1's
    for i in range(2):
        reads.append(dict(name=f"s2p{i}", ref_id=0, start=100, cigar=_SPLICED, mapq=60, xs="+"))
    # unique plus junction donor 350 acceptor 450 (4 reads)
    for i in range(4):
        reads.append(dict(name=f"s2u{i}", ref_id=0, start=300, cigar=_SPLICED, mapq=60, xs="+"))
    return reads


@pytest.fixture
def bam_path(tmp_path):
    p = str(tmp_path / "sample1.bam")
    _build_bam(p, _sample1_reads())
    return p


@pytest.fixture
def bam_path2(tmp_path):
    p = str(tmp_path / "sample2.bam")
    _build_bam(p, _sample2_reads())
    return p


@pytest.fixture
def bam_no_index(tmp_path):
    p = str(tmp_path / "noindex.bam")
    with pysam.AlignmentFile(p, "wb", header=_BAM_HEADER) as af:
        a = pysam.AlignedSegment(af.header)
        a.query_name = "r0"
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigartuples = [(0, 50)]
        a.query_sequence = "A" * 50
        af.write(a)
    # deliberately NOT indexed
    return p


# --- STAR SJ.out.tab ---
# chrom start end strand motif annotated unique multi overhang
# (150,250,'+',8) kept; (450,550,'-',5) kept; strand 0 excluded;
# unique 2 < 3 filtered; overhang 5 < 8 filtered
_STAR_SJ = """\
chr1\t151\t250\t1\t1\t0\t8\t2\t30
chr1\t451\t550\t2\t2\t0\t5\t0\t12
chr1\t701\t800\t0\t1\t0\t10\t0\t20
chr1\t901\t1000\t1\t1\t0\t2\t0\t15
chr1\t100\t200\t1\t1\t0\t9\t0\t5
"""


@pytest.fixture
def star_sj_path(tmp_path):
    p = tmp_path / "SJ.out.tab"
    p.write_text(_STAR_SJ)
    return str(p)
