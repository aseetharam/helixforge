"""Synthetic fixtures for the standalone ``confidence`` scorer (Phase 14).

A tiny combined Helixer HDF5 with known per-base channels and a tiny GFF3 whose
exons land on those channels at literal coordinates. No real data (CLAUDE.md
§12). Both strands present. Channels: [intergenic, UTR, CDS, intron].

chr1 layout (64 bp):  0-9 intergenic, 10-19 CDS, 20-39 intron, 40-49 UTR, 50-63 intergenic
chr2 layout (64 bp):  0-19 intergenic, 20-39 CDS, 40-59 intron, 60-63 UTR
"""

import h5py
import numpy as np
import pytest

_INTERGENIC = [0.85, 0.05, 0.05, 0.05]
_UTR = [0.05, 0.85, 0.05, 0.05]
_CDS = [0.05, 0.05, 0.90, 0.00]
_INTRON = [0.05, 0.05, 0.05, 0.85]


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
def conf_hdf5_path(tmp_path):
    p = tmp_path / "helixer.h5"
    predictions, seqids, start_ends = _build_predictions()
    with h5py.File(p, "w") as f:
        f.create_dataset("predictions", data=predictions)
        f.create_dataset("seqids", data=seqids)
        f.create_dataset("start_ends", data=start_ends)
    return str(p)


@pytest.fixture
def conf_split_paths(tmp_path):
    """(input_h5, predictions_h5) holding the same data as ``conf_hdf5_path``."""
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


# GFF3 (1-based inclusive). Internal (0-based half-open) coords in comments:
#   g1.t1 chr1 +  exons [10,20) CDS + [40,50) UTR ; intron [20,40)  (multi-exon)
#   g2.t1 chr2 -  exons [20,40) CDS + [60,64) UTR ; intron [40,60)  (multi-exon, minus)
#   g3.t1 chr1 +  exon  [10,20) CDS                               (single-exon)
_GFF3 = """\
##gff-version 3
chr1\ttest\tgene\t11\t50\t.\t+\t.\tID=g1
chr1\ttest\tmRNA\t11\t50\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\ttest\texon\t11\t20\t.\t+\t.\tID=g1.t1.e1;Parent=g1.t1
chr1\ttest\texon\t41\t50\t.\t+\t.\tID=g1.t1.e2;Parent=g1.t1
chr2\ttest\tgene\t21\t64\t.\t-\t.\tID=g2
chr2\ttest\tmRNA\t21\t64\t.\t-\t.\tID=g2.t1;Parent=g2
chr2\ttest\texon\t21\t40\t.\t-\t.\tID=g2.t1.e1;Parent=g2.t1
chr2\ttest\texon\t61\t64\t.\t-\t.\tID=g2.t1.e2;Parent=g2.t1
chr1\ttest\tgene\t11\t20\t.\t+\t.\tID=g3
chr1\ttest\tmRNA\t11\t20\t.\t+\t.\tID=g3.t1;Parent=g3
chr1\ttest\texon\t11\t20\t.\t+\t.\tID=g3.t1.e1;Parent=g3.t1
"""


@pytest.fixture
def conf_gff3_path(tmp_path):
    p = tmp_path / "annotation.gff3"
    p.write_text(_GFF3)
    return str(p)


# ==========================================================================
# Phase 15 — RNA-seq evidence fixtures (synthetic BAM / STAR SJ / GFF3 / TPM)
# ==========================================================================
#
# Junctions present in the data (internal 0-based half-open donor/acceptor):
#   chr1 +  (150, 250)   5 reads (BAM) + 8 reads (STAR)  -> merged 13, 2 samples
#   chr1 +  (350, 450)   STAR only, 6 reads              (union, not in any model)
#   chr1 +  (500, 600)   STAR only, 2 reads              (below min_reads=3)
#   chr2 -  (150, 250)   3 reads (BAM)
#
# GFF3 transcripts (internal coords in comments) span both strands and exercise
# supported / contradicted / novel / single-exon cases.

import pysam  # noqa: E402

_BAM_HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}],
}
_QUERY_OPS = {0, 1, 4, 7, 8}
_SPLICED = [(0, 50), (3, 100), (0, 50)]  # 50M 100N 50M -> donor=start+50, acceptor=+100


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


def _ev_reads():
    reads = []
    # chr1 + junction (150, 250): 5 reads
    for i in range(5):
        reads.append(dict(name=f"p{i}", ref_id=0, start=100, cigar=_SPLICED, mapq=60, xs="+"))
    # chr2 - junction (150, 250): 3 reads
    for i in range(3):
        reads.append(dict(name=f"m{i}", ref_id=1, start=100, cigar=_SPLICED, mapq=60, xs="-"))
    return reads


def _ev_reads2():
    # second BAM sample: chr1 + junction (150, 250) with 2 reads
    return [dict(name=f"s2p{i}", ref_id=0, start=100, cigar=_SPLICED, mapq=60, xs="+")
            for i in range(2)]


@pytest.fixture
def ev_bam_path(tmp_path):
    p = str(tmp_path / "ev_sample1.bam")
    _build_bam(p, _ev_reads())
    return p


@pytest.fixture
def ev_bam_path2(tmp_path):
    p = str(tmp_path / "ev_sample2.bam")
    _build_bam(p, _ev_reads2())
    return p


# STAR SJ.out.tab: chrom start(1-based) end(1-based) strand motif annot unique multi overhang
#   chr1 151-250 '+' 8 reads  -> internal (150, 250)  (overlaps BAM junction)
#   chr1 351-450 '+' 6 reads  -> internal (350, 450)
#   chr1 501-600 '+' 2 reads  -> internal (500, 600)  (below min_reads=3, filtered)
_EV_STAR_SJ = """\
chr1\t151\t250\t1\t1\t0\t8\t0\t30
chr1\t351\t450\t1\t1\t0\t6\t0\t30
chr1\t501\t600\t1\t1\t0\t2\t0\t30
"""


@pytest.fixture
def ev_star_sj_path(tmp_path):
    p = tmp_path / "ev_SJ.out.tab"
    p.write_text(_EV_STAR_SJ)
    return str(p)


# GFF3 (1-based inclusive). Internal coords in comments:
#   gA.t1 chr1 +  exons [100,150),[250,300)  intron [150,250)  SUPPORTED
#   gA.t2 chr1 +  exons [100,160),[250,300)  intron [160,250)  CONTRADICTED (shares acceptor 250)
#   gB.t1 chr1 +  exons [600,650),[760,800)  intron [650,760)  NOVEL
#   gS.t1 chr1 +  exon  [100,150)            single-exon (no introns)
#   gM.t1 chr2 -  exons [100,150),[250,300)  intron [150,250)  SUPPORTED (minus)
_EV_GFF3 = """\
##gff-version 3
chr1\ttest\tgene\t101\t300\t.\t+\t.\tID=gA
chr1\ttest\tmRNA\t101\t300\t.\t+\t.\tID=gA.t1;Parent=gA
chr1\ttest\texon\t101\t150\t.\t+\t.\tID=gA.t1.e1;Parent=gA.t1
chr1\ttest\texon\t251\t300\t.\t+\t.\tID=gA.t1.e2;Parent=gA.t1
chr1\ttest\tmRNA\t101\t300\t.\t+\t.\tID=gA.t2;Parent=gA
chr1\ttest\texon\t101\t160\t.\t+\t.\tID=gA.t2.e1;Parent=gA.t2
chr1\ttest\texon\t251\t300\t.\t+\t.\tID=gA.t2.e2;Parent=gA.t2
chr1\ttest\tgene\t601\t800\t.\t+\t.\tID=gB
chr1\ttest\tmRNA\t601\t800\t.\t+\t.\tID=gB.t1;Parent=gB
chr1\ttest\texon\t601\t650\t.\t+\t.\tID=gB.t1.e1;Parent=gB.t1
chr1\ttest\texon\t761\t800\t.\t+\t.\tID=gB.t1.e2;Parent=gB.t1
chr1\ttest\tgene\t101\t150\t.\t+\t.\tID=gS
chr1\ttest\tmRNA\t101\t150\t.\t+\t.\tID=gS.t1;Parent=gS
chr1\ttest\texon\t101\t150\t.\t+\t.\tID=gS.t1.e1;Parent=gS.t1
chr2\ttest\tgene\t101\t300\t.\t-\t.\tID=gM
chr2\ttest\tmRNA\t101\t300\t.\t-\t.\tID=gM.t1;Parent=gM
chr2\ttest\texon\t101\t150\t.\t-\t.\tID=gM.t1.e1;Parent=gM.t1
chr2\ttest\texon\t251\t300\t.\t-\t.\tID=gM.t1.e2;Parent=gM.t1
"""


@pytest.fixture
def ev_gff3_path(tmp_path):
    p = tmp_path / "ev_annotation.gff3"
    p.write_text(_EV_GFF3)
    return str(p)


# StringTie samples. Sample A carries a transcript with the SAME structure as
# gA.t1 (chr1 +, exons 101-150 / 251-300) at TPM 4.0; sample B repeats it at
# TPM 9.0 so the aggregated max TPM (linked by structure hash) is 9.0.
_EV_STRINGTIE_A = """\
chr1\tStringTie\ttranscript\t101\t300\t1000\t+\t.\tgene_id "S.1"; transcript_id "S.1.1"; cov "8.0"; TPM "4.0";
chr1\tStringTie\texon\t101\t150\t1000\t+\t.\tgene_id "S.1"; transcript_id "S.1.1"; exon_number "1";
chr1\tStringTie\texon\t251\t300\t1000\t+\t.\tgene_id "S.1"; transcript_id "S.1.1"; exon_number "2";
"""
_EV_STRINGTIE_B = """\
chr1\tStringTie\ttranscript\t101\t300\t1500\t+\t.\tgene_id "S.9"; transcript_id "S.9.1"; cov "18.0"; TPM "9.0";
chr1\tStringTie\texon\t101\t150\t1500\t+\t.\tgene_id "S.9"; transcript_id "S.9.1"; exon_number "1";
chr1\tStringTie\texon\t251\t300\t1500\t+\t.\tgene_id "S.9"; transcript_id "S.9.1"; exon_number "2";
"""


@pytest.fixture
def ev_stringtie_gtfs(tmp_path):
    """Individual per-sample StringTie GTF paths (the post-FOFN ``--stringtie``
    form ``score_annotation`` now consumes)."""
    a = tmp_path / "ev_sampleA.gtf"
    b = tmp_path / "ev_sampleB.gtf"
    a.write_text(_EV_STRINGTIE_A)
    b.write_text(_EV_STRINGTIE_B)
    return [str(a), str(b)]


# ==========================================================================
# Protein-AED fixtures (synthetic miniprot GFF3; both strands)
# ==========================================================================
#
# miniprot GFF3 is 1-based inclusive. Two perfect alignments, one per strand,
# whose CDS structure matches the gA.t1 (chr1 +) / gM.t1 (chr2 -) models in the
# RNA GFF3 above. Target start 1 → query_coverage = span / max(target_end) = 1.0.
#   MP1 chr1 +  CDS 101-150, 251-300  -> internal CDS [100,150),[250,300); intron (150,250)
#   MP2 chr2 -  CDS 101-150, 251-300  -> internal CDS [100,150),[250,300); intron (150,250)
_EV_MINIPROT_GFF = """\
##gff-version 3
chr1\tminiprot\tmRNA\t101\t300\t500\t+\t.\tID=MP1;Target=P1 1 100;Identity=0.95;Rank=0
chr1\tminiprot\tCDS\t101\t150\t.\t+\t0\tParent=MP1
chr1\tminiprot\tCDS\t251\t300\t.\t+\t0\tParent=MP1
chr2\tminiprot\tmRNA\t101\t300\t500\t-\t.\tID=MP2;Target=P2 1 100;Identity=0.95;Rank=0
chr2\tminiprot\tCDS\t101\t150\t.\t-\t0\tParent=MP2
chr2\tminiprot\tCDS\t251\t300\t.\t-\t0\tParent=MP2
"""


@pytest.fixture
def ev_miniprot_gff(tmp_path):
    p = tmp_path / "ev_miniprot.gff3"
    p.write_text(_EV_MINIPROT_GFF)
    return str(p)
