"""Tests for mikado/emit_gtf.py (Phase 4). Floor: 10."""

import pytest

from helixforge.mikado.emit_gtf import (
    helixer_gff3_to_gtf,
    stringtie_to_labelled_gtf,
    verify_unique_transcript_ids,
)
from helixforge.reconcile.models import Exon, StringTieTranscript


def _parse_gtf(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        rows.append(line.rstrip("\n").split("\t"))
    return rows


def _st(tid, gid, seqid, strand, exons, tpm=5.0, cov=None):
    start = min(e.start for e in exons)
    end = max(e.end for e in exons)
    return StringTieTranscript(tid, gid, seqid, start, end, strand, exons, tpm, "s1", cov)


# --- helixer_gff3_to_gtf ---

def test_helixer_gtf_returns_path(helixer_gff_path, tmp_path):
    out = helixer_gff3_to_gtf(helixer_gff_path, tmp_path / "h.gtf")
    assert out.exists()


def test_helixer_gtf_transcript_coordinate_conversion(helixer_gff_path, tmp_path):
    out = helixer_gff3_to_gtf(helixer_gff_path, tmp_path / "h.gtf")
    rows = _parse_gtf(out)
    tx = next(r for r in rows if r[2] == "transcript" and "gene1" in r[8])
    # internal (100,300) -> GTF 1-based inclusive (101, 300)
    assert tx[3] == "101"
    assert tx[4] == "300"


def test_helixer_gtf_exon_coordinate_conversion(helixer_gff_path, tmp_path):
    out = helixer_gff3_to_gtf(helixer_gff_path, tmp_path / "h.gtf")
    rows = _parse_gtf(out)
    exons = [r for r in rows if r[2] == "exon" and "gene1" in r[8]]
    assert (exons[0][3], exons[0][4]) == ("101", "150")
    assert (exons[1][3], exons[1][4]) == ("201", "300")


def test_helixer_gtf_has_gene_and_transcript_ids(helixer_gff_path, tmp_path):
    out = helixer_gff3_to_gtf(helixer_gff_path, tmp_path / "h.gtf")
    text = out.read_text()
    assert 'gene_id "gene1";' in text
    assert 'transcript_id "gene1.1";' in text


def test_helixer_gtf_cds_lines_with_phase(helixer_gff_path, tmp_path):
    out = helixer_gff3_to_gtf(helixer_gff_path, tmp_path / "h.gtf")
    rows = _parse_gtf(out)
    cds = [r for r in rows if r[2] == "CDS" and "gene1" in r[8]]
    assert len(cds) == 2
    assert cds[0][7] == "0"  # phase/frame column


def test_helixer_gtf_one_transcript_per_gene(helixer_gff_path, tmp_path):
    out = helixer_gff3_to_gtf(helixer_gff_path, tmp_path / "h.gtf")
    rows = _parse_gtf(out)
    tx_ids = {r[8] for r in rows if r[2] == "transcript"}
    assert len(tx_ids) == 2  # gene1, gene2


def test_helixer_gtf_minus_strand_no_cds(helixer_gff_path, tmp_path):
    out = helixer_gff3_to_gtf(helixer_gff_path, tmp_path / "h.gtf")
    rows = _parse_gtf(out)
    gene2 = [r for r in rows if "gene2" in r[8]]
    assert all(r[6] == "-" for r in gene2)
    assert not any(r[2] == "CDS" for r in gene2)


# --- stringtie_to_labelled_gtf ---

def test_stringtie_gtf_label_in_source(tmp_path):
    txs = [_st("t1", "g1", "chr1", "+", [Exon(100, 200)])]
    out = stringtie_to_labelled_gtf(txs, tmp_path / "st.gtf", "sampleA")
    rows = _parse_gtf(out)
    assert all(r[1] == "sampleA" for r in rows)


def test_stringtie_gtf_coordinate_conversion(tmp_path):
    txs = [_st("t1", "g1", "chr1", "+", [Exon(100, 200), Exon(300, 400)])]
    out = stringtie_to_labelled_gtf(txs, tmp_path / "st.gtf", "sampleA")
    rows = _parse_gtf(out)
    tx = next(r for r in rows if r[2] == "transcript")
    assert (tx[3], tx[4]) == ("101", "400")


def test_stringtie_gtf_preserves_tpm(tmp_path):
    txs = [_st("t1", "g1", "chr1", "+", [Exon(100, 200)], tpm=7.5)]
    out = stringtie_to_labelled_gtf(txs, tmp_path / "st.gtf", "sampleA")
    assert 'TPM "7.5";' in out.read_text()


def test_stringtie_gtf_minus_strand(tmp_path):
    txs = [_st("t1", "g1", "chr1", "-", [Exon(100, 200)])]
    out = stringtie_to_labelled_gtf(txs, tmp_path / "st.gtf", "sampleA")
    rows = _parse_gtf(out)
    assert all(r[6] == "-" for r in rows)


def test_verify_unique_transcript_ids_ok(tmp_path):
    txs = [
        _st("t1", "g1", "chr1", "+", [Exon(100, 200)]),
        _st("t2", "g1", "chr1", "+", [Exon(300, 400)]),
    ]
    verify_unique_transcript_ids(txs)  # no raise


def test_verify_unique_transcript_ids_duplicate_raises(tmp_path):
    txs = [
        _st("t1", "g1", "chr1", "+", [Exon(100, 200)]),
        _st("t1", "g1", "chr1", "+", [Exon(300, 400)]),
    ]
    with pytest.raises(ValueError):
        verify_unique_transcript_ids(txs)


def test_stringtie_gtf_duplicate_ids_raises(tmp_path):
    txs = [
        _st("t1", "g1", "chr1", "+", [Exon(100, 200)]),
        _st("t1", "g1", "chr1", "+", [Exon(300, 400)]),
    ]
    with pytest.raises(ValueError):
        stringtie_to_labelled_gtf(txs, tmp_path / "st.gtf", "sampleA")
