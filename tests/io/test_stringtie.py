"""Tests for StringTieParser (Phase 2). Floor: 22."""

import pytest

from helixforge.io.stringtie import StringTieParser


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        StringTieParser().parse_gtf("/no/such.gtf", "s")


def test_parse_filters_tpm_zero(stringtie_gtf_a):
    txs = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")
    # STRG.2.1 has TPM 0.0 -> filtered; STRG.1.1 + STRG.3.1 remain
    ids = {t.transcript_id for t in txs}
    assert ids == {"STRG.1.1", "STRG.3.1"}


def test_parse_count(stringtie_gtf_a):
    txs = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")
    assert len(txs) == 2


def test_parse_sorted_by_seqid_start(stringtie_gtf_a):
    txs = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")
    assert [t.transcript_id for t in txs] == ["STRG.1.1", "STRG.3.1"]


def test_coordinate_conversion(stringtie_gtf_a):
    # GTF 101-300 -> internal start 100; span end 300
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    assert t.start == 100
    assert t.end == 300


def test_exon_coordinates(stringtie_gtf_a):
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    assert [(e.start, e.end) for e in t.exons] == [(100, 150), (200, 300)]


def test_exons_sorted(stringtie_gtf_a):
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    starts = [e.start for e in t.exons]
    assert starts == sorted(starts)


def test_tpm_parsed(stringtie_gtf_a):
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    assert t.tpm == 5.5


def test_coverage_parsed(stringtie_gtf_a):
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    assert t.coverage == 10.0


def test_gene_id_parsed(stringtie_gtf_a):
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    assert t.gene_id == "STRG.1"


def test_sample_id_assigned(stringtie_gtf_a):
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    assert t.sample_id == "sampleA"


def test_strand_parsed(stringtie_gtf_a):
    t = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")[0]
    assert t.strand == "+"


def test_single_exon_transcript(stringtie_gtf_a):
    txs = StringTieParser().parse_gtf(stringtie_gtf_a, "sampleA")
    t3 = next(t for t in txs if t.transcript_id == "STRG.3.1")
    assert [(e.start, e.end) for e in t3.exons] == [(500, 600)]
    assert t3.tpm == 0.3


def test_unstranded_transcript_skipped(tmp_path):
    # StringTie emits '.'-strand single-exon transcripts; they violate the +/-
    # model rule and must be skipped, not crash (real-data regression).
    gtf = tmp_path / "unstranded.gtf"
    gtf.write_text(
        'chr1\tStringTie\ttranscript\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "STRG.1.1"; TPM "5.0";\n'
        'chr1\tStringTie\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "STRG.1.1";\n'
        'chr1\tStringTie\ttranscript\t300\t400\t.\t.\t.\tgene_id "g2"; transcript_id "STRG.2.1"; TPM "9.0";\n'
        'chr1\tStringTie\texon\t300\t400\t.\t.\t.\tgene_id "g2"; transcript_id "STRG.2.1";\n'
    )
    txs = StringTieParser().parse_gtf(str(gtf), "s")
    assert {t.transcript_id for t in txs} == {"STRG.1.1"}
    assert all(t.strand in ("+", "-") for t in txs)


def test_parse_sample_list_count(stringtie_sample_list):
    txs = StringTieParser().parse_sample_list(stringtie_sample_list)
    # sampleA: 2, sampleB: 1
    assert len(txs) == 3


def test_parse_sample_list_sample_ids(stringtie_sample_list):
    txs = StringTieParser().parse_sample_list(stringtie_sample_list)
    assert {t.sample_id for t in txs} == {"sampleA", "sampleB"}


def test_parse_sample_list_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        StringTieParser().parse_sample_list("/no/such/list.txt")


def test_aggregate_groups_matching_structure(stringtie_sample_list):
    parser = StringTieParser()
    txs = parser.parse_sample_list(stringtie_sample_list)
    agg = parser.aggregate_across_samples(txs)
    # default min_tpm=0.5 drops the 0.3-TPM single-exon structure
    assert len(agg) == 1


def test_aggregate_num_samples(stringtie_sample_list):
    parser = StringTieParser()
    txs = parser.parse_sample_list(stringtie_sample_list)
    agg = parser.aggregate_across_samples(txs)
    group = next(iter(agg.values()))
    assert group["num_samples"] == 2


def test_aggregate_max_and_mean_tpm(stringtie_sample_list):
    parser = StringTieParser()
    txs = parser.parse_sample_list(stringtie_sample_list)
    agg = parser.aggregate_across_samples(txs)
    group = next(iter(agg.values()))
    assert group["max_tpm"] == 7.0
    assert group["mean_tpm"] == pytest.approx(6.25)


def test_aggregate_sample_tpms(stringtie_sample_list):
    parser = StringTieParser()
    txs = parser.parse_sample_list(stringtie_sample_list)
    agg = parser.aggregate_across_samples(txs)
    group = next(iter(agg.values()))
    assert group["sample_tpms"] == {"sampleA": 5.5, "sampleB": 7.0}


def test_aggregate_representative_is_highest_tpm(stringtie_sample_list):
    parser = StringTieParser()
    txs = parser.parse_sample_list(stringtie_sample_list)
    agg = parser.aggregate_across_samples(txs)
    group = next(iter(agg.values()))
    assert group["representative"].tpm == 7.0


def test_aggregate_min_tpm_keeps_low_structure(stringtie_sample_list):
    parser = StringTieParser()
    txs = parser.parse_sample_list(stringtie_sample_list)
    agg = parser.aggregate_across_samples(txs, min_tpm=0.1)
    assert len(agg) == 2


def test_aggregate_min_samples_filter(stringtie_sample_list):
    parser = StringTieParser()
    txs = parser.parse_sample_list(stringtie_sample_list)
    agg = parser.aggregate_across_samples(txs, min_tpm=0.1, min_samples=2)
    # only the shared structure has 2 samples
    assert len(agg) == 1


def test_aggregate_empty_input():
    assert StringTieParser().aggregate_across_samples([]) == {}
