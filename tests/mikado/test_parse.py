"""Tests for mikado/parse.py (Phase 5). Floor: 14. Both strands."""

import pytest

from helixforge.mikado.parse import (
    _strip_terminal_stop_codon,
    parse_loci_gff3,
    parse_metrics_tsv,
    parse_scores_tsv,
)
from helixforge.reconcile.models import CDSSegment, MikadoLocus


# --- terminal stop-codon trimming (CLAUDE.md §4.5/§5.2) ---

def test_strip_stop_plus_single_segment():
    out = _strip_terminal_stop_codon([CDSSegment(100, 200, 0)], "+")
    assert [(s.start, s.end) for s in out] == [(100, 197)]


def test_strip_stop_minus_single_segment():
    out = _strip_terminal_stop_codon([CDSSegment(100, 200, 0)], "-")
    # minus strand codes high→low; 3' end is the low coordinate
    assert [(s.start, s.end) for s in out] == [(103, 200)]


def test_strip_stop_plus_multi_segment():
    out = _strip_terminal_stop_codon(
        [CDSSegment(100, 130, 0), CDSSegment(200, 260, 0)], "+"
    )
    assert [(s.start, s.end) for s in out] == [(100, 130), (200, 257)]


def test_strip_stop_minus_multi_segment():
    out = _strip_terminal_stop_codon(
        [CDSSegment(100, 160, 1), CDSSegment(200, 260, 0)], "-"
    )
    # 3' coding end is the low segment; phases preserved
    assert [(s.start, s.end, s.phase) for s in out] == [(103, 160, 1), (200, 260, 0)]


def test_strip_stop_plus_intron_spanning():
    # terminal coding segment is only 2 bp → the stop spans the intron
    out = _strip_terminal_stop_codon(
        [CDSSegment(100, 160, 0), CDSSegment(200, 202, 0)], "+"
    )
    assert [(s.start, s.end) for s in out] == [(100, 159)]


def test_strip_stop_minus_intron_spanning():
    out = _strip_terminal_stop_codon(
        [CDSSegment(100, 102, 0), CDSSegment(200, 260, 0)], "-"
    )
    assert [(s.start, s.end) for s in out] == [(201, 260)]


def test_strip_stop_too_short_returns_none():
    assert _strip_terminal_stop_codon([CDSSegment(100, 102, 0)], "+") is None


@pytest.fixture
def stop_loci(tmp_path):
    gff = tmp_path / "stop.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tMikado\tgene\t101\t300\t.\t+\t.\tID=mikado.sg\n"
        "chr1\tMikado\tmRNA\t101\t300\t.\t+\t.\tID=mikado.sg.1;Parent=mikado.sg\n"
        "chr1\tMikado\texon\t101\t300\t.\t+\t.\tID=e;Parent=mikado.sg.1\n"
        "chr1\tMikado\tCDS\t101\t199\t.\t+\t0\tID=c;Parent=mikado.sg.1\n"
    )
    metrics = tmp_path / "stop.metrics.tsv"
    metrics.write_text(
        "tid\thas_start_codon\thas_stop_codon\nmikado.sg.1\tTrue\tTrue\n"
    )
    return str(gff), str(metrics)


def test_parse_keeps_stop_with_metrics(stop_loci):
    gff, metrics = stop_loci
    t = parse_loci_gff3(gff, metrics_tsv=metrics)[0].transcripts[0]
    # CDS 101-199 -> internal (100,199); stop-inclusive: CDS kept as-is
    assert t.cds[-1].end == 199
    assert t.total_cds_length == 99
    assert t.cds_partial is False


def test_parse_keeps_stop_without_metrics(stop_loci):
    gff, _ = stop_loci
    t = parse_loci_gff3(gff)[0].transcripts[0]
    # stop-inclusive: CDS kept as-is regardless of metrics
    assert t.cds[-1].end == 199
    assert t.total_cds_length == 99


# --- TSV parsers ---

def test_parse_metrics_keyed_by_tid(metrics_tsv_path):
    metrics = parse_metrics_tsv(metrics_tsv_path)
    assert set(metrics) == {"mikado.1G.1", "mikado.1G.2", "mikado.2G.1"}


def test_parse_metrics_numeric_coercion(metrics_tsv_path):
    metrics = parse_metrics_tsv(metrics_tsv_path)
    assert metrics["mikado.1G.1"]["cdna_length"] == 140.0


def test_parse_metrics_bool_coercion(metrics_tsv_path):
    metrics = parse_metrics_tsv(metrics_tsv_path)
    assert metrics["mikado.1G.1"]["is_complete"] is True
    assert metrics["mikado.1G.2"]["is_complete"] is False


def test_parse_scores_keyed_by_tid(scores_tsv_path):
    scores = parse_scores_tsv(scores_tsv_path)
    assert scores["mikado.1G.1"]["score"] == 18.5


def test_parse_tsv_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        parse_metrics_tsv("/no/such.tsv")


# --- loci GFF3 ---

def test_parse_loci_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        parse_loci_gff3("/no/such.gff3")


def test_parse_loci_count(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    assert len(loci) == 2


def test_parse_loci_are_mikado_locus(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    assert all(isinstance(g, MikadoLocus) for g in loci)


def test_parse_loci_sorted_by_seqid_start(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    assert [g.locus_id for g in loci] == ["mikado.1G", "mikado.2G"]


def test_parse_locus_plus_coords(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    g1 = loci[0]
    # GFF 101-300 -> internal (100, 300)
    assert (g1.start, g1.end) == (100, 300)
    assert g1.strand == "+"


def test_parse_locus_minus_coords(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    g2 = loci[1]
    assert (g2.start, g2.end) == (100, 250)
    assert g2.strand == "-"


def test_parse_transcripts_source_mikado(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    for g in loci:
        for t in g.transcripts:
            assert t.source == "mikado"


def test_parse_locus_transcript_count(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    assert len(loci[0].transcripts) == 2  # mikado.1G has .1 and .2
    assert len(loci[1].transcripts) == 1


def test_parse_transcript_exons(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t1 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.1")
    assert [(e.start, e.end) for e in t1.exons] == [(100, 150), (200, 300)]


def test_parse_transcript_cds(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t1 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.1")
    # CDS 112-150 -> (111,150); 201-290 -> (200,290)
    assert [(c.start, c.end) for c in t1.cds] == [(111, 150), (200, 290)]


def test_parse_transcript_no_cds(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t2 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.2")
    assert t2.cds is None


def test_parse_minus_strand_cds(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t = loci[1].transcripts[0]
    assert t.strand == "-"
    assert [(c.start, c.end) for c in t.cds] == [(100, 250)]


def test_parse_protein_id_carried(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t1 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.1")
    assert t1.protein_id == "sp|P12345"


def test_parse_protein_id_absent_is_none(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t2 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.2")
    assert t2.protein_id is None


# --- metrics/scores attached ---

def test_parse_attaches_metrics(loci_gff_path, metrics_tsv_path):
    loci = parse_loci_gff3(loci_gff_path, metrics_tsv=metrics_tsv_path)
    assert loci[0].metrics["mikado.1G.1"]["cdna_length"] == 140.0


def test_parse_attaches_scores(loci_gff_path, scores_tsv_path):
    loci = parse_loci_gff3(loci_gff_path, scores_tsv=scores_tsv_path)
    assert loci[0].scores["mikado.1G.1"]["score"] == 18.5


def test_parse_combined_score_from_scores(loci_gff_path, scores_tsv_path):
    loci = parse_loci_gff3(loci_gff_path, scores_tsv=scores_tsv_path)
    t1 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.1")
    assert t1.combined_score == 18.5


def test_parse_blast_score_from_metrics(loci_gff_path, metrics_tsv_path):
    loci = parse_loci_gff3(loci_gff_path, metrics_tsv=metrics_tsv_path)
    t1 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.1")
    assert t1.blast_score == 816.0
    assert t1.has_homology is True  # drives the Tier-1 gate


def test_parse_blast_score_zero_no_homology(loci_gff_path, metrics_tsv_path):
    loci = parse_loci_gff3(loci_gff_path, metrics_tsv=metrics_tsv_path)
    t2 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.2")
    assert t2.blast_score == 0.0
    assert t2.has_homology is False


def test_parse_blast_score_none_without_metrics(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t1 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.1")
    assert t1.blast_score is None


def test_parse_without_tsv_empty_metrics(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    assert loci[0].metrics == {}
    assert loci[0].scores == {}


def test_parse_without_scores_combined_score_none(loci_gff_path):
    loci = parse_loci_gff3(loci_gff_path)
    t1 = next(t for t in loci[0].transcripts if t.transcript_id == "mikado.1G.1")
    assert t1.combined_score is None


# --- mod-3 crash resilience (Fix 1 + Fix 2) ---


@pytest.fixture
def nonmod3_plus_gff_and_metrics(tmp_path):
    """Plus-strand locus whose CDS is 131 bp (131 % 3 = 2) with metrics
    reporting has_start_codon=True, has_stop_codon=True (the crash scenario).
    A second, valid locus verifies surrounding loci still parse."""
    gff = tmp_path / "nonmod3.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # bad locus: CDS 101-231 → internal (100,231) = 131 bp, 131%3=2
        "chr1\tMikado\tgene\t101\t300\t.\t+\t.\tID=mikado.bad\n"
        "chr1\tMikado\tmRNA\t101\t300\t.\t+\t.\tID=mikado.bad.1;Parent=mikado.bad\n"
        "chr1\tMikado\texon\t101\t300\t.\t+\t.\tParent=mikado.bad.1\n"
        "chr1\tMikado\tCDS\t101\t231\t.\t+\t0\tParent=mikado.bad.1\n"
        # good locus: CDS 501-650 → internal (500,650) = 150 bp, mod-3
        "chr1\tMikado\tgene\t501\t700\t.\t+\t.\tID=mikado.good\n"
        "chr1\tMikado\tmRNA\t501\t700\t.\t+\t.\tID=mikado.good.1;Parent=mikado.good\n"
        "chr1\tMikado\texon\t501\t700\t.\t+\t.\tParent=mikado.good.1\n"
        "chr1\tMikado\tCDS\t501\t650\t.\t+\t0\tParent=mikado.good.1\n"
    )
    metrics = tmp_path / "nonmod3.metrics.tsv"
    metrics.write_text(
        "tid\thas_start_codon\thas_stop_codon\n"
        "mikado.bad.1\tTrue\tTrue\n"
        "mikado.good.1\tTrue\tTrue\n"
    )
    return str(gff), str(metrics)


@pytest.fixture
def nonmod3_minus_gff_and_metrics(tmp_path):
    """Minus-strand locus whose CDS is 131 bp with complete metrics."""
    gff = tmp_path / "nonmod3_minus.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # bad: CDS 101-231 → internal (100,231) = 131 bp, minus strand
        "chr2\tMikado\tgene\t101\t300\t.\t-\t.\tID=mikado.badm\n"
        "chr2\tMikado\tmRNA\t101\t300\t.\t-\t.\tID=mikado.badm.1;Parent=mikado.badm\n"
        "chr2\tMikado\texon\t101\t300\t.\t-\t.\tParent=mikado.badm.1\n"
        "chr2\tMikado\tCDS\t101\t231\t.\t-\t0\tParent=mikado.badm.1\n"
        # good: CDS 501-650 → internal (500,650) = 150 bp, mod-3
        "chr2\tMikado\tgene\t501\t700\t.\t-\t.\tID=mikado.goodm\n"
        "chr2\tMikado\tmRNA\t501\t700\t.\t-\t.\tID=mikado.goodm.1;Parent=mikado.goodm\n"
        "chr2\tMikado\texon\t501\t700\t.\t-\t.\tParent=mikado.goodm.1\n"
        "chr2\tMikado\tCDS\t501\t650\t.\t-\t0\tParent=mikado.goodm.1\n"
    )
    metrics = tmp_path / "nonmod3_minus.metrics.tsv"
    metrics.write_text(
        "tid\thas_start_codon\thas_stop_codon\n"
        "mikado.badm.1\tTrue\tTrue\n"
        "mikado.goodm.1\tTrue\tTrue\n"
    )
    return str(gff), str(metrics)


def test_nonmod3_plus_overridden_to_partial(nonmod3_plus_gff_and_metrics):
    """Non-mod-3 CDS with complete metrics is overridden to partial, not fatal."""
    gff, metrics = nonmod3_plus_gff_and_metrics
    loci = parse_loci_gff3(gff, metrics_tsv=metrics)
    assert len(loci) == 2
    bad = next(g for g in loci if g.locus_id == "mikado.bad")
    t = bad.transcripts[0]
    assert t.total_cds_length == 131
    assert t.cds_partial is True
    assert t.cds_partial_5prime is True
    assert t.cds_partial_3prime is True


def test_nonmod3_minus_overridden_to_partial(nonmod3_minus_gff_and_metrics):
    """Minus-strand: same override logic."""
    gff, metrics = nonmod3_minus_gff_and_metrics
    loci = parse_loci_gff3(gff, metrics_tsv=metrics)
    assert len(loci) == 2
    bad = next(g for g in loci if g.locus_id == "mikado.badm")
    t = bad.transcripts[0]
    assert t.total_cds_length == 131
    assert t.cds_partial is True
    assert t.strand == "-"


def test_nonmod3_surrounding_loci_still_parse(nonmod3_plus_gff_and_metrics):
    """Good locus is unaffected by a sibling's non-mod-3 CDS."""
    gff, metrics = nonmod3_plus_gff_and_metrics
    loci = parse_loci_gff3(gff, metrics_tsv=metrics)
    good = next(g for g in loci if g.locus_id == "mikado.good")
    t = good.transcripts[0]
    assert t.total_cds_length == 150
    assert t.cds_partial is False


def test_nonmod3_without_metrics_uses_heuristic(nonmod3_plus_gff_and_metrics):
    """Without metrics the mod-3 heuristic sets partial=True (no crash)."""
    gff, _ = nonmod3_plus_gff_and_metrics
    loci = parse_loci_gff3(gff)
    bad = next(g for g in loci if g.locus_id == "mikado.bad")
    t = bad.transcripts[0]
    assert t.cds_partial is True
    assert t.total_cds_length == 131


@pytest.fixture
def structerror_gff(tmp_path):
    """A locus with one structurally invalid transcript (CDS outside exon) and
    one valid transcript. A second valid locus follows."""
    gff = tmp_path / "structerror.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # locus with a bad transcript (.1: CDS outside exon) + a good one (.2)
        "chr1\tMikado\tgene\t101\t600\t.\t+\t.\tID=mikado.mixG\n"
        "chr1\tMikado\tmRNA\t101\t300\t.\t+\t.\tID=mikado.mixG.1;Parent=mikado.mixG\n"
        "chr1\tMikado\texon\t101\t200\t.\t+\t.\tParent=mikado.mixG.1\n"
        # CDS 201-300 is outside the exon 101-200 → should fail validation
        "chr1\tMikado\tCDS\t201\t300\t.\t+\t0\tParent=mikado.mixG.1\n"
        "chr1\tMikado\tmRNA\t301\t600\t.\t+\t.\tID=mikado.mixG.2;Parent=mikado.mixG\n"
        "chr1\tMikado\texon\t301\t600\t.\t+\t.\tParent=mikado.mixG.2\n"
        "chr1\tMikado\tCDS\t301\t600\t.\t+\t0\tParent=mikado.mixG.2\n"
        # second valid locus
        "chr1\tMikado\tgene\t801\t1000\t.\t+\t.\tID=mikado.okG\n"
        "chr1\tMikado\tmRNA\t801\t1000\t.\t+\t.\tID=mikado.okG.1;Parent=mikado.okG\n"
        "chr1\tMikado\texon\t801\t1000\t.\t+\t.\tParent=mikado.okG.1\n"
        "chr1\tMikado\tCDS\t801\t1000\t.\t+\t0\tParent=mikado.okG.1\n"
    )
    return str(gff)


def test_structural_error_skips_bad_transcript(structerror_gff):
    """A transcript with CDS outside exon is skipped; sibling and other loci parse."""
    loci = parse_loci_gff3(structerror_gff)
    # mixG keeps the good transcript (.2); bad (.1) is skipped
    mix = next(g for g in loci if g.locus_id == "mikado.mixG")
    assert len(mix.transcripts) == 1
    assert mix.transcripts[0].transcript_id == "mikado.mixG.2"
    # okG is fine
    ok = next(g for g in loci if g.locus_id == "mikado.okG")
    assert len(ok.transcripts) == 1


def test_all_transcripts_fail_skips_locus(tmp_path):
    """If every transcript in a locus fails, the locus is skipped entirely."""
    gff = tmp_path / "allfail.gff3"
    gff.write_text(
        "##gff-version 3\n"
        # every transcript has CDS outside its exon
        "chr1\tMikado\tgene\t101\t600\t.\t+\t.\tID=mikado.failG\n"
        "chr1\tMikado\tmRNA\t101\t300\t.\t+\t.\tID=mikado.failG.1;Parent=mikado.failG\n"
        "chr1\tMikado\texon\t101\t200\t.\t+\t.\tParent=mikado.failG.1\n"
        "chr1\tMikado\tCDS\t201\t300\t.\t+\t0\tParent=mikado.failG.1\n"
        # good locus
        "chr1\tMikado\tgene\t801\t1000\t.\t+\t.\tID=mikado.okG\n"
        "chr1\tMikado\tmRNA\t801\t1000\t.\t+\t.\tID=mikado.okG.1;Parent=mikado.okG\n"
        "chr1\tMikado\texon\t801\t1000\t.\t+\t.\tParent=mikado.okG.1\n"
        "chr1\tMikado\tCDS\t801\t1000\t.\t+\t0\tParent=mikado.okG.1\n"
    )
    loci = parse_loci_gff3(str(gff))
    assert len(loci) == 1
    assert loci[0].locus_id == "mikado.okG"


def test_complete_mod3_not_marked_partial(nonmod3_plus_gff_and_metrics):
    """A truly complete mod-3 CDS is NOT marked partial by the override."""
    gff, metrics = nonmod3_plus_gff_and_metrics
    loci = parse_loci_gff3(gff, metrics_tsv=metrics)
    good = next(g for g in loci if g.locus_id == "mikado.good")
    t = good.transcripts[0]
    assert t.cds_partial is False
    assert t.total_cds_length % 3 == 0
