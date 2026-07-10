"""Phase 7 D1 — backstop CDS projection + cross-check.

Concrete literal coordinates only (CLAUDE.md §12); both strands mandatory.
TransDecoder is mocked (no external tool runs in unit tests).
"""

import os

import pytest

import helixforge.reconcile.cds as cds_mod
from helixforge.reconcile.cds import (
    _include_genomic_stop,
    assign_backstop_cds,
    cds_cross_check,
    extract_transcript_sequence,
    map_transcript_to_genomic,
    project_cds_to_exons,
    _parse_transdecoder_bed,
)
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    MiniprotAlignment,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.utils.sequences import reverse_complement

CONTIG_LEN = 4000


class MockGenome:
    def __init__(self, sequences):
        self.sequences = dict(sequences)

    def get_sequence(self, seqid, start, end, strand="+"):
        seq = self.sequences[seqid][start:end]
        return reverse_complement(seq) if strand == "-" else seq


def genome_with(planted=None):
    arr = ["A"] * CONTIG_LEN
    for pos, seq in (planted or {}).items():
        for i, ch in enumerate(seq):
            arr[pos + i] = ch
    return MockGenome({"chr1": "".join(arr)})


def exons(bounds):
    return [Exon(s, e) for s, e in bounds]


def segs(bounds):
    return [CDSSegment(*b) if len(b) == 3 else CDSSegment(b[0], b[1], 0) for b in bounds]


def backstop_gene(exon_bounds, strand="+", cds=None, tier=4):
    tx = TranscriptCandidate(
        transcript_id="HFG_00001.1",
        locus_id="HFG_00001",
        source="helixer",
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        exons=exons(exon_bounds),
        cds=segs(cds) if cds else None,
    )
    return ReconciledGene(
        gene_id="HFG_00001",
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        tier=tier,
        transcripts=[tx],
        primary_transcript_id="HFG_00001.1",
        classification=LocusClassification("HFG_00001", "SILENT"),
        origin="helixer_backstop",
    )


def mikado_gene(exon_bounds, strand="+", cds=None):
    tx = TranscriptCandidate(
        transcript_id="HFG_00001.1",
        locus_id="HFG_00001",
        source="mikado",
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        exons=exons(exon_bounds),
        cds=segs(cds) if cds else None,
        is_primary=True,
    )
    return ReconciledGene(
        gene_id="HFG_00001",
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        tier=2,
        transcripts=[tx],
        primary_transcript_id="HFG_00001.1",
        classification=LocusClassification("HFG_00001", "EXPRESSED"),
        origin="mikado_1to1",
    )


def alignment(cds_bounds, strand="+", start=None, end=None, pid="sp|TEST"):
    cs = segs(cds_bounds)
    return MiniprotAlignment(
        protein_id=pid,
        seqid="chr1",
        start=start if start is not None else cs[0].start,
        end=end if end is not None else cs[-1].end,
        strand=strand,
        cds_segments=cs,
        query_coverage=0.9,
        identity=0.85,
        score=500.0,
        rank=0,
    )


def bounds(segments):
    return [(s.start, s.end) for s in segments]


def phases(segments):
    return [s.phase for s in segments]


# ---------------------------------------------------------------------------
# project_cds_to_exons
# ---------------------------------------------------------------------------

def test_project_single_segment_within_exon():
    out = project_cds_to_exons(segs([(1000, 1099)]), exons([(1000, 1100)]), "+")
    assert bounds(out) == [(1000, 1099)]
    assert phases(out) == [0]


def test_project_drops_intronic_portion():
    out = project_cds_to_exons(
        segs([(1000, 1400)]), exons([(1000, 1200), (1300, 1500)]), "+"
    )
    assert bounds(out) == [(1000, 1200), (1300, 1400)]  # 200 + 100 = 300
    assert phases(out) == [0, 1]  # (3 - 200 % 3) % 3 = 1


def test_project_returns_none_when_not_mod3():
    out = project_cds_to_exons(segs([(1000, 1101)]), exons([(1000, 1100)]), "+")
    assert out is None  # clipped length 100, not divisible by 3


def test_project_returns_none_when_empty():
    out = project_cds_to_exons(segs([(5000, 5100)]), exons([(1000, 1200)]), "+")
    assert out is None


def test_project_three_segments_plus_phases():
    out = project_cds_to_exons(
        segs([(1000, 1031), (1100, 1131), (1200, 1228)]),
        exons([(1000, 1031), (1100, 1131), (1200, 1228)]),
        "+",
    )
    assert bounds(out) == [(1000, 1031), (1100, 1131), (1200, 1228)]
    assert phases(out) == [0, 2, 1]


def test_project_three_segments_minus_phases():
    out = project_cds_to_exons(
        segs([(1000, 1031), (1100, 1131), (1200, 1228)]),
        exons([(1000, 1031), (1100, 1131), (1200, 1228)]),
        "-",
    )
    # coordinates unchanged (low->high); phases computed in coding (reverse) order
    assert bounds(out) == [(1000, 1031), (1100, 1131), (1200, 1228)]
    assert phases(out) == [1, 2, 0]


def test_project_merges_adjacent_pieces():
    # two CDS segments both inside one exon, abutting -> merged into one
    out = project_cds_to_exons(
        segs([(1000, 1050), (1050, 1099)]), exons([(1000, 1100)]), "+"
    )
    assert bounds(out) == [(1000, 1099)]


def test_project_default_strand_is_plus():
    out = project_cds_to_exons(segs([(1000, 1099)]), exons([(1000, 1100)]))
    assert phases(out) == [0]


# ---------------------------------------------------------------------------
# extract_transcript_sequence — both strands
# ---------------------------------------------------------------------------

def test_extract_transcript_sequence_plus():
    g = genome_with({1000: "ATGAAA", 1100: "GGGTTT"})
    tx = backstop_gene([(1000, 1006), (1100, 1106)], "+").transcripts[0]
    assert extract_transcript_sequence(tx, g) == "ATGAAAGGGTTT"


def test_extract_transcript_sequence_minus():
    g = genome_with({1000: "ATGAAA", 1100: "GGGTTT"})
    tx = backstop_gene([(1000, 1006), (1100, 1106)], "-").transcripts[0]
    # coding order is descending; each piece reverse-complemented
    assert extract_transcript_sequence(tx, g) == reverse_complement(
        "ATGAAA" + "GGGTTT"
    )


# ---------------------------------------------------------------------------
# map_transcript_to_genomic — both strands
# ---------------------------------------------------------------------------

def test_map_transcript_to_genomic_plus():
    tx = backstop_gene([(1000, 1100), (1200, 1300)], "+").transcripts[0]
    out = map_transcript_to_genomic(tx, 0, 150)
    assert bounds(out) == [(1000, 1100), (1200, 1250)]
    assert phases(out) == [0, 2]


def test_map_transcript_to_genomic_minus():
    tx = backstop_gene([(1000, 1100), (1200, 1300)], "-").transcripts[0]
    out = map_transcript_to_genomic(tx, 0, 150)
    assert bounds(out) == [(1050, 1100), (1200, 1300)]
    assert phases(out) == [2, 0]


def test_map_transcript_to_genomic_out_of_range():
    tx = backstop_gene([(1000, 1100), (1200, 1300)], "+").transcripts[0]
    assert map_transcript_to_genomic(tx, 0, 500) is None


def test_map_transcript_to_genomic_not_mod3():
    tx = backstop_gene([(1000, 1100), (1200, 1300)], "+").transcripts[0]
    assert map_transcript_to_genomic(tx, 0, 100) is None


def test_map_transcript_to_genomic_within_single_exon_minus():
    tx = backstop_gene([(1000, 1100), (1200, 1300)], "-").transcripts[0]
    out = map_transcript_to_genomic(tx, 0, 30)  # first 30 coding bases
    assert bounds(out) == [(1270, 1300)]  # high end of the coding-first exon
    assert phases(out) == [0]


# ---------------------------------------------------------------------------
# assign_backstop_cds — miniprot path
# ---------------------------------------------------------------------------

def test_assign_backstop_miniprot_plus_sets_cds_and_tier():
    gene = backstop_gene([(1000, 1200), (1300, 1500)], "+")
    aln = alignment([(1000, 1150), (1300, 1450)], "+", start=1000, end=1500)
    out = assign_backstop_cds(gene, [aln])
    tx = out.transcripts[0]
    assert bounds(tx.cds) == [(1000, 1150), (1300, 1450)]
    assert tx.protein_id == "sp|TEST"
    assert out.tier == 1


def test_assign_backstop_miniprot_minus():
    gene = backstop_gene([(1000, 1200), (1300, 1500)], "-")
    aln = alignment([(1000, 1150), (1300, 1450)], "-", start=1000, end=1500)
    out = assign_backstop_cds(gene, [aln])
    assert out.transcripts[0].cds is not None
    assert out.tier == 1


def test_assign_backstop_non_backstop_unchanged():
    gene = mikado_gene([(1000, 1200)], "+", cds=[(1000, 1198)])
    aln = alignment([(1000, 1150)], "+", start=1000, end=1200)
    out = assign_backstop_cds(gene, [aln])
    assert out is gene
    assert bounds(out.transcripts[0].cds) == [(1000, 1198)]


def test_assign_backstop_no_overlapping_alignment():
    gene = backstop_gene([(1000, 1200)], "+")
    aln = alignment([(5000, 5150)], "+", start=5000, end=5200)
    out = assign_backstop_cds(gene, [aln])
    assert out.transcripts[0].cds is None


def test_assign_backstop_wrong_strand_alignment_ignored():
    gene = backstop_gene([(1000, 1200), (1300, 1500)], "+")
    aln = alignment([(1000, 1150), (1300, 1450)], "-", start=1000, end=1500)
    out = assign_backstop_cds(gene, [aln])
    assert out.transcripts[0].cds is None


def test_assign_backstop_miniprot_mod3_fail_falls_through():
    gene = backstop_gene([(1000, 1200)], "+")
    aln = alignment([(1000, 1100)], "+", start=1000, end=1200)  # 100 bp, not mod3
    out = assign_backstop_cds(gene, [aln])
    assert out.transcripts[0].cds is None


def test_assign_backstop_no_alignments_at_all():
    gene = backstop_gene([(1000, 1200)], "+")
    out = assign_backstop_cds(gene, [])
    assert out.transcripts[0].cds is None


# ---------------------------------------------------------------------------
# assign_backstop_cds — TransDecoder path (mocked)
# ---------------------------------------------------------------------------

def _write_bed(out_dir, strand="+", thick=(0, 150)):
    bed = os.path.join(out_dir, "orf.bed")
    with open(bed, "w") as fh:
        fh.write(
            "tx\t0\t200\tORF\t0\t%s\t%d\t%d\t0\t1\t200\t0\n"
            % (strand, thick[0], thick[1])
        )
    return bed


def test_assign_backstop_transdecoder_path(monkeypatch):
    gene = backstop_gene([(1000, 1100), (1200, 1300)], "+")
    g = genome_with({1000: "ATG"})

    def fake_run_transdecoder(fasta, out_dir, transdecoder_bin_dir=None):
        return _write_bed(out_dir, "+", (0, 150))

    monkeypatch.setattr(cds_mod, "run_transdecoder", fake_run_transdecoder)
    out = assign_backstop_cds(gene, [], genome=g, transdecoder_bin_dir="/fake/bin")
    tx = out.transcripts[0]
    assert bounds(tx.cds) == [(1000, 1100), (1200, 1250)]
    assert tx.protein_id is None  # TransDecoder ORFs carry no homology accession


def test_assign_backstop_transdecoder_not_invoked_without_bindir(monkeypatch):
    gene = backstop_gene([(1000, 1100), (1200, 1300)], "+")
    g = genome_with({1000: "ATG"})
    called = []
    monkeypatch.setattr(
        cds_mod, "run_transdecoder",
        lambda *a, **k: called.append(1) or _write_bed(a[1]),
    )
    out = assign_backstop_cds(gene, [], genome=g, transdecoder_bin_dir=None)
    assert out.transcripts[0].cds is None
    assert called == []


# ---------------------------------------------------------------------------
# _parse_transdecoder_bed
# ---------------------------------------------------------------------------

def test_parse_transdecoder_bed_plus(tmp_path):
    bed = _write_bed(str(tmp_path), "+", (12, 99))
    assert _parse_transdecoder_bed(bed) == (12, 99)


def test_parse_transdecoder_bed_skips_minus(tmp_path):
    bed = _write_bed(str(tmp_path), "-", (12, 99))
    assert _parse_transdecoder_bed(bed) is None


def test_parse_transdecoder_bed_missing_file():
    assert _parse_transdecoder_bed("/nonexistent/orf.bed") is None


# ---------------------------------------------------------------------------
# cds_cross_check
# ---------------------------------------------------------------------------

def test_cross_check_agreement_returns_none():
    gene = mikado_gene([(1000, 1500)], "+", cds=[(1000, 1150)])
    aln = alignment([(1000, 1150)], "+", start=1000, end=1500)
    assert cds_cross_check(gene, [aln]) is None


def test_cross_check_disagreement_flags():
    gene = mikado_gene([(1000, 1500)], "+", cds=[(1000, 1150)])
    aln = alignment([(1300, 1450)], "+", start=1000, end=1500)
    flag = cds_cross_check(gene, [aln])
    assert flag is not None and flag.name == "CDS_DISAGREE"


def test_cross_check_minus_disagreement_flags():
    gene = mikado_gene([(1000, 1500)], "-", cds=[(1000, 1150)])
    aln = alignment([(1300, 1450)], "-", start=1000, end=1500)
    assert cds_cross_check(gene, [aln]).name == "CDS_DISAGREE"


def test_cross_check_backstop_gene_returns_none():
    gene = backstop_gene([(1000, 1500)], "+", cds=[(1000, 1150)])
    aln = alignment([(1300, 1450)], "+", start=1000, end=1500)
    assert cds_cross_check(gene, [aln]) is None


def test_cross_check_no_cds_returns_none():
    gene = mikado_gene([(1000, 1500)], "+", cds=None)
    aln = alignment([(1000, 1150)], "+", start=1000, end=1500)
    assert cds_cross_check(gene, [aln]) is None


def test_cross_check_no_alignment_returns_none():
    gene = mikado_gene([(1000, 1500)], "+", cds=[(1000, 1150)])
    assert cds_cross_check(gene, []) is None


# ---------------------------------------------------------------------------
# _include_genomic_stop — both strands
# ---------------------------------------------------------------------------

def test_include_stop_plus_extends_when_stop_present():
    # CDS [1000,1096) = 96bp (stop-exclusive); genome has TAA at [1096,1099)
    g = genome_with({1096: "TAA"})
    cds = segs([(1000, 1096)])
    ex = exons([(1000, 1200)])
    result = _include_genomic_stop(cds, ex, "+", g, "chr1")
    assert result is not None
    assert [(s.start, s.end) for s in result] == [(1000, 1099)]


def test_include_stop_plus_returns_none_when_no_stop():
    # No stop codon at the 3' boundary
    g = genome_with({1096: "GGG"})
    cds = segs([(1000, 1096)])
    ex = exons([(1000, 1200)])
    assert _include_genomic_stop(cds, ex, "+", g, "chr1") is None


def test_include_stop_minus_extends_when_stop_present():
    # CDS [2004,2100) = 96bp (stop-exclusive); 3' coding end = low genomic
    # genome[2001:2004] with "-" → RC should be a stop; plant "TTA" → RC = "TAA"
    g = genome_with({2001: "TTA"})
    cds = segs([(2004, 2100)])
    ex = exons([(2000, 2100)])
    result = _include_genomic_stop(cds, ex, "-", g, "chr1")
    assert result is not None
    assert [(s.start, s.end) for s in result] == [(2001, 2100)]


def test_include_stop_minus_returns_none_when_no_stop():
    g = genome_with({2001: "GGG"})
    cds = segs([(2004, 2100)])
    ex = exons([(2000, 2100)])
    assert _include_genomic_stop(cds, ex, "-", g, "chr1") is None


def test_include_stop_plus_respects_exon_boundary():
    # CDS fills the exon exactly: extension would overflow
    g = genome_with({1096: "TAA"})
    cds = segs([(1000, 1096)])
    ex = exons([(1000, 1096)])  # exon ends right at CDS end
    assert _include_genomic_stop(cds, ex, "+", g, "chr1") is None


def test_include_stop_minus_respects_exon_boundary():
    g = genome_with({2001: "TTA"})
    cds = segs([(2004, 2100)])
    ex = exons([(2004, 2100)])  # exon starts right at CDS start
    assert _include_genomic_stop(cds, ex, "-", g, "chr1") is None


def test_include_stop_plus_tag_tga():
    # All three stop codons work
    for codon in ("TAA", "TAG", "TGA"):
        g = genome_with({1096: codon})
        cds = segs([(1000, 1096)])
        ex = exons([(1000, 1200)])
        result = _include_genomic_stop(cds, ex, "+", g, "chr1")
        assert result is not None, f"failed for {codon}"
        assert result[-1].end == 1099


def test_include_stop_empty_cds_returns_none():
    g = genome_with()
    assert _include_genomic_stop([], exons([(1000, 1200)]), "+", g, "chr1") is None
