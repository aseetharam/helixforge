"""Phase 18 D1 — miniprot backstop rescue (tier bump + BACKSTOP_RESCUED).

A silent Helixer-only gene that receives a valid projected CDS should not stay
Tier 4. A homology-backed (miniprot accession) CDS re-tiers it to Tier 1; a
CDS-only rescue (no homology) to Tier 2; a no-CDS backstop stays at Tier 3/4.
Every rescue is flagged BACKSTOP_RESCUED and keeps the gene's prior flags.

Concrete literal coordinates only (CLAUDE.md §12); both strands mandatory.
"""

import helixforge.reconcile.cds as cds_mod
from helixforge.qc.flags import HELIXER_ONLY, NO_EXPRESSION
from helixforge.reconcile.cds import _rescue_backstop_cds, assign_backstop_cds
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


def exons(bounds):
    return [Exon(s, e) for s, e in bounds]


def segs(bounds):
    return [CDSSegment(*b) if len(b) == 3 else CDSSegment(b[0], b[1], 0) for b in bounds]


def backstop_gene(exon_bounds, strand="+", status="SILENT", tier=4, flags=None):
    tx = TranscriptCandidate(
        transcript_id="HFG_00001.1",
        locus_id="HFG_00001",
        source="helixer",
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        exons=exons(exon_bounds),
        is_primary=True,
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
        classification=LocusClassification("HFG_00001", status),
        origin="helixer_backstop",
        flags=list(flags) if flags else [HELIXER_ONLY, NO_EXPRESSION],
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


def names(flags):
    return {f.name for f in flags}


def bounds(segments):
    return [(s.start, s.end) for s in segments]


# --- miniprot homology rescue → Tier 1, both strands ---

def test_silent_backstop_homology_rescue_plus_tier1():
    gene = backstop_gene([(1000, 1200), (1300, 1500)], "+", status="SILENT", tier=4)
    aln = alignment([(1000, 1150), (1300, 1450)], "+", start=1000, end=1500)
    out = assign_backstop_cds(gene, [aln])
    assert out.tier == 1
    assert bounds(out.transcripts[0].cds) == [(1000, 1150), (1300, 1450)]
    assert out.transcripts[0].protein_id == "sp|TEST"
    assert "BACKSTOP_RESCUED" in names(out.flags)


def test_silent_backstop_homology_rescue_minus_tier1():
    gene = backstop_gene([(1000, 1200), (1300, 1500)], "-", status="SILENT", tier=4)
    aln = alignment([(1000, 1150), (1300, 1450)], "-", start=1000, end=1500)
    out = assign_backstop_cds(gene, [aln])
    assert out.tier == 1
    assert out.transcripts[0].cds is not None
    assert "BACKSTOP_RESCUED" in names(out.flags)


def test_expressed_backstop_rescue_overrides_tier3():
    # An EXPRESSED backstop would be Tier 3 without a CDS; rescue lifts it to 1.
    gene = backstop_gene([(1000, 1500)], "+", status="EXPRESSED", tier=3)
    aln = alignment([(1000, 1300)], "+", start=1000, end=1500)
    out = assign_backstop_cds(gene, [aln])
    assert out.tier == 1
    assert "BACKSTOP_RESCUED" in names(out.flags)


def test_rescue_preserves_prior_flags():
    gene = backstop_gene([(1000, 1500)], "+", flags=[HELIXER_ONLY, NO_EXPRESSION])
    aln = alignment([(1000, 1300)], "+", start=1000, end=1500)
    out = assign_backstop_cds(gene, [aln])
    assert {"HELIXER_ONLY", "NO_EXPRESSION", "BACKSTOP_RESCUED"} <= names(out.flags)


# --- CDS-only rescue (no homology) → Tier 2 ---

def test_cds_only_rescue_no_homology_tier2_plus():
    gene = backstop_gene([(1000, 1200), (1300, 1500)], "+", status="SILENT", tier=4)
    out = _rescue_backstop_cds(gene, segs([(1000, 1150), (1300, 1450)]), protein_id=None)
    assert out.tier == 2
    assert out.transcripts[0].protein_id is None
    assert "BACKSTOP_RESCUED" in names(out.flags)


def test_cds_only_rescue_no_homology_tier2_minus():
    gene = backstop_gene([(1000, 1200), (1300, 1500)], "-", status="SILENT", tier=4)
    out = _rescue_backstop_cds(gene, segs([(1000, 1150), (1300, 1450)]), protein_id=None)
    assert out.tier == 2
    assert "BACKSTOP_RESCUED" in names(out.flags)


def test_transdecoder_path_rescues_tier2(monkeypatch):
    gene = backstop_gene([(1000, 1500)], "+", status="EXPRESSED", tier=3)
    genome = MockGenome({"chr1": "A" * CONTIG_LEN})
    # No miniprot hit; mock the TransDecoder ORF mapping to a valid CDS.
    monkeypatch.setattr(
        cds_mod, "_transdecoder_cds", lambda *a, **k: segs([(1000, 1300)])
    )
    out = assign_backstop_cds(gene, [], genome=genome, transdecoder_bin_dir="/bin")
    assert out.tier == 2
    assert "BACKSTOP_RESCUED" in names(out.flags)


# --- no rescue: stays Tier 3/4, no flag ---

def test_no_overlapping_alignment_stays_tier4():
    gene = backstop_gene([(1000, 1200)], "+", status="SILENT", tier=4)
    aln = alignment([(5000, 5150)], "+", start=5000, end=5200)
    out = assign_backstop_cds(gene, [aln])
    assert out.tier == 4
    assert out.transcripts[0].cds is None
    assert "BACKSTOP_RESCUED" not in names(out.flags)


def test_wrong_strand_alignment_stays_tier3_minus():
    gene = backstop_gene([(1000, 1500)], "-", status="EXPRESSED", tier=3)
    aln = alignment([(1000, 1300)], "+", start=1000, end=1500)  # plus vs minus gene
    out = assign_backstop_cds(gene, [aln])
    assert out.tier == 3
    assert out.transcripts[0].cds is None
    assert "BACKSTOP_RESCUED" not in names(out.flags)


def test_non_mod3_projection_rejected_stays_tier4():
    # A 149 bp projection (not mod-3) is rejected by project_cds_to_exons; no CDS,
    # tier unchanged, no flag.
    gene = backstop_gene([(1000, 1200)], "+", status="SILENT", tier=4)
    aln = alignment([(1000, 1149)], "+", start=1000, end=1200)
    out = assign_backstop_cds(gene, [aln])
    assert out.tier == 4
    assert out.transcripts[0].cds is None
    assert "BACKSTOP_RESCUED" not in names(out.flags)
