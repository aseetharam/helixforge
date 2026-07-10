"""protein-AED + RNA-AED + parallelism for the standalone ``evidence`` scorer.

Concrete literal coordinates, both strands (CLAUDE.md §12). Covers: v1's exact
RNA-AED formula on known ratios; the protein-AED component ratios + value from a
synthetic miniprot alignment (both strands); each AED column present only when
its evidence is supplied; the per-gene rollup picking the best transcript; and
``--threads N`` producing byte-identical output to serial.
"""

import math

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from helixforge.reconcile.models import CDSSegment, Exon, MiniprotAlignment
from helixforge.score.evidence import (
    PROT_W_CDS,
    PROT_W_PROT,
    PROT_W_STRUCT,
    RNA_W_BOUNDARY,
    RNA_W_COVERAGE,
    RNA_W_JUNCTION,
    protein_aed_from_ratios,
    rna_aed_from_ratios,
    rollup_genes,
    score_annotation,
    score_transcript_protein,
)


def _tx(seqid, strand, exons, cds=None, tid="t"):
    d = {
        "transcript_id": tid, "seqid": seqid, "strand": strand, "exons": exons,
        "start": min(e.start for e in exons), "end": max(e.end for e in exons),
    }
    if cds is not None:
        d["cds"] = cds
    return d


# ---------------------------------------------------------------------------
# RNA-AED — v1's exact formula on known ratios
# ---------------------------------------------------------------------------

def test_rna_aed_perfect_is_zero():
    assert rna_aed_from_ratios(1.0, 1.0, 1.0) == 0.0


def test_rna_aed_no_support_is_one():
    assert rna_aed_from_ratios(0.0, 0.0, 0.0) == 1.0


def test_rna_aed_matches_v1_weighting():
    # v1: 0.5*(1-jr) + 0.3*(1-cr) + 0.2*(1-br)
    jr, cr, br = 0.5, 0.25, 0.0
    expected = (
        RNA_W_JUNCTION * (1 - jr)
        + RNA_W_COVERAGE * (1 - cr)
        + RNA_W_BOUNDARY * (1 - br)
    )
    assert rna_aed_from_ratios(jr, cr, br) == pytest.approx(expected)
    assert rna_aed_from_ratios(jr, cr, br) == pytest.approx(0.25 + 0.225 + 0.2)


def test_rna_aed_clamped_to_unit_interval():
    assert rna_aed_from_ratios(-5.0, 0.0, 0.0) == 1.0      # >1 clamps to 1
    assert rna_aed_from_ratios(2.0, 2.0, 2.0) == 0.0       # negative clamps to 0


# ---------------------------------------------------------------------------
# protein-AED — component ratios + value, both strands
# ---------------------------------------------------------------------------

def test_protein_aed_perfect_match_plus():
    # model introns (150,250); alignment introns (150,250) -> struct 1.0;
    # exon-as-coding fully covered -> cds 1.0; query_coverage 0.8 -> prot 0.8
    exons = [Exon(100, 150), Exon(250, 300)]
    aln = MiniprotAlignment(
        protein_id="P1", seqid="chr1", start=100, end=300, strand="+",
        cds_segments=[CDSSegment(100, 150, 0), CDSSegment(250, 300, 0)],
        query_coverage=0.8, identity=0.95, score=500.0, rank=0,
    )
    m = score_transcript_protein(_tx("chr1", "+", exons), aln)
    assert m["protein_id"] == "P1"
    assert m["protein_struct_ratio"] == 1.0
    assert m["protein_cds_cov_ratio"] == 1.0
    assert m["protein_prot_cov_ratio"] == 0.8
    expected = PROT_W_STRUCT * 0 + PROT_W_CDS * 0 + PROT_W_PROT * (1 - 0.8)
    assert m["protein_aed"] == pytest.approx(expected)
    assert m["protein_aed"] == pytest.approx(0.04)


def test_protein_aed_perfect_match_minus():
    exons = [Exon(100, 150), Exon(250, 300)]
    aln = MiniprotAlignment(
        protein_id="P2", seqid="chr2", start=100, end=300, strand="-",
        cds_segments=[CDSSegment(100, 150, 0), CDSSegment(250, 300, 0)],
        query_coverage=1.0, identity=0.95, score=500.0, rank=0,
    )
    m = score_transcript_protein(_tx("chr2", "-", exons), aln)
    assert m["protein_struct_ratio"] == 1.0
    assert m["protein_cds_cov_ratio"] == 1.0
    assert m["protein_prot_cov_ratio"] == 1.0
    assert m["protein_aed"] == 0.0


def test_protein_aed_structure_mismatch():
    # alignment intron (150,260) != model intron (150,250) -> struct 0;
    # alignment CDS covers [100,150)+[260,300)=90/100 of model coding -> cds 0.9
    exons = [Exon(100, 150), Exon(250, 300)]
    aln = MiniprotAlignment(
        protein_id="P3", seqid="chr1", start=100, end=300, strand="+",
        cds_segments=[CDSSegment(100, 150, 0), CDSSegment(260, 300, 0)],
        query_coverage=1.0, identity=0.9, score=400.0, rank=0,
    )
    m = score_transcript_protein(_tx("chr1", "+", exons), aln)
    assert m["protein_struct_ratio"] == 0.0
    assert m["protein_cds_cov_ratio"] == pytest.approx(0.9)
    assert m["protein_prot_cov_ratio"] == 1.0
    expected = PROT_W_STRUCT * 1 + PROT_W_CDS * (1 - 0.9) + PROT_W_PROT * 0
    assert m["protein_aed"] == pytest.approx(expected)


def test_protein_aed_uses_model_cds_when_present():
    # CDS narrower than exons: coding intervals come from CDS, not exons.
    exons = [Exon(100, 150), Exon(250, 300)]
    cds = [CDSSegment(110, 150, 0), CDSSegment(250, 290, 0)]  # intron (150,250)
    aln = MiniprotAlignment(
        protein_id="P4", seqid="chr1", start=110, end=290, strand="+",
        cds_segments=[CDSSegment(110, 150, 0), CDSSegment(250, 290, 0)],
        query_coverage=1.0, identity=0.95, score=500.0, rank=0,
    )
    m = score_transcript_protein(_tx("chr1", "+", exons, cds=cds), aln)
    assert m["protein_struct_ratio"] == 1.0
    assert m["protein_cds_cov_ratio"] == 1.0


def test_protein_aed_single_exon_struct_ratio_is_blank():
    # Single-exon model: no intron chain to compare. struct_ratio is None (blank),
    # NOT 1.0 — "no introns to disagree" must not read as "introns agree". The AED
    # renormalises over the two present components (cds + prot), both perfect -> 0.
    exons = [Exon(100, 200)]
    aln = MiniprotAlignment(
        protein_id="P5", seqid="chr1", start=100, end=200, strand="+",
        cds_segments=[CDSSegment(100, 200, 0)],
        query_coverage=1.0, identity=0.95, score=500.0, rank=0,
    )
    m = score_transcript_protein(_tx("chr1", "+", exons), aln)
    assert m["protein_struct_ratio"] is None
    assert m["protein_aed"] == 0.0


def test_protein_aed_single_exon_renormalises_over_present_components():
    # Single-exon model, partial protein coverage: struct dropped, AED renormalised
    # over cds (perfect) + prot (0.5). aed = (0.3*0 + 0.2*0.5)/(0.3+0.2) = 0.2.
    exons = [Exon(100, 200)]
    aln = MiniprotAlignment(
        protein_id="P6", seqid="chr5", start=100, end=200, strand="-",
        cds_segments=[CDSSegment(100, 200, 0)],
        query_coverage=0.5, identity=0.95, score=500.0, rank=0,
    )
    m = score_transcript_protein(_tx("chr5", "-", exons), aln)
    assert m["protein_struct_ratio"] is None
    assert m["protein_cds_cov_ratio"] == 1.0
    assert m["protein_prot_cov_ratio"] == 0.5
    expected = (PROT_W_CDS * 0.0 + PROT_W_PROT * 0.5) / (PROT_W_CDS + PROT_W_PROT)
    assert m["protein_aed"] == pytest.approx(expected)
    assert m["protein_aed"] == pytest.approx(0.2)


def test_protein_aed_formula_helper():
    assert protein_aed_from_ratios(1.0, 1.0, 1.0) == 0.0
    assert protein_aed_from_ratios(0.0, 0.0, 0.0) == 1.0


# ---------------------------------------------------------------------------
# score_annotation — protein columns from a synthetic miniprot GFF
# ---------------------------------------------------------------------------

def test_score_annotation_protein_columns_present(ev_gff3_path, ev_miniprot_gff):
    df = score_annotation(ev_gff3_path, miniprot_gff=ev_miniprot_gff)
    rows = {r["transcript_id"]: r for _, r in df.iterrows()}
    # gA.t1 (chr1 +) matches MP1 exactly
    assert rows["gA.t1"]["protein_id"] == "P1"
    assert rows["gA.t1"]["protein_struct_ratio"] == 1.0
    assert rows["gA.t1"]["protein_aed"] == 0.0
    # gM.t1 (chr2 -) matches MP2 exactly
    assert rows["gM.t1"]["protein_id"] == "P2"
    assert rows["gM.t1"]["protein_aed"] == 0.0
    # gA.t2 overlaps MP1 but its intron (160,250) disagrees -> struct 0, aed > 0
    assert rows["gA.t2"]["protein_struct_ratio"] == 0.0
    assert rows["gA.t2"]["protein_aed"] > 0.0


def test_score_annotation_protein_absent_when_no_overlap(ev_gff3_path, ev_miniprot_gff):
    df = score_annotation(ev_gff3_path, miniprot_gff=ev_miniprot_gff)
    rows = {r["transcript_id"]: r for _, r in df.iterrows()}
    # gB.t1 is on chr1 at [600,800) — no alignment overlaps it.
    assert pd.isna(rows["gB.t1"]["protein_aed"])
    assert pd.isna(rows["gB.t1"]["protein_id"])


def test_score_annotation_protein_column_na_without_protein_input(ev_gff3_path, ev_bam_path):
    df = score_annotation(ev_gff3_path, bam_paths=[ev_bam_path])
    assert df["protein_aed"].isna().all()
    assert df["protein_id"].isna().all()


def test_score_annotation_rna_aed_present_with_bam(ev_gff3_path, ev_bam_path, ev_star_sj_path):
    df = score_annotation(
        ev_gff3_path, bam_paths=[ev_bam_path], star_sj_paths=[ev_star_sj_path]
    )
    rows = {r["transcript_id"]: r for _, r in df.iterrows()}
    # gA.t1: junction supported, both exons covered (depth 5 >= 5) -> rna_aed 0
    assert rows["gA.t1"]["rna_junction_ratio"] == 1.0
    assert rows["gA.t1"]["rna_coverage_ratio"] == 1.0
    assert rows["gA.t1"]["rna_boundary_ratio"] == 1.0
    assert rows["gA.t1"]["rna_aed"] == 0.0
    # minus strand model: junction supported (jr=1) but coverage is only 3 reads,
    # below the default min_exon_coverage=5, so cr=br=0 -> rna_aed = 0.3 + 0.2 = 0.5.
    assert rows["gM.t1"]["rna_junction_ratio"] == 1.0
    assert rows["gM.t1"]["rna_coverage_ratio"] == 0.0
    assert rows["gM.t1"]["rna_aed"] == pytest.approx(0.5)


def test_score_annotation_rna_aed_na_without_bam(ev_gff3_path, ev_star_sj_path):
    # SJ only -> junctions but no coverage -> rna_aed is neutral (NA), not a penalty.
    df = score_annotation(ev_gff3_path, star_sj_paths=[ev_star_sj_path])
    assert df["rna_aed"].isna().all()
    assert df["rna_coverage_ratio"].isna().all()


# ---------------------------------------------------------------------------
# per-gene rollup = best transcript
# ---------------------------------------------------------------------------

def test_rollup_picks_best_transcript(ev_gff3_path, ev_bam_path, ev_star_sj_path):
    df = score_annotation(
        ev_gff3_path, bam_paths=[ev_bam_path], star_sj_paths=[ev_star_sj_path]
    )
    gene_df = rollup_genes(df)
    # one row per gene
    assert set(gene_df["gene_id"]) == {"gA", "gB", "gM", "gS"}
    # gene gA has gA.t1 (rna_aed 0, supported) and gA.t2 (contradicted, rna_aed > 0)
    gA = gene_df[gene_df["gene_id"] == "gA"].iloc[0]
    assert gA["transcript_id"] == "gA.t1"
    assert gA["rna_aed"] == 0.0


def test_rollup_empty():
    from helixforge.score.evidence import TSV_COLUMNS
    out = rollup_genes(pd.DataFrame(columns=TSV_COLUMNS))
    assert list(out.columns) == TSV_COLUMNS
    assert len(out) == 0


# ---------------------------------------------------------------------------
# determinism: serial == --threads N
# ---------------------------------------------------------------------------

def test_serial_equals_threaded(ev_gff3_path, ev_bam_path, ev_bam_path2,
                                ev_star_sj_path, ev_stringtie_gtfs, ev_miniprot_gff):
    kwargs = dict(
        bam_paths=[ev_bam_path, ev_bam_path2],
        star_sj_paths=[ev_star_sj_path],
        stringtie_gtfs=ev_stringtie_gtfs,
        miniprot_gff=ev_miniprot_gff,
    )
    serial = score_annotation(ev_gff3_path, threads=1, **kwargs)
    threaded = score_annotation(ev_gff3_path, threads=3, **kwargs)
    assert_frame_equal(serial, threaded)


def test_threaded_extraction_over_multiple_bams(ev_gff3_path, ev_bam_path, ev_bam_path2):
    # chr1 + junction (150,250): 5 (sample1) + 2 (sample2) -> read_count 7, 2 samples.
    from helixforge.score.evidence import collect_junctions
    serial = collect_junctions(bam_paths=[ev_bam_path, ev_bam_path2], threads=1)
    parallel = collect_junctions(bam_paths=[ev_bam_path, ev_bam_path2], threads=2)
    by_key = {(j.seqid, j.donor, j.acceptor): j for j in parallel}
    assert by_key[("chr1", 150, 250)].read_count == 7
    assert by_key[("chr1", 150, 250)].samples == 2
    # identical to serial
    assert [(j.seqid, j.donor, j.acceptor, j.read_count, j.samples) for j in serial] == \
           [(j.seqid, j.donor, j.acceptor, j.read_count, j.samples) for j in parallel]
