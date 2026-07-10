"""Phase 31 D3 — genome-level annotation report.

Concrete literal coordinates / values only (CLAUDE.md §12); both strands. The
report is built from a synthetic gene set and asserted section-by-section, then
serialised to JSON (valid + deterministic) and HTML (renders).
"""

from __future__ import annotations

import json

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.stats.report import (
    build_report,
    write_report_html,
    write_report_json,
)


def _tx(tid, strand="+", cds=True, tpm=5.0, jsf=1.0, exons=None, biotype="protein_coding"):
    exons = exons or [Exon(1000, 1200), Exon(1300, 1500), Exon(1600, 1800)]
    # CDS leaves UTR (starts after exon start, ends before exon end).
    cdsseg = (
        [CDSSegment(1050, 1200, 0), CDSSegment(1300, 1500, 0), CDSSegment(1600, 1649, 1)]
        if cds else None
    )
    return TranscriptCandidate(
        transcript_id=tid, locus_id=tid.rsplit(".", 1)[0], source="mikado",
        seqid="chr1", start=exons[0].start, end=exons[-1].end, strand=strand,
        exons=exons, cds=cdsseg, tpm=tpm, junction_support_fraction=jsf,
        confidence=0.9, combined_score=10.0, is_primary=True, biotype=biotype,
    )


def _gene(gid, strand="+", cds=True, tpm=5.0, jsf=1.0, status="EXPRESSED",
          tier=1, biotype="protein_coding", txs=None):
    txs = txs or [_tx(f"{gid}.1", strand=strand, cds=cds, tpm=tpm, jsf=jsf, biotype=biotype)]
    return ReconciledGene(
        gene_id=gid, seqid="chr1", start=txs[0].start, end=txs[0].end, strand=strand,
        tier=tier, transcripts=txs, primary_transcript_id=txs[0].transcript_id,
        classification=LocusClassification(locus_id=gid, status=status, max_tpm=tpm),
        origin="mikado_1to1", biotype=biotype,
    )


@pytest.fixture
def gene_set():
    # 2 coding (one each strand, EXPRESSED), 1 lncRNA (minus, SILENT, no CDS).
    return [
        _gene("HFG_00001", strand="+"),
        _gene("HFG_00002", strand="-"),
        _gene("HFG_00003", strand="-", cds=False, tpm=0.0, jsf=0.0,
              status="SILENT", tier=3, biotype="lncRNA"),
    ]


def test_report_structure_counts(gene_set):
    r = build_report(gene_set)
    s = r["structure"]
    assert s["gene_count"] == 3
    assert s["transcript_count"] == 3
    assert s["coding_gene_count"] == 2
    assert s["exon_count"] == 9          # 3 transcripts × 3 exons
    assert s["multi_exon_genes"] == 3


def test_report_biotype_and_tier_sections(gene_set):
    r = build_report(gene_set)
    assert r["biotype"] == {"lncRNA": 1, "protein_coding": 2}
    assert r["tier"] == {"1": 2, "3": 1}
    assert r["origin"] == {"mikado_1to1": 3}


def test_report_classification_fractions(gene_set):
    r = build_report(gene_set)
    cls = r["classification"]
    assert cls["counts"] == {"EXPRESSED": 2, "SILENT": 1}
    assert cls["fractions"]["EXPRESSED"] == pytest.approx(2 / 3)
    assert cls["fractions"]["SILENT"] == pytest.approx(1 / 3)


def test_report_pct_with_utr_and_complete(gene_set):
    r = build_report(gene_set)
    s = r["structure"]
    # Both coding genes have UTR (CDS strictly inside the exon span) → 100%.
    assert s["pct_with_utr"] == pytest.approx(100.0)


def test_report_completeness_and_runstats_passthrough(gene_set):
    comp = {"busco": {"complete": 97.5}, "compleasm": {"complete": 98.1}}
    rs = {"merges_accepted": 5, "splits": 2}
    r = build_report(gene_set, completeness=comp, run_stats=rs)
    assert r["completeness"]["busco"]["complete"] == 97.5
    assert r["run_stats"]["merges_accepted"] == 5


def test_report_support_section(gene_set):
    r = build_report(gene_set, bam_stats={"mapped": 90, "unmapped": 10},
                     tpm_threshold=0.5)
    sup = r["support"]
    assert sup["mapping_rate"] == pytest.approx(0.9)
    # 2 of 3 genes are junction-supported (the lncRNA has jsf 0.0).
    assert sup["frac_genes_junction_supported"] == pytest.approx(2 / 3)
    # 2 of 3 genes have tpm >= 0.5.
    assert sup["frac_genes_tpm_pass"] == pytest.approx(2 / 3)


def test_report_multiqc_block(gene_set):
    r = build_report(gene_set, bam_stats={"mapped": 90, "unmapped": 10})
    mq = r["multiqc"]
    assert mq["plot_type"] == "generalstats"
    assert mq["data"]["HelixForge"]["genes"] == 3
    assert mq["data"]["HelixForge"]["coding_genes"] == 2


def test_report_json_valid_and_deterministic(tmp_path, gene_set):
    r = build_report(gene_set, bam_stats={"mapped": 90, "unmapped": 10})
    p1 = tmp_path / "r1.json"
    p2 = tmp_path / "r2.json"
    write_report_json(r, p1)
    write_report_json(build_report(gene_set, bam_stats={"mapped": 90, "unmapped": 10}), p2)
    # valid JSON + byte-identical across runs (sorted keys, NaN→null).
    assert json.loads(p1.read_text())["structure"]["gene_count"] == 3
    assert p1.read_text() == p2.read_text()


def test_report_html_renders(tmp_path, gene_set):
    r = build_report(gene_set, bam_stats={"mapped": 90, "unmapped": 10},
                     completeness={"busco": {"complete": 97.5}},
                     run_stats={"merges_accepted": 5})
    out = tmp_path / "report.html"
    write_report_html(r, out)
    text = out.read_text()
    assert "<html" in text and "HelixForge annotation report" in text
    assert "Biotype" in text and "Completeness" in text and "Evidence support" in text
