"""Tests for ``scripts/validate_evidence.py`` — the evidence-TSV invariant checker.

Concrete literal coordinates / values, both strands (CLAUDE.md §12). Covers: a
correct synthetic TSV passing cleanly; deliberately corrupted TSVs each tripping
the right rule (out-of-range ratio, broken AED reconstruction, wrong blank,
single-exon leaking intron_recall, single-exon struct_ratio populated); the
single-exon boundary/rna_aed staying defined; struct_ratio blank on a zero-intron
gene reconstructing the renormalised protein AED; and the summary-mean denominator
audit reporting n separately for RNA (all) vs protein (hit subset).
"""

import importlib.util
import math
from pathlib import Path

import pandas as pd
import pytest

from helixforge.score.evidence import (
    PROT_W_CDS,
    PROT_W_PROT,
    PROT_W_STRUCT,
    RNA_W_BOUNDARY,
    RNA_W_COVERAGE,
    RNA_W_JUNCTION,
    TSV_COLUMNS,
)

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "validate_evidence.py"


def _load():
    spec = importlib.util.spec_from_file_location("validate_evidence", _SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


VE = _load()


# ---------------------------------------------------------------------------
# Synthetic rows: one multi-exon (+), one multi-exon (-), one single-exon (+),
# one single-exon (-) with a protein hit. All values internally consistent.
# ---------------------------------------------------------------------------


def _rna_aed(jr, cr, br):
    return RNA_W_JUNCTION * (1 - jr) + RNA_W_COVERAGE * (1 - cr) + RNA_W_BOUNDARY * (1 - br)


def _prot_aed_full(s, c, p):
    return PROT_W_STRUCT * (1 - s) + PROT_W_CDS * (1 - c) + PROT_W_PROT * (1 - p)


def _prot_aed_renorm(c, p):
    return (PROT_W_CDS * (1 - c) + PROT_W_PROT * (1 - p)) / (PROT_W_CDS + PROT_W_PROT)


def _good_rows():
    nan = math.nan
    rows = [
        # multi-exon +, 2 introns both supported: jsf=prec=jr=1.0, full intron metrics
        {
            "gene_id": "g1", "transcript_id": "g1.1", "seqid": "chr1", "strand": "+",
            "start": 1000, "end": 2000, "num_exons": 3, "num_introns": 2,
            "supported": 2, "contradicted": 0, "novel_in_data": 0,
            "junction_support_fraction": 1.0, "intron_precision": 1.0,
            "intron_recall": 1.0, "intron_f1": 1.0, "mean_coverage": 50.0, "tpm": 10.0,
            "rna_aed": _rna_aed(1.0, 1.0, 1.0), "rna_junction_ratio": 1.0,
            "rna_coverage_ratio": 1.0, "rna_boundary_ratio": 1.0,
            "protein_id": "P1", "protein_aed": _prot_aed_full(1.0, 1.0, 0.8),
            "protein_struct_ratio": 1.0, "protein_cds_cov_ratio": 1.0,
            "protein_prot_cov_ratio": 0.8,
        },
        # multi-exon -, 4 introns 1 supported, no qualifying-junction recall=None ok? here
        # supported=1 so recall must be defined -> give it 0.5/f1.
        {
            "gene_id": "g2", "transcript_id": "g2.1", "seqid": "chr2", "strand": "-",
            "start": 3000, "end": 4000, "num_exons": 5, "num_introns": 4,
            "supported": 1, "contradicted": 1, "novel_in_data": 2,
            "junction_support_fraction": 0.25, "intron_precision": 0.25,
            "intron_recall": 0.5, "intron_f1": 2 * 0.25 * 0.5 / (0.25 + 0.5),
            "mean_coverage": 5.0, "tpm": nan,
            "rna_aed": _rna_aed(0.25, 1.0, 0.5), "rna_junction_ratio": 0.25,
            "rna_coverage_ratio": 1.0, "rna_boundary_ratio": 0.5,
            "protein_id": nan, "protein_aed": nan, "protein_struct_ratio": nan,
            "protein_cds_cov_ratio": nan, "protein_prot_cov_ratio": nan,
        },
        # single-exon +, no introns: all intron metrics blank; rna_junction_ratio=1.0;
        # boundary defined (1.0); no protein hit.
        {
            "gene_id": "g3", "transcript_id": "g3.1", "seqid": "chr1", "strand": "+",
            "start": 5000, "end": 5500, "num_exons": 1, "num_introns": 0,
            "supported": 0, "contradicted": 0, "novel_in_data": 0,
            "junction_support_fraction": nan, "intron_precision": nan,
            "intron_recall": nan, "intron_f1": nan, "mean_coverage": 30.0, "tpm": 2.0,
            "rna_aed": _rna_aed(1.0, 1.0, 1.0), "rna_junction_ratio": 1.0,
            "rna_coverage_ratio": 1.0, "rna_boundary_ratio": 1.0,
            "protein_id": nan, "protein_aed": nan, "protein_struct_ratio": nan,
            "protein_cds_cov_ratio": nan, "protein_prot_cov_ratio": nan,
        },
        # single-exon -, no introns, WITH a protein hit: struct_ratio blank,
        # protein_aed renormalised over cds+prot. boundary=0.0 (uncovered terminus).
        {
            "gene_id": "g4", "transcript_id": "g4.1", "seqid": "chr2", "strand": "-",
            "start": 7000, "end": 7600, "num_exons": 1, "num_introns": 0,
            "supported": 0, "contradicted": 0, "novel_in_data": 0,
            "junction_support_fraction": nan, "intron_precision": nan,
            "intron_recall": nan, "intron_f1": nan, "mean_coverage": 0.5, "tpm": nan,
            "rna_aed": _rna_aed(1.0, 0.0, 0.0), "rna_junction_ratio": 1.0,
            "rna_coverage_ratio": 0.0, "rna_boundary_ratio": 0.0,
            "protein_id": "P4", "protein_aed": _prot_aed_renorm(1.0, 0.5),
            "protein_struct_ratio": nan, "protein_cds_cov_ratio": 1.0,
            "protein_prot_cov_ratio": 0.5,
        },
    ]
    return [{c: r.get(c, nan) for c in TSV_COLUMNS} for r in rows]


@pytest.fixture
def good_rows():
    return _good_rows()


# ---------------------------------------------------------------------------
# Clean TSV passes
# ---------------------------------------------------------------------------


def test_clean_tsv_passes(good_rows):
    rep = VE.validate(good_rows)
    assert rep.n_categories() == 0, rep.violations


def test_clean_tsv_via_file(good_rows, tmp_path):
    path = tmp_path / "ev.tsv"
    pd.DataFrame(good_rows, columns=TSV_COLUMNS).to_csv(path, sep="\t", index=False)
    rows = VE._read_tsv(path)
    assert VE.validate(rows).n_categories() == 0


# ---------------------------------------------------------------------------
# Corruptions each trip exactly the intended rule
# ---------------------------------------------------------------------------


def test_out_of_range_ratio_flagged(good_rows):
    good_rows[0]["rna_coverage_ratio"] = 1.5  # > 1
    rep = VE.validate(good_rows)
    assert "ratio_out_of_range" in rep.violations
    # the broken value also breaks the AED reconstruction off the stored ratio? no —
    # we recompute aed FROM the (now 1.5) ratio, so reconstruction still matches; only
    # the range rule fires.
    assert any("g1.1" in m for m in rep.violations["ratio_out_of_range"])


def test_broken_rna_aed_reconstruction_flagged(good_rows):
    good_rows[0]["rna_aed"] = 0.42  # inconsistent with its ratios (all 1.0 -> 0.0)
    rep = VE.validate(good_rows)
    assert "rna_aed_reconstruction" in rep.violations
    assert any("g1.1" in m for m in rep.violations["rna_aed_reconstruction"])


def test_broken_protein_aed_reconstruction_flagged(good_rows):
    good_rows[3]["protein_aed"] = 0.99  # renorm case; stored value is wrong
    rep = VE.validate(good_rows)
    assert "protein_aed_reconstruction" in rep.violations


def test_wrong_blank_intron_metric_with_introns_flagged(good_rows):
    good_rows[0]["intron_precision"] = math.nan  # multi-exon must not be blank
    rep = VE.validate(good_rows)
    assert "intron_metric_blank_with_introns" in rep.violations


def test_single_exon_recall_leak_flagged(good_rows):
    # The historical bug: a single-exon model carrying intron_recall=0.0.
    good_rows[2]["intron_recall"] = 0.0
    rep = VE.validate(good_rows)
    assert "intron_metric_populated_no_introns" in rep.violations
    assert any("g3.1" in m for m in rep.violations["intron_metric_populated_no_introns"])


def test_single_exon_struct_ratio_populated_flagged(good_rows):
    # The degenerate struct_ratio=1.0 on a no-intron protein model.
    good_rows[3]["protein_struct_ratio"] = 1.0
    rep = VE.validate(good_rows)
    assert "struct_ratio_populated_no_introns" in rep.violations


def test_protein_metric_blank_with_hit_flagged(good_rows):
    good_rows[0]["protein_cds_cov_ratio"] = math.nan  # has hit -> must be populated
    rep = VE.validate(good_rows)
    assert "protein_metric_blank_with_hit" in rep.violations


def test_count_partition_violation_flagged(good_rows):
    good_rows[1]["novel_in_data"] = 99  # 1+1+99 != 4
    rep = VE.validate(good_rows)
    assert "count_partition" in rep.violations


def test_f1_reconstruction_violation_flagged(good_rows):
    good_rows[1]["intron_f1"] = 0.9  # != harmonic_mean(0.25, 0.5)
    rep = VE.validate(good_rows)
    assert "f1_reconstruction" in rep.violations


# ---------------------------------------------------------------------------
# Single-exon behaviour (Part B.1) + struct_ratio zero-intron (Part B.2)
# ---------------------------------------------------------------------------


def test_single_exon_boundary_and_rna_aed_defined(good_rows):
    # Both single-exon rows carry a defined boundary ratio + rna_aed (not blank),
    # and the AED reconstructs from the stored ratios.
    rep = VE.validate(good_rows)
    assert rep.n_categories() == 0
    g3, g4 = good_rows[2], good_rows[3]
    assert not VE._is_blank(g3["rna_boundary_ratio"])
    assert not VE._is_blank(g3["rna_aed"])
    assert g4["rna_boundary_ratio"] == 0.0  # uncovered terminus -> 0, still defined


def test_zero_intron_struct_ratio_blank_reconstructs_renormalised(good_rows):
    # g4 is single-exon with a hit: struct blank, AED = renorm(cds,prot).
    g4 = good_rows[3]
    assert VE._is_blank(g4["protein_struct_ratio"])
    expected = _prot_aed_renorm(1.0, 0.5)
    assert g4["protein_aed"] == pytest.approx(expected)
    assert VE.validate(good_rows).n_categories() == 0


# ---------------------------------------------------------------------------
# Summary-mean denominator audit (Part B.3)
# ---------------------------------------------------------------------------


def test_summary_denominators_reported_separately(good_rows):
    sd = VE._summary_denominators(good_rows)
    assert sd["n_all"] == 4
    assert sd["n_rna_aed"] == 4          # every row has an RNA AED
    assert sd["n_protein_aed"] == 2      # only g1.1, g4.1 have a protein hit
    assert sd["n_rna_aed"] != sd["n_protein_aed"]
    assert sd["mean_rna_aed"] is not None
    assert sd["mean_protein_aed"] is not None
