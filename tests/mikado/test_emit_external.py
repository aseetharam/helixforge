"""Tests for mikado/emit_external.py — the Helixer prior (Phase 4). Floor: 16.

Synthetic HDF5; metrics must land in [0,1]; the TSV writer must reject
out-of-range values. Both strands exercised.
"""

import pytest

from helixforge.io.hdf5 import HDF5ConfidenceReader
from helixforge.mikado.emit_external import (
    helixer_locus_conf,
    helixer_support,
    helixer_support_and_conf,
    normalize_tpm,
    write_external_scores_tsv,
)
from helixforge.reconcile.models import Exon

# Exon structures over the synthetic HDF5 (see conftest layout).
_TWO_EXON = [Exon(100, 150), Exon(200, 300)]          # intron (150,200) = high
_THREE_EXON = [Exon(100, 150), Exon(200, 300), Exon(450, 500)]  # +(300,450)=low
_SINGLE = [Exon(100, 150)]


@pytest.fixture
def reader(ext_h5_path):
    with HDF5ConfidenceReader(ext_h5_path) as r:
        yield r


# --- helixer_support ---

def test_support_default_weight(reader):
    # 0.7*exon_conf(0.9) + 0.3*intron_frac(1.0) = 0.93
    val = helixer_support(_TWO_EXON, "chr1", "+", reader)
    assert val == pytest.approx(0.93, abs=1e-5)


def test_support_exon_weight_one(reader):
    val = helixer_support(_TWO_EXON, "chr1", "+", reader, exon_weight=1.0)
    assert val == pytest.approx(0.90, abs=1e-5)


def test_support_exon_weight_zero(reader):
    val = helixer_support(_TWO_EXON, "chr1", "+", reader, exon_weight=0.0)
    # intron fraction only (the single intron is high) = 1.0
    assert val == pytest.approx(1.0, abs=1e-5)


def test_support_single_exon_is_exon_conf(reader):
    val = helixer_support(_SINGLE, "chr1", "+", reader)
    assert val == pytest.approx(0.90, abs=1e-5)


def test_support_three_exon_partial_intron_support(reader):
    # introns: (150,200) high, (300,450) low -> frac 0.5 -> 0.7*0.9+0.3*0.5=0.78
    val = helixer_support(_THREE_EXON, "chr1", "+", reader)
    assert val == pytest.approx(0.78, abs=1e-5)


def test_support_intron_cutoff_makes_all_low(reader):
    # cutoff above the high intron's 0.85 -> no introns count -> frac 0
    val = helixer_support(_TWO_EXON, "chr1", "+", reader, exon_weight=0.0,
                          intron_high_cutoff=0.99)
    assert val == pytest.approx(0.0, abs=1e-5)


def test_support_both_strands_equal(reader):
    plus = helixer_support(_TWO_EXON, "chr1", "+", reader)
    minus = helixer_support(_TWO_EXON, "chr1", "-", reader)
    assert plus == pytest.approx(minus)


def test_support_in_unit_range(reader):
    for exons in (_TWO_EXON, _THREE_EXON, _SINGLE):
        val = helixer_support(exons, "chr1", "+", reader)
        assert 0.0 <= val <= 1.0


# --- helixer_locus_conf ---

def test_locus_conf_value(reader):
    # region (100,300): (50*0.9 + 50*0.05 + 100*0.9)/200 = 0.6875
    val = helixer_locus_conf("chr1", 100, 300, reader)
    assert val == pytest.approx(0.6875, abs=1e-5)


def test_locus_conf_in_unit_range(reader):
    assert 0.0 <= helixer_locus_conf("chr1", 0, 512, reader) <= 1.0


def test_locus_conf_both_strands_same(reader):
    # span-based, strand-agnostic
    a = helixer_locus_conf("chr1", 100, 300, reader)
    assert 0.0 <= a <= 1.0


# --- normalize_tpm ---

def test_normalize_tpm_fraction():
    assert normalize_tpm(5.0, 10.0) == pytest.approx(0.5)


def test_normalize_tpm_caps_at_one():
    assert normalize_tpm(20.0, 10.0) == pytest.approx(1.0)


def test_normalize_tpm_zero_max():
    assert normalize_tpm(5.0, 0.0) == 0.0


# --- write_external_scores_tsv ---

def test_tsv_header_and_rows(tmp_path):
    rows = {
        "t1": {"helixer_support": 0.9, "helixer_locus_conf": 0.6},
        "t2": {"helixer_support": 0.1, "helixer_locus_conf": 0.2},
    }
    out = write_external_scores_tsv(rows, tmp_path / "ext.tsv")
    lines = out.read_text().splitlines()
    assert lines[0] == "tid\thelixer_support\thelixer_locus_conf"
    assert len(lines) == 3


def test_tsv_rejects_value_above_one(tmp_path):
    rows = {"t1": {"helixer_support": 1.5}}
    with pytest.raises(ValueError):
        write_external_scores_tsv(rows, tmp_path / "ext.tsv")


def test_tsv_rejects_negative_value(tmp_path):
    rows = {"t1": {"helixer_support": -0.1}}
    with pytest.raises(ValueError):
        write_external_scores_tsv(rows, tmp_path / "ext.tsv")


def test_tsv_accepts_boundary_values(tmp_path):
    rows = {"t1": {"m": 0.0}, "t2": {"m": 1.0}}
    out = write_external_scores_tsv(rows, tmp_path / "ext.tsv")
    assert len(out.read_text().splitlines()) == 3


def test_tsv_rejects_mismatched_metric_keys(tmp_path):
    rows = {"t1": {"a": 0.5, "b": 0.5}, "t2": {"a": 0.5}}
    with pytest.raises(ValueError):
        write_external_scores_tsv(rows, tmp_path / "ext.tsv")


def test_tsv_empty_rows(tmp_path):
    out = write_external_scores_tsv({}, tmp_path / "ext.tsv")
    assert out.read_text() == "tid\n"


# --- Phase 21 D1: single-read scoring == per-exon implementation bit-for-bit ---

def _per_exon_reference(exons, seqid, reader):
    """The pre-Phase-21 multi-read result: (support, locus_conf)."""
    support = helixer_support(exons, seqid, "+", reader)
    locus_conf = helixer_locus_conf(seqid, exons[0].start, exons[-1].end, reader)
    return support, locus_conf


def test_single_read_equals_per_exon_two_exon(reader):
    got = helixer_support_and_conf(_TWO_EXON, "chr1", "+", reader)
    ref = _per_exon_reference(_TWO_EXON, "chr1", reader)
    assert got == ref  # exact float equality (same bytes, same arithmetic)


def test_single_read_equals_per_exon_three_exon(reader):
    got = helixer_support_and_conf(_THREE_EXON, "chr1", "+", reader)
    ref = _per_exon_reference(_THREE_EXON, "chr1", reader)
    assert got == ref


def test_single_read_equals_per_exon_single_exon(reader):
    got = helixer_support_and_conf(_SINGLE, "chr1", "+", reader)
    ref = _per_exon_reference(_SINGLE, "chr1", reader)
    assert got == ref


def test_single_read_both_strands_equal(reader):
    plus = helixer_support_and_conf(_TWO_EXON, "chr1", "+", reader)
    minus = helixer_support_and_conf(_TWO_EXON, "chr1", "-", reader)
    assert plus == minus
    # and minus still equals the per-exon reference (strand-agnostic values)
    assert minus == _per_exon_reference(_TWO_EXON, "chr1", reader)


def test_single_read_custom_exon_weight_equals_per_exon(reader):
    got = helixer_support_and_conf(_THREE_EXON, "chr1", "+", reader, exon_weight=0.4)
    support = helixer_support(_THREE_EXON, "chr1", "+", reader, exon_weight=0.4)
    locus_conf = helixer_locus_conf("chr1", _THREE_EXON[0].start, _THREE_EXON[-1].end, reader)
    assert got == (support, locus_conf)
