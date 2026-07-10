"""Tests for HDF5ConfidenceReader (Phase 1). Floor: 17.

Known high-CDS / high-intron / intergenic regions; cross-chunk stitching;
value ranges in [0,1]; context manager.
"""

import h5py
import numpy as np
import pytest

from helixforge.io.hdf5 import HDF5ConfidenceReader


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        HDF5ConfidenceReader("/no/such/file.h5")


def test_missing_dataset_raises(tmp_path):
    p = tmp_path / "bad.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("seqids", data=np.array([b"chr1"]))
    with pytest.raises(ValueError):
        HDF5ConfidenceReader(str(p))


def test_seqids_property(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.seqids == ["chr1", "chr2"]


def test_region_confidence_high_cds(hdf5_path):
    # chr1 [10,20) is CDS: max(CDS=0.9, UTR=0.05) = 0.9
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_region_confidence("chr1", 10, 20) == pytest.approx(0.90, abs=1e-5)


def test_region_confidence_intergenic_low(hdf5_path):
    # chr1 [0,10) intergenic: max(CDS=0.05, UTR=0.05) = 0.05
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_region_confidence("chr1", 0, 10) == pytest.approx(0.05, abs=1e-5)


def test_region_confidence_utr(hdf5_path):
    # chr1 [40,50) UTR: max(CDS=0.05, UTR=0.85) = 0.85
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_region_confidence("chr1", 40, 50) == pytest.approx(0.85, abs=1e-5)


def test_intron_score_high(hdf5_path):
    # chr1 [20,40) intron channel = 0.85
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_intron_score("chr1", 20, 40) == pytest.approx(0.85, abs=1e-5)


def test_intron_score_low_in_cds(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_intron_score("chr1", 10, 20) == pytest.approx(0.0, abs=1e-5)


def test_per_base_predictions_shape(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        preds = r.get_per_base_predictions("chr1", 10, 20)
        assert preds.shape == (10, 4)


def test_per_base_cross_chunk_stitch(hdf5_path):
    # [15,25) spans chunk0 (0-20) and chunk1 (20-40)
    with HDF5ConfidenceReader(hdf5_path) as r:
        preds = r.get_per_base_predictions("chr1", 15, 25)
        assert preds.shape == (10, 4)
        # positions 15-19 are CDS (channel 2 high)
        assert preds[0, 2] == pytest.approx(0.90, abs=1e-5)
        # positions 20-24 are intron (channel 3 high)
        assert preds[5, 3] == pytest.approx(0.85, abs=1e-5)


def test_per_base_last_partial_chunk(hdf5_path):
    # chr1 [60,64) is the padded final chunk; only 4 real positions
    with HDF5ConfidenceReader(hdf5_path) as r:
        preds = r.get_per_base_predictions("chr1", 60, 64)
        assert preds.shape == (4, 4)


def test_values_in_unit_range(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        preds = r.get_per_base_predictions("chr1", 0, 64)
        assert preds.min() >= 0.0
        assert preds.max() <= 1.0


def test_exon_confidence_single(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_exon_confidence("chr1", [(10, 20)]) == pytest.approx(0.90, abs=1e-5)


def test_exon_confidence_length_weighted(hdf5_path):
    # exon [10,20) CDS=0.9 and [40,50) UTR=0.85 -> (0.9*10 + 0.85*10)/20 = 0.875
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_exon_confidence("chr1", [(10, 20), (40, 50)]) == pytest.approx(0.875, abs=1e-5)


def test_chr2_region_confidence(hdf5_path):
    # chr2 [20,40) is CDS
    with HDF5ConfidenceReader(hdf5_path) as r:
        assert r.get_region_confidence("chr2", 20, 40) == pytest.approx(0.90, abs=1e-5)


def test_unknown_seqid_raises_key_error(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        with pytest.raises(KeyError):
            r.get_per_base_predictions("chrX", 0, 10)


def test_out_of_bounds_raises_index_error(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        with pytest.raises(IndexError):
            r.get_per_base_predictions("chr1", 0, 100)


def test_end_le_start_raises_value_error(hdf5_path):
    with HDF5ConfidenceReader(hdf5_path) as r:
        with pytest.raises(ValueError):
            r.get_per_base_predictions("chr1", 20, 20)


def test_context_manager_double_close(hdf5_path):
    r = HDF5ConfidenceReader(hdf5_path)
    r.close()
    r.close()  # idempotent


def test_reverse_chunk_orientation(hdf5_reverse_path):
    # chrR is a reverse chunk; CDS placed at genomic 0-4
    with HDF5ConfidenceReader(hdf5_reverse_path) as r:
        preds = r.get_per_base_predictions("chrR", 0, 5)
        assert preds.shape == (5, 4)
        assert np.allclose(preds[:, 2], 0.90, atol=1e-5)
        # genomic 5-9 is intergenic
        tail = r.get_per_base_predictions("chrR", 5, 10)
        assert np.allclose(tail[:, 0], 0.85, atol=1e-5)


# --- Phase 21 D1: the single-read scorer slices ONE block bit-identically to the
#     per-exon get_exon_confidence / get_intron_score / get_region_confidence reads.

from helixforge.mikado.emit_external import _score_components  # noqa: E402
from helixforge.reconcile.models import Exon  # noqa: E402


def _reference_components(exons, seqid, r, exon_weight=0.7, intron_high_cutoff=0.5):
    """Recompute the metric components via the per-window reader calls."""
    exon_conf = max(0.0, min(1.0, r.get_exon_confidence(seqid, exons)))
    ordered = sorted(exons, key=lambda e: e.start)
    introns = [(ordered[i].end, ordered[i + 1].start) for i in range(len(ordered) - 1)]
    if introns:
        high = sum(
            1 for s, e in introns
            if r.get_intron_score(seqid, s, e) >= intron_high_cutoff
        )
        frac = high / len(introns)
        support = max(0.0, min(1.0, exon_weight * exon_conf + (1 - exon_weight) * frac))
    else:
        frac = 0.0
        support = exon_conf
    lo = min(e.start for e in exons)
    hi = max(e.end for e in exons)
    locus_conf = max(0.0, min(1.0, r.get_region_confidence(seqid, lo, hi)))
    return support, locus_conf, exon_conf, frac, len(introns)


def test_single_read_components_multi_exon_bit_identical(hdf5_path):
    exons = [Exon(20, 30), Exon(40, 50)]  # chr1: intron (30,40) is in the INTRON band
    with HDF5ConfidenceReader(hdf5_path) as r:
        c = _score_components(exons, "chr1", r, 0.7, 0.5)
        s, lc, ec, fr, ni = _reference_components(exons, "chr1", r)
        assert (c["helixer_support"], c["helixer_locus_conf"]) == (s, lc)
        assert (c["exon_confidence"], c["intron_coincidence_fraction"]) == (ec, fr)
        assert c["num_introns"] == ni == 1


def test_single_read_components_single_exon_bit_identical(hdf5_path):
    exons = [Exon(0, 20)]  # chr1: single exon, no introns
    with HDF5ConfidenceReader(hdf5_path) as r:
        c = _score_components(exons, "chr1", r, 0.7, 0.5)
        s, lc, ec, fr, ni = _reference_components(exons, "chr1", r)
        assert (c["helixer_support"], c["helixer_locus_conf"], c["num_introns"]) == (s, lc, 0)
        assert c["helixer_support"] == c["exon_confidence"]


def test_single_read_components_reverse_chunk_bit_identical(hdf5_reverse_path):
    # chrR is a reverse-orientation chunk; exons + intron all inside genomic [0,10).
    exons = [Exon(0, 4), Exon(6, 9)]  # intron (4,6)
    with HDF5ConfidenceReader(hdf5_reverse_path) as r:
        c = _score_components(exons, "chrR", r, 0.7, 0.5)
        s, lc, ec, fr, ni = _reference_components(exons, "chrR", r)
        assert (c["helixer_support"], c["helixer_locus_conf"]) == (s, lc)
        assert (c["exon_confidence"], c["intron_coincidence_fraction"]) == (ec, fr)
        assert c["num_introns"] == ni == 1
