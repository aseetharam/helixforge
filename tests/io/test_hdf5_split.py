"""Tests for native split-Helixer HDF5 support (Phase 13 D1). Floor: 12.

Synthetic split fixtures (a metadata ``*_input.h5`` + a ``*_predictions.h5``)
carry the *same* values as the combined fixture, so ``from_helixer_outputs``
must return identical confidence numbers. ``detect_helixer_layout`` is checked
for all three layouts. The core confidence assertions are parametrized across
both layouts (combined + split) rather than forked. Reverse-strand covered.
"""

import h5py
import numpy as np
import pytest

from helixforge.io.hdf5 import (
    LAYOUT_COMBINED,
    LAYOUT_SPLIT_METADATA,
    LAYOUT_SPLIT_PREDICTIONS,
    HDF5ConfidenceReader,
    detect_helixer_layout,
)


# --- layout detection (all three cases) ---


def test_detect_combined(hdf5_path):
    assert detect_helixer_layout(hdf5_path) == LAYOUT_COMBINED


def test_detect_split_metadata(hdf5_split_paths):
    inp, _pred = hdf5_split_paths
    assert detect_helixer_layout(inp) == LAYOUT_SPLIT_METADATA


def test_detect_split_predictions(hdf5_split_paths):
    _inp, pred = hdf5_split_paths
    assert detect_helixer_layout(pred) == LAYOUT_SPLIT_PREDICTIONS


def test_detect_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        detect_helixer_layout("/no/such/file.h5")


def test_detect_malformed_raises(tmp_path):
    p = tmp_path / "bad.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("seqids", data=np.array([b"chr1"]))
    with pytest.raises(ValueError):
        detect_helixer_layout(str(p))


# --- construction from a split pair ---


def test_from_helixer_outputs_explicit(hdf5_split_paths):
    inp, pred = hdf5_split_paths
    with HDF5ConfidenceReader.from_helixer_outputs(inp, pred) as r:
        assert r.seqids == ["chr1", "chr2"]


def test_split_sibling_inference(hdf5_split_paths):
    # Passing only the *_input.h5 must auto-resolve the *_predictions.h5 sibling.
    inp, _pred = hdf5_split_paths
    with HDF5ConfidenceReader(inp) as r:
        assert r.get_region_confidence("chr1", 10, 20) == pytest.approx(0.90, abs=1e-5)


def test_split_missing_predictions_raises(tmp_path):
    inp = tmp_path / "Lonely_input.h5"
    with h5py.File(inp, "w") as f:
        g = f.create_group("data")
        g.create_dataset("seqids", data=np.array([b"chr1"]))
        g.create_dataset("start_ends", data=np.array([[0, 10]], dtype=np.int64))
    # No sibling *_predictions.h5 on disk → clear error.
    with pytest.raises(FileNotFoundError):
        HDF5ConfidenceReader(str(inp))


def test_predictions_half_as_primary_raises(hdf5_split_paths):
    # Handing the predictions half in as the primary path is a usage error.
    _inp, pred = hdf5_split_paths
    with pytest.raises(ValueError):
        HDF5ConfidenceReader(pred)


def test_split_double_close_idempotent(hdf5_split_paths):
    inp, pred = hdf5_split_paths
    r = HDF5ConfidenceReader.from_helixer_outputs(inp, pred)
    r.close()
    r.close()  # both handles closed; idempotent


# --- combined vs split must agree, value-for-value (parametrized) ---


@pytest.fixture
def reader(request):
    """Yield a reader over the requested layout: 'combined' or 'split'."""
    if request.param == "combined":
        r = HDF5ConfidenceReader(request.getfixturevalue("hdf5_path"))
    else:
        inp, pred = request.getfixturevalue("hdf5_split_paths")
        r = HDF5ConfidenceReader.from_helixer_outputs(inp, pred)
    yield r
    r.close()


@pytest.mark.parametrize("reader", ["combined", "split"], indirect=True)
def test_region_confidence_cds(reader):
    # chr1 [10,20) is CDS: max(CDS=0.9, UTR=0.05) = 0.9
    assert reader.get_region_confidence("chr1", 10, 20) == pytest.approx(0.90, abs=1e-5)


@pytest.mark.parametrize("reader", ["combined", "split"], indirect=True)
def test_region_confidence_intergenic(reader):
    assert reader.get_region_confidence("chr1", 0, 10) == pytest.approx(0.05, abs=1e-5)


@pytest.mark.parametrize("reader", ["combined", "split"], indirect=True)
def test_intron_score(reader):
    assert reader.get_intron_score("chr1", 20, 40) == pytest.approx(0.85, abs=1e-5)


@pytest.mark.parametrize("reader", ["combined", "split"], indirect=True)
def test_exon_confidence_length_weighted(reader):
    # exon [10,20) CDS=0.9 and [40,50) UTR=0.85 -> (0.9*10 + 0.85*10)/20 = 0.875
    assert reader.get_exon_confidence("chr1", [(10, 20), (40, 50)]) == pytest.approx(
        0.875, abs=1e-5
    )


@pytest.mark.parametrize("reader", ["combined", "split"], indirect=True)
def test_cross_chunk_stitch(reader):
    preds = reader.get_per_base_predictions("chr1", 15, 25)
    assert preds.shape == (10, 4)
    assert preds[0, 2] == pytest.approx(0.90, abs=1e-5)  # CDS
    assert preds[5, 3] == pytest.approx(0.85, abs=1e-5)  # intron


@pytest.mark.parametrize("reader", ["combined", "split"], indirect=True)
def test_chr2_region(reader):
    assert reader.get_region_confidence("chr2", 20, 40) == pytest.approx(0.90, abs=1e-5)


# --- reverse strand on the split layout ---


def test_split_reverse_chunk(hdf5_reverse_split_paths):
    inp, pred = hdf5_reverse_split_paths
    with HDF5ConfidenceReader.from_helixer_outputs(inp, pred) as r:
        preds = r.get_per_base_predictions("chrR", 0, 5)
        assert preds.shape == (5, 4)
        assert np.allclose(preds[:, 2], 0.90, atol=1e-5)  # CDS at genomic 0-4
        tail = r.get_per_base_predictions("chrR", 5, 10)
        assert np.allclose(tail[:, 0], 0.85, atol=1e-5)  # intergenic at 5-9


def test_split_matches_combined_exactly(hdf5_path, hdf5_split_paths):
    inp, pred = hdf5_split_paths
    with HDF5ConfidenceReader(hdf5_path) as c, \
            HDF5ConfidenceReader.from_helixer_outputs(inp, pred) as s:
        a = c.get_per_base_predictions("chr1", 0, 64)
        b = s.get_per_base_predictions("chr1", 0, 64)
        assert np.array_equal(a, b)
