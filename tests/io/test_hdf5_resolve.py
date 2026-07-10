"""Tests for ``open_confidence_reader`` — the reconcile-side half resolver.

``reconcile`` is handed whichever Helixer HDF5 the user has (a combined file, the
``*_input.h5`` metadata half, or the bare ``*_predictions.h5`` softmax half) and
must resolve the missing partner exactly the way ``confidence`` does. These tests
lock that resolution: sibling auto-detection, explicit pairing (either order), the
combined passthrough, and a clear error when the predictions half is orphaned.
Both strands are covered via the reverse-orientation split fixture. Floor: 8.
"""

import h5py
import numpy as np
import pytest

from helixforge.io.hdf5 import (
    HDF5ConfidenceReader,
    _sibling_metadata,
    open_confidence_reader,
)


# --- sibling-name inference (reverse of _sibling_predictions) ---


def test_sibling_metadata_from_predictions():
    assert _sibling_metadata("/d/Zea_predictions.h5") == "/d/Zea_input.h5"


def test_sibling_metadata_non_conforming_name():
    assert _sibling_metadata("/d/something.h5") is None


# --- predictions half auto-resolves its *_input.h5 sibling ---


def test_predictions_half_resolves_sibling(hdf5_split_paths):
    # The exact failure from the bug report: a bare *_predictions.h5 with a
    # co-located *_input.h5 must open, not raise.
    _inp, pred = hdf5_split_paths
    with open_confidence_reader(pred) as r:
        assert r.seqids == ["chr1", "chr2"]
        # chr1 [10,20) is CDS: max(CDS=0.9, UTR=0.05) = 0.9
        assert r.get_region_confidence("chr1", 10, 20) == pytest.approx(0.90, abs=1e-5)


def test_predictions_half_resolves_sibling_minus_strand(hdf5_reverse_split_paths):
    # Reverse-orientation chunk: CDS at genomic 0-4 -> max(CDS) confidence.
    _inp, pred = hdf5_reverse_split_paths
    with open_confidence_reader(pred) as r:
        assert r.seqids == ["chrR"]
        assert r.get_region_confidence("chrR", 0, 5) == pytest.approx(0.90, abs=1e-5)


# --- both halves named explicitly (order-independent) ---


def test_explicit_pair_predictions_primary(hdf5_split_paths):
    inp, pred = hdf5_split_paths
    with open_confidence_reader(pred, input_h5=inp) as r:
        assert r.get_region_confidence("chr1", 10, 20) == pytest.approx(0.90, abs=1e-5)


def test_explicit_pair_order_independent(hdf5_split_paths):
    # Passing the metadata half as the primary and predictions via input_h5
    # still pairs correctly (resolver keys on detected layout, not argument order).
    inp, pred = hdf5_split_paths
    with open_confidence_reader(inp, input_h5=pred) as r:
        assert r.get_region_confidence("chr1", 10, 20) == pytest.approx(0.90, abs=1e-5)


# --- metadata half + combined passthrough still work ---


def test_metadata_half_primary_autofinds_predictions(hdf5_split_paths):
    inp, _pred = hdf5_split_paths
    with open_confidence_reader(inp) as r:
        assert r.get_region_confidence("chr1", 0, 10) == pytest.approx(0.05, abs=1e-5)


def test_combined_file_passthrough(hdf5_path):
    with open_confidence_reader(hdf5_path) as r:
        assert r.get_region_confidence("chr1", 10, 20) == pytest.approx(0.90, abs=1e-5)


# --- orphaned predictions half: clear, actionable error (not the raw ValueError) ---


def test_orphan_predictions_half_clear_error(tmp_path):
    pred = tmp_path / "Orphan_predictions.h5"
    with h5py.File(pred, "w") as f:
        f.create_dataset("predictions", data=np.zeros((1, 10, 4), dtype=np.float32))
    # No sibling *_input.h5 and none supplied -> FileNotFoundError naming the fix.
    with pytest.raises(FileNotFoundError, match="helixer-input-h5"):
        open_confidence_reader(str(pred))


def test_orphan_predictions_half_supplied_input_recovers(tmp_path):
    pred = tmp_path / "Orphan_predictions.h5"
    inp = tmp_path / "Elsewhere_input.h5"  # deliberately NOT the sibling name
    with h5py.File(pred, "w") as f:
        f.create_dataset(
            "predictions",
            data=np.tile(np.array([0.05, 0.05, 0.90, 0.0], dtype=np.float32), (1, 10, 1)),
        )
    with h5py.File(inp, "w") as f:
        g = f.create_group("data")
        g.create_dataset("seqids", data=np.array([b"chr1"]))
        g.create_dataset("start_ends", data=np.array([[0, 10]], dtype=np.int64))
    with open_confidence_reader(str(pred), input_h5=str(inp)) as r:
        assert r.get_region_confidence("chr1", 0, 10) == pytest.approx(0.90, abs=1e-5)


# --- the resolver returns a real reader (not a wrapper) ---


def test_returns_reader_instance(hdf5_split_paths):
    _inp, pred = hdf5_split_paths
    r = open_confidence_reader(pred)
    try:
        assert isinstance(r, HDF5ConfidenceReader)
    finally:
        r.close()
