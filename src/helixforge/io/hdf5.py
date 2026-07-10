"""Helixer HDF5 confidence reader."""

from __future__ import annotations

import os
from typing import Any

import h5py
import numpy as np

# Channel indices in /predictions.
INTERGENIC, UTR, CDS, INTRON = 0, 1, 2, 3

# Layout tags returned by detect_helixer_layout.
LAYOUT_COMBINED = "combined"
LAYOUT_SPLIT_METADATA = "split_metadata"
LAYOUT_SPLIT_PREDICTIONS = "split_predictions"


def _layout_of(f: h5py.File) -> str:
    """Classify an open ``h5py.File`` as one of the three Helixer layouts.

    Raises ``ValueError`` if it matches none (e.g. a malformed file missing the
    required datasets).
    """
    has_pred = "predictions" in f
    has_top_meta = "seqids" in f and "start_ends" in f
    has_data_meta = "data" in f and "seqids" in f["data"] and "start_ends" in f["data"]
    if has_pred and has_top_meta:
        return LAYOUT_COMBINED
    if has_data_meta and not has_pred:
        return LAYOUT_SPLIT_METADATA
    if has_pred and not has_top_meta and not has_data_meta:
        return LAYOUT_SPLIT_PREDICTIONS
    raise ValueError(
        "unrecognised Helixer HDF5 layout: expected a combined file "
        "(/predictions + /seqids + /start_ends), a metadata half "
        "(data/seqids + data/start_ends), or a predictions half (/predictions)"
    )


def detect_helixer_layout(path: str | os.PathLike[str]) -> str:
    """Return the layout tag for ``path`` without constructing a reader.

    One of :data:`LAYOUT_COMBINED`, :data:`LAYOUT_SPLIT_METADATA`,
    :data:`LAYOUT_SPLIT_PREDICTIONS`.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"HDF5 not found: {path}")
    with h5py.File(str(path), "r") as f:
        return _layout_of(f)


def _sibling_predictions(meta_path: str | os.PathLike[str]) -> str | None:
    """Infer the ``*_predictions.h5`` path next to a ``*_input.h5`` metadata file.

    Returns ``None`` when the name does not follow Helixer's ``*_input.h5``
    convention (the caller then requires an explicit predictions path).
    """
    p = str(meta_path)
    if p.endswith("_input.h5"):
        return p[: -len("_input.h5")] + "_predictions.h5"
    return None


def _sibling_metadata(pred_path: str | os.PathLike[str]) -> str | None:
    """Infer the ``*_input.h5`` metadata path next to a ``*_predictions.h5`` file.

    The reverse of :func:`_sibling_predictions`. Returns ``None`` when the name
    does not follow Helixer's ``*_predictions.h5`` convention (the caller then
    requires an explicit metadata path).
    """
    p = str(pred_path)
    if p.endswith("_predictions.h5"):
        return p[: -len("_predictions.h5")] + "_input.h5"
    return None


def open_confidence_reader(
    h5_path: str | os.PathLike[str],
    input_h5: str | os.PathLike[str] | None = None,
) -> HDF5ConfidenceReader:
    """Open a reader, resolving Helixer's split input/predictions halves.

    Helixer emits two files: a metadata half (``*_input.h5``, coordinates) and a
    predictions half (``*_predictions.h5``, softmax). The reader needs both. This
    resolves the pair the way the ``confidence`` command does, so callers may pass
    whichever half they have:

    - a **combined** single-file HDF5 → used directly;
    - a **metadata half** (``*_input.h5``) → predictions sibling auto-found;
    - a **predictions half** (``*_predictions.h5``) → metadata taken from
      ``input_h5`` if given, else the ``*_input.h5`` sibling next to it;
    - both halves passed explicitly (``h5_path`` + ``input_h5``, any order) →
      paired via :meth:`HDF5ConfidenceReader.from_helixer_outputs`.

    Raises ``FileNotFoundError`` with an actionable message when a predictions
    half is given but its metadata partner can neither be located nor was
    supplied.
    """
    layout = detect_helixer_layout(h5_path)

    if input_h5 is not None:
        # Both halves named explicitly (order-independent): pair them up by layout.
        other = detect_helixer_layout(input_h5)
        by_layout = {layout: str(h5_path), other: str(input_h5)}
        meta = by_layout.get(LAYOUT_SPLIT_METADATA)
        pred = by_layout.get(LAYOUT_SPLIT_PREDICTIONS)
        if meta is not None and pred is not None:
            return HDF5ConfidenceReader.from_helixer_outputs(meta, pred)
        # Not a clean metadata+predictions pair — fall through to single-path logic
        # (e.g. h5_path is already a combined file and input_h5 is redundant).

    if layout == LAYOUT_SPLIT_PREDICTIONS:
        meta = str(input_h5) if input_h5 is not None else _sibling_metadata(h5_path)
        if meta is None or not os.path.exists(meta):
            raise FileNotFoundError(
                f"{str(h5_path)!r} is a Helixer predictions half; its '*_input.h5' "
                f"metadata partner was not found (looked for {meta!r}). Supply it "
                "with --helixer-input-h5, or place the '*_input.h5' sibling "
                "alongside the predictions file."
            )
        return HDF5ConfidenceReader.from_helixer_outputs(meta, str(h5_path))

    # Combined file, or a metadata half whose predictions sibling __init__ finds.
    return HDF5ConfidenceReader(h5_path)


class HDF5ConfidenceReader:
    """Random-access reader over Helixer per-base predictions."""

    def __init__(
        self,
        h5_path: str | os.PathLike[str],
        predictions_h5: str | os.PathLike[str] | None = None,
    ) -> None:
        """Open a confidence reader, auto-detecting the on-disk layout.

        ``h5_path`` is either a combined file or the metadata half
        (``*_input.h5``). For the split layout, ``predictions_h5`` names the
        softmax half; if omitted it is inferred as the ``*_predictions.h5``
        sibling. Public method signatures and returns are identical for both
        layouts.
        """
        if not os.path.exists(h5_path):
            raise FileNotFoundError(f"HDF5 not found: {h5_path}")
        self.h5_path = str(h5_path)
        self._meta_h5: h5py.File | None = None
        self._pred_h5: h5py.File | None = None

        meta = h5py.File(self.h5_path, "r")
        self._meta_h5 = meta
        try:
            layout = _layout_of(meta)
        except ValueError:
            meta.close()
            self._meta_h5 = None
            raise

        if layout == LAYOUT_SPLIT_PREDICTIONS:
            meta.close()
            self._meta_h5 = None
            raise ValueError(
                f"{self.h5_path!r} is a Helixer predictions half (no metadata). "
                "Pass the matching '*_input.h5' as the primary path, or use "
                "HDF5ConfidenceReader.from_helixer_outputs(input_h5, predictions_h5)."
            )

        self._pred: Any  # h5py Dataset; typed as Any (h5py lacks stubs)
        if layout == LAYOUT_COMBINED:
            self._pred_h5 = meta  # same handle
            self._pred = meta["predictions"]  # left on disk; sliced reads only
            raw_seqids = meta["seqids"][:]
            start_ends = meta["start_ends"][:]
        else:  # LAYOUT_SPLIT_METADATA
            pred_path = predictions_h5 or _sibling_predictions(self.h5_path)
            if pred_path is None or not os.path.exists(str(pred_path)):
                meta.close()
                self._meta_h5 = None
                raise FileNotFoundError(
                    "split-Helixer metadata given but the predictions half was "
                    f"not found: {pred_path!r}. Pass predictions_h5=... explicitly."
                )
            pred = h5py.File(str(pred_path), "r")
            self._pred_h5 = pred
            if "predictions" not in pred:
                self.close()
                raise ValueError(
                    f"predictions file {pred_path!r} missing dataset /predictions"
                )
            self._pred = pred["predictions"]
            raw_seqids = meta["data/seqids"][:]
            start_ends = meta["data/start_ends"][:]

        # Build a per-scaffold chunk index: seqid -> list of (row, low, high, reversed).
        self._index: dict[str, list[tuple[int, int, int, bool]]] = {}
        self._scaffold_len: dict[str, int] = {}
        for row, (sid_raw, (s, e)) in enumerate(zip(raw_seqids, start_ends)):
            sid = sid_raw.decode() if isinstance(sid_raw, bytes) else str(sid_raw)
            reverse = s > e
            low, high = (int(e), int(s)) if reverse else (int(s), int(e))
            self._index.setdefault(sid, []).append((row, low, high, reverse))
            self._scaffold_len[sid] = max(self._scaffold_len.get(sid, 0), high)
        for sid in self._index:
            self._index[sid].sort(key=lambda t: t[1])

    @classmethod
    def from_helixer_outputs(
        cls,
        input_h5: str | os.PathLike[str],
        predictions_h5: str | os.PathLike[str],
    ) -> HDF5ConfidenceReader:
        """Build a reader from a native Helixer output pair (zero-copy).

        ``input_h5`` holds ``data/seqids`` + ``data/start_ends``; ``predictions_h5``
        holds the softmax ``/predictions``. Equivalent to ``cls(input_h5,
        predictions_h5=predictions_h5)`` and returns the same confidence values
        as the equivalent combined file.
        """
        return cls(input_h5, predictions_h5=predictions_h5)

    # --- context manager / lifecycle ---
    def __enter__(self) -> HDF5ConfidenceReader:
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.close()

    def close(self) -> None:
        # In the combined layout both attributes are the same handle; closing an
        # already-closed h5py.File is a no-op, so the dedupe is just tidiness.
        handles: list[h5py.File] = []
        for h in (self._pred_h5, self._meta_h5):
            if h is not None and h not in handles:
                handles.append(h)
        for h in handles:
            h.close()
        self._pred_h5 = None
        self._meta_h5 = None

    @property
    def seqids(self) -> list[str]:
        return sorted(self._index.keys())

    # --- core sliced access ---
    def get_per_base_predictions(self, seqid: str, start: int, end: int) -> np.ndarray:
        """Return an ``(L, 4)`` array for ``[start, end)`` in genomic low→high order.

        Stitches across chunks and trims chunk overlap via numpy slicing (no
        per-base Python loops). Positions not covered by any chunk are zero.
        """
        if seqid not in self._index:
            raise KeyError(f"unknown seqid: {seqid!r}")
        if end <= start:
            raise ValueError(f"require start < end, got start={start} end={end}")
        scaffold_len = self._scaffold_len[seqid]
        if start < 0 or end > scaffold_len:
            raise IndexError(
                f"region {start}-{end} out of bounds for {seqid!r} (len {scaffold_len})"
            )

        out = np.zeros((end - start, 4), dtype=np.float64)
        for row, low, high, reverse in self._index[seqid]:
            if high <= start or low >= end:
                continue
            ov_low = max(low, start)
            ov_high = min(high, end)
            o0 = ov_low - start
            o1 = ov_high - start
            if not reverse:
                c0 = ov_low - low
                c1 = ov_high - low
                out[o0:o1, :] = self._pred[row, c0:c1, :]
            else:
                # genomic position p maps to chunk index (high - 1 - p)
                c0 = high - ov_high
                c1 = high - ov_low
                out[o0:o1, :] = self._pred[row, c0:c1, :][::-1]
        return out

    # --- confidence summaries ---
    def get_region_confidence(self, seqid: str, start: int, end: int) -> float:
        """Mean over the region of ``max(CDS, UTR)`` per-base probability."""
        preds = self.get_per_base_predictions(seqid, start, end)
        return float(np.maximum(preds[:, CDS], preds[:, UTR]).mean())

    def get_exon_confidence(self, seqid: str, exons: Any) -> float:
        """Exon-length-weighted mean of ``max(CDS, UTR)`` over the given exons.

        ``exons`` is an iterable of either ``Exon`` objects or ``(start, end)``
        tuples (0-based half-open).
        """
        total_len = 0
        weighted = 0.0
        for ex in exons:
            s, e = (ex.start, ex.end) if hasattr(ex, "start") else (ex[0], ex[1])
            n = e - s
            if n <= 0:
                continue
            preds = self.get_per_base_predictions(seqid, s, e)
            weighted += float(np.maximum(preds[:, CDS], preds[:, UTR]).sum())
            total_len += n
        if total_len == 0:
            return 0.0
        return weighted / total_len

    def get_intron_score(self, seqid: str, start: int, end: int) -> float:
        """Mean of the intron channel over the region."""
        preds = self.get_per_base_predictions(seqid, start, end)
        return float(preds[:, INTRON].mean())

    def get_cds_channel_confidence(self, seqid: str, exons: Any) -> float:
        """Exon-length-weighted mean of the **CDS channel alone** over the exons.

        Distinct from :meth:`get_exon_confidence` (which is ``max(CDS, UTR)``):
        this returns the coding-channel probability on its own, the signal used
        to discriminate a genuine non-coding RNA (Helixer's
        CDS channel is quiet) from a fragmentary/failed coding locus (Helixer still
        predicts coding signal but no ORF was admissible). ``exons`` is an iterable
        of ``Exon`` objects or ``(start, end)`` tuples (0-based half-open).
        """
        total_len = 0
        weighted = 0.0
        for ex in exons:
            s, e = (ex.start, ex.end) if hasattr(ex, "start") else (ex[0], ex[1])
            n = e - s
            if n <= 0:
                continue
            preds = self.get_per_base_predictions(seqid, s, e)
            weighted += float(preds[:, CDS].sum())
            total_len += n
        if total_len == 0:
            return 0.0
        return weighted / total_len
