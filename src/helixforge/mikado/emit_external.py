"""External-scores emitter: the Helixer prior."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

# Re-exported from helixforge.constants (Phase 19, §1.6). DEFAULT_EXON_WEIGHT is
# the weighting between the exon-confidence component (i) and the
# junction-coincidence component (ii) of helixer_support; DEFAULT_INTRON_HIGH_CUTOFF
# (biology-bearing) is the Helixer intron-probability cutoff for component (ii).
from helixforge.constants import DEFAULT_EXON_WEIGHT, DEFAULT_INTRON_HIGH_CUTOFF
from helixforge.io.hdf5 import CDS, INTRON, UTR

if TYPE_CHECKING:
    from helixforge.io.hdf5 import HDF5ConfidenceReader


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def _exon_bounds(ex: Any) -> tuple[int, int]:
    """``(start, end)`` from an ``Exon`` (``.start``/``.end``) or a 2-tuple."""
    return (ex.start, ex.end) if hasattr(ex, "start") else (ex[0], ex[1])


def _introns_from_exons(exons: list[Any]) -> list[tuple[int, int]]:
    """Return intron ``(start, end)`` tuples between consecutive sorted exons."""
    ordered = sorted(exons, key=lambda e: e.start)
    return [(ordered[i].end, ordered[i + 1].start) for i in range(len(ordered) - 1)]


def _introns_from_pairs(pairs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Intron ``(start, end)`` tuples from ``(start, end)`` pairs (sorted by start).

    Identical to :func:`_introns_from_exons` but over already-extracted bounds, so
    the single-read scorer below need not re-touch exon objects.
    """
    ordered = sorted(pairs, key=lambda p: p[0])
    return [(ordered[i][1], ordered[i + 1][0]) for i in range(len(ordered) - 1)]


# ---------------------------------------------------------------------------
# Single-read scoring
# ---------------------------------------------------------------------------
#
# The per-exon path issues ~2E+1 sliced HDF5 reads per transcript (one per exon
# via get_exon_confidence, one per intron via get_intron_score, plus a full-span
# get_region_confidence that re-reads the locus). Here we fetch the transcript's
# full ``[min_exon_start, max_exon_end)`` prediction block **once** and compute
# every metric by numpy-slicing that in-memory array: collapsing it to one
# h5py read per transcript. This is the *same arithmetic over the same bytes*
# reorganized: a sliced view of the big block at ``[s-block_start, e-block_start)``
# is element-identical to ``get_per_base_predictions(seqid, s, e)`` (the reader's
# stitching is position-wise, independent of the query window), and the reductions
# run in the same order, so the floats are bit-identical (the golden gate holds).


def _score_components(
    exons: list[Any],
    seqid: str,
    h5_reader: HDF5ConfidenceReader,
    exon_weight: float,
    intron_high_cutoff: float,
    conf_span: tuple[int, int] | None = None,
) -> dict[str, float | int]:
    """All Helixer metric components for one transcript from **one** HDF5 read.

    Returns ``{helixer_support, helixer_locus_conf, exon_confidence,
    intron_coincidence_fraction, num_introns}``. Values are bit-identical to the
    per-exon implementations (:meth:`HDF5ConfidenceReader.get_exon_confidence` /
    ``get_intron_score`` / ``get_region_confidence`` composed in ``helixer_support``
    and ``helixer_locus_conf``). The exon sum runs in **passed exon order** (to
    match ``get_exon_confidence``'s running-sum reduction order); the intron loop
    runs in genomic-ascending order (to match ``_introns_from_exons``).

    ``conf_span`` is the ``(start, end)`` over which ``helixer_locus_conf`` is
    taken; default ``None`` uses the exon envelope (``min start``..``max end``),
    which is exactly what the pipeline emitter passes. A caller wanting the legacy
    explicit-span behavior (``score/confidence``) passes the transcript's own
    span; the single read is widened to cover it so the slice stays valid.
    """
    pairs = [_exon_bounds(ex) for ex in exons]
    if not pairs:
        return {
            "helixer_support": 0.0,
            "helixer_locus_conf": 0.0,
            "exon_confidence": 0.0,
            "intron_coincidence_fraction": 0.0,
            "num_introns": 0,
        }

    ex_lo = min(s for s, _ in pairs)
    ex_hi = max(e for _, e in pairs)
    conf_lo, conf_hi = conf_span if conf_span is not None else (ex_lo, ex_hi)
    block_start = min(ex_lo, conf_lo)
    block_end = max(ex_hi, conf_hi)
    # One read for the whole transcript span; every metric slices this in-memory.
    block = h5_reader.get_per_base_predictions(seqid, block_start, block_end)

    # (i) exon-length-weighted mean of max(CDS, UTR), passed exon order.
    total_len = 0
    weighted = 0.0
    for s, e in pairs:
        n = e - s
        if n <= 0:
            continue
        sl = block[s - block_start : e - block_start]
        weighted += float(np.maximum(sl[:, CDS], sl[:, UTR]).sum())
        total_len += n
    exon_conf = _clamp01(weighted / total_len if total_len else 0.0)

    # (ii) fraction of introns whose intron-channel mean >= cutoff.
    introns = _introns_from_pairs(pairs)
    num_introns = len(introns)
    if num_introns:
        high = sum(
            1
            for (s, e) in introns
            if float(block[s - block_start : e - block_start][:, INTRON].mean())
            >= intron_high_cutoff
        )
        intron_frac = high / num_introns
        support = _clamp01(exon_weight * exon_conf + (1.0 - exon_weight) * intron_frac)
    else:
        intron_frac = 0.0
        support = (
            exon_conf  # single-exon: component (i) alone (matches helixer_support)
        )

    # locus confidence == region mean of max(CDS, UTR) over the conf span.
    region = block[conf_lo - block_start : conf_hi - block_start]
    locus_conf = _clamp01(float(np.maximum(region[:, CDS], region[:, UTR]).mean()))

    return {
        "helixer_support": support,
        "helixer_locus_conf": locus_conf,
        "exon_confidence": exon_conf,
        "intron_coincidence_fraction": intron_frac,
        "num_introns": num_introns,
    }


def helixer_support_and_conf(
    exons: list[Any],
    seqid: str,
    strand: str,
    h5_reader: HDF5ConfidenceReader,
    exon_weight: float = DEFAULT_EXON_WEIGHT,
    intron_high_cutoff: float = DEFAULT_INTRON_HIGH_CUTOFF,
) -> tuple[float, float]:
    """``(helixer_support, helixer_locus_conf)`` for a transcript in **one** read.

    Bit-identical to calling :func:`helixer_support` and :func:`helixer_locus_conf`
    separately (the latter over ``exons[0].start..exons[-1].end``), but collapses
    the ~2E+1 sliced HDF5 reads to a single ``get_per_base_predictions`` block.
    ``strand`` is part of the API for symmetry but does not affect the value.
    """
    c = _score_components(exons, seqid, h5_reader, exon_weight, intron_high_cutoff)
    return float(c["helixer_support"]), float(c["helixer_locus_conf"])


def helixer_support(
    exons: list[Any],
    seqid: str,
    strand: str,
    h5_reader: HDF5ConfidenceReader,
    exon_weight: float = DEFAULT_EXON_WEIGHT,
    intron_high_cutoff: float = DEFAULT_INTRON_HIGH_CUTOFF,
) -> float:
    """Per-transcript Helixer structural support in [0, 1].

    Combines two components:
      (i)  exon-length-weighted mean of ``max(CDS, UTR)`` Helixer probability
           over the transcript's exons (``get_exon_confidence``);
      (ii) the fraction of the transcript's introns whose intron-channel score
           (``get_intron_score``) is ``>= intron_high_cutoff``.

    The blend is ``exon_weight*(i) + (1 - exon_weight)*(ii)``. For single-exon
    transcripts there are no introns, so the result is component (i) alone.
    ``strand`` is part of the API for symmetry but does not affect the value
    (the Helixer channels are strand-agnostic per genomic position).
    """
    exon_conf = _clamp01(h5_reader.get_exon_confidence(seqid, exons))

    introns = _introns_from_exons(exons)
    if not introns:
        return exon_conf

    high = sum(
        1
        for (s, e) in introns
        if h5_reader.get_intron_score(seqid, s, e) >= intron_high_cutoff
    )
    intron_frac = high / len(introns)

    support = exon_weight * exon_conf + (1.0 - exon_weight) * intron_frac
    return _clamp01(support)


def helixer_locus_conf(
    seqid: str, start: int, end: int, h5_reader: HDF5ConfidenceReader
) -> float:
    """Helixer region confidence over a span (transcript or locus), clamped [0,1]."""
    return _clamp01(h5_reader.get_region_confidence(seqid, start, end))


def normalize_tpm(tpm: float, locus_max_tpm: float) -> float:
    """Normalize a TPM to [0, 1] as ``min(1, tpm / locus_max_tpm)`` (0 if max<=0)."""
    if locus_max_tpm <= 0:
        return 0.0
    return _clamp01(tpm / locus_max_tpm)


def write_external_scores_tsv(
    rows: dict[str, dict[str, float]], out_path: str | Path
) -> Path:
    """Write the external-scores TSV (``tid`` + metric columns); return ``Path``.

    ``rows`` is ``{tid: {metric: value}}``. Every value must be in [0, 1], a
    violation raises ``ValueError`` (Mikado rejects out-of-range external
    metrics; we never silently clamp here). All rows must share the same metric
    keys.
    """
    out_path = Path(out_path)
    if not rows:
        out_path.write_text("tid\n")
        return out_path

    first_tid = next(iter(rows))
    metric_names = list(rows[first_tid].keys())
    metric_set = set(metric_names)

    with out_path.open("w") as fh:
        fh.write("\t".join(["tid", *metric_names]) + "\n")
        for tid, metrics in rows.items():
            if set(metrics.keys()) != metric_set:
                raise ValueError(
                    f"transcript {tid!r} metrics {sorted(metrics)} differ from "
                    f"expected {sorted(metric_set)}"
                )
            values = []
            for name in metric_names:
                v = metrics[name]
                if not 0.0 <= v <= 1.0:
                    raise ValueError(
                        f"external metric {name!r} for {tid!r} is {v}, "
                        "must be in [0, 1] (Mikado rejects out-of-range values)"
                    )
                values.append(f"{v:.6g}")
            fh.write("\t".join([tid, *values]) + "\n")
    return out_path
