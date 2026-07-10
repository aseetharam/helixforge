"""Standalone Helixer-confidence scoring of any GFF3."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from helixforge.io.gff import GFF3Parser
from helixforge.io.hdf5 import HDF5ConfidenceReader
from helixforge.mikado.emit_external import (  # type: ignore[attr-defined]  # re-exported from constants; mypy --follow-imports=silent can't trace the re-export
    DEFAULT_EXON_WEIGHT,
    DEFAULT_INTRON_HIGH_CUTOFF,
    _score_components,
)
from helixforge.utils.regions import gff3_to_internal, parse_region

if TYPE_CHECKING:
    import pandas as pd

# Below this helixer_support a transcript is counted "low confidence" in the
# summary. Documented, overridable parameter.
DEFAULT_LOW_CONFIDENCE_CUTOFF = 0.5

# Stable column order for the per-transcript table (so an empty DataFrame still
# carries the schema and ``summarize_confidence`` works on it).
TSV_COLUMNS = [
    "gene_id",
    "transcript_id",
    "seqid",
    "strand",
    "start",
    "end",
    "helixer_support",
    "helixer_locus_conf",
    "exon_confidence",
    "intron_coincidence_fraction",
    "num_exons",
    "num_introns",
]


def _adapt_transcript(
    transcript_like: Any,
) -> tuple[str, str, list[Any], int, int]:
    """Normalise a transcript to ``(seqid, strand, exons, start, end)``.

    Accepts both a :class:`~helixforge.reconcile.models.TranscriptCandidate`
    (attributes) and a ``parse_genes_generic`` transcript dict enriched with the
    parent gene's ``seqid``/``strand`` (mapping). The span (``start``/``end``) is
    taken verbatim when present, else derived as the exon envelope.
    """
    if isinstance(transcript_like, dict):
        seqid: str = transcript_like["seqid"]
        strand: str = transcript_like["strand"]
        exons: list[Any] = transcript_like["exons"]
        start: int | None = transcript_like.get("start")
        end: int | None = transcript_like.get("end")
    else:
        seqid = transcript_like.seqid
        strand = transcript_like.strand
        exons = transcript_like.exons
        start = getattr(transcript_like, "start", None)
        end = getattr(transcript_like, "end", None)

    if not exons:
        raise ValueError("transcript has no exons; cannot score confidence")
    if start is None or end is None:
        start = min(e.start for e in exons)
        end = max(e.end for e in exons)
    return seqid, strand, exons, start, end


def score_transcript_confidence(
    transcript_like: Any,
    h5_reader: HDF5ConfidenceReader,
    exon_weight: float = DEFAULT_EXON_WEIGHT,
    intron_high_cutoff: float = DEFAULT_INTRON_HIGH_CUTOFF,
) -> dict[str, Any]:
    """Helixer-confidence metrics for one transcript, all values in [0, 1].

    Returns ``{helixer_support, helixer_locus_conf, exon_confidence,
    intron_coincidence_fraction, num_exons, num_introns}``. The first two are the
    exact Mikado external metrics; ``exon_confidence`` and
    ``intron_coincidence_fraction`` are the two components of ``helixer_support``
    exposed for inspection. A single-exon transcript has no introns, so
    ``intron_coincidence_fraction`` is ``0.0`` and ``helixer_support ==
    exon_confidence``.
    """
    seqid, strand, exons, start, end = _adapt_transcript(transcript_like)

    # One HDF5 read per transcript: helixer_support,
    # helixer_locus_conf and both of its components come from a single prediction
    # block, bit-identical to the per-exon get_exon_confidence/get_intron_score/
    # get_region_confidence path. The span used for locus_conf is the exon
    # envelope (== start..end here, since _adapt_transcript derives them so).
    c = _score_components(
        exons,
        seqid,
        h5_reader,
        exon_weight,
        intron_high_cutoff,
        conf_span=(start, end),
    )
    return {
        "helixer_support": c["helixer_support"],
        "helixer_locus_conf": c["helixer_locus_conf"],
        "exon_confidence": c["exon_confidence"],
        "intron_coincidence_fraction": c["intron_coincidence_fraction"],
        "num_exons": len(exons),
        "num_introns": c["num_introns"],
    }


def _region_to_internal(
    region: str | None,
) -> tuple[str | None, int | None, int | None]:
    """Parse a region string to ``(seqid, lo, hi)`` internal coords (or ``None``).

    ``None`` region → ``(None, None, None)`` (no filter). A bare ``seqid`` →
    ``(seqid, None, None)``. A ``seqid:start-end`` span is 1-based inclusive (the
    user-facing convention) and is converted here, the single boundary, to
    0-based half-open ``[lo, hi)``.
    """
    if region is None:
        return None, None, None
    seqid, start, end = parse_region(region)
    if start is None:
        return seqid, None, None
    assert end is not None  # parse_region returns (None,None) or (int,int) together
    lo, hi = gff3_to_internal(start, end)
    return seqid, lo, hi


def score_annotation(
    gff3_path: str | Path,
    h5_path: str | Path,
    region: str | None = None,
    predictions_h5: str | Path | None = None,
    exon_weight: float = DEFAULT_EXON_WEIGHT,
    intron_high_cutoff: float = DEFAULT_INTRON_HIGH_CUTOFF,
) -> pd.DataFrame:
    """Score every transcript of any GFF3 against the Helixer HDF5 prior.

    Returns a ``pandas.DataFrame`` with one row per transcript (columns
    :data:`TSV_COLUMNS`). ``gff3_path`` need not be HelixForge output — it is
    parsed with ``GFF3Parser.parse_genes_generic`` so any GFF3 works. ``h5_path``
    may be a combined or a split (``*_input.h5``) Helixer file; the predictions
    half is auto-detected or given via ``predictions_h5``. ``region`` (1-based
    ``seqid`` or ``seqid:start-end``) keeps only transcripts on that scaffold (and
    overlapping the span, if given). The reader is always closed.
    """
    import pandas as pd

    r_seqid, r_lo, r_hi = _region_to_internal(region)

    reader = HDF5ConfidenceReader(h5_path, predictions_h5=predictions_h5)
    try:
        genes = GFF3Parser(gff3_path).parse_genes_generic()
        rows = []
        for gene in genes:
            if r_seqid is not None and gene["seqid"] != r_seqid:
                continue
            for tx in gene["transcripts"]:
                exons = tx["exons"]
                if not exons:
                    continue
                t_start = min(e.start for e in exons)
                t_end = max(e.end for e in exons)
                if r_lo is not None and not (t_start < r_hi and t_end > r_lo):
                    continue
                transcript_like = {
                    "seqid": gene["seqid"],
                    "strand": gene["strand"],
                    "exons": exons,
                    "start": t_start,
                    "end": t_end,
                }
                metrics = score_transcript_confidence(
                    transcript_like,
                    reader,
                    exon_weight=exon_weight,
                    intron_high_cutoff=intron_high_cutoff,
                )
                rows.append(
                    {
                        "gene_id": gene["gene_id"],
                        "transcript_id": tx["transcript_id"],
                        "seqid": gene["seqid"],
                        "strand": gene["strand"],
                        "start": t_start,
                        "end": t_end,
                        **metrics,
                    }
                )
        return pd.DataFrame(rows, columns=TSV_COLUMNS)
    finally:
        reader.close()


def summarize_confidence(
    df: pd.DataFrame,
    low_confidence_cutoff: float = DEFAULT_LOW_CONFIDENCE_CUTOFF,
) -> dict[str, Any]:
    """Distribution summary of ``helixer_support`` over a scored DataFrame.

    Returns ``{n_transcripts, mean, median, q25, q75, min, max,
    low_confidence_cutoff, n_below_cutoff}``. An empty DataFrame yields
    ``n_transcripts == 0`` and ``None`` for the distribution statistics.
    """
    support = df["helixer_support"]
    n = int(len(support))
    if n == 0:
        return {
            "n_transcripts": 0,
            "mean_helixer_support": None,
            "median_helixer_support": None,
            "q25_helixer_support": None,
            "q75_helixer_support": None,
            "min_helixer_support": None,
            "max_helixer_support": None,
            "low_confidence_cutoff": low_confidence_cutoff,
            "n_below_cutoff": 0,
        }
    return {
        "n_transcripts": n,
        "mean_helixer_support": float(support.mean()),
        "median_helixer_support": float(support.median()),
        "q25_helixer_support": float(support.quantile(0.25)),
        "q75_helixer_support": float(support.quantile(0.75)),
        "min_helixer_support": float(support.min()),
        "max_helixer_support": float(support.max()),
        "low_confidence_cutoff": low_confidence_cutoff,
        "n_below_cutoff": int((support < low_confidence_cutoff).sum()),
    }


def write_confidence_tsv(df: pd.DataFrame, path: str | Path) -> Path:
    """Write the per-transcript scored DataFrame to a TSV; return ``Path``."""
    path = Path(path)
    df.to_csv(path, sep="\t", index=False)
    return path


def write_confidence_bigwig(
    h5_path: str | Path,
    chrom_sizes: dict[str, int],
    path: str | Path,
    channel: str = "genic",
    predictions_h5: str | Path | None = None,
) -> Path:
    """Write Helixer per-base confidence as a bigWig; return ``Path``.

    Thin wrapper over :func:`helixforge.viz.tracks.write_confidence_bigwig` that
    opens the HDF5 reader (combined or split) from ``h5_path`` and closes it
    afterwards. ``viz/tracks`` is imported lazily so this module's import surface
    stays minimal (and ``pyBigWig``/plotting deps are only touched on demand).
    """
    from helixforge.viz.tracks import write_confidence_bigwig as _write_bigwig

    reader = HDF5ConfidenceReader(h5_path, predictions_h5=predictions_h5)
    try:
        return _write_bigwig(reader, chrom_sizes, path, channel=channel)
    finally:
        reader.close()
