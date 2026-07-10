"""Helixer master locus list: the gene-set anchor."""

from __future__ import annotations

from typing import TypedDict, overload

import attrs

from helixforge.io.gff import GFF3Parser
from helixforge.io.hdf5 import HDF5ConfidenceReader, open_confidence_reader
from helixforge.reconcile.models import Exon, HelixerLocus
from helixforge.utils.intervals import IntervalIndex
from helixforge.utils.logging import get_logger
from helixforge.utils.regions import parse_region

_log = get_logger(__name__)


class _LocusGroup(TypedDict):
    seqid: str
    strand: str
    end: int
    members: list[HelixerLocus]


def load_helixer_loci(
    gff3_path: str,
    h5_path: str | None = None,
    region: str | tuple[str, int, int] | None = None,
    input_h5: str | None = None,
) -> list[HelixerLocus]:
    """Parse → (enrich) → (region-filter) → merge → sorted by (seqid, start).

    ``region`` may be a ``(seqid, start, end)`` tuple or a ``"seqid:start-end"``
    string, in internal 0-based half-open coordinates.

    ``h5_path`` accepts any Helixer HDF5 half (combined, ``*_input.h5``, or
    ``*_predictions.h5``); ``input_h5`` names the metadata partner when
    ``h5_path`` is a bare predictions half, see :func:`open_confidence_reader`.
    """
    loci = GFF3Parser(gff3_path).parse_helixer_genes()

    if h5_path is not None:
        with open_confidence_reader(h5_path, input_h5=input_h5) as reader:
            loci = enrich_with_confidence(loci, reader)

    if region is not None:
        if isinstance(region, str):
            seqid, start, end = parse_region(region)
        else:
            seqid, start, end = region
        loci = filter_loci_by_region(loci, seqid, start, end)

    loci = merge_overlapping_loci(loci)
    loci.sort(key=lambda g: (g.seqid, g.start))
    return loci


def enrich_with_confidence(
    loci: list[HelixerLocus],
    h5_reader: HDF5ConfidenceReader,
) -> list[HelixerLocus]:
    """Return loci with ``confidence`` set from ``get_region_confidence``.

    On scaffold mismatch (seqid absent from the HDF5, or an out-of-bounds
    region) a warning is logged and the locus is left unchanged (``None``).
    """
    available = set(h5_reader.seqids)
    out: list[HelixerLocus] = []
    for locus in loci:
        if locus.seqid not in available:
            _log.warning(
                "scaffold %s not in HDF5; leaving confidence=None for %s",
                locus.seqid,
                locus.gene_id,
            )
            out.append(locus)
            continue
        try:
            conf = h5_reader.get_region_confidence(locus.seqid, locus.start, locus.end)
            out.append(attrs.evolve(locus, confidence=conf))
        except (KeyError, IndexError) as exc:
            _log.warning(
                "could not read confidence for %s (%s:%d-%d): %s",
                locus.gene_id,
                locus.seqid,
                locus.start,
                locus.end,
                exc,
            )
            out.append(locus)
    return out


def _mean_confidence(loci: list[HelixerLocus]) -> float | None:
    vals = [g.confidence for g in loci if g.confidence is not None]
    if not vals:
        return None
    return sum(vals) / len(vals)


@overload
def merge_overlapping_loci(
    loci: list[HelixerLocus],
    return_provenance: bool,
) -> list[HelixerLocus] | tuple[list[HelixerLocus], dict[str, list[str]]]: ...
@overload
def merge_overlapping_loci(
    loci: list[HelixerLocus],
) -> list[HelixerLocus]: ...
def merge_overlapping_loci(
    loci: list[HelixerLocus],
    return_provenance: bool = False,
) -> list[HelixerLocus] | tuple[list[HelixerLocus], dict[str, list[str]]]:
    """Merge overlapping same-strand loci (greedy, single pass).

    Only loci sharing seqid + strand and with genuinely overlapping spans
    (``next.start < current.end``) are merged; touching spans and opposite
    strands are not. A merged locus keeps the **first** gene_id, ``min`` start /
    ``max`` end, the mean of available confidences, and the union of exons
    (CDS is dropped on merge). Non-merged loci are returned unchanged.

    With ``return_provenance=True`` also returns ``{kept_gene_id:
    [constituent_gene_ids...]}`` (every output locus has an entry; singletons
    map to ``[gene_id]``).
    """
    ordered = sorted(loci, key=lambda g: (g.seqid, g.strand, g.start))
    open_groups: list[
        _LocusGroup
    ] = []  # each: {"seqid", "strand", "end" (running max), "members"}

    for locus in ordered:
        cur = open_groups[-1] if open_groups else None
        if (
            cur is not None
            and locus.seqid == cur["seqid"]
            and locus.strand == cur["strand"]
            and locus.start < cur["end"]
        ):
            cur["members"].append(locus)
            cur["end"] = max(cur["end"], locus.end)
        else:
            open_groups.append(
                {
                    "seqid": locus.seqid,
                    "strand": locus.strand,
                    "end": locus.end,
                    "members": [locus],
                }
            )

    groups: list[list[HelixerLocus]] = [g["members"] for g in open_groups]
    result: list[HelixerLocus] = []
    provenance: dict[str, list[str]] = {}
    for constituents in groups:
        if len(constituents) == 1:
            locus = constituents[0]
            result.append(locus)
            provenance[locus.gene_id] = [locus.gene_id]
            continue
        first = constituents[0]
        start = min(g.start for g in constituents)
        end = max(g.end for g in constituents)
        exons: list[Exon] = []
        for g in constituents:
            exons = _merge_exon_lists(exons, g.exons)
        fused = HelixerLocus(
            gene_id=first.gene_id,
            seqid=first.seqid,
            start=start,
            end=end,
            strand=first.strand,
            confidence=_mean_confidence(constituents),
            exons=exons,
            cds=None,
        )
        result.append(fused)
        provenance[first.gene_id] = [g.gene_id for g in constituents]

    if return_provenance:
        return result, provenance
    return result


def _merge_exon_lists(a: list[Exon], b: list[Exon]) -> list[Exon]:
    """Standard interval merge of two exon lists; touching exons stay separate."""
    combined = sorted(list(a) + list(b), key=lambda e: e.start)
    if not combined:
        return []
    merged = [Exon(combined[0].start, combined[0].end)]
    for ex in combined[1:]:
        last = merged[-1]
        if ex.start < last.end:  # strict overlap; end==start (touching) not merged
            if ex.end > last.end:
                merged[-1] = Exon(last.start, ex.end)
        else:
            merged.append(Exon(ex.start, ex.end))
    return merged


def filter_loci_by_region(
    loci: list[HelixerLocus],
    seqid: str | None,
    start: int | None,
    end: int | None,
) -> list[HelixerLocus]:
    """Return loci on ``seqid`` overlapping ``[start, end)``."""
    return [
        g
        for g in loci
        if g.seqid == seqid and g.start < end and g.end > start  # type: ignore[operator]  # None checked by caller
    ]


def build_locus_index(loci: list[HelixerLocus]) -> IntervalIndex:
    """Flat ``IntervalIndex`` over all loci (data = original-list index).

    Coordinate-only, the caller is responsible for single-scaffold use.
    """
    idx = IntervalIndex()
    idx.add_intervals([(g.start, g.end, i) for i, g in enumerate(loci)])
    return idx


def build_locus_index_by_scaffold(loci: list[HelixerLocus]) -> dict[str, IntervalIndex]:
    """Return ``{seqid: IntervalIndex}``; each index's data is the original index.

    Reused by Phase 6 reconciliation to map query hits back to ``loci``.
    """
    by_scaffold: dict[str, list[tuple[int, int, int]]] = {}
    for i, g in enumerate(loci):
        by_scaffold.setdefault(g.seqid, []).append((g.start, g.end, i))
    indices: dict[str, IntervalIndex] = {}
    for seqid, intervals in by_scaffold.items():
        idx = IntervalIndex()
        idx.add_intervals(intervals)
        indices[seqid] = idx
    return indices
