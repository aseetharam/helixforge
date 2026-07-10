"""Optional EDTA transposable-element gating.

EDTA is the *only* signal that calls a feature a transposon, no TE is ever
inferred from any other source. With ``--te-annotation`` the pipeline computes
each model's overlap with EDTA TE features and:

* always emits :data:`~helixforge.qc.flags.TE_OVERLAP` when any overlap occurs
  (flag-only inspection, see what *would* be reclassified), and, separately,
* **gates** the coding call: a good-ORF gene whose model-fraction TE overlap is
  at or above the threshold is reclassified ``transposable_element`` (TE-encoded
  transposase / gag-pol ORFs are real ORFs, but they are not host genes).

The gate uses the EDTA ``Classification`` attribute, not bare intersection: EDTA
also annotates non-TE repeats (knob/satellite, centromere, rDNA, low-complexity),
which must never count toward the coding gate. Only the configured TE orders do.

Without ``--te-annotation`` none of this runs and behavior is identical to today.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterable

import attrs

from helixforge.qc.flags import TE_OVERLAP, dedup_flags
from helixforge.utils.intervals import IntervalIndex
from helixforge.utils.logging import get_logger
from helixforge.utils.regions import gff3_to_internal

if TYPE_CHECKING:
    from helixforge.reconcile.models import Exon, ReconciledGene, TranscriptCandidate

_log = get_logger(__name__)

# Top-level EDTA ``Classification`` order tokens (the part before ``/``,
# lowercased) that are genuine transposable elements. Satellites (knob),
# centromeric/subtelomeric/rDNA repeats, simple/low-complexity repeats are
# deliberately excluded: they are repeats, not transposons.
DEFAULT_TE_CLASSES: frozenset[str] = frozenset(
    {"ltr", "dna", "mite", "tir", "helitron", "line", "sine"}
)

# Default model-fraction overlap at/above which a good-ORF gene is reclassified a
# transposable element. High enough that a small incidental overlap (a TE abutting
# a real gene's UTR) does not flag a real gene.
DEFAULT_TE_OVERLAP_THRESHOLD = 0.5

# EDTA biotype assigned to a TE-gated gene. Matches the biotype the publication
# filter already excludes (``utils/filters.publication_ready``).
TE_BIOTYPE = "transposable_element"


def _classification_order(attrs_field: str) -> str | None:
    """Return the lowercased order token of a GFF3 ``Classification=Order/Super``.

    EDTA writes ``Classification=LTR/Gypsy`` etc. The order (before ``/``) is the
    gate key. Returns ``None`` if the attribute is absent.
    """
    for chunk in attrs_field.split(";"):
        chunk = chunk.strip()
        if chunk.startswith("Classification="):
            value = chunk[len("Classification=") :].strip()
            return value.split("/", 1)[0].strip().lower()
    return None


def parse_edta_te_intervals(
    gff3_path: str,
    te_classes: Iterable[str] | None = None,
) -> dict[str, IntervalIndex]:
    """Parse an EDTA TE GFF3 → ``{seqid: IntervalIndex}`` of TE-class intervals.

    Only features whose ``Classification`` order is in ``te_classes`` (default
    :data:`DEFAULT_TE_CLASSES`) are kept; knob/satellite/centromere/rDNA/
    low-complexity rows are dropped. Coordinates are converted from 1-based
    inclusive GFF3 to internal 0-based half-open (I/O boundary). Overlapping TE
    intervals are merged per seqid so the overlap fraction never double-counts a
    base.
    """
    keep = (
        DEFAULT_TE_CLASSES
        if te_classes is None
        else frozenset(c.strip().lower() for c in te_classes if c.strip())
    )
    raw: dict[str, list[tuple[int, int]]] = {}
    kept = dropped = 0
    with open(gff3_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            order = _classification_order(cols[8])
            if order is None or order not in keep:
                dropped += 1
                continue
            try:
                start, end = gff3_to_internal(int(cols[3]), int(cols[4]))
            except ValueError:
                continue
            if end <= start:
                continue
            raw.setdefault(cols[0], []).append((start, end))
            kept += 1
    index: dict[str, IntervalIndex] = {}
    for seqid, ivals in raw.items():
        merged = _merge_intervals(ivals)
        idx = IntervalIndex()
        idx.add_intervals([(s, e) for s, e in merged])
        index[seqid] = idx
    _log.info(
        "EDTA TE annotation: kept %d TE-class features on %d seqids "
        "(%d non-TE features dropped); classes=%s",
        kept,
        len(index),
        dropped,
        sorted(keep),
    )
    return index


def _merge_intervals(ivals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/touching ``[start, end)`` intervals (sorted ascending)."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(ivals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def te_overlap_fraction(exons: list["Exon"], te_index: IntervalIndex) -> float:
    """Fraction of the model's exonic bases that overlap a TE interval.

    Model-fraction (not reciprocal): TE-overlapped exonic bases / total exonic
    bases. The TE intervals are pre-merged, so the per-exon clipped overlaps sum
    without double-counting. Returns ``0.0`` for an empty model.
    """
    total = sum(e.end - e.start for e in exons)
    if total == 0:
        return 0.0
    overlap = 0
    for ex in exons:
        for s, e, *_ in te_index.query_with_data(ex.start, ex.end):
            lo, hi = max(ex.start, s), min(ex.end, e)
            if hi > lo:
                overlap += hi - lo
    return overlap / total


def _primary(gene: "ReconciledGene") -> "TranscriptCandidate":
    for t in gene.transcripts:
        if t.transcript_id == gene.primary_transcript_id:
            return t
    return gene.transcripts[0]


def gate_te(
    genes: list["ReconciledGene"],
    te_index: dict[str, IntervalIndex],
    *,
    threshold: float = DEFAULT_TE_OVERLAP_THRESHOLD,
) -> tuple[list["ReconciledGene"], int, int]:
    """Flag TE overlaps and gate the coding call on TE-class overlap.

    Returns ``(genes, n_flagged, n_reclassified)``:

    * ``n_flagged``, genes that got :data:`TE_OVERLAP` (any TE overlap > 0).
    * ``n_reclassified``, good-ORF (``protein_coding``) genes whose overlap was
      at/above ``threshold`` and were reclassified ``transposable_element`` and
      demoted to Tier 4.

    Flag emission and the gating action are kept distinct: TE_OVERLAP marks every
    overlapping gene regardless of size, so a user can inspect what *would* be
    reclassified before it changes output. Never drops a gene or changes structure.
    """
    out: list["ReconciledGene"] = []
    n_flagged = n_reclassified = 0
    for gene in genes:
        idx = te_index.get(gene.seqid)
        if idx is None or idx.is_empty:
            out.append(gene)
            continue
        frac = te_overlap_fraction(_primary(gene).exons, idx)
        if frac <= 0.0:
            out.append(gene)
            continue
        n_flagged += 1
        gene = attrs.evolve(gene, flags=dedup_flags([*gene.flags, TE_OVERLAP]))
        if frac >= threshold and gene.biotype == "protein_coding":
            transcripts = [
                attrs.evolve(t, biotype=TE_BIOTYPE) for t in gene.transcripts
            ]
            gene = attrs.evolve(
                gene, biotype=TE_BIOTYPE, tier=4, transcripts=transcripts
            )
            n_reclassified += 1
        out.append(gene)
    return out, n_flagged, n_reclassified
