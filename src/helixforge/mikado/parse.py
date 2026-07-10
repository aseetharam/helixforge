"""Mikado output parser."""

from __future__ import annotations

import logging
import os
from typing import Any, cast

_log = logging.getLogger(__name__)

from helixforge.io.gff import build_gffutils_db
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    MikadoLocus,
    TranscriptCandidate,
)
from helixforge.utils.regions import gff3_to_internal

_GENE_TYPES = ("gene", "ncRNA_gene")
_TRANSCRIPT_TYPES = ("mRNA", "transcript", "ncRNA")
# mRNA attribute keys that may carry Mikado's chosen protein target.
_PROTEIN_ATTR_KEYS = ("protein_id", "blast_target", "best_blast_target")


def _phase_from_frame(frame: str) -> int:
    return int(frame) if frame in ("0", "1", "2") else 0


def _coerce(value: str) -> float | bool | str:
    """Coerce a TSV cell to float/bool where possible, else keep the string."""
    try:
        return float(value)
    except ValueError:
        if value == "True":
            return True
        if value == "False":
            return False
        return value


def _read_keyed_tsv(path: str) -> dict[str, dict[str, float | bool | str]]:
    """Read a TSV whose key column is ``tid`` (or the first column)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"TSV not found: {path}")
    out: dict[str, dict[str, float | bool | str]] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        key_idx = header.index("tid") if "tid" in header else 0
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            tid = cols[key_idx]
            out[tid] = {
                header[i]: _coerce(cols[i])
                for i in range(len(header))
                if i < len(cols) and i != key_idx
            }
    return out


def parse_metrics_tsv(path: str) -> dict[str, dict[str, float | bool | str]]:
    """Parse a Mikado metrics TSV → ``{transcript_id: {metric: value}}``."""
    return _read_keyed_tsv(path)


def parse_scores_tsv(path: str) -> dict[str, dict[str, float | bool | str]]:
    """Parse a Mikado scores TSV → ``{transcript_id: {metric: value}}``."""
    return _read_keyed_tsv(path)


def _first_attr(feature: Any, keys: tuple[str, ...]) -> str | None:
    for key in keys:
        vals = feature.attributes.get(key)
        if vals:
            return cast(str, vals[0])
    return None


def parse_loci_gff3(
    loci_gff3: str,
    metrics_tsv: str | None = None,
    scores_tsv: str | None = None,
    dbfn: str | None = None,
    keep_db: bool = False,
) -> list[MikadoLocus]:
    """Parse ``mikado.loci.gff3`` into ``MikadoLocus`` objects.

    Each gene becomes a locus; each child transcript becomes a
    ``TranscriptCandidate(source='mikado')`` with exons + CDS. If the metrics /
    scores TSVs are given, per-transcript dicts are attached to the locus
    (keyed by transcript id) and ``combined_score`` is taken from the scores
    ``score`` column when present.
    """
    if not os.path.exists(loci_gff3):
        raise FileNotFoundError(f"loci GFF3 not found: {loci_gff3}")

    metrics_by_tid = parse_metrics_tsv(metrics_tsv) if metrics_tsv else {}
    scores_by_tid = parse_scores_tsv(scores_tsv) if scores_tsv else {}

    # Default None keeps the in-memory DB; a persistent path avoids re-parsing.
    db: Any = build_gffutils_db(loci_gff3, dbfn=dbfn, keep_db=keep_db)

    loci: list[MikadoLocus] = []
    n_skipped = 0
    for gene_type in _GENE_TYPES:
        for gene in db.features_of_type(gene_type, order_by=("seqid", "start")):
            transcripts: list[TranscriptCandidate] = []
            locus_metrics: dict[str, dict[str, float | bool | str]] = {}
            locus_scores: dict[str, dict[str, float | bool | str]] = {}
            for ttype in _TRANSCRIPT_TYPES:
                for tx in db.children(gene, featuretype=ttype, order_by="start"):
                    try:
                        candidate = _build_candidate(
                            db, tx, gene.id, metrics_by_tid, scores_by_tid
                        )
                    except (ValueError, TypeError) as exc:
                        _log.warning(
                            "skipping transcript %s in locus %s: %s",
                            tx.id,
                            gene.id,
                            exc,
                        )
                        n_skipped += 1
                        continue
                    transcripts.append(candidate)
                    if tx.id in metrics_by_tid:
                        locus_metrics[tx.id] = metrics_by_tid[tx.id]
                    if tx.id in scores_by_tid:
                        locus_scores[tx.id] = scores_by_tid[tx.id]

            if not transcripts:
                _log.warning(
                    "skipping locus %s: all transcripts failed validation",
                    gene.id,
                )
                continue

            transcripts.sort(key=lambda t: t.transcript_id)
            g_start, g_end = gff3_to_internal(gene.start, gene.end)
            loci.append(
                MikadoLocus(
                    locus_id=gene.id,
                    seqid=gene.seqid,
                    start=g_start,
                    end=g_end,
                    strand=gene.strand,
                    transcripts=transcripts,
                    metrics=cast(dict[str, object], locus_metrics),
                    scores=cast(dict[str, object], locus_scores),
                )
            )

    if n_skipped:
        _log.warning(
            "skipped %d transcript(s) due to structural violations", n_skipped
        )

    loci.sort(key=lambda g: (g.seqid, g.start))
    return loci


def _strip_terminal_stop_codon(
    cds: list[CDSSegment], strand: str, n: int = 3
) -> list[CDSSegment] | None:
    """Trim the terminal stop codon (``n`` nt) from the coding 3' end of a CDS.

    Mikado's loci GFF3 stores CDS **including** the stop codon, but the internal
    convention is CDS-excludes-stop. When Mikado's ``has_stop_codon`` metric is
    true this removes the last ``n`` coding nucleotides, walking across CDS
    segments when the stop spans an intron, and returns a new genomic-ascending
    ``CDSSegment`` list, or ``None`` if the CDS is too short to trim (the caller
    then keeps the original).

    Only the 3' (downstream) coding end is shortened, so every upstream phase,
    including the 5'-most segment's frame-fixing phase, is preserved. ``+`` codes
    ascending (trim the high-coordinate end); ``-`` codes descending (trim the
    low-coordinate end). This is a deterministic convention normalisation, **not**
    heuristic frame trimming: mod-3 is unaffected (n=3).
    """
    bounds: list[tuple[int, int, int]] = [(s.start, s.end, s.phase) for s in cds]
    remaining = n
    if strand == "+":
        i = len(bounds) - 1
        while remaining > 0 and i >= 0:
            s, e, ph = bounds[i]
            take = min(remaining, e - s)
            e -= take
            remaining -= take
            if e > s:
                bounds[i] = (s, e, ph)
            else:
                bounds.pop(i)
            i -= 1
    else:
        while remaining > 0 and bounds:
            s, e, ph = bounds[0]
            take = min(remaining, e - s)
            s += take
            remaining -= take
            if e > s:
                bounds[0] = (s, e, ph)
                break
            bounds.pop(0)
    if remaining > 0 or not bounds:
        return None
    return [CDSSegment(s, e, ph) for (s, e, ph) in bounds]


def _build_candidate(
    db: Any,
    tx: Any,
    locus_id: str,
    metrics_by_tid: dict[str, dict[str, float | bool | str]],
    scores_by_tid: dict[str, dict[str, float | bool | str]],
) -> TranscriptCandidate:
    exons = [
        Exon(*gff3_to_internal(e.start, e.end))
        for e in db.children(tx, featuretype="exon", order_by="start")
    ]
    cds_feats = list(db.children(tx, featuretype="CDS", order_by="start"))
    cds: list[CDSSegment] | None = (
        [
            CDSSegment(*gff3_to_internal(c.start, c.end), _phase_from_frame(c.frame))
            for c in cds_feats
        ]
        if cds_feats
        else None
    )

    score_row = scores_by_tid.get(tx.id, {})
    combined_score_raw = score_row.get("score")
    combined_score: float | None = (
        float(combined_score_raw) if combined_score_raw is not None else None
    )

    # Partial ORF? Prefer Mikado's authoritative completeness metrics
    # (has_start_codon / has_stop_codon); fall back to the mod-3 heuristic when
    # the metrics TSV is absent. A partial CDS is exempt from the mod-3 model
    # invariant (reading frame preserved by the 5' segment phase).
    metric_row = metrics_by_tid.get(tx.id, {})
    has_start = metric_row.get("has_start_codon")
    has_stop = metric_row.get("has_stop_codon")

    # CDS is stop-inclusive (§4.5): the stop codon is the last 3 bases of the
    # CDS. Mikado's loci GFF3 already stores CDS including the stop codon, so
    # no trimming is needed: the CDS is kept as-is.

    # 5'/3' partiality are independent: a 5'-partial-but-3'-complete
    # transcript must still have its stop codon verified, and vice versa. Mikado's
    # has_start_codon/has_stop_codon give each end directly; without metrics we
    # fall back to the mod-3 heuristic and treat both ends the same.
    cds_partial_5prime: bool
    cds_partial_3prime: bool
    if has_start is not None and has_stop is not None:
        cds_partial_5prime = not bool(has_start)
        cds_partial_3prime = not bool(has_stop)
    elif cds is not None:
        both = (sum(seg.end - seg.start for seg in cds) % 3) != 0
        cds_partial_5prime = both
        cds_partial_3prime = both
    else:
        cds_partial_5prime = False
        cds_partial_3prime = False

    # Safety: if metrics say complete but the CDS is empirically non-mod-3, the
    # metrics are stale (computed on the prepared transcript, not the picked
    # locus). Override to partial, empirical CDS is the ground truth.
    if cds is not None and not cds_partial_5prime and not cds_partial_3prime:
        total_cds = sum(seg.end - seg.start for seg in cds)
        if total_cds % 3 != 0:
            _log.warning(
                "%s: metrics report complete ORF but CDS total=%d is not "
                "mod-3; overriding to partial (metrics stale after pick)",
                tx.id,
                total_cds,
            )
            cds_partial_5prime = True
            cds_partial_3prime = True

    # Protein-homology support: Mikado's loci GFF3 has no hit accession, but the
    # metrics TSV carries the best-hit BLAST/DIAMOND score, the Tier-1 homology
    # signal (assign_tier / TranscriptCandidate.has_homology).
    blast_score_raw = metric_row.get("blast_score")
    blast_score: float | None = (
        float(blast_score_raw) if blast_score_raw is not None else None
    )

    start, end = gff3_to_internal(tx.start, tx.end)
    return TranscriptCandidate(
        transcript_id=tx.id,
        locus_id=locus_id,
        source="mikado",
        seqid=tx.seqid,
        start=start,
        end=end,
        strand=tx.strand,
        exons=exons,
        cds=cds,
        cds_partial_5prime=cds_partial_5prime,
        cds_partial_3prime=cds_partial_3prime,
        protein_id=_first_attr(tx, _PROTEIN_ATTR_KEYS),
        blast_score=blast_score,
        combined_score=combined_score,
    )
