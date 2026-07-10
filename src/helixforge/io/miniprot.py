"""miniprot GFF3 parsing."""

from __future__ import annotations

import os
from typing import Any

from helixforge.io.gff import build_gffutils_db
from helixforge.reconcile.models import CDSSegment, MiniprotAlignment
from helixforge.utils.regions import gff3_to_internal, internal_to_gff3


def _phase_from_frame(frame: str) -> int:
    return int(frame) if frame in ("0", "1", "2") else 0


def _first_attr(feature: Any, key: str, default: Any = None) -> Any:
    vals = feature.attributes.get(key)
    if not vals:
        return default
    return " ".join(vals)


def _parse_target(feature: Any) -> tuple[str, int, int] | None:
    """Return ``(protein_id, target_start, target_end)`` (1-based) or None."""
    target = _first_attr(feature, "Target")
    if not target:
        return None
    parts = target.split()
    if len(parts) < 3:
        return None
    return parts[0], int(parts[1]), int(parts[2])


def _normalize_identity(raw: Any) -> float:
    """Identity attribute → fraction in [0,1] (handles 0.95, '95', '88.0%')."""
    if raw is None:
        return 0.0
    s = str(raw).strip().rstrip("%")
    val = float(s)
    if val > 1.0:
        val = val / 100.0
    return val


class MiniprotParser:
    """Parser for miniprot protein-alignment GFF3.

    ``dbfn``/``keep_db``: pass an on-disk path to build the gffutils
    DB once and reuse it instead of rebuilding ``:memory:`` each construction; the
    default (``None``) is the historical in-memory behavior.
    """

    def __init__(
        self,
        gff_path: str | os.PathLike[str],
        dbfn: str | None = None,
        keep_db: bool = False,
    ) -> None:
        if not os.path.exists(gff_path):
            raise FileNotFoundError(f"miniprot GFF3 not found: {gff_path}")
        self.gff_path = str(gff_path)
        self._db: Any = build_gffutils_db(self.gff_path, dbfn=dbfn, keep_db=keep_db)
        # protein_id -> max target_end (used as the protein-length denominator).
        self._protein_len: dict[str, int] = {}
        for feat in self._alignment_features():
            tgt = _parse_target(feat)
            if tgt is None:
                continue
            pid, _ts, te = tgt
            self._protein_len[pid] = max(self._protein_len.get(pid, 0), te)

    def _alignment_features(self) -> Any:
        """Yield the per-alignment features (mRNA, or gene with no mRNA child)."""
        seen_mrna = False
        for mrna in self._db.features_of_type("mRNA"):
            seen_mrna = True
            yield mrna
        if not seen_mrna:
            for gene in self._db.features_of_type("gene"):
                yield gene

    def _build_alignment(self, feat: Any) -> MiniprotAlignment | None:
        tgt = _parse_target(feat)
        if tgt is None:
            return None
        pid, t_start, t_end = tgt

        cds_feats = sorted(
            self._db.children(feat, featuretype="CDS"), key=lambda c: c.start
        )
        if not cds_feats:
            return None
        cds_segments = [
            CDSSegment(*gff3_to_internal(c.start, c.end), _phase_from_frame(c.frame))
            for c in cds_feats
        ]

        denom = self._protein_len.get(pid, t_end)
        target_span = t_end - t_start + 1
        coverage = target_span / denom if denom else 0.0
        coverage = min(1.0, max(0.0, coverage))

        identity = _normalize_identity(_first_attr(feat, "Identity"))
        rank = int(float(_first_attr(feat, "Rank", "0")))
        try:
            score = float(feat.score)
        except (TypeError, ValueError):
            score = 0.0

        start, end = gff3_to_internal(feat.start, feat.end)
        return MiniprotAlignment(
            protein_id=pid,
            seqid=feat.seqid,
            start=start,
            end=end,
            strand=feat.strand,
            cds_segments=cds_segments,
            query_coverage=coverage,
            identity=identity,
            score=score,
            rank=rank,
        )

    def parse(
        self, min_coverage: float = 0.0, min_identity: float = 0.0
    ) -> list[MiniprotAlignment]:
        """Parse all alignments, filtered and sorted by (seqid, start)."""
        out: list[MiniprotAlignment] = []
        for feat in self._alignment_features():
            aln = self._build_alignment(feat)
            if aln is None:
                continue
            if aln.query_coverage < min_coverage or aln.identity < min_identity:
                continue
            out.append(aln)
        out.sort(key=lambda a: (a.seqid, a.start))
        return out

    def parse_for_region(
        self,
        seqid: str,
        start: int,
        end: int,
        min_coverage: float = 0.0,
        min_identity: float = 0.0,
    ) -> list[MiniprotAlignment]:
        """Parse alignments overlapping ``[start, end)`` (internal coords)."""
        g_start, g_end = internal_to_gff3(start, end)
        out: list[MiniprotAlignment] = []
        for feat in self._db.region(
            seqid=seqid, start=g_start, end=g_end, featuretype="mRNA"
        ):
            aln = self._build_alignment(feat)
            if aln is None:
                continue
            if aln.query_coverage < min_coverage or aln.identity < min_identity:
                continue
            out.append(aln)
        out.sort(key=lambda a: (a.seqid, a.start))
        return out

    @staticmethod
    def get_best_per_locus(
        alignments: list[MiniprotAlignment],
        seqid: str,
        start: int,
        end: int,
        strand: str,
    ) -> list[MiniprotAlignment]:
        """Same-strand alignments overlapping the locus, sorted by rank↑, score↓."""
        hits = [
            a
            for a in alignments
            if a.seqid == seqid
            and a.strand == strand
            and a.start < end
            and a.end > start
        ]
        hits.sort(key=lambda a: (a.rank, -a.score))
        return hits
