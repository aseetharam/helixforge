"""Optional VCF / haplotype-awareness scaffolding."""

from __future__ import annotations

import gzip
import os
from typing import TYPE_CHECKING, Iterable

import attrs

from helixforge.qc.flags import VARIANT_IMPACTED, dedup_flags
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene

_log = get_logger(__name__)

# SnpEff/VEP sequence-ontology consequence terms that disable or truncate an ORF
# (the "high impact" class). A CDS overlapping any of these is worth flagging.
HIGH_IMPACT_TERMS = frozenset(
    {
        "stop_gained",
        "start_lost",
        "stop_lost",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "transcript_ablation",
        "exon_loss_variant",
    }
)


@attrs.frozen
class Variant:
    """One decomposed (biallelic) VCF record in internal 0-based coordinates.

    ``start``/``end`` is the half-open REF span (``[pos-1, pos-1+len(ref))``), so a
    SNV spans one base and an indel spans its REF length, used for CDS overlap.
    ``high_impact`` is precomputed from the INFO consequence annotation.
    """

    seqid: str
    start: int
    end: int
    ref: str
    alt: str
    high_impact: bool = False


def variant_is_high_impact(info: dict[str, str]) -> bool:
    """Classify a variant's INFO as high-impact (premature stop / splice / frameshift).

    Recognised signals, in order: an explicit ``IMPACT=HIGH``; a SnpEff ``ANN``,
    VEP ``CSQ``, or bcftools ``BCSQ`` field carrying a ``|HIGH|`` impact token or
    one of :data:`HIGH_IMPACT_TERMS`. No annotation → not high-impact (we never
    guess an impact we cannot read).
    """
    if info.get("IMPACT", "").upper() == "HIGH":
        return True
    for key in ("ANN", "CSQ", "BCSQ"):
        ann = info.get(key)
        if not ann:
            continue
        low = ann.lower()
        if "|high|" in low:
            return True
        if any(term in low for term in HIGH_IMPACT_TERMS):
            return True
    return False


def _parse_info(field: str) -> dict[str, str]:
    """Parse a VCF column-8 INFO string into a ``{key: value}`` dict (flags → '')."""
    info: dict[str, str] = {}
    if field in (".", ""):
        return info
    for item in field.split(";"):
        if not item:
            continue
        key, _, value = item.partition("=")
        info[key] = value
    return info


def parse_vcf_line(line: str, *, require_decomposed: bool = True) -> Variant | None:
    """Parse one VCF data line into a :class:`Variant` (``None`` for header/blank).

    A multiallelic record (comma-separated ALT) raises ``ValueError`` when
    ``require_decomposed``, decompose upstream with ``bcftools norm -m -`` first.
    """
    line = line.rstrip("\n")
    if not line or line.startswith("#"):
        return None
    cols = line.split("\t")
    if len(cols) < 8:
        raise ValueError(f"malformed VCF record (need >=8 columns): {line[:60]!r}")
    seqid, pos_s, _id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
    try:
        pos = int(pos_s)
    except ValueError as exc:
        raise ValueError(f"non-numeric VCF POS {pos_s!r} on {seqid}") from exc
    alts = [a for a in alt.split(",") if a not in ("", ".")]
    if require_decomposed and len(alts) > 1:
        raise ValueError(
            f"multiallelic VCF site {seqid}:{pos} has {len(alts)} ALT alleles, "
            "decompose to biallelic records first (`bcftools norm -m -`); "
            "HelixForge will not reason over an un-split multiallelic site."
        )
    info = _parse_info(cols[7])
    start = pos - 1  # 1-based VCF POS -> internal 0-based
    end = start + max(len(ref), 1)
    return Variant(
        seqid=seqid,
        start=start,
        end=end,
        ref=ref,
        alt=alts[0] if alts else "",
        high_impact=variant_is_high_impact(info),
    )


def load_vcf(
    path: str | os.PathLike[str],
    *,
    require_decomposed: bool = True,
) -> list[Variant]:
    """Load a (optionally ``.gz``) VCF into :class:`Variant` records.

    Raises ``FileNotFoundError`` for a missing file and ``ValueError`` for a
    malformed or non-decomposed-multiallelic record (see :func:`parse_vcf_line`).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"VCF not found: {path}")
    opener = gzip.open if str(path).endswith(".gz") else open
    variants: list[Variant] = []
    with opener(path, "rt") as fh:
        for line in fh:
            v = parse_vcf_line(line, require_decomposed=require_decomposed)
            if v is not None:
                variants.append(v)
    return variants


def _index_high_impact(variants: Iterable[Variant]) -> dict[str, list[tuple[int, int]]]:
    """``{seqid: sorted [(start, end), ...]}`` over only the high-impact variants."""
    by_seqid: dict[str, list[tuple[int, int]]] = {}
    for v in variants:
        if v.high_impact:
            by_seqid.setdefault(v.seqid, []).append((v.start, v.end))
    for spans in by_seqid.values():
        spans.sort()
    return by_seqid


def _cds_overlaps_variant(gene: ReconciledGene, spans: list[tuple[int, int]]) -> bool:
    """True if any CDS interval of any transcript overlaps a variant span."""
    for tx in gene.transcripts:
        for seg in tx.cds or ():
            for v_start, v_end in spans:
                if v_start < seg.end and v_end > seg.start:  # half-open overlap
                    return True
    return False


def flag_variant_impacted_genes(
    genes: list[ReconciledGene],
    variants: Iterable[Variant],
    *,
    enabled: bool = False,
) -> list[ReconciledGene]:
    """Attach ``VARIANT_IMPACTED`` to genes whose CDS overlaps a high-impact variant.

    **Disabled by default** → returns ``genes`` unchanged (the identity), so the
    standard / golden path is untouched. When enabled, a gene whose CDS overlaps a
    high-impact variant gets the flag added (deduped); the gene is otherwise
    unchanged, coordinates, tier, biotype, CDS all preserved. The flag is
    INFO-grade interpretation, not a rejection.
    """
    if not enabled:
        return list(genes)
    index = _index_high_impact(variants)
    if not index:
        return list(genes)
    out: list[ReconciledGene] = []
    for gene in genes:
        spans = index.get(gene.seqid)
        if spans and _cds_overlaps_variant(gene, spans):
            out.append(
                attrs.evolve(gene, flags=dedup_flags([*gene.flags, VARIANT_IMPACTED]))
            )
        else:
            out.append(gene)
    return out
