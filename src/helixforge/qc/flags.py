"""QC flag registry: the single source of truth for flag constants."""

from __future__ import annotations

from collections.abc import Iterable

from helixforge.reconcile.models import AS_EVENT_KINDS, QCFlag


def dedup_flags(flags: Iterable[QCFlag]) -> list[QCFlag]:
    """De-duplicate ``flags`` by ``QCFlag.name``, order-preserving, first-wins.

    Canonical dedup logic shared by ``reconcile/{pipeline,mikado_integrate,cds,
    fallback}.py``. Idempotent: ``dedup_flags(dedup_flags(x)) == dedup_flags(x)``.
    """
    seen: set[str] = set()
    out: list[QCFlag] = []
    for f in flags:
        if f.name not in seen:
            seen.add(f.name)
            out.append(f)
    return out


# Internal accumulator, exposed as the public ALL_FLAGS registry below.
_REGISTRY: dict[str, QCFlag] = {}


def _register(flag: QCFlag) -> QCFlag:
    if flag.name in _REGISTRY:
        raise ValueError(f"duplicate flag name {flag.name!r}")
    _REGISTRY[flag.name] = flag
    return flag


# --- evidence ---
HELIXER_ONLY = _register(
    QCFlag(
        "HELIXER_ONLY",
        "evidence",
        "INFO",
        "Locus supported by Helixer only; no expression/protein evidence.",
    )
)
NO_EXPRESSION = _register(
    QCFlag(
        "NO_EXPRESSION",
        "evidence",
        "WARNING",
        "No RNA-seq expression evidence at this locus.",
    )
)
PUTATIVE_CODING = _register(
    QCFlag(
        "PUTATIVE_CODING",
        "evidence",
        "INFO",
        "Protein-coding on the strength of the model's own complete ORF alone, "
        "no RNA-seq expression and no protein homology corroborate it. The coding "
        "call stands (absence of corroboration is neutral, not a non-coding "
        "signal); this flag marks it as putative rather than well-supported.",
    )
)

# --- confidence ---
LOW_CONF = _register(
    QCFlag(
        "LOW_CONF",
        "confidence",
        "WARNING",
        "Low Helixer confidence at this locus.",
    )
)
VERY_LOW_CONF = _register(
    QCFlag(
        "VERY_LOW_CONF",
        "confidence",
        "ERROR",
        "Very low Helixer confidence at this locus.",
    )
)

# --- structure ---
NO_START = _register(
    QCFlag(
        "NO_START",
        "structure",
        "ERROR",
        "CDS does not begin with a start codon (ATG).",
    )
)
NO_STOP = _register(
    QCFlag(
        "NO_STOP",
        "structure",
        "ERROR",
        "CDS does not end with a stop codon.",
    )
)
INTERNAL_STOP = _register(
    QCFlag(
        "INTERNAL_STOP",
        "structure",
        "CRITICAL",
        "Premature (internal) stop codon in the CDS translation.",
    )
)
AMBIGUOUS_CODON = _register(
    QCFlag(
        "AMBIGUOUS_CODON",
        "structure",
        "INFO",
        "A start/stop codon overlaps an N/IUPAC ambiguous base; the codon is "
        "indeterminate (neither confirmed nor rejected), emitted instead of "
        "NO_START/NO_STOP.",
    )
)
PARTIAL_ORF = _register(
    QCFlag(
        "PARTIAL_ORF",
        "structure",
        "WARNING",
        "CDS is a partial ORF (missing start and/or stop); total CDS length need "
        "not be divisible by 3 (frame preserved by 5' segment phase).",
    )
)
BACKSTOP_RESCUED = _register(
    QCFlag(
        "BACKSTOP_RESCUED",
        "structure",
        "INFO",
        "Silent Helixer-only backstop gene rescued with a valid projected CDS "
        "(re-tiered out of Tier 3/4).",
    )
)
SHORT_CDS = _register(
    QCFlag(
        "SHORT_CDS",
        "structure",
        "WARNING",
        "CDS is shorter than the minimum expected length.",
    )
)
SHORT_EXON = _register(
    QCFlag(
        "SHORT_EXON",
        "structure",
        "WARNING",
        "Exon shorter than the minimum expected length.",
    )
)
LONG_INTRON = _register(
    QCFlag(
        "LONG_INTRON",
        "structure",
        "WARNING",
        "Intron longer than the maximum expected length.",
    )
)
VARIANT_IMPACTED = _register(
    QCFlag(
        "VARIANT_IMPACTED",
        "structure",
        "WARNING",
        "This gene's CDS overlaps a high-impact / originally-multiallelic variant "
        "(premature stop, splice-disrupting, frameshift) from an optional VCF input "
        "EXPERIMENTAL: input plumbing + flag only; off by default "
        "and never re-types/re-tiers the gene. Multiallelic sites must be decomposed "
        "(`bcftools norm -m -`) upstream.",
    )
)

# --- homology ---
NO_HOMOL = _register(
    QCFlag(
        "NO_HOMOL",
        "homology",
        "INFO",
        "No protein homology support for this transcript.",
    )
)
DOMAIN_COMPLETE = _register(
    QCFlag(
        "DOMAIN_COMPLETE",
        "homology",
        "INFO",
        "The primary ORF spans a complete recognized protein domain (Pfam/InterPro "
        "via the functional-annotation hook), a soft ORF-credibility signal "
        "stronger than a bare ORF of equal length. Opt-in: emitted "
        "only when functional annotation ran.",
    )
)

# --- splice ---
ALL_JUNCTIONS_SUPPORTED = _register(
    QCFlag(
        "ALL_JUNCTIONS_SUPPORTED",
        "splice",
        "INFO",
        "All introns are supported by filtered splice junctions.",
    )
)
PARTIAL_JUNCTION_SUPPORT = _register(
    QCFlag(
        "PARTIAL_JUNCTION_SUPPORT",
        "splice",
        "WARNING",
        "Some but not all introns are junction-supported.",
    )
)
NO_JUNCTION_SUPPORT = _register(
    QCFlag(
        "NO_JUNCTION_SUPPORT",
        "splice",
        "WARNING",
        "No introns are supported by filtered splice junctions.",
    )
)
NON_CANONICAL_SPLICE = _register(
    QCFlag(
        "NON_CANONICAL_SPLICE",
        "splice",
        "INFO",
        "One or more of this gene's introns use a non-canonical splice motif "
        "(not GT-AG / GC-AG / AT-AC). In backstop junction correction a "
        "non-canonical motif is a hard reject.",
    )
)

# --- locus (new in v3) ---
LOCUS_SPLIT = _register(
    QCFlag(
        "LOCUS_SPLIT",
        "locus",
        "INFO",
        "One Helixer locus split into multiple reconciled genes.",
    )
)
LOCUS_MERGE = _register(
    QCFlag(
        "LOCUS_MERGE",
        "locus",
        "INFO",
        "Multiple Helixer loci merged into one reconciled gene.",
    )
)
MERGE_REJECTED = _register(
    QCFlag(
        "MERGE_REJECTED",
        "locus",
        "WARNING",
        "A proposed locus merge was rejected during reconciliation.",
    )
)
NOVEL_LOCUS = _register(
    QCFlag(
        "NOVEL_LOCUS",
        "locus",
        "INFO",
        "Reconciled gene with no corresponding Helixer locus (evidence-only).",
    )
)
FROM_TANGLED_LOCUS = _register(
    QCFlag(
        "FROM_TANGLED_LOCUS",
        "locus",
        "INFO",
        "Novel-style gene retained from a tangled many-to-many component where a "
        "Mikado locus lost all its Helixer claimants.",
    )
)
CDS_DISAGREE = _register(
    QCFlag(
        "CDS_DISAGREE",
        "locus",
        "WARNING",
        "miniprot and Mikado CDS projections disagree.",
    )
)
PSEUDOGENE_CANDIDATE = _register(
    QCFlag(
        "PSEUDOGENE_CANDIDATE",
        "locus",
        "WARNING",
        "Homology-backed CDS carries a disabling lesion (premature stop / mod-3 "
        "frameshift), a processed/unitary pseudogene candidate, not a failed "
        "coding gene. Gene biotype is set to 'pseudogene' and the "
        "gene is demoted out of Tier 1.",
    )
)


# --- transposable element (optional EDTA gating) ---
TE_OVERLAP = _register(
    QCFlag(
        "TE_OVERLAP",
        "locus",
        "INFO",
        "The model overlaps a transposable-element feature from an optional EDTA "
        "TE annotation (only configured TE classes count, satellites/knobs/"
        "low-complexity are excluded). Flag only: emitted whenever any TE overlap "
        "occurs, independent of whether the overlap is large enough to reclassify "
        "the gene as a transposable element.",
    )
)


# --- factories for dynamic flags ---
def tier_flag(n: int) -> QCFlag:
    """Return the ``TIER_<n>`` flag (n ∈ {1, 2, 3, 4})."""
    if n not in (1, 2, 3, 4):
        raise ValueError(f"tier must be 1, 2, 3 or 4, got {n}")
    return QCFlag(f"TIER_{n}", "confidence", "INFO", f"Quality tier {n}.")


def as_event_flag(kind: str) -> QCFlag:
    """Return the ``AS_<kind>`` flag for an alternative-splicing event kind."""
    if kind not in AS_EVENT_KINDS:
        raise ValueError(f"kind must be one of {AS_EVENT_KINDS}, got {kind!r}")
    return QCFlag(
        f"AS_{kind}", "splice", "INFO", f"Alternative-splicing event: {kind}."
    )


# Public registry of the static, named flag constants.
ALL_FLAGS: dict[str, QCFlag] = dict(_REGISTRY)


def get_flag(name: str) -> QCFlag:
    """Return the registered :class:`QCFlag` for ``name``; KeyError if missing."""
    return ALL_FLAGS[name]
