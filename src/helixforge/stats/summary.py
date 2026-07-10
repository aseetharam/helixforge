"""Genome-level summary counts over a reconciled gene set."""

from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING, Any, Iterable

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene


def summarize_genes(genes: Iterable[ReconciledGene]) -> dict[str, Any]:
    """Return the canonical golden-count summary for a reconciled gene set.

    Parameters
    ----------
    genes:
        Iterable of ``reconcile.models.ReconciledGene``.

    Returns
    -------
    dict with keys:
        ``genes``     total gene count
        ``total_tx``  total transcripts across all genes
        ``multi``     genes with >1 transcript (isoform-bearing)
        ``coding``    genes with at least one CDS-bearing transcript
        ``tier``      {tier_int: count}
        ``origin``    {origin_str: count}
        ``biotype``   {biotype_str_or_"none": count}
        ``as_kinds``  {as_event_kind: count}
        ``flags``     {qc_flag_name: count}

    The nested count dicts are plain ``dict``s with deterministic key order
    (sorted) so two summaries compare equal iff the counts match.
    """
    genes = list(genes)

    tier = Counter(g.tier for g in genes)
    origin = Counter(g.origin for g in genes)
    # ``biotype`` is None until Phase 29/30 type a gene; bucket None as "none" so
    # the summary key set is stable (JSON has no None key) and golden-comparable.
    biotype = Counter(g.biotype or "none" for g in genes)
    flags = Counter(f.name for g in genes for f in g.flags)
    as_kinds = Counter(e.kind for g in genes for e in g.as_events)
    total_tx = sum(len(g.transcripts) for g in genes)
    multi = sum(1 for g in genes if len(g.transcripts) > 1)
    coding = sum(1 for g in genes if any(t.cds for t in g.transcripts))

    return {
        "genes": len(genes),
        "total_tx": total_tx,
        "multi": multi,
        "coding": coding,
        "tier": dict(sorted(tier.items())),
        "origin": dict(sorted(origin.items())),
        "biotype": dict(sorted(biotype.items())),
        "as_kinds": dict(sorted(as_kinds.items())),
        "flags": dict(sorted(flags.items())),
    }
