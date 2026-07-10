"""Per-run decision telemetry."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene


@dataclass
class RunStats:
    """Mutable per-run decision counters; emitted once at run end."""

    # locus correspondence (reconcile)
    merges_accepted: int = 0
    merges_rejected: int = 0
    splits: int = 0
    # backstop CDS rescue source (cds.assign_backstop_cds)
    backstop_rescued_miniprot: int = 0
    backstop_rescued_transdecoder: int = 0
    backstop_rescued_none: int = 0
    # backstop junction correction (fallback.refine_backstop_gene)
    junction_corrections_applied: int = 0
    junction_corrections_reverted: int = 0
    # isoform admission (reconcile._adopt)
    isoforms_dropped_redundant: int = 0
    # derived from the final gene set
    partial_orfs: int = 0
    tier_counts: dict[int, int] = field(default_factory=dict)
    origin_counts: dict[str, int] = field(default_factory=dict)
    # Biotype tally (None bucketed as "none").
    biotype_counts: dict[str, int] = field(default_factory=dict)

    # The decision counters summed by :meth:`merge` (the derived tallies —
    # ``partial_orfs`` / ``tier_counts`` / ``origin_counts`` — are recomputed from
    # the final gene set via ``populate_from_genes``, never merged).
    _DECISION_COUNTERS = (
        "merges_accepted",
        "merges_rejected",
        "splits",
        "backstop_rescued_miniprot",
        "backstop_rescued_transdecoder",
        "backstop_rescued_none",
        "junction_corrections_applied",
        "junction_corrections_reverted",
        "isoforms_dropped_redundant",
    )

    def bump(self, name: str, n: int = 1) -> None:
        """Add ``n`` to integer counter ``name`` (no-op-safe single chokepoint)."""
        setattr(self, name, getattr(self, name) + n)

    def merge(self, other: "RunStats") -> None:
        """Fold another ``RunStats``'s decision counters into this one.

        Used by the parallel finalize stage: each worker bumps its
        own ``RunStats`` (a process can't mutate the parent's), and the deltas are
        summed back here so the run telemetry is identical to the serial path. Only
        the decision counters are summed; the derived tallies are recomputed from
        the final gene set, so merging them would double-count.
        """
        for name in self._DECISION_COUNTERS:
            self.bump(name, getattr(other, name))

    def populate_from_genes(self, genes: list[ReconciledGene]) -> None:
        """Fill the derived tallies (tier / origin / partial-ORF) from ``genes``.

        Idempotent given a fixed gene set — recomputes the derived counters from
        scratch each call, so it never double-counts when re-invoked.
        """
        from helixforge.qc.flags import PARTIAL_ORF

        self.tier_counts = dict(sorted(Counter(g.tier for g in genes).items()))
        self.origin_counts = dict(sorted(Counter(g.origin for g in genes).items()))
        self.biotype_counts = dict(
            sorted(Counter(g.biotype or "none" for g in genes).items())
        )
        self.partial_orfs = sum(1 for g in genes if PARTIAL_ORF in g.flags)

    def headline(self) -> dict[str, int]:
        """The flat headline numbers (for the report sidecar / aggregate summary)."""
        return {
            "merges_accepted": self.merges_accepted,
            "merges_rejected": self.merges_rejected,
            "splits": self.splits,
            "backstop_rescued_miniprot": self.backstop_rescued_miniprot,
            "backstop_rescued_transdecoder": self.backstop_rescued_transdecoder,
            "backstop_rescued_none": self.backstop_rescued_none,
            "junction_corrections_applied": self.junction_corrections_applied,
            "junction_corrections_reverted": self.junction_corrections_reverted,
            "isoforms_dropped_redundant": self.isoforms_dropped_redundant,
            "partial_orfs": self.partial_orfs,
        }

    def as_dict(self) -> dict[str, object]:
        """Full telemetry dict (headline + per-tier / per-origin breakdowns)."""
        out: dict[str, object] = dict(self.headline())
        out["tier_counts"] = dict(self.tier_counts)
        out["origin_counts"] = dict(self.origin_counts)
        out["biotype_counts"] = dict(self.biotype_counts)
        return out
