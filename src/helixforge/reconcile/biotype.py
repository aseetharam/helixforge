"""Gene/transcript biotype assignment."""

from __future__ import annotations

from typing import TYPE_CHECKING

import attrs

from helixforge.constants import LNCRNA_MIN_LENGTH, NONCODING_MAX_CDS_CONF
from helixforge.qc.flags import PUTATIVE_CODING, dedup_flags

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene, TranscriptCandidate


def _primary(gene: ReconciledGene) -> TranscriptCandidate:
    for t in gene.transcripts:
        if t.transcript_id == gene.primary_transcript_id:
            return t
    return gene.transcripts[0]


def classify_biotype(
    gene: ReconciledGene,
    *,
    cds_channel_conf: float | None = None,
    min_lncrna_length: int = LNCRNA_MIN_LENGTH,
    max_noncoding_cds_conf: float = NONCODING_MAX_CDS_CONF,
) -> str:
    """Return the biotype string for ``gene`` from its existing signals.

    Pure decision logic (no mutation). ``cds_channel_conf`` is the exon-weighted
    Helixer CDS-channel probability of the primary transcript (from
    :meth:`io.hdf5.HDF5ConfidenceReader.get_cds_channel_confidence`); when ``None``
    the signal is treated as absent (it cannot block a lncRNA call, but the other
    four criteria still must hold). A CDS-channel confidence *above*
    ``max_noncoding_cds_conf`` vetoes a lncRNA call (Helixer still saw coding
    signal → indeterminate, not a confident lncRNA).
    """
    primary = _primary(gene)
    if primary.cds:
        return "protein_coding"
    # ORF-less: a non-coding locus. Decide lncRNA vs ncRNA_undetermined.
    status = gene.classification.status
    expressed = status == "EXPRESSED"
    spliced = primary.num_introns >= 1
    long_enough = primary.total_exon_length >= min_lncrna_length
    cds_quiet = cds_channel_conf is None or cds_channel_conf <= max_noncoding_cds_conf
    if expressed and spliced and long_enough and cds_quiet:
        return "lncRNA"
    return "ncRNA_undetermined"


def _coherent_noncoding_tier(gene: ReconciledGene) -> int:
    """The tier a non-coding gene should carry (never Tier 1/2; D3).

    A lncRNA/ncRNA is not a CDS+homology gene, so it must not sit in Tier 1/2. It
    tiers by expression like a backstop: EXPRESSED/LOW → 3, SILENT → 4. A gene
    already at Tier 3/4 (a CDS-less backstop) keeps its tier.
    """
    if gene.tier >= 3:
        return gene.tier
    status = gene.classification.status
    return 3 if status in ("EXPRESSED", "LOW") else 4


def assign_biotype(
    gene: ReconciledGene,
    *,
    cds_channel_conf: float | None = None,
    min_lncrna_length: int = LNCRNA_MIN_LENGTH,
    max_noncoding_cds_conf: float = NONCODING_MAX_CDS_CONF,
) -> ReconciledGene:
    """Return ``gene`` with its ``biotype`` (and a coherent tier) assigned.

    A gene that already carries a non-None ``biotype`` (e.g. a Phase-29
    ``pseudogene``) is returned unchanged, an explicit biotype set upstream wins.
    Otherwise the gene biotype is derived via :func:`classify_biotype`, mirrored
    onto every transcript (so the GFF3/GTF can emit ``transcript_biotype`` per
    mRNA), and a non-coding biotype demotes a Tier-1/2 gene to its coherent
    non-coding tier. Never drops or restructures the gene.
    """
    if gene.biotype is not None:
        return gene
    biotype = classify_biotype(
        gene,
        cds_channel_conf=cds_channel_conf,
        min_lncrna_length=min_lncrna_length,
        max_noncoding_cds_conf=max_noncoding_cds_conf,
    )
    new_tier = gene.tier
    if biotype in ("lncRNA", "ncRNA_undetermined"):
        new_tier = _coherent_noncoding_tier(gene)
    transcripts = [attrs.evolve(t, biotype=biotype) for t in gene.transcripts]
    flags = gene.flags
    # A coding gene resting solely on its own ORF: no expression, no homology,
    # is honestly *putative* coding (the call stands; we only mark its support
    # level). Expression/homology adjust confidence within coding, never biotype.
    if biotype == "protein_coding":
        primary = _primary(gene)
        expressed = gene.classification.status == "EXPRESSED"
        if not expressed and not primary.has_homology:
            flags = dedup_flags([*gene.flags, PUTATIVE_CODING])
    return attrs.evolve(
        gene, biotype=biotype, tier=new_tier, transcripts=transcripts, flags=flags
    )
