"""Pseudogene-candidate typing."""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterable

import attrs

from helixforge.qc.flags import INTERNAL_STOP, PSEUDOGENE_CANDIDATE, dedup_flags

if TYPE_CHECKING:
    from helixforge.reconcile.models import QCFlag, ReconciledGene, TranscriptCandidate


def _primary(gene: ReconciledGene) -> TranscriptCandidate:
    for t in gene.transcripts:
        if t.transcript_id == gene.primary_transcript_id:
            return t
    return gene.transcripts[0]


def _homology_backed(transcript: TranscriptCandidate, min_homology: float) -> bool:
    """True when the ORF has protein-homology coverage above ``min_homology``.

    A real hit accession (``protein_id``) always counts. Otherwise the BLAST/
    DIAMOND score must be present **and** at/above ``min_homology`` — so a stronger
    threshold demands stronger homology evidence than the bare ``has_homology``
    gate. With the default ``min_homology=0.0`` this reduces to ``has_homology``.
    """
    if transcript.protein_id is not None:
        return True
    return (
        transcript.blast_score is not None
        and transcript.blast_score > 0
        and transcript.blast_score >= min_homology
    )


def _has_disabling_lesion(gene: ReconciledGene, gate_flags: Iterable[QCFlag]) -> bool:
    """True if the primary CDS carries a premature stop or a mod-3 frameshift.

    The premature-stop signal is taken from the structural codon gate
    (``INTERNAL_STOP``) so the two cannot disagree. A complete CDS whose total
    length is not divisible by 3 is a frameshift lesion (defense in depth — the
    model normally blocks this at construction). A *partial* ORF (missing start/
    stop) is **not** a lesion: it is a truncated-but-valid ORF, not a disabled one.
    """
    t = _primary(gene)
    if not t.cds:
        return False
    if any(f.name == INTERNAL_STOP.name for f in gate_flags):
        return True
    if not t.cds_partial and t.total_cds_length % 3 != 0:
        return True
    return False


def is_pseudogene_candidate(
    gene: ReconciledGene,
    gate_flags: Iterable[QCFlag],
    *,
    min_homology: float = 0.0,
) -> bool:
    """Whether ``gene`` is a homology-backed disabled ORF (pseudogene candidate).

    Requires (1) a CDS, (2) homology coverage above ``min_homology``, and (3) a
    disabling lesion (premature stop / mod-3 frameshift). All three are needed —
    a broken ORF *without* homology stays a normal flagged coding gene.
    """
    gate_flags = list(gate_flags)
    if not _primary(gene).cds:
        return False
    if not _homology_backed(_primary(gene), min_homology):
        return False
    return _has_disabling_lesion(gene, gate_flags)


def apply_pseudogene_typing(
    gene: ReconciledGene,
    gate_flags: Iterable[QCFlag],
    *,
    min_homology: float = 0.0,
) -> ReconciledGene:
    """Return ``gene`` re-typed as a pseudogene candidate when it qualifies.

    Sets ``biotype='pseudogene'``, adds ``PSEUDOGENE_CANDIDATE``, and demotes a
    Tier-1 gene to Tier 2 (a pseudogene is not a Tier-1 coding gene; Tiers 2–4 are
    left as-is). A gene that already carries a non-None ``biotype`` is returned
    unchanged (an explicit biotype set elsewhere wins). Non-qualifying genes are
    returned unchanged. Never drops or restructures the gene.
    """
    if gene.biotype is not None:
        return gene
    gate_flags = list(gate_flags)
    if not is_pseudogene_candidate(gene, gate_flags, min_homology=min_homology):
        return gene
    new_tier = 2 if gene.tier == 1 else gene.tier
    return attrs.evolve(
        gene,
        biotype="pseudogene",
        tier=new_tier,
        flags=dedup_flags([*gene.flags, PSEUDOGENE_CANDIDATE]),
    )
