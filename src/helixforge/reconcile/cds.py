"""Backstop CDS projection + cross-check."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any

import attrs

from helixforge.constants import CROSS_CHECK_OVERLAP as _CROSS_CHECK_OVERLAP
from helixforge.io.miniprot import MiniprotParser
from helixforge.mikado.run import run_transdecoder
from helixforge.qc.flags import BACKSTOP_RESCUED, CDS_DISAGREE, dedup_flags
from helixforge.reconcile.as_events import reciprocal_overlap
from helixforge.reconcile.mikado_integrate import assign_tier
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    MiniprotAlignment,
    QCFlag,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.utils.logging import get_logger
from helixforge.utils.sequences import is_stop_codon

_log = get_logger(__name__)

MIKADO_ORIGINS = ("mikado_1to1", "split", "merge", "novel")

# _CROSS_CHECK_OVERLAP re-exported from helixforge.constants (Phase 19, §1.6):
# below this reciprocal CDS overlap, an independent miniprot ORF is considered to
# materially disagree with the chosen Mikado ORF (cross-check only — never edits).


# ---------------------------------------------------------------------------
# Phase assignment (coding order; first segment = 0)
# ---------------------------------------------------------------------------


def _phased_segments(
    bounds_ascending: list[tuple[int, int]], strand: str
) -> list[CDSSegment]:
    """Build ``CDSSegment``s from ascending ``(start, end)`` bounds with phases.

    Phase of a segment = ``(3 - preceding_cds_len % 3) % 3`` where lengths
    accumulate in **coding** order (ascending on ``+``, descending on ``-``).
    The returned list is in genomic-ascending order.
    """
    coding = bounds_ascending if strand == "+" else list(reversed(bounds_ascending))
    phase_by_bound: dict[tuple[int, int], int] = {}
    running = 0
    for s, e in coding:
        phase_by_bound[(s, e)] = (3 - running % 3) % 3
        running += e - s
    return [CDSSegment(s, e, phase_by_bound[(s, e)]) for (s, e) in bounds_ascending]


def _merge_adjacent(bounds: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge touching/overlapping ``(start, end)`` intervals (sorted ascending)."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(bounds):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


# ---------------------------------------------------------------------------
# D1: projection
# ---------------------------------------------------------------------------


def project_cds_to_exons(
    cds_segments: list[CDSSegment],
    exons: list[Exon],
    strand: str = "+",
) -> list[CDSSegment] | None:
    """Clip CDS segments to exon boundaries; return phased ``CDSSegment``s or None.

    Each input CDS segment is intersected with every exon (pairwise overlap);
    portions falling in introns are dropped. Surviving pieces are merged where
    adjacent, sorted, and phased in coding order. Returns ``None`` if nothing
    survives or the total length is not a multiple of 3 — **no trimming** is done
    to force frame.

    ``strand`` (default ``'+'``) selects the coding order for phase computation;
    it does not change coordinate storage.
    """
    pieces: list[tuple[int, int]] = []
    for seg in cds_segments:
        s0, e0 = seg.start, seg.end
        for ex in exons:
            lo = max(s0, ex.start)
            hi = min(e0, ex.end)
            if hi > lo:
                pieces.append((lo, hi))
    if not pieces:
        return None
    merged = _merge_adjacent(pieces)
    total = sum(e - s for s, e in merged)
    if total == 0 or total % 3 != 0:
        return None
    return _phased_segments(merged, strand)


# ---------------------------------------------------------------------------
# Transcript ↔ genome sequence helpers (both strands)
# ---------------------------------------------------------------------------


def extract_transcript_sequence(transcript: TranscriptCandidate, genome: Any) -> str:
    """Return the spliced transcript sequence in coding (5'→3') orientation.

    Concatenates exon sequences; ascending genomic order on ``+`` and descending
    (reverse-complemented) on ``-`` (``RC(a+b) == RC(b)+RC(a)``).
    """
    exons = list(transcript.exons)  # sorted ascending by model invariants
    if transcript.strand == "-":
        exons = list(reversed(exons))
    return "".join(
        genome.get_sequence(transcript.seqid, ex.start, ex.end, transcript.strand)
        for ex in exons
    )


def map_transcript_to_genomic(
    transcript: TranscriptCandidate,
    tx_start: int,
    tx_end: int,
) -> list[CDSSegment] | None:
    """Map a half-open transcript-coordinate ORF ``[tx_start, tx_end)`` to genomic
    ``CDSSegment``s clipped to exons, phased in coding order (both strands).

    Transcript coordinates run 5'→3' in coding direction (position 0 = the first
    base of :func:`extract_transcript_sequence`). Returns ``None`` if the window
    is empty/out of range or the mapped CDS is not mod-3.
    """
    length = transcript.total_exon_length
    if tx_start < 0 or tx_end > length or tx_end <= tx_start:
        return None

    exons = list(transcript.exons)
    ordered = exons if transcript.strand == "+" else list(reversed(exons))

    bounds: list[tuple[int, int]] = []
    cursor = 0  # transcript-space offset at the 5' edge of the current exon
    for ex in ordered:
        ex_len = len(ex)
        ex_lo_tx, ex_hi_tx = cursor, cursor + ex_len
        # overlap of [tx_start, tx_end) with this exon in transcript space
        ov_lo = max(tx_start, ex_lo_tx)
        ov_hi = min(tx_end, ex_hi_tx)
        if ov_hi > ov_lo:
            off_lo = ov_lo - ex_lo_tx  # bases into the exon (from its 5' edge)
            off_hi = ov_hi - ex_lo_tx
            if transcript.strand == "+":
                g_lo = ex.start + off_lo
                g_hi = ex.start + off_hi
            else:
                # 5' edge of a minus-strand exon is its high genomic coordinate
                g_hi = ex.end - off_lo
                g_lo = ex.end - off_hi
            bounds.append((g_lo, g_hi))
        cursor += ex_len

    if not bounds:
        return None
    merged = _merge_adjacent(bounds)
    total = sum(e - s for s, e in merged)
    if total == 0 or total % 3 != 0:
        return None
    return _phased_segments(merged, transcript.strand)


# ---------------------------------------------------------------------------
# Stop-inclusive extension (miniprot CDS is stop-exclusive by construction)
# ---------------------------------------------------------------------------


def _include_genomic_stop(
    cds_segments: list[CDSSegment],
    exons: list[Exon],
    strand: str,
    genome: Any,
    seqid: str,
    transl_table: int = 1,
) -> list[CDSSegment] | None:
    """Extend stop-exclusive CDS by 3 bp to include the genomic in-frame stop.

    Protein alignments (miniprot) produce CDS that ends at the last coding
    codon, excluding the stop. This checks the single deterministic in-frame
    position (the 3 nt immediately 3' of the CDS in the genome) and, if it is
    a valid stop codon AND the extension fits within an exon, returns the
    extended + re-phased CDS. Returns ``None`` when no valid stop is present or
    the extension overflows the exon — the caller keeps the original CDS and
    ``validate`` flags it.
    """
    if not cds_segments:
        return None
    bounds = [(s.start, s.end) for s in cds_segments]
    if strand == "+":
        last_s, last_e = bounds[-1]
        stop_seq = genome.get_sequence(seqid, last_e, last_e + 3, "+")
        if not stop_seq or len(stop_seq) != 3 or not is_stop_codon(stop_seq, transl_table):
            return None
        new_end = last_e + 3
        if not any(ex.start <= last_s and new_end <= ex.end for ex in exons):
            return None
        bounds[-1] = (last_s, new_end)
    else:
        first_s, first_e = bounds[0]
        if first_s < 3:
            return None
        stop_seq = genome.get_sequence(seqid, first_s - 3, first_s, "-")
        if not stop_seq or len(stop_seq) != 3 or not is_stop_codon(stop_seq, transl_table):
            return None
        new_start = first_s - 3
        if not any(ex.start <= new_start and first_e <= ex.end for ex in exons):
            return None
        bounds[0] = (new_start, first_e)
    return _phased_segments(bounds, strand)


# ---------------------------------------------------------------------------
# D1: backstop CDS assignment
# ---------------------------------------------------------------------------


def _rescue_backstop_cds(
    gene: ReconciledGene,
    cds_segments: list[CDSSegment],
    protein_id: str | None = None,
) -> ReconciledGene:
    """Return ``gene`` with its backstop transcript given ``cds_segments``,
    re-tiered and flagged ``BACKSTOP_RESCUED``.

    Uses ``attrs.evolve`` so model validation re-runs — a CDS
    that violates any invariant raises and is treated by the caller as a failed
    projection. The tier is recomputed from the rescued primary via
    :func:`assign_tier`: a now-CDS-bearing backstop tiers like a Mikado-origin
    gene — Tier 1 if the CDS is homology-backed (a real miniprot ``protein_id``),
    Tier 2 if it is CDS-only (e.g. a TransDecoder ORF with no homology) — instead
    of staying at the silent-backstop Tier 3/4.
    """
    t = gene.transcripts[0]
    new_t = attrs.evolve(t, cds=list(cds_segments), protein_id=protein_id)
    tier = assign_tier(new_t, gene.classification, gene.origin)
    new_gene = attrs.evolve(
        gene,
        transcripts=[new_t],
        tier=tier,
        flags=dedup_flags([*gene.flags, BACKSTOP_RESCUED]),
    )
    return new_gene


def _transdecoder_cds(
    transcript: TranscriptCandidate,
    genome: Any,
    transdecoder_bin_dir: str,
) -> list[CDSSegment] | None:
    """Predict an ORF on the spliced transcript with TransDecoder → CDS or None.

    Writes the single coding-orientation transcript to a FASTA, runs TransDecoder
    (subprocess via :func:`run_transdecoder`), reads the ``thickStart``/
    ``thickEnd`` (ORF span, transcript coordinates) from the resulting BED12, and
    maps it back to genomic CDS. Mocked in unit tests.
    """
    seq = extract_transcript_sequence(transcript, genome)
    with tempfile.TemporaryDirectory() as work:
        fasta = os.path.join(work, "transcript.fasta")
        with open(fasta, "w") as fh:
            fh.write(f">{transcript.transcript_id}\n{seq}\n")
        try:
            bed = run_transdecoder(
                fasta, work, transdecoder_bin_dir=transdecoder_bin_dir
            )
        except Exception as exc:  # noqa: BLE001 - external tool; never fatal here
            _log.info("TransDecoder failed for %s: %s", transcript.transcript_id, exc)
            return None
        orf = _parse_transdecoder_bed(bed)
    if orf is None:
        return None
    tx_start, tx_end = orf
    return map_transcript_to_genomic(transcript, tx_start, tx_end)


def _parse_transdecoder_bed(bed_path: str | Path) -> tuple[int, int] | None:
    """Return ``(thickStart, thickEnd)`` of the first ``+``-strand ORF, or None.

    TransDecoder BED12 reports the ORF on the sequence we supplied; since that
    sequence is already in coding orientation we keep only ``+``-strand ORFs
    (``-`` would be antisense to the gene). ``thickStart`` is 0-based,
    ``thickEnd`` half-open — i.e. already transcript-space coordinates.
    """
    if not os.path.exists(bed_path):
        return None
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            strand = cols[5]
            if strand != "+":
                continue
            return int(cols[6]), int(cols[7])
    return None


def _parse_transdecoder_bed_multi(bed_path: str | Path) -> dict[str, tuple[int, int]]:
    """Map each input sequence id → its first ``+``-strand ORF ``(start, end)``.

    The multi-FASTA TransDecoder run emits one BED12 with up to one
    line per supplied sequence (``--single_best_only``); column 0 is the sequence
    id we wrote (the transcript id). Keeps the first ``+``-strand ORF per id (a
    ``-`` ORF would be antisense to the coding-oriented sequence we supplied).
    """
    orfs: dict[str, tuple[int, int]] = {}
    if not os.path.exists(bed_path):
        return orfs
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8 or cols[5] != "+":
                continue
            seq_id = cols[0]
            if seq_id not in orfs:
                orfs[seq_id] = (int(cols[6]), int(cols[7]))
    return orfs


def batch_backstop_transdecoder(
    genes: list[ReconciledGene],
    genome: Any,
    transdecoder_bin_dir: str | None,
    stats: Any = None,
) -> list[ReconciledGene]:
    """Rescue all CDS-less backstop genes with **one** TransDecoder invocation.

    Replaces the per-gene ``cds._transdecoder_cds`` subprocess (one spawn per
    CDS-less backstop gene — hundreds on a real genome) with a
    single multi-FASTA ``TransDecoder.LongOrfs``/``Predict`` over every CDS-less
    backstop transcript, parsing the resulting ORFs back per gene. Genes already
    carrying a CDS (e.g. miniprot-rescued) and Mikado-origin genes are untouched.

    Identical, order-independent result to the per-gene path: each transcript's
    spliced coding-orientation sequence is the same input TransDecoder would have
    seen alone, and the per-record ORF is mapped back with the same
    :func:`map_transcript_to_genomic`. Returns the (possibly) updated gene list in
    the **same order** it was given.

    ``stats`` (optional :class:`~helixforge.reconcile.runstats.RunStats`) moves a
    rescued gene from the ``backstop_rescued_none`` bucket (where the miniprot-only
    :func:`assign_backstop_cds` left it) into ``backstop_rescued_transdecoder`` so
    the telemetry matches the old single-pass semantics.
    """
    if genome is None or transdecoder_bin_dir is None:
        return list(genes)

    pending: list[tuple[int, ReconciledGene, TranscriptCandidate]] = []
    for i, gene in enumerate(genes):
        if gene.origin != "helixer_backstop":
            continue
        transcript = gene.transcripts[0]
        if transcript.cds:
            continue
        pending.append((i, gene, transcript))
    if not pending:
        return list(genes)

    with tempfile.TemporaryDirectory() as work:
        fasta = os.path.join(work, "backstop_transcripts.fasta")
        with open(fasta, "w") as fh:
            for _i, _gene, transcript in pending:
                seq = extract_transcript_sequence(transcript, genome)
                fh.write(f">{transcript.transcript_id}\n{seq}\n")
        try:
            bed = run_transdecoder(
                fasta, work, transdecoder_bin_dir=transdecoder_bin_dir
            )
        except Exception as exc:  # noqa: BLE001 - external tool; never fatal here
            _log.info(
                "batched TransDecoder failed (%d transcripts): %s", len(pending), exc
            )
            return list(genes)
        orfs = _parse_transdecoder_bed_multi(bed)

    out = list(genes)
    for i, gene, transcript in pending:
        orf = orfs.get(transcript.transcript_id)
        if orf is None:
            continue
        cds = map_transcript_to_genomic(transcript, orf[0], orf[1])
        if cds is None:
            continue
        try:
            out[i] = _rescue_backstop_cds(gene, cds)
        except ValueError as exc:
            _log.info(
                "batched TransDecoder CDS for %s rejected by validation: %s",
                gene.gene_id,
                exc,
            )
            continue
        if stats is not None:
            stats.bump("backstop_rescued_none", -1)
            stats.bump("backstop_rescued_transdecoder", 1)
    return out


def _confirm_backstop_homology(
    gene: ReconciledGene,
    transcript: TranscriptCandidate,
    alignments: list[MiniprotAlignment],
    stats: Any = None,
) -> ReconciledGene:
    """Confirm a backstop gene's intrinsic CDS against miniprot homology.

    The gene already has a (Helixer-intrinsic) CDS, so this never changes the ORF
    coordinates. If a same-locus, same-strand miniprot ORF projects onto the
    transcript and reciprocally overlaps the intrinsic CDS above
    :data:`_CROSS_CHECK_OVERLAP`, the gene is marked homology-backed (the hit
    accession is recorded as ``protein_id``) and re-tiered via
    :func:`assign_tier` — typically Tier 1. Genes already carrying homology, and
    genes with no agreeing hit, are returned unchanged.
    """
    if transcript.has_homology or transcript.cds is None:
        return gene
    hits = MiniprotParser.get_best_per_locus(
        alignments, gene.seqid, gene.start, gene.end, gene.strand
    )
    for aln in hits:
        projected = project_cds_to_exons(
            aln.cds_segments, transcript.exons, gene.strand
        )
        if projected is None:
            continue
        if reciprocal_overlap(transcript.cds, projected) < _CROSS_CHECK_OVERLAP:
            continue
        new_t = attrs.evolve(transcript, protein_id=aln.protein_id)
        tier = assign_tier(new_t, gene.classification, gene.origin)
        if stats is not None:
            stats.bump("backstop_rescued_miniprot")
        return attrs.evolve(gene, transcripts=[new_t], tier=tier)
    return gene


def assign_backstop_cds(
    gene: ReconciledGene,
    alignments: list[MiniprotAlignment],
    genome: Any = None,
    transdecoder_bin_dir: str | None = None,
    stats: Any = None,
) -> ReconciledGene:
    """Give a ``helixer_backstop`` gene a CDS; return the (possibly) updated gene.

    Tries the best same-strand miniprot alignment first (homology-backed →
    re-tier to 1), then optionally TransDecoder. Leaves ``cds=None`` (gene
    unchanged) if neither yields a valid, mod-3, in-exon CDS. Mikado-origin genes
    are returned untouched (their CDS comes from ``mikado pick``).

    ``stats`` (optional :class:`~helixforge.reconcile.runstats.RunStats`) records
    the rescue source (miniprot / TransDecoder / none) for each backstop gene;
    observation only — it never changes which source wins.
    """
    if gene.origin != "helixer_backstop":
        return gene
    transcript = gene.transcripts[0]

    # The model may already carry its intrinsic (Helixer) CDS — that ORF makes the
    # gene coding on its own. miniprot then only *confirms homology* (re-tier to 1
    # + record the hit accession); it never replaces the model's ORF coordinates.
    if transcript.cds is not None:
        return _confirm_backstop_homology(gene, transcript, alignments, stats)

    # --- 1. miniprot projection ---
    hits = MiniprotParser.get_best_per_locus(
        alignments, gene.seqid, gene.start, gene.end, gene.strand
    )
    for aln in hits:
        projected = project_cds_to_exons(
            aln.cds_segments, transcript.exons, gene.strand
        )
        if projected is None:
            continue
        # miniprot CDS is stop-exclusive; extend to include the genomic stop
        if genome is not None:
            extended = _include_genomic_stop(
                projected, transcript.exons, gene.strand,
                genome, gene.seqid,
            )
            if extended is not None:
                projected = extended
        try:
            rescued = _rescue_backstop_cds(gene, projected, protein_id=aln.protein_id)
            if stats is not None:
                stats.bump("backstop_rescued_miniprot")
            return rescued
        except ValueError as exc:
            _log.info(
                "miniprot CDS for %s rejected by model validation: %s",
                gene.gene_id,
                exc,
            )
            continue

    # --- 2. TransDecoder fallback (optional) ---
    if genome is not None and transdecoder_bin_dir is not None:
        td = _transdecoder_cds(transcript, genome, transdecoder_bin_dir)
        if td is not None:
            try:
                # CDS-only rescue (no homology accession) → Tier 2 via assign_tier.
                rescued = _rescue_backstop_cds(gene, td)
                if stats is not None:
                    stats.bump("backstop_rescued_transdecoder")
                return rescued
            except ValueError as exc:
                _log.info(
                    "TransDecoder CDS for %s rejected by model validation: %s",
                    gene.gene_id,
                    exc,
                )

    if stats is not None:
        stats.bump("backstop_rescued_none")
    return gene


# ---------------------------------------------------------------------------
# D1: cross-check (Mikado-origin genes; read-only)
# ---------------------------------------------------------------------------


def cds_cross_check(
    gene: ReconciledGene, alignments: list[MiniprotAlignment]
) -> QCFlag | None:
    """Compare a Mikado-origin gene's ORF to the best independent miniprot CDS.

    Returns ``CDS_DISAGREE`` if an overlapping same-strand miniprot ORF, once
    projected onto the gene's exons, has reciprocal CDS overlap below
    :data:`_CROSS_CHECK_OVERLAP` with the chosen primary ORF. **Never changes the
    CDS** — this is a flag only. Returns ``None`` when the gene is
    not a Mikado-origin coding gene or there is nothing to compare against.
    """
    if gene.origin not in MIKADO_ORIGINS:
        return None
    primary = next(
        (t for t in gene.transcripts if t.transcript_id == gene.primary_transcript_id),
        gene.transcripts[0],
    )
    if not primary.cds:
        return None

    hits = MiniprotParser.get_best_per_locus(
        alignments, gene.seqid, gene.start, gene.end, gene.strand
    )
    for aln in hits:
        projected = project_cds_to_exons(aln.cds_segments, primary.exons, gene.strand)
        if projected is None:
            continue
        overlap = reciprocal_overlap(primary.cds, projected)
        return CDS_DISAGREE if overlap < _CROSS_CHECK_OVERLAP else None
    return None
