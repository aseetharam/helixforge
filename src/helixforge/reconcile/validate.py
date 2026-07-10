"""Structural codon gate + QC flags (read-only)."""

from __future__ import annotations

from typing import TYPE_CHECKING

from helixforge.constants import DEFAULT_TRANSL_TABLE
from helixforge.qc.flags import (
    AMBIGUOUS_CODON,
    INTERNAL_STOP,
    LONG_INTRON,
    NO_START,
    NO_STOP,
    PARTIAL_ORF,
    SHORT_CDS,
    SHORT_EXON,
)
from helixforge.utils.logging import get_logger
from helixforge.utils.sequences import (
    has_ambiguous_base,
    is_start_codon,
    is_stop_codon,
)
from helixforge.utils.sequences import check_internal_stops as _seq_internal_stops

if TYPE_CHECKING:
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.reconcile.models import (
        CDSSegment,
        QCFlag,
        ReconciledGene,
        TranscriptCandidate,
    )

_log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Sequence extraction (coding direction, both strands)
# ---------------------------------------------------------------------------


def _safe_get_sequence(
    genome: GenomeAccessor,
    seqid: str,
    start: int,
    end: int,
    strand: str,
) -> str | None:
    """``genome.get_sequence`` guarded against off-contig / bad coordinates.

    Returns ``None`` (rather than raising) when the requested window runs off
    the start of the contig, has non-positive width, or the backend rejects it
    — the caller then skips that sequence-dependent check instead of crashing on
    an edge-of-contig gene.
    """
    if start < 0 or end <= start:
        return None
    try:
        return genome.get_sequence(seqid, start, end, strand)
    except (KeyError, IndexError, ValueError):
        return None


def extract_cds_sequence(
    transcript: TranscriptCandidate,
    genome: GenomeAccessor,
) -> str | None:
    """Return the CDS nucleotide sequence in coding (5'→3') orientation.

    CDS segments are stored low→high genomically. On ``+`` the
    coding order is ascending; on ``-`` it is descending and each piece is
    reverse-complemented (``RC(a+b) == RC(b)+RC(a)``). Returns ``None`` if the
    transcript has no CDS or any segment cannot be read.
    """
    if not transcript.cds:
        return None
    segs = list(transcript.cds)  # already sorted ascending by model invariants
    if transcript.strand == "-":
        segs = list(reversed(segs))
    pieces: list[str] = []
    for seg in segs:
        piece = _safe_get_sequence(
            genome, transcript.seqid, seg.start, seg.end, transcript.strand
        )
        if piece is None:
            return None
        pieces.append(piece)
    return "".join(pieces)


def _coding_first_segment(transcript: TranscriptCandidate) -> CDSSegment | None:
    """The 5'-most CDS segment in coding direction (``+`` low, ``-`` high)."""
    if not transcript.cds:
        return None
    return transcript.cds[0] if transcript.strand == "+" else transcript.cds[-1]


# ---------------------------------------------------------------------------
# Per-rule checks (each returns a QCFlag or None; never mutates)
# ---------------------------------------------------------------------------


def check_start_codon(
    transcript: TranscriptCandidate,
    genome: GenomeAccessor | None,
    transl_table: int = DEFAULT_TRANSL_TABLE,
) -> QCFlag | None:
    """``NO_START`` if the CDS does not begin with ATG.

    Skipped only for **5'-partial** ORFs (a 5'-partial CDS legitimately has no
    start) and when the codon cannot be read. A 3'-partial-but-5'-complete CDS
    still gets its start verified.

    If the start codon overlaps an ``N``/IUPAC ambiguous base the codon is
    *indeterminate* → ``AMBIGUOUS_CODON`` (INFO), never ``NO_START``: a masked
    assembly gap must not masquerade as a broken start.
    """
    if not transcript.cds or genome is None or transcript.cds_partial_5prime:
        return None
    first = transcript.cds[0]
    last = transcript.cds[-1]
    if transcript.strand == "+":
        codon = _safe_get_sequence(
            genome, transcript.seqid, first.start, first.start + 3, "+"
        )
    else:
        codon = _safe_get_sequence(
            genome, transcript.seqid, last.end - 3, last.end, "-"
        )
    if codon is None or len(codon) != 3:
        return None
    if has_ambiguous_base(codon):
        return AMBIGUOUS_CODON
    return None if is_start_codon(codon) else NO_START


def check_stop_codon(
    transcript: TranscriptCandidate,
    genome: GenomeAccessor | None,
    transl_table: int = DEFAULT_TRANSL_TABLE,
) -> QCFlag | None:
    """``NO_STOP`` if the terminal CDS codon is not a stop (§5.2, stop-inclusive).

    The CDS includes the stop codon as its last 3 bases (Helixer/Ensembl
    convention). This checks the last codon of the CDS in coding direction.

    Skipped only for **3'-partial** ORFs (a 3'-partial CDS legitimately has no
    stop) and when the codon cannot be read. A 5'-partial-but-3'-complete CDS
    still gets its stop verified.

    A stop codon overlapping an ``N``/IUPAC base is indeterminate →
    ``AMBIGUOUS_CODON`` (INFO), never ``NO_STOP``. The stop set is
    table-aware (``transl_table``) so organellar contigs use the right code.
    """
    if not transcript.cds or genome is None or transcript.cds_partial_3prime:
        return None
    first = transcript.cds[0]
    last = transcript.cds[-1]
    if transcript.strand == "+":
        codon = _safe_get_sequence(
            genome, transcript.seqid, last.end - 3, last.end, "+"
        )
    else:
        codon = _safe_get_sequence(
            genome, transcript.seqid, first.start, first.start + 3, "-"
        )
    if codon is None or len(codon) != 3:
        return None
    if has_ambiguous_base(codon):
        return AMBIGUOUS_CODON
    return None if is_stop_codon(codon, transl_table) else NO_STOP


def check_internal_stops(
    transcript: TranscriptCandidate,
    genome: GenomeAccessor | None,
    transl_table: int = DEFAULT_TRANSL_TABLE,
) -> QCFlag | None:
    """``INTERNAL_STOP`` if the CDS translation has a premature stop (§5.5).

    An ``N``/IUPAC internal codon is "unknown" (translates to ``X``), not a
    premature stop — so a masked base no longer silently suppresses
    *or* fabricates the check. ``transl_table`` selects the genetic code.
    """
    if not transcript.cds or genome is None:
        return None
    seq = extract_cds_sequence(transcript, genome)
    if seq is None:
        return None
    first = _coding_first_segment(transcript)
    phase = first.phase if first is not None else 0
    try:
        positions = _seq_internal_stops(seq, phase, transl_table)
    except ValueError:
        # Truly non-IUPAC character in the extracted window — cannot evaluate.
        return None
    return INTERNAL_STOP if positions else None


def check_cds_mod3(transcript: TranscriptCandidate) -> bool:
    """Defense in depth: total CDS length divisible by 3 for **complete** ORFs.

    The model already enforces this at construction; this
    helper re-asserts it so a CDS built by a path that bypassed validation is
    still caught. Returns ``True`` when the invariant holds (or the CDS is
    partial / absent), ``False`` otherwise. Emits no flag of its own.
    """
    if not transcript.cds or transcript.cds_partial:
        return True
    return transcript.total_cds_length % 3 == 0


def check_cds_within_exons(transcript: TranscriptCandidate) -> bool:
    """Defense in depth: every CDS segment contained within an exon (§4.7).

    Enforced at construction; re-asserted here. Returns ``True``/``False``;
    emits no flag of its own.
    """
    if not transcript.cds:
        return True
    return all(
        any(seg.start >= ex.start and seg.end <= ex.end for ex in transcript.exons)
        for seg in transcript.cds
    )


def check_short_cds(
    transcript: TranscriptCandidate, threshold: int = 300
) -> QCFlag | None:
    """``SHORT_CDS`` if the CDS is shorter than ``threshold`` nucleotides."""
    if not transcript.cds:
        return None
    return SHORT_CDS if transcript.total_cds_length < threshold else None


def check_short_exons(
    transcript: TranscriptCandidate, threshold: int = 10
) -> QCFlag | None:
    """``SHORT_EXON`` if any exon is shorter than ``threshold`` nucleotides."""
    return SHORT_EXON if any(len(ex) < threshold for ex in transcript.exons) else None


def check_long_introns(
    transcript: TranscriptCandidate, threshold: int = 100_000
) -> QCFlag | None:
    """``LONG_INTRON`` if any intron is longer than ``threshold`` nucleotides."""
    return LONG_INTRON if any(len(i) > threshold for i in transcript.introns) else None


# ---------------------------------------------------------------------------
# Gene-level entry points
# ---------------------------------------------------------------------------


def _primary_transcript(gene: ReconciledGene) -> TranscriptCandidate:
    for t in gene.transcripts:
        if t.transcript_id == gene.primary_transcript_id:
            return t
    return gene.transcripts[0]


def _resolve_table(
    seqid: str,
    transl_table: int,
    transl_table_map: dict[str, int] | None,
) -> int:
    """Per-scaffold genetic code: ``transl_table_map[seqid]`` else the default.

    Lets organellar (plastid/mito) contigs use the right NCBI code while every
    unmapped (nuclear) scaffold keeps the run default.
    """
    if transl_table_map:
        return transl_table_map.get(seqid, transl_table)
    return transl_table


def validate_gene(
    gene: ReconciledGene,
    genome: GenomeAccessor | None = None,
    short_cds_threshold: int = 300,
    short_exon_threshold: int = 10,
    long_intron_threshold: int = 100_000,
    transl_table: int = DEFAULT_TRANSL_TABLE,
    transl_table_map: dict[str, int] | None = None,
) -> list[QCFlag]:
    """Validate a gene's primary transcript; return the list of ``QCFlag``s.

    Read-only. Sequence checks need ``genome`` (skipped if ``None``); CDS checks
    need a CDS (skipped if absent). The genetic code is ``transl_table``,
    overridable per-scaffold via ``transl_table_map`` (``{seqid: table_id}``) so
    organellar contigs translate correctly. Flags are de-duplicated by name.
    """
    t = _primary_transcript(gene)
    table = _resolve_table(t.seqid, transl_table, transl_table_map)
    flags: list[QCFlag] = []

    if t.cds and t.cds_partial:
        flags.append(PARTIAL_ORF)

    for flag in (
        check_start_codon(t, genome, table),
        check_stop_codon(t, genome, table),
        check_internal_stops(t, genome, table),
        check_short_cds(t, short_cds_threshold),
        check_short_exons(t, short_exon_threshold),
        check_long_introns(t, long_intron_threshold),
    ):
        if flag is not None:
            flags.append(flag)

    # Defense in depth — should never fire given model invariants, but a CDS
    # built outside the models would be caught here.
    if not check_cds_mod3(t):
        _log.warning(
            "gene %s: CDS mod-3 invariant violated post-construction", gene.gene_id
        )
    if not check_cds_within_exons(t):
        _log.warning(
            "gene %s: CDS-within-exons invariant violated post-construction",
            gene.gene_id,
        )

    seen: set[str] = set()
    out: list[QCFlag] = []
    for f in flags:
        if f.name not in seen:
            seen.add(f.name)
            out.append(f)
    return out


def validate_all(
    genes: list[ReconciledGene],
    genome: GenomeAccessor | None = None,
    short_cds_threshold: int = 300,
    short_exon_threshold: int = 10,
    long_intron_threshold: int = 100_000,
    transl_table: int = DEFAULT_TRANSL_TABLE,
    transl_table_map: dict[str, int] | None = None,
) -> dict[str, list[QCFlag]]:
    """Validate every gene; return ``{gene_id: [QCFlag, ...]}``."""
    return {
        g.gene_id: validate_gene(
            g,
            genome=genome,
            short_cds_threshold=short_cds_threshold,
            short_exon_threshold=short_exon_threshold,
            long_intron_threshold=long_intron_threshold,
            transl_table=transl_table,
            transl_table_map=transl_table_map,
        )
        for g in genes
    }
