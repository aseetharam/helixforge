"""All HelixForge data models."""

from __future__ import annotations

from typing import Any, Sequence

import attrs

# ---------------------------------------------------------------------------
# Allowed value sets (single source of truth for model-level enums)
# ---------------------------------------------------------------------------

VALID_STRANDS = ("+", "-")
VALID_PHASES = (0, 1, 2)

CLASSIFICATION_STATUSES = ("EXPRESSED", "LOW", "SILENT")
EVIDENCE_SOURCES = ("stringtie", "bam_coverage", "bigwig", "none")
CANDIDATE_SOURCES = ("stringtie", "helixer", "miniprot", "mikado", "ncrna")

FLAG_CATEGORIES = ("evidence", "confidence", "structure", "homology", "splice", "locus")
FLAG_SEVERITIES = ("INFO", "WARNING", "ERROR", "CRITICAL")

AS_EVENT_KINDS = ("ES", "IR", "A5", "A3", "ALT_TSS", "ALT_TES", "MX")

RECONCILED_TIERS = (1, 2, 3, 4)
RECONCILED_ORIGINS = (
    "mikado_1to1",
    "split",
    "merge",
    "helixer_backstop",
    "novel",
)

# Gene/transcript biotype taxonomy (Ensembl convention). Introduced in Phase 29
# (locus biology) for ``pseudogene``; the non-coding members (``lncRNA`` /
# ``ncRNA_undetermined``) are populated by Phase 30 (``reconcile/biotype.py``),
# which also sets the transcript-level carrier ``TranscriptCandidate.biotype``.
# The structured-ncRNA members (``tRNA`` / ``rRNA`` / ``snoRNA`` / ``snRNA`` /
# ``miRNA``) are only produced by the optional, off-by-default tRNAscan-SE /
# Infernal hook (``prep/ncrna.py``). Defined *once* here so every
# module shares the same carrier (default None = not yet classified). ``None`` is
# always allowed.
#
# ``CORE_GENE_BIOTYPES`` are the four values ``assign_biotype`` derives from the
# pipeline's own signals; ``STRUCTURED_NCRNA_BIOTYPES`` come only from the
# external structured-ncRNA hook. ``GENE_BIOTYPES`` is their union, the full set
# the models accept.
CORE_GENE_BIOTYPES = (
    "protein_coding",
    "lncRNA",
    "ncRNA_undetermined",
    "pseudogene",
)
STRUCTURED_NCRNA_BIOTYPES = (
    "tRNA",
    "rRNA",
    "snoRNA",
    "snRNA",
    "miRNA",
)
# Assigned only by the optional EDTA TE gate (off by default) to a good-ORF gene
# whose model overlaps a transposable element above the configured threshold.
TE_BIOTYPES = ("transposable_element",)
GENE_BIOTYPES = CORE_GENE_BIOTYPES + STRUCTURED_NCRNA_BIOTYPES + TE_BIOTYPES


# ---------------------------------------------------------------------------
# Shared private validators (reuse; never duplicate validation logic)
# ---------------------------------------------------------------------------


def _validate_coordinates(
    start: object, end: object, *, label: str = "coordinates"
) -> None:
    """0-based half-open: integers, ``start >= 0`` and ``start < end``."""
    if not isinstance(start, int) or isinstance(start, bool):
        raise ValueError(f"{label}: start must be an int, got {start!r}")
    if not isinstance(end, int) or isinstance(end, bool):
        raise ValueError(f"{label}: end must be an int, got {end!r}")
    if start < 0:
        raise ValueError(f"{label}: start must be >= 0, got {start}")
    if end <= start:
        raise ValueError(f"{label}: require start < end, got start={start} end={end}")


def _validate_strand(strand: object) -> None:
    if strand not in VALID_STRANDS:
        raise ValueError(f"strand must be one of {VALID_STRANDS}, got {strand!r}")


def _validate_phase(phase: object) -> None:
    if phase not in VALID_PHASES:
        raise ValueError(f"phase must be one of {VALID_PHASES}, got {phase!r}")


def _validate_fraction(value: object, name: str) -> None:
    """Value in [0, 1] (or ``None``, caller decides whether None is allowed)."""
    if value is None:
        return
    if not 0.0 <= float(value) <= 1.0:  # type: ignore[arg-type]
        raise ValueError(f"{name} must be in [0, 1], got {value}")


def _validate_intervals_sorted_nonoverlapping(
    items: Sequence[Any], *, label: str = "intervals"
) -> None:
    """Each item exposes ``.start``/``.end``; ascending, disjoint, ``start < end``."""
    prev_end = None
    for i, it in enumerate(items):
        if it.end <= it.start:
            raise ValueError(
                f"{label}[{i}]: require start < end, got {it.start}-{it.end}"
            )
        if prev_end is not None and it.start < prev_end:
            raise ValueError(
                f"{label} must be sorted ascending and non-overlapping; "
                f"item {i} start={it.start} overlaps previous end={prev_end}"
            )
        prev_end = it.end


def _validate_within_span(
    items: Sequence[Any], start: int, end: int, *, label: str = "interval"
) -> None:
    """Every item lies within the genomic span ``[start, end)``."""
    for i, it in enumerate(items):
        if it.start < start or it.end > end:
            raise ValueError(
                f"{label}[{i}] ({it.start}-{it.end}) is outside span ({start}-{end})"
            )


def _validate_exon_structure(
    exons: Sequence[Any],
    start: int,
    end: int,
    *,
    require_nonempty: bool = True,
    label: str = "exon",
) -> None:
    if require_nonempty and len(exons) == 0:
        raise ValueError(f"at least one {label} required")
    _validate_intervals_sorted_nonoverlapping(exons, label=label)
    _validate_within_span(exons, start, end, label=label)


def _validate_cds_within_exons(
    cds_segments: Sequence[Any], exons: Sequence[Any]
) -> None:
    """Every CDS segment must be contained within a single exon."""
    for i, seg in enumerate(cds_segments):
        if not any(seg.start >= ex.start and seg.end <= ex.end for ex in exons):
            raise ValueError(
                f"CDS segment[{i}] ({seg.start}-{seg.end}) is not contained within any exon"
            )


def _validate_cds_total_mod3(
    cds_segments: Sequence[Any], partial: bool = False
) -> None:
    """Total CDS length divisible by 3.

    The mod-3 invariant holds for **complete** ORFs (start + stop present). Real
    Mikado/TransDecoder output also contains *partial* ORFs (5'/3'-partial, at
    contig edges or genuinely incomplete genes) whose total CDS length need not
    be a multiple of 3, reading frame is preserved by the GFF3 ``phase`` on the
    5'-most coding segment. When ``partial`` is set, the mod-3 check is skipped;
    all other CDS invariants (within-exon, sorted, non-overlapping) still apply.
    """
    if partial:
        return
    total = sum(len(seg) for seg in cds_segments)
    if total % 3 != 0:
        raise ValueError(f"total CDS length must be divisible by 3, got {total}")


# ---------------------------------------------------------------------------
# Primitive coordinate models
# ---------------------------------------------------------------------------


@attrs.define
class Interval:
    """A generic 0-based half-open genomic interval ``[start, end)``."""

    start: int
    end: int

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="Interval")

    def __len__(self) -> int:
        return self.end - self.start


@attrs.define
class Exon:
    """A single exon, 0-based half-open."""

    start: int
    end: int

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="Exon")

    def __len__(self) -> int:
        return self.end - self.start


@attrs.define
class CDSSegment:
    """A CDS segment with reading-frame phase ∈ {0, 1, 2}."""

    start: int
    end: int
    phase: int = 0

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="CDSSegment")
        _validate_phase(self.phase)

    def __len__(self) -> int:
        return self.end - self.start


# --- splice-site motif classes ---
# A junction's donor/acceptor dinucleotides classify (strand-aware) into one
# canonical class or 'non-canonical'. ``None`` means the motif was never
# evaluated (no genome / no STAR motif column), distinct from 'non-canonical'.
# AT-AC is the minor (U12) spliceosome motif and counts as canonical here.
SPLICE_GT_AG = "GT-AG"
SPLICE_GC_AG = "GC-AG"
SPLICE_AT_AC = "AT-AC"
SPLICE_NON_CANONICAL = "non-canonical"
CANONICAL_SPLICE_MOTIFS = frozenset({SPLICE_GT_AG, SPLICE_GC_AG, SPLICE_AT_AC})
VALID_SPLICE_MOTIFS = CANONICAL_SPLICE_MOTIFS | {SPLICE_NON_CANONICAL}


@attrs.define
class SpliceJunction:
    """An intron defined by donor/acceptor genomic positions (low, high)."""

    seqid: str
    donor: int
    acceptor: int
    strand: str
    read_count: int
    samples: int = 1
    # Splice-site motif class (one of VALID_SPLICE_MOTIFS) or None
    # when never evaluated. Default None preserves all existing construction.
    canonical: str | None = None
    # STAR multi-mapping read count (col 7 of SJ.out.tab); surfaced
    # so a junction supported only by multi-mappers stays visible.
    multimap_reads: int = 0

    def __attrs_post_init__(self) -> None:
        # donor < acceptor (stored low->high regardless of strand)
        _validate_coordinates(self.donor, self.acceptor, label="SpliceJunction")
        _validate_strand(self.strand)
        if self.read_count < 0:
            raise ValueError(f"read_count must be >= 0, got {self.read_count}")
        if self.samples < 1:
            raise ValueError(f"samples must be >= 1, got {self.samples}")
        if self.canonical is not None and self.canonical not in VALID_SPLICE_MOTIFS:
            raise ValueError(
                f"canonical must be one of {sorted(VALID_SPLICE_MOTIFS)} or None, "
                f"got {self.canonical!r}"
            )
        if self.multimap_reads < 0:
            raise ValueError(f"multimap_reads must be >= 0, got {self.multimap_reads}")

    @property
    def intron_length(self) -> int:
        return self.acceptor - self.donor

    @property
    def is_canonical(self) -> bool:
        """True iff the motif was evaluated and is a canonical class.

        ``None`` (never evaluated) returns False, callers that must treat an
        unevaluated junction as "not known non-canonical" should test
        ``canonical in CANONICAL_SPLICE_MOTIFS or canonical is None``.
        """
        return self.canonical in CANONICAL_SPLICE_MOTIFS


# ---------------------------------------------------------------------------
# Evidence / source models
# ---------------------------------------------------------------------------


@attrs.define
class HelixerLocus:
    """A Helixer gene locus (one transcript per gene by design)."""

    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    confidence: float | None = None
    exons: list[Exon] = attrs.field(factory=list)
    cds: list[CDSSegment] | None = None

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="HelixerLocus")
        _validate_strand(self.strand)
        _validate_fraction(self.confidence, "confidence")
        _validate_exon_structure(
            self.exons, self.start, self.end, require_nonempty=False
        )
        if self.cds is not None:
            _validate_intervals_sorted_nonoverlapping(self.cds, label="CDS")
            _validate_cds_within_exons(self.cds, self.exons)

    @property
    def span(self) -> Interval:
        return Interval(self.start, self.end)


@attrs.define
class StringTieTranscript:
    """An assembled transcript from StringTie with TPM/coverage."""

    transcript_id: str
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    exons: list[Exon]
    tpm: float
    sample_id: str
    coverage: float | None = None

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="StringTieTranscript")
        _validate_strand(self.strand)
        _validate_exon_structure(
            self.exons, self.start, self.end, require_nonempty=True
        )
        if self.tpm < 0:
            raise ValueError(f"tpm must be >= 0, got {self.tpm}")
        if self.coverage is not None and self.coverage < 0:
            raise ValueError(f"coverage must be >= 0, got {self.coverage}")


@attrs.define
class MiniprotAlignment:
    """A protein-to-genome alignment from miniprot."""

    protein_id: str
    seqid: str
    start: int
    end: int
    strand: str
    cds_segments: list[CDSSegment]
    query_coverage: float
    identity: float
    score: float
    rank: int = 0

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="MiniprotAlignment")
        _validate_strand(self.strand)
        if len(self.cds_segments) == 0:
            raise ValueError("at least one CDS segment required")
        _validate_intervals_sorted_nonoverlapping(self.cds_segments, label="CDS")
        _validate_within_span(self.cds_segments, self.start, self.end, label="CDS")
        _validate_fraction(self.query_coverage, "query_coverage")
        _validate_fraction(self.identity, "identity")
        if self.query_coverage is None or self.identity is None:
            raise ValueError("query_coverage and identity are required (not None)")


# ---------------------------------------------------------------------------
# Classification + candidate models
# ---------------------------------------------------------------------------


@attrs.define
class LocusClassification:
    """Expression-based classification of a Helixer locus (Phase 3)."""

    locus_id: str
    status: str
    max_tpm: float | None = None
    mean_coverage: float | None = None
    evidence_source: str = "none"
    num_samples_expressed: int = 0

    def __attrs_post_init__(self) -> None:
        if self.status not in CLASSIFICATION_STATUSES:
            raise ValueError(
                f"status must be one of {CLASSIFICATION_STATUSES}, got {self.status!r}"
            )
        if self.evidence_source not in EVIDENCE_SOURCES:
            raise ValueError(
                f"evidence_source must be one of {EVIDENCE_SOURCES}, "
                f"got {self.evidence_source!r}"
            )
        if self.num_samples_expressed < 0:
            raise ValueError(
                f"num_samples_expressed must be >= 0, got {self.num_samples_expressed}"
            )


@attrs.define
class TranscriptCandidate:
    """A candidate transcript from any source, anchored to a locus."""

    transcript_id: str
    locus_id: str
    source: str
    seqid: str
    start: int
    end: int
    strand: str
    exons: list[Exon]
    cds: list[CDSSegment] | None = None
    # Partial-ORF granularity: 5'-partial means no start codon,
    # 3'-partial means no stop codon. They are independent, a transcript can be
    # 5'-partial-but-3'-complete, in which case the stop is still verifiable. The
    # mod-3 exemption keys on *either* end being partial; the codon gate keys on
    # the matching end (start↔5', stop↔3'). ``cds_partial`` is the derived OR.
    cds_partial_5prime: bool = False
    cds_partial_3prime: bool = False
    tpm: float | None = None
    protein_id: str | None = None
    blast_score: float | None = None
    junction_support_fraction: float | None = None
    confidence: float | None = None
    combined_score: float | None = None
    is_primary: bool = False
    # Transcript-level biotype (Ensembl ``transcript_biotype``), one of
    # GENE_BIOTYPES or None. The carrier deferred from Phase 29, set by Phase 30's
    # ``assign_biotype`` (mirrors the gene biotype onto the transcripts) so the
    # GFF3/GTF can emit ``transcript_biotype`` per mRNA. None keeps every existing
    # construction valid and the golden output unchanged until a phase sets it.
    biotype: str | None = None
    # TRaCE election rank (Phase 33b): the 1-based position of this transcript in
    # the Transcript Ranking and Canonical Election ordering (1 = elected
    # canonical/primary). ``None`` means TRaCE was not run (the default, off-by-
    # default path), the primary is then the highest-``combined_score`` isoform.
    # An optional carrier only: a None default leaves every existing construction
    # and the golden output byte-for-byte unchanged.
    trace_rank: int | None = None
    # Backward-compat init alias for the old single ``cds_partial`` bool. When
    # supplied (not None) it sets *both* granular flags; ``cds_partial`` itself is
    # now a read-only property (= 5'-partial OR 3'-partial). Not stored long-term:
    # reset to None in post-init so it never clobbers a later ``attrs.evolve`` of
    # the granular flags (preserves attrs.evolve immutable-update semantics).
    _cds_partial_compat: bool | None = attrs.field(
        default=None, alias="cds_partial", repr=False, eq=False
    )

    def __attrs_post_init__(self) -> None:
        if self._cds_partial_compat is not None:
            both = bool(self._cds_partial_compat)
            self.cds_partial_5prime = both
            self.cds_partial_3prime = both
            self._cds_partial_compat = None
        if self.source not in CANDIDATE_SOURCES:
            raise ValueError(
                f"source must be one of {CANDIDATE_SOURCES}, got {self.source!r}"
            )
        if self.biotype is not None and self.biotype not in GENE_BIOTYPES:
            raise ValueError(
                f"biotype must be one of {GENE_BIOTYPES} or None, got {self.biotype!r}"
            )
        if self.trace_rank is not None and (
            not isinstance(self.trace_rank, int)
            or isinstance(self.trace_rank, bool)
            or self.trace_rank < 1
        ):
            raise ValueError(
                f"trace_rank must be a positive int or None, got {self.trace_rank!r}"
            )
        _validate_coordinates(self.start, self.end, label="TranscriptCandidate")
        _validate_strand(self.strand)
        _validate_exon_structure(
            self.exons, self.start, self.end, require_nonempty=True
        )
        if self.cds is not None:
            _validate_intervals_sorted_nonoverlapping(self.cds, label="CDS")
            _validate_cds_within_exons(self.cds, self.exons)
            # Partial ORFs: mod-3 enforced only for complete CDS.
            # Either end being partial exempts the whole CDS from mod-3.
            _validate_cds_total_mod3(self.cds, partial=self.cds_partial)
        _validate_fraction(self.junction_support_fraction, "junction_support_fraction")
        _validate_fraction(self.confidence, "confidence")

    # --- computed structural properties ---
    @property
    def cds_partial(self) -> bool:
        """True if the CDS is partial at either end (no start and/or no stop).

        Derived OR of :attr:`cds_partial_5prime` / :attr:`cds_partial_3prime`,
        kept for backward compatibility with the pre-granularity call sites."""
        return self.cds_partial_5prime or self.cds_partial_3prime

    @property
    def num_exons(self) -> int:
        return len(self.exons)

    @property
    def num_introns(self) -> int:
        return max(0, len(self.exons) - 1)

    @property
    def total_exon_length(self) -> int:
        return sum(len(e) for e in self.exons)

    @property
    def total_cds_length(self) -> int:
        if self.cds is None:
            return 0
        return sum(len(c) for c in self.cds)

    @property
    def introns(self) -> list[Interval]:
        return [
            Interval(self.exons[i].end, self.exons[i + 1].start)
            for i in range(len(self.exons) - 1)
        ]

    @property
    def has_homology(self) -> bool:
        """True if the ORF is protein-homology-backed.

        Either a hit accession (``protein_id``) or a positive BLAST/DIAMOND score
        (``blast_score``). Mikado's loci GFF3 carries no accession, so in practice
        the BLAST score from the metrics TSV is the homology signal (used for the
        Tier-1 gate)."""
        return self.protein_id is not None or (self.blast_score or 0) > 0


# ---------------------------------------------------------------------------
# QC flag record (frozen constant)
# ---------------------------------------------------------------------------


@attrs.frozen
class QCFlag:
    """A QC flag definition. Constants live in ``qc/flags.py``."""

    name: str
    category: str
    severity: str
    description: str = ""

    def __attrs_post_init__(self) -> None:
        if self.category not in FLAG_CATEGORIES:
            raise ValueError(
                f"category must be one of {FLAG_CATEGORIES}, got {self.category!r}"
            )
        if self.severity not in FLAG_SEVERITIES:
            raise ValueError(
                f"severity must be one of {FLAG_SEVERITIES}, got {self.severity!r}"
            )


# ---------------------------------------------------------------------------
# v3-new models
# ---------------------------------------------------------------------------


@attrs.frozen
class ASEvent:
    """An alternative-splicing event (frozen record).

    kinds: ES (exon skip), IR (intron retention), A5/A3 (alt 5'/3' splice
    site), ALT_TSS / ALT_TES (alt transcription start/end), MX (mutually
    exclusive exons).
    """

    kind: str
    seqid: str
    start: int
    end: int
    strand: str
    support_read_count: int = 0

    def __attrs_post_init__(self) -> None:
        if self.kind not in AS_EVENT_KINDS:
            raise ValueError(f"kind must be one of {AS_EVENT_KINDS}, got {self.kind!r}")
        _validate_coordinates(self.start, self.end, label="ASEvent")
        _validate_strand(self.strand)
        if self.support_read_count < 0:
            raise ValueError(
                f"support_read_count must be >= 0, got {self.support_read_count}"
            )


@attrs.define
class MikadoLocus:
    """A locus parsed from ``mikado.loci.gff3`` + metrics/scores TSVs (Phase 5)."""

    locus_id: str
    seqid: str
    start: int
    end: int
    strand: str
    transcripts: list[TranscriptCandidate]
    metrics: dict[str, object] = attrs.field(factory=dict)
    scores: dict[str, object] = attrs.field(factory=dict)

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="MikadoLocus")
        _validate_strand(self.strand)
        for t in self.transcripts:
            if t.source != "mikado":
                raise ValueError(
                    f"MikadoLocus transcripts must have source='mikado', "
                    f"got {t.source!r} for {t.transcript_id!r}"
                )
            if t.seqid != self.seqid:
                raise ValueError(
                    f"transcript {t.transcript_id!r} seqid {t.seqid!r} != locus seqid "
                    f"{self.seqid!r}"
                )
            if t.strand != self.strand:
                raise ValueError(
                    f"transcript {t.transcript_id!r} strand {t.strand!r} != locus strand "
                    f"{self.strand!r}"
                )


@attrs.define
class ReconciledGene:
    """The final reconciled gene (supersedes v2 ``GeneResult``)."""

    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    tier: int
    transcripts: list[TranscriptCandidate]
    primary_transcript_id: str
    classification: LocusClassification
    origin: str
    as_events: list[ASEvent] = attrs.field(factory=list)
    flags: list[QCFlag] = attrs.field(factory=list)
    merged_from: list[str] = attrs.field(factory=list)
    # Gene biotype (Ensembl convention), one of GENE_BIOTYPES or None when not
    # yet classified. Introduced in Phase 29 (set to ``pseudogene`` for a
    # homology-backed disabled ORF); the full protein_coding/lncRNA/ncRNA
    # taxonomy is assigned in Phase 30. None default keeps every existing
    # construction valid and the golden output unchanged until a phase sets it.
    biotype: str | None = None

    def __attrs_post_init__(self) -> None:
        _validate_coordinates(self.start, self.end, label="ReconciledGene")
        _validate_strand(self.strand)
        if self.tier not in RECONCILED_TIERS:
            raise ValueError(
                f"tier must be one of {RECONCILED_TIERS}, got {self.tier!r}"
            )
        if self.biotype is not None and self.biotype not in GENE_BIOTYPES:
            raise ValueError(
                f"biotype must be one of {GENE_BIOTYPES} or None, got {self.biotype!r}"
            )
        if self.origin not in RECONCILED_ORIGINS:
            raise ValueError(
                f"origin must be one of {RECONCILED_ORIGINS}, got {self.origin!r}"
            )
        if len(self.transcripts) == 0:
            raise ValueError("ReconciledGene requires at least one transcript")
        tids = [t.transcript_id for t in self.transcripts]
        if self.primary_transcript_id not in tids:
            raise ValueError(
                f"primary_transcript_id {self.primary_transcript_id!r} not among "
                f"transcripts {tids}"
            )
        for t in self.transcripts:
            if t.seqid != self.seqid:
                raise ValueError(
                    f"transcript {t.transcript_id!r} seqid {t.seqid!r} != gene seqid "
                    f"{self.seqid!r}"
                )
            if t.strand != self.strand:
                raise ValueError(
                    f"transcript {t.transcript_id!r} strand {t.strand!r} != gene strand "
                    f"{self.strand!r}"
                )


@attrs.define
class IsoformAdmission:
    """Audit record for an isoform admission decision (drives Phase 9 stats)."""

    transcript_id: str
    gene_id: str
    admitted: bool
    reason: str
    novel_as_event: ASEvent | None = None
    score: float | None = None
    redundant_with: str | None = None

    def __attrs_post_init__(self) -> None:
        if not isinstance(self.admitted, bool):
            raise ValueError(f"admitted must be a bool, got {self.admitted!r}")
        if not isinstance(self.reason, str):
            raise ValueError(f"reason must be a str, got {self.reason!r}")
        if self.novel_as_event is not None and not isinstance(
            self.novel_as_event, ASEvent
        ):
            raise ValueError("novel_as_event must be an ASEvent or None")
