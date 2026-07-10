"""BAM / STAR evidence I/O."""

from __future__ import annotations

import os
from contextlib import ExitStack
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import pysam

if TYPE_CHECKING:
    import numpy as np

from helixforge.constants import COVERAGE_MAX_DEPTH, UNIQUE_MIN_MAPQ
from helixforge.io.validate import FormatError, _star_int
from helixforge.reconcile.models import (
    CANONICAL_SPLICE_MOTIFS,
    SPLICE_AT_AC,
    SPLICE_GC_AG,
    SPLICE_GT_AG,
    SPLICE_NON_CANONICAL,
    SpliceJunction,
)
from helixforge.utils.logging import get_logger
from helixforge.utils.sequences import reverse_complement

_log = get_logger(__name__)

# CIGAR op codes (pysam): M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8
_REF_CONSUMING = frozenset({0, 2, 3, 7, 8})
_MATCH_OPS = frozenset({0, 7, 8})
_CIGAR_N = 3

# STAR SJ.out.tab strand codes: 0=undefined, 1='+', 2='-'.
_STAR_STRAND_CODES_LOCAL = frozenset({"0", "1", "2"})

# STAR SJ.out.tab intron-motif codes (col 4) -> canonical class.
# 0=non-canonical; 1/2=GT/AG (+/-); 3/4=GC/AG (+/-); 5/6=AT/AC (+/-).
# A free, high-confidence signal STAR already computed against the genome.
_STAR_MOTIF_CANONICAL: dict[str, str] = {
    "0": SPLICE_NON_CANONICAL,
    "1": SPLICE_GT_AG,
    "2": SPLICE_GT_AG,
    "3": SPLICE_GC_AG,
    "4": SPLICE_GC_AG,
    "5": SPLICE_AT_AC,
    "6": SPLICE_AT_AC,
}

# Junction-strand sourcing policy:
#   "xs"            — strand from the XS tag only; drop reads without XS (legacy).
#   "motif"         — strand from the canonical splice motif only (needs a genome).
#   "xs_then_motif" — XS when present, else infer from the motif (default). With
#                     no genome supplied this degrades to legacy "xs" behaviour.
STRAND_SOURCES = frozenset({"xs", "motif", "xs_then_motif"})
DEFAULT_STRAND_SOURCE = "xs_then_motif"


def classify_splice_motif(donor_dinuc: str, acceptor_dinuc: str, strand: str) -> str:
    """Classify a splice junction motif from its donor/acceptor dinucleotides.

    ``donor_dinuc`` is the genomic-forward sequence at ``[donor, donor + 2)``
    (the intron's low-coordinate end) and ``acceptor_dinuc`` is ``[acceptor - 2,
    acceptor)`` (the high-coordinate end) — both read on the **forward** genomic
    strand regardless of ``strand``. The classification is strand-aware: on the
    minus strand the coding-orientation donor/acceptor are the reverse
    complements of the genomic acceptor/donor. Returns one of ``'GT-AG'``,
    ``'GC-AG'``, ``'AT-AC'`` (canonical) or ``'non-canonical'``. Never raises:
    an ``N``/IUPAC or short dinucleotide simply classifies as non-canonical.
    """
    if strand == "-":
        # Coding 5' donor sits at the genomic high end; coding 3' acceptor at the
        # low end. Reverse-complement both so the canonical table below applies.
        cod_donor = reverse_complement(acceptor_dinuc)
        cod_acceptor = reverse_complement(donor_dinuc)
    else:
        cod_donor = donor_dinuc.upper()
        cod_acceptor = acceptor_dinuc.upper()
    pair = (cod_donor, cod_acceptor)
    if pair == ("GT", "AG"):
        return SPLICE_GT_AG
    if pair == ("GC", "AG"):
        return SPLICE_GC_AG
    if pair == ("AT", "AC"):
        return SPLICE_AT_AC
    return SPLICE_NON_CANONICAL


def infer_strand_from_motif(donor_dinuc: str, acceptor_dinuc: str) -> str | None:
    """Infer junction strand from the canonical splice-site motif.

    A canonical motif read in the forward orientation implies ``'+'``; one
    canonical only when read in the reverse orientation implies ``'-'``. Returns
    ``None`` when the motif is non-canonical in both orientations (cannot orient
    — the caller then drops the junction unless XS provides the strand).
    """
    if (
        classify_splice_motif(donor_dinuc, acceptor_dinuc, "+")
        in CANONICAL_SPLICE_MOTIFS
    ):
        return "+"
    if (
        classify_splice_motif(donor_dinuc, acceptor_dinuc, "-")
        in CANONICAL_SPLICE_MOTIFS
    ):
        return "-"
    return None


# A .bai index uses 16-bit (R-tree) bins capped at 2^29 - 1 bp; any contig
# longer than 2^29 bp cannot be addressed by .bai and a pysam fetch over it
# *silently returns nothing*. Wheat/barley/maize
# chromosomes exceed this — the .csi index is mandatory for them.
CSI_REQUIRED_THRESHOLD = 2**29  # 536_870_912 bp (~512 Mb)


def _index_path_exists(bam_path: str) -> bool:
    return (
        os.path.exists(bam_path + ".bai")
        or os.path.exists(bam_path + ".csi")
        or os.path.exists(bam_path + ".crai")
        or os.path.exists(os.path.splitext(bam_path)[0] + ".bai")
    )


def _has_csi_index(bam_path: str) -> bool:
    """True if a ``.csi`` index exists for ``bam_path`` (``<bam>.csi`` form)."""
    return os.path.exists(str(bam_path) + ".csi") or os.path.exists(
        os.path.splitext(str(bam_path))[0] + ".csi"
    )


# --- CRAM support -----------------------------------------------------------
# CRAM is BAM's reference-compressed sibling: alignments store only the diff
# against the genome, so the FASTA must be supplied (``reference_filename``) to
# decode them. We accept it transparently — detected by extension or the 4-byte
# ``CRAM`` magic — so junctions/coverage work on a CRAM exactly as on a BAM.


def _is_cram(path: str | os.PathLike[str]) -> bool:
    """True if ``path`` is a CRAM (``.cram`` extension or the ``CRAM`` magic)."""
    if str(path).lower().endswith(".cram"):
        return True
    try:
        with open(path, "rb") as fh:
            return fh.read(4) == b"CRAM"
    except OSError:
        return False


def _has_crai_index(path: str | os.PathLike[str]) -> bool:
    """True if a ``.crai`` index exists for a CRAM (``<cram>.crai`` form)."""
    return os.path.exists(str(path) + ".crai") or os.path.exists(
        os.path.splitext(str(path))[0] + ".crai"
    )


def _open_alignment(
    path: str | os.PathLike[str],
    *,
    reference_filename: str | os.PathLike[str] | None = None,
) -> pysam.AlignmentFile:
    """Open a BAM **or** CRAM, picking the mode + threading the reference.

    A CRAM (``.cram`` / ``CRAM`` magic) is opened ``"rc"`` with the genome FASTA
    as ``reference_filename`` so its reference-compressed records decode; a BAM is
    opened ``"rb"`` exactly as before. ``reference_filename`` is ignored for BAM
    (and may be ``None`` for a CRAM that embeds/locates its own reference).
    """
    if _is_cram(path):
        ref = str(reference_filename) if reference_filename is not None else None
        return pysam.AlignmentFile(str(path), "rc", reference_filename=ref)
    return pysam.AlignmentFile(str(path), "rb")


def _pileup_kwargs(handle: pysam.AlignmentFile) -> dict[str, Any]:
    """Extra ``pileup()`` kwargs that depend on the file type.

    ``AlignmentFile.pileup`` defaults ``multiple_iterators=True``, which htslib
    does **not** implement for CRAM — pysam then emits a ``UserWarning`` and
    silently falls back to a single iterator. We pre-empt that for CRAM by
    passing ``multiple_iterators=False`` explicitly (the same effective
    behaviour, no warning). A BAM is left untouched so its pileup semantics are
    byte-for-byte unchanged (the audit's ``io/bam.py:612`` item).
    """
    if getattr(handle, "is_cram", False):
        return {"multiple_iterators": False}
    return {}


def require_csi_for_large_contigs(bam_path: str, handle: pysam.AlignmentFile) -> None:
    """Raise if ``handle`` has a contig > 2^29 bp but only a ``.bai`` index.

    A ``.bai`` cannot address coordinates beyond 2^29 bp, so pysam fetch over
    such a contig silently returns nothing — a catastrophic, silent evidence
    loss on the largest plant chromosomes. When a
    large contig is present a ``.csi`` index is required.

    A CRAM is exempt: the ``.crai`` index addresses large coordinates natively
    (it is not subject to the ``.bai`` 2^29 bin limit), so a CRAM with a
    ``.crai`` over a >512 Mb chromosome is fine.
    """
    if _is_cram(bam_path):
        return
    try:
        refs = list(handle.references or ())
        lengths = list(handle.lengths or ())
    except (ValueError, OSError):
        return  # no @SQ dictionary (check_sq=False) — nothing to enforce here
    big = [r for r, ln in zip(refs, lengths) if ln > CSI_REQUIRED_THRESHOLD]
    if big and not _has_csi_index(bam_path):
        raise FormatError(
            f"BAM has contig(s) > 2^29 bp ({', '.join(big[:3])}"
            f"{'…' if len(big) > 3 else ''}) but no .csi index — a .bai cannot "
            "address these coordinates (pysam fetch would silently return "
            "nothing). Re-index with `samtools index -c`.",
            path=bam_path,
        )


def bam_mapping_stats(
    bam_path: str | os.PathLike[str],
    *,
    reference_filename: str | os.PathLike[str] | None = None,
) -> dict[str, float | int]:
    """Overall mapped/unmapped read counts + mapping rate for a BAM/CRAM.

    Reads the index statistics (``idxstats``: per-reference mapped + the
    unmapped-read tally) so the genome-level report can answer "how well did the
    RNA-seq map" without a full pass over the alignments. A CRAM input is opened
    with the genome FASTA as
    ``reference_filename``.

    Returns ``{mapped, unmapped, total, mapping_rate}`` where ``mapping_rate`` =
    ``mapped / total`` (0.0 when empty — no div-by-zero). Requires an
    index; raises ``FileNotFoundError`` for a missing file and ``ValueError`` for
    a missing index, mirroring :class:`JunctionExtractor`.
    """
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    with _open_alignment(bam_path, reference_filename=reference_filename) as af:
        if not af.has_index():
            raise ValueError(f"BAM index (.bai/.csi) missing for {bam_path}")
        mapped = 0
        unmapped = 0
        for st in af.get_index_statistics():
            mapped += int(st.mapped)
            unmapped += int(st.unmapped)
        # Index also tracks reads with no reference (unmapped, unplaced).
        unmapped += int(af.nocoordinate)
    total = mapped + unmapped
    rate = (mapped / total) if total else 0.0
    return {
        "mapped": mapped,
        "unmapped": unmapped,
        "total": total,
        "mapping_rate": rate,
    }


def overall_mapping_rate(
    bam_paths: Sequence[str | os.PathLike[str]],
    *,
    reference_filename: str | os.PathLike[str] | None = None,
) -> dict[str, float | int] | None:
    """Aggregate :func:`bam_mapping_stats` across ``bam_paths`` (None if empty/all fail).

    Sums mapped/unmapped over every BAM that opens cleanly; a BAM that is missing
    or unindexed is logged and skipped so one bad sample never sinks the summary.
    Returns ``None`` when no BAM yields stats. ``reference_filename`` (the genome
    FASTA) is threaded so CRAM inputs decode offline; it is ignored
    for BAM.
    """
    if not bam_paths:
        return None
    mapped = unmapped = 0
    n_ok = 0
    for p in bam_paths:
        try:
            s = bam_mapping_stats(p, reference_filename=reference_filename)
        except (FileNotFoundError, ValueError, OSError) as exc:
            _log.warning("mapping-rate: skipping %s (%s)", p, exc)
            continue
        mapped += int(s["mapped"])
        unmapped += int(s["unmapped"])
        n_ok += 1
    if n_ok == 0:
        return None
    total = mapped + unmapped
    rate = (mapped / total) if total else 0.0
    return {
        "mapped": mapped,
        "unmapped": unmapped,
        "total": total,
        "mapping_rate": rate,
        "num_bams": n_ok,
    }


class JunctionExtractor:
    """Extract splice junctions from a coordinate-sorted, indexed BAM **or CRAM**.

    A CRAM input needs its genome FASTA passed as
    ``reference_filename`` to decode the reference-compressed records; detection
    is automatic (extension / magic) and junction extraction is otherwise
    identical. The index check covers ``.bai``/``.csi``/``.crai``.
    """

    def __init__(
        self,
        bam_path: str | os.PathLike[str],
        *,
        reference_filename: str | os.PathLike[str] | None = None,
    ) -> None:
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM not found: {bam_path}")
        self.bam_path = str(bam_path)
        self._bam: pysam.AlignmentFile | None = _open_alignment(
            self.bam_path, reference_filename=reference_filename
        )
        if not self._bam.has_index():
            self._bam.close()
            self._bam = None
            raise ValueError(
                f"alignment index (.bai/.csi/.crai) missing for {bam_path}"
            )
        try:
            require_csi_for_large_contigs(self.bam_path, self._bam)
        except FormatError:
            self._bam.close()
            self._bam = None
            raise

    # --- context manager ---
    def __enter__(self) -> JunctionExtractor:
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.close()

    def close(self) -> None:
        if self._bam is not None:
            self._bam.close()
            self._bam = None

    @staticmethod
    def _read_strand(read: pysam.AlignedSegment) -> str | None:
        """Strand from the XS tag; None if absent/unknown."""
        try:
            xs: Any = read.get_tag("XS")
        except KeyError:
            return None
        return xs if xs in ("+", "-") else None

    def extract_junctions(
        self,
        seqid: str,
        start: int,
        end: int,
        min_mapq: int = 10,
        min_overhang: int = 8,
        *,
        genome: Any | None = None,
        strand_source: str = DEFAULT_STRAND_SOURCE,
    ) -> list[SpliceJunction]:
        """Return junctions in ``[start, end)`` (0-based) aggregated by identity.

        Walks each read's CIGAR; every ``N`` op is an intron with
        ``donor = ref position`` and ``acceptor = donor + N``. Both flanking
        exon blocks must be ``>= min_overhang`` aligned bases. Unmapped,
        secondary, supplementary, and duplicate reads are skipped, and reads
        below ``min_mapq`` are dropped.

        Strand sourcing follows ``strand_source``:

        * ``"xs"`` — strand from the XS tag only; a read without XS is dropped.
        * ``"motif"`` — strand inferred from the canonical splice motif read off
          ``genome`` (a :class:`~helixforge.io.fasta.GenomeAccessor`); dropped
          when the motif is non-canonical or no genome is given.
        * ``"xs_then_motif"`` (default) — XS when present, otherwise the motif.
          **With no genome supplied this is identical to legacy ``"xs"``**, so
          existing call sites are byte-for-byte unchanged.

        When ``genome`` is supplied every emitted junction also carries its
        ``canonical`` motif class.
        """
        if strand_source not in STRAND_SOURCES:
            raise ValueError(
                f"strand_source must be one of {sorted(STRAND_SOURCES)}, "
                f"got {strand_source!r}"
            )
        if strand_source == "motif" and genome is None:
            _log.warning(
                "extract_junctions(strand_source='motif') with no genome: every "
                "junction will be dropped (no way to orient). Supply a genome or "
                "use 'xs'/'xs_then_motif'."
            )

        # Cache donor/acceptor dinucleotides per (donor, acceptor) so the genome
        # is read at most once per distinct junction, not once per read.
        dinuc_cache: dict[tuple[int, int], tuple[str, str] | None] = {}

        def _dinuc(donor: int, acceptor: int) -> tuple[str, str] | None:
            if genome is None:
                return None
            k = (donor, acceptor)
            if k not in dinuc_cache:
                try:
                    d = genome.get_sequence(seqid, donor, donor + 2, "+")
                    a = genome.get_sequence(seqid, acceptor - 2, acceptor, "+")
                    dinuc_cache[k] = (d, a)
                except Exception:  # noqa: BLE001 - never fail extraction on a fetch
                    dinuc_cache[k] = None
            return dinuc_cache[k]

        # key (donor, acceptor, strand) -> read_count
        counts: dict[tuple[int, int, str], int] = {}
        for read in self._bam.fetch(seqid, start, end):  # type: ignore[union-attr]
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
            ):
                continue
            if read.mapping_quality < min_mapq:
                continue
            cigar = read.cigartuples
            if not cigar or not any(op == _CIGAR_N for op, _ in cigar):
                continue
            xs_strand = self._read_strand(read)

            ref_pos = read.reference_start
            for i, (op, length) in enumerate(cigar):
                if op == _CIGAR_N:
                    donor = ref_pos
                    acceptor = ref_pos + length
                    left = self._flank_overhang(cigar, i, direction=-1)
                    right = self._flank_overhang(cigar, i, direction=+1)
                    if left >= min_overhang and right >= min_overhang:
                        strand = self._resolve_strand(
                            xs_strand, donor, acceptor, strand_source, _dinuc
                        )
                        if strand is not None:
                            key = (donor, acceptor, strand)
                            counts[key] = counts.get(key, 0) + 1
                    ref_pos += length
                elif op in _REF_CONSUMING:
                    ref_pos += length

        junctions: list[SpliceJunction] = []
        for (donor, acceptor, strand), count in counts.items():
            canonical: str | None = None
            dn = _dinuc(donor, acceptor)
            if dn is not None:
                canonical = classify_splice_motif(dn[0], dn[1], strand)
            junctions.append(
                SpliceJunction(
                    seqid,
                    donor,
                    acceptor,
                    strand,
                    read_count=count,
                    canonical=canonical,
                )
            )
        junctions.sort(key=lambda j: (j.donor, j.acceptor))
        return junctions

    @staticmethod
    def _resolve_strand(
        xs_strand: str | None,
        donor: int,
        acceptor: int,
        strand_source: str,
        dinuc_fn: Any,
    ) -> str | None:
        """Resolve a junction's strand per the ``strand_source`` policy.

        Returns ``'+'``/``'-'`` or ``None`` (drop the junction). See
        :meth:`extract_junctions` for the policy semantics.
        """
        if strand_source == "xs":
            return xs_strand
        # "motif" / "xs_then_motif"
        if strand_source == "xs_then_motif" and xs_strand is not None:
            return xs_strand
        dn = dinuc_fn(donor, acceptor)
        if dn is None:
            # No genome (or fetch failed): motif inference impossible. Fall back
            # to XS for the hybrid policy; drop for the strict motif policy.
            return xs_strand if strand_source == "xs_then_motif" else None
        return infer_strand_from_motif(dn[0], dn[1])

    @staticmethod
    def _flank_overhang(
        cigar: list[tuple[int, int]], n_index: int, direction: int
    ) -> int:
        """Sum aligned (M/=/X) bases in the exon block adjacent to cigar[n_index].

        ``direction`` is -1 (left flank) or +1 (right flank). The block stops at
        the next ``N`` op or the read end.
        """
        total = 0
        i = n_index + direction
        while 0 <= i < len(cigar):
            op, length = cigar[i]
            if op == _CIGAR_N:
                break
            if op in _MATCH_OPS:
                total += length
            i += direction
        return total

    def extract_junctions_multisample(
        self,
        bam_paths: list[str],
        seqid: str,
        start: int,
        end: int,
        min_mapq: int = 10,
        min_overhang: int = 8,
        *,
        genome: Any | None = None,
        strand_source: str = DEFAULT_STRAND_SOURCE,
    ) -> list[SpliceJunction]:
        """Merge junctions across BAMs: read_count summed, samples counted.

        ``genome`` / ``strand_source`` are passed through to each per-BAM
        :meth:`extract_junctions`; the merged junction keeps the
        canonical motif class (identical across samples for the same intron).
        """
        # key (donor, acceptor, strand) -> [read_count, sample_count, canonical]
        merged: dict[tuple[int, int, str], list[Any]] = {}
        for path in bam_paths:
            with JunctionExtractor(path) as je:
                for j in je.extract_junctions(
                    seqid,
                    start,
                    end,
                    min_mapq=min_mapq,
                    min_overhang=min_overhang,
                    genome=genome,
                    strand_source=strand_source,
                ):
                    key = (j.donor, j.acceptor, j.strand)
                    if key not in merged:
                        merged[key] = [0, 0, None]
                    merged[key][0] += j.read_count
                    merged[key][1] += 1
                    if merged[key][2] is None:
                        merged[key][2] = j.canonical
        out = [
            SpliceJunction(
                seqid,
                donor,
                acceptor,
                strand,
                read_count=rc,
                samples=ns,
                canonical=canon,
            )
            for (donor, acceptor, strand), (rc, ns, canon) in merged.items()
        ]
        out.sort(key=lambda j: (j.donor, j.acceptor))
        return out


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Coalesce overlapping / adjacent ``[start, end)`` intervals into maximal spans.

    Used to turn a per-exon region list into the minimal set of contiguous spans
    to scan, so one htslib pass covers a whole gene locus instead of one indexed
    seek per exon.
    """
    out: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1] = (out[-1][0], e)
        else:
            out.append((s, e))
    return out


class CoverageCalculator:
    """Per-base / mean coverage from a BAM or bigWig source."""

    def __init__(self, handle: Any, source_type: str) -> None:
        self._handle: Any = handle
        self.source_type = source_type

    @classmethod
    def from_bam(
        cls,
        bam_path: str | os.PathLike[str],
        *,
        reference_filename: str | os.PathLike[str] | None = None,
    ) -> CoverageCalculator:
        """Open a BAM **or CRAM** coverage source.

        A CRAM (``.cram`` / ``CRAM`` magic) is opened with the genome FASTA as
        ``reference_filename`` so coverage decodes; a BAM ignores it. The
        large-chromosome CSI gate applies to BAM only (a ``.crai`` addresses
        >2^29 bp natively).
        """
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM not found: {bam_path}")
        handle = _open_alignment(bam_path, reference_filename=reference_filename)
        try:
            require_csi_for_large_contigs(str(bam_path), handle)
        except FormatError:
            handle.close()
            raise
        return cls(handle, "bam")

    @classmethod
    def from_bigwig(cls, bw_path: str | os.PathLike[str]) -> CoverageCalculator:
        try:
            import pyBigWig
        except ImportError as exc:  # pragma: no cover - exercised via monkeypatch
            raise ImportError(
                "bigWig support requires pyBigWig: pip install helixforge[bigwig]"
            ) from exc
        if not os.path.exists(bw_path):
            raise FileNotFoundError(f"bigWig not found: {bw_path}")
        return cls(pyBigWig.open(bw_path), "bigwig")

    # --- context manager ---
    def __enter__(self) -> CoverageCalculator:
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.close()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def region_coverage_array(
        self,
        seqid: str,
        start: int,
        end: int,
        *,
        min_mapq: int = 0,
        max_depth: int = COVERAGE_MAX_DEPTH,
    ) -> list[float]:
        """Return a length-``(end-start)`` list of per-base depth.

        Read-filtering policy: the pileup uses
        ``stepper="all"``, which already excludes unmapped / **secondary** /
        QC-fail / **duplicate** reads, and ``max_depth`` is set **explicitly**
        (default :data:`~helixforge.constants.COVERAGE_MAX_DEPTH`) so a deep
        pileup is no longer silently capped at htslib's ~8000 default. Reads
        below ``min_mapq`` are excluded (default 0 = count every primary read —
        the legacy behaviour; pass ``min_mapq=UNIQUE_MIN_MAPQ`` for unique-only
        coverage, consistent with junction extraction).
        """
        if end <= start:
            raise ValueError(f"require start < end, got start={start} end={end}")
        length = end - start
        if self.source_type == "bam":
            depth: list[float] = [0] * length
            for col in self._handle.pileup(
                seqid,
                start,
                end,
                truncate=True,
                stepper="all",
                max_depth=max_depth,
                **_pileup_kwargs(self._handle),
            ):
                pos = col.reference_pos
                if start <= pos < end:
                    # Count only truly aligned bases: exclude intron skips (N)
                    # and deletions (D), which htslib otherwise reports; and
                    # reads below the MAPQ floor.
                    depth[pos - start] = sum(
                        1
                        for p in col.pileups
                        if not p.is_refskip
                        and not p.is_del
                        and p.alignment.mapping_quality >= min_mapq
                    )
            return depth
        # bigwig — precomputed signal, no per-read MAPQ/depth notion.
        import math

        vals = self._handle.values(seqid, start, end)
        return [0.0 if (v is None or math.isnan(v)) else float(v) for v in vals]

    def region_coverage_arrays(
        self,
        regions: Sequence[tuple[str, int, int]],
        *,
        min_mapq: int = 0,
        max_depth: int = COVERAGE_MAX_DEPTH,
    ) -> dict[tuple[str, int, int], "np.ndarray"]:
        """Per-base depth for many regions in **one pass per merged span**.

        Returns ``{(seqid, start, end): depth_array}`` whose values are **exactly**
        what :meth:`region_coverage_array` returns for each region individually,
        but issues one pileup per merged span (overlapping / adjacent regions
        coalesced by :func:`_merge_intervals`) instead of one pileup per region.
        Because a base's pileup depth does not depend on the query window, slicing
        a span's depth to a contained region is identical to querying that region
        directly — so this is byte-for-byte the old path, just with the per-exon
        indexed seeks collapsed to one per gene locus.

        The reason ``--threads`` did nothing on the old per-exon path was twofold:
        (1) O(num_regions) indexed seeks per BAM dominated wall time, and (2) the
        ``pileup`` column loop is pure-Python and GIL-bound, so a *thread* pool
        over BAMs could not use more than one core. (1) is fixed here (per-locus
        seeks); (2) is fixed by the caller scanning BAMs in a **process** pool.
        """
        import numpy as np

        for seqid, s, e in regions:
            if e <= s:
                raise ValueError(f"require start < end, got start={s} end={e}")

        by_seqid: dict[str, list[tuple[int, int]]] = {}
        for seqid, s, e in regions:
            by_seqid.setdefault(seqid, []).append((s, e))

        out: dict[tuple[str, int, int], "np.ndarray"] = {}
        for seqid, ivs in by_seqid.items():
            uniq = sorted(set(ivs))  # dedup shared exons; sorted for the walk below
            spans = _merge_intervals(uniq)
            # One pileup per merged span (reusing the exact per-region counting),
            # then a single sorted walk assigns each region to its containing span
            # (both lists ascending) — O(n), not O(spans x regions).
            depths = [
                np.asarray(
                    self.region_coverage_array(
                        seqid, lo, hi, min_mapq=min_mapq, max_depth=max_depth
                    ),
                    dtype=np.float64,
                )
                for lo, hi in spans
            ]
            si = 0
            for s, e in uniq:
                while si < len(spans) and not (spans[si][0] <= s and e <= spans[si][1]):
                    si += 1
                lo, _hi = spans[si]
                out[(seqid, s, e)] = depths[si][s - lo : e - lo]
        return out

    def region_coverage_both(
        self,
        seqid: str,
        start: int,
        end: int,
        *,
        unique_min_mapq: int = UNIQUE_MIN_MAPQ,
        max_depth: int = COVERAGE_MAX_DEPTH,
    ) -> tuple[list[float], list[float]]:
        """Return ``(all_depth, unique_depth)`` per-base arrays in one pileup pass.

        ``all_depth`` counts every primary read (``min_mapq=0``); ``unique_depth``
        counts only reads with ``mapping_quality >= unique_min_mapq``. Reporting
        both lets a caller judge paralog/polyploid (multi-mapper-inflated)
        regions — where the two diverge — cheaply.
        For a bigWig source both arrays are identical (no per-read MAPQ).
        """
        if end <= start:
            raise ValueError(f"require start < end, got start={start} end={end}")
        if self.source_type != "bam":
            arr = self.region_coverage_array(seqid, start, end)
            return arr, list(arr)
        length = end - start
        all_d: list[float] = [0] * length
        uniq_d: list[float] = [0] * length
        for col in self._handle.pileup(
            seqid,
            start,
            end,
            truncate=True,
            stepper="all",
            max_depth=max_depth,
            **_pileup_kwargs(self._handle),
        ):
            pos = col.reference_pos
            if start <= pos < end:
                a = 0
                u = 0
                for p in col.pileups:
                    if p.is_refskip or p.is_del:
                        continue
                    a += 1
                    if p.alignment.mapping_quality >= unique_min_mapq:
                        u += 1
                all_d[pos - start] = a
                uniq_d[pos - start] = u
        return all_d, uniq_d

    def mean_coverage(
        self,
        seqid: str,
        start: int,
        end: int,
        *,
        genome: Any | None = None,
        n_mask: list[bool] | None = None,
        min_mapq: int = 0,
        max_depth: int = COVERAGE_MAX_DEPTH,
    ) -> float:
        """Mean depth over the region; 0.0 if no reads/values.

        Assembly gaps (reference ``N`` runs) carry depth 0 by construction, so a
        gene spanning a draft-scaffold gap is biased toward LOW/SILENT
        When the reference ``N``-mask is supplied — either a
        ``genome`` accessor exposing ``get_sequence(seqid, start, end, '+')`` or
        an explicit ``n_mask`` (one bool per base, ``True`` = reference ``N``) —
        those positions are excluded from the **denominator**, so the mean is
        taken over non-``N`` (mappable) bases only. With no mask the behaviour is
        byte-identical to before (mean over the whole region).

        ``min_mapq`` / ``max_depth`` are passed through to
        :meth:`region_coverage_array`; the defaults (``min_mapq=0``)
        keep the legacy "count every primary read" semantics.
        """
        arr = self.region_coverage_array(
            seqid, start, end, min_mapq=min_mapq, max_depth=max_depth
        )
        if not arr:
            return 0.0
        mask = self._resolve_n_mask(seqid, start, end, genome, n_mask)
        if mask is not None:
            kept = [d for d, is_n in zip(arr, mask) if not is_n]
            if not kept:
                return 0.0  # whole region is an assembly gap
            return float(sum(kept)) / len(kept)
        return float(sum(arr)) / len(arr)

    def mean_coverage_both(
        self,
        seqid: str,
        start: int,
        end: int,
        *,
        genome: Any | None = None,
        n_mask: list[bool] | None = None,
        unique_min_mapq: int = UNIQUE_MIN_MAPQ,
        max_depth: int = COVERAGE_MAX_DEPTH,
    ) -> tuple[float, float]:
        """Return ``(all_mean, unique_mean)`` depth, N-masked if a mask is given.

        Mirrors :meth:`mean_coverage` but reports the all-reads and unique-only
        means together (one pileup pass via :meth:`region_coverage_both`) so a
        caller can spot multi-mapper-inflated paralog/polyploid regions
        The same reference-``N`` denominator masking applies to
        both.
        """
        all_arr, uniq_arr = self.region_coverage_both(
            seqid, start, end, unique_min_mapq=unique_min_mapq, max_depth=max_depth
        )
        if not all_arr:
            return 0.0, 0.0
        mask = self._resolve_n_mask(seqid, start, end, genome, n_mask)

        def _mean(arr: list[float]) -> float:
            if mask is not None:
                kept = [d for d, is_n in zip(arr, mask) if not is_n]
                if not kept:
                    return 0.0
                return float(sum(kept)) / len(kept)
            return float(sum(arr)) / len(arr)

        return _mean(all_arr), _mean(uniq_arr)

    @staticmethod
    def _resolve_n_mask(
        seqid: str,
        start: int,
        end: int,
        genome: Any | None,
        n_mask: list[bool] | None,
    ) -> list[bool] | None:
        """Return a per-base ``True``=reference-``N`` mask, or ``None`` if no source.

        An explicit ``n_mask`` wins; otherwise the reference sequence is fetched
        from ``genome`` and each ``N`` (case-insensitive) marked. A fetch failure
        falls back to ``None`` (no masking) rather than crashing the coverage pass.
        """
        if n_mask is not None:
            return n_mask
        if genome is None:
            return None
        try:
            ref = genome.get_sequence(seqid, start, end, "+")
        except Exception:  # noqa: BLE001 - never fail coverage on a mask fetch
            return None
        return [ch in ("N", "n") for ch in ref]

    @classmethod
    def mean_coverage_multisample(
        cls,
        sources: list[str],
        seqid: str,
        start: int,
        end: int,
        source_type: str = "bam",
        *,
        genome: Any | None = None,
    ) -> float:
        """Mean of per-source mean coverage across ``sources`` (0.0 if empty).

        ``genome`` (optional) supplies the reference ``N``-mask passed through to
        each per-source :meth:`mean_coverage`.
        """
        opener = cls.from_bam if source_type == "bam" else cls.from_bigwig
        means: list[float] = []
        for src in sources:
            with opener(src) as cov:
                means.append(cov.mean_coverage(seqid, start, end, genome=genome))
        if not means:
            return 0.0
        return float(sum(means)) / len(means)


class CoveragePool:
    """Open each coverage source **once** for a whole classification/scoring pass.

    The per-locus ``from_bam``/``from_bigwig`` open in
    ``classify._get_exonic_coverage`` would re-open every source for every locus
    (100 K loci × 10 BAMs ⇒ 1 M opens). This context manager opens one
    :class:`CoverageCalculator` per source on ``__enter__`` and closes them all
    deterministically on ``__exit__``; callers query the held handles per locus.
    Mean values are unchanged — only the open/close lifecycle moves.

    ``calculator_cls`` is injectable so a caller (and the unit tests) can supply a
    fake; it defaults to :class:`CoverageCalculator`. Handles are managed via the
    context-manager protocol (``ExitStack``), so a fake need only implement
    ``from_bam``/``from_bigwig`` + ``__enter__``/``__exit__`` (no ``close``).
    """

    def __init__(
        self,
        bam_paths: list[str] | None = None,
        bigwig_paths: list[str] | None = None,
        calculator_cls: type[CoverageCalculator] | None = None,
    ) -> None:
        self._calculator_cls: type[CoverageCalculator] = (
            calculator_cls or CoverageCalculator
        )
        self.bam_paths: list[str] = list(bam_paths or [])
        self.bigwig_paths: list[str] = list(bigwig_paths or [])
        self._stack: ExitStack | None = None
        self._bam: list[CoverageCalculator] = []
        self._bigwig: list[CoverageCalculator] = []

    def __enter__(self) -> CoveragePool:
        self._stack = ExitStack()
        cls = self._calculator_cls
        self._bam = [self._stack.enter_context(cls.from_bam(p)) for p in self.bam_paths]
        self._bigwig = [
            self._stack.enter_context(cls.from_bigwig(p)) for p in self.bigwig_paths
        ]
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        if self._stack is not None:
            self._stack.close()
            self._stack = None
        self._bam = []
        self._bigwig = []

    def calculators(self, source_type: str) -> list[CoverageCalculator]:
        """Open :class:`CoverageCalculator` handles for ``'bam'`` or ``'bigwig'``."""
        return self._bigwig if source_type == "bigwig" else self._bam


def parse_star_sj_tab(
    sj_path: str | os.PathLike[str],
    min_unique_reads: int = 3,
    min_overhang: int = 8,
    *,
    genome: Any | None = None,
) -> list[SpliceJunction]:
    """Parse a STAR ``SJ.out.tab`` → ``list[SpliceJunction]``.

    Columns: chrom, intron_start (1-based), intron_end (1-based inclusive),
    strand (0=undef, 1='+', 2='-'), motif, annotated, unique reads, multi-map
    reads, max spliced overhang. Converts ``donor = start - 1``,
    ``acceptor = end``. Strand 0 and junctions below the read/overhang
    thresholds are excluded.

    The **motif column (col 4)** is no longer
    discarded — it is mapped to the canonical class and carried on each
    junction's ``canonical`` field, and the **multi-map column (col 7)** is
    carried on ``multimap_reads``. When ``genome`` (a
    :class:`~helixforge.io.fasta.GenomeAccessor`) is supplied the genome-derived
    motif is computed and **cross-checked** against STAR's; a disagreement is
    logged and the genome call (ground truth) wins.
    """
    if not os.path.exists(sj_path):
        raise FileNotFoundError(f"STAR SJ tab not found: {sj_path}")

    strand_map: dict[str, str] = {"1": "+", "2": "-"}
    out: list[SpliceJunction] = []
    with open(sj_path) as fh:
        for idx, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                raise FormatError(
                    f"STAR SJ line has {len(cols)} columns, need 9 (truncated?)",
                    path=sj_path,
                    line=idx,
                    value=line.rstrip("\n")[:60],
                )
            chrom = cols[0]
            # Contextual numeric parsing: a header line or a
            # comma-decimal now raises a located FormatError, not a bare
            # uncontextualised ValueError.
            intron_start = _star_int(cols[1], path=sj_path, line=idx, col=2)
            intron_end = _star_int(cols[2], path=sj_path, line=idx, col=3)
            strand_code = cols[3]
            motif_code = cols[4]
            unique_reads = _star_int(cols[6], path=sj_path, line=idx, col=7)
            multimap_reads = _star_int(cols[7], path=sj_path, line=idx, col=8)
            overhang = _star_int(cols[8], path=sj_path, line=idx, col=9)

            if strand_code not in _STAR_STRAND_CODES_LOCAL:
                raise FormatError(
                    "invalid STAR strand code (expected 0, 1, or 2)",
                    path=sj_path,
                    line=idx,
                    value=strand_code,
                )
            strand = strand_map.get(strand_code)
            if strand is None:
                continue  # undefined strand (code 0) — preserved skip
            if unique_reads < min_unique_reads or overhang < min_overhang:
                continue

            donor = intron_start - 1
            acceptor = intron_end
            # STAR-derived canonical class (col 4); unknown codes -> None.
            canonical = _STAR_MOTIF_CANONICAL.get(motif_code)
            if genome is not None:
                geno = _genome_motif(genome, chrom, donor, acceptor, strand)
                if geno is not None:
                    if canonical is not None and geno != canonical:
                        _log.warning(
                            "STAR/genome splice-motif disagreement at %s:%d-%d "
                            "(%s): STAR=%s genome=%s — using genome",
                            chrom,
                            donor,
                            acceptor,
                            strand,
                            canonical,
                            geno,
                        )
                    canonical = geno  # genome is ground truth
            out.append(
                SpliceJunction(
                    seqid=chrom,
                    donor=donor,
                    acceptor=acceptor,
                    strand=strand,
                    read_count=unique_reads,
                    canonical=canonical,
                    multimap_reads=multimap_reads,
                )
            )
    out.sort(key=lambda j: (j.seqid, j.donor, j.acceptor))
    return out


def _genome_motif(
    genome: Any, seqid: str, donor: int, acceptor: int, strand: str
) -> str | None:
    """Genome-derived canonical class for ``[donor, acceptor)``, or None on fetch fail."""
    try:
        d = genome.get_sequence(seqid, donor, donor + 2, "+")
        a = genome.get_sequence(seqid, acceptor - 2, acceptor, "+")
    except Exception:  # noqa: BLE001 - never fail parsing on a mask fetch
        return None
    return classify_splice_motif(d, a, strand)
