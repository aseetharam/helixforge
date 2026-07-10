"""Standalone RNA-seq + protein scoring of any GFF3."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import attrs

from helixforge.io.bam import (
    CoverageCalculator,
    JunctionExtractor,
    parse_star_sj_tab,
)
from helixforge.io.gff import GFF3Parser
from helixforge.io.stringtie import (
    TpmRecord as _StRecord,
    best_overlapping_tpm as _best_overlapping_tpm,
    build_tpm_overlap_index,
)
from helixforge.stats.evidence_concordance import intron_concordance
from helixforge.utils.regions import gff3_to_internal, parse_region

if TYPE_CHECKING:
    import pandas as pd
    from helixforge.io.bam import SpliceJunction  # type: ignore[attr-defined]  # re-exported from reconcile.models; no __all__ in io.bam
    from helixforge.reconcile.models import MiniprotAlignment

# Default minimum read support for a junction to qualify as evidence (mirrors
# stats.evidence_concordance.intron_concordance and io.bam.parse_star_sj_tab).
DEFAULT_MIN_READS = 3
DEFAULT_MIN_MAPQ = 10
DEFAULT_MIN_OVERHANG = 8
# Coverage to count an exon as "expressed" / a boundary as supported (v1 default).
DEFAULT_MIN_EXON_COVERAGE = 5

# RNA-AED weights: v1 ``EvidenceScorerConfig`` defaults, ported verbatim. They
# sum to 1.0 (junction is the primary evidence for multi-exon genes).
RNA_W_JUNCTION = 0.5
RNA_W_COVERAGE = 0.3
RNA_W_BOUNDARY = 0.2

# protein-AED weights: mirror the RNA weights (structure primary, like junctions).
PROT_W_STRUCT = 0.5
PROT_W_CDS = 0.3
PROT_W_PROT = 0.2

# Stable column order for the per-transcript table (so an empty DataFrame still
# carries the schema and ``summarize_evidence`` works on it). The two AED blocks
# are appended after the historical RNA columns; each is ``NA`` when its evidence
# was not supplied.
TSV_COLUMNS = [
    "gene_id",
    "transcript_id",
    "seqid",
    "strand",
    "start",
    "end",
    "num_exons",
    "num_introns",
    "supported",
    "contradicted",
    "novel_in_data",
    "junction_support_fraction",
    "intron_precision",
    "intron_recall",
    "intron_f1",
    "mean_coverage",
    "tpm",
    # RNA-AED block (present only when BAM coverage was supplied)
    "rna_aed",
    "rna_junction_ratio",
    "rna_coverage_ratio",
    "rna_boundary_ratio",
    # protein-AED block (present only when a proteome / miniprot GFF was supplied)
    "protein_id",
    "protein_aed",
    "protein_struct_ratio",
    "protein_cds_cov_ratio",
    "protein_prot_cov_ratio",
]


# ---------------------------------------------------------------------------
# Lightweight shims so the reused gene-oriented stats helpers can score a single
# transcript. They expose only the attributes those helpers read; deliberately
# no ``introns`` attribute, so ``_intron_bounds`` derives introns from exons.
# ---------------------------------------------------------------------------


class _TxShim:
    def __init__(
        self,
        transcript_id: str | None,
        seqid: str,
        strand: str,
        exons: list[Any],
    ) -> None:
        self.transcript_id = transcript_id
        self.seqid = seqid
        self.strand = strand
        self.exons = exons


class _GeneShim:
    def __init__(
        self,
        gene_id: str,
        seqid: str,
        strand: str,
        start: int,
        end: int,
        transcripts: list[_TxShim],
    ) -> None:
        self.gene_id = gene_id
        self.seqid = seqid
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts = transcripts


def _adapt_transcript(
    transcript_like: Any,
) -> tuple[str | None, str, str, list[Any], int, int]:
    """Normalise to ``(transcript_id, seqid, strand, exons, start, end)``.

    Accepts a ``parse_genes_generic``-style transcript dict (enriched with the
    parent gene's ``seqid``/``strand``) or any object with the matching
    attributes. The span is taken verbatim when present, else the exon envelope.
    """
    if isinstance(transcript_like, dict):
        tid: str | None = transcript_like.get("transcript_id")
        seqid: str = transcript_like["seqid"]
        strand: str = transcript_like["strand"]
        exons: list[Any] = transcript_like["exons"]
        start: int | None = transcript_like.get("start")
        end: int | None = transcript_like.get("end")
    else:
        tid = getattr(transcript_like, "transcript_id", None)
        seqid = transcript_like.seqid
        strand = transcript_like.strand
        exons = transcript_like.exons
        start = getattr(transcript_like, "start", None)
        end = getattr(transcript_like, "end", None)

    if not exons:
        raise ValueError("transcript has no exons; cannot score evidence")
    if start is None or end is None:
        start = min(e.start for e in exons)
        end = max(e.end for e in exons)
    return tid, seqid, strand, exons, start, end


# The TPM overlap index + per-model lookup live in ``io.stringtie`` so the
# ``evidence`` scorer and the ``reconcile`` pipeline share one implementation
# (imported above as ``_StRecord`` / ``_best_overlapping_tpm``).


def _region_to_internal(
    region: str | None,
) -> tuple[str | None, int | None, int | None]:
    """Parse a region string to ``(seqid, lo, hi)`` internal coords (or ``None``).

    ``None`` → ``(None, None, None)`` (no filter). A bare ``seqid`` →
    ``(seqid, None, None)``. A ``seqid:start-end`` span is 1-based inclusive and
    is converted here, the single boundary, to 0-based half-open ``[lo, hi)``.
    """
    if region is None:
        return None, None, None
    seqid, start, end = parse_region(region)
    if start is None:
        return seqid, None, None
    assert end is not None  # parse_region returns (None,None) or (int,int) together
    lo, hi = gff3_to_internal(start, end)
    return seqid, lo, hi


# ---------------------------------------------------------------------------
# Junction collection (process-parallel across files)
# ---------------------------------------------------------------------------


def _bam_targets(
    bam_paths: list[str | Path],
    r_seqid: str | None,
    r_lo: int | None,
    r_hi: int | None,
    *,
    reference_filename: str | Path | None = None,
) -> list[tuple[str, int, int]]:
    """Per-contig ``(seqid, start, end)`` fetch windows for the BAM/CRAM scan.

    With no region, every reference in the (first) BAM/CRAM header is scanned end
    to end. With a region, only the matching contig (whole, or the given span).
    A CRAM source needs its genome FASTA as ``reference_filename`` to decode
    a BAM ignores it.
    """
    from helixforge.io.bam import _open_alignment

    with _open_alignment(bam_paths[0], reference_filename=reference_filename) as af:
        refs = list(zip(af.references, af.lengths))
    if r_seqid is not None:
        for name, length in refs:
            if name == r_seqid:
                lo = r_lo if r_lo is not None else 0
                hi = r_hi if r_hi is not None else length
                return [(name, lo, hi)]
        return []  # region seqid not present in the BAM
    return [(name, 0, length) for name, length in refs]


def _extract_bam_junctions(
    path: str | Path,
    targets: list[tuple[str, int, int]],
    min_mapq: int,
    min_overhang: int,
    reference_filename: str | Path | None,
) -> list[SpliceJunction]:
    """All splice junctions from one BAM/CRAM over ``targets`` (one file = a worker)."""
    out: list[SpliceJunction] = []
    with JunctionExtractor(path, reference_filename=reference_filename) as je:
        for seqid, lo, hi in targets:
            out.extend(
                je.extract_junctions(
                    seqid, lo, hi, min_mapq=min_mapq, min_overhang=min_overhang
                )
            )
    return out


def _extract_star_junctions(
    path: str | Path,
    min_reads: int,
    min_overhang: int,
    r_seqid: str | None,
    r_lo: int | None,
    r_hi: int | None,
) -> list[SpliceJunction]:
    """STAR SJ.out.tab junctions from one file, region-filtered (one file = a worker)."""
    out: list[SpliceJunction] = []
    for j in parse_star_sj_tab(
        path, min_unique_reads=min_reads, min_overhang=min_overhang
    ):
        if r_seqid is not None and j.seqid != r_seqid:
            continue
        if r_lo is not None and not (j.acceptor > r_lo and j.donor < r_hi):  # type: ignore[operator]  # r_hi is int when r_lo is not None
            continue
        out.append(j)
    return out


def collect_junctions(
    bam_paths: list[str | Path] | None = None,
    star_sj_paths: list[str | Path] | None = None,
    region: str | None = None,
    min_reads: int = DEFAULT_MIN_READS,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    min_overhang: int = DEFAULT_MIN_OVERHANG,
    reference_filename: str | Path | None = None,
    threads: int = 1,
) -> list[SpliceJunction]:
    """Union of splice junctions from BAM and/or STAR SJ files.

    Identical ``(seqid, donor, acceptor, strand)`` junctions are merged: read
    counts summed and ``samples`` set to the number of distinct source files that
    observed it. At least one of ``bam_paths`` / ``star_sj_paths`` is required.

    Per-file extraction is independent, so with ``threads > 1`` the files are
    scanned in a **process** pool (one file per worker); the per-file read scan +
    CIGAR walk is GIL-bound, so processes, not threads, engage multiple cores.
    The merge is then done deterministically in this process, so the result is
    identical to the serial path. ``min_mapq`` / ``min_overhang`` gate BAM
    extraction; ``min_reads`` is the per-file STAR unique-read floor. Returns a
    list sorted by ``(seqid, donor, acceptor)``.
    """
    if not bam_paths and not star_sj_paths:
        raise ValueError(
            "collect_junctions requires at least one of bam_paths / star_sj_paths"
        )
    r_seqid, r_lo, r_hi = _region_to_internal(region)
    bam_paths = list(bam_paths or [])
    star_sj_paths = list(star_sj_paths or [])

    targets = (
        _bam_targets(
            bam_paths, r_seqid, r_lo, r_hi, reference_filename=reference_filename
        )
        if bam_paths
        else []
    )

    # One task per source file; source_id keeps the per-file identity for the
    # ``samples`` count regardless of completion order. ``partial`` (not a lambda)
    # so each task is picklable for the process pool: the per-file read scan +
    # CIGAR walk is GIL-bound, so processes (not threads) are what engage cores.
    tasks: list[tuple[tuple[str, int], Any]] = []
    for i, path in enumerate(bam_paths):
        tasks.append(
            (
                ("bam", i),
                partial(
                    _extract_bam_junctions,
                    path,
                    targets,
                    min_mapq,
                    min_overhang,
                    reference_filename,
                ),
            )
        )
    for i, path in enumerate(star_sj_paths):
        tasks.append(
            (
                ("sj", i),
                partial(
                    _extract_star_junctions,
                    path,
                    min_reads,
                    min_overhang,
                    r_seqid,
                    r_lo,
                    r_hi,
                ),
            )
        )

    per_source: list[tuple[tuple[str, int], list[SpliceJunction]]] = []
    n_workers = max(1, min(threads, len(tasks)))
    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futs = {ex.submit(fn): sid for sid, fn in tasks}
            for fut in as_completed(futs):
                per_source.append((futs[fut], fut.result()))
    else:
        per_source = [(sid, fn()) for sid, fn in tasks]

    # Deterministic merge: order by source_id so accumulation is reproducible.
    merged: dict[tuple[str, int, int, str], dict[str, Any]] = {}
    for source_id, juncs in sorted(per_source, key=lambda t: t[0]):
        for j in juncs:
            key = (j.seqid, j.donor, j.acceptor, j.strand)
            rec = merged.get(key)
            if rec is None:
                merged[key] = {"rep": j, "reads": j.read_count, "sources": {source_id}}
            else:
                rec["reads"] += j.read_count
                rec["sources"].add(source_id)

    out = [
        attrs.evolve(rec["rep"], read_count=rec["reads"], samples=len(rec["sources"]))
        for rec in merged.values()
    ]
    out.sort(key=lambda j: (j.seqid, j.donor, j.acceptor))
    return out


# ---------------------------------------------------------------------------
# Coverage collection (precomputed once, never re-read per gene)
# ---------------------------------------------------------------------------


def _bam_exon_arrays(
    path: str | Path,
    regions: list[tuple[str, int, int]],
    min_mapq: int,
    reference_filename: str | Path | None,
) -> dict[tuple[str, int, int], Any]:
    """Per-base depth array for each ``(seqid, start, end)`` region in one BAM.

    One pass over the BAM: :meth:`CoverageCalculator.region_coverage_arrays` issues
    a single htslib ``count_coverage`` per merged span (gene locus), not a pileup
    per exon, so the per-BAM scan is C-level and seek-light. The returned depth
    arrays are byte-identical to the old per-exon pileup path (same Q13 base
    floor, same read filter).
    """
    cov = CoverageCalculator.from_bam(path, reference_filename=reference_filename)
    try:
        return cov.region_coverage_arrays(regions, min_mapq=min_mapq)
    finally:
        cov.close()


def collect_coverage(
    bam_paths: list[str | Path],
    exon_regions: Sequence[tuple[str, int, int]],
    *,
    min_mapq: int = 0,
    reference_filename: str | Path | None = None,
    threads: int = 1,
) -> dict[tuple[str, int, int], tuple[float, float, float]]:
    """Combine per-exon coverage across BAMs into ``region -> (mean, median, min)``.

    Each BAM is scanned independently (process-parallel; one htslib
    ``count_coverage`` pass per gene locus, not a pileup per exon) and the
    per-base depth arrays are combined by **element-wise mean across
    samples**, so a single BAM reproduces that BAM's coverage exactly and adding
    samples averages depth rather than inflating it. The combined per-exon array
    yields ``(mean, median, min)``; the per-base ``median`` feeds both the
    expressed-exon ratio and the boundary ratio (v1 semantics), ``mean`` the
    exon-weighted coverage. The cache is keyed by the exon region so the scoring
    phase never touches a BAM again.

    Memory is bounded by the total exonic length (one combined array per distinct
    exon), independent of BAM count; a transient per-BAM array set is held only
    while a worker is in flight.
    """
    import numpy as np

    regions = sorted(set(exon_regions))
    if not regions:
        return {}

    accum: dict[tuple[str, int, int], Any] = {}

    def _merge(per_bam: dict[tuple[str, int, int], Any]) -> None:
        for key, arr in per_bam.items():
            cur = accum.get(key)
            accum[key] = arr if cur is None else cur + arr

    n_workers = max(1, min(threads, len(bam_paths)))
    if n_workers > 1:
        # Process pool, not threads: count_coverage releases the GIL during the
        # htslib scan but the per-span A+C+G+T summation is GIL-bound, so threads
        # plateau at ~1 core; processes scale.
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futs = [
                ex.submit(_bam_exon_arrays, p, regions, min_mapq, reference_filename)
                for p in bam_paths
            ]
            for fut in as_completed(futs):
                _merge(fut.result())
    else:
        for p in bam_paths:
            _merge(_bam_exon_arrays(p, regions, min_mapq, reference_filename))

    n = len(bam_paths)
    cache: dict[tuple[str, int, int], tuple[float, float, float]] = {}
    for key, summed in accum.items():
        combined = summed / n if n else summed
        if combined.size == 0:
            cache[key] = (0.0, 0.0, 0.0)
        else:
            cache[key] = (
                float(combined.mean()),
                float(np.median(combined)),
                float(combined.min()),
            )
    return cache


def _transcript_coverage_metrics(
    seqid: str,
    strand: str,
    exons: list[Any],
    cache: dict[tuple[str, int, int], tuple[float, float, float]],
    min_exon_coverage: int,
) -> tuple[float, float, float]:
    """``(weighted_mean, coverage_ratio, boundary_ratio)`` for one transcript.

    ``coverage_ratio`` is the fraction of exons whose per-base **median** depth ≥
    ``min_exon_coverage`` (v1 ``EvidenceScorer.score_gene`` semantics).

    ``boundary_ratio`` is ``(start_supported + stop_supported) / 2`` where a
    boundary is supported when its **terminal exon is expressed**, the same
    per-base median ≥ threshold test, applied to the lowest- and highest-
    coordinate exons (plus strand: the low-coordinate exon holds the start, the
    high-coordinate exon the stop; minus strand swaps the labels, the ratio is
    symmetric, but the labelling is kept correct). This replaces the old per-exon
    *minimum* depth: a single low-coverage base anywhere in a terminal exon, a
    near-universal coverage dip, or a Helixer UTR annotated past where reads
    actually reach, drove the minimum to ~0 and pinned ``boundary_ratio`` to 0
    even for deeply covered genes. The median is robust to both.
    """
    ordered = sorted(exons, key=lambda e: e.start)
    stats = [cache.get((seqid, e.start, e.end), (0.0, 0.0, 0.0)) for e in ordered]
    total = sum(e.end - e.start for e in ordered)
    weighted = (
        sum(st[0] * (e.end - e.start) for st, e in zip(stats, ordered)) / total
        if total
        else 0.0
    )
    n = len(ordered)
    n_expressed = sum(1 for st in stats if st[1] >= min_exon_coverage)
    coverage_ratio = n_expressed / n if n else 0.0

    # Terminal-exon expression (median depth) decides boundary support.
    first_median, last_median = stats[0][1], stats[-1][1]
    if strand == "+":
        start_sup = first_median >= min_exon_coverage
        stop_sup = last_median >= min_exon_coverage
    else:
        start_sup = last_median >= min_exon_coverage
        stop_sup = first_median >= min_exon_coverage
    boundary_ratio = (int(start_sup) + int(stop_sup)) / 2.0
    return weighted, coverage_ratio, boundary_ratio


# ---------------------------------------------------------------------------
# AED formulas
# ---------------------------------------------------------------------------


def rna_aed_from_ratios(
    junction_ratio: float,
    coverage_ratio: float,
    boundary_ratio: float,
    weights: tuple[float, float, float] = (
        RNA_W_JUNCTION,
        RNA_W_COVERAGE,
        RNA_W_BOUNDARY,
    ),
) -> float:
    """v1 ``_calculate_aed`` ported verbatim: weighted ``1 - ratio`` distance.

    ``aed = w_j*(1-jr) + w_c*(1-cr) + w_b*(1-br)``, clamped to ``[0, 1]``. Lower
    is better (0 = perfect agreement, 1 = no support).
    """
    wj, wc, wb = weights
    aed = (
        wj * (1.0 - junction_ratio)
        + wc * (1.0 - coverage_ratio)
        + wb * (1.0 - boundary_ratio)
    )
    return min(1.0, max(0.0, aed))


def protein_aed_from_ratios(
    struct_ratio: float | None,
    cds_cov_ratio: float,
    prot_cov_ratio: float,
    weights: tuple[float, float, float] = (PROT_W_STRUCT, PROT_W_CDS, PROT_W_PROT),
) -> float:
    """Protein analogue of :func:`rna_aed_from_ratios`, identical distance shape.

    ``struct_ratio is None`` means the model has no intron chain to compare
    (single-exon coding). The structural term is then **dropped** and the distance
    is renormalised over the two remaining components (``w_cds`` + ``w_prot``), so a
    missing structural comparison neither credits nor penalises the AED, it never
    silently becomes ``struct_ratio = 1.0`` ("introns agree" when there are none).
    """
    ws, wc, wp = weights
    if struct_ratio is None:
        denom = wc + wp
        if denom == 0:
            return 0.0
        aed = (wc * (1.0 - cds_cov_ratio) + wp * (1.0 - prot_cov_ratio)) / denom
        return min(1.0, max(0.0, aed))
    aed = (
        ws * (1.0 - struct_ratio)
        + wc * (1.0 - cds_cov_ratio)
        + wp * (1.0 - prot_cov_ratio)
    )
    return min(1.0, max(0.0, aed))


def _introns_from_intervals(intervals: list[Any]) -> list[tuple[int, int]]:
    """Gaps between sorted, non-overlapping ``[start, end)`` intervals (introns)."""
    s = sorted(intervals, key=lambda x: x.start)
    return [
        (s[i].end, s[i + 1].start)
        for i in range(len(s) - 1)
        if s[i + 1].start > s[i].end
    ]


def _covered_fraction(
    model_intervals: list[tuple[int, int]],
    cover_intervals: list[tuple[int, int]],
) -> float:
    """Fraction of ``model_intervals`` bases overlapped by ``cover_intervals``.

    Both lists are non-overlapping (model exons/CDS and miniprot CDS segments
    respectively), so per-model-interval overlap sums without double counting.
    """
    total = sum(e - s for s, e in model_intervals)
    if total == 0:
        return 0.0
    covered = 0
    for ms, me in model_intervals:
        for cs, ce in cover_intervals:
            lo, hi = max(ms, cs), min(me, ce)
            if hi > lo:
                covered += hi - lo
    return covered / total


def score_transcript_protein(
    transcript_like: Any,
    alignment: MiniprotAlignment,
) -> dict[str, Any]:
    """protein-AED + its three component ratios for one transcript vs one alignment.

    The model's *coding intervals* are its CDS when present, else its exons (so a
    CDS-less GFF still gets a structural comparison). ``struct_ratio`` is the
    fraction of the model's (coding) introns whose donor/acceptor exactly match a
    miniprot-alignment intron; ``cds_cov_ratio`` is the fraction of the model's
    coding bases covered by the alignment's CDS segments; ``prot_cov_ratio`` is
    the alignment's reference-protein coverage (``query_coverage``).
    """
    tid, seqid, strand, exons, _start, _end = _adapt_transcript(transcript_like)
    cds = (
        transcript_like.get("cds")
        if isinstance(transcript_like, dict)
        else getattr(transcript_like, "cds", None)
    )
    coding = cds if cds else exons

    model_introns = _introns_from_intervals(coding)
    aln_introns = set(_introns_from_intervals(alignment.cds_segments))
    struct_ratio: float | None
    if not model_introns:
        # No model intron chain: structural agreement is undefined, not perfect.
        # Report it blank (``None``) rather than 1.0, which would mean "introns
        # agree"; the protein AED renormalises over the remaining components.
        struct_ratio = None
    else:
        struct_ratio = sum(1 for i in model_introns if i in aln_introns) / len(
            model_introns
        )

    model_cov_intervals = [(c.start, c.end) for c in coding]
    aln_cov_intervals = [(c.start, c.end) for c in alignment.cds_segments]
    cds_cov_ratio = _covered_fraction(model_cov_intervals, aln_cov_intervals)
    prot_cov_ratio = float(alignment.query_coverage)

    return {
        "protein_id": alignment.protein_id,
        "protein_struct_ratio": struct_ratio,
        "protein_cds_cov_ratio": cds_cov_ratio,
        "protein_prot_cov_ratio": prot_cov_ratio,
        "protein_aed": protein_aed_from_ratios(
            struct_ratio, cds_cov_ratio, prot_cov_ratio
        ),
    }


# ---------------------------------------------------------------------------
# Protein alignment loading + indexing
# ---------------------------------------------------------------------------


def _index_alignments(
    alignments: list[MiniprotAlignment],
) -> dict[str, list[MiniprotAlignment]]:
    """Bucket alignments by seqid so per-transcript best-overlap is a short scan."""
    index: dict[str, list[MiniprotAlignment]] = {}
    for a in alignments:
        index.setdefault(a.seqid, []).append(a)
    return index


def _best_alignment(
    index: dict[str, list[MiniprotAlignment]],
    seqid: str,
    start: int,
    end: int,
    strand: str,
) -> MiniprotAlignment | None:
    """Best same-strand alignment overlapping ``[start, end)`` (rank↑, then score↓)."""
    hits = [
        a
        for a in index.get(seqid, ())
        if a.strand == strand and a.start < end and a.end > start
    ]
    if not hits:
        return None
    hits.sort(key=lambda a: (a.rank, -a.score))
    return hits[0]


def load_alignments(
    *,
    proteins: Sequence[str | Path] | None = None,
    miniprot_gff: str | Path | None = None,
    genome: str | Path | None = None,
    threads: int = 1,
    miniprot_bin: str = "miniprot",
) -> list[MiniprotAlignment]:
    """Protein alignments, either parsed from a precomputed GFF or aligned fresh.

    A precomputed ``miniprot_gff`` is parsed directly (the re-run-friendly path
    and the path unit tests exercise). Otherwise ``proteins`` are aligned to
    ``genome`` with the existing :func:`helixforge.prep.protein_align.run_miniprot`
    wrapper (miniprot's own ``-t threads``); the wrapper, never reimplemented
    alignment, owns the subprocess. Multiple ``proteins`` files are concatenated
    into one query FASTA before alignment.
    """
    # Imports are local so the RNA-only path never pulls the protein I/O.
    from helixforge.io.miniprot import MiniprotParser

    if miniprot_gff is not None:
        return MiniprotParser(miniprot_gff).parse()

    proteins = list(proteins or [])
    if not proteins:
        return []
    if genome is None:
        raise ValueError("aligning --proteins from scratch requires a genome FASTA")

    import shutil
    import tempfile

    from helixforge.prep.protein_align import run_miniprot

    tmpdir = Path(tempfile.mkdtemp(prefix="helixforge_evidence_miniprot_"))
    try:
        if len(proteins) == 1:
            query = proteins[0]
        else:
            query = tmpdir / "proteins.faa"
            with open(query, "w") as out:
                for p in proteins:
                    out.write(Path(p).read_text())
        out_gff = run_miniprot(
            genome,
            query,
            tmpdir / "miniprot.gff",
            threads=threads,
            miniprot_bin=miniprot_bin,
        )
        return MiniprotParser(out_gff).parse()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Per-transcript RNA scoring (junction concordance; unchanged public behaviour)
# ---------------------------------------------------------------------------


def score_transcript_evidence(
    transcript_like: Any,
    junctions: list[SpliceJunction],
    coverage: float | None = None,
    tpm_index: dict[str, list[_StRecord]] | None = None,
    min_reads: int = DEFAULT_MIN_READS,
) -> dict[str, Any]:
    """RNA-seq junction-concordance metrics for one transcript.

    Reuses :func:`intron_concordance` (so junction logic is never duplicated) and
    returns ``{num_introns, supported, contradicted, novel_in_data,
    junction_support_fraction, intron_precision, intron_recall, intron_f1}``,
    plus ``mean_coverage`` (only when ``coverage`` is given) and ``tpm`` (only
    when ``tpm_index`` has a StringTie transcript overlapping this model, see
    :func:`_best_overlapping_tpm`).

    A single-exon transcript has no introns: counts are ``0`` and the
    fraction/precision/recall/f1 fields are ``None`` (recorded as "no introns",
    never a divide-by-zero).
    """
    tid, seqid, strand, exons, start, end = _adapt_transcript(transcript_like)

    tx = _TxShim(tid or "tx", seqid, strand, exons)
    gene = _GeneShim(tid or "gene", seqid, strand, start, end, [tx])
    ic = intron_concordance(gene, junctions, min_reads=min_reads)["isoforms"][
        tx.transcript_id
    ]

    num = ic["num_introns"]
    jsf = ic["supported"] / num if num else None
    # A transcript with no introns has no intron chain to score: every intron
    # metric must be blank. ``intron_concordance`` already returns precision/f1 as
    # ``None`` here, but it returns ``recall == 0.0`` (not ``None``) whenever a
    # qualifying junction merely overlaps the locus: matched/qualifying = 0/k.
    # That leaks a meaningless 0.0 into ``intron_recall`` for single-exon models.
    # Force all four to ``None`` so "blank iff num_introns == 0" holds exactly.
    if num:
        precision, recall, f1 = ic["precision"], ic["recall"], ic["f1"]
    else:
        precision = recall = f1 = None
    result = {
        "num_introns": num,
        "supported": ic["supported"],
        "contradicted": ic["contradicted"],
        "novel_in_data": ic["novel"],
        "junction_support_fraction": jsf,
        "intron_precision": precision,
        "intron_recall": recall,
        "intron_f1": f1,
    }

    if coverage is not None:
        result["mean_coverage"] = coverage
    if tpm_index:
        tpm = _best_overlapping_tpm(tpm_index, seqid, strand, exons)
        if tpm is not None:
            result["tpm"] = tpm
    return result


# ---------------------------------------------------------------------------
# Process-pool worker (merged evidence broadcast once via the initializer)
# ---------------------------------------------------------------------------

_W_JUNCTIONS: list[Any] | None = None
_W_COVERAGE: dict[tuple[str, int, int], tuple[float, float, float]] | None = None
_W_ALN_INDEX: dict[str, list[Any]] | None = None
_W_TPM: dict[str, list[_StRecord]] | None = None
_W_PARAMS: dict[str, Any] = {}


def _init_worker(
    junctions: list[Any],
    coverage: dict[tuple[str, int, int], tuple[float, float, float]] | None,
    aln_index: dict[str, list[Any]] | None,
    tpm: dict[str, list[_StRecord]] | None,
    params: dict[str, Any],
) -> None:
    global _W_JUNCTIONS, _W_COVERAGE, _W_ALN_INDEX, _W_TPM, _W_PARAMS
    _W_JUNCTIONS = junctions
    _W_COVERAGE = coverage
    _W_ALN_INDEX = aln_index
    _W_TPM = tpm
    _W_PARAMS = params


def _score_job(job: dict[str, Any]) -> dict[str, Any]:
    """Score one transcript into a full TSV row using the broadcast evidence."""
    seqid = job["seqid"]
    strand = job["strand"]
    exons = job["exons"]
    t_start, t_end = job["start"], job["end"]

    has_cov = bool(_W_PARAMS["has_coverage"]) and _W_COVERAGE is not None
    coverage_mean: float | None = None
    cov_ratio: float | None = None
    bound_ratio: float | None = None
    if has_cov:
        assert _W_COVERAGE is not None  # narrowed by has_cov
        coverage_mean, cov_ratio, bound_ratio = _transcript_coverage_metrics(
            seqid, strand, exons, _W_COVERAGE, _W_PARAMS["min_exon_coverage"]
        )

    transcript_like = {
        "transcript_id": job["transcript_id"],
        "seqid": seqid,
        "strand": strand,
        "exons": exons,
        "start": t_start,
        "end": t_end,
    }
    metrics = score_transcript_evidence(
        transcript_like,
        _W_JUNCTIONS or [],
        coverage=coverage_mean,
        tpm_index=_W_TPM,
        min_reads=_W_PARAMS["min_reads"],
    )

    row: dict[str, Any] = {
        "gene_id": job["gene_id"],
        "transcript_id": job["transcript_id"],
        "seqid": seqid,
        "strand": strand,
        "start": t_start,
        "end": t_end,
        "num_exons": len(exons),
        "mean_coverage": metrics.get("mean_coverage"),
        "tpm": metrics.get("tpm"),
        **{k: v for k, v in metrics.items() if k not in ("mean_coverage", "tpm")},
    }

    # RNA-AED: only meaningful with BAM coverage (needs all three ratios).
    if has_cov and cov_ratio is not None and bound_ratio is not None:
        num = metrics["num_introns"]
        jr = metrics["supported"] / num if num else 1.0  # v1: single-exon -> 1.0
        row["rna_junction_ratio"] = jr
        row["rna_coverage_ratio"] = cov_ratio
        row["rna_boundary_ratio"] = bound_ratio
        row["rna_aed"] = rna_aed_from_ratios(
            jr, cov_ratio, bound_ratio, _W_PARAMS["rna_weights"]
        )

    # protein-AED: only when a proteome / miniprot GFF was supplied.
    if _W_ALN_INDEX is not None:
        best = _best_alignment(_W_ALN_INDEX, seqid, t_start, t_end, strand)
        if best is not None:
            row.update(
                score_transcript_protein(
                    {**transcript_like, "cds": job.get("cds")}, best
                )
            )
    return row


# ---------------------------------------------------------------------------
# Whole-annotation scoring
# ---------------------------------------------------------------------------


def _build_tpm_index(
    stringtie_gtfs: Sequence[str | Path],
) -> dict[str, list[_StRecord]]:
    """``seqid -> [(strand, start, end, exon_tuples, max_tpm), ...]`` for overlap lookup.

    Thin wrapper over :func:`io.stringtie.build_tpm_overlap_index` (the shared
    implementation used by ``reconcile`` too); see :func:`_best_overlapping_tpm`
    for the matching rule.
    """
    return build_tpm_overlap_index(stringtie_gtfs)


def score_annotation(
    gff3_path: str | Path,
    *,
    bam_paths: Sequence[str | Path] | None = None,
    star_sj_paths: Sequence[str | Path] | None = None,
    stringtie_gtfs: Sequence[str | Path] | None = None,
    proteins: Sequence[str | Path] | None = None,
    miniprot_gff: str | Path | None = None,
    genome: str | Path | None = None,
    region: str | None = None,
    min_reads: int = DEFAULT_MIN_READS,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    min_overhang: int = DEFAULT_MIN_OVERHANG,
    min_exon_coverage: int = DEFAULT_MIN_EXON_COVERAGE,
    rna_weights: tuple[float, float, float] = (
        RNA_W_JUNCTION,
        RNA_W_COVERAGE,
        RNA_W_BOUNDARY,
    ),
    protein_weights: tuple[float, float, float] = (
        PROT_W_STRUCT,
        PROT_W_CDS,
        PROT_W_PROT,
    ),
    reference_filename: str | Path | None = None,
    miniprot_bin: str = "miniprot",
    threads: int = 1,
) -> pd.DataFrame:
    """Score every transcript of any GFF3 against the supplied evidence.

    Returns a per-transcript ``pandas.DataFrame`` (columns :data:`TSV_COLUMNS`).
    ``gff3_path`` need not be HelixForge output, it is parsed with
    ``GFF3Parser.parse_genes_generic`` so any GFF3 works. Junctions are unioned
    from ``bam_paths`` / ``star_sj_paths``; ``tpm`` from ``stringtie_gtfs`` by
    exonic overlap; per-transcript coverage + ``rna_aed`` from BAMs;
    ``protein_aed`` from ``miniprot_gff`` (precomputed) or ``proteins`` aligned to
    ``genome``. Each AED block is ``NA`` when its evidence was not supplied
    (missing evidence is neutral). ``region`` (1-based ``seqid`` or
    ``seqid:start-end``) restricts the scored transcripts.

    ``threads`` is the single parallelism dial: extraction (one BAM/SJ per worker:
    a single count_coverage pass per locus + the junction scan; plus the miniprot
    run) runs in a process pool and the merged evidence is built once; per-gene
    scoring then runs in a process pool over gene batches. No BAM is opened or
    queried inside the per-transcript loop. Output is byte-identical to
    ``threads=1`` (rows are sorted by gene/transcript).
    """
    import pandas as pd

    bam_paths = list(bam_paths or [])
    star_sj_paths = list(star_sj_paths or [])
    threads = max(1, int(threads))

    # ---- 1. Parse the annotation + build the (region-filtered) job list. ----
    r_seqid, r_lo, r_hi = _region_to_internal(region)
    genes = GFF3Parser(gff3_path).parse_genes_generic()
    jobs: list[dict[str, Any]] = []
    exon_region_set: set[tuple[str, int, int]] = set()
    for gene in genes:
        if r_seqid is not None and gene["seqid"] != r_seqid:
            continue
        for tx in gene["transcripts"]:
            exons = tx["exons"]
            if not exons:
                continue
            t_start = min(e.start for e in exons)
            t_end = max(e.end for e in exons)
            if r_lo is not None and not (t_start < r_hi and t_end > r_lo):
                continue
            jobs.append(
                {
                    "gene_id": gene["gene_id"],
                    "transcript_id": tx["transcript_id"],
                    "seqid": gene["seqid"],
                    "strand": gene["strand"],
                    "exons": exons,
                    "cds": tx.get("cds"),
                    "start": t_start,
                    "end": t_end,
                }
            )
            if bam_paths:
                for e in exons:
                    exon_region_set.add((gene["seqid"], e.start, e.end))

    # ---- 2. Extraction (process-parallel; merged evidence built once). ----
    junctions: list[Any] = []
    if bam_paths or star_sj_paths:
        junctions = collect_junctions(
            bam_paths=bam_paths or None,
            star_sj_paths=star_sj_paths or None,
            region=region,
            min_reads=min_reads,
            min_mapq=min_mapq,
            min_overhang=min_overhang,
            reference_filename=reference_filename,
            threads=threads,
        )

    coverage_cache: dict[tuple[str, int, int], tuple[float, float, float]] | None = None
    if bam_paths:
        coverage_cache = collect_coverage(
            bam_paths,
            sorted(exon_region_set),
            reference_filename=reference_filename,
            threads=threads,
        )

    tpm_index = _build_tpm_index(stringtie_gtfs) if stringtie_gtfs else None

    aln_index: dict[str, list[Any]] | None = None
    if proteins or miniprot_gff:
        alignments = load_alignments(
            proteins=proteins,
            miniprot_gff=miniprot_gff,
            genome=genome,
            threads=threads,
            miniprot_bin=miniprot_bin,
        )
        aln_index = _index_alignments(alignments)

    # ---- 3. Scoring (process-parallel over gene batches). ----
    params = {
        "has_coverage": bool(bam_paths),
        "min_exon_coverage": min_exon_coverage,
        "min_reads": min_reads,
        "rna_weights": rna_weights,
        "protein_weights": protein_weights,
    }
    if threads > 1 and len(jobs) > 1:
        chunksize = max(1, len(jobs) // (threads * 4))
        with ProcessPoolExecutor(
            max_workers=threads,
            initializer=_init_worker,
            initargs=(junctions, coverage_cache, aln_index, tpm_index, params),
        ) as ex:
            rows = list(ex.map(_score_job, jobs, chunksize=chunksize))
    else:
        _init_worker(junctions, coverage_cache, aln_index, tpm_index, params)
        rows = [_score_job(job) for job in jobs]

    # ---- 4. Deterministic output (sort so threads==serial). ----
    rows.sort(key=lambda r: (r["gene_id"], r["transcript_id"]))
    return pd.DataFrame(rows, columns=TSV_COLUMNS)


# ---------------------------------------------------------------------------
# Per-gene rollup
# ---------------------------------------------------------------------------


def _best_transcript_key(row: "pd.Series") -> tuple[float, float, float, str]:
    """Sort key picking a gene's best-supported transcript (smallest wins).

    Ranking, in order: lowest ``rna_aed`` (best RNA agreement), then lowest
    ``protein_aed``, then highest ``junction_support_fraction`` (more confirmed
    introns), then ``transcript_id`` for a deterministic tie-break. Absent AEDs
    sort last so a transcript with evidence always beats one without.
    """
    import math

    def _na(v: Any, worst: float) -> float:
        return (
            worst if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)
        )

    rna = _na(row.get("rna_aed"), math.inf)
    prot = _na(row.get("protein_aed"), math.inf)
    jsf = _na(row.get("junction_support_fraction"), -1.0)
    return (rna, prot, -jsf, str(row["transcript_id"]))


def rollup_genes(df: pd.DataFrame) -> pd.DataFrame:
    """One row per gene: its best-supported transcript (see :func:`_best_transcript_key`).

    The schema is identical to the per-transcript table (:data:`TSV_COLUMNS`); the
    chosen transcript's row is carried verbatim. Rows are sorted by ``gene_id``.
    """
    import pandas as pd

    if df.empty:
        return pd.DataFrame(columns=TSV_COLUMNS)
    best_rows = []
    for _gene_id, group in df.groupby("gene_id", sort=True):
        ranked = sorted((row for _, row in group.iterrows()), key=_best_transcript_key)
        best_rows.append(ranked[0])
    out = pd.DataFrame(best_rows, columns=TSV_COLUMNS).reset_index(drop=True)
    return out.sort_values("gene_id").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Summary + writers
# ---------------------------------------------------------------------------


def summarize_evidence(df: pd.DataFrame) -> dict[str, Any]:
    """Distribution summary over a scored DataFrame.

    Returns ``{n_transcripts, n_multi_exon, n_single_exon, n_fully_supported,
    fraction_fully_supported, mean_junction_support_fraction, mean_intron_f1,
    median_intron_f1}`` plus ``mean_rna_aed`` / ``mean_protein_aed``.

    Each AED mean is averaged over a **different** denominator and that count is
    reported alongside it so the two are never silently compared on unlike sets:

    - ``mean_rna_aed`` is over the ``n_rna_aed`` transcripts that have an RNA AED
      (every transcript, once BAM coverage was supplied), and
    - ``mean_protein_aed`` is over only the ``n_protein_aed`` transcripts with a
      protein hit, a strict subset.

    Comparing the two means directly is comparing all transcripts against the
    protein-hit subset; read each next to its ``n_…`` count. "Fully supported" and
    the junction/F1 distributions are over **multi-exon** transcripts only.
    """
    n = int(len(df))
    multi = df[df["num_introns"] > 0]
    n_multi = int(len(multi))
    n_single = n - n_multi

    def _present(col: str) -> "pd.Series":
        import pandas as pd

        if col not in df.columns:
            return pd.Series(dtype=float)
        return df[col].dropna()

    def _mean(vals: "pd.Series") -> float | None:
        return float(vals.mean()) if len(vals) else None

    rna_vals = _present("rna_aed")
    prot_vals = _present("protein_aed")
    base: dict[str, Any] = {
        "n_transcripts": n,
        "n_multi_exon": n_multi,
        "n_single_exon": n_single,
        "mean_rna_aed": _mean(rna_vals),
        "n_rna_aed": int(len(rna_vals)),
        "mean_protein_aed": _mean(prot_vals),
        "n_protein_aed": int(len(prot_vals)),
    }

    if n_multi == 0:
        base.update(
            {
                "n_fully_supported": 0,
                "fraction_fully_supported": None,
                "mean_junction_support_fraction": None,
                "mean_intron_f1": None,
                "median_intron_f1": None,
            }
        )
        return base

    jsf = multi["junction_support_fraction"]
    f1 = multi["intron_f1"].dropna()
    n_full = int((jsf >= 1.0).sum())
    base.update(
        {
            "n_fully_supported": n_full,
            "fraction_fully_supported": n_full / n_multi,
            "mean_junction_support_fraction": float(jsf.mean()),
            "mean_intron_f1": float(f1.mean()) if len(f1) else None,
            "median_intron_f1": float(f1.median()) if len(f1) else None,
        }
    )
    return base


def write_evidence_tsv(df: pd.DataFrame, path: str | Path) -> Path:
    """Write the scored DataFrame to a TSV; return ``Path``."""
    path = Path(path)
    df.to_csv(path, sep="\t", index=False)
    return path
