"""Expression classification of Helixer loci."""

from __future__ import annotations

from contextlib import nullcontext
from typing import TYPE_CHECKING, cast

from helixforge.io.bam import CoverageCalculator, CoveragePool
from helixforge.io.stringtie import StringTieParser
from helixforge.reconcile.models import LocusClassification
from helixforge.utils.intervals import IntervalIndex

if TYPE_CHECKING:
    from helixforge.io.bam import CoverageCalculator
    from helixforge.reconcile.models import HelixerLocus, StringTieTranscript


def classify_loci(
    loci: list[HelixerLocus],
    stringtie_aggregated: dict[str, dict[str, object]] | None = None,
    stringtie_transcripts: list[StringTieTranscript] | None = None,
    bam_paths: list[str] | None = None,
    bigwig_paths: list[str] | None = None,
    min_tpm: float = 0.5,
    min_samples: int = 1,
    coverage_threshold: float = 2.0,
    near_zero_coverage: float = 0.1,
) -> list[LocusClassification]:
    """Classify each locus EXPRESSED / LOW / SILENT (same order as input).

    StringTie is tried first; if a locus has no same-strand StringTie overlap,
    coverage is used when BAM/bigWig sources are given, otherwise the locus is
    SILENT with ``evidence_source='none'``.
    """
    results: list[LocusClassification] = []
    st_indices: dict[str, IntervalIndex] = {}  # seqid -> IntervalIndex (built lazily)

    # Open each coverage source ONCE for the whole pass instead of per locus.
    # bigWig is preferred over BAM, so the
    # pool only opens the type that will actually be queried. The context manager
    # closes every handle deterministically when the pass ends.
    use_bigwig = bool(bigwig_paths)
    pool_cm = (
        CoveragePool(
            bam_paths=None if use_bigwig else bam_paths,
            bigwig_paths=bigwig_paths if use_bigwig else None,
            calculator_cls=CoverageCalculator,
        )
        if (bam_paths or bigwig_paths)
        else nullcontext()
    )

    with pool_cm as pool:
        for locus in loci:
            classification: LocusClassification | None = None

            if stringtie_transcripts:
                idx = st_indices.get(locus.seqid)
                if idx is None:
                    idx = _build_stringtie_index(stringtie_transcripts, locus.seqid)
                    st_indices[locus.seqid] = idx
                classification = _classify_by_stringtie(
                    locus,
                    stringtie_transcripts,
                    idx,
                    min_tpm,
                    min_samples,
                    stringtie_aggregated,
                )

            if classification is None:
                if bam_paths or bigwig_paths:
                    assert isinstance(pool, CoveragePool)
                    classification = _classify_by_coverage(
                        locus,
                        pool,
                        coverage_threshold,
                        near_zero_coverage,
                    )
                else:
                    classification = LocusClassification(
                        locus.gene_id, "SILENT", evidence_source="none"
                    )

            results.append(classification)
    return results


def _build_stringtie_index(
    transcripts: list[StringTieTranscript],
    seqid: str,
) -> IntervalIndex:
    """Index transcripts on ``seqid``; payload is the global transcript index."""
    idx = IntervalIndex()
    idx.add_intervals(
        [(t.start, t.end, i) for i, t in enumerate(transcripts) if t.seqid == seqid]
    )
    return idx


def _classify_by_stringtie(
    locus: HelixerLocus,
    transcripts: list[StringTieTranscript],
    st_index: IntervalIndex,
    min_tpm: float,
    min_samples: int,
    aggregated: dict[str, dict[str, object]] | None,
) -> LocusClassification | None:
    """Classify via same-strand StringTie overlap, or ``None`` to fall back.

    Returns EXPRESSED if any overlapping same-strand structure meets the TPM and
    sample thresholds, LOW if there is overlap but none qualify, or ``None`` if
    there is no same-strand overlap at all.
    """
    hits = st_index.query_with_data(locus.start, locus.end)
    overlappers = [
        transcripts[data[2]]
        for data in hits
        if transcripts[data[2]].strand == locus.strand
    ]
    if not overlappers:
        return None

    expressed = False
    best_tpm = -1.0
    best_samples = 0
    for t in overlappers:
        key = StringTieParser._structure_key(t)
        if aggregated and key in aggregated:
            tpm = cast(float, aggregated[key]["max_tpm"])
            num_samples = cast(int, aggregated[key]["num_samples"])
        else:
            tpm = t.tpm
            num_samples = 1
        if tpm >= min_tpm and num_samples >= min_samples:
            expressed = True
        if tpm > best_tpm:
            best_tpm = tpm
            best_samples = num_samples

    status = "EXPRESSED" if expressed else "LOW"
    return LocusClassification(
        locus_id=locus.gene_id,
        status=status,
        max_tpm=best_tpm,
        evidence_source="stringtie",
        num_samples_expressed=best_samples,
    )


def _classify_by_coverage(
    locus: HelixerLocus,
    pool: CoveragePool,
    coverage_threshold: float,
    near_zero_coverage: float,
) -> LocusClassification:
    """Classify via weighted exonic coverage; bigWig is preferred over BAM.

    ``pool`` is the pass-scoped :class:`~helixforge.io.bam.CoveragePool` of already
    open coverage handles, no per-locus open happens here.
    """
    if pool.bigwig_paths:
        calculators, evidence = pool.calculators("bigwig"), "bigwig"
    else:
        calculators, evidence = pool.calculators("bam"), "bam_coverage"

    coverage = _get_exonic_coverage(locus, calculators)
    if coverage >= coverage_threshold:
        status = "EXPRESSED"
    elif coverage >= near_zero_coverage:
        status = "LOW"
    else:
        status = "SILENT"
    return LocusClassification(
        locus_id=locus.gene_id,
        status=status,
        mean_coverage=coverage,
        evidence_source=evidence,
    )


def _get_exonic_coverage(
    locus: HelixerLocus,
    calculators: list[CoverageCalculator],
) -> float:
    """Exon-length-weighted mean coverage across open handles (span if no exons).

    ``calculators`` is the list of already-open
    :class:`~helixforge.io.bam.CoverageCalculator` handles for the chosen source
    type; this function never opens or closes a source.
    """
    regions: list[tuple[int, int]] = (
        [(e.start, e.end) for e in locus.exons]
        if locus.exons
        else [(locus.start, locus.end)]
    )
    total_len = sum(e - s for s, e in regions)
    if total_len == 0:
        return 0.0

    per_source: list[float] = []
    for cov in calculators:
        weighted = 0.0
        for s, e in regions:
            weighted += cov.mean_coverage(locus.seqid, s, e) * (e - s)
        per_source.append(weighted / total_len)

    if not per_source:
        return 0.0
    return sum(per_source) / len(per_source)
