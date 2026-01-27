"""Output writers for RNA-seq evidence scoring.

This module provides functions for exporting evidence scores in various formats:

- TSV reports with detailed per-gene metrics
- Summary statistics
- GFF3 attribute updates

Example:
    >>> from helixforge.core.evidence import EvidenceScorer
    >>> from helixforge.core.evidence_output import (
    ...     write_evidence_report_tsv,
    ...     write_gene_evidence_summary_tsv,
    ...     update_gff_with_evidence,
    ... )
    >>>
    >>> scores = list(scorer.score_genes(genes, junctions))
    >>> write_evidence_report_tsv(scores, "evidence_report.tsv")
    >>> write_gene_evidence_summary_tsv(scores, "evidence_summary.tsv")
    >>> update_gff_with_evidence(genes, scores, "annotated.gff3")
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from helixforge.core.evidence import EvidenceScore
    from helixforge.io.gff import GeneModel

logger = logging.getLogger(__name__)


# =============================================================================
# TSV Report Writers
# =============================================================================


def write_evidence_report_tsv(
    scores: list["EvidenceScore"],
    output_path: Path | str,
    include_details: bool = False,
) -> None:
    """Write detailed evidence report to TSV.

    Args:
        scores: List of EvidenceScore objects.
        output_path: Output file path.
        include_details: Include per-junction and per-exon details.
    """
    output_path = Path(output_path)

    headers = [
        "gene_id",
        "transcript_id",
        "seqid",
        "strand",
        "n_introns",
        "n_introns_supported",
        "n_introns_exact",
        "junction_support_ratio",
        "n_exons",
        "n_exons_expressed",
        "mean_exon_coverage",
        "exon_coverage_ratio",
        "start_supported",
        "stop_supported",
        "strand_consistent",
        "evidence_level",
        "aed",
        "flags",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)

        for score in scores:
            row = [
                score.gene_id,
                score.transcript_id,
                score.seqid,
                score.strand,
                score.n_introns,
                score.n_introns_supported,
                score.n_introns_exact,
                f"{score.junction_support_ratio:.4f}",
                score.n_exons,
                score.n_exons_expressed,
                f"{score.mean_exon_coverage:.2f}",
                f"{score.exon_coverage_ratio:.4f}",
                "yes" if score.start_supported else "no",
                "yes" if score.stop_supported else "no",
                "yes" if score.strand_consistent else "no",
                score.evidence_level.value,
                f"{score.aed:.4f}",
                ",".join(score.flags) if score.flags else "",
            ]
            writer.writerow(row)

    logger.info(f"Wrote evidence report with {len(scores)} genes to {output_path}")


def write_gene_evidence_summary_tsv(
    scores: list["EvidenceScore"],
    output_path: Path | str,
) -> None:
    """Write simplified summary TSV for downstream processing.

    This format is designed for compatibility with QC aggregation.

    Args:
        scores: List of EvidenceScore objects.
        output_path: Output file path.
    """
    output_path = Path(output_path)

    headers = [
        "gene_id",
        "evidence_level",
        "aed",
        "junction_support_ratio",
        "exon_coverage_ratio",
        "n_flags",
        "flags",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)

        for score in scores:
            row = [
                score.gene_id,
                score.evidence_level.value,
                f"{score.aed:.4f}",
                f"{score.junction_support_ratio:.4f}",
                f"{score.exon_coverage_ratio:.4f}",
                len(score.flags),
                ",".join(score.flags) if score.flags else "",
            ]
            writer.writerow(row)

    logger.info(f"Wrote evidence summary with {len(scores)} genes to {output_path}")


def write_junction_details_tsv(
    scores: list["EvidenceScore"],
    output_path: Path | str,
) -> None:
    """Write per-junction details to TSV.

    Args:
        scores: List of EvidenceScore objects.
        output_path: Output file path.
    """
    output_path = Path(output_path)

    headers = [
        "gene_id",
        "transcript_id",
        "intron_start",
        "intron_end",
        "intron_length",
        "observed",
        "read_count",
        "unique_count",
        "exact_match",
        "shift_donor",
        "shift_acceptor",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)

        for score in scores:
            for je in score.junction_evidence:
                row = [
                    score.gene_id,
                    score.transcript_id,
                    je.intron_start,
                    je.intron_end,
                    je.intron_length,
                    "yes" if je.observed else "no",
                    je.read_count,
                    je.unique_count,
                    "yes" if je.exact_match else "no",
                    je.shift_donor,
                    je.shift_acceptor,
                ]
                writer.writerow(row)

    n_junctions = sum(len(s.junction_evidence) for s in scores)
    logger.info(f"Wrote junction details with {n_junctions} junctions to {output_path}")


def write_exon_details_tsv(
    scores: list["EvidenceScore"],
    output_path: Path | str,
) -> None:
    """Write per-exon details to TSV.

    Args:
        scores: List of EvidenceScore objects.
        output_path: Output file path.
    """
    output_path = Path(output_path)

    headers = [
        "gene_id",
        "transcript_id",
        "exon_start",
        "exon_end",
        "exon_length",
        "mean_coverage",
        "median_coverage",
        "min_coverage",
        "max_coverage",
        "fraction_covered",
        "coverage_uniformity",
        "is_expressed",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)

        for score in scores:
            for ee in score.exon_evidence:
                row = [
                    score.gene_id,
                    score.transcript_id,
                    ee.exon_start,
                    ee.exon_end,
                    ee.exon_length,
                    f"{ee.mean_coverage:.2f}",
                    f"{ee.median_coverage:.2f}",
                    ee.min_coverage,
                    ee.max_coverage,
                    f"{ee.fraction_covered:.4f}",
                    f"{ee.coverage_uniformity:.4f}",
                    "yes" if ee.is_expressed else "no",
                ]
                writer.writerow(row)

    n_exons = sum(len(s.exon_evidence) for s in scores)
    logger.info(f"Wrote exon details with {n_exons} exons to {output_path}")


# =============================================================================
# GFF3 Output
# =============================================================================


def update_gff_with_evidence(
    genes: list["GeneModel"],
    scores: list["EvidenceScore"],
    output_path: Path | str,
    source: str = "HelixForge",
) -> None:
    """Write GFF3 with evidence attributes added.

    Adds the following attributes to gene features:
    - evidence_level: full/partial/minimal/none
    - aed: Annotation Edit Distance (0-1)
    - junction_support: fraction of supported junctions
    - exon_coverage: fraction of expressed exons
    - evidence_flags: comma-separated list of flags

    Args:
        genes: List of GeneModel objects.
        scores: List of EvidenceScore objects (must match genes).
        output_path: Output GFF3 file path.
        source: Source field for GFF3.
    """
    from helixforge.io.gff import GFF3Writer

    # Build score lookup
    score_lookup = {s.gene_id: s for s in scores}

    with GFF3Writer(output_path, source=source) as writer:
        writer.write_header()

        for gene in genes:
            score = score_lookup.get(gene.gene_id)

            if score is not None:
                # Add evidence attributes
                gene.attributes["evidence_level"] = score.evidence_level.value
                gene.attributes["aed"] = f"{score.aed:.4f}"
                gene.attributes["junction_support"] = f"{score.junction_support_ratio:.4f}"
                gene.attributes["exon_coverage"] = f"{score.exon_coverage_ratio:.4f}"
                if score.flags:
                    gene.attributes["evidence_flags"] = ",".join(score.flags)

            writer.write_gene(gene)

    logger.info(f"Wrote GFF3 with evidence for {len(genes)} genes to {output_path}")


def update_genes_with_evidence(
    genes: list["GeneModel"],
    scores: list["EvidenceScore"],
) -> list["GeneModel"]:
    """Update gene models in-place with evidence attributes.

    This modifies the gene objects directly rather than writing to file.

    Args:
        genes: List of GeneModel objects.
        scores: List of EvidenceScore objects.

    Returns:
        The same list of genes with updated attributes.
    """
    score_lookup = {s.gene_id: s for s in scores}

    for gene in genes:
        score = score_lookup.get(gene.gene_id)
        if score is not None:
            gene.attributes["evidence_level"] = score.evidence_level.value
            gene.attributes["aed"] = f"{score.aed:.4f}"
            gene.attributes["junction_support"] = f"{score.junction_support_ratio:.4f}"
            gene.attributes["exon_coverage"] = f"{score.exon_coverage_ratio:.4f}"
            if score.flags:
                gene.attributes["evidence_flags"] = ",".join(score.flags)

    return genes


# =============================================================================
# Summary Statistics
# =============================================================================


def generate_evidence_summary_text(
    scores: list["EvidenceScore"],
) -> str:
    """Generate human-readable summary of evidence scores.

    Args:
        scores: List of EvidenceScore objects.

    Returns:
        Formatted summary text.
    """
    from helixforge.core.evidence import EvidenceLevel, summarize_evidence_scores

    if not scores:
        return "No genes to summarize."

    summary = summarize_evidence_scores(scores)

    lines = [
        "Evidence Score Summary",
        "=" * 40,
        "",
        f"Total genes:           {summary['n_genes']:,}",
        "",
        "Evidence Level Distribution:",
        f"  Full support:        {summary['n_full_support']:,} ({summary['n_full_support']/summary['n_genes']*100:.1f}%)",
        f"  Partial support:     {summary['n_partial_support']:,} ({summary['n_partial_support']/summary['n_genes']*100:.1f}%)",
        f"  Minimal support:     {summary['n_minimal_support']:,} ({summary['n_minimal_support']/summary['n_genes']*100:.1f}%)",
        f"  No support:          {summary['n_no_support']:,} ({summary['n_no_support']/summary['n_genes']*100:.1f}%)",
        "",
        "AED Statistics:",
        f"  Mean AED:            {summary['mean_aed']:.4f}",
        f"  Median AED:          {summary['median_aed']:.4f}",
        "",
    ]

    # Additional stats
    n_single_exon = sum(1 for s in scores if s.is_single_exon)
    n_full_junction = sum(1 for s in scores if s.has_full_junction_support and not s.is_single_exon)
    n_multi_exon = sum(1 for s in scores if not s.is_single_exon)

    lines.extend([
        "Gene Structure:",
        f"  Single-exon genes:   {n_single_exon:,} ({n_single_exon/summary['n_genes']*100:.1f}%)",
        f"  Multi-exon genes:    {n_multi_exon:,} ({n_multi_exon/summary['n_genes']*100:.1f}%)",
        "",
    ])

    if n_multi_exon > 0:
        lines.append(
            f"  Full junction support (multi-exon): {n_full_junction:,} "
            f"({n_full_junction/n_multi_exon*100:.1f}%)"
        )

    return "\n".join(lines)


def write_evidence_summary_txt(
    scores: list["EvidenceScore"],
    output_path: Path | str,
) -> None:
    """Write human-readable summary to text file.

    Args:
        scores: List of EvidenceScore objects.
        output_path: Output file path.
    """
    output_path = Path(output_path)
    text = generate_evidence_summary_text(scores)

    with open(output_path, "w") as f:
        f.write(text)
        f.write("\n")

    logger.info(f"Wrote evidence summary to {output_path}")


# =============================================================================
# Loading Functions
# =============================================================================


def load_evidence_report_tsv(
    input_path: Path | str,
) -> list["EvidenceScore"]:
    """Load evidence scores from TSV report.

    Args:
        input_path: Path to evidence report TSV.

    Returns:
        List of EvidenceScore objects.
    """
    from helixforge.core.evidence import EvidenceLevel, EvidenceScore

    input_path = Path(input_path)
    scores = []

    with open(input_path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            flags = row.get("flags", "").split(",") if row.get("flags") else []

            score = EvidenceScore(
                gene_id=row["gene_id"],
                transcript_id=row["transcript_id"],
                seqid=row["seqid"],
                strand=row["strand"],
                n_introns=int(row["n_introns"]),
                n_introns_supported=int(row["n_introns_supported"]),
                n_introns_exact=int(row["n_introns_exact"]),
                junction_support_ratio=float(row["junction_support_ratio"]),
                n_exons=int(row["n_exons"]),
                n_exons_expressed=int(row["n_exons_expressed"]),
                mean_exon_coverage=float(row["mean_exon_coverage"]),
                exon_coverage_ratio=float(row["exon_coverage_ratio"]),
                start_supported=row["start_supported"].lower() == "yes",
                stop_supported=row["stop_supported"].lower() == "yes",
                strand_consistent=row["strand_consistent"].lower() == "yes",
                evidence_level=EvidenceLevel(row["evidence_level"]),
                aed=float(row["aed"]),
                flags=flags,
            )
            scores.append(score)

    logger.info(f"Loaded {len(scores)} evidence scores from {input_path}")
    return scores
