"""Export writers for downstream tools."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Any

from helixforge.io.gff import GFF3Writer
from helixforge.reconcile.validate import extract_cds_sequence
from helixforge.utils.atomic import atomic_write
from helixforge.utils.logging import get_logger
from helixforge.utils.regions import internal_to_gff3
from helixforge.utils.sequences import translate

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene

_log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Sequence helpers (coding direction, both strands)
# ---------------------------------------------------------------------------


def _cdna_sequence(tx: Any, genome: Any) -> str | None:
    """Spliced cDNA in coding (5'→3') orientation.

    Exons are stored low→high genomically; on ``-`` the coding order is
    descending and each piece is reverse-complemented by ``get_sequence``. Returns
    ``None`` if any exon cannot be read.
    """
    exons = list(tx.exons)
    if tx.strand == "-":
        exons = list(reversed(exons))
    pieces: list[str] = []
    for ex in exons:
        try:
            piece = genome.get_sequence(tx.seqid, ex.start, ex.end, tx.strand)
        except (KeyError, IndexError, ValueError) as exc:
            # A missing/renamed scaffold or an out-of-range exon previously
            # yielded a silently shorter FASTA. Log the context so a dropped
            # transcript is visible, then skip it (behaviour preserved).
            _log.warning(
                "skipping %s: cannot read exon %s:%d-%d (%s)",
                tx.transcript_id,
                tx.seqid,
                ex.start,
                ex.end,
                exc,
            )
            return None
        pieces.append(piece)
    return "".join(pieces)


def _protein_sequence(tx: Any, genome: Any) -> str | None:
    """Translate the CDS (coding direction), dropping a trailing stop. ``None`` if no CDS."""
    cds_seq = extract_cds_sequence(tx, genome)
    if cds_seq is None:
        return None
    first = tx.cds[0] if tx.strand == "+" else tx.cds[-1]
    return translate(cds_seq, first.phase).rstrip("*")


def _write_fasta(records: list[tuple[str, str]], path: str, width: int = 60) -> str:
    with atomic_write(path) as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
    return str(path)


def write_protein_fasta(genes: list[ReconciledGene], genome: Any, path: str) -> str:
    """Protein FASTA for every coding isoform (BUSCO/compleasm/OMArk/DIAMOND)."""
    records: list[tuple[str, str]] = []
    for gene in genes:
        for tx in gene.transcripts:
            if not tx.cds:
                continue
            prot = _protein_sequence(tx, genome)
            if prot:
                records.append((tx.transcript_id, prot))
    return _write_fasta(records, path)


def write_cds_fasta(genes: list[ReconciledGene], genome: Any, path: str) -> str:
    """Nucleotide CDS FASTA (coding direction) for every coding isoform."""
    records: list[tuple[str, str]] = []
    for gene in genes:
        for tx in gene.transcripts:
            if not tx.cds:
                continue
            cds_seq = extract_cds_sequence(tx, genome)
            if cds_seq:
                records.append((tx.transcript_id, cds_seq))
    return _write_fasta(records, path)


def write_cdna_fasta(genes: list[ReconciledGene], genome: Any, path: str) -> str:
    """Spliced cDNA FASTA (coding direction) for every isoform."""
    records: list[tuple[str, str]] = []
    for gene in genes:
        for tx in gene.transcripts:
            seq = _cdna_sequence(tx, genome)
            if seq:
                records.append((tx.transcript_id, seq))
    return _write_fasta(records, path)


# ---------------------------------------------------------------------------
# GTF (legacy tools)
# ---------------------------------------------------------------------------


def _gtf_attrs(
    gene_id: str,
    transcript_id: str,
    gene_biotype: str | None = None,
    transcript_biotype: str | None = None,
) -> str:
    out = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
    if gene_biotype:
        out += f' gene_biotype "{gene_biotype}";'
    if transcript_biotype:
        out += f' transcript_biotype "{transcript_biotype}";'
    return out


def write_gtf(
    genes: list[ReconciledGene], path: str, source: str = "HelixForge"
) -> str:
    """Write a GTF (1-based inclusive) for legacy tools.

    Emits ``transcript`` / ``exon`` / ``CDS`` records with ``gene_id`` +
    ``transcript_id`` (+ Ensembl ``gene_biotype`` / ``transcript_biotype`` when
    assigned) attributes. CDS frame is the model phase. Coordinate conversion
    (internal → GTF) happens here, at the I/O boundary.
    """
    with atomic_write(path) as fh:
        for gene in genes:
            strand = gene.strand
            g_biotype = getattr(gene, "biotype", None)
            for tx in gene.transcripts:
                attrs = _gtf_attrs(
                    gene.gene_id,
                    tx.transcript_id,
                    g_biotype,
                    getattr(tx, "biotype", None) or g_biotype,
                )
                t_start, t_end = internal_to_gff3(tx.start, tx.end)
                fh.write(
                    f"{tx.seqid}\t{source}\ttranscript\t{t_start}\t{t_end}\t.\t"
                    f"{strand}\t.\t{attrs}\n"
                )
                for ex in tx.exons:
                    e_start, e_end = internal_to_gff3(ex.start, ex.end)
                    fh.write(
                        f"{tx.seqid}\t{source}\texon\t{e_start}\t{e_end}\t.\t"
                        f"{strand}\t.\t{attrs}\n"
                    )
                for c in tx.cds or []:
                    c_start, c_end = internal_to_gff3(c.start, c.end)
                    fh.write(
                        f"{tx.seqid}\t{source}\tCDS\t{c_start}\t{c_end}\t.\t"
                        f"{strand}\t{c.phase}\t{attrs}\n"
                    )
    return str(path)


# ---------------------------------------------------------------------------
# AGAT-clean GFF3
# ---------------------------------------------------------------------------


def write_agat_clean_gff3(
    genes: list[ReconciledGene],
    path: str,
    source: str = "HelixForge",
    functional: dict[str, Any] | None = None,
) -> str:
    """Write a coordinate-sorted GFF3 that ``agat_sp_statistics`` accepts.

    Wraps :class:`io.gff.GFF3Writer` (consistent ID/Parent, valid phases, the
    full gene→mRNA→exon/CDS hierarchy) after sorting genes by ``(seqid, start)``
    so AGAT/gffcompare see records in genomic order. ``functional``
    optionally adds ``Ontology_term``/``Dbxref`` per transcript; the GFF3 stays
    AGAT-clean because those are reserved GFF3 attributes. Returns the path.
    """
    ordered = sorted(genes, key=lambda g: (g.seqid, g.start, g.end))
    GFF3Writer(str(path)).write_genes(ordered, source=source, functional=functional)
    return str(path)


# ---------------------------------------------------------------------------
# Per-gene JSON (Phase-10 viz / dashboards)
# ---------------------------------------------------------------------------


def _exon_json(ex: Any) -> dict[str, int]:
    return {"start": ex.start, "end": ex.end}


def _cds_json(seg: Any) -> dict[str, int]:
    return {"start": seg.start, "end": seg.end, "phase": seg.phase}


def _transcript_json(tx: Any, primary_id: str) -> dict[str, Any]:
    return {
        "transcript_id": tx.transcript_id,
        "source": tx.source,
        "seqid": tx.seqid,
        "start": tx.start,
        "end": tx.end,
        "strand": tx.strand,
        "is_primary": tx.transcript_id == primary_id,
        "biotype": getattr(tx, "biotype", None),
        "exons": [_exon_json(e) for e in tx.exons],
        "cds": [_cds_json(c) for c in tx.cds] if tx.cds else None,
        "cds_partial": tx.cds_partial,
        "protein_id": tx.protein_id,
        "blast_score": tx.blast_score,
        "tpm": tx.tpm,
        "junction_support": tx.junction_support_fraction,
        "helixer_support": tx.confidence,
        "combined_score": tx.combined_score,
        "has_homology": tx.has_homology,
    }


def _gene_json(gene: ReconciledGene) -> dict[str, Any]:
    return {
        "gene_id": gene.gene_id,
        "seqid": gene.seqid,
        "start": gene.start,
        "end": gene.end,
        "strand": gene.strand,
        "tier": gene.tier,
        "origin": gene.origin,
        "biotype": getattr(gene, "biotype", None),
        "merged_from": list(gene.merged_from),
        "classification": {
            "status": gene.classification.status,
            "max_tpm": gene.classification.max_tpm,
            "mean_coverage": gene.classification.mean_coverage,
            "evidence_source": gene.classification.evidence_source,
            "num_samples_expressed": gene.classification.num_samples_expressed,
        },
        "primary_transcript_id": gene.primary_transcript_id,
        "num_isoforms": len(gene.transcripts),
        "transcripts": [
            _transcript_json(t, gene.primary_transcript_id) for t in gene.transcripts
        ],
        "as_events": [
            {
                "kind": e.kind,
                "seqid": e.seqid,
                "start": e.start,
                "end": e.end,
                "strand": e.strand,
                "support_read_count": e.support_read_count,
            }
            for e in gene.as_events
        ],
        "flags": [f.name for f in gene.flags],
    }


def build_gene_records(genes: list[ReconciledGene]) -> list[dict[str, Any]]:
    """Return the list of per-gene plain-dict records (also used by the viz)."""
    return [_gene_json(g) for g in genes]


def write_per_gene_json(genes: list[ReconciledGene], path: str) -> str:
    """Write per-gene JSON (structure, isoforms, AS events, tier, flags, scores).

    Coordinates are the internal 0-based half-open convention; the viz layer
    converts at its own boundary. Round-trips via ``json``.
    """
    records = build_gene_records(genes)
    with atomic_write(path) as fh:
        json.dump(records, fh, indent=2)
    return str(path)
