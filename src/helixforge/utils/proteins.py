"""CDS-to-protein extraction from GFF3 + genome FASTA."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from helixforge.io.fasta import GenomeAccessor
from helixforge.io.gff import GFF3Parser
from helixforge.utils.sequences import reverse_complement, translate

log = logging.getLogger(__name__)


def _cds_sequence(
    genome: GenomeAccessor,
    seqid: str,
    strand: str,
    cds_segments: list[Any],
) -> str:
    """Concatenate CDS segments into a coding sequence in reading order.

    Segments are always stored low-to-high; for minus-strand genes the
    concatenation is reversed and reverse-complemented.
    """
    parts = [genome.get_sequence(seqid, seg.start, seg.end) for seg in cds_segments]
    if strand == "-":
        parts.reverse()
        return reverse_complement("".join(parts))
    return "".join(parts)


def _write_fasta_record(
    fh: Any,
    header: str,
    sequence: str,
    line_width: int = 60,
) -> None:
    fh.write(f">{header}\n")
    for i in range(0, len(sequence), line_width):
        fh.write(sequence[i : i + line_width] + "\n")


def extract_proteins(
    gff3_path: str | Path,
    genome_path: str | Path,
    output_fasta: str | Path,
    longest_only: bool = True,
    transl_table: int = 1,
    min_length: int = 30,
) -> dict[str, int]:
    """Translate CDS regions from gene models to protein FASTA.

    Args:
        gff3_path: GFF3 with gene models (any valid GFF3).
        genome_path: Reference genome FASTA with .fai index.
        output_fasta: Output protein FASTA path.
        longest_only: Keep only the longest isoform per gene.
        transl_table: NCBI genetic code table id.
        min_length: Minimum protein length in amino acids.

    Returns:
        {protein_id: length_aa} dict.
    """
    parser = GFF3Parser(str(gff3_path))
    genes = parser.parse_genes_generic()

    collected: list[tuple[str, str]] = []
    skipped = 0

    with GenomeAccessor(str(genome_path)) as genome:
        for gene in genes:
            gene_id = gene["gene_id"]
            strand = gene["strand"]
            seqid = gene["seqid"]
            candidates: list[tuple[str, str]] = []

            for tx in gene["transcripts"]:
                if tx["cds"] is None:
                    continue
                try:
                    cds_seq = _cds_sequence(genome, seqid, strand, tx["cds"])
                except (KeyError, IndexError, ValueError) as exc:
                    log.warning("skipping %s: %s", tx["transcript_id"], exc)
                    skipped += 1
                    continue

                protein = translate(cds_seq, transl_table=transl_table)
                if len(protein) < min_length:
                    skipped += 1
                    continue

                tid = tx["transcript_id"]
                pid = f"{gene_id}.{tid}"
                candidates.append((pid, protein))

            if not candidates:
                continue

            if longest_only:
                candidates.sort(key=lambda c: len(c[1]), reverse=True)
                candidates = candidates[:1]

            collected.extend(candidates)

    results: dict[str, int] = {}
    with open(output_fasta, "w") as fh:
        for pid, protein in collected:
            _write_fasta_record(fh, f"{pid} length={len(protein)}", protein)
            results[pid] = len(protein)

    log.info(
        "extracted %d proteins from %d genes (%d skipped)",
        len(results),
        len(genes),
        skipped,
    )
    return results
