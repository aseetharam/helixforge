"""Mikado GTF emitters."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from helixforge.io.gff import GFF3Parser
from helixforge.utils.regions import internal_to_gff3

if TYPE_CHECKING:
    from helixforge.reconcile.models import StringTieTranscript


def _attrs(
    gene_id: str, transcript_id: str, extra: dict[str, str] | None = None
) -> str:
    parts = [f'gene_id "{gene_id}";', f' transcript_id "{transcript_id}";']
    if extra:
        for key, value in extra.items():
            parts.append(f' {key} "{value}";')
    return "".join(parts)


def _gtf_line(
    seqid: str,
    source: str,
    feature: str,
    start: int,
    end: int,
    strand: str,
    frame: str,
    attributes: str,
) -> str:
    """One GTF row; ``start``/``end`` are internal 0-based half-open and are
    converted to 1-based inclusive here."""
    g_start, g_end = internal_to_gff3(start, end)
    return (
        f"{seqid}\t{source}\t{feature}\t{g_start}\t{g_end}\t.\t"
        f"{strand}\t{frame}\t{attributes}\n"
    )


def helixer_gff3_to_gtf(gff3_path: str | Path, out_path: str | Path) -> Path:
    """Convert a Helixer GFF3 to GTF for ``mikado prepare`` (one transcript/gene).

    Each gene becomes a ``transcript`` + ``exon`` (+ optional ``CDS``) block with
    ``gene_id``/``transcript_id`` attributes. Returns the output ``Path``.
    """
    loci = GFF3Parser(gff3_path).parse_helixer_genes()
    out_path = Path(out_path)
    with out_path.open("w") as fh:
        for locus in loci:
            gene_id = locus.gene_id
            tid = f"{gene_id}.1"
            attrs = _attrs(gene_id, tid)
            fh.write(
                _gtf_line(
                    locus.seqid,
                    "Helixer",
                    "transcript",
                    locus.start,
                    locus.end,
                    locus.strand,
                    ".",
                    attrs,
                )
            )
            for ex in locus.exons:
                fh.write(
                    _gtf_line(
                        locus.seqid,
                        "Helixer",
                        "exon",
                        ex.start,
                        ex.end,
                        locus.strand,
                        ".",
                        attrs,
                    )
                )
            if locus.cds:
                for seg in locus.cds:
                    fh.write(
                        _gtf_line(
                            locus.seqid,
                            "Helixer",
                            "CDS",
                            seg.start,
                            seg.end,
                            locus.strand,
                            str(seg.phase),
                            attrs,
                        )
                    )
    return out_path


def verify_unique_transcript_ids(transcripts: list[StringTieTranscript]) -> None:
    """Raise ``ValueError`` if any transcript_id repeats (Mikado needs them unique)."""
    seen: set[str] = set()
    for t in transcripts:
        if t.transcript_id in seen:
            raise ValueError(f"duplicate transcript_id: {t.transcript_id!r}")
        seen.add(t.transcript_id)


def stringtie_to_labelled_gtf(
    transcripts: list[StringTieTranscript],
    out_path: str | Path,
    label: str,
) -> Path:
    """Write ``StringTieTranscript`` objects to GTF, ``label`` in the source column.

    transcript_ids are verified unique first. TPM/coverage are preserved as
    attributes. Returns the output ``Path``.
    """
    verify_unique_transcript_ids(transcripts)
    out_path = Path(out_path)
    with out_path.open("w") as fh:
        for t in transcripts:
            extra: dict[str, str] = {}
            if t.tpm is not None:
                extra["TPM"] = f"{t.tpm}"
            if t.coverage is not None:
                extra["cov"] = f"{t.coverage}"
            attrs = _attrs(t.gene_id, t.transcript_id, extra)
            fh.write(
                _gtf_line(
                    t.seqid, label, "transcript", t.start, t.end, t.strand, ".", attrs
                )
            )
            for ex in t.exons:
                fh.write(
                    _gtf_line(
                        t.seqid, label, "exon", ex.start, ex.end, t.strand, ".", attrs
                    )
                )
    return out_path
