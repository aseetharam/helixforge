"""GFF3 <-> GTF format conversion."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from helixforge.utils.logging import get_logger

_log = get_logger(__name__)


def _infer_format(path: str | Path) -> str:
    """Return ``'gff3'`` or ``'gtf'`` from the file extension, or raise."""
    suffix = Path(path).suffix.lower()
    if suffix in (".gff3", ".gff"):
        return "gff3"
    if suffix == ".gtf":
        return "gtf"
    raise ValueError(
        f"cannot infer format from extension {suffix!r} "
        f"(expected .gff3, .gff, or .gtf): {path}"
    )


def gff3_to_gtf(
    input_path: str | Path,
    output_path: str | Path,
    source: str = "HelixForge",
) -> int:
    """Convert a GFF3 to GTF. Returns the number of genes written."""
    from helixforge.io.gff import GFF3Parser
    from helixforge.utils.atomic import atomic_write
    from helixforge.utils.regions import internal_to_gff3

    parser = GFF3Parser(str(input_path))
    genes = parser.parse_genes_generic()

    with atomic_write(str(output_path)) as fh:
        for gene in genes:
            strand = gene["strand"]
            for tx in gene["transcripts"]:
                tid = tx["transcript_id"]
                gid = gene["gene_id"]
                attrs = f'gene_id "{gid}"; transcript_id "{tid}";'

                exons = tx["exons"]
                if not exons:
                    continue
                tx_start = min(e.start for e in exons)
                tx_end = max(e.end for e in exons)
                t_start, t_end = internal_to_gff3(tx_start, tx_end)
                fh.write(
                    f"{gene['seqid']}\t{source}\ttranscript\t{t_start}\t{t_end}"
                    f"\t.\t{strand}\t.\t{attrs}\n"
                )
                for ex in exons:
                    e_start, e_end = internal_to_gff3(ex.start, ex.end)
                    fh.write(
                        f"{gene['seqid']}\t{source}\texon\t{e_start}\t{e_end}"
                        f"\t.\t{strand}\t.\t{attrs}\n"
                    )
                for c in tx["cds"] or []:
                    c_start, c_end = internal_to_gff3(c.start, c.end)
                    fh.write(
                        f"{gene['seqid']}\t{source}\tCDS\t{c_start}\t{c_end}"
                        f"\t.\t{strand}\t{c.phase}\t{attrs}\n"
                    )

    _log.info("gff3_to_gtf: %d genes → %s", len(genes), output_path)
    return len(genes)


def gtf_to_gff3(
    input_path: str | Path,
    output_path: str | Path,
    source: str = "HelixForge",
) -> int:
    """Convert a GTF to GFF3. Returns the number of genes written.

    Parses the GTF by extracting ``gene_id`` and ``transcript_id`` from the
    attribute column, groups features into a gene->transcript hierarchy, and
    writes a valid GFF3 with proper ``ID``/``Parent`` attributes. Coordinates
    are converted at the I/O boundary (GTF is 1-based inclusive; internal is
    0-based half-open).
    """
    from helixforge.utils.atomic import atomic_write
    from helixforge.utils.regions import gff3_to_internal, internal_to_gff3

    genes: dict[str, dict[str, Any]] = {}
    tx_order: dict[str, list[str]] = {}

    with open(input_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, src, ftype, start_s, end_s, _, strand, phase_s, attrs_s = cols[:9]
            if ftype not in ("transcript", "exon", "CDS", "mRNA"):
                continue

            gene_id = _parse_gtf_attr(attrs_s, "gene_id")
            transcript_id = _parse_gtf_attr(attrs_s, "transcript_id")
            if not gene_id or not transcript_id:
                continue

            start, end = gff3_to_internal(int(start_s), int(end_s))

            if gene_id not in genes:
                genes[gene_id] = {
                    "seqid": seqid,
                    "strand": strand,
                    "transcripts": {},
                    "start": start,
                    "end": end,
                }
                tx_order[gene_id] = []
            g = genes[gene_id]
            g["start"] = min(g["start"], start)
            g["end"] = max(g["end"], end)

            if transcript_id not in g["transcripts"]:
                g["transcripts"][transcript_id] = {"exons": [], "cds": []}
                tx_order[gene_id].append(transcript_id)
            tx = g["transcripts"][transcript_id]

            if ftype == "exon":
                tx["exons"].append((start, end))
            elif ftype == "CDS":
                phase = int(phase_s) if phase_s in ("0", "1", "2") else 0
                tx["cds"].append((start, end, phase))

    with atomic_write(str(output_path)) as fh:
        fh.write("##gff-version 3\n")
        for gene_id, g in genes.items():
            g_start, g_end = internal_to_gff3(g["start"], g["end"])
            fh.write(
                f"{g['seqid']}\t{source}\tgene\t{g_start}\t{g_end}\t.\t"
                f"{g['strand']}\t.\tID={gene_id}\n"
            )
            for tid in tx_order[gene_id]:
                tx = g["transcripts"][tid]
                exons = sorted(tx["exons"])
                if not exons:
                    continue
                tx_start, tx_end = internal_to_gff3(exons[0][0], exons[-1][1])
                fh.write(
                    f"{g['seqid']}\t{source}\tmRNA\t{tx_start}\t{tx_end}\t.\t"
                    f"{g['strand']}\t.\tID={tid};Parent={gene_id}\n"
                )
                for i, (es, ee) in enumerate(exons, 1):
                    e_start, e_end = internal_to_gff3(es, ee)
                    fh.write(
                        f"{g['seqid']}\t{source}\texon\t{e_start}\t{e_end}\t.\t"
                        f"{g['strand']}\t.\tID={tid}.exon{i};Parent={tid}\n"
                    )
                for i, (cs, ce, phase) in enumerate(sorted(tx["cds"]), 1):
                    c_start, c_end = internal_to_gff3(cs, ce)
                    fh.write(
                        f"{g['seqid']}\t{source}\tCDS\t{c_start}\t{c_end}\t.\t"
                        f"{g['strand']}\t{phase}\tID={tid}.CDS{i};Parent={tid}\n"
                    )

    _log.info("gtf_to_gff3: %d genes → %s", len(genes), output_path)
    return len(genes)


def _parse_gtf_attr(attrs: str, key: str) -> str | None:
    """Extract a GTF attribute value by key (e.g. ``gene_id "g1"`` → ``g1``)."""
    for part in attrs.split(";"):
        part = part.strip()
        if part.startswith(key):
            rest = part[len(key) :].strip()
            return rest.strip('"').strip("'")
    return None
