"""``helixforge utils``: standalone data-prep and post-run utilities."""

from __future__ import annotations

import sys

import click


def _stub(name: str) -> None:
    click.echo(f"Not yet implemented: {name}", err=True)
    raise SystemExit(0)


@click.group()
def utils() -> None:
    """Standalone data-prep and post-run utilities."""


# ---------------------------------------------------------------------------
# fetch-db
# ---------------------------------------------------------------------------


@utils.command(
    "fetch-db",
    epilog="""\b
Examples:
  helixforge utils fetch-db --list
\b
  helixforge utils fetch-db --db uniprot_sprot --cache-dir /scratch/dbs
\b
  helixforge utils fetch-db --db /local/proteins.fa --format-tool diamond --force
""",
)
@click.option("--db", default=None, help="Database name or local FASTA path.")
@click.option(
    "--cache-dir",
    default=None,
    type=click.Path(file_okay=False),
    help="Cache directory [default: ~/.helixforge/databases].",
)
@click.option(
    "--force", is_flag=True, default=False, help="Re-download even if cached."
)
@click.option(
    "--list",
    "list_dbs",
    is_flag=True,
    default=False,
    help="List available databases and exit.",
)
@click.option(
    "--format-tool", default="diamond", help="Formatting tool [default: diamond]."
)
def fetch_db(
    db: str | None,
    cache_dir: str | None,
    force: bool,
    list_dbs: bool,
    format_tool: str,
) -> None:
    """Download, decompress, and Diamond-format reference protein databases."""
    from helixforge.utils.databases import (
        DatabaseManager,
        PREDEFINED_DATABASES,
    )

    if list_dbs:
        click.echo(f"{'Name':<20} {'Description'}")
        click.echo("-" * 60)
        for key, cfg in PREDEFINED_DATABASES.items():
            click.echo(f"{key:<20} {cfg['description']}")
        return

    if db is None:
        raise click.UsageError("--db is required unless --list is given.")

    manager = DatabaseManager(
        cache_dir=cache_dir,
        diamond_bin=format_tool,
    )
    try:
        info = manager.get_database(db, force_download=force)
    except FileNotFoundError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)
    except ValueError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    click.echo(f"Database: {info.name}")
    click.echo(f"  FASTA:     {info.path}")
    if info.formatted_path:
        click.echo(f"  Diamond:   {info.formatted_path}")
    if info.n_sequences is not None:
        click.echo(f"  Sequences: {info.n_sequences:,}")
    if info.checksum:
        click.echo(f"  MD5:       {info.checksum}")


# ---------------------------------------------------------------------------
# extract-proteins
# ---------------------------------------------------------------------------


@utils.command(
    "extract-proteins",
    epilog="""\b
Examples:
  helixforge utils extract-proteins --gff3 genes.gff3 --genome genome.fa --out proteins.fa
\b
  helixforge utils extract-proteins --gff3 genes.gff3 --genome genome.fa \\
      --out proteins.fa --all-isoforms --min-length 50 --transl-table 11
""",
)
@click.option(
    "--gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Input GFF3.",
)
@click.option(
    "--genome",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Reference genome FASTA.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output protein FASTA.",
)
@click.option(
    "--longest-only/--all-isoforms",
    default=True,
    help="Keep only longest isoform per gene.",
)
@click.option(
    "--transl-table", default=1, type=int, help="NCBI translation table [default: 1]."
)
@click.option(
    "--min-length",
    default=30,
    type=int,
    help="Minimum protein length in aa [default: 30].",
)
def extract_proteins_cmd(
    gff3: str,
    genome: str,
    out: str,
    longest_only: bool,
    transl_table: int,
    min_length: int,
) -> None:
    """Translate CDS from GFF3 + genome to protein FASTA."""
    from helixforge.utils.proteins import extract_proteins

    results = extract_proteins(
        gff3_path=gff3,
        genome_path=genome,
        output_fasta=out,
        longest_only=longest_only,
        transl_table=transl_table,
        min_length=min_length,
    )
    n_genes = len({pid.rsplit(".", 1)[0] for pid in results})
    click.echo(f"Extracted {len(results)} proteins from {n_genes} genes to {out}")


# ---------------------------------------------------------------------------
# align
# ---------------------------------------------------------------------------


@utils.command(
    "align",
    epilog="""\b
Examples:
  helixforge utils align --genome genome.fa --proteins uniprot.fa --out miniprot.gff
\b
  helixforge utils align --genome genome.fa --proteins uniprot.fa \\
      --out miniprot.gff --threads 16 --miniprot-bin /opt/bin/miniprot
""",
)
@click.option(
    "--genome",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Reference genome FASTA.",
)
@click.option(
    "--proteins",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Protein FASTA.",
)
@click.option(
    "--out", required=True, type=click.Path(dir_okay=False), help="Output GFF."
)
@click.option("--threads", default=4, type=int, help="Number of threads [default: 4].")
@click.option(
    "--miniprot-bin", default="miniprot", help="miniprot binary [default: miniprot]."
)
def align(
    genome: str,
    proteins: str,
    out: str,
    threads: int,
    miniprot_bin: str,
) -> None:
    """Run miniprot protein-to-genome alignment."""
    from helixforge.prep.protein_align import run_miniprot

    result = run_miniprot(
        genome_fasta=genome,
        proteome_fasta=proteins,
        out_gff=out,
        threads=threads,
        miniprot_bin=miniprot_bin,
        force=True,
    )
    click.echo(f"miniprot alignment written to {result}")


# ---------------------------------------------------------------------------
# qc
# ---------------------------------------------------------------------------


@utils.command(
    "qc",
    epilog="""\b
Examples:
  helixforge utils qc --gff3 helixforge.gff3 --out report.html
\b
  helixforge utils qc --gff3 helixforge.gff3 --genome genome.fa \\
      --helixer-h5 predictions.h5 --out report.json --format json
""",
)
@click.option(
    "--gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Reconciled GFF3.",
)
@click.option(
    "--genome",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Reference genome FASTA.",
)
@click.option(
    "--helixer-h5",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Helixer HDF5.",
)
@click.option(
    "--out", required=True, type=click.Path(dir_okay=False), help="Output report path."
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["html", "json", "tsv"]),
    default="html",
    help="Output format [default: html].",
)
def qc(
    gff3: str,
    genome: str | None,
    helixer_h5: str | None,
    out: str,
    fmt: str,
) -> None:
    """Genome-wide QC report (HTML/JSON/TSV)."""
    from helixforge.stats.genome_report import generate_qc_report

    stats = generate_qc_report(
        gff3_path=gff3,
        output_path=out,
        genome_path=genome,
        helixer_h5_path=helixer_h5,
        fmt=fmt,
    )
    click.echo(
        f"QC report written to {out} ({stats['total_genes']} genes, format={fmt})"
    )


# ---------------------------------------------------------------------------
# convert
# ---------------------------------------------------------------------------


@utils.command(
    "convert",
    epilog="""\b
Examples:
  helixforge utils convert --input genes.gff3 --out genes.gtf
\b
  helixforge utils convert --input genes.gtf --out genes.gff3 \\
      --from-format gtf --to-format gff3
""",
)
@click.option(
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Input file.",
)
@click.option(
    "--out", required=True, type=click.Path(dir_okay=False), help="Output file."
)
@click.option(
    "--from-format",
    default=None,
    type=str,
    help="Source format: gff3 | gtf (inferred from extension).",
)
@click.option(
    "--to-format",
    default=None,
    type=str,
    help="Target format: gff3 | gtf (inferred from extension).",
)
def convert(
    input_path: str,
    out: str,
    from_format: str | None,
    to_format: str | None,
) -> None:
    """GFF3 <-> GTF format conversion."""
    from helixforge.utils.convert import (
        _infer_format,
        gff3_to_gtf,
        gtf_to_gff3,
    )

    src_fmt = from_format or _infer_format(input_path)
    dst_fmt = to_format or _infer_format(out)

    if src_fmt == dst_fmt:
        raise click.UsageError(
            f"source and target formats are the same ({src_fmt}); "
            "use --from-format / --to-format to override"
        )

    if src_fmt == "gff3" and dst_fmt == "gtf":
        n = gff3_to_gtf(input_path, out)
    elif src_fmt == "gtf" and dst_fmt == "gff3":
        n = gtf_to_gff3(input_path, out)
    else:
        raise click.UsageError(f"unsupported conversion: {src_fmt} → {dst_fmt}")
    click.echo(f"Converted {n} genes: {input_path} ({src_fmt}) → {out} ({dst_fmt})")


# ---------------------------------------------------------------------------
# filter
# ---------------------------------------------------------------------------


@utils.command(
    "filter",
    epilog="""\b
Examples:
  helixforge utils filter --gff3 helixforge.gff3 --out tier1.gff3 --max-tier 1
\b
  helixforge utils filter --gff3 helixforge.gff3 --out pub.gff3 \\
      --preset publication_ready --exclude-biotype transposable_element
""",
)
@click.option(
    "--gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Input GFF3.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(dir_okay=False),
    help="Output filtered GFF3.",
)
@click.option(
    "--preset",
    type=click.Choice(["high_confidence", "publication_ready", "custom"]),
    default="custom",
    help="Filter preset [default: custom].",
)
@click.option(
    "--min-confidence", default=None, type=float, help="Minimum confidence score."
)
@click.option("--min-tier", default=None, type=int, help="Minimum tier (inclusive).")
@click.option("--max-tier", default=None, type=int, help="Maximum tier (inclusive).")
@click.option(
    "--exclude-biotype",
    multiple=True,
    type=str,
    help="Biotype(s) to exclude, repeatable.",
)
def filter_cmd(
    gff3: str,
    out: str,
    preset: str,
    min_confidence: float | None,
    min_tier: int | None,
    max_tier: int | None,
    exclude_biotype: tuple[str, ...],
) -> None:
    """Tiered output filtering by confidence/evidence/biotype."""
    from helixforge.utils.filters import (
        FilterCriteria,
        GeneFilter,
        parse_gff3_genes,
        write_filtered_gff3,
    )

    genes = parse_gff3_genes(gff3)

    if preset == "high_confidence":
        filt = GeneFilter.high_confidence()
    elif preset == "publication_ready":
        filt = GeneFilter.publication_ready()
    else:
        criteria = FilterCriteria(
            min_tier=min_tier,
            max_tier=max_tier,
            exclude_biotypes=list(exclude_biotype),
            min_confidence=min_confidence,
        )
        filt = GeneFilter(criteria)

    kept, excluded = filt.apply(genes)
    write_filtered_gff3(kept, out, header_from=gff3)
    click.echo(
        f"Kept {len(kept)} / {len(genes)} genes → {out} (excluded {len(excluded)})"
    )


# ---------------------------------------------------------------------------
# summarize
# ---------------------------------------------------------------------------


@utils.command(
    "summarize",
    epilog="""\b
Examples:
  helixforge utils summarize --gff3 helixforge.gff3
\b
  helixforge utils summarize --gff3 helixforge.gff3 --out stats.json --format json
""",
)
@click.option(
    "--gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Input GFF3.",
)
@click.option(
    "--out",
    default=None,
    type=click.Path(dir_okay=False),
    help="Output file (default: stdout).",
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["tsv", "json", "markdown"]),
    default="tsv",
    help="Output format [default: tsv].",
)
def summarize(
    gff3: str,
    out: str | None,
    fmt: str,
) -> None:
    """Annotation statistics table."""
    import json as _json

    from helixforge.io.gff import GFF3Parser

    parser = GFF3Parser(gff3)
    genes = parser.parse_genes_generic()

    gene_count = len(genes)
    transcript_count = sum(len(g["transcripts"]) for g in genes)
    multi_isoform = sum(1 for g in genes if len(g["transcripts"]) > 1)
    coding = sum(1 for g in genes if any(t.get("cds") for t in g["transcripts"]))
    mono_exon = sum(
        1 for g in genes if all(len(t.get("exons", [])) == 1 for t in g["transcripts"])
    )

    stats = {
        "gene_count": gene_count,
        "transcript_count": transcript_count,
        "multi_isoform_genes": multi_isoform,
        "coding_genes": coding,
        "mono_exon_genes": mono_exon,
    }

    if fmt == "json":
        text = _json.dumps(stats, indent=2)
    elif fmt == "markdown":
        lines = ["| Metric | Value |", "| --- | --- |"]
        for k, v in stats.items():
            lines.append(f"| {k} | {v} |")
        text = "\n".join(lines)
    else:
        lines = ["metric\tvalue"]
        for k, v in stats.items():
            lines.append(f"{k}\t{v}")
        text = "\n".join(lines)

    if out:
        from pathlib import Path

        Path(out).write_text(text + "\n")
        click.echo(f"Summary written to {out}")
    else:
        click.echo(text)
