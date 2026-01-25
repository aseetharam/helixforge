"""Command-line interface for HelixForge.

This module provides the main entry point for the helixforge CLI tool.
It uses Click to define commands and subcommands for various operations.

Commands:
    refine: Main refinement pipeline
    add-evidence: Add RNA-seq evidence to predictions
    qc: Generate quality control reports
    validate: Homology-based validation
    confidence: Calculate confidence scores for gene predictions
    viz: Generate visualizations

Example:
    $ helixforge --help
    $ helixforge refine --predictions helixer.h5 --genome genome.fa
    $ helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa -o scores.tsv
    $ helixforge qc --gff refined.gff3 --output report.html
"""

from pathlib import Path
from typing import Optional

import click
from rich.console import Console

# Initialize rich console for pretty output
console = Console()


@click.group()
@click.version_option(prog_name="helixforge")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
@click.option("-q", "--quiet", is_flag=True, help="Suppress non-error output.")
@click.pass_context
def main(ctx: click.Context, verbose: bool, quiet: bool) -> None:
    """HelixForge: Refine Helixer gene predictions into publication-quality annotations.

    HelixForge integrates RNA-seq evidence, homology validation, and confidence
    scoring to transform raw Helixer predictions into high-quality genome
    annotations.
    """
    # Ensure context object exists
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet


# =============================================================================
# refine command
# =============================================================================


@main.command()
@click.option(
    "-p",
    "--predictions",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Helixer HDF5 predictions file.",
)
@click.option(
    "-g",
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Reference genome FASTA file.",
)
@click.option(
    "-b",
    "--bam",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="RNA-seq BAM file(s) for evidence. Can be specified multiple times.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output GFF3 file path.",
)
@click.option(
    "--min-confidence",
    type=float,
    default=0.5,
    show_default=True,
    help="Minimum confidence score to retain a gene.",
)
@click.option(
    "-j",
    "--threads",
    type=int,
    default=1,
    show_default=True,
    help="Number of parallel threads.",
)
@click.option(
    "--chunk-size",
    type=int,
    default=1000000,
    show_default=True,
    help="Genomic chunk size for parallel processing.",
)
@click.pass_context
def refine(
    ctx: click.Context,
    predictions: Path,
    genome: Path,
    bam: tuple[Path, ...],
    output: Path,
    min_confidence: float,
    threads: int,
    chunk_size: int,
) -> None:
    """Refine Helixer predictions using evidence and confidence scoring.

    This is the main HelixForge pipeline that takes raw Helixer predictions
    and produces refined, high-quality gene annotations.

    \b
    Steps performed:
    1. Load Helixer HDF5 predictions
    2. Extract RNA-seq evidence (if BAM files provided)
    3. Calculate confidence scores
    4. Refine gene boundaries
    5. Detect and resolve merge/split errors
    6. Apply quality filters
    7. Write refined GFF3

    Example:
        $ helixforge refine -p helixer.h5 -g genome.fa -b rnaseq.bam -o refined.gff3
    """
    # TODO: Implement refinement pipeline
    console.print("[yellow]refine command not yet implemented[/yellow]")
    console.print(f"  Predictions: {predictions}")
    console.print(f"  Genome: {genome}")
    console.print(f"  BAM files: {bam if bam else 'None'}")
    console.print(f"  Output: {output}")
    console.print(f"  Min confidence: {min_confidence}")
    console.print(f"  Threads: {threads}")
    raise SystemExit(1)


# =============================================================================
# add-evidence command
# =============================================================================


@main.command("add-evidence")
@click.option(
    "-g",
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input GFF3 file with gene predictions.",
)
@click.option(
    "-b",
    "--bam",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    multiple=True,
    help="RNA-seq BAM file(s). Can be specified multiple times.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output GFF3 file with evidence annotations.",
)
@click.option(
    "--min-reads",
    type=int,
    default=3,
    show_default=True,
    help="Minimum reads to support a splice junction.",
)
@click.pass_context
def add_evidence(
    ctx: click.Context,
    gff: Path,
    bam: tuple[Path, ...],
    output: Path,
    min_reads: int,
) -> None:
    """Add RNA-seq evidence to existing gene predictions.

    Annotates gene models with supporting evidence from RNA-seq data,
    including splice junction support and coverage metrics.

    Example:
        $ helixforge add-evidence -g predictions.gff3 -b rnaseq.bam -o annotated.gff3
    """
    # TODO: Implement evidence addition
    console.print("[yellow]add-evidence command not yet implemented[/yellow]")
    raise SystemExit(1)


# =============================================================================
# qc command
# =============================================================================


@main.command()
@click.option(
    "-g",
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input GFF3 file to analyze.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output report file (HTML or PDF based on extension).",
)
@click.option(
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome for sequence-based QC metrics.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["html", "pdf", "json"]),
    default="html",
    show_default=True,
    help="Output format for the report.",
)
@click.pass_context
def qc(
    ctx: click.Context,
    gff: Path,
    output: Path,
    genome: Optional[Path],
    output_format: str,
) -> None:
    """Generate quality control report for gene predictions.

    Analyzes gene predictions and produces a comprehensive QC report
    including statistics, flag summaries, and visualizations.

    Example:
        $ helixforge qc -g refined.gff3 -o report.html
    """
    # TODO: Implement QC report generation
    console.print("[yellow]qc command not yet implemented[/yellow]")
    raise SystemExit(1)


# =============================================================================
# validate command
# =============================================================================


@main.command()
@click.option(
    "-g",
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input GFF3 file with gene predictions.",
)
@click.option(
    "-p",
    "--proteins",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Protein database for homology search (FASTA format).",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output GFF3 with validation annotations.",
)
@click.option(
    "--tool",
    type=click.Choice(["diamond", "blast"]),
    default="diamond",
    show_default=True,
    help="Homology search tool to use.",
)
@click.option(
    "-e",
    "--evalue",
    type=float,
    default=1e-5,
    show_default=True,
    help="E-value threshold for homology matches.",
)
@click.pass_context
def validate(
    ctx: click.Context,
    gff: Path,
    proteins: Path,
    output: Path,
    tool: str,
    evalue: float,
) -> None:
    """Validate gene predictions using protein homology.

    Searches predicted proteins against a reference database and
    annotates genes with homology-based validation scores.

    Example:
        $ helixforge validate -g refined.gff3 -p uniprot.fa -o validated.gff3
    """
    # TODO: Implement homology validation
    console.print("[yellow]validate command not yet implemented[/yellow]")
    raise SystemExit(1)


# =============================================================================
# confidence command
# =============================================================================


@main.command()
@click.option(
    "-p",
    "--predictions",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Helixer HDF5 predictions file.",
)
@click.option(
    "-g",
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="GFF3 file with gene predictions.",
)
@click.option(
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Reference genome FASTA file.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output TSV file for confidence scores.",
)
@click.option(
    "--bed",
    type=click.Path(path_type=Path),
    help="Output BED file for genome browser visualization.",
)
@click.option(
    "--low-conf-bed",
    type=click.Path(path_type=Path),
    help="Output BED file for low-confidence regions within genes.",
)
@click.option(
    "--plot-dir",
    type=click.Path(path_type=Path),
    help="Directory for per-gene confidence plots.",
)
@click.option(
    "--distribution-plot",
    type=click.Path(path_type=Path),
    help="Output path for genome-wide distribution plot.",
)
@click.option(
    "--threshold",
    type=float,
    default=0.7,
    show_default=True,
    help="Threshold for low-confidence regions.",
)
@click.option(
    "-j",
    "--threads",
    type=int,
    default=1,
    show_default=True,
    help="Number of parallel threads.",
)
@click.option(
    "--format",
    "plot_format",
    type=click.Choice(["html", "png", "pdf"]),
    default="html",
    show_default=True,
    help="Format for visualization output.",
)
@click.pass_context
def confidence(
    ctx: click.Context,
    predictions: Path,
    gff: Path,
    genome: Path,
    output: Path,
    bed: Optional[Path],
    low_conf_bed: Optional[Path],
    plot_dir: Optional[Path],
    distribution_plot: Optional[Path],
    threshold: float,
    threads: int,
    plot_format: str,
) -> None:
    """Calculate confidence scores for gene predictions.

    Computes multi-factor confidence metrics for each gene based on
    Helixer HDF5 predictions, including:

    \b
    - Mean/min/median class probabilities
    - Shannon entropy (prediction uncertainty)
    - Exon/intron boundary sharpness
    - CDS coding consistency
    - Per-exon confidence scores

    Genes are classified as high (>=0.85), medium (>=0.70), or low (<0.70)
    confidence, with specific flags for problematic regions.

    Example:
        $ helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa -o scores.tsv
        $ helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa -o scores.tsv --bed scores.bed --distribution-plot dist.html
    """
    from helixforge.core.confidence import (
        ConfidenceCalculator,
        ConfidenceWriter,
    )
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser
    from helixforge.io.hdf5 import HelixerHDF5Reader

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    # Get FAI path
    fai_path = genome.with_suffix(genome.suffix + ".fai")
    if not fai_path.exists():
        console.print(f"[red]Error:[/red] FAI index not found: {fai_path}")
        console.print("Run 'samtools faidx' to create the index.")
        raise SystemExit(1)

    if not quiet:
        console.print(f"[blue]Loading predictions from:[/blue] {predictions}")
        console.print(f"[blue]Loading genes from:[/blue] {gff}")
        console.print(f"[blue]Reference genome:[/blue] {genome}")

    try:
        # Load data
        with HelixerHDF5Reader(predictions, fai_path) as reader:
            with GenomeAccessor(genome) as genome_accessor:
                parser = GFF3Parser(gff)
                genes = list(parser.iter_genes())

                if not quiet:
                    console.print(f"[green]Loaded {len(genes)} genes[/green]")

                # Create calculator
                calc = ConfidenceCalculator(
                    reader,
                    genome_accessor,
                    low_conf_threshold=threshold,
                )

                # Score genes
                if not quiet:
                    console.print(
                        f"[blue]Scoring genes with {threads} thread(s)...[/blue]"
                    )

                scores = list(
                    calc.score_genes_parallel(genes, n_workers=threads)
                )

                # Write TSV output
                ConfidenceWriter.to_tsv(scores, output)
                if not quiet:
                    console.print(f"[green]Wrote TSV to:[/green] {output}")

                # Write BED if requested
                if bed:
                    ConfidenceWriter.to_bed(scores, bed)
                    if not quiet:
                        console.print(f"[green]Wrote BED to:[/green] {bed}")

                # Write low-confidence regions BED if requested
                if low_conf_bed:
                    ConfidenceWriter.low_confidence_regions_bed(scores, low_conf_bed)
                    if not quiet:
                        console.print(
                            f"[green]Wrote low-conf regions to:[/green] {low_conf_bed}"
                        )

                # Generate distribution plot if requested
                if distribution_plot:
                    from helixforge.viz.genome import plot_confidence_distribution

                    plot_confidence_distribution(
                        scores,
                        output_path=distribution_plot,
                        format=plot_format,
                    )
                    if not quiet:
                        console.print(
                            f"[green]Wrote distribution plot to:[/green] {distribution_plot}"
                        )

                # Generate per-gene plots if requested
                if plot_dir:
                    from helixforge.viz.locus import plot_gene_confidence_batch

                    plot_dir.mkdir(parents=True, exist_ok=True)
                    paths = plot_gene_confidence_batch(
                        genes,
                        scores,
                        plot_dir,
                        format=plot_format,
                        calc=calc,
                    )
                    if not quiet:
                        console.print(
                            f"[green]Generated {len(paths)} gene plots in:[/green] {plot_dir}"
                        )

                # Print summary
                if not quiet:
                    high_count = sum(1 for s in scores if s.confidence_class == "high")
                    medium_count = sum(
                        1 for s in scores if s.confidence_class == "medium"
                    )
                    low_count = sum(1 for s in scores if s.confidence_class == "low")

                    console.print("\n[bold]Summary:[/bold]")
                    console.print(f"  [green]High confidence:[/green] {high_count}")
                    console.print(f"  [yellow]Medium confidence:[/yellow] {medium_count}")
                    console.print(f"  [red]Low confidence:[/red] {low_count}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise SystemExit(1)


# =============================================================================
# splice command
# =============================================================================


@main.command()
@click.option(
    "--helixer-gff",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Helixer GFF3 output file.",
)
@click.option(
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Genome FASTA file (indexed).",
)
@click.option(
    "--rnaseq-bam",
    type=click.Path(exists=True, path_type=Path),
    help="RNA-seq BAM file (sorted, indexed).",
)
@click.option(
    "--junctions-bed",
    type=click.Path(exists=True, path_type=Path),
    help="Pre-extracted junctions BED (alternative to BAM).",
)
@click.option(
    "-o",
    "--output-gff",
    type=click.Path(path_type=Path),
    required=True,
    help="Output refined GFF3.",
)
@click.option(
    "-r",
    "--report",
    type=click.Path(path_type=Path),
    help="Output splice refinement report TSV.",
)
@click.option(
    "--corrections-detail",
    type=click.Path(path_type=Path),
    help="Output detailed corrections TSV.",
)
@click.option(
    "--unsupported-bed",
    type=click.Path(path_type=Path),
    help="Output BED of unsupported introns.",
)
@click.option(
    "--max-shift",
    default=15,
    type=int,
    show_default=True,
    help="Maximum splice site correction distance.",
)
@click.option(
    "--min-reads",
    default=3,
    type=int,
    show_default=True,
    help="Minimum junction read support.",
)
@click.option(
    "--adjust-boundaries",
    is_flag=True,
    help="Also adjust start/stop codons.",
)
@click.option(
    "-j",
    "--workers",
    default=1,
    type=int,
    show_default=True,
    help="Number of parallel workers.",
)
@click.option(
    "-v",
    "--verbose",
    "verbose_flag",
    is_flag=True,
    help="Verbose output.",
)
@click.pass_context
def splice(
    ctx: click.Context,
    helixer_gff: Path,
    genome: Path,
    rnaseq_bam: Optional[Path],
    junctions_bed: Optional[Path],
    output_gff: Path,
    report: Optional[Path],
    corrections_detail: Optional[Path],
    unsupported_bed: Optional[Path],
    max_shift: int,
    min_reads: int,
    adjust_boundaries: bool,
    workers: int,
    verbose_flag: bool,
) -> None:
    """Refine splice sites using RNA-seq evidence.

    This command refines splice sites in Helixer gene predictions using
    RNA-seq junction evidence and position weight matrix (PWM) scoring.

    \b
    Features:
    - Correct splice sites to match empirical junctions
    - Score candidate sites using PWMs for canonical splice motifs
    - Handle canonical (GT-AG) and non-canonical splice sites
    - Flag genes with no RNA-seq support or conflicting evidence
    - Optionally adjust start/stop codons

    Example:
        $ helixforge splice --helixer-gff helixer.gff3 --genome genome.fa \\
            --rnaseq-bam rnaseq.bam -o refined.gff3 -r splice_report.tsv
    """
    from helixforge.core.boundaries import BoundaryAdjuster
    from helixforge.core.splice import (
        SpliceRefiner,
        SpliceReportWriter,
    )
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser, GFF3Writer

    verbose = ctx.obj.get("verbose", False) or verbose_flag
    quiet = ctx.obj.get("quiet", False)

    # Validate inputs
    if not rnaseq_bam and not junctions_bed:
        console.print(
            "[red]Error:[/red] Provide either --rnaseq-bam or --junctions-bed"
        )
        raise SystemExit(1)

    if not quiet:
        console.print(f"[blue]Loading GFF from:[/blue] {helixer_gff}")
        console.print(f"[blue]Reference genome:[/blue] {genome}")
        if rnaseq_bam:
            console.print(f"[blue]RNA-seq BAM:[/blue] {rnaseq_bam}")
        if junctions_bed:
            console.print(f"[blue]Junctions BED:[/blue] {junctions_bed}")

    try:
        # Load genome
        genome_accessor = GenomeAccessor(genome)

        # Extract junctions
        if rnaseq_bam:
            from helixforge.io.bam import JunctionExtractor

            if not quiet:
                console.print("[blue]Extracting junctions from BAM...[/blue]")
            extractor = JunctionExtractor(rnaseq_bam)
            junctions = extractor.extract_all(min_reads=min_reads)
            if not quiet:
                n_junctions = sum(len(j) for j in junctions.values())
                console.print(f"[green]Extracted {n_junctions} junctions[/green]")
        else:
            # Load from BED
            from helixforge.io.bam import load_junctions_from_bed

            if not quiet:
                console.print("[blue]Loading junctions from BED...[/blue]")
            junctions = load_junctions_from_bed(junctions_bed)
            if not quiet:
                n_junctions = sum(len(j) for j in junctions.values())
                console.print(f"[green]Loaded {n_junctions} junctions[/green]")

        # Load genes
        parser = GFF3Parser(helixer_gff)
        genes = list(parser.iter_genes())
        if not quiet:
            console.print(f"[green]Loaded {len(genes)} genes[/green]")

        # Create refiner
        refiner = SpliceRefiner(
            genome_accessor,
            junctions,
            max_shift=max_shift,
            min_junction_reads=min_reads,
        )

        # Refine genes
        if not quiet:
            console.print(
                f"[blue]Refining splice sites with {workers} worker(s)...[/blue]"
            )

        refined_genes = []
        reports = []

        for refined_gene, splice_report in refiner.refine_genes_parallel(
            genes, n_workers=workers
        ):
            refined_genes.append(refined_gene)
            reports.append(splice_report)

        # Optionally adjust boundaries
        if adjust_boundaries:
            if not quiet:
                console.print("[blue]Adjusting start/stop codons...[/blue]")

            adjuster = BoundaryAdjuster(genome_accessor)
            adjusted_genes = []
            for gene in refined_genes:
                gene, _ = adjuster.adjust_start_codon(gene)
                gene, _ = adjuster.adjust_stop_codon(gene)
                adjusted_genes.append(gene)
            refined_genes = adjusted_genes

        # Write output GFF
        writer = GFF3Writer(output_gff)
        writer.write_genes(refined_genes)
        if not quiet:
            console.print(f"[green]Wrote refined GFF to:[/green] {output_gff}")

        # Write reports if requested
        if report:
            SpliceReportWriter.write_tsv(reports, report)
            if not quiet:
                console.print(f"[green]Wrote splice report to:[/green] {report}")

        if corrections_detail:
            SpliceReportWriter.write_corrections_detail(reports, corrections_detail)
            if not quiet:
                console.print(
                    f"[green]Wrote corrections detail to:[/green] {corrections_detail}"
                )

        if unsupported_bed:
            genes_dict = {g.gene_id: g for g in refined_genes}
            SpliceReportWriter.write_unsupported_introns(
                reports, genes_dict, unsupported_bed
            )
            if not quiet:
                console.print(
                    f"[green]Wrote unsupported introns to:[/green] {unsupported_bed}"
                )

        # Print summary
        if not quiet:
            stats = SpliceReportWriter.summary_statistics(reports)
            console.print("\n[bold]Summary:[/bold]")
            console.print(f"  Total genes: {stats['total_genes']}")
            console.print(f"  Total introns: {stats['total_introns']}")
            console.print(
                f"  Supported introns: {stats['supported_introns']} "
                f"({stats['support_rate']*100:.1f}%)"
            )
            console.print(
                f"  Corrections made: {stats['corrections_made']} "
                f"({stats['correction_rate']*100:.1f}%)"
            )
            console.print(
                f"  Canonical introns: {stats['canonical_introns']} "
                f"({stats['canonical_rate']*100:.1f}%)"
            )

        genome_accessor.close()

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise SystemExit(1)


# =============================================================================
# viz command
# =============================================================================


@main.command()
@click.option(
    "-g",
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input GFF3 file with gene predictions.",
)
@click.option(
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA for sequence display.",
)
@click.option(
    "--bam",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="RNA-seq BAM file(s) for coverage tracks.",
)
@click.option(
    "--region",
    type=str,
    help="Genomic region to visualize (e.g., chr1:1000-2000).",
)
@click.option(
    "--gene",
    type=str,
    help="Gene ID to visualize.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file for static visualization.",
)
@click.option(
    "--interactive",
    is_flag=True,
    help="Launch interactive visualization server.",
)
@click.option(
    "--port",
    type=int,
    default=8050,
    show_default=True,
    help="Port for interactive server.",
)
@click.pass_context
def viz(
    ctx: click.Context,
    gff: Path,
    genome: Optional[Path],
    bam: tuple[Path, ...],
    region: Optional[str],
    gene: Optional[str],
    output: Optional[Path],
    interactive: bool,
    port: int,
) -> None:
    """Visualize gene predictions and evidence.

    Generate static plots or launch an interactive browser for
    exploring gene models and supporting evidence.

    Example:
        $ helixforge viz -g refined.gff3 --gene GENE001 -o gene_plot.png
        $ helixforge viz -g refined.gff3 --interactive
    """
    # TODO: Implement visualization
    console.print("[yellow]viz command not yet implemented[/yellow]")
    raise SystemExit(1)


if __name__ == "__main__":
    main()
