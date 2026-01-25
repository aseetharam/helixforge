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
# Chunk-aware processing options
@click.option(
    "--region",
    type=str,
    default=None,
    help="Process only this region (format: seqid:start-end, 1-based inclusive).",
)
@click.option(
    "--chunk-id",
    type=str,
    default=None,
    help="Chunk identifier for logging and output naming in parallel mode.",
)
@click.option(
    "--scaffold",
    type=str,
    default=None,
    help="Process only this scaffold (simpler alternative to --region).",
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
    region: Optional[str],
    chunk_id: Optional[str],
    scaffold: Optional[str],
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

    For parallel execution, use --region or --scaffold to process a subset
    of genes, and --chunk-id for logging and output organization.

    Example:
        $ helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa -o scores.tsv
        $ helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa -o scores.tsv --bed scores.bed --distribution-plot dist.html
        $ helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa --region chr1:1-1000000 --chunk-id chunk_001 -o chunk_001.tsv
    """
    from helixforge.core.confidence import (
        ConfidenceCalculator,
        ConfidenceWriter,
    )
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser
    from helixforge.io.hdf5 import HelixerHDF5Reader
    from helixforge.utils.regions import (
        GenomicRegion,
        parse_region,
        region_from_scaffold,
        validate_region,
    )

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    # Log chunk ID if provided
    if chunk_id and not quiet:
        console.print(f"[blue]Chunk ID:[/blue] {chunk_id}")

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
                # Parse region constraints
                target_region: GenomicRegion | None = None

                if region:
                    try:
                        target_region = parse_region(region)
                        validate_region(target_region, genome_accessor)
                        if not quiet:
                            console.print(f"[blue]Processing region:[/blue] {target_region}")
                    except ValueError as e:
                        console.print(f"[red]Error:[/red] {e}")
                        raise SystemExit(1)
                elif scaffold:
                    # Scaffold-only mode: process entire scaffold
                    if scaffold not in genome_accessor.scaffold_lengths:
                        console.print(
                            f"[red]Error:[/red] Scaffold '{scaffold}' not found in genome"
                        )
                        raise SystemExit(1)
                    scaffold_len = genome_accessor.scaffold_lengths[scaffold]
                    target_region = region_from_scaffold(scaffold, scaffold_len)
                    if not quiet:
                        console.print(
                            f"[blue]Processing scaffold:[/blue] {scaffold} "
                            f"(length: {scaffold_len:,})"
                        )

                # Load and filter genes
                parser = GFF3Parser(gff)

                if target_region:
                    genes = parser.get_genes_in_region(
                        target_region.seqid,
                        target_region.start,
                        target_region.end,
                    )
                    if not quiet:
                        console.print(f"[green]Found {len(genes)} genes in region[/green]")
                else:
                    genes = list(parser.iter_genes())
                    if not quiet:
                        console.print(f"[green]Loaded {len(genes)} genes[/green]")

                # Handle empty gene list
                if not genes:
                    if not quiet:
                        console.print("[yellow]No genes found in specified region[/yellow]")
                    # Write empty output with header
                    ConfidenceWriter.to_tsv([], output)
                    if not quiet:
                        console.print(f"[green]Wrote empty TSV to:[/green] {output}")
                    return

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
# Chunk-aware processing options
@click.option(
    "--region",
    type=str,
    default=None,
    help="Process only this region (format: seqid:start-end, 1-based inclusive).",
)
@click.option(
    "--chunk-id",
    type=str,
    default=None,
    help="Chunk identifier for logging and output naming in parallel mode.",
)
@click.option(
    "--scaffold",
    type=str,
    default=None,
    help="Process only this scaffold (simpler alternative to --region).",
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
    region: Optional[str],
    chunk_id: Optional[str],
    scaffold: Optional[str],
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

    For parallel execution, use --region or --scaffold to process a subset
    of genes, and --chunk-id for logging and output organization.

    Example:
        $ helixforge splice --helixer-gff helixer.gff3 --genome genome.fa \\
            --rnaseq-bam rnaseq.bam -o refined.gff3 -r splice_report.tsv
        $ helixforge splice --helixer-gff helixer.gff3 --genome genome.fa \\
            --rnaseq-bam rnaseq.bam --region chr1:1-1000000 --chunk-id chunk_001 -o chunk_001.gff3
    """
    from helixforge.core.boundaries import BoundaryAdjuster
    from helixforge.core.splice import (
        SpliceRefiner,
        SpliceReportWriter,
    )
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser, GFF3Writer
    from helixforge.utils.regions import (
        GenomicRegion,
        parse_region,
        region_from_scaffold,
        validate_region,
    )

    verbose = ctx.obj.get("verbose", False) or verbose_flag
    quiet = ctx.obj.get("quiet", False)

    # Log chunk ID if provided
    if chunk_id and not quiet:
        console.print(f"[blue]Chunk ID:[/blue] {chunk_id}")

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

        # Parse region constraints
        target_region: GenomicRegion | None = None

        if region:
            try:
                target_region = parse_region(region)
                validate_region(target_region, genome_accessor)
                if not quiet:
                    console.print(f"[blue]Processing region:[/blue] {target_region}")
            except ValueError as e:
                console.print(f"[red]Error:[/red] {e}")
                genome_accessor.close()
                raise SystemExit(1)
        elif scaffold:
            # Scaffold-only mode: process entire scaffold
            if scaffold not in genome_accessor.scaffold_lengths:
                console.print(
                    f"[red]Error:[/red] Scaffold '{scaffold}' not found in genome"
                )
                genome_accessor.close()
                raise SystemExit(1)
            scaffold_len = genome_accessor.scaffold_lengths[scaffold]
            target_region = region_from_scaffold(scaffold, scaffold_len)
            if not quiet:
                console.print(
                    f"[blue]Processing scaffold:[/blue] {scaffold} "
                    f"(length: {scaffold_len:,})"
                )

        # Extract junctions
        if rnaseq_bam:
            from helixforge.io.bam import JunctionExtractor

            if not quiet:
                console.print("[blue]Extracting junctions from BAM...[/blue]")
            extractor = JunctionExtractor(rnaseq_bam)

            if target_region:
                # Extract only for the target region
                junctions = {
                    target_region.seqid: extractor.extract_region(
                        target_region.seqid,
                        target_region.start,
                        target_region.end,
                        min_reads=min_reads,
                    )
                }
            else:
                junctions = extractor.extract_all(min_reads=min_reads)

            if not quiet:
                n_junctions = sum(len(j) for j in junctions.values())
                console.print(f"[green]Extracted {n_junctions} junctions[/green]")
        else:
            # Load from BED
            from helixforge.io.bam import load_junctions_from_bed

            if not quiet:
                console.print("[blue]Loading junctions from BED...[/blue]")
            all_junctions = load_junctions_from_bed(junctions_bed)

            if target_region:
                # Filter junctions to region
                junctions = {}
                for seqid, seqid_junctions in all_junctions.items():
                    if seqid == target_region.seqid:
                        filtered = [
                            j for j in seqid_junctions
                            if j.start >= target_region.start
                            and j.end <= target_region.end
                        ]
                        if filtered:
                            junctions[seqid] = filtered
            else:
                junctions = all_junctions

            if not quiet:
                n_junctions = sum(len(j) for j in junctions.values())
                console.print(f"[green]Loaded {n_junctions} junctions[/green]")

        # Load and filter genes
        parser = GFF3Parser(helixer_gff)

        if target_region:
            genes = parser.get_genes_in_region(
                target_region.seqid,
                target_region.start,
                target_region.end,
            )
            if not quiet:
                console.print(f"[green]Found {len(genes)} genes in region[/green]")
        else:
            genes = list(parser.iter_genes())
            if not quiet:
                console.print(f"[green]Loaded {len(genes)} genes[/green]")

        # Handle empty gene list
        if not genes:
            if not quiet:
                console.print("[yellow]No genes found in specified region[/yellow]")
            # Write empty GFF
            with GFF3Writer(output_gff) as writer:
                writer.write_header()
            if not quiet:
                console.print(f"[green]Wrote empty GFF to:[/green] {output_gff}")
            genome_accessor.close()
            return

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


# =============================================================================
# parallel command group
# =============================================================================


@main.group()
def parallel():
    """Parallel execution for large genome processing.

    Commands for managing parallel processing of large genomes,
    including chunking strategies and task file generation.

    Recommended workflow:
        1. Create a chunk plan: helixforge parallel plan
        2. Generate task file: helixforge parallel tasks
        3. Execute with HyperShell: hs launch --parallelism 32 < tasks.txt
        4. Aggregate results: helixforge parallel aggregate
    """
    pass


@parallel.command("plan")
@click.option(
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Genome FASTA file.",
)
@click.option(
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    help="GFF3 file for gene-aware chunking.",
)
@click.option(
    "--strategy",
    type=click.Choice(["scaffold", "size", "genes", "adaptive"]),
    default="scaffold",
    show_default=True,
    help="Chunking strategy.",
)
@click.option(
    "--chunk-size",
    type=int,
    help="Chunk size (bases for 'size', gene count for 'genes').",
)
@click.option(
    "--min-chunk-size",
    type=int,
    default=100000,
    show_default=True,
    help="Minimum chunk size in bases.",
)
@click.option(
    "--max-chunk-size",
    type=int,
    help="Maximum chunk size (for splitting large scaffolds).",
)
@click.option(
    "--target-chunks",
    type=int,
    help="Target number of chunks for 'adaptive' strategy.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output JSON file for chunk plan.",
)
def create_chunk_plan(
    genome: Path,
    gff: Optional[Path],
    strategy: str,
    chunk_size: Optional[int],
    min_chunk_size: int,
    max_chunk_size: Optional[int],
    target_chunks: Optional[int],
    output: Path,
) -> None:
    """Create a chunking plan for parallel execution.

    Generates a JSON file containing genomic chunk definitions that
    can be used with HyperShell, GNU Parallel, or local multiprocessing.

    Example:
        $ helixforge parallel plan --genome genome.fa -o chunks.json
        $ helixforge parallel plan --genome genome.fa --gff genes.gff3 --strategy genes -o chunks.json
    """
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.parallel.chunker import GenomeChunker, ChunkStrategy

    try:
        genome_accessor = GenomeAccessor(genome)

        gff_parser = None
        if gff:
            from helixforge.io.gff import GFF3Parser
            gff_parser = GFF3Parser(gff)

        chunker = GenomeChunker(genome_accessor, gff_parser)

        plan = chunker.create_plan(
            strategy=ChunkStrategy(strategy),
            chunk_size=chunk_size,
            min_chunk_size=min_chunk_size,
            max_chunk_size=max_chunk_size,
            target_chunks=target_chunks,
        )

        plan.save(output)

        # Print summary
        summary = plan.summary()
        console.print(f"[green]Created chunk plan:[/green] {output}")
        console.print(f"  Strategy: {summary['strategy']}")
        console.print(f"  Number of chunks: {summary['n_chunks']}")
        console.print(f"  Total bases: {summary['total_bases']:,}")
        console.print(f"  Mean chunk size: {summary['mean_chunk_size']:,.0f}")
        console.print(f"  Min chunk size: {summary['min_chunk_size']:,}")
        console.print(f"  Max chunk size: {summary['max_chunk_size']:,}")

        genome_accessor.close()

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise SystemExit(1)


@parallel.command("tasks")
@click.option(
    "--chunk-plan",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Chunk plan JSON file.",
)
@click.option(
    "--command",
    required=True,
    help="Command template with placeholders: {chunk_id}, {seqid}, {start}, {end} (1-based for CLI), {start_0}, {end_0} (0-based), {size}, {output_dir}.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output task file.",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    help="Output directory for results (creates {output_dir} placeholder).",
)
@click.option(
    "--wrapper",
    type=click.Path(path_type=Path),
    help="Generate wrapper script for complex workflows.",
)
@click.option(
    "--wrapper-setup",
    multiple=True,
    help="Setup commands for wrapper (module loads, conda activate, etc.).",
)
@click.option(
    "--include-logging",
    is_flag=True,
    help="Add stdout/stderr redirection per task.",
)
def generate_tasks(
    chunk_plan: Path,
    command: str,
    output: Path,
    output_dir: Optional[Path],
    wrapper: Optional[Path],
    wrapper_setup: tuple[str, ...],
    include_logging: bool,
) -> None:
    """Generate task file for parallel execution.

    Output can be used with HyperShell, GNU Parallel, or xargs.

    \b
    The command template can use these placeholders:
    - {chunk_id}: Unique chunk identifier
    - {seqid}: Scaffold/chromosome name
    - {start}: Start coordinate (1-based, for CLI --region flags)
    - {end}: End coordinate (1-based inclusive, for CLI --region flags)
    - {start_0}: Start coordinate (0-based, for internal use)
    - {end_0}: End coordinate (0-based exclusive, for internal use)
    - {size}: Chunk size in bases
    - {output_dir}: Output directory (if --output-dir provided)

    Note: {start} and {end} output 1-based coordinates suitable for CLI
    commands like --region chr1:1000-2000. The internal 0-based half-open
    coordinates are converted automatically.

    Examples:
        # Simple task file
        $ helixforge parallel tasks --chunk-plan chunks.json \\
            --command 'helixforge confidence --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o {output_dir}/{chunk_id}.tsv' \\
            --output tasks.txt --output-dir outputs/

        # Execute with HyperShell
        $ hs launch --parallelism 32 < tasks.txt

        # Execute with GNU Parallel
        $ parallel -j 32 < tasks.txt

        # With wrapper script for complex workflows
        $ helixforge parallel tasks --chunk-plan chunks.json \\
            --wrapper wrapper.sh --wrapper-setup 'module load python' \\
            --output tasks.txt
    """
    from helixforge.parallel.chunker import ChunkPlan
    from helixforge.parallel.taskgen import TaskGenerator, generate_hypershell_command

    try:
        plan = ChunkPlan.load(chunk_plan)
        gen = TaskGenerator(plan)

        if wrapper:
            # Generate wrapper script and task file
            wrapper_cmds = [command]
            task_file = gen.generate_with_wrapper(
                output_path=output,
                wrapper_script=wrapper,
                chunk_plan_output=chunk_plan,  # Keep using original plan
            )

            # Also generate the wrapper script
            TaskGenerator.generate_wrapper_script(
                output_path=wrapper,
                chunk_plan_path=chunk_plan,
                commands=wrapper_cmds,
                setup_commands=list(wrapper_setup) if wrapper_setup else None,
                output_dir=output_dir,
            )
            console.print(f"[green]Generated wrapper script:[/green] {wrapper}")
        else:
            # Generate direct task file
            task_file = gen.generate(
                command_template=command,
                output_path=output,
                output_dir=output_dir,
                include_logging=include_logging,
            )

        console.print(f"[green]Generated task file:[/green] {output}")
        console.print(f"  Number of tasks: {task_file.n_tasks}")

        # Show preview
        preview = task_file.preview(3)
        if preview:
            console.print("\n[bold]Preview (first 3 tasks):[/bold]")
            for i, line in enumerate(preview, 1):
                console.print(f"  {i}. {line[:80]}{'...' if len(line) > 80 else ''}")

        # Show execution hint
        console.print(f"\n[blue]Execute with HyperShell:[/blue]")
        console.print(f"  {generate_hypershell_command(output, parallelism=8)}")
        console.print(f"\n[blue]Or with GNU Parallel:[/blue]")
        console.print(f"  parallel -j 8 < {output}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise SystemExit(1)


@parallel.command("aggregate")
@click.option(
    "--input-dir",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Directory containing chunk outputs.",
)
@click.option(
    "--pattern",
    required=True,
    help="Glob pattern for chunk output files (e.g., '*.gff3').",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output file for aggregated results.",
)
@click.option(
    "--type",
    "agg_type",
    type=click.Choice(["concat", "merge_gff", "merge_tsv"]),
    default="concat",
    show_default=True,
    help="Aggregation method.",
)
@click.option(
    "--sort-by",
    type=click.Choice(["name", "position"]),
    help="Sort order for merged output (for GFF3).",
)
def aggregate_outputs(
    input_dir: Path,
    pattern: str,
    output: Path,
    agg_type: str,
    sort_by: Optional[str],
) -> None:
    """Aggregate outputs from parallel chunk processing.

    Combines output files from HyperShell, GNU Parallel, or local
    parallel processing into a single output file.

    Example:
        $ helixforge parallel aggregate --input-dir outputs/ \\
            --pattern '*.gff3' -o combined.gff3 --type merge_gff

        $ helixforge parallel aggregate --input-dir outputs/ \\
            --pattern '*.tsv' -o combined.tsv --type merge_tsv
    """
    try:
        # Find matching files
        files = sorted(input_dir.glob(pattern))

        if not files:
            console.print(f"[red]No files matching '{pattern}' found in {input_dir}[/red]")
            raise SystemExit(1)

        console.print(f"Found {len(files)} files matching '{pattern}'")

        output.parent.mkdir(parents=True, exist_ok=True)

        if agg_type == "concat":
            # Simple concatenation
            with open(output, "w") as out_f:
                for f in files:
                    out_f.write(f.read_text())
                    if not f.read_text().endswith("\n"):
                        out_f.write("\n")

        elif agg_type == "merge_gff":
            # Merge GFF3 files with proper header handling
            seen_header = False
            with open(output, "w") as out_f:
                for f in files:
                    for line in f.read_text().splitlines():
                        if line.startswith("##gff-version"):
                            if not seen_header:
                                out_f.write(line + "\n")
                                seen_header = True
                        elif line.startswith("###"):
                            continue  # Skip section separators
                        else:
                            out_f.write(line + "\n")
                out_f.write("###\n")  # Final terminator

        elif agg_type == "merge_tsv":
            # Merge TSV files keeping only first header
            seen_header = False
            with open(output, "w") as out_f:
                for f in files:
                    lines = f.read_text().splitlines()
                    for i, line in enumerate(lines):
                        if i == 0 and line.startswith("#"):
                            if not seen_header:
                                out_f.write(line + "\n")
                                seen_header = True
                        elif i == 0 and not seen_header:
                            # First file, assume header
                            out_f.write(line + "\n")
                            seen_header = True
                        elif i == 0:
                            continue  # Skip header in subsequent files
                        else:
                            out_f.write(line + "\n")

        console.print(f"[green]Aggregated {len(files)} files to:[/green] {output}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise SystemExit(1)


@parallel.command("example-sbatch")
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output path for SBATCH script.",
)
@click.option(
    "--executor",
    type=click.Choice(["hypershell", "parallel"]),
    default="hypershell",
    show_default=True,
    help="Which parallel executor to use.",
)
def write_example_sbatch(
    output: Path,
    executor: str,
) -> None:
    """Write example SLURM SBATCH script.

    Generates a template SBATCH script for running HyperShell or
    GNU Parallel on a SLURM cluster. Customize the script for your
    cluster's configuration.

    Example:
        $ helixforge parallel example-sbatch -o run_helixforge.sbatch
        $ helixforge parallel example-sbatch -o run_helixforge.sbatch --executor parallel
    """
    from helixforge.parallel.slurm import write_example_sbatch as _write_sbatch

    try:
        _write_sbatch(output, executor=executor)
        console.print(f"[green]Wrote example SBATCH script:[/green] {output}")
        console.print("\n[yellow]Remember to customize for your cluster:[/yellow]")
        console.print("  - Update partition name")
        console.print("  - Adjust time, memory, and CPU allocation")
        console.print("  - Add module loads for your environment")
        console.print("  - Set up conda/virtualenv activation")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise SystemExit(1)


@parallel.command("suggest")
@click.option(
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Genome FASTA file.",
)
@click.option(
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    help="GFF3 file for gene count.",
)
@click.option(
    "--memory",
    type=float,
    default=64.0,
    show_default=True,
    help="Available memory in GB.",
)
@click.option(
    "--workers",
    type=int,
    default=8,
    show_default=True,
    help="Number of parallel workers.",
)
def suggest_parameters(
    genome: Path,
    gff: Optional[Path],
    memory: float,
    workers: int,
) -> None:
    """Suggest optimal chunking parameters.

    Analyzes genome size and gene count to recommend chunking
    strategy and parameters based on available resources.

    Example:
        $ helixforge parallel suggest --genome genome.fa --memory 128 --workers 16
    """
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.parallel.chunker import suggest_chunk_parameters

    try:
        genome_accessor = GenomeAccessor(genome)
        genome_size = genome_accessor.total_length

        n_genes = 0
        if gff:
            from helixforge.io.gff import GFF3Parser
            parser = GFF3Parser(gff)
            n_genes = sum(1 for _ in parser.iter_genes())

        suggestion = suggest_chunk_parameters(
            genome_size=genome_size,
            n_genes=n_genes,
            available_memory_gb=memory,
            n_workers=workers,
        )

        console.print("[bold]Suggested Chunking Parameters:[/bold]")
        console.print(f"  Genome size: {genome_size:,} bp")
        console.print(f"  Gene count: {n_genes:,}")
        console.print(f"  Strategy: {suggestion['strategy'].value}")
        if suggestion.get('chunk_size'):
            console.print(f"  Chunk size: {suggestion['chunk_size']:,}")
        if suggestion.get('target_chunks'):
            console.print(f"  Target chunks: {suggestion['target_chunks']}")
        console.print(f"  Memory per worker: {suggestion['memory_per_worker_mb']:.0f} MB")
        console.print(f"  Max chunk bases: {suggestion['max_chunk_bases']:,}")
        console.print(f"\n  [blue]Rationale:[/blue] {suggestion['rationale']}")

        genome_accessor.close()

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise SystemExit(1)


# =============================================================================
# homology command group
# =============================================================================


@main.group()
def homology():
    """Homology-based validation of gene predictions.

    Commands for searching predicted proteins against reference databases,
    validating gene models, and detecting issues like chimeras and fragments.

    Recommended workflow:
        1. Extract proteins: helixforge homology extract-proteins
        2. Search database: helixforge homology search
        3. Validate genes: helixforge homology validate
    """
    pass


@homology.command("extract-proteins")
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
    required=True,
    help="Reference genome FASTA file.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output protein FASTA file.",
)
@click.option(
    "--longest-isoform",
    is_flag=True,
    default=True,
    show_default=True,
    help="Output only longest isoform per gene.",
)
@click.pass_context
def extract_proteins(
    ctx: click.Context,
    gff: Path,
    genome: Path,
    output: Path,
    longest_isoform: bool,
) -> None:
    """Extract protein sequences from gene models.

    Translates CDS features from GFF3 gene models and writes
    protein sequences in FASTA format for homology searching.

    Example:
        $ helixforge homology extract-proteins -g genes.gff3 --genome genome.fa -o proteins.fa
    """
    from helixforge.homology.search import extract_proteins_from_gff
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]Loading GFF from:[/blue] {gff}")
        console.print(f"[blue]Reference genome:[/blue] {genome}")

    try:
        with GenomeAccessor(genome) as genome_accessor:
            parser = GFF3Parser(gff)

            protein_lengths = extract_proteins_from_gff(
                parser,
                genome_accessor,
                output,
                longest_isoform=longest_isoform,
            )

            if not quiet:
                console.print(f"[green]Extracted {len(protein_lengths)} proteins to:[/green] {output}")
                total_aa = sum(protein_lengths.values())
                console.print(f"  Total amino acids: {total_aa:,}")
                avg_len = total_aa / len(protein_lengths) if protein_lengths else 0
                console.print(f"  Average length: {avg_len:.1f} aa")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@homology.command("format-db")
@click.option(
    "-i",
    "--input",
    "input_fasta",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input protein FASTA file.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output database path (auto-generated if omitted).",
)
@click.option(
    "--tool",
    type=click.Choice(["diamond", "mmseqs2"]),
    default="diamond",
    show_default=True,
    help="Search tool to format database for.",
)
@click.option(
    "-j",
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Number of threads.",
)
@click.pass_context
def format_db(
    ctx: click.Context,
    input_fasta: Path,
    output: Optional[Path],
    tool: str,
    threads: int,
) -> None:
    """Format protein database for homology searching.

    Creates a Diamond or MMseqs2 database from a FASTA file.

    Example:
        $ helixforge homology format-db -i swissprot.fa -o swissprot
        $ helixforge homology format-db -i proteins.fa --tool mmseqs2 -o proteins_db
    """
    from helixforge.homology.search import HomologySearch, SearchTool

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]Formatting database from:[/blue] {input_fasta}")
        console.print(f"[blue]Tool:[/blue] {tool}")

    try:
        searcher = HomologySearch(
            tool=SearchTool(tool),
            threads=threads,
        )

        formatted_path = searcher.format_database(input_fasta, output)

        if not quiet:
            console.print(f"[green]Created database:[/green] {formatted_path}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@homology.command("download-db")
@click.option(
    "--database",
    type=click.Choice([
        "swissprot",
        "swissprot_plants",
        "swissprot_fungi",
        "uniref90",
        "uniref50",
    ]),
    required=True,
    help="Database to download.",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(path_type=Path),
    help="Output directory (default: ~/.helixforge/databases).",
)
@click.option(
    "--format/--no-format",
    default=True,
    show_default=True,
    help="Format database for Diamond after download.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Force re-download even if cached.",
)
@click.pass_context
def download_db(
    ctx: click.Context,
    database: str,
    output_dir: Optional[Path],
    format: bool,
    force: bool,
) -> None:
    """Download reference protein database.

    Downloads protein databases from UniProt for homology validation.

    Available databases:
        - swissprot: Swiss-Prot complete (high-quality, curated)
        - swissprot_plants: Swiss-Prot plant proteins only
        - swissprot_fungi: Swiss-Prot fungal proteins only
        - uniref90: UniRef90 clustered (larger coverage)
        - uniref50: UniRef50 clustered (smaller, faster)

    Example:
        $ helixforge homology download-db --database swissprot_plants
        $ helixforge homology download-db --database uniref90 -o /data/databases
    """
    from helixforge.homology.databases import DatabaseManager

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]Downloading database:[/blue] {database}")

    try:
        manager = DatabaseManager(
            cache_dir=output_dir,
            auto_format=format,
        )

        db_info = manager.get_database(database, force_download=force)

        if not quiet:
            console.print(f"[green]Downloaded:[/green] {db_info.path}")
            if db_info.n_sequences:
                console.print(f"  Sequences: {db_info.n_sequences:,}")
            if db_info.formatted_path:
                console.print(f"  Formatted database: {db_info.formatted_path}")
            if db_info.download_date:
                console.print(f"  Download date: {db_info.download_date.isoformat()}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@homology.command("search")
@click.option(
    "-q",
    "--query",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Query protein FASTA file.",
)
@click.option(
    "-d",
    "--database",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Diamond database (.dmnd) or MMseqs2 database.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output file for search results (TSV format).",
)
@click.option(
    "--tool",
    type=click.Choice(["diamond", "mmseqs2"]),
    default="diamond",
    show_default=True,
    help="Search tool to use.",
)
@click.option(
    "-e",
    "--evalue",
    type=float,
    default=1e-5,
    show_default=True,
    help="E-value threshold.",
)
@click.option(
    "--max-hits",
    type=int,
    default=10,
    show_default=True,
    help="Maximum hits per query.",
)
@click.option(
    "--sensitive",
    is_flag=True,
    help="Use sensitive search mode (slower but more hits).",
)
@click.option(
    "-j",
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Number of threads.",
)
@click.pass_context
def search(
    ctx: click.Context,
    query: Path,
    database: Path,
    output: Path,
    tool: str,
    evalue: float,
    max_hits: int,
    sensitive: bool,
    threads: int,
) -> None:
    """Run homology search against protein database.

    Searches query proteins against a reference database using
    Diamond or MMseqs2 and outputs results in tabular format.

    Example:
        $ helixforge homology search -q proteins.fa -d swissprot.dmnd -o hits.tsv
        $ helixforge homology search -q proteins.fa -d swissprot.dmnd -o hits.tsv --sensitive
    """
    from helixforge.homology.search import HomologySearch, SearchTool

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]Query:[/blue] {query}")
        console.print(f"[blue]Database:[/blue] {database}")
        console.print(f"[blue]Tool:[/blue] {tool}")
        console.print(f"[blue]E-value:[/blue] {evalue}")
        console.print(f"[blue]Threads:[/blue] {threads}")

    try:
        searcher = HomologySearch(
            tool=SearchTool(tool),
            database=database,
            threads=threads,
            evalue=evalue,
            max_target_seqs=max_hits,
            sensitive=sensitive,
        )

        result_path = searcher.search(query, output)

        # Parse to count hits
        hits_by_query = searcher.parse_results(result_path)
        n_queries_with_hits = len(hits_by_query)
        total_hits = sum(len(h) for h in hits_by_query.values())

        if not quiet:
            console.print(f"[green]Search complete:[/green] {result_path}")
            console.print(f"  Queries with hits: {n_queries_with_hits}")
            console.print(f"  Total hits: {total_hits}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@homology.command("validate")
@click.option(
    "-g",
    "--gff",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input GFF3 file with gene predictions.",
)
@click.option(
    "-s",
    "--search-results",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Homology search results (TSV from 'homology search').",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output validation report (TSV).",
)
@click.option(
    "--te-bed",
    type=click.Path(exists=True, path_type=Path),
    help="BED file with transposable element annotations.",
)
@click.option(
    "--thresholds",
    type=click.Choice(["default", "strict", "relaxed"]),
    default="default",
    show_default=True,
    help="Validation threshold preset.",
)
@click.option(
    "--chimera-report",
    type=click.Path(path_type=Path),
    help="Output report for chimeric genes.",
)
@click.option(
    "--fragment-report",
    type=click.Path(path_type=Path),
    help="Output report for fragment groups.",
)
@click.option(
    "--output-gff",
    type=click.Path(path_type=Path),
    help="Output GFF3 with validation status attributes.",
)
@click.pass_context
def validate_homology(
    ctx: click.Context,
    gff: Path,
    search_results: Path,
    output: Path,
    te_bed: Optional[Path],
    thresholds: str,
    chimera_report: Optional[Path],
    fragment_report: Optional[Path],
    output_gff: Optional[Path],
) -> None:
    """Validate gene predictions using homology evidence.

    Analyzes homology search results to classify genes and detect
    issues like chimeras, fragments, and TE overlaps.

    \b
    Validation status:
    - complete: Good coverage of known protein
    - partial: Incomplete coverage
    - no_hit: No significant homology
    - chimeric: Appears to be fusion of multiple genes
    - fragmented: Appears to be part of split gene
    - te_overlap: Overlaps transposable element

    Example:
        $ helixforge homology validate -g genes.gff3 -s hits.tsv -o validation.tsv
        $ helixforge homology validate -g genes.gff3 -s hits.tsv -o validation.tsv --te-bed te.bed --chimera-report chimeras.tsv
    """
    from helixforge.homology.search import HomologySearch, SearchTool, get_sequence_lengths
    from helixforge.homology.validate import (
        HomologyValidator,
        HomologyStatus,
        ValidationThresholds,
        load_te_annotations,
        summarize_validation,
    )
    from helixforge.io.gff import GFF3Parser, GFF3Writer

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]GFF:[/blue] {gff}")
        console.print(f"[blue]Search results:[/blue] {search_results}")
        console.print(f"[blue]Thresholds:[/blue] {thresholds}")

    try:
        # Load gene models
        parser = GFF3Parser(gff)
        genes = list(parser.iter_genes())
        if not quiet:
            console.print(f"[green]Loaded {len(genes)} genes[/green]")

        # Create gene info mapping for TE overlap
        gene_info = {
            g.gene_id: (g.seqid, g.start, g.end)
            for g in genes
        }

        # Load TE annotations if provided
        te_intervals = {}
        if te_bed:
            if not quiet:
                console.print(f"[blue]Loading TE annotations from:[/blue] {te_bed}")
            te_intervals = load_te_annotations(te_bed)

        # Select thresholds
        if thresholds == "strict":
            threshold_obj = ValidationThresholds.strict()
        elif thresholds == "relaxed":
            threshold_obj = ValidationThresholds.relaxed()
        else:
            threshold_obj = ValidationThresholds.default()

        # Parse search results
        searcher = HomologySearch(tool=SearchTool.DIAMOND)
        hits_by_gene = searcher.parse_results(search_results)
        if not quiet:
            console.print(f"[green]Loaded hits for {len(hits_by_gene)} genes[/green]")

        # Create validator
        validator = HomologyValidator(
            thresholds=threshold_obj,
            te_annotations=te_intervals,
        )

        # Validate
        results = validator.validate_from_search(hits_by_gene, gene_info)

        # Write main report
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w") as f:
            # Write header
            f.write("\t".join([
                "gene_id", "status", "best_hit_id", "n_hits",
                "query_coverage", "subject_coverage", "identity", "evalue",
                "is_chimeric", "fragment_group", "te_overlap_fraction", "flags", "notes"
            ]) + "\n")

            for gene_id, result in sorted(results.items()):
                row = [
                    gene_id,
                    result.status.value,
                    result.best_hit.subject_id if result.best_hit else "",
                    str(result.n_hits),
                    f"{result.query_coverage:.3f}" if result.query_coverage else "",
                    f"{result.subject_coverage:.3f}" if result.subject_coverage else "",
                    f"{result.identity:.1f}" if result.identity else "",
                    f"{result.evalue:.2e}" if result.evalue else "",
                    "yes" if result.chimeric_evidence else "no",
                    result.fragment_group or "",
                    f"{result.te_overlap_fraction:.3f}",
                    ";".join(result.flags),
                    ";".join(result.notes),
                ]
                f.write("\t".join(row) + "\n")

        if not quiet:
            console.print(f"[green]Wrote validation report to:[/green] {output}")

        # Write chimera report if requested
        if chimera_report:
            chimeric = [r for r in results.values() if r.chimeric_evidence]
            with open(chimera_report, "w") as f:
                f.write("gene_id\tsubject_a\tsubject_b\tgap\toverlap\tconfidence\n")
                for r in chimeric:
                    ev = r.chimeric_evidence
                    f.write(f"{ev.gene_id}\t{ev.subject_a}\t{ev.subject_b}\t{ev.gap}\t{ev.overlap}\t{ev.confidence:.3f}\n")
            if not quiet:
                console.print(f"[green]Wrote chimera report ({len(chimeric)} chimeras):[/green] {chimera_report}")

        # Write fragment report if requested
        if fragment_report:
            fragmented = [r for r in results.values() if r.fragment_group]
            # Group by fragment group
            groups = {}
            for r in fragmented:
                if r.fragment_group not in groups:
                    groups[r.fragment_group] = []
                groups[r.fragment_group].append(r.gene_id)

            with open(fragment_report, "w") as f:
                f.write("fragment_group\tgene_ids\tn_fragments\n")
                for group_id, gene_ids in sorted(groups.items()):
                    f.write(f"{group_id}\t{','.join(gene_ids)}\t{len(gene_ids)}\n")
            if not quiet:
                console.print(f"[green]Wrote fragment report ({len(groups)} groups):[/green] {fragment_report}")

        # Write annotated GFF if requested
        if output_gff:
            # Add validation attributes to genes
            gene_map = {g.gene_id: g for g in genes}
            for gene_id, result in results.items():
                if gene_id in gene_map:
                    gene = gene_map[gene_id]
                    gene.attributes["homology_status"] = result.status.value
                    if result.best_hit:
                        gene.attributes["best_hit"] = result.best_hit.subject_id
                    if result.flags:
                        gene.attributes["homology_flags"] = ",".join(result.flags)

            writer = GFF3Writer(output_gff)
            writer.write_genes(genes)
            if not quiet:
                console.print(f"[green]Wrote annotated GFF to:[/green] {output_gff}")

        # Print summary
        if not quiet:
            summary = summarize_validation(results)
            console.print("\n[bold]Summary:[/bold]")
            console.print(f"  Total genes: {summary['total']}")
            console.print(f"  [green]Complete:[/green] {summary['complete']} ({100*summary['complete']/summary['total']:.1f}%)")
            console.print(f"  [yellow]Partial:[/yellow] {summary['partial']} ({100*summary['partial']/summary['total']:.1f}%)")
            console.print(f"  [red]No hit:[/red] {summary['no_hit']} ({100*summary['no_hit']/summary['total']:.1f}%)")
            console.print(f"  [red]Chimeric:[/red] {summary['chimeric']}")
            console.print(f"  [red]Fragmented:[/red] {summary['fragmented']}")
            console.print(f"  [red]TE overlap:[/red] {summary['te_overlap']}")
            console.print(f"  With homology: {summary['pct_with_homology']:.1f}%")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@homology.command("list-databases")
def list_available_databases() -> None:
    """List available databases for download.

    Shows predefined databases that can be downloaded with
    'helixforge homology download-db'.

    Example:
        $ helixforge homology list-databases
    """
    from helixforge.homology.databases import list_databases

    descriptions = list_databases()

    console.print("[bold]Available Databases:[/bold]\n")
    for name, description in descriptions.items():
        console.print(f"  [blue]{name}[/blue]")
        console.print(f"    {description}\n")

    console.print("[bold]Usage:[/bold]")
    console.print("  $ helixforge homology download-db --database swissprot_plants")


if __name__ == "__main__":
    main()
