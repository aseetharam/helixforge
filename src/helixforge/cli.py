"""Command-line interface for HelixForge.

This module provides the main entry point for the helixforge CLI tool.
It uses Click to define commands and subcommands for various operations.

Commands:
    refine: Main refinement pipeline (requires HDF5 + RNA-seq)
    evidence: Score genes with RNA-seq evidence (scoring only)
    qc: Generate quality control reports
    validate: Homology-based validation
    confidence: Calculate confidence scores for gene predictions
    viz: Generate visualizations

Example:
    $ helixforge --help
    $ helixforge refine -p helixer.h5 -g helixer.gff3 --genome genome.fa --rnaseq-bam rnaseq.bam -o refined.gff3
    $ helixforge evidence -g predictions.gff3 -b rnaseq.bam -o annotated.gff3
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
    "--helixer-h5",
    "-p",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Helixer HDF5 predictions file (required for confidence scoring).",
)
@click.option(
    "--helixer-gff",
    "-g",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Helixer GFF3 predictions file.",
)
@click.option(
    "--genome",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Reference genome FASTA file.",
)
@click.option(
    "--rnaseq-bam",
    type=str,
    multiple=True,
    help="RNA-seq BAM file(s). Comma-separated or repeated. Required unless --junctions-bed.",
)
@click.option(
    "--rnaseq-bam-list",
    type=click.Path(exists=True, path_type=Path),
    help="File containing BAM paths (one per line).",
)
@click.option(
    "--junctions-bed",
    type=str,
    multiple=True,
    help="Junction file(s) in BED or STAR SJ.out.tab format. Alternative to BAM.",
)
@click.option(
    "--junctions-list",
    type=click.Path(exists=True, path_type=Path),
    help="File containing junction file paths (one per line).",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output refined GFF3 file.",
)
@click.option(
    "-r",
    "--report",
    type=click.Path(path_type=Path),
    help="Output refine report TSV with all scores and corrections.",
)
@click.option(
    "--splice-details",
    type=click.Path(path_type=Path),
    help="Output detailed splice corrections TSV.",
)
@click.option(
    "--evidence-details",
    type=click.Path(path_type=Path),
    help="Output per-gene evidence details TSV.",
)
@click.option(
    "--unsupported-bed",
    type=click.Path(path_type=Path),
    help="Output BED of introns without RNA-seq support.",
)
# Splice refinement options
@click.option(
    "--max-shift",
    type=int,
    default=15,
    show_default=True,
    help="Maximum splice site correction distance (bp).",
)
@click.option(
    "--min-reads",
    type=int,
    default=3,
    show_default=True,
    help="Minimum reads to support a splice junction.",
)
@click.option(
    "--min-tissues",
    type=int,
    default=1,
    show_default=True,
    help="Minimum samples/tissues supporting a junction. When >1, --min-reads applies per-sample.",
)
# Boundary adjustment options
@click.option(
    "--adjust-boundaries/--no-adjust-boundaries",
    default=True,
    show_default=True,
    help="Adjust start/stop codon boundaries.",
)
@click.option(
    "--boundary-window",
    type=int,
    default=30,
    show_default=True,
    help="Search window for boundary adjustment (bp).",
)
# Evidence scoring options
@click.option(
    "--min-coverage",
    type=int,
    default=5,
    show_default=True,
    help="Minimum coverage to count exon as expressed.",
)
@click.option(
    "--boundary-tolerance",
    type=int,
    default=10,
    show_default=True,
    help="Maximum bp shift for near-match junction comparison.",
)
@click.option(
    "--no-coverage",
    is_flag=True,
    help="Skip exon coverage analysis (junction-only evidence).",
)
# Confidence options
@click.option(
    "--confidence-threshold",
    type=float,
    default=0.5,
    show_default=True,
    help="Minimum confidence to flag LOW_CONF.",
)
# Execution options
@click.option(
    "-j",
    "--workers",
    type=int,
    default=1,
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
    "--scaffold",
    type=str,
    default=None,
    help="Process only this scaffold.",
)
@click.option(
    "--chunk-id",
    type=str,
    default=None,
    help="Chunk identifier for logging in parallel mode.",
)
@click.pass_context
def refine(
    ctx: click.Context,
    helixer_h5: Path,
    helixer_gff: Path,
    genome: Path,
    rnaseq_bam: tuple[str, ...],
    rnaseq_bam_list: Optional[Path],
    junctions_bed: tuple[str, ...],
    junctions_list: Optional[Path],
    output: Path,
    report: Optional[Path],
    splice_details: Optional[Path],
    evidence_details: Optional[Path],
    unsupported_bed: Optional[Path],
    max_shift: int,
    min_reads: int,
    min_tissues: int,
    adjust_boundaries: bool,
    boundary_window: int,
    min_coverage: int,
    boundary_tolerance: int,
    no_coverage: bool,
    confidence_threshold: float,
    workers: int,
    verbose_flag: bool,
    region: Optional[str],
    scaffold: Optional[str],
    chunk_id: Optional[str],
) -> None:
    """Refine Helixer predictions using RNA-seq evidence and confidence scoring.

    This is the main HelixForge pipeline that combines splice site correction,
    boundary adjustment, confidence scoring, and evidence scoring into a single
    workflow. RNA-seq evidence is REQUIRED for refinement.

    \b
    Steps performed:
    1. Correct splice sites using RNA-seq junctions
    2. Adjust start/stop codon boundaries (optional)
    3. Calculate confidence scores from Helixer HDF5
    4. Score evidence support from RNA-seq
    5. Write refined GFF3 with all scores as attributes

    \b
    Output GFF3 attributes:
    - confidence: Gene confidence score (0-1)
    - evidence_level: full/partial/minimal/none
    - aed: Annotation Edit Distance (0-1, MAKER-compatible)
    - junction_support: Fraction of supported junctions
    - splice_corrections: Number of corrected splice sites

    \b
    Examples:
        # Basic refinement with single BAM
        $ helixforge refine -p helixer.h5 -g helixer.gff3 --genome genome.fa \\
            --rnaseq-bam rnaseq.bam -o refined.gff3

        # Multi-tissue with reports
        $ helixforge refine -p helixer.h5 -g helixer.gff3 --genome genome.fa \\
            --rnaseq-bam liver.bam,brain.bam --min-tissues 2 \\
            -o refined.gff3 -r refine_report.tsv

        # Using STAR junction files
        $ helixforge refine -p helixer.h5 -g helixer.gff3 --genome genome.fa \\
            --junctions-bed sample1_SJ.out.tab,sample2_SJ.out.tab \\
            -o refined.gff3
    """
    from helixforge.core.refine import RefineConfig, RefinePipeline, RefineReportWriter
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser, GFF3Writer
    from helixforge.io.hdf5 import HelixerHDF5Reader
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

    # Collect BAM paths
    bam_paths: list[Path] = []
    for bam_arg in rnaseq_bam:
        for bam_str in bam_arg.split(","):
            bam_str = bam_str.strip()
            if bam_str:
                bam_path = Path(bam_str)
                if not bam_path.exists():
                    console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
                    raise SystemExit(1)
                bam_paths.append(bam_path)

    if rnaseq_bam_list:
        with open(rnaseq_bam_list) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    bam_path = Path(line)
                    if not bam_path.exists():
                        console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
                        raise SystemExit(1)
                    bam_paths.append(bam_path)

    # Collect junction file paths
    junction_paths: list[Path] = []
    for junc_arg in junctions_bed:
        for junc_str in junc_arg.split(","):
            junc_str = junc_str.strip()
            if junc_str:
                junc_path = Path(junc_str)
                if not junc_path.exists():
                    console.print(f"[red]Error:[/red] Junction file not found: {junc_path}")
                    raise SystemExit(1)
                junction_paths.append(junc_path)

    if junctions_list:
        with open(junctions_list) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    junc_path = Path(line)
                    if not junc_path.exists():
                        console.print(f"[red]Error:[/red] Junction file not found: {junc_path}")
                        raise SystemExit(1)
                    junction_paths.append(junc_path)

    # Validate - require at least one RNA-seq source
    if not bam_paths and not junction_paths:
        console.print(
            "[red]Error:[/red] RNA-seq evidence is required for refinement.\n"
            "Provide --rnaseq-bam, --rnaseq-bam-list, --junctions-bed, or --junctions-list"
        )
        raise SystemExit(1)

    # Validate min-tissues
    n_samples = len(bam_paths) if bam_paths else len(junction_paths)
    if min_tissues > 1 and n_samples < 2:
        console.print(
            "[red]Error:[/red] --min-tissues > 1 requires multiple input files"
        )
        raise SystemExit(1)

    if not quiet:
        console.print(f"[blue]Helixer H5:[/blue] {helixer_h5}")
        console.print(f"[blue]Helixer GFF:[/blue] {helixer_gff}")
        console.print(f"[blue]Genome:[/blue] {genome}")
        if bam_paths:
            console.print(f"[blue]RNA-seq BAMs:[/blue] {len(bam_paths)} file(s)")
            for bp in bam_paths:
                console.print(f"  - {bp.name}")
        if junction_paths:
            console.print(f"[blue]Junction files:[/blue] {len(junction_paths)} file(s)")
            for jp in junction_paths:
                console.print(f"  - {jp.name}")
        if min_tissues > 1:
            console.print(f"[blue]Min tissues:[/blue] {min_tissues}")
        console.print(f"[blue]Output:[/blue] {output}")

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
            if scaffold not in genome_accessor.scaffold_lengths:
                console.print(f"[red]Error:[/red] Scaffold '{scaffold}' not found in genome")
                genome_accessor.close()
                raise SystemExit(1)
            scaffold_len = genome_accessor.scaffold_lengths[scaffold]
            target_region = region_from_scaffold(scaffold, scaffold_len)
            if not quiet:
                console.print(f"[blue]Processing scaffold:[/blue] {scaffold} ({scaffold_len:,} bp)")

        # Load HDF5 predictions
        if not quiet:
            console.print("[dim]Loading HDF5 predictions...[/dim]")
        hdf5_reader = HelixerHDF5Reader(helixer_h5)

        # Load genes from GFF
        if not quiet:
            console.print("[dim]Loading gene predictions...[/dim]")
        parser = GFF3Parser(helixer_gff)
        all_genes = list(parser.iter_genes())

        # Filter to region if specified
        if target_region:
            genes = [
                g for g in all_genes
                if g.seqid == target_region.seqid
                and g.start >= target_region.start
                and g.end <= target_region.end
            ]
        else:
            genes = all_genes

        if not quiet:
            console.print(f"[dim]Loaded {len(genes)} genes[/dim]")

        # Configure pipeline
        config = RefineConfig(
            max_shift=max_shift,
            min_junction_reads=min_reads,
            min_tissues=min_tissues,
            adjust_boundaries=adjust_boundaries,
            boundary_search_window=boundary_window,
            confidence_threshold=confidence_threshold,
            min_exon_coverage=min_coverage,
            boundary_tolerance=boundary_tolerance,
            skip_coverage=no_coverage,
        )

        # Create pipeline
        if not quiet:
            console.print("[dim]Initializing refinement pipeline...[/dim]")

        pipeline = RefinePipeline(
            genome=genome_accessor,
            hdf5_reader=hdf5_reader,
            bam_files=bam_paths if bam_paths else None,
            junction_files=junction_paths if junction_paths else None,
            config=config,
            verbose=verbose,
            quiet=quiet,
        )

        # Refine genes
        if not quiet:
            console.print("[dim]Refining genes...[/dim]")

        refined_genes = []
        for gene in genes:
            refined = pipeline.refine_gene(gene)
            refined_genes.append(refined)

        # Write output GFF
        if not quiet:
            console.print("[dim]Writing refined GFF...[/dim]")

        writer = GFF3Writer(output)
        for refined in refined_genes:
            # Add additional attributes not already set by the pipeline
            # (pipeline already adds: confidence_score, evidence_score, aed,
            #  junction_support, mean_coverage, flags)
            attrs = refined.gene.attributes.copy()
            if refined.evidence:
                attrs["evidence_level"] = refined.evidence.evidence_level.value
            if refined.splice_corrections > 0:
                attrs["splice_corrections"] = str(refined.splice_corrections)
            if refined.boundary_adjusted:
                attrs["boundary_adjusted"] = "yes"

            # Write gene with updated attributes
            refined.gene.attributes = attrs
            writer.write_gene(refined.gene)
        writer.close()

        # Write reports
        if report:
            report_writer = RefineReportWriter()
            report_writer.write_summary_tsv(refined_genes, report)
            if not quiet:
                console.print(f"[green]Wrote refine report:[/green] {report}")

        if splice_details:
            report_writer = RefineReportWriter()
            report_writer.write_splice_details_tsv(refined_genes, splice_details)
            if not quiet:
                console.print(f"[green]Wrote splice details:[/green] {splice_details}")

        if evidence_details:
            from helixforge.core.evidence_output import write_evidence_report_tsv
            scores = [r.evidence for r in refined_genes if r.evidence is not None]
            write_evidence_report_tsv(scores, evidence_details)
            if not quiet:
                console.print(f"[green]Wrote evidence details:[/green] {evidence_details}")

        if unsupported_bed:
            report_writer = RefineReportWriter()
            report_writer.write_unsupported_introns_bed(refined_genes, unsupported_bed)
            if not quiet:
                console.print(f"[green]Wrote unsupported introns:[/green] {unsupported_bed}")

        # Clean up
        hdf5_reader.close()
        genome_accessor.close()

        # Print summary
        if not quiet:
            n_corrected = sum(1 for r in refined_genes if r.splice_corrections > 0)
            n_boundary = sum(1 for r in refined_genes if r.boundary_adjusted)
            n_with_evidence = sum(1 for r in refined_genes if r.evidence is not None)

            mean_conf = 0.0
            if refined_genes:
                confs = [r.confidence.overall_score for r in refined_genes if r.confidence]
                if confs:
                    mean_conf = sum(confs) / len(confs)

            mean_aed = 0.0
            if n_with_evidence > 0:
                aeds = [r.evidence.aed for r in refined_genes if r.evidence]
                mean_aed = sum(aeds) / len(aeds) if aeds else 0.0

            console.print("")
            console.print("[bold]Refinement Summary:[/bold]")
            console.print(f"  Total genes:           {len(refined_genes):,}")
            console.print(f"  Splice corrections:    {n_corrected:,} genes modified")
            console.print(f"  Boundary adjustments:  {n_boundary:,} genes")
            console.print(f"  Mean confidence:       {mean_conf:.4f}")
            console.print(f"  Mean AED:              {mean_aed:.4f}")
            console.print("")
            console.print(f"[green]Wrote refined GFF:[/green] {output}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


# =============================================================================
# evidence command
# =============================================================================


@main.command("evidence")
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
    type=str,
    multiple=True,
    help="RNA-seq BAM file(s). Comma-separated or specify multiple times.",
)
@click.option(
    "--bam-list",
    type=click.Path(exists=True, path_type=Path),
    help="File containing BAM paths (one per line).",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output GFF3 file with evidence annotations.",
)
@click.option(
    "--report",
    "-r",
    type=click.Path(path_type=Path),
    help="Output detailed evidence report TSV.",
)
@click.option(
    "--summary",
    type=click.Path(path_type=Path),
    help="Output evidence summary TSV for QC aggregation.",
)
@click.option(
    "--junction-details",
    type=click.Path(path_type=Path),
    help="Output per-junction details TSV.",
)
@click.option(
    "--exon-details",
    type=click.Path(path_type=Path),
    help="Output per-exon details TSV.",
)
@click.option(
    "--min-reads",
    type=int,
    default=3,
    show_default=True,
    help="Minimum reads to support a splice junction.",
)
@click.option(
    "--min-coverage",
    type=int,
    default=5,
    show_default=True,
    help="Minimum coverage to count exon as expressed.",
)
@click.option(
    "--boundary-tolerance",
    type=int,
    default=10,
    show_default=True,
    help="Maximum bp shift for near-match junctions.",
)
@click.option(
    "--no-coverage",
    is_flag=True,
    help="Skip exon coverage analysis (junction-only mode).",
)
@click.pass_context
def add_evidence(
    ctx: click.Context,
    gff: Path,
    bam: tuple[str, ...],
    bam_list: Optional[Path],
    output: Path,
    report: Optional[Path],
    summary: Optional[Path],
    junction_details: Optional[Path],
    exon_details: Optional[Path],
    min_reads: int,
    min_coverage: int,
    boundary_tolerance: int,
    no_coverage: bool,
) -> None:
    """Score gene predictions with RNA-seq evidence.

    Calculates evidence support for each gene based on:
    - Splice junction support from RNA-seq alignments
    - Exon coverage profiles (unless --no-coverage)
    - Boundary agreement at start/stop codons

    Outputs an Annotation Edit Distance (AED) score for each gene,
    where 0 = perfect support and 1 = no support.

    This command performs SCORING ONLY - it does not modify splice sites
    or boundaries. Use 'helixforge refine' for full refinement.

    \b
    Output attributes added to GFF3:
    - evidence_level: full/partial/minimal/none
    - aed: Annotation Edit Distance (0-1)
    - junction_support: fraction of supported junctions
    - exon_coverage: fraction of expressed exons
    - evidence_flags: any warning flags

    \b
    Examples:
        # Basic usage
        $ helixforge evidence -g predictions.gff3 -b rnaseq.bam -o annotated.gff3

        # Multiple BAMs
        $ helixforge evidence -g predictions.gff3 -b tissue1.bam,tissue2.bam -o annotated.gff3

        # With detailed reports
        $ helixforge evidence -g predictions.gff3 -b rnaseq.bam -o annotated.gff3 \\
            --report evidence_report.tsv --summary evidence_summary.tsv
    """
    from helixforge.core.evidence import EvidenceScorer, EvidenceScorerConfig
    from helixforge.core.evidence_output import (
        update_gff_with_evidence,
        write_evidence_report_tsv,
        write_exon_details_tsv,
        write_evidence_summary_txt,
        write_gene_evidence_summary_tsv,
        write_junction_details_tsv,
    )
    from helixforge.io.bam import JunctionExtractor
    from helixforge.io.gff import GFF3Parser

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    # Collect BAM files
    bam_paths: list[Path] = []

    # From --bam (may be comma-separated)
    for bam_arg in bam:
        for path_str in bam_arg.split(","):
            path_str = path_str.strip()
            if path_str:
                bam_path = Path(path_str)
                if not bam_path.exists():
                    console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
                    raise SystemExit(1)
                bam_paths.append(bam_path)

    # From --bam-list
    if bam_list:
        with open(bam_list) as f:
            for line in f:
                path_str = line.strip()
                if path_str and not path_str.startswith("#"):
                    bam_path = Path(path_str)
                    if not bam_path.exists():
                        console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
                        raise SystemExit(1)
                    bam_paths.append(bam_path)

    if not bam_paths:
        console.print("[red]Error:[/red] At least one BAM file is required")
        raise SystemExit(1)

    if not quiet:
        console.print(f"[blue]Input GFF:[/blue] {gff}")
        console.print(f"[blue]BAM files:[/blue] {len(bam_paths)}")
        for bam_path in bam_paths:
            console.print(f"  - {bam_path}")
        console.print(f"[blue]Output GFF:[/blue] {output}")

    try:
        # Parse GFF
        if not quiet:
            console.print("[dim]Loading gene predictions...[/dim]")
        parser = GFF3Parser(gff)
        genes = list(parser.iter_genes())
        if not quiet:
            console.print(f"[dim]Loaded {len(genes)} genes[/dim]")

        # Extract junctions from all BAM files
        if not quiet:
            console.print("[dim]Extracting splice junctions...[/dim]")

        all_junctions: dict[str, list] = {}
        primary_extractor = None

        for bam_path in bam_paths:
            extractor = JunctionExtractor(bam_path, min_reads=1)
            junctions = extractor.extract_all(min_reads=1)

            # Merge junctions
            for seqid, juncs in junctions.items():
                if seqid not in all_junctions:
                    all_junctions[seqid] = []
                all_junctions[seqid].extend(juncs)

            # Keep first extractor for coverage analysis
            if primary_extractor is None and not no_coverage:
                primary_extractor = extractor
            else:
                extractor.close()

        total_junctions = sum(len(j) for j in all_junctions.values())
        if not quiet:
            console.print(f"[dim]Extracted {total_junctions:,} junctions[/dim]")

        # Configure scorer
        config = EvidenceScorerConfig(
            min_junction_reads=min_reads,
            min_exon_coverage=min_coverage,
            boundary_tolerance=boundary_tolerance,
        )
        scorer = EvidenceScorer(config)

        # Score genes
        if not quiet:
            console.print("[dim]Scoring genes...[/dim]")

        scores = list(scorer.score_genes(
            genes,
            all_junctions,
            extractor=primary_extractor if not no_coverage else None,
        ))

        # Close extractor if still open
        if primary_extractor is not None:
            primary_extractor.close()

        # Write outputs
        if not quiet:
            console.print("[dim]Writing outputs...[/dim]")

        # Main GFF output
        update_gff_with_evidence(genes, scores, output)

        # Optional reports
        if report:
            write_evidence_report_tsv(scores, report)
            if not quiet:
                console.print(f"[green]Wrote detailed report:[/green] {report}")

        if summary:
            write_gene_evidence_summary_tsv(scores, summary)
            if not quiet:
                console.print(f"[green]Wrote summary:[/green] {summary}")

        if junction_details:
            write_junction_details_tsv(scores, junction_details)
            if not quiet:
                console.print(f"[green]Wrote junction details:[/green] {junction_details}")

        if exon_details:
            write_exon_details_tsv(scores, exon_details)
            if not quiet:
                console.print(f"[green]Wrote exon details:[/green] {exon_details}")

        # Print summary
        if not quiet:
            from helixforge.core.evidence import EvidenceLevel, summarize_evidence_scores

            stats = summarize_evidence_scores(scores)
            console.print("")
            console.print("[bold]Evidence Summary:[/bold]")
            console.print(f"  Total genes:         {stats['n_genes']:,}")
            console.print(f"  Mean AED:            {stats['mean_aed']:.4f}")
            console.print(f"  Full support:        {stats['n_full_support']:,} ({stats['n_full_support']/stats['n_genes']*100:.1f}%)")
            console.print(f"  Partial support:     {stats['n_partial_support']:,} ({stats['n_partial_support']/stats['n_genes']*100:.1f}%)")
            console.print(f"  Minimal support:     {stats['n_minimal_support']:,} ({stats['n_minimal_support']/stats['n_genes']*100:.1f}%)")
            console.print(f"  No support:          {stats['n_no_support']:,} ({stats['n_no_support']/stats['n_genes']*100:.1f}%)")
            console.print("")
            console.print(f"[green]Wrote annotated GFF:[/green] {output}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


# =============================================================================
# qc command group
# =============================================================================


@main.group()
def qc():
    """Quality control and reporting.

    Commands for aggregating module results, generating QC reports,
    filtering genes by criteria, and creating tiered output sets.

    Recommended workflow:
        1. Aggregate results: helixforge qc aggregate
        2. Generate report: helixforge qc report
        3. Filter genes: helixforge qc filter
        4. Create tiered outputs: helixforge qc tiered-output
    """
    pass


@qc.command("aggregate")
@click.option(
    "--refine-tsv",
    type=click.Path(exists=True, path_type=Path),
    help="Refine report TSV from 'helixforge refine'. Contains confidence, splice, and evidence scores.",
)
@click.option(
    "--confidence-tsv",
    type=click.Path(exists=True, path_type=Path),
    help="Confidence scores TSV from 'helixforge confidence'. Alternative if not using --refine-tsv.",
)
@click.option(
    "--splice-tsv",
    type=click.Path(exists=True, path_type=Path),
    help="Splice report TSV from 'helixforge splice' (deprecated). Use --refine-tsv instead.",
)
@click.option(
    "--homology-tsv",
    type=click.Path(exists=True, path_type=Path),
    help="Homology validation TSV from 'helixforge homology validate'.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output aggregated QC TSV file.",
)
@click.option(
    "--high-threshold",
    type=float,
    default=0.85,
    show_default=True,
    help="Minimum confidence for high tier.",
)
@click.option(
    "--medium-threshold",
    type=float,
    default=0.70,
    show_default=True,
    help="Minimum confidence for medium tier.",
)
@click.option(
    "--low-threshold",
    type=float,
    default=0.50,
    show_default=True,
    help="Minimum confidence for low tier (below = reject).",
)
@click.pass_context
def qc_aggregate(
    ctx: click.Context,
    refine_tsv: Optional[Path],
    confidence_tsv: Optional[Path],
    splice_tsv: Optional[Path],
    homology_tsv: Optional[Path],
    output: Path,
    high_threshold: float,
    medium_threshold: float,
    low_threshold: float,
) -> None:
    """Aggregate QC results from refine, confidence, splice, and homology modules.

    Combines results from individual analysis modules into unified GeneQC
    objects with tier classifications and appropriate flags.

    The preferred workflow is to use --refine-tsv from 'helixforge refine' output,
    which includes confidence, splice, and evidence scores in a single file.

    \b
    Examples:
        # Using refine output (recommended)
        $ helixforge qc aggregate --refine-tsv refine_report.tsv --homology-tsv validation.tsv -o qc_results.tsv

        # Using separate module outputs
        $ helixforge qc aggregate --confidence-tsv scores.tsv --splice-tsv splice_report.tsv -o qc_results.tsv
    """
    from helixforge.qc import QCAggregator, QCAggregatorConfig, export_qc_tsv

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not refine_tsv and not confidence_tsv and not splice_tsv and not homology_tsv:
        console.print("[red]Error:[/red] At least one input TSV is required")
        raise SystemExit(1)

    if not quiet:
        if refine_tsv:
            console.print(f"[blue]Refine TSV:[/blue] {refine_tsv}")
        if confidence_tsv:
            console.print(f"[blue]Confidence TSV:[/blue] {confidence_tsv}")
        if splice_tsv:
            console.print(f"[blue]Splice TSV:[/blue] {splice_tsv}")
        if homology_tsv:
            console.print(f"[blue]Homology TSV:[/blue] {homology_tsv}")

    try:
        config = QCAggregatorConfig(
            high_confidence_threshold=high_threshold,
            medium_confidence_threshold=medium_threshold,
            low_confidence_threshold=low_threshold,
        )

        aggregator = QCAggregator(config)
        gene_qcs = aggregator.aggregate_from_files(
            refine_tsv=refine_tsv,
            confidence_tsv=confidence_tsv,
            splice_tsv=splice_tsv,
            homology_tsv=homology_tsv,
        )

        export_qc_tsv(gene_qcs, output)

        if not quiet:
            console.print(f"[green]Aggregated {len(gene_qcs)} genes to:[/green] {output}")

            # Print tier summary
            from helixforge.qc import summarize_qc_results
            summary = summarize_qc_results(gene_qcs)
            tier_counts = summary.get("tier_counts", {})
            console.print("\n[bold]Tier Distribution:[/bold]")
            for tier in ["high", "medium", "low", "reject"]:
                count = tier_counts.get(tier, 0)
                pct = count / len(gene_qcs) * 100 if gene_qcs else 0
                console.print(f"  {tier.capitalize()}: {count} ({pct:.1f}%)")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@qc.command("report")
@click.option(
    "--qc-tsv",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Aggregated QC TSV from 'helixforge qc aggregate'.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output report file (HTML, TXT, or JSON based on extension).",
)
@click.option(
    "--title",
    type=str,
    default="HelixForge QC Report",
    show_default=True,
    help="Report title.",
)
@click.option(
    "--description",
    type=str,
    default="",
    help="Report description.",
)
@click.pass_context
def qc_report(
    ctx: click.Context,
    qc_tsv: Path,
    output: Path,
    title: str,
    description: str,
) -> None:
    """Generate QC report from aggregated results.

    Creates comprehensive reports in HTML (with Chart.js visualizations),
    plain text, or JSON formats.

    Example:
        $ helixforge qc report --qc-tsv qc_results.tsv -o report.html
        $ helixforge qc report --qc-tsv qc_results.tsv -o summary.txt
    """
    from helixforge.qc import (
        GeneQC,
        QCReportGenerator,
        generate_json_report,
        generate_summary_report,
    )

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]Loading QC data from:[/blue] {qc_tsv}")

    try:
        # Load QC data
        gene_qcs = _load_qc_tsv(qc_tsv)

        if not quiet:
            console.print(f"[green]Loaded {len(gene_qcs)} genes[/green]")

        # Generate report based on extension
        suffix = output.suffix.lower()

        if suffix == ".html":
            generator = QCReportGenerator()
            generator.generate(gene_qcs, output, title=title, description=description)
        elif suffix == ".txt":
            generate_summary_report(gene_qcs, output)
        elif suffix == ".json":
            generate_json_report(gene_qcs, output)
        else:
            console.print(f"[red]Error:[/red] Unsupported format: {suffix}")
            console.print("  Supported: .html, .txt, .json")
            raise SystemExit(1)

        if not quiet:
            console.print(f"[green]Wrote report to:[/green] {output}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@qc.command("filter")
@click.option(
    "--qc-tsv",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Aggregated QC TSV from 'helixforge qc aggregate'.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    required=True,
    help="Output filtered gene list (TXT) or filtered QC TSV.",
)
@click.option(
    "--preset",
    type=click.Choice([
        "high_confidence",
        "publication_ready",
        "has_homology",
        "has_rnaseq_support",
        "no_structural_issues",
    ]),
    help="Use a preset filter profile.",
)
@click.option(
    "--tiers",
    type=str,
    help="Comma-separated list of tiers to include (e.g., 'high,medium').",
)
@click.option(
    "--min-confidence",
    type=float,
    help="Minimum confidence score.",
)
@click.option(
    "--exclude-flags",
    type=str,
    help="Comma-separated flag codes to exclude.",
)
@click.option(
    "--gene-list-only",
    is_flag=True,
    help="Output gene IDs only (one per line).",
)
@click.pass_context
def qc_filter(
    ctx: click.Context,
    qc_tsv: Path,
    output: Path,
    preset: Optional[str],
    tiers: Optional[str],
    min_confidence: Optional[float],
    exclude_flags: Optional[str],
    gene_list_only: bool,
) -> None:
    """Filter genes based on QC criteria.

    Apply preset filter profiles or custom criteria to select genes.

    \b
    Preset profiles:
    - high_confidence: High tier, confidence >= 0.85, no warnings
    - publication_ready: High/medium tier, no errors, no critical issues
    - has_homology: Genes with protein homology support
    - has_rnaseq_support: Genes with RNA-seq splice junction support
    - no_structural_issues: No internal stops, frameshifts, or missing codons

    Example:
        $ helixforge qc filter --qc-tsv qc_results.tsv --preset publication_ready -o filtered.txt --gene-list-only
        $ helixforge qc filter --qc-tsv qc_results.tsv --tiers high,medium --min-confidence 0.7 -o filtered_qc.tsv
    """
    from helixforge.qc import (
        FilterCriteria,
        Flags,
        GeneFilter,
        export_qc_tsv,
        summarize_filter_results,
    )

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]Loading QC data from:[/blue] {qc_tsv}")

    try:
        gene_qcs = _load_qc_tsv(qc_tsv)

        if not quiet:
            console.print(f"[green]Loaded {len(gene_qcs)} genes[/green]")

        # Create filter based on options
        if preset:
            if preset == "high_confidence":
                gene_filter = GeneFilter.high_confidence()
            elif preset == "publication_ready":
                gene_filter = GeneFilter.publication_ready()
            elif preset == "has_homology":
                gene_filter = GeneFilter.has_homology()
            elif preset == "has_rnaseq_support":
                gene_filter = GeneFilter.has_rnaseq_support()
            elif preset == "no_structural_issues":
                gene_filter = GeneFilter.no_structural_issues()
            else:
                raise ValueError(f"Unknown preset: {preset}")

            if not quiet:
                console.print(f"[blue]Using preset:[/blue] {preset}")

        else:
            # Build custom criteria
            allowed_tiers = None
            if tiers:
                allowed_tiers = [t.strip() for t in tiers.split(",")]

            exclude_flag_list = []
            if exclude_flags:
                for code in exclude_flags.split(","):
                    flag = Flags.get_by_code(code.strip())
                    if flag:
                        exclude_flag_list.append(flag)
                    else:
                        console.print(f"[yellow]Warning:[/yellow] Unknown flag code: {code.strip()}")

            criteria = FilterCriteria(
                name="custom",
                description="Custom filter from CLI",
                min_confidence=min_confidence,
                allowed_tiers=allowed_tiers,
                exclude_flags=exclude_flag_list,
            )
            gene_filter = GeneFilter(criteria)

        # Apply filter
        result = gene_filter.apply(gene_qcs)

        if not quiet:
            console.print(f"\n[bold]Filter Results:[/bold]")
            console.print(summarize_filter_results(result))

        # Write output
        if gene_list_only:
            with open(output, "w") as f:
                for gene_id in result.passed_ids():
                    f.write(f"{gene_id}\n")
        else:
            passed_qcs = {qc.gene_id: qc for qc in result.passed}
            export_qc_tsv(passed_qcs, output)

        if not quiet:
            console.print(f"\n[green]Wrote {result.pass_count} genes to:[/green] {output}")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@qc.command("tiered-output")
@click.option(
    "--qc-tsv",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Aggregated QC TSV from 'helixforge qc aggregate'.",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Output directory for tiered gene lists.",
)
@click.option(
    "--prefix",
    type=str,
    default="genes",
    show_default=True,
    help="Prefix for output file names.",
)
@click.option(
    "--input-gff",
    type=click.Path(exists=True, path_type=Path),
    help="Input GFF3 to create tiered GFF outputs.",
)
@click.pass_context
def qc_tiered_output(
    ctx: click.Context,
    qc_tsv: Path,
    output_dir: Path,
    prefix: str,
    input_gff: Optional[Path],
) -> None:
    """Generate tiered gene lists and optionally GFF files.

    Creates separate files for each quality tier (high, medium, low, reject),
    and optionally creates filtered GFF files.

    Example:
        $ helixforge qc tiered-output --qc-tsv qc_results.tsv -o tiered_output/
        $ helixforge qc tiered-output --qc-tsv qc_results.tsv -o tiered_output/ --input-gff genes.gff3
    """
    from helixforge.qc import write_filtered_gff, write_tiered_gene_lists

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    if not quiet:
        console.print(f"[blue]Loading QC data from:[/blue] {qc_tsv}")

    try:
        gene_qcs = _load_qc_tsv(qc_tsv)

        if not quiet:
            console.print(f"[green]Loaded {len(gene_qcs)} genes[/green]")

        # Write tiered gene lists
        output_files = write_tiered_gene_lists(gene_qcs, output_dir, prefix)

        if not quiet:
            console.print("\n[bold]Created tiered gene lists:[/bold]")
            for tier, path in sorted(output_files.items()):
                n_genes = sum(1 for _ in open(path))
                console.print(f"  {tier}: {path} ({n_genes} genes)")

        # Write tiered GFFs if requested
        if input_gff:
            for tier in ["high", "medium", "low"]:
                gff_path = output_dir / f"{prefix}_{tier}.gff3"
                n_written = write_filtered_gff(
                    gene_qcs, input_gff, gff_path, tiers=[tier]
                )
                if not quiet:
                    console.print(f"  {tier} GFF: {gff_path} ({n_written} genes)")

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        raise SystemExit(1)


@qc.command("list-flags")
@click.option(
    "--category",
    type=click.Choice(["confidence", "splice", "homology", "structure", "annotation"]),
    help="Filter flags by category.",
)
@click.option(
    "--severity",
    type=click.Choice(["info", "warning", "error", "critical"]),
    help="Filter flags by severity.",
)
def qc_list_flags(
    category: Optional[str],
    severity: Optional[str],
) -> None:
    """List all available QC flags.

    Shows flag codes, names, descriptions, categories, and severity levels.

    Example:
        $ helixforge qc list-flags
        $ helixforge qc list-flags --category homology
        $ helixforge qc list-flags --severity error
    """
    from helixforge.qc import FlagCategory, FlagSeverity, Flags

    flags = Flags.get_all()

    # Filter if requested
    if category:
        cat = FlagCategory(category)
        flags = [f for f in flags if f.category == cat]

    if severity:
        sev = FlagSeverity(severity)
        flags = [f for f in flags if f.severity == sev]

    # Sort by category, then severity, then code
    flags = sorted(flags, key=lambda f: (f.category.value, f.severity.value, f.code))

    console.print("[bold]QC Flags:[/bold]\n")

    current_category = None
    for flag in flags:
        if flag.category != current_category:
            current_category = flag.category
            console.print(f"\n[bold blue]{current_category.value.upper()}[/bold blue]")

        severity_color = {
            "info": "cyan",
            "warning": "yellow",
            "error": "red",
            "critical": "bold red",
        }.get(flag.severity.value, "white")

        console.print(f"  [{severity_color}]{flag.code}[/{severity_color}]")
        console.print(f"    Name: {flag.name}")
        console.print(f"    Description: {flag.description}")
        console.print(f"    Severity: {flag.severity.value}")
        console.print()


def _load_qc_tsv(path: Path) -> dict:
    """Load GeneQC objects from aggregated TSV file."""
    import csv
    from helixforge.qc import Flags, GeneQC

    gene_qcs = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row.get("gene_id", "")
            if not gene_id:
                continue

            gene_qc = GeneQC(
                gene_id=gene_id,
                tier=row.get("tier", "unclassified"),
            )

            # Parse scores
            if row.get("confidence_score"):
                try:
                    gene_qc.confidence_score = float(row["confidence_score"])
                except ValueError:
                    pass

            if row.get("splice_score"):
                try:
                    gene_qc.splice_score = float(row["splice_score"])
                except ValueError:
                    pass

            if row.get("homology_score"):
                try:
                    gene_qc.homology_score = float(row["homology_score"])
                except ValueError:
                    pass

            # Parse flags
            flag_codes = row.get("flag_codes", "")
            if flag_codes:
                for code in flag_codes.split(","):
                    code = code.strip()
                    flag = Flags.get_by_code(code)
                    if flag:
                        gene_qc.add_flag(flag)

            gene_qcs[gene_id] = gene_qc

    return gene_qcs


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
    required=False,
    default=None,
    help="Reference genome FASTA file. Optional if --input-h5 is provided or auto-detected.",
)
@click.option(
    "--input-h5",
    "input_h5",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    default=None,
    help="Helixer input HDF5 file (contains coordinate mapping). "
    "If not provided, will auto-detect from predictions file location. "
    "Required for strand-aware prediction retrieval.",
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
    genome: Optional[Path],
    input_h5: Optional[Path],
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

    The command requires coordinate mapping information, which can come from:
    \b
    - --input-h5: Helixer input HDF5 file (preferred, enables strand-aware mapping)
    - --genome: Reference FASTA with .fai index (fallback)
    - Auto-detection: Will look for *_input.h5 alongside *_predictions.h5

    For parallel execution, use --region or --scaffold to process a subset
    of genes, and --chunk-id for logging and output organization.

    Example:
        $ helixforge confidence -p predictions.h5 -g genes.gff3 --input-h5 input.h5 -o scores.tsv
        $ helixforge confidence -p predictions.h5 -g genes.gff3 -o scores.tsv  # auto-detect input.h5
        $ helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa -o scores.tsv  # fallback
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
    )

    verbose = ctx.obj.get("verbose", False)
    quiet = ctx.obj.get("quiet", False)

    # Log chunk ID if provided
    if chunk_id and not quiet:
        console.print(f"[blue]Chunk ID:[/blue] {chunk_id}")

    # Determine coordinate mapping source
    fai_path = None
    if genome is not None:
        fai_path = genome.with_suffix(genome.suffix + ".fai")
        if not fai_path.exists():
            console.print(f"[red]Error:[/red] FAI index not found: {fai_path}")
            console.print("Run 'samtools faidx' to create the index.")
            raise SystemExit(1)

    if not quiet:
        console.print(f"[blue]Loading predictions from:[/blue] {predictions}")
        console.print(f"[blue]Loading genes from:[/blue] {gff}")
        if input_h5:
            console.print(f"[blue]Helixer input file:[/blue] {input_h5}")
        elif genome:
            console.print(f"[blue]Reference genome:[/blue] {genome}")
        else:
            console.print("[blue]Coordinate mapping:[/blue] auto-detect from predictions path")

    try:
        # Load HDF5 reader with appropriate coordinate mapping
        reader = HelixerHDF5Reader(
            predictions,
            fasta_index=fai_path,
            input_h5_path=input_h5,
        )

        # Get scaffold info from reader's coordinate index
        scaffold_lengths = reader.coord_index.scaffold_lengths

        if not quiet and reader.has_dual_strand:
            console.print("[green]Using strand-aware coordinate mapping[/green]")

        # Optionally load genome accessor for validation
        genome_accessor = None
        if genome is not None:
            genome_accessor = GenomeAccessor(genome)

        try:
            # Parse region constraints
            target_region: GenomicRegion | None = None

            if region:
                try:
                    target_region = parse_region(region)
                    # Validate against scaffold lengths from reader
                    if target_region.seqid not in scaffold_lengths:
                        raise ValueError(f"Scaffold '{target_region.seqid}' not found")
                    scaffold_len = scaffold_lengths[target_region.seqid]
                    if target_region.end > scaffold_len:
                        raise ValueError(
                            f"Region end ({target_region.end}) exceeds scaffold length ({scaffold_len})"
                        )
                    if not quiet:
                        console.print(f"[blue]Processing region:[/blue] {target_region}")
                except ValueError as e:
                    console.print(f"[red]Error:[/red] {e}")
                    raise SystemExit(1)
            elif scaffold:
                # Scaffold-only mode: process entire scaffold
                if scaffold not in scaffold_lengths:
                    console.print(
                        f"[red]Error:[/red] Scaffold '{scaffold}' not found. "
                        f"Available: {list(scaffold_lengths.keys())}"
                    )
                    raise SystemExit(1)
                scaffold_len = scaffold_lengths[scaffold]
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

                # Check if a file exists with the directory name
                if plot_dir.exists() and not plot_dir.is_dir():
                    console.print(
                        f"[red]Error:[/red] '{plot_dir}' exists but is not a directory. "
                        "Please remove it or use a different --plot-dir path."
                    )
                    raise SystemExit(1)
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

        finally:
            # Clean up resources
            if genome_accessor is not None:
                genome_accessor.close()
            reader.close()

    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise SystemExit(1)


# =============================================================================
# splice command (DEPRECATED - use refine instead)
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
    type=str,
    multiple=True,
    help="RNA-seq BAM file(s) (sorted, indexed). Accepts comma-separated list or can be specified multiple times.",
)
@click.option(
    "--rnaseq-bam-list",
    type=click.Path(exists=True, path_type=Path),
    help="File containing BAM paths (one per line).",
)
@click.option(
    "--junctions-bed",
    type=str,
    multiple=True,
    help="Junction file(s) in BED or STAR SJ.out.tab format. Comma-separated or repeated.",
)
@click.option(
    "--junctions-list",
    type=click.Path(exists=True, path_type=Path),
    help="File containing junction file paths (one per line). Supports BED and STAR SJ.out.tab.",
)
@click.option(
    "--min-tissues",
    default=1,
    type=int,
    show_default=True,
    help="Minimum number of samples/tissues supporting a junction. When >1, --min-reads applies per-sample.",
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
    rnaseq_bam: tuple[str, ...],
    rnaseq_bam_list: Optional[Path],
    junctions_bed: tuple[str, ...],
    junctions_list: Optional[Path],
    output_gff: Path,
    report: Optional[Path],
    corrections_detail: Optional[Path],
    unsupported_bed: Optional[Path],
    max_shift: int,
    min_reads: int,
    min_tissues: int,
    adjust_boundaries: bool,
    workers: int,
    verbose_flag: bool,
    region: Optional[str],
    chunk_id: Optional[str],
    scaffold: Optional[str],
) -> None:
    """[DEPRECATED] Refine splice sites using RNA-seq evidence.

    WARNING: This command is deprecated. Use 'helixforge refine' instead,
    which combines splice correction, boundary adjustment, confidence scoring,
    and evidence scoring in a single pipeline.

    This command refines splice sites in Helixer gene predictions using
    RNA-seq junction evidence and position weight matrix (PWM) scoring.

    \b
    Features:
    - Correct splice sites to match empirical junctions
    - Score candidate sites using PWMs for canonical splice motifs
    - Handle canonical (GT-AG) and non-canonical splice sites
    - Flag genes with no RNA-seq support or conflicting evidence
    - Optionally adjust start/stop codons
    - Multi-tissue/sample support for improved junction confidence

    \b
    Multi-BAM Mode:
    When multiple BAM files are provided (e.g., different tissues), junctions
    are aggregated across samples. Use --min-tissues to require junction
    support from multiple samples. In multi-BAM mode, --min-reads applies
    per-sample when --min-tissues > 1.

    For parallel execution, use --region or --scaffold to process a subset
    of genes, and --chunk-id for logging and output organization.

    \b
    Examples:
        # Single BAM
        $ helixforge splice --helixer-gff helixer.gff3 --genome genome.fa \\
            --rnaseq-bam rnaseq.bam -o refined.gff3 -r splice_report.tsv

        # Multiple BAMs (multi-tissue)
        $ helixforge splice --helixer-gff helixer.gff3 --genome genome.fa \\
            --rnaseq-bam liver.bam --rnaseq-bam brain.bam --rnaseq-bam heart.bam \\
            --min-tissues 2 --min-reads 3 -o refined.gff3 -r splice_report.tsv

        # Region-based parallel processing
        $ helixforge splice --helixer-gff helixer.gff3 --genome genome.fa \\
            --rnaseq-bam rnaseq.bam --region chr1:1-1000000 --chunk-id chunk_001 -o chunk_001.gff3
    """
    import warnings
    warnings.warn(
        "The 'splice' command is deprecated. Use 'helixforge refine' instead, "
        "which combines splice correction, confidence scoring, and evidence scoring.",
        DeprecationWarning,
        stacklevel=2,
    )
    console.print(
        "[yellow]Warning:[/yellow] The 'splice' command is deprecated. "
        "Use 'helixforge refine' for full refinement pipeline."
    )

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

    # Expand and validate BAM paths
    # Supports: repeated flags, comma-separated, and file list
    bam_paths: list[Path] = []

    # Process --rnaseq-bam arguments (may contain comma-separated values)
    for bam_arg in rnaseq_bam:
        # Split by comma, handling potential whitespace
        for bam_str in bam_arg.split(","):
            bam_str = bam_str.strip()
            if bam_str:
                bam_path = Path(bam_str)
                if not bam_path.exists():
                    console.print(f"[red]Error:[/red] BAM file not found: {bam_path}")
                    raise SystemExit(1)
                bam_paths.append(bam_path)

    # Process --rnaseq-bam-list file
    if rnaseq_bam_list:
        with open(rnaseq_bam_list) as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if line and not line.startswith("#"):
                    bam_path = Path(line)
                    if not bam_path.exists():
                        console.print(
                            f"[red]Error:[/red] BAM file not found: {bam_path} "
                            f"(from {rnaseq_bam_list})"
                        )
                        raise SystemExit(1)
                    bam_paths.append(bam_path)

    # Expand and validate junction file paths
    # Supports: repeated flags, comma-separated, and file list
    junction_paths: list[Path] = []

    # Process --junctions-bed arguments (may contain comma-separated values)
    for junc_arg in junctions_bed:
        for junc_str in junc_arg.split(","):
            junc_str = junc_str.strip()
            if junc_str:
                junc_path = Path(junc_str)
                if not junc_path.exists():
                    console.print(f"[red]Error:[/red] Junction file not found: {junc_path}")
                    raise SystemExit(1)
                junction_paths.append(junc_path)

    # Process --junctions-list file
    if junctions_list:
        with open(junctions_list) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    junc_path = Path(line)
                    if not junc_path.exists():
                        console.print(
                            f"[red]Error:[/red] Junction file not found: {junc_path} "
                            f"(from {junctions_list})"
                        )
                        raise SystemExit(1)
                    junction_paths.append(junc_path)

    # Validate inputs - need at least one source
    if not bam_paths and not junction_paths:
        console.print(
            "[red]Error:[/red] Provide --rnaseq-bam, --rnaseq-bam-list, "
            "--junctions-bed, or --junctions-list"
        )
        raise SystemExit(1)

    # Validate min-tissues (applies to whichever input type is provided)
    n_samples = len(bam_paths) if bam_paths else len(junction_paths)
    if min_tissues > 1 and n_samples < 2:
        console.print(
            "[red]Error:[/red] --min-tissues > 1 requires multiple input files"
        )
        raise SystemExit(1)

    if not quiet:
        console.print(f"[blue]Loading GFF from:[/blue] {helixer_gff}")
        console.print(f"[blue]Reference genome:[/blue] {genome}")
        if bam_paths:
            if len(bam_paths) == 1:
                console.print(f"[blue]RNA-seq BAM:[/blue] {bam_paths[0]}")
            else:
                console.print(f"[blue]RNA-seq BAMs:[/blue] {len(bam_paths)} files")
                for bam in bam_paths:
                    console.print(f"  - {bam.name}")
                if min_tissues > 1:
                    console.print(
                        f"[blue]Min tissues:[/blue] {min_tissues} "
                        f"(--min-reads={min_reads} applied per-sample)"
                    )
        if junction_paths:
            if len(junction_paths) == 1:
                console.print(f"[blue]Junctions file:[/blue] {junction_paths[0]}")
            else:
                console.print(f"[blue]Junction files:[/blue] {len(junction_paths)} files")
                for jp in junction_paths:
                    console.print(f"  - {jp.name}")
                if min_tissues > 1:
                    console.print(
                        f"[blue]Min tissues:[/blue] {min_tissues} "
                        f"(--min-reads={min_reads} applied per-sample)"
                    )

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
        if bam_paths:
            if len(bam_paths) == 1:
                # Single BAM mode - original behavior
                from helixforge.io.bam import JunctionExtractor

                if not quiet:
                    console.print("[blue]Extracting junctions from BAM...[/blue]")
                extractor = JunctionExtractor(bam_paths[0])

                if target_region:
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

                extractor.close()

                if not quiet:
                    n_junctions = sum(len(j) for j in junctions.values())
                    console.print(f"[green]Extracted {n_junctions} junctions[/green]")
            else:
                # Multi-BAM mode - aggregate across samples
                from helixforge.io.bam import (
                    aggregate_junctions_multi_sample,
                    multi_sample_to_standard_junctions,
                )

                if not quiet:
                    console.print(
                        f"[blue]Extracting and aggregating junctions from "
                        f"{len(bam_paths)} BAM files...[/blue]"
                    )

                # Prepare region tuple if needed
                region_tuple = None
                if target_region:
                    region_tuple = (
                        target_region.seqid,
                        target_region.start,
                        target_region.end,
                    )

                # When min_tissues > 1, min_reads applies per-sample
                # Otherwise, we use min_reads=1 per-sample and filter on total
                if min_tissues > 1:
                    min_reads_per_sample = min_reads
                    min_reads_total = 1  # Will be filtered by min_samples
                else:
                    min_reads_per_sample = 1
                    min_reads_total = min_reads

                multi_junctions = aggregate_junctions_multi_sample(
                    bam_paths=bam_paths,
                    min_reads_per_sample=min_reads_per_sample,
                    min_samples=min_tissues,
                    region=region_tuple,
                )

                # Convert to standard junctions for SpliceRefiner
                junctions = multi_sample_to_standard_junctions(multi_junctions)

                # Apply total read filter if min_tissues == 1
                if min_tissues == 1 and min_reads_total > 1:
                    junctions = {
                        seqid: [j for j in juncs if j.read_count >= min_reads_total]
                        for seqid, juncs in junctions.items()
                    }
                    junctions = {k: v for k, v in junctions.items() if v}

                if not quiet:
                    n_junctions = sum(len(j) for j in junctions.values())
                    n_samples = len(bam_paths)
                    console.print(
                        f"[green]Aggregated {n_junctions} junctions from "
                        f"{n_samples} samples[/green]"
                    )
                    if min_tissues > 1:
                        console.print(
                            f"[green]  (requiring >= {min_tissues} tissues with "
                            f">= {min_reads} reads each)[/green]"
                        )
        else:
            # Load from junction files (BED or STAR SJ.out.tab)
            from helixforge.io.bam import (
                aggregate_junctions_from_files,
                load_junctions_auto,
                multi_sample_to_standard_junctions,
            )

            if len(junction_paths) == 1:
                # Single file mode
                if not quiet:
                    console.print("[blue]Loading junctions from file...[/blue]")
                all_junctions = load_junctions_auto(junction_paths[0], min_reads=min_reads)
            else:
                # Multi-file mode - aggregate across samples
                if not quiet:
                    console.print(
                        f"[blue]Loading and aggregating junctions from "
                        f"{len(junction_paths)} files...[/blue]"
                    )

                if min_tissues > 1:
                    min_reads_per_sample = min_reads
                    min_reads_total = 1
                else:
                    min_reads_per_sample = 1
                    min_reads_total = min_reads

                multi_junctions = aggregate_junctions_from_files(
                    file_paths=junction_paths,
                    min_reads_per_sample=min_reads_per_sample,
                    min_samples=min_tissues,
                )
                all_junctions = multi_sample_to_standard_junctions(multi_junctions)

                # Apply total read filter if min_tissues == 1
                if min_tissues == 1 and min_reads_total > 1:
                    all_junctions = {
                        seqid: [j for j in juncs if j.read_count >= min_reads_total]
                        for seqid, juncs in all_junctions.items()
                    }
                    all_junctions = {k: v for k, v in all_junctions.items() if v}

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
                if len(junction_paths) == 1:
                    console.print(f"[green]Loaded {n_junctions} junctions[/green]")
                else:
                    console.print(
                        f"[green]Aggregated {n_junctions} junctions from "
                        f"{len(junction_paths)} files[/green]"
                    )
                    if min_tissues > 1:
                        console.print(
                            f"[green]  (requiring >= {min_tissues} samples with "
                            f">= {min_reads} reads each)[/green]"
                        )

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
        3. Execute with HyperShell: hs cluster tasks.txt --num-tasks 32
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
        # Simple task file (include all required options for your command)
        $ helixforge parallel tasks --chunk-plan chunks.json \\
            --command 'helixforge confidence -p predictions.h5 -g predictions.gff3 --genome genome.fa --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o {output_dir}/{chunk_id}.tsv' \\
            --output tasks.txt --output-dir outputs/

        # Execute with HyperShell
        $ hs cluster tasks.txt --num-tasks 32

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
