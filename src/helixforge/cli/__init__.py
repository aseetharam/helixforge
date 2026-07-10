"""HelixForge command-line interface."""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any, Callable

import click

from helixforge.bench.ablation import run_ablation
from helixforge.bench.wrappers import benchmark_all
from helixforge.reconcile.pipeline import PipelineConfig, run_pipeline
from helixforge.stats.before_after import before_after_table, write_summary
from helixforge.viz.interactive import interactive_index
from helixforge.viz.locus_plot import plot_loci, plot_locus
from helixforge.viz.tracks import write_bed12


def _console() -> Any:
    """A ``rich`` console (lazy import keeps the dependency at the entry point)."""
    from rich.console import Console

    return Console()


# ---------------------------------------------------------------------------
# Shared PipelineConfig options
# ---------------------------------------------------------------------------

# Each entry is a click.option whose dest name matches a PipelineConfig field
# (list fields are renamed in _build_config). Applied via @pipeline_options to
# the commands that build a full pipeline run (reconcile, viz, benchmark ablation).
_PIPELINE_OPTIONS: list[Callable[[Any], Any]] = [
    click.option(
        "--genome",
        "genome_fasta",
        required=True,
        type=click.Path(exists=True, dir_okay=False),
        help="Genome FASTA (required).",
    ),
    click.option(
        "--helixer",
        "helixer_gff3",
        required=True,
        type=click.Path(exists=True, dir_okay=False),
        help="Helixer GFF3 (required).",
    ),
    click.option(
        "--helixer-h5",
        "helixer_h5",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="Helixer confidence HDF5: a combined file, the '*_input.h5' "
        "metadata half, or the '*_predictions.h5' softmax half (the partner "
        "half is auto-detected as the sibling, or named via --helixer-input-h5).",
    ),
    click.option(
        "--helixer-input-h5",
        "helixer_input_h5",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="Helixer '*_input.h5' metadata half (coordinate mapping). Needed "
        "only when --helixer-h5 is a bare '*_predictions.h5' whose '*_input.h5' "
        "sibling is not co-located. Mirrors 'confidence --input-h5'.",
    ),
    # Multi-file evidence inputs: repeatable + comma-separated, each with a
    # ``-list`` FOFN companion (v1 convention; merged by expand_file_args in
    # _build_config). Direct flags are plain strings (NOT click.Path(exists)) so
    # a comma-joined value is not mis-validated as one path; existence is checked
    # during expansion instead.
    click.option(
        "--stringtie",
        "stringtie",
        multiple=True,
        type=str,
        help="StringTie GTF, one per sample. Comma-separated or repeated.",
    ),
    click.option(
        "--stringtie-list",
        "stringtie_list",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="File of StringTie GTF paths, one per line "
        "(blank lines and # comments ignored).",
    ),
    click.option(
        "--bam",
        "bam",
        multiple=True,
        type=str,
        help="RNA-seq BAM. Comma-separated or repeated.",
    ),
    click.option(
        "--bam-list",
        "bam_list",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="File of BAM paths, one per line (blank lines and # comments ignored).",
    ),
    click.option(
        "--bigwig",
        "bigwig",
        multiple=True,
        type=str,
        help="Coverage bigWig. Comma-separated or repeated.",
    ),
    click.option(
        "--bigwig-list",
        "bigwig_list",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="File of bigWig paths, one per line "
        "(blank lines and # comments ignored).",
    ),
    click.option(
        "--star-sj",
        "star_sj",
        multiple=True,
        type=str,
        help="STAR SJ.out.tab. Comma-separated or repeated.",
    ),
    click.option(
        "--star-sj-list",
        "star_sj_list",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="File of SJ.out.tab paths, one per line "
        "(blank lines and # comments ignored).",
    ),
    click.option(
        "--miniprot",
        "miniprot_gff",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="miniprot GFF for backstop CDS projection. "
        "Not a substitute for --protein-db; does not enable Mikado.",
    ),
    click.option(
        "--protein-db",
        "protein_db",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="Protein FASTA for DIAMOND homology. Required (with StringTie) "
        "to enable Mikado reconciliation and isoform discovery.",
    ),
    # Optional EDTA TE gating. EDTA is the only TE signal; off unless given.
    click.option(
        "--te-annotation",
        "te_annotation",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="EDTA TE GFF3. Flags TE-overlapping models (TE_OVERLAP) and "
        "reclassifies good-ORF genes above --te-overlap-threshold as "
        "transposable_element. Uses the EDTA Classification order; "
        "knob/satellite/low-complexity never gate.",
    ),
    click.option(
        "--te-overlap-threshold",
        "te_overlap_threshold",
        default=0.5,
        type=float,
        help="Model-fraction TE overlap at/above which a good-ORF gene is "
        "reclassified transposable_element (default 0.5).",
    ),
    click.option(
        "--te-class",
        "te_class",
        multiple=True,
        type=str,
        help="EDTA Classification order to treat as a TE (repeatable). Default "
        "covers true TE orders (LTR, DNA, MITE, TIR, Helitron, LINE, SINE) and "
        "excludes knob/satellite/low-complexity.",
    ),
    # tool binaries
    click.option("--transdecoder-bin-dir", "transdecoder_bin_dir", default=None),
    click.option(
        "--backstop-transdecoder/--no-backstop-transdecoder",
        "backstop_transdecoder",
        default=False,
    ),
    click.option("--portcullis-bin", "portcullis_bin", default=None),
    click.option("--mikado-bin", "mikado_bin", default="mikado"),
    click.option("--diamond-bin", "diamond_bin", default="diamond"),
    click.option(
        "--mikado-configure/--no-mikado-configure", "use_mikado_configure", default=True
    ),
    # Functional annotation (opt-in "completeness extra"). Off by
    # default → no subprocess, golden path untouched. The tool selector matches
    # FUNCTION_TOOLS in prep/function.py; emapper.py is eggNOG's binary, exposed
    # via --eggnog-bin.
    click.option(
        "--functional-annotation/--no-functional-annotation",
        "functional_annotation",
        default=False,
        help="Run the opt-in functional-annotation hook (InterProScan / eggNOG).",
    ),
    click.option(
        "--functional-tool",
        "functional_tool",
        default="interproscan",
        type=click.Choice(["interproscan", "eggnog", "both"]),
        help="Annotation tool(s) to run for --functional-annotation.",
    ),
    click.option(
        "--functional-db",
        "functional_db",
        default=None,
        help="eggNOG data dir (required for --functional-tool eggnog/both).",
    ),
    click.option(
        "--interproscan-bin",
        "interproscan_bin",
        default="interproscan.sh",
        help="InterProScan binary (bare command or absolute path).",
    ),
    click.option(
        "--eggnog-bin",
        "eggnog_bin",
        default="emapper.py",
        help="eggNOG-mapper binary (bare command or absolute path).",
    ),
    # organellar / non-standard genetic code (M-REALIZE-2). Default table 1 is
    # the standard nuclear code; --transl-table-map overrides it per-seqid so
    # organellar contigs (mito/plastid) validate against the right code. The
    # code path already exists in validate.py; these flags wire it to the CLI.
    click.option(
        "--transl-table",
        "transl_table",
        default=1,
        type=int,
        help="NCBI genetic-code table id for CDS validation "
        "(default 1 = standard nuclear code).",
    ),
    click.option(
        "--transl-table-map",
        "transl_table_map_path",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="Per-seqid genetic-code overrides; lines 'seqid=table' "
        "(e.g. chrMt=1, chrPt=11). 'seqid<TAB/space>table' also "
        "accepted. Overrides --transl-table for the named seqids.",
    ),
    # scoring + Helixer coupling levers
    click.option(
        "--scoring-profile",
        "scoring_profile",
        default="strict",
        type=click.Choice(["strict", "permissive"]),
    ),
    click.option(
        "--helixer-reference/--no-helixer-reference",
        "helixer_is_reference",
        default=True,
        help="Inject Helixer as is_reference=true.",
    ),
    click.option(
        "--helixer-support-weight",
        "helixer_support_weight",
        default=1.0,
        type=float,
        help="Scale external helixer_support (0 disables coupling).",
    ),
    # output
    click.option("--output-prefix", "output_prefix", default="helixforge"),
    click.option("--report-path", "report_path", default=None),
    click.option("--id-map-path", "id_map_path", default=None),
    click.option("--work-dir", "work_dir", default=None),
    # region / chunk
    click.option("--region", "region", default=None, help="seqid or seqid:start-end."),
    click.option("--chunk-id", "chunk_id", default=None),
    # HFG id-range allocation (Phase 16; defaults reproduce global numbering)
    click.option(
        "--id-base",
        "id_base",
        default=1,
        type=int,
        help="Lowest HFG number for this chunk's reserved range.",
    ),
    click.option(
        "--novel-base",
        "novel_base",
        default=90000,
        type=int,
        help="Lowest HFG number for this chunk's novel sub-range.",
    ),
    # resources
    click.option("--procs", "procs", default=1, type=int),
    click.option("--threads", "threads", default=4, type=int),
    # classification
    click.option("--min-tpm", "min_tpm", default=0.5, type=float),
    click.option("--min-samples", "min_samples", default=1, type=int),
    click.option("--coverage-threshold", "coverage_threshold", default=2.0, type=float),
    click.option("--near-zero-coverage", "near_zero_coverage", default=0.1, type=float),
    # AS / pick knobs
    click.option("--as-report/--no-as-report", "as_report", default=True),
    click.option(
        "--only-confirmed-introns/--allow-unconfirmed-introns",
        "only_confirmed_introns",
        default=True,
    ),
    click.option("--max-isoforms", "max_isoforms", default=5, type=int),
    click.option(
        "--keep-retained-introns/--drop-retained-introns",
        "keep_retained_introns",
        default=False,
    ),
    click.option(
        "--pad/--no-pad",
        "pad",
        default=True,
        help="Unify isoform termini (the strongest consistency lever).",
    ),
    click.option("--chimera-split/--no-chimera-split", "chimera_split", default=True),
    click.option("--flank", "flank", default=200, type=int),
    # reconciliation
    click.option("--reciprocal-overlap", "reciprocal_overlap", default=0.5, type=float),
    click.option("--min-cds-overlap", "min_cds_overlap", default=0.6, type=float),
    click.option("--min-cdna-overlap", "min_cdna_overlap", default=0.6, type=float),
    click.option("--admit-novel/--no-admit-novel", "admit_novel", default=False),
    click.option(
        "--novel-evidence-floor", "novel_evidence_floor", default=None, type=float
    ),
    # TRaCE canonical-transcript election (Phase 33b)
    click.option(
        "--trace-primary/--no-trace-primary",
        "trace_primary",
        default=False,
        help="Elect the canonical/primary isoform per gene by ranked-choice "
        "voting (TRaCE) instead of highest combined_score.",
    ),
    click.option(
        "--trace-max-aed",
        "trace_max_aed",
        default=0.5,
        type=float,
        help="TRaCE: max AED for a sample to vote for a candidate.",
    ),
    click.option(
        "--trace-min-tpm",
        "trace_min_tpm",
        default=0.5,
        type=float,
        help="TRaCE: min TPM for a sample's assembled transcript to vote.",
    ),
    click.option(
        "--trace-min-overlap",
        "trace_min_overlap",
        default=0.5,
        type=float,
        help="TRaCE: min proportion-overlap for a sample transcript to vote.",
    ),
    click.option(
        "--trace-weight-domain",
        "trace_weight_domain",
        default=9.0,
        type=float,
        help="TRaCE: domain-coverage voter weight.",
    ),
    click.option(
        "--trace-weight-protein",
        "trace_weight_protein",
        default=6.0,
        type=float,
        help="TRaCE: protein (CDS) length voter weight.",
    ),
    click.option(
        "--trace-weight-cdna",
        "trace_weight_cdna",
        default=3.0,
        type=float,
        help="TRaCE: transcript (cDNA) length voter weight.",
    ),
    click.option(
        "--trace-use-domain/--no-trace-use-domain",
        "trace_use_domain",
        default=True,
        help="TRaCE: include the domain-coverage voter when per-isoform "
        "domain coverage is available.",
    ),
    # structural validation
    click.option("--short-cds-threshold", "short_cds_threshold", default=300, type=int),
    click.option(
        "--short-exon-threshold", "short_exon_threshold", default=10, type=int
    ),
    click.option(
        "--long-intron-threshold", "long_intron_threshold", default=100_000, type=int
    ),
    # backstop junction correction
    click.option("--junction-tolerance", "junction_tolerance", default=0, type=int),
    click.option("--junction-min-reads", "junction_min_reads", default=3, type=int),
]

# click dest names for the multi-file evidence inputs (direct flag + -list
# FOFN). Handled explicitly in _build_config via expand_file_args, so they are
# skipped by the generic field copy below.
_EVIDENCE_KEYS: frozenset[str] = frozenset(
    {
        "bam",
        "bam_list",
        "star_sj",
        "star_sj_list",
        "stringtie",
        "stringtie_list",
        "bigwig",
        "bigwig_list",
    }
)


def pipeline_options(func: Callable[..., Any]) -> Callable[..., Any]:
    """Stack the shared ``PipelineConfig`` options onto a command."""
    for option in reversed(_PIPELINE_OPTIONS):
        func = option(func)
    return func


def _load_alias_map(path: str | None) -> dict[str, str] | None:
    """Load a ``{alias: canonical}`` seqid map from a JSON file or 2-column TSV."""
    if not path:
        return None
    import json
    import os

    text = open(path).read()
    if os.path.splitext(path)[1].lower() == ".json" or text.lstrip().startswith("{"):
        return {str(k): str(v) for k, v in json.loads(text).items()}
    alias: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            alias[parts[0]] = parts[1]
    return alias


def _validate_transl_table(table: int) -> int:
    """Return ``table`` if it is a known NCBI genetic-code id, else error out."""
    from helixforge.utils.sequences import CODON_TABLES

    if table not in CODON_TABLES:
        raise click.ClickException(
            f"unknown NCBI transl_table id {table}; "
            f"known tables: {sorted(CODON_TABLES)}"
        )
    return table


def _load_transl_table_map(
    path: str, *, genome_fasta: str | None = None
) -> dict[str, int]:
    """Parse a per-seqid genetic-code map file → ``{seqid: table_id}``.

    Accepts ``seqid=table`` or ``seqid<whitespace>table`` lines; ``#`` comments
    and blank lines are skipped. Every table id must be a known NCBI code, and —
    when ``genome_fasta`` is supplied and readable — every seqid must be present
    in the genome (a typo like ``ChrMt`` vs ``chrMt`` is a clear error, not a
    silent no-op). Raises :class:`click.ClickException` on any problem.
    """
    genome_seqids: set[str] | None = None
    if genome_fasta:
        try:
            from helixforge.io.fasta import GenomeAccessor

            with GenomeAccessor(genome_fasta) as genome:
                genome_seqids = set(genome.get_seqids())
        except Exception:  # noqa: BLE001 - genome unreadable → skip the seqid cross-check
            genome_seqids = None

    table_map: dict[str, int] = {}
    with open(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                seqid, _, value = line.partition("=")
            else:
                parts = line.split()
                if len(parts) < 2:
                    raise click.ClickException(
                        f"--transl-table-map line {lineno}: expected "
                        f"'seqid=table' or 'seqid table', got {line!r}"
                    )
                seqid, value = parts[0], parts[1]
            seqid = seqid.strip()
            value = value.strip()
            if not seqid:
                raise click.ClickException(
                    f"--transl-table-map line {lineno}: empty seqid in {line!r}"
                )
            try:
                table = int(value)
            except ValueError:
                raise click.ClickException(
                    f"--transl-table-map line {lineno}: table id {value!r} for "
                    f"seqid {seqid!r} is not an integer"
                ) from None
            _validate_transl_table(table)
            if genome_seqids is not None and seqid not in genome_seqids:
                preview = ", ".join(sorted(genome_seqids)[:8])
                raise click.ClickException(
                    f"--transl-table-map line {lineno}: seqid {seqid!r} is not in "
                    f"the genome FASTA (have: {preview}…)"
                )
            table_map[seqid] = table
    return table_map


def _build_config(kwargs: dict[str, Any]) -> PipelineConfig:
    """Map collected click kwargs to a :class:`PipelineConfig` (ignores extras).

    Two keys get bespoke handling before the generic field copy: ``transl_table``
    is validated against the known NCBI codes, and ``transl_table_map_path`` (a
    file path, not a config field) is parsed + validated into the
    ``transl_table_map`` dict the pipeline expects.
    """
    from helixforge.io.fofn import expand_file_args, expand_stringtie_args

    fields = {f.name for f in dataclasses.fields(PipelineConfig)}
    cfg_kwargs: dict[str, Any] = {}

    # Multi-file evidence inputs: merge the repeatable/comma flag with its
    # ``-list`` FOFN companion (v1 convention) before the generic copy. Each
    # downstream PipelineConfig field already takes a plain list of paths.
    if "bam" in kwargs or "bam_list" in kwargs:
        cfg_kwargs["bam_paths"] = expand_file_args(
            kwargs.get("bam", ()), kwargs.get("bam_list"), label="BAM file"
        )
    if "star_sj" in kwargs or "star_sj_list" in kwargs:
        cfg_kwargs["star_sj_paths"] = expand_file_args(
            kwargs.get("star_sj", ()),
            kwargs.get("star_sj_list"),
            label="SJ.out.tab file",
        )
    if "stringtie" in kwargs or "stringtie_list" in kwargs:
        cfg_kwargs["stringtie_list"] = expand_stringtie_args(
            kwargs.get("stringtie", ()),
            kwargs.get("stringtie_list"),
            label="StringTie GTF",
        )
    if "bigwig" in kwargs or "bigwig_list" in kwargs:
        cfg_kwargs["bigwig_paths"] = expand_file_args(
            kwargs.get("bigwig", ()),
            kwargs.get("bigwig_list"),
            label="bigWig file",
        )
    # --te-class is repeatable; an empty selection means "use the default TE
    # class set" (None), not "treat nothing as a TE" ([]).
    if "te_class" in kwargs:
        te_classes = list(kwargs.get("te_class", ()))
        cfg_kwargs["te_classes"] = te_classes or None

    for key, value in kwargs.items():
        if key in _EVIDENCE_KEYS:
            continue
        field = key
        if field not in fields:
            continue
        if isinstance(value, tuple):  # click multiple=True → tuple → list
            value = list(value)
        cfg_kwargs[field] = value

    if "transl_table" in cfg_kwargs and cfg_kwargs["transl_table"] is not None:
        _validate_transl_table(int(cfg_kwargs["transl_table"]))
    map_path = kwargs.get("transl_table_map_path")
    if map_path:
        cfg_kwargs["transl_table_map"] = _load_transl_table_map(
            map_path, genome_fasta=kwargs.get("genome_fasta")
        )
    return PipelineConfig(**cfg_kwargs)


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


class IntentGroup(click.Group):
    """Top-level group that lists commands by user intent, not alphabetically.

    The grouping is purely a help-presentation concern (titled sections +
    essentials-first ordering); it adds no nesting. Hidden commands (e.g. the
    dev-only ``benchmark``) and any command not assigned to a section are
    omitted from the sections and collected under "Other commands".
    """

    #: section title → ordered command names (essentials first).
    SECTIONS: list[tuple[str, list[str]]] = [
        ("Annotate", ["reconcile", "parallel"]),
        ("Score & inspect", ["confidence", "evidence", "stats"]),
        ("Preflight & utilities", ["doctor", "viz", "utils"]),
    ]

    def format_commands(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        listed: set[str] = set()
        for title, names in self.SECTIONS:
            rows: list[tuple[str, str]] = []
            for name in names:
                cmd = self.get_command(ctx, name)
                if cmd is None or cmd.hidden:
                    continue
                rows.append((name, cmd.get_short_help_str(limit=78)))
                listed.add(name)
            if rows:
                with formatter.section(title):
                    formatter.write_dl(rows)
        extra: list[tuple[str, str]] = []
        for name in self.list_commands(ctx):
            if name in listed:
                continue
            cmd = self.get_command(ctx, name)
            if cmd is None or cmd.hidden:
                continue
            extra.append((name, cmd.get_short_help_str(limit=78)))
        if extra:
            with formatter.section("Other commands"):
                formatter.write_dl(extra)


@click.group(cls=IntentGroup)
@click.version_option(package_name="helixforge", message="%(prog)s %(version)s")
def main() -> None:
    """HelixForge v3 — isoform-aware refinement of Helixer annotations.

    \b
    Annotate:            reconcile (one region, whole genome, or --scatter
                         chunked) and parallel (distributed scatter-gather).
    Score & inspect:     confidence, evidence, stats (read-only; never edit models).
    Preflight & utils:   doctor (validate inputs + tools), viz, utils.
    """


@main.command(
    epilog="""\b
Examples:
  # RNA-seq (StringTie + BAM) + protein homology, strict scoring (single pass)
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \\
      --helixer-h5 helixer.h5 --stringtie sampleA.gtf,sampleB.gtf \\
      --bam sampleA.bam,sampleB.bam --protein-db proteins.fa \\
      --output-prefix run1
\b
  # Protein-only backstop CDS from a precomputed miniprot GFF, permissive scoring
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \\
      --miniprot miniprot.gff --scoring-profile permissive \\
      --output-prefix run1
\b
  # Chunked whole genome: partition into gene-safe chunks, run locally, aggregate
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \\
      --helixer-h5 helixer.h5 --stringtie-list stringtie.list \\
      --scatter auto --workers 16 --output-prefix maize
\b
  # Chunked on a cluster: emit a Slurm array (submit it, then `parallel aggregate`)
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \\
      --scatter 64 --hpc slurm --output-prefix maize
\b
  # One region only with a reserved id range (the shape distributed runs emit)
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \\
      --helixer-h5 helixer.h5 --region chr1:1-2000000 \\
      --id-base 1 --output-prefix chunk_chr1
"""
)
@pipeline_options
@click.option(
    "--scatter",
    "scatter",
    default="off",
    show_default=True,
    help="Chunked execution: off | auto | <N chunks>. 'off' is a single "
    "in-process pass; 'auto'/<N> partition the genome into gene-safe chunks "
    "(same result, less memory). For a scheduler-driven distributed run drive "
    "`parallel plan/tasks/aggregate` yourself.",
)
@click.option(
    "--hpc",
    "hpc",
    type=click.Choice(["local", "hypershell", "slurm"]),
    default="local",
    show_default=True,
    help="Chunked mode only: run chunks locally now, run them now via HyperShell "
    "(`hs cluster`), or emit a Slurm array to submit.",
)
@click.option(
    "--workers",
    "workers",
    default=4,
    type=int,
    show_default=True,
    help="Chunked --hpc local: process-pool size (also the default HyperShell "
    "--num-tasks concurrency).",
)
@click.option(
    "--hs-bin",
    "hs_bin",
    default="hs",
    show_default=True,
    help="Chunked --hpc hypershell: the HyperShell executable (install via "
    "`pip install helixforge[hpc]`).",
)
@click.option(
    "--target-loci-per-chunk", "target_loci_per_chunk", default=None, type=int
)
@click.option("--min-boundary-gap", "min_boundary_gap", default=None, type=int)
@click.option(
    "--out-prefix",
    "out_prefix",
    default=None,
    help="Chunked mode: prefix for the merged outputs + manifest "
    "(defaults to --output-prefix).",
)
@click.option(
    "--master-id-map",
    "master_id_map_path",
    default=None,
    help="Chunked mode: genome-wide master id_map.json (defaults to --id-map-path).",
)
@click.option(
    "--stitch/--no-stitch",
    "stitch",
    default=False,
    help="Chunked mode: report boundary-recoverable cross-chunk merges after aggregate.",
)
@click.option(
    "--finalize-workers",
    "finalize_workers",
    default=1,
    type=int,
    help="Worker processes for the per-gene finalize stage (intra-chunk).",
)
@click.option(
    "--seqid-aliases",
    "seqid_aliases",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Chunked mode: seqid alias map (JSON {alias: canonical} or 2-col TSV) "
    "for the integrity preflight.",
)
@click.option(
    "--skip-preflight",
    "skip_preflight",
    is_flag=True,
    default=False,
    help="Chunked mode: skip the input-integrity preflight (not recommended).",
)
@click.option(
    "--resume/--no-resume",
    "resume",
    default=True,
    show_default=True,
    help="Reuse completed stage outputs (PREP, Mikado) from the work dir "
    "instead of re-running them. A stage is skipped only if its outputs "
    "exist and validate (non-empty + content hash). --no-resume forces a "
    "fresh run even when valid outputs exist.",
)
@click.option(
    "--mikado-loci",
    "mikado_loci_gff3",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to an externally-run Mikado loci GFF3 (e.g. mikado.loci.gff3). "
    "Skips the entire PREP + Mikado stage; companion files "
    "(*.metrics.tsv, *.scores.tsv) are auto-detected next to the GFF3 "
    "and used when present. The rest of the pipeline (backstop, "
    "reconcile, output) runs normally on the supplied loci.",
)
def reconcile(
    scatter: str,
    hpc: str,
    workers: int,
    hs_bin: str,
    target_loci_per_chunk: int | None,
    min_boundary_gap: int | None,
    out_prefix: str | None,
    master_id_map_path: str | None,
    stitch: bool,
    seqid_aliases: str | None,
    skip_preflight: bool,
    **kwargs: Any,
) -> None:
    """Reconcile Helixer models with evidence → tiered GFF3 + per-gene report.

    The main HelixForge command. It runs PREP → MIKADO → RECONCILE → OUTPUT end
    to end: prepares Mikado inputs from the Helixer models plus RNA-seq and
    protein evidence, reconciles Mikado loci against the Helixer gene set
    (rescuing silent genes, stabilising HFG ids, calling isoforms), applies the
    structural codon gate, and writes the tiered annotation. Unlike ``evidence``
    (which only *scores* an existing GFF3) this command *builds and fixes* models.

    \b
    Execution mode:
    - default (``--scatter off``): one in-process pass over the whole genome
      or a single ``--region``.
    - ``--scatter auto|<N>``: partition the genome into gene-safe chunks, run
      them (``--hpc local`` now, or ``--hpc slurm`` to emit an array to submit),
      and aggregate into one verified annotation (+ a run manifest). Chunking
      changes memory/speed, not the result. For a fully scheduler-driven
      distributed run, drive ``parallel plan/tasks/aggregate`` directly.

    \b
    Outputs (``--output-prefix`` is the stem):
    - <prefix>.gff3                       full reconciled annotation (1-based GFF3)
    - <prefix>.tier1/2/3.gff3             cumulative per-tier GFF3 subsets
    - <prefix>.report.tsv                 per-gene report (one row per gene)
    - <prefix>.id_map.json                Helixer-locus → HFG id map (stable reruns)
    - <prefix>.report.json / .report.html / .run_stats.json   run summaries
    - (chunked) <out-prefix>.manifest.json   plan + chunk prefixes + id_map

    \b
    Per-gene report (<prefix>.report.tsv) columns:
    gene_id, seqid, start, end, strand, tier, origin, biotype,
    primary_transcript_id, num_isoforms, num_as_events, has_cds, protein_id,
    max_tpm, junction_support, helixer_support, combined_score, aed, flags.
    """
    console = _console()
    # finalize_workers stays in kwargs (a PipelineConfig field) → _build_config.
    config = _build_config(kwargs)

    # --- single in-process pass: identical to the historical `reconcile` ---
    if scatter == "off":
        console.print(
            f"[bold]Reconciling[/bold] → prefix [cyan]{config.output_prefix}[/cyan] "
            f"(profile: {config.scoring_profile})"
        )
        genes = run_pipeline(config)
        console.print(
            f"[green]Done[/green]: {len(genes)} genes → {config.output_prefix}.gff3 "
            f"(+ tier1/2/3), report {config.report_path}"
        )
        return

    # --- chunked scatter-gather path (plan → run chunks → aggregate) ---
    from helixforge.parallel.run import run_genome
    from helixforge.prep.preflight import PreflightError, require_preflight

    # Required input-integrity gate before a whole-genome scatter run: refuse to
    # run on inconsistent/malformed inputs. Non-mutating, so a passing gate is
    # count-neutral.
    if not skip_preflight:
        try:
            report = require_preflight(
                config, alias_map=_load_alias_map(seqid_aliases)
            )
        except PreflightError as exc:
            raise click.ClickException(str(exc)) from exc
        for warning in report.warnings():
            console.print(f"[yellow]preflight warning[/yellow]: {warning}")

    # Resolve --scatter: "auto" | an integer chunk count.
    scatter_val: str | int = scatter
    if scatter != "auto":
        try:
            scatter_val = int(scatter)
        except ValueError:
            raise click.ClickException(
                f"--scatter must be off | auto | <int>, got {scatter!r}"
            )
    console.print(
        f"[bold]Reconciling (chunked)[/bold] scatter=[cyan]{scatter}[/cyan] "
        f"hpc={hpc} → prefix [cyan]{out_prefix or config.output_prefix}[/cyan]"
    )
    result = run_genome(
        config,
        scatter=scatter_val,
        hpc=hpc,
        workers=workers,
        hs_bin=hs_bin,
        target_loci_per_chunk=target_loci_per_chunk,
        min_boundary_gap=min_boundary_gap,
        out_prefix=out_prefix,
        master_id_map_path=master_id_map_path,
        stitch=stitch,
    )
    if result.mode == "single":
        assert result.genes is not None
        console.print(
            f"[green]Done[/green]: {len(result.genes)} genes → "
            f"manifest {result.manifest_path}"
        )
    elif result.mode == "scatter-slurm":
        assert result.plan is not None
        console.print(
            f"[green]Wrote[/green] Slurm array {result.script_path} "
            f"({len(result.plan)} chunks) → submit, then aggregate. "
            f"Manifest {result.manifest_path}"
        )
    else:
        assert result.aggregate is not None
        agg = result.aggregate
        console.print(
            f"[green]Done[/green]: {agg.num_genes} genes, {agg.num_loci} "
            f"loci over {len(result.chunk_prefixes)} chunks; "
            f"recovered merges={result.boundary_merges_recovered}. "
            f"Manifest {result.manifest_path}"
        )


@main.command(
    epilog="""\b
Examples:
  # just the external-tool preflight (mikado, diamond, …)
  helixforge doctor
\b
  # tools + input-integrity (seqid concordance / format sniff / BAM index)
  helixforge doctor --genome genome.fa --helixer helixer.gff3 \\
      --helixer-h5 helixer.h5 --bam sampleA.bam --star-sj sampleA.SJ.out.tab
\b
  # verify an emitted Mikado config/scoring against the detected mikado version
  helixforge doctor --check-config configuration.yaml --check-scoring scoring.yaml
"""
)
@click.option(
    "--mikado-bin", "mikado_bin", default=None, help="Override the mikado binary."
)
@click.option(
    "--diamond-bin", "diamond_bin", default=None, help="Override the diamond binary."
)
@click.option(
    "--check-config",
    "check_config",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Verify an emitted Mikado configuration.yaml against the detected version.",
)
@click.option(
    "--check-scoring",
    "check_scoring",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Scoring YAML to verify alongside --check-config.",
)
# Input-integrity preflight inputs: when --genome + --helixer are
# given, doctor also runs the seqid-concordance / format-sniffer / CSI gate.
@click.option(
    "--genome",
    "genome_fasta",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Genome FASTA (integrity check).",
)
@click.option(
    "--helixer",
    "helixer_gff3",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Helixer GFF3 (integrity check).",
)
@click.option(
    "--helixer-h5",
    "helixer_h5",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Helixer HDF5 (integrity check).",
)
@click.option(
    "--stringtie",
    "stringtie",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False),
    help="StringTie GTF (repeatable).",
)
@click.option(
    "--bam",
    "bam",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False),
    help="RNA-seq BAM (repeatable).",
)
@click.option(
    "--star-sj",
    "star_sj",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False),
    help="STAR SJ.out.tab (repeatable).",
)
@click.option(
    "--miniprot",
    "miniprot_gff",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="miniprot GFF (integrity check).",
)
@click.option(
    "--reference",
    "reference",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Genome FASTA for CRAM decode — doctor reports whether CRAM "
    "BAM inputs have an offline reference (--reference or REF_CACHE).",
)
@click.option(
    "--seqid-aliases",
    "seqid_aliases",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Seqid alias map (JSON {alias: canonical} or 2-col TSV).",
)
def doctor(
    mikado_bin: str | None,
    diamond_bin: str | None,
    check_config: str | None,
    check_scoring: str | None,
    genome_fasta: str | None,
    helixer_gff3: str | None,
    helixer_h5: str | None,
    stringtie: tuple[str, ...],
    bam: tuple[str, ...],
    star_sj: tuple[str, ...],
    miniprot_gff: str | None,
    reference: str | None,
    seqid_aliases: str | None,
) -> None:
    """Preflight: validate inputs + resolve external tools before a run.

    The single input-validation entry point — run it before a real
    ``reconcile``. With no inputs it resolves and versions every external tool
    the pipeline shells out to (mikado, diamond, transdecoder, portcullis, …)
    and checks environment hygiene (console-script shim, CRAM reference). Given
    ``--genome`` + ``--helixer`` (and optional evidence) it also runs the
    input-integrity preflight: seqid concordance across inputs, GFF3/GTF/FASTA/BAM
    format sniffing, and BAM/CRAM index (CSI) presence. ``--check-config``
    verifies an emitted Mikado config/scoring file against the detected Mikado
    version.

    Nothing is modified and nothing is annotated — this command only inspects.
    It exits non-zero if a required tool is missing, an emitted config fails
    schema verification, or the input-integrity gate fails.

    \b
    Reports (printed to the console; advisory checks warn, gates fail):
    - external-tool table: tool, resolved path, detected version, status
    - environment hygiene: console-script shim, CRAM-reference readiness
    - (optional) Mikado config/scoring verification verdict
    - (optional) input-integrity report: seqid concordance, format, index
    """
    from types import SimpleNamespace

    from helixforge.prep.doctor import (
        check_cram_reference,
        check_scoring_profiles,
        check_shim,
        check_tools,
        verify_emitted_config,
    )
    from helixforge.prep.preflight import run_preflight

    console = _console()
    overrides: dict[str, str] = {}
    if mikado_bin:
        overrides["mikado"] = mikado_bin
        overrides["mikado_compare"] = mikado_bin
    if diamond_bin:
        overrides["diamond"] = diamond_bin
    report = check_tools(overrides)
    console.print(report.render())

    scoring_status = check_scoring_profiles()
    console.print(scoring_status.render())

    # Environment hygiene: a mis-pointed console-script shim and
    # CRAM-reference readiness. Both are advisory (warn, never block).
    shim = check_shim()
    if shim.foreign:
        console.print(f"[yellow]warning[/yellow]: {shim.render()}")
    else:
        console.print(shim.render())
    cram = check_cram_reference(list(bam), reference=reference)
    if cram.has_cram and not cram.ok:
        console.print(f"[yellow]warning[/yellow]: {cram.render()}")
    elif cram.has_cram:
        console.print(cram.render())

    config_bad = False
    if check_config:
        verification = verify_emitted_config(
            check_config,
            scoring_path=check_scoring,
            mikado_bin=(mikado_bin or "mikado"),
        )
        console.print(verification.render())
        config_bad = not verification.ok

    preflight_bad = False
    if genome_fasta and helixer_gff3:
        cfg = SimpleNamespace(
            genome_fasta=genome_fasta,
            helixer_gff3=helixer_gff3,
            helixer_h5=helixer_h5,
            stringtie_list=list(stringtie),
            bam_paths=list(bam),
            star_sj_paths=list(star_sj),
            miniprot_gff=miniprot_gff,
        )
        preflight = run_preflight(cfg, alias_map=_load_alias_map(seqid_aliases))
        console.print(preflight.render())
        preflight_bad = not preflight.ok()

    if not report.ok():
        missing = ", ".join(s.key for s in report.missing_required())
        raise click.ClickException(f"required tool(s) missing: {missing}")
    if config_bad:
        raise click.ClickException("emitted Mikado config failed schema verification")
    if preflight_bad:
        raise click.ClickException("input-integrity preflight failed")


def _require_viz(plot_format: str, console: Any) -> bool:
    """Return True if the plotting backend for ``plot_format`` is importable.

    The standalone ``confidence`` plot options need the optional ``viz`` extra
    (``html`` → plotly, ``png``/``pdf`` → matplotlib). When the backend is
    absent we print a clear install hint and return False so the command still
    completes (the TSV/BED outputs were already written) — visualization
    degrades gracefully instead of crashing with an ImportError traceback.
    """
    import importlib

    backend = "plotly" if plot_format == "html" else "matplotlib"
    try:
        importlib.import_module(backend)
    except ImportError:
        console.print(
            f"[yellow]Skipping plots:[/yellow] '{backend}' is not installed "
            f"(needed for --format {plot_format}). "
            "Install it with: [bold]pip install 'helixforge[viz]'[/bold]"
        )
        return False
    return True


@main.command(
    epilog="""\b
Examples:
  # explicit input HDF5 (preferred — enables strand-aware coordinate mapping)
  helixforge confidence -p predictions.h5 -g genes.gff3 \\
      --input-h5 input.h5 -o scores.tsv
\b
  # auto-detect *_input.h5 next to the predictions file
  helixforge confidence -p predictions.h5 -g genes.gff3 -o scores.tsv
\b
  # FASTA-index fallback + a genome-wide distribution plot and summary TSV
  helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa \\
      -o scores.tsv --distribution-plot dist.html --summary-tsv summary.tsv
\b
  # one scaffold only, low-confidence plots, for a parallel chunk
  helixforge confidence -p predictions.h5 -g genes.gff3 --input-h5 input.h5 \\
      -o chr1.scores.tsv --scaffold chr1 --plot-dir plots --plot-threshold 0.7 \\
      --chunk-id chr1
"""
)
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
    help="Directory for per-gene confidence plots. Use with --plot-threshold or --max-plots to limit output for large genomes.",
)
@click.option(
    "--plot-threshold",
    type=float,
    default=None,
    help="Only plot genes with confidence below this threshold (e.g., 0.7). Recommended for large genomes.",
)
@click.option(
    "--max-plots",
    type=int,
    default=None,
    help="Maximum number of gene plots to generate. Plots lowest-confidence genes first.",
)
@click.option(
    "--distribution-plot",
    type=click.Path(path_type=Path),
    help="Output path for genome-wide distribution plot.",
)
@click.option(
    "--summary-tsv",
    type=click.Path(path_type=Path),
    help="Write the on-screen distribution summary (metric/value, long format) to this TSV.",
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
def confidence(
    predictions: Path,
    gff: Path,
    genome: Path | None,
    input_h5: Path | None,
    output: Path,
    bed: Path | None,
    low_conf_bed: Path | None,
    plot_dir: Path | None,
    plot_threshold: float | None,
    max_plots: int | None,
    distribution_plot: Path | None,
    summary_tsv: Path | None,
    threshold: float,
    threads: int,
    plot_format: str,
    region: str | None,
    chunk_id: str | None,
    scaffold: str | None,
) -> None:
    """Inspect (read-only): score genes against the Helixer HDF5 confidence track.

    A standalone, HDF5-only scorer: it reads the Helixer softmax predictions and
    computes multi-factor confidence metrics per gene (class probabilities,
    Shannon entropy, boundary sharpness, CDS coding consistency, per-exon
    scores). It needs *only* the Helixer HDF5 + a GFF3 — no Mikado, no RNA-seq,
    no external toolchain. Use ``evidence`` instead to score against RNA-seq /
    protein evidence, and ``reconcile`` to actually build/fix models; this
    command never modifies a model, it only annotates confidence.

    Genes are classed high (>=0.85), medium (>=0.70), or low (<0.70), but the
    scores are genome-relative — prefer a cutoff from the printed distribution
    (or ``--summary-tsv``) over the fixed class thresholds.

    Coordinate mapping comes from (in order of preference): ``--input-h5``
    (strand-aware), auto-detected ``*_input.h5`` next to the predictions, or
    ``--genome`` FASTA + ``.fai``. For parallel runs use ``--region`` /
    ``--scaffold`` to subset and ``--chunk-id`` for output organization.

    \b
    Outputs:
    - -o/--output TSV, one row per gene, columns:
      gene_id, seqid, start, end, strand, mean_prob, min_prob, median_prob,
      entropy, boundary_sharpness, coding_consistency, worst_exon_score,
      overall_score, confidence_class, flags, n_low_conf_regions, n_exons
    - --bed: BED12 of genes coloured by confidence; --low-conf-bed: BED of
      low-confidence sub-regions
    - --summary-tsv: long-format metric/value distribution (n, mean, median,
      std, min, max, p5/p25/p50/p75/p95; per-component and per-exon stats;
      n_high/n_medium/n_low)
    - --distribution-plot / --plot-dir: genome-wide and per-gene plots
    """
    from helixforge.core.confidence import (
        ConfidenceCalculator,
        ConfidenceWriter,
    )
    from helixforge.core.gff import GFF3Parser
    from helixforge.core.hdf5 import HelixerHDF5Reader
    from helixforge.core.regions import (
        GenomicRegion,
        parse_region,
        region_from_scaffold,
    )
    from helixforge.io.fasta import GenomeAccessor

    console = _console()
    # v3's ``main`` group has no global --verbose/--quiet; the standalone
    # confidence command prints its progress lines unconditionally.
    verbose = False
    quiet = False

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
            console.print(
                "[blue]Coordinate mapping:[/blue] auto-detect from predictions path"
            )

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
                        console.print(
                            f"[blue]Processing region:[/blue] {target_region}"
                        )
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
                console.print(f"[blue]Scoring genes with {threads} thread(s)...[/blue]")

            scores = list(calc.score_genes_parallel(genes, n_workers=threads))

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
            if distribution_plot and _require_viz(plot_format, console):
                from helixforge.core.plots import plot_confidence_distribution

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
            if plot_dir and _require_viz(plot_format, console):
                from helixforge.core.plots import plot_gene_confidence_batch

                # Check if a file exists with the directory name
                if plot_dir.exists() and not plot_dir.is_dir():
                    console.print(
                        f"[red]Error:[/red] '{plot_dir}' exists but is not a directory. "
                        "Please remove it or use a different --plot-dir path."
                    )
                    raise SystemExit(1)

                # Filter genes/scores for plotting
                plot_genes = genes
                plot_scores = scores

                # Apply threshold filter if specified
                if plot_threshold is not None:
                    filtered = [
                        (g, s)
                        for g, s in zip(genes, scores)
                        if s.overall_score < plot_threshold
                    ]
                    plot_genes = [g for g, _ in filtered]
                    plot_scores = [s for _, s in filtered]
                    if not quiet:
                        console.print(
                            f"[dim]Filtering plots to {len(plot_genes)} genes "
                            f"with confidence < {plot_threshold}[/dim]"
                        )

                # Sort by confidence (lowest first) and apply max_plots limit
                if max_plots is not None and len(plot_genes) > max_plots:
                    # Sort by confidence score ascending
                    sorted_pairs = sorted(
                        zip(plot_genes, plot_scores),
                        key=lambda x: x[1].overall_score,
                    )
                    plot_genes = [g for g, _ in sorted_pairs[:max_plots]]
                    plot_scores = [s for _, s in sorted_pairs[:max_plots]]
                    if not quiet:
                        console.print(
                            f"[dim]Limiting to {max_plots} lowest-confidence genes[/dim]"
                        )

                # Warn if generating many plots
                if (
                    len(plot_genes) > 500
                    and plot_threshold is None
                    and max_plots is None
                ):
                    console.print(
                        f"[yellow]Warning:[/yellow] About to generate {len(plot_genes)} "
                        "individual plot files. This may take a while and use significant "
                        "disk space.\n"
                        "  Consider using --plot-threshold to only plot low-confidence genes,\n"
                        "  or --max-plots to limit the number of plots."
                    )

                if len(plot_genes) == 0:
                    if not quiet:
                        console.print(
                            "[dim]No genes match the plotting criteria. "
                            "No plots generated.[/dim]"
                        )
                else:
                    plot_dir.mkdir(parents=True, exist_ok=True)
                    paths = plot_gene_confidence_batch(
                        plot_genes,
                        plot_scores,
                        plot_dir,
                        format=plot_format,
                        calc=calc,
                    )
                    if not quiet:
                        console.print(
                            f"[green]Generated {len(paths)} gene plots in:[/green] {plot_dir}"
                        )

            # Distribution summary (stats over the overall score + components).
            # Always computed so --summary-tsv works under --quiet; printed only
            # when not quiet. Additive: the high/medium/low counts below stay.
            from helixforge.core.confidence import (
                summarize_confidence_scores,
                write_confidence_summary_tsv,
            )

            summary = summarize_confidence_scores(scores)

            if summary_tsv:
                write_confidence_summary_tsv(summary, summary_tsv)
                if not quiet:
                    console.print(
                        f"[green]Wrote summary stats to:[/green] {summary_tsv}"
                    )

            if not quiet:
                from rich.table import Table

                table = Table(title="Confidence distribution (overall score)")
                table.add_column("metric")
                table.add_column("value", justify="right")
                for metric, value in summary:
                    if value is None:
                        rendered = "—"
                    elif isinstance(value, float):
                        rendered = f"{value:.4f}"
                    else:
                        rendered = str(value)
                    table.add_row(metric, rendered)
                console.print(table)
                console.print(
                    "[dim]Scores are genome-relative: pick a cutoff from this "
                    "distribution rather than the fixed 0.85/0.70 class thresholds.[/dim]"
                )

            # Print class-count summary (existing output; unchanged).
            if not quiet:
                high_count = sum(1 for s in scores if s.confidence_class == "high")
                medium_count = sum(1 for s in scores if s.confidence_class == "medium")
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


@main.command(
    epilog="""\b
Examples:
  # RNA-seq only (junctions + coverage) → rna_aed
  helixforge evidence --gff3 annotation.gff3 \\
      --bam sampleA.bam,sampleB.bam --sj sampleA.SJ.out.tab \\
      --out evidence.tsv
\b
  # add protein evidence by aligning a TE-filtered proteome with miniprot
  helixforge evidence --gff3 annotation.gff3 --bam sampleA.bam \\
      --proteins refprot.fa --genome genome.fa --out evidence.tsv
\b
  # reuse a precomputed miniprot GFF (no realignment) + StringTie TPM
  helixforge evidence --gff3 annotation.gff3 --stringtie-list stringtie.list \\
      --proteins refprot.fa --miniprot-gff miniprot.gff \\
      --out evidence.tsv --gene-out evidence.gene.tsv
"""
)
@click.option(
    "--gff3",
    "gff3_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Annotation GFF3 to score (any GFF3, not just HelixForge output).",
)
@click.option(
    "--bam",
    "bam",
    multiple=True,
    type=str,
    help="RNA-seq BAM (sorted + indexed). Comma-separated or repeated.",
)
@click.option(
    "--bam-list",
    "bam_list",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="File of BAM paths, one per line (blank lines and # comments ignored).",
)
@click.option(
    "--sj",
    "star_sj",
    multiple=True,
    type=str,
    help="STAR SJ.out.tab. Comma-separated or repeated.",
)
@click.option(
    "--sj-list",
    "sj_list",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="File of SJ.out.tab paths, one per line (blank lines and # comments ignored).",
)
@click.option(
    "--stringtie",
    "stringtie",
    multiple=True,
    type=str,
    help="StringTie GTF, one per sample (TPM). Comma-separated or repeated.",
)
@click.option(
    "--stringtie-list",
    "stringtie_list",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="File of StringTie GTF paths, one per line "
    "(blank lines and # comments ignored).",
)
@click.option(
    "--proteins",
    "proteins",
    multiple=True,
    type=str,
    help="Reference proteome FASTA for the protein-AED axis. "
    "Comma-separated or repeated. Aligned with miniprot (needs "
    "--genome) unless --miniprot-gff is supplied. NOTE: use a "
    "TE-filtered proteome — TE proteins align well and would give "
    "TE models a deceptively low protein_aed.",
)
@click.option(
    "--proteins-list",
    "proteins_list",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="File of proteome FASTA paths, one per line "
    "(blank lines and # comments ignored).",
)
@click.option(
    "--miniprot-gff",
    "miniprot_gff",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Precomputed miniprot GFF3 for --proteins; skips realigning on "
    "re-runs (parsed directly, no miniprot needed).",
)
@click.option(
    "--genome",
    "genome",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Genome FASTA — the miniprot target when aligning --proteins "
    "from scratch. Not needed with --miniprot-gff.",
)
@click.option(
    "--reference",
    "reference",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Genome FASTA for CRAM decode (passed as reference_filename to "
    "pysam, so CRAM never triggers a remote ENA fetch). Ignored for "
    "BAM. If omitted, htslib honours a REF_CACHE / REF_PATH env cache.",
)
@click.option(
    "--region",
    "region",
    default=None,
    help="Restrict to seqid or seqid:start-end (1-based).",
)
@click.option(
    "--out", "out_path", default="evidence.tsv", help="Per-transcript TSV output path."
)
@click.option(
    "--gene-out",
    "gene_out_path",
    default=None,
    help="Per-gene TSV (best-supported transcript per gene). "
    "Defaults to '<out>.gene.tsv'.",
)
@click.option(
    "--min-reads",
    "min_reads",
    default=3,
    type=int,
    help="Minimum junction read support to qualify as evidence.",
)
@click.option(
    "--min-mapq",
    "min_mapq",
    default=10,
    type=int,
    help="Minimum read MAPQ for BAM junction extraction.",
)
@click.option(
    "--min-overhang",
    "min_overhang",
    default=8,
    type=int,
    help="Minimum spliced-read overhang (bp) on each side of a junction.",
)
@click.option(
    "-j",
    "--threads",
    "threads",
    default=1,
    type=int,
    show_default=True,
    help="Parallelism dial: process-parallel evidence extraction "
    "(one BAM/SJ per worker — a single coverage pass per locus, not "
    "per exon — plus the miniprot run) then process-parallel per-gene "
    "scoring. Output is identical to -j 1.",
)
def evidence(
    gff3_path: str,
    bam: tuple[str, ...],
    bam_list: str | None,
    star_sj: tuple[str, ...],
    sj_list: str | None,
    stringtie: tuple[str, ...],
    stringtie_list: str | None,
    proteins: tuple[str, ...],
    proteins_list: str | None,
    miniprot_gff: str | None,
    genome: str | None,
    reference: str | None,
    region: str | None,
    out_path: str,
    gene_out_path: str | None,
    min_reads: int,
    min_mapq: int,
    min_overhang: int,
    threads: int,
) -> None:
    """Inspect (read-only): score any GFF3 against RNA-seq + protein evidence.

    A standalone, evidence-only scorer. It works on *any* GFF3 (not just
    HelixForge output) and needs only RNA-seq (BAM / STAR SJ / StringTie) and/or
    protein evidence — no Mikado, no Helixer HDF5, no external toolchain beyond
    an optional miniprot run. This command only *scores*; it never modifies a
    model (use ``reconcile`` for that), and it scores against evidence rather
    than the Helixer track (use ``confidence`` for that).

    Reports two **separate** AED scores (never fused): ``rna_aed`` (junctions +
    coverage + boundary, from BAM/SJ) and ``protein_aed`` (intron structure + CDS
    coverage + reference-protein coverage, from miniprot). Each is in [0, 1],
    lower = better, and is populated only when its evidence was supplied —
    missing evidence is neutral, not a penalty.

    \b
    Outputs:
    - --out: per-transcript TSV, one row per transcript, columns:
      gene_id, transcript_id, seqid, strand, start, end, num_exons, num_introns,
      supported, contradicted, novel_in_data, junction_support_fraction,
      intron_precision, intron_recall, intron_f1, mean_coverage, tpm,
      then the RNA-AED block (when BAM coverage given): rna_aed,
      rna_junction_ratio, rna_coverage_ratio, rna_boundary_ratio,
      and the protein-AED block (when proteins/miniprot given): protein_id,
      protein_aed, protein_struct_ratio, protein_cds_cov_ratio,
      protein_prot_cov_ratio
    - --gene-out: per-gene TSV (best-supported transcript per gene), same columns;
      defaults to '<out>.gene.tsv'
    - a printed summary table (n_transcripts, fraction_fully_supported, mean
      intron F1, mean rna_aed / protein_aed, …)
    """
    from helixforge.io.fofn import expand_file_args, expand_stringtie_args
    from helixforge.score.evidence import (
        rollup_genes,
        score_annotation,
        summarize_evidence,
        write_evidence_tsv,
    )

    console = _console()
    bam_paths = expand_file_args(bam, bam_list, label="BAM file")
    sj_paths = expand_file_args(star_sj, sj_list, label="SJ.out.tab file")
    stringtie_gtfs = expand_stringtie_args(
        stringtie,
        stringtie_list,
        label="StringTie GTF",
        warn=lambda msg: console.print(f"[yellow]deprecated[/yellow]: {msg}"),
    )
    protein_paths = expand_file_args(proteins, proteins_list, label="proteome FASTA")

    if (
        not bam_paths
        and not sj_paths
        and not stringtie_gtfs
        and not protein_paths
        and not miniprot_gff
    ):
        raise click.ClickException(
            "at least one evidence source required: --bam, --sj, --stringtie, "
            "--proteins, or --miniprot-gff"
        )
    if protein_paths and not miniprot_gff and not genome:
        raise click.ClickException(
            "--proteins needs either --miniprot-gff (precomputed) or --genome "
            "(to align with miniprot)"
        )

    df = score_annotation(
        gff3_path,
        bam_paths=bam_paths or None,
        star_sj_paths=sj_paths or None,
        stringtie_gtfs=stringtie_gtfs or None,
        proteins=protein_paths or None,
        miniprot_gff=miniprot_gff,
        genome=genome,
        region=region,
        min_reads=min_reads,
        min_mapq=min_mapq,
        min_overhang=min_overhang,
        reference_filename=reference,
        threads=threads,
    )
    write_evidence_tsv(df, out_path)

    gene_out = gene_out_path or f"{out_path}.gene.tsv"
    gene_df = rollup_genes(df)
    write_evidence_tsv(gene_df, gene_out)

    summary = summarize_evidence(df)

    from rich.table import Table

    table = Table(title="Evidence summary (RNA + protein AED)")
    table.add_column("metric")
    table.add_column("value", justify="right")
    for key, value in summary.items():
        rendered = (
            "—"
            if value is None
            else (f"{value:.4f}" if isinstance(value, float) else str(value))
        )
        table.add_row(key, rendered)
    console.print(table)
    console.print(
        f"[green]Wrote[/green] {len(df)} transcript rows → {out_path}  "
        f"({len(gene_df)} genes → {gene_out})"
    )


@main.command(
    epilog="""\b
Examples:
  # before/after on two GFF3s
  helixforge stats --helixer helixer.gff3 --reconciled helixforge.gff3 \\
      --out before_after.md
\b
  # add the Helixer-support column (needs the HDF5)
  helixforge stats --helixer helixer.gff3 --reconciled helixforge.gff3 \\
      --helixer-h5 helixer.h5 --out before_after.md
"""
)
@click.option(
    "--helixer",
    "helixer_gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Helixer GFF3 (before).",
)
@click.option(
    "--reconciled",
    "reconciled_gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="HelixForge GFF3 (after).",
)
@click.option(
    "--genome",
    "genome_fasta",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Genome FASTA (optional).",
)
@click.option(
    "--helixer-h5",
    "helixer_h5",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Helixer HDF5 → populates the Helixer-support before/after column.",
)
@click.option(
    "--out", "out_path", default="before_after.md", help="Markdown summary output path."
)
def stats(
    helixer_gff3: str,
    reconciled_gff3: str,
    genome_fasta: str | None,
    helixer_h5: str | None,
    out_path: str,
) -> None:
    """Inspect (read-only): before/after table of a Helixer vs HelixForge GFF3.

    Works directly off two GFF3 paths (no rerun, no Mikado): it summarises the
    Helixer input and the reconciled output and tabulates the change. Pass
    ``--helixer-h5`` to add the mean-Helixer-support row. Use this to show what
    reconciliation did; use ``benchmark`` for external accuracy metrics against a
    reference / BUSCO / OMArk.

    \b
    Output (--out Markdown), one row per metric, columns:
    metric, helixer, helixforge, delta. Metric rows include gene and transcript
    counts, coding genes, isoforms/gene (mean/max), mono/multi-exon genes, CDS
    length (mean/median), % complete ORFs, mean AED, and mean Helixer support.
    """
    console = _console()
    df = before_after_table(
        helixer_gff3, reconciled_gff3, genome=None, h5_path=helixer_h5
    )
    write_summary(df, out_path)
    console.print(f"[green]Wrote[/green] before/after summary → {out_path}")


@main.command(
    epilog="""\b
Examples:
  # static SVG of the 50 most AS-complex loci
  helixforge viz --genome genome.fa --helixer helixer.gff3 \\
      --helixer-h5 helixer.h5 --stringtie-list stringtie.list \\
      --out-dir figures --mode static --top-n 50
\b
  # one gene, as a PDF
  helixforge viz --genome genome.fa --helixer helixer.gff3 \\
      --out-dir figures --mode static --gene HFG_00042 --fmt pdf
\b
  # BED12 browser track for the whole run
  helixforge viz --genome genome.fa --helixer helixer.gff3 \\
      --out-dir tracks --mode tracks
"""
)
@pipeline_options
@click.option(
    "--out-dir",
    "out_dir",
    required=True,
    type=click.Path(file_okay=False),
    help="Directory for the figures / tracks.",
)
@click.option(
    "--mode",
    "mode",
    type=click.Choice(["static", "interactive", "tracks"]),
    default="static",
    help="Output kind.",
)
@click.option("--gene", "gene_id", default=None, help="Restrict to a single gene id.")
@click.option(
    "--top-n",
    "top_n",
    default=None,
    type=int,
    help="Static mode: only the N most AS-complex loci.",
)
@click.option("--fmt", "fmt", default="svg", help="Static figure format (svg/pdf/png).")
def viz(
    out_dir: str,
    mode: str,
    gene_id: str | None,
    top_n: int | None,
    fmt: str,
    **kwargs: Any,
) -> None:
    """Render per-locus plots, interactive pages, or browser tracks for a run.

    Visualizes a reconciled run. Because the rich ``ReconciledGene`` set is not
    serialised to disk, ``viz`` takes the same inputs as ``reconcile`` and
    rebuilds the gene set via the pipeline before plotting — so it needs the full
    pipeline option set, not a finished GFF3.

    \b
    Output (written under --out-dir), by --mode:
    - static:      one figure per locus (--fmt svg/pdf/png); --top-n keeps the N
                   most AS-complex loci, --gene restricts to a single gene id
    - interactive: an HTML index of interactive per-locus pages
    - tracks:      a BED12 browser track (<output-prefix>.bed)
    """
    import os

    console = _console()
    config = _build_config(kwargs)
    genes = run_pipeline(config)

    if gene_id:
        genes = [g for g in genes if g.gene_id == gene_id]
        if not genes:
            raise click.ClickException(
                f"gene id {gene_id!r} not found in the reconciled set"
            )

    os.makedirs(out_dir, exist_ok=True)
    if mode == "static":
        if gene_id:
            plot_locus(genes[0], out_path=os.path.join(out_dir, f"{gene_id}.{fmt}"))
        else:
            plot_loci(genes, out_dir, top_n_by_as_complexity=top_n, fmt=fmt)
    elif mode == "interactive":
        interactive_index(genes, out_dir)
    else:  # tracks
        write_bed12(genes, os.path.join(out_dir, f"{config.output_prefix}.bed"))
    console.print(f"[green]Wrote[/green] {mode} output → {out_dir}")


@main.group(hidden=True)
def benchmark() -> None:
    """Benchmark an annotation (all tools) or run the ablation figure engine.

    Developer / paper-figure tooling — hidden from the main command list (it is
    not part of the annotation workflow), but still fully runnable.
    """


@benchmark.command(
    "all",
    epilog="""\b
Examples:
  # structural + completeness metrics (whatever tools are installed)
  helixforge benchmark all --reconciled helixforge.gff3 \\
      --proteins proteins.fa --out-dir bench_out
\b
  # add a trusted reference (gffcompare/mikado compare) + BUSCO lineage + OMArk
  helixforge benchmark all --reconciled helixforge.gff3 --proteins proteins.fa \\
      --reference reference.gff3 --lineage eudicots_odb10 --omadb LUCA.h5 \\
      --out-dir bench_out --threads 16
""",
)
@click.option(
    "--reconciled",
    "reconciled_gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--proteins",
    "proteins_fa",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option("--out-dir", "out_dir", required=True, type=click.Path(file_okay=False))
@click.option(
    "--reference",
    "reference_gff3",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Trusted reference GFF3.",
)
@click.option("--lineage", "lineage", default=None, help="compleasm/BUSCO lineage.")
@click.option(
    "--omadb",
    "omadb",
    default=None,
    type=click.Path(exists=True),
    help="OMArk OMA database.",
)
@click.option("--threads", "threads", default=4, type=int)
def benchmark_all_cmd(
    reconciled_gff3: str,
    proteins_fa: str,
    out_dir: str,
    reference_gff3: str | None,
    lineage: str | None,
    omadb: str | None,
    threads: int,
) -> None:
    """Run the applicable benchmarking tools → one comparable TSV table.

    Runs each external benchmarking tool whose inputs are present (AGAT structural
    counts; gffcompare / ``mikado compare`` against ``--reference``; BUSCO or
    compleasm for ``--lineage``; OMArk for ``--omadb``) and flattens every result
    into one long table. Unlike ``stats`` (which compares before/after intrinsic
    counts) this reports *external* accuracy/completeness metrics.

    \b
    Output: <out-dir>/benchmark.tsv, one row per metric, columns:
    tool, metric, value, status.
    """
    import os

    console = _console()
    df = benchmark_all(
        reconciled_gff3,
        proteins_fa,
        out_dir,
        reference_gff3=reference_gff3,
        lineage=lineage,
        omadb=omadb,
        threads=threads,
    )
    os.makedirs(out_dir, exist_ok=True)
    out_tsv = os.path.join(out_dir, "benchmark.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    console.print(f"[green]Wrote[/green] {len(df)} metrics → {out_tsv}")


@benchmark.command(
    "ablation",
    epilog="""\b
Examples:
  # default variant set, benchmarked against a reference
  helixforge benchmark ablation --genome genome.fa --helixer helixer.gff3 \\
      --helixer-h5 helixer.h5 --stringtie-list stringtie.list \\
      --reference reference.gff3 --out-dir ablation_out
\b
  # only the Helixer-coupling levers
  helixforge benchmark ablation --genome genome.fa --helixer helixer.gff3 \\
      --helixer-h5 helixer.h5 --stringtie-list stringtie.list \\
      --variants full,no_helixer_support,no_reference_flag --out-dir ablation_out
""",
)
@pipeline_options
@click.option("--out-dir", "out_dir", required=True, type=click.Path(file_okay=False))
@click.option(
    "--variants",
    "variants",
    default="full,no_helixer_support,no_reference_flag,no_pad,strict_vs_permissive",
    help="Comma-separated ablation variants.",
)
@click.option(
    "--reference",
    "reference_gff3",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option("--lineage", "lineage", default=None)
@click.option("--omadb", "omadb", default=None, type=click.Path(exists=True))
def ablation_cmd(
    out_dir: str,
    variants: str,
    reference_gff3: str | None,
    lineage: str | None,
    omadb: str | None,
    **kwargs: Any,
) -> None:
    """Rerun the pipeline under each Helixer-coupling lever → ablation table TSV.

    Reruns the full reconciliation once per ``--variants`` entry (e.g. ``full``,
    ``no_helixer_support``, ``no_reference_flag``, ``no_pad``,
    ``strict_vs_permissive``), benchmarks each output, and tabulates the metrics
    side by side — the engine behind the manuscript ablation figure. Takes the
    full ``reconcile`` option set since it re-runs the pipeline for each variant.

    \b
    Output: <out-dir>/ablation.tsv, one row per variant, columns:
    variant, num_genes, and one flattened benchmark column per metric
    (``<tool>.<metric>``, e.g. agat.gene_count, busco.complete).
    """
    import os

    console = _console()
    config = _build_config(kwargs)
    variant_list = [v.strip() for v in variants.split(",") if v.strip()]
    df = run_ablation(
        config,
        variant_list,
        out_dir,
        reference_gff3=reference_gff3,
        lineage=lineage,
        omadb=omadb,
    )
    os.makedirs(out_dir, exist_ok=True)
    out_tsv = os.path.join(out_dir, "ablation.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    console.print(
        f"[green]Wrote[/green] ablation table ({len(df)} variants) → {out_tsv}"
    )


@main.group()
def parallel() -> None:
    """Annotate at scale: chunk the genome into regions for distributed reconcile.

    The scheduler-driven scatter-gather primitives (plan / tasks / aggregate /
    suggest): partition the genome into gene-safe chunks, emit per-chunk
    ``reconcile`` tasks for any executor, and merge the outputs into one
    verified annotation. For a single-host run, ``reconcile --scatter`` does the
    same plan→run→aggregate in one command.
    """


@parallel.command(
    "plan",
    epilog="""\b
Examples:
  # one chunk per scaffold (default), stable ids seeded from a prior run
  helixforge parallel plan --genome maize.fa --gff helixer.gff3 \\
      --id-map master_id_map.json -o plan.json
\b
  # ~2000 loci per chunk
  helixforge parallel plan --genome maize.fa --gff helixer.gff3 \\
      --strategy genes --chunk-size 2000 -o plan.json
""",
)
@click.option(
    "--genome",
    "genome_fai",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Genome FASTA or .fai index (for scaffold sizes).",
)
@click.option(
    "--helixer",
    "--gff",
    "helixer_gff3",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Helixer GFF3 (master loci; drives gene-aware cuts + id ranges).",
)
@click.option(
    "--strategy",
    "strategy",
    type=click.Choice(["scaffold", "size", "genes", "adaptive"]),
    default="scaffold",
    show_default=True,
    help="Chunking strategy (v1 vocabulary).",
)
@click.option(
    "--chunk-size",
    "chunk_size",
    default=None,
    type=int,
    help="Bases per chunk for 'size'; loci per chunk for 'genes'.",
)
@click.option(
    "--min-chunk-size",
    "min_chunk_size",
    default=100_000,
    type=int,
    show_default=True,
    help="Smallest chunk (bp) for the size-based strategies.",
)
@click.option(
    "--max-chunk-size",
    "max_chunk_size",
    default=None,
    type=int,
    help="Split a scaffold longer than this (bp) under 'scaffold'.",
)
@click.option(
    "--target-chunks",
    "target_chunks",
    default=None,
    type=int,
    help="Target chunk count for 'adaptive' (from `parallel suggest`).",
)
@click.option(
    "--min-boundary-gap",
    "min_boundary_gap",
    default=None,
    type=int,
    help="Only cut in inter-locus gaps at least this wide (bp), so a "
    "gene is never split. Default max(--flank, 1000).",
)
@click.option("--flank", "flank", default=200, type=int)
@click.option(
    "--id-map",
    "id_map_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Master id_map.json to seed stable HFG ranges from.",
)
@click.option(
    "-o", "--output", "out_path", default="plan.json", help="Output plan path."
)
def parallel_plan(
    genome_fai: str,
    helixer_gff3: str,
    strategy: str,
    chunk_size: int | None,
    min_chunk_size: int,
    max_chunk_size: int | None,
    target_chunks: int | None,
    min_boundary_gap: int | None,
    flank: int,
    id_map_path: str | None,
    out_path: str,
) -> None:
    """Partition the genome (v1 strategies) + reserve disjoint HFG ranges → plan.json.

    Step 1 of the chunked path (plan → tasks → run → aggregate). Cuts the genome
    only in inter-locus gaps wider than ``--min-boundary-gap`` so no gene — and no
    Mikado merge — is ever split, then reserves each chunk a disjoint, contiguous
    HFG number range (seeded from ``--id-map`` for run-stable ids). The plan is
    consumed by ``parallel tasks``.

    \b
    Output: plan.json with top-level {strategy, min_boundary_gap, flank,
    total_loci} and a "chunks" list; each chunk has: chunk_id, regions,
    num_loci, locus_ids, id_base, id_range, novel_base, novel_range,
    est_resources.
    """
    import json

    from helixforge.parallel.plan import (
        partition_by_strategy,
        reserve_id_ranges,
        write_plan,
    )

    console = _console()
    plan = partition_by_strategy(
        genome_fai,
        helixer_gff3,
        strategy=strategy,
        chunk_size=chunk_size,
        min_chunk_size=min_chunk_size,
        max_chunk_size=max_chunk_size,
        target_chunks=target_chunks,
        min_boundary_gap=min_boundary_gap,
        flank=flank,
    )
    id_map = json.loads(open(id_map_path).read()) if id_map_path else None
    reserve_id_ranges(plan, id_map=id_map)
    write_plan(plan, out_path)
    console.print(
        f"[green]Wrote[/green] {len(plan)} chunks ({plan.total_loci} loci, "
        f"strategy={plan.strategy}) → {out_path}\n"
        f"  next: [cyan]helixforge parallel tasks --plan {out_path} "
        f"--genome … --helixer … -o tasks.txt[/cyan]"
    )


@parallel.command(
    "tasks",
    epilog="""\b
Examples:
  # default template (per-chunk reconcile from the attached inputs)
  helixforge parallel tasks --plan plan.json \\
      --genome maize.fa --helixer helixer.gff3 --helixer-h5 preds.h5 \\
      --stringtie-list stringtie.list -o tasks.txt
  parallel -j 16 < tasks.txt          # or: hs cluster tasks.txt --num-tasks 16
\b
  # full override + cluster setup wrapper
  helixforge parallel tasks --plan plan.json --genome g.fa --helixer h.gff3 \\
      --command 'helixforge reconcile --genome g.fa --helixer h.gff3 \\
                 --region {region} --id-base {id_start} \\
                 --output-prefix {output_dir}/{chunk_id}' \\
      --wrapper run_chunk.sh --wrapper-setup 'module load helixforge' -o tasks.txt
""",
)
@pipeline_options
@click.option(
    "--plan",
    "plan_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="plan.json from `parallel plan`.",
)
@click.option(
    "--command",
    "command",
    default=None,
    help="Command template with chunk placeholders ({chunk_id} {region} "
    "{seqid} {start} {end} {start_0} {end_0} {size} {id_start} "
    "{novel_start} {output_dir}). Default: a per-chunk "
    "`helixforge reconcile` built from the attached inputs.",
)
@click.option(
    "-o",
    "--output",
    "out_path",
    default="tasks.txt",
    help="Task file to write (one command per chunk).",
)
@click.option(
    "--output-dir",
    "output_dir",
    default="chunks",
    help="Per-chunk output directory ({output_dir} placeholder).",
)
@click.option(
    "--wrapper",
    "wrapper",
    default=None,
    help="Also write a reusable per-task wrapper script (setup once).",
)
@click.option(
    "--wrapper-setup",
    "wrapper_setup",
    multiple=True,
    help="Setup line for the wrapper (module load / conda activate); repeatable.",
)
@click.option(
    "--include-logging",
    "include_logging",
    is_flag=True,
    help="Redirect each task's stdout/stderr to <output-dir>/logs/<chunk>.log.",
)
@click.option(
    "--helixforge-bin",
    "helixforge_bin",
    default="helixforge",
    help="helixforge binary used in the default command template.",
)
def parallel_tasks(
    plan_path: str,
    command: str | None,
    out_path: str,
    output_dir: str,
    wrapper: str | None,
    wrapper_setup: tuple[str, ...],
    include_logging: bool,
    helixforge_bin: str,
    **kwargs: Any,
) -> None:
    """Expand a command template over the plan → an executor-agnostic task file.

    Step 2 of the chunked path. Expands a command template over every chunk in
    plan.json; the default template is a per-chunk ``helixforge reconcile`` built
    from the inputs you attach, so it works out of the box. Run the resulting file
    with any executor — GNU parallel, a Slurm array, xargs, HyperShell — then
    ``parallel aggregate`` the per-chunk outputs.

    \b
    Template placeholders, substituted per chunk:
    {chunk_id} {region} {seqid} {start} {end} {start_0} {end_0} {size}
    {id_start} {novel_start} {output_dir}
    ({region}/{start_0}/{end_0} are internal 0-based; {seqid}/{start}/{end} are
    1-based inclusive for the CLI.)

    \b
    Output:
    - -o/--output: the task file, one command per chunk (one line each)
    - --wrapper: an optional reusable per-task wrapper script (with
      --wrapper-setup lines) so cluster setup runs once per task
    """
    from helixforge.parallel.plan import read_plan
    from helixforge.parallel.taskgen import generate_task_file
    from helixforge.parallel.tasks import default_reconcile_template

    console = _console()
    base_config = _build_config(kwargs)
    plan = read_plan(plan_path)
    template = command or default_reconcile_template(
        base_config, helixforge_bin=helixforge_bin
    )
    console.print(f"[bold]Command template[/bold]:\n  [cyan]{template}[/cyan]")

    task_file = generate_task_file(
        plan,
        template,
        out_path,
        output_dir=output_dir,
        include_logging=include_logging,
        wrapper=wrapper,
        wrapper_setup=tuple(wrapper_setup),
    )
    console.print(
        f"[green]Wrote[/green] {task_file.n_tasks} tasks → {out_path}"
        + (f" (wrapper {task_file.wrapper_path})" if task_file.wrapper_path else "")
    )
    for i, line in enumerate(task_file.preview(3), 1):
        console.print(f"  {i}. {line[:100]}{'…' if len(line) > 100 else ''}")
    console.print(
        "  run it: [cyan]parallel -j 16 < "
        + out_path
        + "[/cyan]  (or `hs cluster`, `xargs`, a Slurm array)"
    )


@parallel.command(
    "aggregate",
    epilog="""\b
Example:
  helixforge parallel aggregate --input-dir chunks/ --pattern '*.gff3' \\
      -o maize_helixforge.gff3 --master-id-map master_id_map.json
""",
)
@click.option(
    "--input-dir",
    "input_dir",
    default=None,
    type=click.Path(exists=True, file_okay=False),
    help="Directory of per-chunk outputs (use with --pattern).",
)
@click.option(
    "--pattern",
    "pattern",
    default="*.gff3",
    show_default=True,
    help="Glob for the per-chunk GFF3s under --input-dir.",
)
@click.option(
    "--chunk-output",
    "chunk_outputs",
    multiple=True,
    help="Explicit per-chunk output prefix (repeatable; alternative to --input-dir).",
)
@click.option(
    "-o",
    "--output",
    "output",
    default=None,
    help="Merged GFF3 path, e.g. combined.gff3 ('.gff3' stripped for the prefix).",
)
@click.option(
    "--out-prefix",
    "out_prefix",
    default=None,
    help="Output prefix (alternative to -o/--output).",
)
@click.option(
    "--master-id-map",
    "master_id_map_path",
    default=None,
    help="Master id_map.json to fold chunk maps into (created if absent).",
)
def parallel_aggregate(
    input_dir: str | None,
    pattern: str,
    chunk_outputs: tuple[str, ...],
    output: str | None,
    out_prefix: str | None,
    master_id_map_path: str | None,
) -> None:
    """Collect per-chunk outputs by pattern → one annotation; verify + fail closed.

    Step 4 of the chunked path. Merges the per-chunk outputs (by ``--input-dir`` +
    ``--pattern``, or explicit ``--chunk-output`` prefixes) into one genome-wide
    annotation and folds the chunk id_maps into ``--master-id-map``. The merge is
    verifying: it asserts globally-unique gene AND transcript ids and that every
    Helixer locus is covered exactly once, raising on any violation.

    \b
    Outputs (under the resolved prefix):
    - <prefix>.gff3 (+ .tier1/2/3.gff3) and <prefix>.report.tsv — the merged
      annotation + per-gene report (same columns as ``reconcile``)
    - the master id_map.json (created if absent)
    - a printed summary: num_genes, num_loci, tier_counts, origin_counts
    """
    from helixforge.parallel.aggregate import aggregate, prefixes_from_pattern

    console = _console()
    if input_dir:
        prefixes = prefixes_from_pattern(input_dir, pattern)
    elif chunk_outputs:
        prefixes = list(chunk_outputs)
    else:
        raise click.UsageError(
            "give --input-dir (+ --pattern) or --chunk-output prefixes."
        )

    if out_prefix is None:
        if output is None:
            raise click.UsageError(
                "give -o/--output (e.g. combined.gff3) or --out-prefix."
            )
        out_prefix = output[: -len(".gff3")] if output.endswith(".gff3") else output

    result = aggregate(prefixes, out_prefix, master_id_map_path=master_id_map_path)
    console.print(
        f"[green]Aggregated[/green] {len(prefixes)} chunks → "
        f"{result.num_genes} genes, {result.num_loci} loci\n"
        f"  GFF3: {result.gff3_path}  report: {result.report_path}\n"
        f"  tier={result.tier_counts}  origin={result.origin_counts}"
    )


@parallel.command(
    "suggest",
    epilog="""\b
Examples:
  # default node profile
  helixforge parallel suggest --genome maize.fa
\b
  # tune to a specific node + array/walltime caps
  helixforge parallel suggest --genome maize.fa --cores-per-node 32 \\
      --mem-per-node 128 --max-array-size 500 --walltime-cap 12
""",
)
@click.option(
    "--genome",
    "genome_fai",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Genome FASTA or .fai index.",
)
@click.option("--cores-per-node", "cores_per_node", default=16, type=int)
@click.option(
    "--mem-per-node",
    "mem_gb_per_node",
    default=64,
    type=float,
    help="RAM per node, GB.",
)
@click.option("--max-array-size", "max_array_size", default=1000, type=int)
@click.option(
    "--walltime-cap",
    "walltime_cap_hours",
    default=24,
    type=float,
    help="Walltime cap, hours.",
)
def parallel_suggest(
    genome_fai: str,
    cores_per_node: int,
    mem_gb_per_node: float,
    max_array_size: int,
    walltime_cap_hours: float,
) -> None:
    """Recommend chunk count + per-chunk resources (heuristic; prints trade-offs).

    Step 0 of the chunked path. From the genome size/scaffold profile and your
    node spec (``--cores-per-node``, ``--mem-per-node``, ``--max-array-size``,
    ``--walltime-cap``) it suggests a granularity and per-chunk resources to feed
    into ``parallel plan`` / ``run --scatter``. Heuristic only — it reads sizes,
    runs nothing, and modifies nothing.

    \b
    Output (printed recommendation): target_chunks, target_loci_per_chunk,
    min_boundary_gap, and per-chunk procs / threads / mem_gb / walltime_min,
    with a rationale explaining the trade-offs.
    """
    from helixforge.parallel.suggest import suggest_plan

    console = _console()
    suggestion = suggest_plan(
        genome_fai,
        hpc_profile={
            "cores_per_node": cores_per_node,
            "mem_gb_per_node": mem_gb_per_node,
            "max_array_size": max_array_size,
            "walltime_cap_hours": walltime_cap_hours,
        },
    )
    console.print(suggestion.render())


@parallel.command(
    "example-sbatch",
    epilog="""\b
Examples:
  # HyperShell wrapper (default)
  helixforge parallel example-sbatch -o run.sbatch
\b
  # GNU Parallel wrapper
  helixforge parallel example-sbatch -o run.sbatch --executor parallel
""",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False),
    help="Where to write the example SBATCH script.",
)
@click.option(
    "--executor",
    type=click.Choice(["hypershell", "parallel"]),
    default="hypershell",
    show_default=True,
    help="Run the task file with HyperShell or GNU Parallel.",
)
def parallel_example_sbatch(output_path: str, executor: str) -> None:
    """Write a copy-paste SBATCH wrapper that runs a ``tasks`` file on one node.

    Convenience starting point for the scheduler-agnostic path: it runs the whole
    ``tasks.txt`` (from ``parallel tasks``) inside a single allocation, fanning
    out across ``$SLURM_CPUS_PER_TASK`` cores via HyperShell or GNU Parallel. Edit
    the partition / time / resources for your cluster. For a true one-task-per-
    chunk Slurm *array* driven off a plan, use ``reconcile --hpc slurm`` instead.
    """
    from helixforge.parallel.taskgen import write_example_sbatch

    path = write_example_sbatch(output_path, executor=executor)
    _console().print(f"Wrote example sbatch ({executor}) → {path}")


from helixforge.cli.utils import utils as utils_group  # noqa: E402

main.add_command(utils_group, "utils")


if __name__ == "__main__":
    main()
