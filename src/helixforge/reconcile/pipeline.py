"""Pipeline orchestration: the HelixForge Python API."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import attrs

from helixforge.constants import (
    LNCRNA_MIN_LENGTH,
    MERGE_MIN_GAP_READS,
    NONCODING_MAX_CDS_CONF,
    NOVEL_ID_BASE,
    PARALOG_IDENTITY_K,
)
from helixforge.io.bam import parse_star_sj_tab
from helixforge.io.fasta import GenomeAccessor
from helixforge.io.gff import GFF3Writer
from helixforge.io.hdf5 import HDF5ConfidenceReader, open_confidence_reader
from helixforge.io.miniprot import MiniprotParser
from helixforge.io.stringtie import (
    StringTieParser,
    best_overlapping_tpm,
    tpm_overlap_index_from_transcripts,
)
from helixforge.mikado.config import (
    install_scoring_profile,
    validate_scoring_profile,
    write_configuration,
    write_input_list,
)
from helixforge.mikado.emit_external import (
    helixer_support_and_conf,
    normalize_tpm,
    write_external_scores_tsv,
)
from helixforge.mikado.emit_gtf import helixer_gff3_to_gtf, stringtie_to_labelled_gtf
from helixforge.mikado.emit_junctions import junctions_to_portcullis_tab, run_portcullis
from helixforge.mikado.parse import parse_loci_gff3
from helixforge.mikado.run import (
    run_diamond,
    run_pick,
    run_prepare,
    run_serialise,
    run_transdecoder,
)
from helixforge.qc.flags import dedup_flags
from helixforge.reconcile.cds import (
    assign_backstop_cds,
    batch_backstop_transdecoder,
    cds_cross_check,
)
from helixforge.reconcile.classify import classify_loci
from helixforge.reconcile.fallback import JunctionIndex, refine_backstop_gene
from helixforge.reconcile.pseudogene import apply_pseudogene_typing
from helixforge.reconcile.biotype import assign_biotype
from helixforge.reconcile.locus import load_helixer_loci
from helixforge.reconcile.mikado_integrate import IdAllocator, reconcile
from helixforge.reconcile.models import Exon, SpliceJunction, _validate_fraction
from helixforge.reconcile.runstats import RunStats
from helixforge.reconcile.validate import validate_all
from helixforge.utils.atomic import atomic_write
from helixforge.utils.checkpoint import Checkpoint
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.reconcile.models import (
        HelixerLocus,
        LocusClassification,
        MikadoLocus,
        MiniprotAlignment,
        QCFlag,
        ReconciledGene,
        TranscriptCandidate,
    )

_log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
#
# The ~60-field flat ``PipelineConfig`` is grouped
# into five validated sub-configs (``ClassificationConfig``, ``ASConfig``,
# ``ReconcileConfig``, ``ValidationConfig``, ``ResourceConfig``). To keep every
# existing construction path working byte-for-byte: ``PipelineConfig(min_tpm=…)``
# in tests, ``dataclasses.replace(cfg, pad=…)`` in the ablation runner / chunk
# builder, and the CLI's flat-kwargs mapping, **the flat fields stay on
# PipelineConfig**. The sub-configs are exposed as validated *views* (read
# properties), their per-field validators run at construction via
# ``__post_init__``, and ``from_flat`` / the YAML round-trip (D2) move between the
# flat and nested representations. So the grouping adds validation + a nested
# config surface without breaking a single flat caller.


def _require_nonneg(value: float | int | None, name: str) -> None:
    if value is not None and value < 0:
        raise ValueError(f"{name} must be >= 0, got {value}")


def _require_positive_int(value: float | int, name: str) -> None:
    if int(value) < 1:
        raise ValueError(f"{name} must be >= 1, got {value}")


@dataclass
class ClassificationConfig:
    """Expression/coverage classification thresholds (Phase 3)."""

    min_tpm: float = 0.5
    min_samples: int = 1
    coverage_threshold: float = 2.0
    near_zero_coverage: float = 0.1

    def __post_init__(self) -> None:
        _require_nonneg(self.min_tpm, "min_tpm")
        _require_positive_int(self.min_samples, "min_samples")
        _require_nonneg(self.coverage_threshold, "coverage_threshold")
        _validate_fraction(self.near_zero_coverage, "near_zero_coverage")


@dataclass
class ASConfig:
    """Mikado pick / alternative-splicing knobs."""

    as_report: bool = True
    only_confirmed_introns: bool = True
    max_isoforms: int = 5
    keep_retained_introns: bool = False
    pad: bool = True
    chimera_split: bool = True
    flank: int = 200

    def __post_init__(self) -> None:
        _require_positive_int(self.max_isoforms, "max_isoforms")
        _require_nonneg(self.flank, "flank")


@dataclass
class ReconcileConfig:
    """Locus-correspondence + isoform-admission thresholds.

    Paralog/tandem-array merge guards:
    ``merge_min_gap_reads`` / ``merge_require_canonical`` gate an accepted merge on
    a canonical, sufficiently-read-supported bridging junction; the optional
    ``paralog_identity_threshold`` (None = disabled) rejects a merge of two
    high-k-mer-identity adjacent loci (a recent-duplication / homeolog signal).
    """

    reciprocal_overlap: float = 0.5
    min_cds_overlap: float = 0.6
    min_cdna_overlap: float = 0.6
    admit_novel: bool = False
    novel_evidence_floor: float | None = None
    merge_min_gap_reads: int = MERGE_MIN_GAP_READS
    merge_require_canonical: bool = True
    paralog_identity_threshold: float | None = None
    paralog_kmer_k: int = PARALOG_IDENTITY_K

    def __post_init__(self) -> None:
        _validate_fraction(self.reciprocal_overlap, "reciprocal_overlap")
        _validate_fraction(self.min_cds_overlap, "min_cds_overlap")
        _validate_fraction(self.min_cdna_overlap, "min_cdna_overlap")
        _require_nonneg(self.merge_min_gap_reads, "merge_min_gap_reads")
        _require_positive_int(self.paralog_kmer_k, "paralog_kmer_k")
        if self.paralog_identity_threshold is not None:
            _validate_fraction(
                self.paralog_identity_threshold, "paralog_identity_threshold"
            )


@dataclass
class ValidationConfig:
    """Structural codon-gate + backstop junction-correction thresholds (Phase 7)."""

    short_cds_threshold: int = 300
    short_exon_threshold: int = 10
    long_intron_threshold: int = 100_000
    junction_tolerance: int = 0
    junction_min_reads: int = 3
    # blast_score floor for treating a disabled, homology-backed ORF
    # as a pseudogene candidate. 0.0 = any positive homology.
    pseudogene_min_homology: float = 0.0
    # lncRNA minimum spliced-transcript length (nt)
    # and the Helixer CDS-channel confidence ceiling for a non-coding call.
    lncrna_min_length: int = LNCRNA_MIN_LENGTH
    noncoding_max_cds_conf: float = NONCODING_MAX_CDS_CONF

    def __post_init__(self) -> None:
        for name in (
            "short_cds_threshold",
            "short_exon_threshold",
            "long_intron_threshold",
            "junction_tolerance",
            "junction_min_reads",
            "pseudogene_min_homology",
            "lncrna_min_length",
            "noncoding_max_cds_conf",
        ):
            _require_nonneg(getattr(self, name), name)
        _validate_fraction(self.noncoding_max_cds_conf, "noncoding_max_cds_conf")


@dataclass
class ResourceConfig:
    """Process / thread budget for external tool invocations."""

    procs: int = 1
    threads: int = 4

    def __post_init__(self) -> None:
        _require_positive_int(self.procs, "procs")
        _require_positive_int(self.threads, "threads")


# (attr-name on PipelineConfig view → sub-config class). The order is the nested
# YAML group order. Used by from_flat / to_nested_dict / __post_init__ validation.
_SUBCONFIG_VIEWS: tuple[
    tuple[
        str,
        type[
            ClassificationConfig
            | ASConfig
            | ReconcileConfig
            | ValidationConfig
            | ResourceConfig
        ],
    ],
    ...,
] = (
    ("classification", ClassificationConfig),
    ("as_config", ASConfig),
    ("reconciliation", ReconcileConfig),
    ("validation", ValidationConfig),
    ("resources", ResourceConfig),
)


@dataclass
class PipelineConfig:
    """All inputs, tool paths, and tuning knobs for one ``run_pipeline`` call.

    Defaults mirror Mikado's recommended settings and the downstream functions' own defaults
    (classification, AS knobs, reconciliation, validation), they are *reused*,
    not reinvented. Only ``genome_fasta`` and ``helixer_gff3`` are required.
    """

    # --- required inputs ---
    genome_fasta: str
    helixer_gff3: str

    # --- optional evidence inputs ---
    helixer_h5: str | None = None
    # Metadata half (*_input.h5) when helixer_h5 is a bare *_predictions.h5 and the
    # sibling cannot be auto-located. Mirrors the `confidence` command's --input-h5.
    helixer_input_h5: str | None = None
    stringtie_list: list[str] = field(default_factory=list)  # StringTie GTF paths
    bam_paths: list[str] = field(default_factory=list)
    bigwig_paths: list[str] = field(default_factory=list)
    star_sj_paths: list[str] = field(default_factory=list)
    miniprot_gff: str | None = None
    protein_db: str | None = None

    # --- external tool binaries ---
    transdecoder_bin_dir: str | None = None
    # Per-backstop-gene TransDecoder ORF fallback (cds.py). Off by default: on a
    # real genome it spawns one TransDecoder subprocess per CDS-less backstop
    # gene (hundreds), which is rarely worth it, miniprot is the primary source.
    backstop_transdecoder: bool = False
    portcullis_bin: str | None = None
    mikado_bin: str = "mikado"
    diamond_bin: str = "diamond"
    use_mikado_configure: bool = True  # real `mikado configure` vs. templated YAML

    # --- scoring ---
    scoring_profile: str = "strict"

    # --- Helixer coupling levers (the project's novel contribution).
    #     Defaults reproduce the standard run; the ablation
    #     engine toggles these to quantify each lever's contribution.
    #     ``helixer_is_reference`` drops Helixer's ``is_reference=true`` flag.
    #     ``helixer_support_weight`` scales the external helixer_support /
    #     helixer_locus_conf metrics (0.0 == external metric weight off). ---
    helixer_is_reference: bool = True
    helixer_support_weight: float = 1.0

    # --- output ---
    output_prefix: str = "helixforge"
    report_path: str | None = None
    id_map_path: str | None = None
    work_dir: str | None = None

    # --- genome-level report. On by default: a run
    #     emits ``<prefix>.report.json`` (MultiQC-compatible) + ``.report.html``
    #     summarising structure / biotype / tier / completeness / evidence
    #     support. Reporting only, count-neutral (never changes the gene set);
    #     best-effort, so a report failure never sinks a finished run. ---
    write_report: bool = True

    # --- GFF3 spec-completeness + provenance. On by default:
    #     the emitted GFF3 carries ``##sequence-region`` directives (one per
    #     contig) and a ``#!`` provenance preamble (HelixForge version, resolved
    #     external-tool versions, the parameter hash, input MD5s) so the output is
    #     self-describing and gt-gff3validator / AGAT clean. Header-only, purely
    #     count-neutral (the parsed gene/tier counts never change). ``embed_fasta``
    #     additionally appends the genome after ``##FASTA`` (off by default; large).
    gff3_sequence_regions: bool = True
    gff3_provenance: bool = True
    gff3_embed_fasta: bool = False

    # --- functional annotation hook. Off by default: an
    #     opt-in InterProScan/eggNOG-mapper stage that attaches GO/Pfam/InterPro
    #     attributes (Ontology_term/Dbxref) + a domain-completeness credibility
    #     flag to the protein_coding set. Needs the external tool on PATH; never
    #     on the default/golden path. ---
    functional_annotation: bool = False
    functional_tool: str = "interproscan"  # interproscan | eggnog | both
    functional_db: str | None = None  # eggNOG data dir (eggnog/both only)
    interproscan_bin: str = "interproscan.sh"
    eggnog_bin: str = "emapper.py"

    # --- region / chunk (single subset per call; no chunk engine here) ---
    region: object | None = None  # "seqid" | "seqid:start-end" | (seqid,start,end)
    chunk_id: str | None = None

    # --- HFG id-range allocation (Phase 16 §C3). Defaults reproduce the global
    #     single-process numbering; ``parallel/plan.py`` reserves a disjoint
    #     range per chunk so a scattered run stays collision-free without a
    #     shared run-time lock. Defaults keep the M3 golden counts/IDs identical.
    id_base: int = 1
    novel_base: int = NOVEL_ID_BASE

    # --- shared-Mikado reuse (Phase 17 ablations). When set to a prior run's
    #     ``mikado_run`` dir, the heavy prepare/external/TransDecoder/DIAMOND/
    #     serialise steps are skipped (their artifacts reused) and only
    #     ``mikado pick`` reruns. Valid only for knobs that do NOT change the
    #     prepared transcripts or the serialise DB (pad / scoring_profile /
    #     AS knobs), the ablation runner sets it solely for those variants.
    reuse_mikado_dir: str | None = None

    # --- stage checkpoint/resume (Phase 22 §3.1). Off by default reproduces the
    #     cold-run behavior and counts exactly (no manifest written). When True,
    #     ``run_pipeline`` records completed stages in ``work_dir/checkpoint.json``
    #     and, on a re-run, skips any stage whose outputs still exist and validate
    #     (the MIKADO loci are re-parsed from disk instead of re-running the
    #     external chain; OUTPUT files are not rewritten). This is the in-package
    #     generalization of ``reuse_mikado_dir`` and the m3_run ``--resume``. ---
    resume: bool = False

    # --- external Mikado loci (skip the entire Mikado stage). When set to a path
    #     to a ``mikado.loci.gff3`` file (from an externally-run Mikado), the PREP
    #     and MIKADO stages are skipped entirely and the pipeline goes straight to
    #     consume/reconcile. Companion files (``*.metrics.tsv``, ``*.scores.tsv``)
    #     are auto-detected in the same directory; a warning is logged if absent.
    mikado_loci_gff3: str | None = None

    # --- resources ---
    procs: int = 1
    threads: int = 4

    # --- intra-node parallel finalize. Number of
    #     worker processes for the embarrassingly-parallel per-gene finalize loop
    #     (junction correction + CDS projection + cross-check), dispatched by
    #     scaffold. Default 1 runs the serial path **byte-identically** (the M3
    #     golden config never sets it), so this is opt-in scaling only. The result
    #     is sorted back into gene order, so it is independent of worker count. ---
    finalize_workers: int = 1

    # --- classification (Phase 3 defaults) ---
    min_tpm: float = 0.5
    min_samples: int = 1
    coverage_threshold: float = 2.0
    near_zero_coverage: float = 0.1

    # --- Mikado pick / AS knobs ---
    as_report: bool = True
    only_confirmed_introns: bool = True
    max_isoforms: int = 5
    keep_retained_introns: bool = False
    pad: bool = True
    chimera_split: bool = True
    flank: int = 200

    # --- reconciliation (Phase 6 defaults) ---
    reciprocal_overlap: float = 0.5
    min_cds_overlap: float = 0.6
    min_cdna_overlap: float = 0.6
    admit_novel: bool = False
    novel_evidence_floor: float | None = None

    # --- TRaCE canonical-transcript election (Phase 33b). Off by default: the
    #     primary isoform of a multi-transcript gene is the highest-combined_score
    #     one (the historical rule), so the M3 golden ids/primary are untouched.
    #     When on, ``reconcile/trace.py`` holds a ranked-choice election (RNA-seq
    #     sample + length voters) to elect the canonical/primary isoform and order
    #     the rest; this reorders + renumbers isoforms only (never adds/drops genes,
    #     changes tiers/origins, or the AS-event set). The Arabidopsis poster + maize
    #     runs enable it. Defaults mirror ``trace.TraceParams`` (the TRaCE paper). ---
    trace_primary: bool = False
    trace_max_aed: float = 0.5
    trace_min_tpm: float = 0.5
    trace_min_overlap: float = 0.5
    trace_weight_domain: float = 9.0
    trace_weight_protein: float = 6.0
    trace_weight_cdna: float = 3.0
    trace_use_domain: bool = True

    # --- paralog/tandem-array merge guards ---
    merge_min_gap_reads: int = MERGE_MIN_GAP_READS
    merge_require_canonical: bool = True
    paralog_identity_threshold: float | None = None
    paralog_kmer_k: int = PARALOG_IDENTITY_K

    # --- structural validation thresholds (Phase 7 defaults) ---
    short_cds_threshold: int = 300
    short_exon_threshold: int = 10
    long_intron_threshold: int = 100_000
    # blast_score floor for pseudogene-candidate typing.
    pseudogene_min_homology: float = 0.0
    # non-coding biotype thresholds.
    lncrna_min_length: int = LNCRNA_MIN_LENGTH
    noncoding_max_cds_conf: float = NONCODING_MAX_CDS_CONF

    # --- genetic code ---
    # NCBI transl_table id (default 1 = standard code; nuclear plant genes).
    # ``transl_table_map`` ({seqid: table_id}) overrides it per-scaffold so
    # organellar (plastid/mito) contigs translate with the right code; an
    # unmapped scaffold keeps the default, so the default run is unchanged.
    transl_table: int = 1
    transl_table_map: dict[str, int] | None = None

    # --- backstop junction correction (Phase 7 defaults) ---
    junction_tolerance: int = 0
    junction_min_reads: int = 3

    # --- structured-ncRNA hook ---
    # Off by default: an opt-in completeness step that runs tRNAscan-SE / Infernal
    # and merges tRNA/rRNA/snoRNA loci into the final GFF3. Needs the external tool
    # on PATH (or an absolute bin); never on the default/golden path.
    ncrna_scan: bool = False
    ncrna_tool: str = "trnascan"
    ncrna_rfam_cm: str | None = None

    # --- experimental VCF / haplotype awareness ---
    # Off by default. An optional **decomposed** VCF (multiallelic sites split via
    # `bcftools norm -m -`) flags genes whose CDS overlaps a high-impact variant
    # (premature stop / splice-disrupting / frameshift) with VARIANT_IMPACTED.
    # Input plumbing + flag only: never re-types/re-tiers; with no vcf_path (the
    # default/golden path) it is a no-op. Full pangenome projection is future work.
    vcf_path: str | None = None
    vcf_flag_impact: bool = False

    # --- optional EDTA TE gating. Off unless ``te_annotation`` is set. EDTA is
    # the only signal that calls a TE: with it, each model's overlap with TE-class
    # features is flagged (TE_OVERLAP) and a good-ORF gene above
    # ``te_overlap_threshold`` is reclassified ``transposable_element``. Only the
    # configured ``te_classes`` (EDTA Classification orders; None = the true-TE
    # default set) count, knob/satellite/centromere/rDNA/low-complexity are
    # excluded. With no te_annotation (the default/golden path) it is a no-op.
    te_annotation: str | None = None
    te_overlap_threshold: float = 0.5
    te_classes: list[str] | None = None

    def __post_init__(self) -> None:
        # Derive output paths from the prefix when not given explicitly.
        if self.report_path is None:
            self.report_path = f"{self.output_prefix}.report.tsv"
        if self.id_map_path is None:
            self.id_map_path = f"{self.output_prefix}.id_map.json"
        if self.work_dir is None:
            base = f"{self.output_prefix}_work"
            self.work_dir = f"{base}/{self.chunk_id}" if self.chunk_id else base
        # A region round-tripped through YAML (D2) comes back as a list; restore
        # the tuple form so an explicit-span region compares equal after reload.
        if isinstance(self.region, list):
            self.region = tuple(self.region)
        # Run the grouped sub-config validators. Building each view
        # invokes its __post_init__, rejecting an out-of-range knob at
        # construction (e.g. min_cds_overlap > 1) regardless of how it was set.
        for name, _klass in _SUBCONFIG_VIEWS:
            getattr(self, name)

    # --- grouped sub-config views -----------------------------
    # Each returns a freshly-built, validated sub-config from the flat fields.
    # Read-only by design: mutate the flat field (or use ``from_flat``); the view
    # is a projection, not stored state.

    @property
    def classification(self) -> ClassificationConfig:
        return ClassificationConfig(
            self.min_tpm,
            self.min_samples,
            self.coverage_threshold,
            self.near_zero_coverage,
        )

    @property
    def as_config(self) -> ASConfig:
        return ASConfig(
            self.as_report,
            self.only_confirmed_introns,
            self.max_isoforms,
            self.keep_retained_introns,
            self.pad,
            self.chimera_split,
            self.flank,
        )

    @property
    def reconciliation(self) -> ReconcileConfig:
        return ReconcileConfig(
            self.reciprocal_overlap,
            self.min_cds_overlap,
            self.min_cdna_overlap,
            self.admit_novel,
            self.novel_evidence_floor,
            self.merge_min_gap_reads,
            self.merge_require_canonical,
            self.paralog_identity_threshold,
            self.paralog_kmer_k,
        )

    @property
    def validation(self) -> ValidationConfig:
        return ValidationConfig(
            self.short_cds_threshold,
            self.short_exon_threshold,
            self.long_intron_threshold,
            self.junction_tolerance,
            self.junction_min_reads,
            self.pseudogene_min_homology,
            self.lncrna_min_length,
            self.noncoding_max_cds_conf,
        )

    @property
    def resources(self) -> ResourceConfig:
        return ResourceConfig(self.procs, self.threads)

    # --- flat ↔ nested bridges ---------------------------

    @classmethod
    def from_flat(cls, **kwargs: Any) -> "PipelineConfig":
        """Build a :class:`PipelineConfig` from flat kwargs **or** nested groups.

        Any key matching a sub-config view name (``classification``, ``as_config``,
        ``reconciliation``, ``validation``, ``resources``) whose value is a dict or
        a sub-config instance is expanded into its flat fields, so this accepts
        both the historical flat CLI kwargs and the nested YAML structure (D2).
        Unknown nested keys flow straight through to ``cls(**flat)``.
        """
        flat: dict[str, Any] = dict(kwargs)
        for name, klass in _SUBCONFIG_VIEWS:
            if name not in flat:
                continue
            sub = flat.pop(name)
            if sub is None:
                continue
            if isinstance(sub, klass):
                import dataclasses as _dc

                sub = _dc.asdict(sub)
            flat.update(sub)
        return cls(**flat)

    def to_nested_dict(self) -> dict[str, Any]:
        """Serialise to a nested dict: flat top-level fields + the five groups."""
        import dataclasses as _dc

        grouped_fields: set[str] = set()
        for _name, klass in _SUBCONFIG_VIEWS:
            grouped_fields.update(f.name for f in _dc.fields(klass))
        out: dict[str, Any] = {}
        for f in _dc.fields(self):
            if f.name in grouped_fields:
                continue
            out[f.name] = getattr(self, f.name)
        for name, _klass in _SUBCONFIG_VIEWS:
            out[name] = _dc.asdict(getattr(self, name))
        # A tuple region cannot be safe-dumped to YAML; emit it as a list. The
        # __post_init__ list→tuple restore makes the reload compare equal.
        if isinstance(out.get("region"), tuple):
            out["region"] = list(out["region"])
        return out

    def to_yaml(self, path: str | Path) -> str | None:
        """Write the fully-resolved nested config to ``path`` (PyYAML optional)."""
        try:
            import yaml
        except ImportError:
            _log.warning(
                "PyYAML unavailable: skipping resolved-config dump to %s", path
            )
            return None
        with atomic_write(str(path)) as fh:
            yaml.safe_dump(
                self.to_nested_dict(), fh, default_flow_style=False, sort_keys=False
            )
        return str(path)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "PipelineConfig":
        """Deserialise a nested run-config YAML into a :class:`PipelineConfig`."""
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - environment-dependent
            raise RuntimeError(
                "PyYAML is required to read a --config run.yaml; install pyyaml "
                "or pass the configuration via CLI flags."
            ) from exc
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        return cls.from_flat(**data)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


class _Timer:
    """Log a stage's start/end with elapsed wall-clock seconds."""

    def __init__(self, stage: str) -> None:
        self.stage = stage

    def __enter__(self) -> "_Timer":
        self._t0 = time.perf_counter()
        _log.info("[%s] start", self.stage)
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> None:
        _log.info("[%s] done in %.2fs", self.stage, time.perf_counter() - self._t0)


def _primary(gene: "ReconciledGene") -> "TranscriptCandidate":
    for t in gene.transcripts:
        if t.transcript_id == gene.primary_transcript_id:
            return t
    return gene.transcripts[0]


def _structure_key(seqid: str, strand: str, exons: list[Exon]) -> str:
    """Match ``StringTieParser._structure_key`` so prepared tx can look up TPM."""
    exon_tuples = tuple((e.start, e.end) for e in exons)
    return f"{seqid}:{strand}:{exon_tuples}"


def _merge_star_junctions(
    sj_paths: list[str], min_reads: int = 3
) -> list[SpliceJunction]:
    """Union STAR ``SJ.out.tab`` junctions across samples (sum reads, count samples).

    Carries the STAR-derived canonical motif class and the summed
    multi-mapping read count through the union, so the backstop
    junction-correction canonicity gate and downstream consumers see them.
    """
    agg: dict[tuple[str, int, int, str], list[Any]] = {}
    for p in sj_paths:
        for j in parse_star_sj_tab(str(p), min_unique_reads=min_reads):
            key = (j.seqid, j.donor, j.acceptor, j.strand)
            if key in agg:
                agg[key][0] += j.read_count
                agg[key][1] += 1
                agg[key][3] += j.multimap_reads
                if agg[key][2] is None:
                    agg[key][2] = j.canonical
            else:
                agg[key] = [j.read_count, 1, j.canonical, j.multimap_reads]
    out = [
        SpliceJunction(
            seqid=s,
            donor=d,
            acceptor=a,
            strand=st,
            read_count=reads,
            samples=samples,
            canonical=canon,
            multimap_reads=multi,
        )
        for (s, d, a, st), (reads, samples, canon, multi) in agg.items()
    ]
    out.sort(key=lambda j: (j.seqid, j.donor, j.acceptor))
    return out


def _parse_prepared_gtf(
    gtf_path: str,
) -> dict[str, tuple[str, str, list[Exon]]]:
    """Minimal reader for ``mikado_prepared.gtf`` → ``{tid: (seqid, strand, [Exon])}``.

    Unlike ``StringTieParser`` this keeps **every** prepared transcript (no TPM
    filter), each needs an external-scores row. Coordinates are converted from
    1-based GTF to internal 0-based half-open here (I/O boundary).
    """
    records: dict[str, dict[str, Any]] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            seqid, strand, attrs_field = cols[0], cols[6], cols[8]
            tid: str | None = None
            for field_str in attrs_field.split(";"):
                field_str = field_str.strip()
                if field_str.startswith("transcript_id"):
                    tid = (
                        field_str.split('"')[1]
                        if '"' in field_str
                        else field_str.split()[-1]
                    )
                    break
            if tid is None:
                continue
            rec = records.setdefault(
                tid, {"seqid": seqid, "strand": strand, "exons": []}
            )
            rec["exons"].append(Exon(int(cols[3]) - 1, int(cols[4])))
    out: dict[str, tuple[str, str, list[Exon]]] = {}
    for tid, rec in records.items():
        exons = sorted(rec["exons"], key=lambda e: e.start)
        out[tid] = (rec["seqid"], rec["strand"], exons)
    return out


def _build_external_scores(
    prepared_gtf: str,
    h5_reader: HDF5ConfidenceReader | None,
    struct_tpm: dict[str, float],
    global_max_tpm: float,
    helixer_weight: float = 1.0,
) -> dict[str, dict[str, float]]:
    """One external-scores row per prepared transcript (helixer + tpm, all in [0,1]).

    When no HDF5 is available the Helixer-derived metrics default to ``0.0``
    (neutral; logged once by the caller) so Mikado still gets a complete TSV.
    ``helixer_weight`` scales the two Helixer metrics (the ablation lever, 0.0
    switches the Helixer↔evidence coupling off; result stays in [0, 1]).
    """
    prepared = _parse_prepared_gtf(prepared_gtf)
    rows: dict[str, dict[str, float]] = {}
    for tid, (seqid, strand, exons) in prepared.items():
        if h5_reader is not None:
            # One HDF5 read per transcript: both
            # Helixer metrics come from a single prediction block. Bit-identical
            # to the per-exon helixer_support / helixer_locus_conf calls.
            support, locus_conf = helixer_support_and_conf(
                exons, seqid, strand, h5_reader
            )
            hs = helixer_weight * support
            hlc = helixer_weight * locus_conf
        else:
            hs = hlc = 0.0
        raw_tpm = struct_tpm.get(_structure_key(seqid, strand, exons), 0.0)
        rows[tid] = {
            "helixer_support": hs,
            "helixer_locus_conf": hlc,
            "tpm": normalize_tpm(raw_tpm, global_max_tpm),
        }
    return rows


def _patch_as_knobs(config_path: str | Path, knobs: dict[str, Any]) -> None:
    """Inject ``pick.alternative_splicing`` knobs into a ``mikado configure`` YAML.

    ``mikado configure`` writes a complete config but with default AS settings;
    the §7 knobs (esp. ``pad``) must be set afterward. Uses PyYAML if available
    (Mikado's own env ships it); warns and skips otherwise.
    """
    try:
        import yaml
    except ImportError:
        _log.warning(
            "PyYAML unavailable: using mikado configure AS defaults; set "
            "pick.alternative_splicing manually for the strict benchmark."
        )
        return
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    cfg.setdefault("pick", {}).setdefault("alternative_splicing", {}).update(knobs)
    with open(config_path, "w") as fh:
        yaml.safe_dump(cfg, fh, default_flow_style=False, sort_keys=False)
    _log.info("patched AS knobs into %s", Path(config_path).name)


def _should_run_mikado(config: PipelineConfig) -> bool:
    """Mikado runs only with expression evidence *and* a protein DB.

    The MIKADO chain includes DIAMOND
    homology and a ``serialise`` that needs ``--blast_targets``; without a
    protein DB or StringTie assemblies there is nothing to reconcile, so the
    pipeline degrades to the Helixer-backstop path (basic mode), every gene is
    carried through as ``helixer_backstop``.
    """
    return bool(config.stringtie_list) and bool(config.protein_db)


# ---------------------------------------------------------------------------
# Stage A: input loading
# ---------------------------------------------------------------------------


def _load_inputs(config: PipelineConfig) -> dict[str, Any]:
    """Open/parse every input; return a dict of in-memory evidence + open handles.

    The ``genome`` and ``h5`` handles are returned **open** and must be closed by
    the caller (``run_pipeline`` does so in ``finally``).
    """
    for label, path in (
        ("genome_fasta", config.genome_fasta),
        ("helixer_gff3", config.helixer_gff3),
    ):
        if not path or not Path(path).exists():
            raise FileNotFoundError(f"required input {label} not found: {path!r}")

    genome = GenomeAccessor(config.genome_fasta)
    # Resolve Helixer's split input/predictions halves the same way `confidence`
    # does: helixer_h5 may be a combined file, a *_input.h5, or a *_predictions.h5
    # (sibling/helixer_input_h5 supplies the metadata partner).
    h5 = (
        open_confidence_reader(config.helixer_h5, input_h5=config.helixer_input_h5)
        if config.helixer_h5
        else None
    )

    # --- Helixer loci (+ HDF5 enrichment, region filter) ---
    region = config.region
    if isinstance(region, str) and ":" not in region:
        loci = load_helixer_loci(
            config.helixer_gff3,
            h5_path=config.helixer_h5,
            input_h5=config.helixer_input_h5,
        )
        loci = [g for g in loci if g.seqid == region]
    else:
        loci = load_helixer_loci(
            config.helixer_gff3,
            h5_path=config.helixer_h5,
            input_h5=config.helixer_input_h5,
            region=region,  # type: ignore[arg-type]  # object|None; load_helixer_loci untyped
        )

    # --- StringTie assemblies ---
    st_parser = StringTieParser()
    per_sample: dict[str, list[Any]] = {}
    all_st: list[Any] = []
    for gtf in config.stringtie_list:
        sid = Path(gtf).stem
        tx = st_parser.parse_gtf(str(gtf), sid)
        per_sample[sid] = tx
        all_st.extend(tx)
    st_agg: dict[str, Any] = (
        st_parser.aggregate_across_samples(all_st) if all_st else {}
    )
    struct_tpm: dict[str, float] = {k: v["max_tpm"] for k, v in st_agg.items()}
    global_max_tpm = max(struct_tpm.values(), default=0.0)
    # Per-model TPM overlap index (min_tpm=0 keeps every structure; the threshold
    # comparison happens later). Built from the already-parsed transcripts so the
    # GTFs are not read a second time. Shared with the `evidence` scorer.
    tpm_overlap_index = (
        tpm_overlap_index_from_transcripts(all_st) if all_st else {}
    )

    # --- splice junctions (STAR SJ union) for serialise + backstop refinement ---
    junctions: list[SpliceJunction] = (
        _merge_star_junctions(config.star_sj_paths, config.junction_min_reads)
        if config.star_sj_paths
        else []
    )
    # Build the junction index once; threaded into
    # backstop refinement so every intron is an O(1)/O(log J) lookup, not a scan
    # of the whole set. Identical match decisions, so the golden gate holds.
    junction_index = JunctionIndex.from_junctions(junctions)

    # --- miniprot alignments (backstop CDS + cross-check) ---
    alignments: list[MiniprotAlignment] = (
        MiniprotParser(config.miniprot_gff).parse()
        if config.miniprot_gff and Path(config.miniprot_gff).exists()
        else []
    )

    _log.info(
        "loaded: %d helixer loci, %d StringTie tx / %d samples, %d junctions, "
        "%d miniprot alignments",
        len(loci),
        len(all_st),
        len(per_sample),
        len(junctions),
        len(alignments),
    )
    return {
        "genome": genome,
        "h5": h5,
        "helixer_loci": loci,
        "per_sample": per_sample,
        "stringtie_all": all_st,
        "stringtie_agg": st_agg,
        "struct_tpm": struct_tpm,
        "global_max_tpm": global_max_tpm,
        "tpm_overlap_index": tpm_overlap_index,
        "junctions": junctions,
        "junction_index": junction_index,
        "miniprot_alignments": alignments,
    }


# ---------------------------------------------------------------------------
# Stage B: Mikado
# ---------------------------------------------------------------------------


def _run_mikado_stage(
    config: PipelineConfig,
    inputs: dict[str, Any],
    classifications: list[LocusClassification],
) -> list[MikadoLocus]:
    """PREP + MIKADO: emit inputs, run the chain, parse loci → ``[MikadoLocus]``.

    Returns ``[]`` in basic mode (no Mikado) so reconciliation carries every
    Helixer locus through as a backstop gene.
    """
    if not _should_run_mikado(config):
        if config.stringtie_list and not config.protein_db:
            _log.warning(
                "basic mode: StringTie evidence present but --protein-db not "
                "supplied: Mikado requires both; skipping Mikado. "
                "All Helixer loci become backstop genes (no isoform discovery). "
                "Supply --protein-db to enable the full reconciliation chain."
            )
        else:
            _log.info(
                "basic mode: no expression+protein evidence → skipping Mikado; "
                "all Helixer loci become backstop genes"
            )
        return []

    assert config.work_dir is not None  # set by __post_init__
    work = Path(config.work_dir).resolve()
    mik_in = work / "mikado_inputs"
    mik_run = work / "mikado_run"
    mik_in.mkdir(parents=True, exist_ok=True)
    mik_run.mkdir(parents=True, exist_ok=True)

    # --- PREP: GTFs, junctions, list, scoring, configuration ---
    # On resume, if the PREP outputs already exist, reuse them instead of
    # regenerating.  The config_path and scoring_path are the critical outputs
    # that the MIKADO chain consumes.
    config_path_candidate = mik_in / "configuration.yaml"
    prep_reused = False
    if (
        config.resume
        and config_path_candidate.exists()
        and config_path_candidate.stat().st_size > 0
        and (mik_in / "list.txt").exists()
    ):
        config_path: str | Path = config_path_candidate
        scoring_candidates = sorted(mik_in.glob("*.yaml"))
        scoring_candidates = [
            p for p in scoring_candidates if p.name != "configuration.yaml"
        ]
        if scoring_candidates:
            scoring_path: str | Path = scoring_candidates[0]
            junctions_bed: str | Path = mik_in / "junctions.bed"
            if not Path(junctions_bed).exists():
                junctions_bed = _emit_junctions_bed(config, inputs, mik_in)
            _log.info(
                "[PREP] resume: reusing existing PREP outputs in %s", mik_in
            )
            prep_reused = True

    if not prep_reused:
        with _Timer("PREP"):
            helixer_gtf = helixer_gff3_to_gtf(
                config.helixer_gff3, mik_in / "helixer.gtf"
            )
            list_entries = [
                {
                    "file": str(helixer_gtf),
                    "label": "helixer",
                    "is_reference": config.helixer_is_reference,
                    "score": 0,
                }
            ]
            for sid, tx in inputs["per_sample"].items():
                if not tx:
                    continue
                out_gtf = mik_in / f"stringtie_{sid}.gtf"
                stringtie_to_labelled_gtf(tx, str(out_gtf), sid)
                list_entries.append(
                    {"file": str(out_gtf), "label": sid, "score": 0}
                )
            list_path = write_input_list(list_entries, mik_in / "list.txt")

            junctions_bed = _emit_junctions_bed(config, inputs, mik_in)
            scoring_path = install_scoring_profile(config.scoring_profile, mik_in)

            as_knobs: dict[str, Any] = {
                "report": config.as_report,
                "only_confirmed_introns": config.only_confirmed_introns,
                "min_cds_overlap": config.min_cds_overlap,
                "min_cdna_overlap": config.min_cdna_overlap,
                "max_isoforms": config.max_isoforms,
                "keep_retained_introns": config.keep_retained_introns,
                "pad": config.pad,
            }
            config_path = write_configuration(
                genome_fa=config.genome_fasta,
                list_path=str(list_path),
                scoring_path=str(scoring_path),
                junctions_path=str(junctions_bed),
                out_path=str(mik_in / "configuration.yaml"),
                chimera_split=config.chimera_split,
                flank=config.flank,
                use_subprocess=config.use_mikado_configure,
                mikado_bin=config.mikado_bin,
                **as_knobs,
            )
            if config.use_mikado_configure:
                _patch_as_knobs(config_path, as_knobs)

    # --- MIKADO: prepare → external → transdecoder/diamond → serialise → pick ---
    # ``reuse_mikado_dir`` (Phase 17 ablations) lets a pick-only variant reuse the
    # prepared transcripts + serialise DB of a prior run, re-running only pick.
    reuse_dir = Path(config.reuse_mikado_dir) if config.reuse_mikado_dir else None
    with _Timer("MIKADO"):
        prepared_gtf = _reuse_link(reuse_dir, "mikado_prepared.gtf", mik_run)
        prepared_fasta = _reuse_link(reuse_dir, "mikado_prepared.fasta", mik_run)
        if prepared_gtf is None or prepared_fasta is None:
            prepared_gtf, prepared_fasta = run_prepare(
                str(config_path),
                mik_run,
                procs=config.procs,
                mikado_bin=config.mikado_bin,
            )

        db = _reuse_link(reuse_dir, "mikado.db", mik_run)
        if db is None:
            if inputs["h5"] is None:
                _log.warning(
                    "no HDF5: helixer_support/helixer_locus_conf default to 0.0"
                )
            rows = _build_external_scores(
                str(prepared_gtf),
                inputs["h5"],
                inputs["struct_tpm"],
                inputs["global_max_tpm"],
                helixer_weight=config.helixer_support_weight,
            )
            external_tsv = write_external_scores_tsv(
                rows, mik_in / "external_scores.tsv"
            )

            assert config.protein_db is not None  # guaranteed by _should_run_mikado
            orfs_bed = run_transdecoder(
                prepared_fasta,
                mik_run,
                transdecoder_bin_dir=config.transdecoder_bin_dir,
            )
            blast_out = run_diamond(
                prepared_fasta,
                config.protein_db,
                mik_run,
                threads=config.threads,
                diamond_bin=config.diamond_bin,
            )
            run_serialise(
                str(config_path),
                prepared_fasta,
                orfs_bed,
                blast_out,
                config.protein_db,
                str(junctions_bed),
                str(external_tsv),
                config.genome_fasta,
                mik_run,
                mikado_bin=config.mikado_bin,
            )
        else:
            _log.info(
                "reuse: skipping prepare/external/transdecoder/diamond/serialise, "
                "reusing Mikado DB + prepared transcripts from %s (re-running pick only)",
                reuse_dir,
            )
        loci_gff3 = run_pick(
            str(config_path),
            str(scoring_path),
            prepared_gtf,
            mik_run,
            procs=config.procs,
            mikado_bin=config.mikado_bin,
        )

    mikado_loci: list[MikadoLocus] = parse_loci_gff3(
        str(loci_gff3),
        metrics_tsv=str(mik_run / "mikado.loci.metrics.tsv"),
        scores_tsv=str(mik_run / "mikado.loci.scores.tsv"),
    )
    _log.info("parsed %d Mikado loci", len(mikado_loci))
    return mikado_loci


def _mikado_output_paths(config: PipelineConfig) -> list[Path]:
    """The three files ``parse_loci_gff3`` needs: the MIKADO stage's checkpoint."""
    assert config.work_dir is not None  # set by __post_init__
    mik_run = Path(config.work_dir) / "mikado_run"
    return [
        mik_run / "mikado.loci.gff3",
        mik_run / "mikado.loci.metrics.tsv",
        mik_run / "mikado.loci.scores.tsv",
    ]


def _reparse_mikado_loci(config: PipelineConfig) -> list[MikadoLocus]:
    """Re-parse the existing Mikado loci on resume (skips the external chain)."""
    assert config.work_dir is not None  # set by __post_init__
    mik_run = Path(config.work_dir) / "mikado_run"
    loci: list[MikadoLocus] = parse_loci_gff3(
        str(mik_run / "mikado.loci.gff3"),
        metrics_tsv=str(mik_run / "mikado.loci.metrics.tsv"),
        scores_tsv=str(mik_run / "mikado.loci.scores.tsv"),
    )
    _log.info("[MIKADO] resume: re-parsed %d loci from checkpoint", len(loci))
    return loci


def _parse_external_mikado_loci(config: PipelineConfig) -> list["MikadoLocus"]:
    """Parse externally-supplied Mikado loci (``--mikado-loci``).

    Auto-detects companion ``*.metrics.tsv`` and ``*.scores.tsv`` next to the GFF3
    using Mikado's naming convention (``<stem>.metrics.tsv``). Warns if absent,
    they carry combined_score and blast_score but are not strictly required.
    """
    assert config.mikado_loci_gff3 is not None
    gff3 = Path(config.mikado_loci_gff3)
    if not gff3.exists():
        raise FileNotFoundError(f"--mikado-loci GFF3 not found: {gff3}")
    if gff3.stat().st_size == 0:
        raise ValueError(f"--mikado-loci GFF3 is empty: {gff3}")

    stem = gff3.with_suffix("")  # strip .gff3
    metrics = Path(f"{stem}.metrics.tsv")
    scores = Path(f"{stem}.scores.tsv")

    metrics_str: str | None = str(metrics) if metrics.exists() else None
    scores_str: str | None = str(scores) if scores.exists() else None

    if metrics_str is None:
        _log.warning(
            "--mikado-loci: metrics TSV not found (%s), blast_score and "
            "partial-ORF detection will fall back to heuristics",
            metrics,
        )
    if scores_str is None:
        _log.warning(
            "--mikado-loci: scores TSV not found (%s), combined_score will "
            "be None for all transcripts",
            scores,
        )

    loci: list[MikadoLocus] = parse_loci_gff3(
        str(gff3), metrics_tsv=metrics_str, scores_tsv=scores_str
    )
    _log.info(
        "--mikado-loci: parsed %d external loci from %s "
        "(metrics=%s, scores=%s)",
        len(loci),
        gff3,
        "yes" if metrics_str else "no",
        "yes" if scores_str else "no",
    )
    return loci


def _reuse_link(reuse_dir: Path | None, name: str, dest_dir: Path) -> Path | None:
    """Symlink ``reuse_dir/name`` into ``dest_dir`` and return it; ``None`` if absent.

    Used by the Phase-17 ablation reuse path so ``mikado pick`` finds the reused
    prepared transcripts + serialise DB in its own ``--output-dir`` (mirrors the
    M3 ``--resume`` symlink approach). Falls back to a copy if symlinks are
    unsupported.
    """
    if reuse_dir is None:
        return None
    src = Path(reuse_dir) / name
    if not src.exists():
        return None
    dest = Path(dest_dir) / name
    if not dest.exists():
        try:
            dest.symlink_to(src.resolve())
        except OSError:
            import shutil

            shutil.copy2(src, dest)
    return dest


def _emit_junctions_bed(
    config: PipelineConfig, inputs: dict[str, Any], mik_in: Path
) -> Path:
    """Produce the Portcullis-format junctions BED for Mikado.

    Prefers Portcullis (if a binary + BAMs are given), else writes the STAR-SJ
    union (possibly empty when no SJ inputs exist).
    """
    if config.portcullis_bin and config.bam_paths:
        bam_paths_union: list[str | Path] = list(config.bam_paths)
        return run_portcullis(
            config.genome_fasta,
            bam_paths_union,
            mik_in,
            threads=config.threads,
            portcullis_bin=config.portcullis_bin,
        )
    return junctions_to_portcullis_tab(inputs["junctions"], mik_in / "junctions.bed")


# ---------------------------------------------------------------------------
# Stage C: per-gene finalisation (Phase 7: backstop CDS, junction fix, gate)
# ---------------------------------------------------------------------------


def _finalize_one_gene(
    gene: "ReconciledGene",
    config: PipelineConfig,
    genome: GenomeAccessor,
    alignments: list[MiniprotAlignment],
    junctions: JunctionIndex | list[SpliceJunction],
    stats: RunStats | None = None,
) -> "ReconciledGene":
    """Finalize a single gene (backstop CDS + junction fix, or Mikado cross-check).

    Backstop genes are junction-corrected **before** CDS projection so the CDS is
    built on the verified structure (a deviation from the phase-spec listing order
    that strictly dominates it: correcting after a CDS is set risks reverting a
    valid junction fix on CDS-containment, leaving the structure unfixed). The
    per-gene TransDecoder fallback is **not** run here, CDS-less backstop genes
    are rescued in one batched invocation afterward, so this function
    is a pure, picklable, per-gene transform safe to dispatch to a worker.
    """
    extra: list[QCFlag] = []
    if gene.origin == "helixer_backstop":
        gene = refine_backstop_gene(
            gene,
            junctions,
            tolerance=config.junction_tolerance,
            min_reads=config.junction_min_reads,
            stats=stats,
        )
        # transdecoder_bin_dir=None: miniprot-only here; the TransDecoder fallback
        # for the still-CDS-less backstop genes is batched (one subprocess) below.
        gene = assign_backstop_cds(
            gene,
            alignments,
            genome=genome,
            transdecoder_bin_dir=None,
            stats=stats,
        )
    else:
        flag = cds_cross_check(gene, alignments)
        if flag is not None:
            extra.append(flag)
    if extra:
        gene = attrs.evolve(gene, flags=dedup_flags([*gene.flags, *extra]))
    return gene


def _finalize_scaffold_worker(
    payload: tuple[
        list[ReconciledGene],
        PipelineConfig,
        list[MiniprotAlignment],
        JunctionIndex | list[SpliceJunction],
    ],
) -> tuple[list[ReconciledGene], RunStats]:
    """Process one scaffold's genes in a worker; re-opens the genome per process.

    ``payload`` is ``(genes, config, alignments, junctions)``. The genome FASTA
    handle is **not** picklable, so each worker opens its own ``GenomeAccessor``
    from ``config.genome_fasta`` (Phase 21's ``JunctionIndex`` / miniprot
    alignments *are* picklable and are passed in). Returns ``(refined_genes,
    worker_stats)`` so the parent can sum the decision counters back (a process
    can't mutate the parent's ``RunStats``).
    """
    genes, config, alignments, junctions = payload
    genome = GenomeAccessor(config.genome_fasta)
    stats = RunStats()
    try:
        refined = [
            _finalize_one_gene(g, config, genome, alignments, junctions, stats)
            for g in genes
        ]
    finally:
        genome.close()
    return refined, stats


def _finalize_genes_parallel(
    genes: list[ReconciledGene],
    config: PipelineConfig,
    inputs: dict[str, Any],
    workers: int,
    stats: RunStats | None = None,
    _executor_cls: Any = None,
) -> list[ReconciledGene]:
    """Run :func:`_finalize_one_gene` over a bounded pool, grouped by scaffold.

    Genes are independent ``attrs`` objects; partitioning by scaffold keeps each
    worker's locality and lets it open one genome handle. The result is rebuilt in
    the **original gene order** (worker scheduling never changes counts/IDs/flags).
    ``_executor_cls`` is a test seam (the
    suite injects a ``ThreadPoolExecutor`` so the per-worker genome re-open is
    observed in-process); production uses ``ProcessPoolExecutor``.
    """
    from concurrent.futures import ProcessPoolExecutor

    executor_cls = _executor_cls or ProcessPoolExecutor
    alignments: list[MiniprotAlignment] = inputs["miniprot_alignments"]
    junctions: JunctionIndex | list[SpliceJunction] = (
        inputs.get("junction_index") or inputs["junctions"]
    )

    by_scaffold: dict[str, list[ReconciledGene]] = {}
    for g in genes:
        by_scaffold.setdefault(g.seqid, []).append(g)
    payloads = [(grp, config, alignments, junctions) for grp in by_scaffold.values()]

    refined_by_id: dict[str, ReconciledGene] = {}
    merged = RunStats()
    with executor_cls(max_workers=workers) as executor:
        for refined, worker_stats in executor.map(_finalize_scaffold_worker, payloads):
            for g in refined:
                refined_by_id[g.gene_id] = g
            merged.merge(worker_stats)
    if stats is not None:
        stats.merge(merged)
    # Reassemble in the order reconcile produced (sorted by seqid/start).
    return [refined_by_id[g.gene_id] for g in genes]


def _finalize_genes(
    genes: list[ReconciledGene],
    config: PipelineConfig,
    inputs: dict[str, Any],
    stats: RunStats | None = None,
    _executor_cls: Any = None,
) -> list[ReconciledGene]:
    """Backstop CDS + junction correction, Mikado-origin cross-check, codon gate.

    The per-gene refine loop runs serially (``finalize_workers <= 1``) or across a
    bounded process pool by scaffold; either path
    yields the identical gene set. CDS-less backstop genes are then rescued in one
    batched TransDecoder invocation (replacing the per-gene subprocess), and the
    read-only structural codon gate merges its flags into every gene.
    """
    genome: GenomeAccessor = inputs["genome"]
    alignments: list[MiniprotAlignment] = inputs["miniprot_alignments"]
    # Prefer the pre-built junction index; fall back to the raw list
    # if an older caller's inputs dict predates it (both give identical results).
    junctions: JunctionIndex | list[SpliceJunction] = (
        inputs["junction_index"] if "junction_index" in inputs else inputs["junctions"]
    )

    workers = getattr(config, "finalize_workers", 1) or 1
    if workers > 1 and len(genes) > 1:
        refined = _finalize_genes_parallel(
            genes, config, inputs, workers, stats=stats, _executor_cls=_executor_cls
        )
    else:
        refined = [
            _finalize_one_gene(g, config, genome, alignments, junctions, stats)
            for g in genes
        ]

    # Batched backstop TransDecoder: one multi-FASTA invocation over
    # every still-CDS-less backstop transcript, gated by ``backstop_transdecoder``.
    if config.backstop_transdecoder and config.transdecoder_bin_dir:
        refined = batch_backstop_transdecoder(
            refined, genome, config.transdecoder_bin_dir, stats=stats
        )

    # Structural codon gate (read-only) → merge flags into each gene.
    flag_map: dict[str, list[QCFlag]] = validate_all(
        refined,
        genome=genome,
        short_cds_threshold=config.short_cds_threshold,
        short_exon_threshold=config.short_exon_threshold,
        long_intron_threshold=config.long_intron_threshold,
        transl_table=config.transl_table,
        transl_table_map=config.transl_table_map,
    )
    # The Helixer CDS-channel confidence signal used to
    # discriminate lncRNA (channel quiet) from a fragmentary/failed coding locus.
    # Computed only for ORF-less genes (irrelevant to coding/pseudogene biotypes)
    # and only when an HDF5 reader is present; a fetch failure degrades to None.
    h5: HDF5ConfidenceReader | None = inputs.get("h5")

    def _cds_channel_conf(gene: ReconciledGene) -> float | None:
        if h5 is None:
            return None
        p = _primary(gene)
        if p.cds:
            return None
        try:
            return h5.get_cds_channel_confidence(gene.seqid, p.exons)
        except Exception:  # noqa: BLE001 - never let a confidence lookup sink the run
            return None

    out: list[ReconciledGene] = []
    for g in refined:
        gate_flags = flag_map.get(g.gene_id, [])
        merged = attrs.evolve(g, flags=dedup_flags([*g.flags, *gate_flags]))
        # A homology-backed disabled ORF (premature stop / mod-3
        # frameshift) becomes a pseudogene candidate, biotype 'pseudogene',
        # PSEUDOGENE_CANDIDATE flag, Tier-1 demotion. Uses the gate flags so the
        # INTERNAL_STOP signal is shared, never recomputed.
        typed = apply_pseudogene_typing(
            merged, gate_flags, min_homology=config.pseudogene_min_homology
        )
        # Assign protein_coding / lncRNA / ncRNA_undetermined from
        # the existing signals (runs AFTER pseudogene typing, which it never
        # overrides: pseudogene wins for a disabled homolog). Non-coding biotypes
        # get a coherent tier (never Tier 1/2).
        out.append(
            assign_biotype(
                typed,
                cds_channel_conf=_cds_channel_conf(typed),
                min_lncrna_length=config.lncrna_min_length,
                max_noncoding_cds_conf=config.noncoding_max_cds_conf,
            )
        )
    return out


def _assign_stringtie_tpm(
    genes: list[ReconciledGene], inputs: dict[str, Any]
) -> list[ReconciledGene]:
    """Attach per-transcript StringTie TPM by exonic overlap (count-neutral).

    Mikado consumes the StringTie TPM only as an external scoring metric and does
    not emit it back, so without this step every reconciled transcript carries
    ``tpm=None`` and the report's ``frac_genes_tpm_pass`` is always 0, even with
    valid StringTie input. Here each transcript is assigned the TPM of the
    same-strand StringTie structure with the greatest exonic overlap (the exact
    rule the ``evidence`` scorer uses; shared via
    :func:`io.stringtie.best_overlapping_tpm`). Sets ``tpm`` only, never changes
    structure, CDS, tier, biotype, or the gene set.
    """
    index = inputs.get("tpm_overlap_index") or {}
    if not index:
        return genes
    out: list[ReconciledGene] = []
    for gene in genes:
        new_tx = []
        changed = False
        for t in gene.transcripts:
            tpm = best_overlapping_tpm(index, t.seqid, t.strand, t.exons)
            if tpm is not None and tpm != t.tpm:
                new_tx.append(attrs.evolve(t, tpm=tpm))
                changed = True
            else:
                new_tx.append(t)
        out.append(attrs.evolve(gene, transcripts=new_tx) if changed else gene)
    return out


def _maybe_add_structured_ncrna(
    genes: list[ReconciledGene], config: PipelineConfig
) -> list[ReconciledGene]:
    """Append structured-ncRNA loci from the opt-in hook (no-op when disabled).

    Off by default (``config.ncrna_scan`` False) → returns ``genes`` unchanged and
    runs no subprocess, so the default/golden path is untouched. When enabled, runs
    tRNAscan-SE / Infernal and merges the resulting tRNA/rRNA/snoRNA loci, keeping
    the gene list sorted by ``(seqid, start)``.
    """
    if not config.ncrna_scan:
        return genes
    from helixforge.prep._subprocess import ToolError
    from helixforge.prep.ncrna import scan_structured_ncrna

    assert config.work_dir is not None
    try:
        ncrna = scan_structured_ncrna(
            config.genome_fasta,
            Path(config.work_dir) / "ncrna",
            enabled=True,
            tool=config.ncrna_tool,
            rfam_cm=config.ncrna_rfam_cm,
            threads=config.threads,
        )
    except (ToolError, FileNotFoundError, OSError, ValueError) as exc:
        _log.warning(
            "structured-ncRNA hook skipped: tool missing or failed (%s): %s",
            type(exc).__name__,
            exc,
        )
        return genes
    if not ncrna:
        return genes
    _log.info("structured-ncRNA hook: merging %d loci into the gene set", len(ncrna))
    return sorted([*genes, *ncrna], key=lambda g: (g.seqid, g.start, g.end))


def _maybe_flag_variants(
    genes: list[ReconciledGene], config: PipelineConfig
) -> list[ReconciledGene]:
    """Flag CDS-overlapping high-impact VCF variants (experimental; no-op by default).

    Off unless both ``config.vcf_flag_impact`` and ``config.vcf_path`` are set
    Returns ``genes`` unchanged when disabled, so the default/golden path is
    untouched. When enabled, loads the **decomposed** VCF and attaches
    ``VARIANT_IMPACTED`` to genes whose CDS overlaps a high-impact variant, flag
    only, never re-typing or re-tiering.
    """
    if not config.vcf_flag_impact or not config.vcf_path:
        return genes
    from helixforge.reconcile.vcf import flag_variant_impacted_genes, load_vcf

    variants = load_vcf(config.vcf_path)
    flagged = flag_variant_impacted_genes(genes, variants, enabled=True)
    n = sum(1 for g in flagged if any(f.name == "VARIANT_IMPACTED" for f in g.flags))
    _log.info(
        "VCF impact (experimental): %d gene(s) flagged VARIANT_IMPACTED "
        "from %d variant(s)",
        n,
        len(variants),
    )
    return flagged


def _maybe_gate_te(
    genes: list[ReconciledGene], config: PipelineConfig
) -> list[ReconciledGene]:
    """Flag + gate transposable-element-encoded ORFs from an optional EDTA GFF3.

    Off unless ``config.te_annotation`` is set → returns ``genes`` unchanged, so
    the default/golden path is untouched (EDTA is the only TE signal). When set,
    flags every TE-overlapping model (``TE_OVERLAP``) and reclassifies good-ORF
    genes whose model-fraction TE overlap is at/above ``te_overlap_threshold`` as
    ``transposable_element``, using the EDTA ``Classification`` order, so
    knob/satellite/low-complexity never gate.
    """
    if not config.te_annotation:
        return genes
    if not Path(config.te_annotation).exists():
        _log.warning("TE annotation not found: %s, skipping TE gating", config.te_annotation)
        return genes
    from helixforge.reconcile.te import gate_te, parse_edta_te_intervals

    te_index = parse_edta_te_intervals(config.te_annotation, te_classes=config.te_classes)
    genes, n_flagged, n_reclassified = gate_te(
        genes, te_index, threshold=config.te_overlap_threshold
    )
    _log.info(
        "EDTA TE gating: %d gene(s) overlap a TE feature (TE_OVERLAP); "
        "%d good-ORF gene(s) reclassified transposable_element "
        "(overlap >= %.2f)",
        n_flagged,
        n_reclassified,
        config.te_overlap_threshold,
    )
    return genes


# ---------------------------------------------------------------------------
# Phase 31: functional annotation (D1/D2) + genome-level report (D3/D4)
# ---------------------------------------------------------------------------


def _maybe_annotate_function(
    genes: list[ReconciledGene],
    config: PipelineConfig,
    genome: GenomeAccessor,
) -> tuple[list[ReconciledGene], dict[str, Any]]:
    """Run the opt-in functional-annotation hook.

    Off by default (``config.functional_annotation`` False) → returns ``(genes,
    {})`` and runs no subprocess, so the default/golden path is untouched. When
    enabled, writes a protein FASTA for the ``protein_coding``/``pseudogene`` set,
    runs InterProScan/eggNOG, and returns ``(genes_with_DOMAIN_COMPLETE_flag,
    functional_records)``. The records feed the GFF3 ``Ontology_term``/``Dbxref``
    emission; the credibility flag is added by ``apply_domain_credibility`` (D2,
    count-neutral, an INFO flag only).
    """
    if not config.functional_annotation:
        return genes, {}
    from helixforge.export.writers import write_protein_fasta
    from helixforge.prep._subprocess import ToolError
    from helixforge.prep.function import (
        ANNOTATABLE_BIOTYPES,
        annotate_function,
        apply_domain_credibility,
    )

    assert config.work_dir is not None
    out_dir = Path(config.work_dir) / "function"
    out_dir.mkdir(parents=True, exist_ok=True)
    coding = [
        g for g in genes if (g.biotype or "protein_coding") in ANNOTATABLE_BIOTYPES
    ]
    proteins = out_dir / "proteins.fa"
    write_protein_fasta(coding, genome, str(proteins))
    try:
        records = annotate_function(
            coding,
            proteins,
            out_dir=out_dir,
            tool=config.functional_tool,
            db=config.functional_db,
            enabled=True,
            threads=config.threads,
            interproscan_bin=config.interproscan_bin,
            eggnog_bin=config.eggnog_bin,
        )
    except (ToolError, FileNotFoundError, OSError, ValueError) as exc:
        _log.warning("functional annotation skipped (%s): %s", type(exc).__name__, exc)
        return genes, {}
    genes = apply_domain_credibility(genes, records)
    return genes, records


def _maybe_trace_reorder(
    genes: list[ReconciledGene],
    config: PipelineConfig,
    inputs: dict[str, Any],
    functional: dict[str, Any],
) -> tuple[list[ReconciledGene], dict[str, Any]]:
    """Re-elect each multi-transcript gene's canonical/primary isoform via TRaCE.

    Off by default (``config.trace_primary`` False) → returns ``(genes,
    functional)`` unchanged, so the golden ids/primary are preserved (Phase 33b).
    When on, runs :func:`reconcile.trace.trace_order` per multi-transcript gene
    and renumbers via :func:`mikado_integrate._renumber` with the elected order.

    **Reorder-only and count-neutral**: it permutes + renumbers a gene's
    transcripts (changing ``.N`` ids, ``is_primary``, ``primary_transcript_id``,
    and ``trace_rank``) but keeps the gene's ``as_events`` and ``flags`` exactly as
    derived during reconciliation, so the gene / tier / origin / AS-event set are
    unchanged (the golden gate for the TRaCE-on variant). Functional records,
    keyed by transcript id, are re-keyed onto the new ids so the GFF3
    ``Ontology_term``/``Dbxref`` still attach to the right isoform.
    """
    if not config.trace_primary:
        return genes, functional

    from helixforge.reconcile.mikado_integrate import _renumber
    from helixforge.reconcile.trace import (
        TraceParams,
        domain_coverage_by_transcript,
        trace_order,
    )

    params = TraceParams(
        max_aed=config.trace_max_aed,
        min_tpm=config.trace_min_tpm,
        min_overlap=config.trace_min_overlap,
        weight_domain=config.trace_weight_domain,
        weight_protein=config.trace_weight_protein,
        weight_cdna=config.trace_weight_cdna,
        use_domain=config.trace_use_domain,
    )
    per_sample: dict[str, list[Any]] = inputs.get("per_sample", {})
    domain_cov = (
        domain_coverage_by_transcript(genes, functional)
        if (config.trace_use_domain and functional)
        else None
    )

    remap: dict[str, str] = {}
    n_reordered = 0
    out: list[ReconciledGene] = []
    for g in genes:
        if len(g.transcripts) < 2:
            out.append(g)
            continue
        ordered = trace_order(
            g.transcripts, per_sample, domain_coverage=domain_cov, params=params
        )
        order_ids = [t.transcript_id for t in ordered]
        new_tx = _renumber(g.transcripts, g.gene_id, order=order_ids)
        new_tx = [attrs.evolve(t, trace_rank=i + 1) for i, t in enumerate(new_tx)]
        # old structural id (ordered[i]) → new positional id (new_tx[i]).
        for old, new in zip(order_ids, (t.transcript_id for t in new_tx)):
            remap[old] = new
        if order_ids != [t.transcript_id for t in g.transcripts]:
            n_reordered += 1
        out.append(
            attrs.evolve(
                g,
                transcripts=new_tx,
                primary_transcript_id=new_tx[0].transcript_id,
            )
        )

    if functional and remap:
        functional = {remap.get(k, k): v for k, v in functional.items()}
    _log.info(
        "TRaCE: elected canonical transcript for %d multi-isoform gene(s) "
        "(%d reordered from the combined_score primary)",
        sum(1 for g in genes if len(g.transcripts) > 1),
        n_reordered,
    )
    return out, functional


def _maybe_write_report(
    genes: list[ReconciledGene],
    config: PipelineConfig,
    inputs: dict[str, Any],
    stats: RunStats,
) -> None:
    """Write the default end-of-run ``report.json`` + ``report.html``.

    On by default (``config.write_report``). **Best-effort + count-neutral**: it
    only reads the finished gene set + run telemetry, and any failure is logged
    and swallowed so it never sinks a completed run. The bench
    completeness tools are *not* run here (off the default path); the structural
    + support sections come entirely from in-memory data.
    """
    if not config.write_report:
        return
    try:
        from helixforge.io.bam import overall_mapping_rate
        from helixforge.stats.report import (
            build_report,
            write_report_html,
            write_report_json,
        )

        # CRAM inputs decode against the (always-present) genome FASTA, so the
        # mapping summary never needs a remote reference fetch; a
        # BAM ignores reference_filename, so this is count-neutral on the golden
        # path.
        bam_stats = (
            overall_mapping_rate(
                config.bam_paths, reference_filename=config.genome_fasta
            )
            if config.bam_paths
            else None
        )
        # The default report is built from in-memory data only (no per-gene HDF5
        # re-read): mean Helixer support is left to the richer `helixforge stats`
        # path so the end-of-run report stays cheap on a whole-genome gene set.
        report = build_report(
            genes,
            junctions=inputs.get("junctions"),
            min_reads=config.junction_min_reads,
            tpm_threshold=config.min_tpm,
            bam_stats=bam_stats,
            run_stats=stats.as_dict(),
        )
        write_report_json(report, f"{config.output_prefix}.report.json")
        write_report_html(report, f"{config.output_prefix}.report.html")
        _log.info(
            "wrote genome-level report: %s.report.json + .report.html",
            config.output_prefix,
        )
    except Exception as exc:  # noqa: BLE001 - reporting never sinks a finished run
        _log.warning("could not write genome-level report: %s", exc)


# ---------------------------------------------------------------------------
# Stage D: output
# ---------------------------------------------------------------------------


def _gff3_metadata(
    config: PipelineConfig,
) -> tuple[dict[str, int] | None, Any | None, str | None]:
    """Resolve ``(sequence_regions, provenance, embed_fasta)`` for the GFF3 header.

    Best-effort and count-neutral: scaffold lengths come from
    the genome ``.fai`` and the provenance preamble from the resolved config +
    configured tool binaries + input MD5s. Any failure degrades to ``None`` so a
    finished run is never sunk by a header-metadata problem.
    """
    seq_regions: dict[str, int] | None = None
    if config.gff3_sequence_regions:
        try:
            with GenomeAccessor(config.genome_fasta) as genome:
                seq_regions = genome.get_scaffold_lengths()
        except (OSError, ValueError) as exc:  # pragma: no cover - defensive
            _log.warning("could not read scaffold lengths for GFF3 header: %s", exc)

    provenance: Any | None = None
    if config.gff3_provenance:
        try:
            from helixforge.provenance import build_provenance

            td = config.transdecoder_bin_dir
            tool_bins = {
                "mikado": config.mikado_bin,
                "diamond": config.diamond_bin,
                "transdecoder": (
                    str(Path(td) / "TransDecoder.LongOrfs")
                    if td
                    else "TransDecoder.LongOrfs"
                ),
            }
            input_files = {
                "genome_fasta": config.genome_fasta,
                "helixer_gff3": config.helixer_gff3,
                "helixer_h5": config.helixer_h5,
                "miniprot_gff": config.miniprot_gff,
                "protein_db": config.protein_db,
            }
            provenance = build_provenance(
                params=config.to_nested_dict(),
                tool_bins=tool_bins,
                input_files=input_files,
            )
        except Exception as exc:  # pragma: no cover - never sink a run on provenance
            _log.warning("could not build GFF3 provenance preamble: %s", exc)

    embed = config.genome_fasta if config.gff3_embed_fasta else None
    return seq_regions, provenance, embed


def _write_outputs(
    genes: list[ReconciledGene],
    config: PipelineConfig,
    functional: dict[str, Any] | None = None,
) -> None:
    """Write the full GFF3 + cumulative tier1/tier2/tier3 GFF3s.

    ``functional`` optionally threads the InterProScan/eggNOG
    records into the GFF3 so each mRNA carries ``Ontology_term``/``Dbxref``;
    ``None`` (the default path) keeps the per-feature output byte-identical
    (count-neutral). The GFF3 header additionally carries ``##sequence-region``
    directives + a ``#!`` provenance preamble unless disabled
    on the config; these are comments/directives and never change parsed counts.
    """
    seq_regions, provenance, embed = _gff3_metadata(config)

    def _writer(path: str) -> GFF3Writer:
        return GFF3Writer(
            path,
            sequence_regions=seq_regions,
            provenance=provenance,
            embed_fasta=embed,
        )

    _writer(f"{config.output_prefix}.gff3").write_genes(genes, functional=functional)
    for n in (1, 2, 3):
        subset = [g for g in genes if g.tier <= n]
        _writer(f"{config.output_prefix}.tier{n}.gff3").write_genes(
            subset, functional=functional
        )
    _log.info("wrote GFF3: %s.gff3 (+ tier1/2/3)", config.output_prefix)


_REPORT_COLUMNS = (
    "gene_id",
    "seqid",
    "start",
    "end",
    "strand",
    "tier",
    "origin",
    "biotype",
    "primary_transcript_id",
    "num_isoforms",
    "num_as_events",
    "has_cds",
    "protein_id",
    "max_tpm",
    "junction_support",
    "helixer_support",
    "combined_score",
    "aed",
    "flags",
)


def _fmt(value: object) -> str:
    return "" if value is None else str(value)


def _write_report(genes: list[ReconciledGene], report_path: str) -> None:
    """One TSV row per gene with the Phase 8 report columns."""
    with atomic_write(report_path) as fh:
        fh.write("\t".join(_REPORT_COLUMNS) + "\n")
        for g in genes:
            p = _primary(g)
            row = [
                g.gene_id,
                g.seqid,
                g.start,
                g.end,
                g.strand,
                g.tier,
                g.origin,
                _fmt(g.biotype),
                g.primary_transcript_id,
                len(g.transcripts),
                len(g.as_events),
                "true" if p.cds else "false",
                _fmt(p.protein_id),
                _fmt(p.tpm),
                _fmt(p.junction_support_fraction),
                _fmt(p.confidence),  # helixer_support (None unless enriched)
                _fmt(p.combined_score),
                "",  # aed, computed in Phase 9 stats
                ",".join(f.name for f in g.flags),
            ]
            fh.write("\t".join(_fmt(v) for v in row) + "\n")
    _log.info("wrote report: %s (%d genes)", report_path, len(genes))


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def run_pipeline(config: PipelineConfig) -> list[ReconciledGene]:
    """Run the full A→D pipeline; return the list of ``ReconciledGene``.

    Opens the genome FASTA and (optional) Helixer HDF5 once and closes them in a
    ``finally``. The Helixer id → HFG map is loaded from
    ``config.id_map_path`` if present and re-saved, giving stable gene IDs across
    reruns.
    """
    assert config.work_dir is not None  # set by __post_init__
    validate_scoring_profile(config.scoring_profile)

    # Fail-fast preflight: verify work-dir writability and tool-chain readiness
    # before the expensive input-loading and Mikado stages.  A missing tool or
    # an unwritable dir is caught in seconds, not after hours.
    from helixforge.prep.preflight import check_tool_chain, check_work_dir

    workdir_errors = check_work_dir(config)
    if workdir_errors:
        raise RuntimeError(
            "work-dir preflight failed: fix before running:\n  "
            + "\n  ".join(workdir_errors)
        )
    tool_errors = check_tool_chain(config)
    if tool_errors:
        raise FileNotFoundError(
            "tool-chain preflight failed: install missing tools:\n  "
            + "\n  ".join(tool_errors)
        )

    if config.stringtie_list and not config.protein_db:
        import warnings

        warnings.warn(
            "--protein-db not supplied: Mikado will not run. Output will be "
            "backstop-only with no isoform discovery. Supply --protein-db to "
            "enable the full reconciliation chain.",
            UserWarning,
            stacklevel=2,
        )
    inputs = _load_inputs(config)
    genome: GenomeAccessor = inputs["genome"]
    h5: HDF5ConfidenceReader | None = inputs["h5"]
    checkpoint = Checkpoint(
        Path(config.work_dir) / "checkpoint.json", enabled=config.resume
    )
    stats = RunStats()
    try:
        loci: list[HelixerLocus] = inputs["helixer_loci"]

        # Dump the fully-resolved nested run-config for reproducibility (D2,
        # §3.6). Best-effort: never let a config-dump failure sink the run.
        _dump_resolved_config(config)

        with _Timer("CLASSIFY"):
            classifications: list[LocusClassification] = classify_loci(
                loci,
                stringtie_aggregated=inputs["stringtie_agg"] or None,
                stringtie_transcripts=inputs["stringtie_all"] or None,
                bam_paths=config.bam_paths or None,
                bigwig_paths=config.bigwig_paths or None,
                min_tpm=config.min_tpm,
                min_samples=config.min_samples,
                coverage_threshold=config.coverage_threshold,
                near_zero_coverage=config.near_zero_coverage,
            )

        # MIKADO is the expensive stage. Three paths, in priority order:
        #
        # 1. --mikado-loci: externally-supplied Mikado result → skip the
        #    entire PREP + MIKADO chain and parse the supplied files.
        # 2. resume + valid checkpoint: re-parse the on-disk loci instead
        #    of re-running the external chain.
        # 3. Fresh run: PREP → MIKADO → parse.
        #
        # Basic mode (no StringTie + protein_db) produces no loci files and
        # is never checkpointed: path 3 always applies.
        if config.mikado_loci_gff3:
            _log.info(
                "--mikado-loci: skipping PREP + MIKADO, using external loci"
            )
            mikado_loci = _parse_external_mikado_loci(config)
        else:
            mikado_outputs = _mikado_output_paths(config)
            if (
                config.resume
                and _should_run_mikado(config)
                and checkpoint.is_complete("MIKADO", mikado_outputs)
            ):
                _log.info(
                    "[MIKADO] resume: checkpoint valid, skipping external chain"
                )
                mikado_loci = _reparse_mikado_loci(config)
            else:
                mikado_loci = _run_mikado_stage(config, inputs, classifications)
                if _should_run_mikado(config):
                    checkpoint.mark("MIKADO", mikado_outputs)

        with _Timer("RECONCILE"):
            id_map = _load_id_map(config.id_map_path)
            genes, id_map, _admissions = reconcile(
                loci,
                classifications,
                mikado_loci,
                id_map=id_map,
                admit_novel=config.admit_novel,
                novel_evidence_floor=config.novel_evidence_floor,
                reciprocal_overlap=config.reciprocal_overlap,
                min_cds_overlap=config.min_cds_overlap,
                min_cdna_overlap=config.min_cdna_overlap,
                allocator=IdAllocator(
                    base=config.id_base, novel_base=config.novel_base
                ),
                stats=stats,
                # Paralog/tandem-array merge guards. The pre-built
                # junction index (canonical motif + read counts from STAR) gates an
                # accepted merge on a canonical, >= merge_min_gap_reads bridge; the
                # genome + paralog threshold (default None=off) gate high-identity
                # adjacent loci.
                junctions=inputs.get("junction_index") or inputs.get("junctions"),
                merge_min_gap_reads=config.merge_min_gap_reads,
                merge_require_canonical=config.merge_require_canonical,
                merge_tolerance=config.junction_tolerance,
                genome=(
                    genome if config.paralog_identity_threshold is not None else None
                ),
                paralog_identity_threshold=config.paralog_identity_threshold,
                paralog_kmer_k=config.paralog_kmer_k,
            )
            _save_id_map(config.id_map_path, id_map)
            genes = _finalize_genes(genes, config, inputs, stats=stats)
            # Attach per-transcript StringTie TPM by exonic overlap (the metric
            # Mikado never emits back). Count-neutral, sets tpm only.
            genes = _assign_stringtie_tpm(genes, inputs)
            genes = _maybe_add_structured_ncrna(genes, config)
            genes = _maybe_flag_variants(genes, config)
            genes = _maybe_gate_te(genes, config)

        # Opt-in functional annotation (off by default → no-op,
        # genes + functional map unchanged). Runs before OUTPUT so the GFF3 carries
        # Ontology_term/Dbxref and the DOMAIN_COMPLETE credibility flag is present.
        genes, functional = _maybe_annotate_function(genes, config, genome)

        # TRaCE canonical-transcript election (off by default → no-op).
        # Runs after isoforms + any functional domains exist so the domain voter
        # can use per-isoform coverage when available; re-keys the functional
        # records onto the elected ids. Reorder-only, counts/tiers/AS unchanged.
        genes, functional = _maybe_trace_reorder(genes, config, inputs, functional)

        with _Timer("OUTPUT"):
            assert config.report_path is not None  # set by __post_init__
            output_files = [
                Path(f"{config.output_prefix}.gff3"),
                Path(config.report_path),
            ]
            # Functional annotation augments the GFF3; force a rewrite when it ran
            # so a resumed run picks up the new attributes.
            if (
                config.resume
                and not functional
                and checkpoint.is_complete("OUTPUT", output_files)
            ):
                _log.info("[OUTPUT] resume: checkpoint valid, outputs already written")
            else:
                _write_outputs(genes, config, functional=functional or None)
                _write_report(genes, config.report_path)
                checkpoint.mark("OUTPUT", output_files)

        stats.populate_from_genes(genes)
        _write_run_stats(config, stats)
        # Default end-of-run genome-level report (best-effort).
        _maybe_write_report(genes, config, inputs, stats)
        _log_summary(loci, mikado_loci, genes, stats)
        return genes
    finally:
        if genome is not None:
            genome.close()
        if h5 is not None:
            h5.close()


def _load_id_map(id_map_path: str | None) -> dict[str, str] | None:
    if id_map_path is None:
        return None
    p = Path(id_map_path)
    if p.exists():
        data: dict[str, str] = json.loads(p.read_text())
        return data
    return None


def _save_id_map(id_map_path: str | None, id_map: dict[str, str]) -> None:
    # Atomic (Phase 22, §3.2): a truncated id_map.json silently breaks HFG ID
    # stability, so it is never written in place.
    if id_map_path is None:
        return
    with atomic_write(id_map_path) as fh:
        fh.write(json.dumps(id_map, indent=0, sort_keys=True))


def _dump_resolved_config(config: PipelineConfig) -> None:
    """Write ``work_dir/run.resolved.yaml`` (D2, §3.6); never fatal."""
    try:
        if config.work_dir is None:
            return
        work = Path(config.work_dir)
        work.mkdir(parents=True, exist_ok=True)
        config.to_yaml(work / "run.resolved.yaml")
    except Exception as exc:  # noqa: BLE001 - observability never sinks a run
        _log.warning("could not write resolved config: %s", exc)


def _write_run_stats(config: PipelineConfig, stats: RunStats) -> None:
    """Write the run-telemetry sidecar ``<output_prefix>.run_stats.json`` (D3)."""
    path = f"{config.output_prefix}.run_stats.json"
    try:
        with atomic_write(path) as fh:
            fh.write(json.dumps(stats.as_dict(), indent=2, sort_keys=True))
    except OSError as exc:
        _log.warning("could not write run stats %s: %s", path, exc)


def _log_summary(
    loci: list[HelixerLocus],
    mikado_loci: list[MikadoLocus],
    genes: list[ReconciledGene],
    stats: RunStats | None = None,
) -> None:
    from collections import Counter

    tier: Counter[int] = Counter(g.tier for g in genes)
    origin: Counter[str] = Counter(g.origin for g in genes)
    total_tx = sum(len(g.transcripts) for g in genes)
    multi = sum(1 for g in genes if len(g.transcripts) > 1)
    _log.info(
        "pipeline complete: %d helixer loci, %d mikado loci → %d genes "
        "(%d isoforms, %d multi-isoform); tier=%s origin=%s",
        len(loci),
        len(mikado_loci),
        len(genes),
        total_tx,
        multi,
        dict(sorted(tier.items())),
        dict(sorted(origin.items())),
    )
    if stats is not None:
        # One structured decision-telemetry block: exactly the
        # numbers a bioinformatician needs to trust the run.
        h = stats.headline()
        _log.info(
            "run telemetry: merges accepted=%d rejected=%d; splits=%d; "
            "backstop rescued miniprot=%d transdecoder=%d none=%d; "
            "junction corrections applied=%d reverted=%d; "
            "isoforms dropped redundant=%d; partial ORFs=%d",
            h["merges_accepted"],
            h["merges_rejected"],
            h["splits"],
            h["backstop_rescued_miniprot"],
            h["backstop_rescued_transdecoder"],
            h["backstop_rescued_none"],
            h["junction_corrections_applied"],
            h["junction_corrections_reverted"],
            h["isoforms_dropped_redundant"],
            h["partial_orfs"],
        )
