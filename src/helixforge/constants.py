"""Central home for HelixForge's tuning constants."""

from __future__ import annotations

# --- HFG id allocation (reconcile/mikado_integrate.py) ---
# Novel (Mikado-only) genes draw HFG numbers from this high range, outside the
# Helixer-anchored block, so they never collide with Helixer-derived ids.
NOVEL_ID_BASE = 90000
# Structured-ncRNA loci added by the optional tRNAscan-SE / Infernal hook
# draw HFG numbers from a separate high block above NOVEL_ID_BASE,
# so they never collide with Helixer-anchored ids or Mikado-novel ids.
NCRNA_ID_BASE = 95000

# --- genetic code (utils/sequences.py) --- [BIOLOGY: per-species/organelle override]
# Default NCBI ``transl_table`` id. 1 = the standard code used by nuclear plant
# genes; organellar (plastid/mito) contigs can be remapped per-scaffold without
# changing this default.
DEFAULT_TRANSL_TABLE = 1

# --- CDS cross-check (reconcile/cds.py) --- [BIOLOGY: per-species override candidate]
# Below this reciprocal CDS overlap, an independent miniprot ORF is judged to
# materially disagree with the chosen Mikado ORF (cross-check flag only; never
# edits the CDS). A homology/ORF-coherence threshold, not a universal constant.
CROSS_CHECK_OVERLAP = 0.8

# --- backstop structural floors (reconcile/fallback.py) ---
# Minimum exon length the junction-correction transaction will accept; below it
# the correction reverts (a flagged gene always beats a destroyed one).
MIN_EXON_BP = 3
# [BIOLOGY: per-species override candidate] Minimum intron length — a splicing
# assumption. A "correction" producing a sub-20 bp gap is rejected as implausible.
MIN_INTRON_BP = 20

# --- read-filtering policy (io/bam.py) --- [BIOLOGY: paralog/polyploid sensitive]
# Junction extraction keeps only uniquely-mapped
# reads (min_mapq>=10 with STAR's 255/3/1/0 scheme), but coverage historically
# applied no MAPQ filter and the default htslib max_depth (~8000) silently capped
# deep pileups. These knobs make the policy consistent and configurable.
# MAPQ threshold above which a read counts as uniquely mapped (STAR: 255=unique,
# 3=2-mapper, 1, 0). Used for the "unique-only" coverage mode and to mirror the
# junction min_mapq default.
UNIQUE_MIN_MAPQ = 10
# Explicit pileup depth ceiling (replaces the silent htslib ~8000 default) so a
# highly-expressed gene's coverage is not truncated. Large but bounded to keep a
# pathological pileup from exhausting memory; raise per-run if needed.
COVERAGE_MAX_DEPTH = 1_000_000

# --- paralog/tandem-array merge guards (reconcile/mikado_integrate.py) ---
# [BIOLOGY: paralog/tandem-array sensitive] A Mikado locus that bridges two
# adjacent Helixer loci is only accepted as a genuine fusion when the bridging
# intron is (a) canonical and (b)
# backed by a splice junction with >= this many reads spanning the *inter-genic*
# gap specifically. Guards against a chimeric/read-through transcript fusing two
# real genes in NBS-LRR / F-box tandem arrays. Matches the junction-extraction
# min-reads default so a merge needs the same evidence floor as any other intron.
MERGE_MIN_GAP_READS = 3
# k-mer width for the alignment-free pairwise CDS-identity proxy used to detect
# recent duplications / homeologs (a high-identity adjacent pair is unlikely to
# be a single split gene and is *not* fused). Pure engineering knob.
PARALOG_IDENTITY_K = 12
# [BIOLOGY: per-species override candidate] Default k-mer Jaccard identity above
# which two adjacent Helixer loci are treated as recent paralogs/homeologs and
# their merge is rejected. Only consulted when the per-run threshold is set
# (config default None = identity guard disabled, so the default run is unchanged).
PARALOG_IDENTITY_THRESHOLD = 0.9

# --- non-coding biotype classification (reconcile/biotype.py) ---
# [BIOLOGY] Minimum spliced-transcript length (nt) for
# a lncRNA call — the GENCODE/Ensembl operational lncRNA floor. An expressed,
# multi-exonic, ORF-less locus shorter than this is left ncRNA_undetermined rather
# than asserted to be a long non-coding RNA.
LNCRNA_MIN_LENGTH = 200
# [BIOLOGY: per-species override candidate] Helixer CDS-channel probability
# (channel 2, exon-length-weighted) at/below which an ORF-less locus is treated as
# coding-channel-quiet — i.e. consistent with a genuine non-coding RNA rather than
# a fragmentary/failed coding gene. An ORF-less locus with CDS-channel confidence
# *above* this (Helixer still "saw" coding signal but no ORF was admissible) is
# left ncRNA_undetermined, not asserted to be a lncRNA. Consulted only when an HDF5
# reader is available; with no HDF5 the signal is absent and does not block.
NONCODING_MAX_CDS_CONF = 0.5

# --- scatter-gather boundary gap (parallel/plan.py) ---
# Heuristic floor for a cut-eligible inter-locus gap when the caller gives none;
# wide enough to clear a typical gene / Mikado-merge bridging span.
PLAN_DEFAULT_MIN_BOUNDARY_GAP = 1000

# --- small-scaffold bin-packing (parallel/plan.py) ---
# A fragmented draft assembly can have 10^5-10^6 tiny contigs; one chunk per
# scaffold would blow past Slurm MaxArraySize. When packing is enabled, a
# scaffold at or below this length is "small" and eligible to share a chunk with
# other small scaffolds (a gene is never split — only whole, uncut scaffolds are
# packed). 1 Mb is a generous ceiling for a draft contig vs a real chromosome arm.
PLAN_SMALL_SCAFFOLD_BP = 1_000_000
# Cap on how many small scaffolds a single packed chunk may hold, so a combined
# chunk stays a reasonable unit of work even when target_loci_per_chunk is unset.
PLAN_MAX_SCAFFOLDS_PER_CHUNK = 200

# --- Helixer-support weighting (mikado/emit_external.py) ---
# Default weight of the exon-confidence component (i) vs the junction-coincidence
# component (ii) when combining them into the [0, 1] ``helixer_support`` metric.
DEFAULT_EXON_WEIGHT = 0.7
# [BIOLOGY: per-species override candidate] Helixer intron-probability cutoff
# above which a transcript intron counts as "Helixer-supported" in component (ii).
DEFAULT_INTRON_HIGH_CUTOFF = 0.5

# --- subprocess error reporting (prep/_subprocess.py) ---
# Number of trailing stderr characters kept in a failed-tool exception message.
STDERR_TAIL = 2000

# --- granularity / resource heuristics (parallel/suggest.py) ---
# Labelled heuristics, not promises (tune per cluster / genome).
TARGET_CHUNK_BP = 40_000_000  # ~40 Mb of sequence per chunk
MEM_GB_PER_MB_SEQUENCE = 0.05  # ~50 MB RAM per Mb of chunk sequence
MIN_CHUNK_MEM_GB = 4  # never request less than this per chunk
WALLTIME_MIN_PER_MB = 1.5  # ~1.5 wall-minutes per Mb of chunk sequence
MIN_WALLTIME_MIN = 60  # never request less than this walltime
# Conservative cut-eligible gap floor used by ``suggest`` (clears most Mikado
# merges); distinct from the ``plan`` floor above, which is the planner's own.
SUGGEST_DEFAULT_MIN_BOUNDARY_GAP = 2000
