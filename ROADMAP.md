# HelixForge Development Roadmap

## Overview
Track development progress across all phases. Update status after completing each milestone.

## Phase Status Legend
- [ ] Not started
- [~] In progress
- [x] Complete

---

## Phase 0: Project Scaffold
Status: [x] Complete

Deliverables:
- [x] pyproject.toml with dependencies
- [x] Directory structure with module stubs
- [x] CLAUDE.md project context
- [x] README.md
- [x] Pre-commit configuration

---

## Phase 1: I/O Foundation
Status: [x] Complete

Deliverables:
- [x] io/hdf5.py - Chunked parallel HDF5 reader
- [x] io/gff.py - GFF3 parser and writer
- [x] io/fasta.py - Indexed FASTA access
- [x] io/bam.py - RNA-seq junction extraction
- [x] Unit tests for all I/O modules
- [x] Integration test with synthetic data

Key requirements:
- Configurable chunk sizes for HPC memory management
- Thread-safe parallel reading
- Memory-mapped options for large files
- Streaming/iterator patterns to avoid loading full files

Implementation notes:
- All modules use 0-based, half-open coordinates internally
- GFF3 coordinates converted from 1-based on read/write
- Thread-safe HDF5 access with threading.Lock
- CoordinateIndex class maps HDF5 array positions to genomic coords
- attrs used for data classes (GeneModel, TranscriptModel)
- NamedTuple used for immutable data (SpliceJunction, CoverageProfile)
- pyfaidx wrapper for indexed FASTA access
- pysam for BAM parsing with CIGAR-based junction extraction

---

## Phase 2: Confidence Scoring
Status: [x] Complete

Deliverables:
- [x] core/confidence.py - Per-gene confidence metrics
- [x] Confidence output TSV format
- [x] Confidence output BED format (for genome browser)
- [x] Visualization of confidence distributions
- [x] CLI confidence subcommand
- [x] Unit tests

Metrics implemented:
- mean_prob: Mean probability of called class across gene
- min_prob: Minimum probability (weakest point)
- median_prob: Median probability across gene
- entropy: Shannon entropy (prediction uncertainty)
- boundary_sharpness: Transition sharpness at exon/intron boundaries
- coding_consistency: CDS frame consistency score
- exon_min: Worst-scoring exon

Implementation notes:
- GeneConfidence attrs class with all metrics and overall weighted score
- RegionConfidence for per-base visualization data
- ConfidenceCalculator supports parallel processing via ThreadPoolExecutor
- ConfidenceWriter outputs TSV, BED, DataFrame, low-confidence regions BED
- Confidence classification: high (>=0.85), medium (>=0.70), low (<0.70)
- Flags identify specific issues: weak_exon, high_entropy, uncertain_boundary, etc.
- viz/locus.py plot_gene_confidence() for per-gene interactive/static plots
- viz/genome.py plot_confidence_distribution(), plot_confidence_by_scaffold()
- CLI: helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa -o scores.tsv
- Supports --bed, --low-conf-bed, --plot-dir, --distribution-plot options

---

## Phase 3: Splice Refinement
Status: [x] Complete

Deliverables:
- [x] core/splice.py - Junction matching and correction
- [x] Position weight matrix scoring
- [x] Canonical/non-canonical site handling (GT-AG, GC-AG, AT-AC)
- [x] Unit tests with synthetic junctions
- [x] core/boundaries.py - Start/stop codon adjustment
- [x] data/pwm_plant.json - Plant splice site PWM data
- [x] CLI splice subcommand
- [x] viz/locus.py - splice refinement visualization
- [x] viz/genome.py - genome-wide splice summary
- [x] tests/test_splice.py - unit tests
- [x] tests/test_splice_integration.py - integration tests
- [x] docs/splice.md - documentation

Data structures implemented:
- SpliceSiteType enum (GT_AG, GC_AG, AT_AC, OTHER)
- SpliceSite attrs class with canonical detection
- IntronModel attrs class with splice type and RNA-seq support
- SpliceCorrection attrs class for tracking changes
- GeneSpliceReport attrs class with per-gene statistics
- PositionWeightMatrix class with from_sequences(), score(), load_plant_defaults()

Implementation notes:
- SpliceRefiner matches introns to RNA-seq junctions within configurable max_shift
- PWM scoring uses log-likelihood ratio against background model
- BoundaryAdjuster finds optimal start/stop codons with Kozak scoring
- Parallel gene processing via ThreadPoolExecutor
- SpliceReportWriter outputs TSV and computes summary statistics
- Plant PWMs derived from Arabidopsis/rice/maize consensus sequences

CLI usage:
```
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    --output-gff refined.gff3 \
    --report splice_report.tsv \
    --max-shift 15 \
    --min-reads 3 \
    --adjust-boundaries \
    --workers 4
```

---

## Phase 4: Parallel Execution
Status: [x] Complete (Revised)

Deliverables:
- [x] parallel/chunker.py - Genome partitioning strategies
- [x] parallel/executor.py - Multiprocessing wrapper with memory monitoring
- [x] parallel/taskgen.py - Task file generation for HyperShell/GNU Parallel
- [x] parallel/slurm.py - Minimal SLURM utilities (environment detection, example scripts)
- [x] templates/chunk_wrapper.sh - Generic wrapper script template
- [x] CLI parallel subcommand with plan, tasks, aggregate, example-sbatch, suggest
- [x] tests/test_parallel_*.py - Comprehensive test suite
- [x] docs/parallel.md - Documentation
- [x] **Chunk-aware flags in core commands (--region, --chunk-id, --scaffold)**
- [x] utils/regions.py - Region parsing and validation utilities
- [x] tests/test_utils_regions.py - Region parsing tests
- [x] tests/test_cli_chunk_aware.py - Chunk-aware CLI tests

Data structures implemented:
- ChunkStrategy enum (BY_SCAFFOLD, BY_SIZE, BY_GENES, ADAPTIVE)
- GenomicChunk attrs class with overlap/containment methods
- ChunkPlan attrs class with JSON serialization
- GenomeChunker class with strategy-specific partitioning
- ExecutorBackend enum (SERIAL, THREADS, PROCESSES)
- MemoryStats, TaskResult, ExecutionStats attrs classes
- MemoryMonitor class for background memory tracking
- ParallelExecutor class with multi-backend support
- ChunkProcessor high-level processing interface
- TaskFile attrs class for task file metadata
- TaskGenerator class for HyperShell/GNU Parallel task files

Chunking strategies:
- BY_SCAFFOLD: One chunk per scaffold, with optional max size splitting
- BY_SIZE: Fixed base-pair windows with scaffold boundary handling
- BY_GENES: Fixed gene count per chunk, respects scaffold boundaries
- ADAPTIVE: Automatic balancing based on gene density and target chunk count

Execution features:
- Serial, threaded, and process-based backends
- Progress callback support
- Error handling with continue_on_error option
- Memory monitoring with thresholds and logging
- Optimal worker count calculation

HPC parallelization approach:
- Task file generation (one command per line) for scheduler-agnostic execution
- Works with HyperShell, GNU Parallel, xargs, or custom scripts
- Wrapper script generation for complex multi-command workflows
- Minimal SLURM utilities for users who need direct SLURM integration
- Example SBATCH scripts for HyperShell and GNU Parallel

CLI usage:
```
# Create chunking plan
helixforge parallel plan --genome genome.fa --strategy scaffold -o chunks.json

# Generate task file
helixforge parallel tasks --chunk-plan chunks.json \
    --command 'helixforge confidence --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o {output_dir}/{chunk_id}.tsv' \
    --output tasks.txt --output-dir outputs/

# Execute with HyperShell (recommended)
hs launch --parallelism 32 < tasks.txt

# Or execute with GNU Parallel
parallel -j 32 < tasks.txt

# Aggregate results
helixforge parallel aggregate --input-dir outputs/ --pattern '*.tsv' -o combined.tsv --type merge_tsv

# Get example SLURM script
helixforge parallel example-sbatch -o run.sbatch --executor hypershell

# Get parameter suggestions
helixforge parallel suggest --genome genome.fa --memory 64 --workers 16
```

Design rationale:
- SLURM configurations vary significantly across HPC sites
- HyperShell provides scheduler-agnostic parallelization (SLURM, PBS, SGE, LSF)
- Task files are simple, inspectable, and portable
- Users configure HyperShell once; HelixForge generates task files

Chunk-aware command notes:
- Commands (`confidence`, `splice`) accept --region, --chunk-id, --scaffold flags
- Region format: seqid:start-end (1-based inclusive, e.g., chr1:1000-2000)
- Internal coordinates: 0-based half-open (Python convention)
- Task template {start}/{end}: 1-based for CLI commands
- Task template {start_0}/{end_0}: 0-based for internal use
- parse_region() converts CLI input to internal coords
- format_command() converts internal coords back to CLI format

---

## Phase 5: Homology Pipeline
Status: [x] Complete

Deliverables:
- [x] homology/search.py - Diamond/MMseqs2 wrapper
- [x] homology/validate.py - Coverage and identity scoring
- [x] homology/databases.py - Database download and management
- [x] TE overlap detection
- [x] Chimera detection (fused genes)
- [x] Fragment detection (split genes)
- [x] CLI homology subcommand group
- [x] Unit tests (tests/test_homology.py)
- [ ] homology/domains.py - InterProScan integration (optional, deferred)

Data structures implemented:
- SearchTool enum (DIAMOND, MMSEQS2)
- HomologyHit attrs class with coverage calculations
- GeneHomology attrs class for aggregated results
- HomologyStatus enum (COMPLETE, PARTIAL, NO_HIT, TE_OVERLAP, CHIMERIC, FRAGMENTED)
- ValidationThresholds class with presets (default, strict, relaxed)
- ValidationResult attrs class with detailed metrics
- ChimericEvidence NamedTuple for fusion detection
- FragmentGroup attrs class for split gene detection
- TEInterval attrs class for TE overlap
- DatabaseInfo dataclass for database management
- DatabaseType enum for supported databases

Implementation notes:
- HomologySearch class supports Diamond and MMseqs2 with format_database(), search(), parse_results()
- HomologyValidator detects chimeras (non-overlapping hits to different proteins) and fragments
- Fragment detection identifies multiple genes hitting same subject in non-overlapping regions
- TE overlap calculated from BED file annotations
- Protein extraction and translation from GFF3/genome
- Database download utilities for Swiss-Prot and UniRef

CLI usage:
```
# Extract proteins from gene models
helixforge homology extract-proteins -g genes.gff3 --genome genome.fa -o proteins.fa

# Format database for Diamond
helixforge homology format-db -i swissprot.fa -o swissprot

# Download Swiss-Prot plant proteins
helixforge homology download-db --database swissprot_plants

# Run homology search
helixforge homology search -q proteins.fa -d swissprot.dmnd -o hits.tsv

# Validate genes using homology
helixforge homology validate -g genes.gff3 -s hits.tsv -o validation.tsv \
    --te-bed te.bed --chimera-report chimeras.tsv --fragment-report fragments.tsv

# List available databases
helixforge homology list-databases
```

---

## Phase 6: QC System
Status: [x] Complete

Deliverables:
- [x] qc/flags.py - QC flag definitions (QCFlag, Flags registry, GeneQC)
- [x] qc/aggregate.py - QCAggregator for combining module results
- [x] qc/filters.py - GeneFilter with preset profiles
- [x] qc/report.py - HTML report generation with Chart.js
- [x] CLI qc subcommand group (aggregate, report, filter, tiered-output, list-flags)
- [x] tests/test_qc.py - Comprehensive test suite

Data structures implemented:
- FlagSeverity enum (INFO, WARNING, ERROR, CRITICAL) with comparison operators
- FlagCategory enum (CONFIDENCE, SPLICE, HOMOLOGY, STRUCTURE, ANNOTATION)
- QCFlag attrs class (frozen, hashable for use in sets)
- Flags registry class with 23 predefined flags
- GeneQC attrs class with tier classification and flag management
- QCAggregatorConfig for threshold configuration
- QCAggregator for combining confidence/splice/homology results
- FilterCriteria for custom filter definitions
- FilterResult for filter operation results
- GeneFilter with preset profiles (high_confidence, publication_ready, etc.)
- QCReportGenerator with embedded HTML/Chart.js template
- ReportConfig and ReportData for report customization

Tier classification logic:
- reject: CRITICAL flags or confidence < 0.50
- low: ERROR flags or confidence < 0.70
- medium: WARNING flags or confidence < 0.85
- high: No significant flags and confidence >= 0.85

CLI usage:
```
# Aggregate results from modules
helixforge qc aggregate --confidence-tsv scores.tsv --splice-tsv splice.tsv --homology-tsv validation.tsv -o qc_results.tsv

# Generate HTML report
helixforge qc report --qc-tsv qc_results.tsv -o report.html

# Filter genes with presets
helixforge qc filter --qc-tsv qc_results.tsv --preset publication_ready -o filtered.txt --gene-list-only

# Create tiered output
helixforge qc tiered-output --qc-tsv qc_results.tsv -o tiered/ --input-gff genes.gff3

# List all flags
helixforge qc list-flags --category homology
```

---

## Phase 7: Visualization
Status: [ ] Not started

Deliverables:
- [ ] viz/locus.py - Single-gene plots
- [ ] viz/genome.py - Genome-wide summaries
- [ ] viz/interactive.py - Optional dashboard (stretch goal)

---

## Phase 8: CLI Integration
Status: [ ] Not started

Deliverables:
- [ ] Full subcommand implementation
- [ ] Config file support (YAML)
- [ ] Progress bars and logging
- [ ] Shell completion

---

## Phase 9: Testing and Documentation
Status: [ ] Not started

Deliverables:
- [ ] Comprehensive pytest suite
- [ ] Test coverage >80%
- [ ] mkdocs documentation site
- [ ] Tutorial notebooks
- [ ] Example workflows

---

## Phase 10: Isoform Support (v0.2)
Status: [ ] Not started

Deliverables:
- [ ] isoforms/evidence.py - Alternative splicing detection
- [ ] isoforms/reconstruct.py - StringTie integration
- [ ] isoforms/select.py - Primary transcript selection
- [ ] Updated GFF3 output with isoforms

---

## Notes
- Update this file after completing each deliverable
- Add discovered tasks as sub-items
- Record blockers or design decisions in phase notes
