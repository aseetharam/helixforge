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
Status: [ ] Not started

Deliverables:
- [ ] parallel/chunker.py - Genome partitioning strategies
- [ ] parallel/executor.py - Multiprocessing wrapper
- [ ] parallel/slurm.py - SLURM array job generator
- [ ] Scaling tests on multi-core systems

---

## Phase 5: Homology Pipeline
Status: [ ] Not started

Deliverables:
- [ ] homology/search.py - Diamond/MMseqs2 wrapper
- [ ] homology/validate.py - Coverage and identity scoring
- [ ] homology/domains.py - InterProScan integration (optional)
- [ ] TE overlap detection
- [ ] Unit tests

---

## Phase 6: QC System
Status: [ ] Not started

Deliverables:
- [ ] qc/flags.py - QC flag definitions
- [ ] qc/filters.py - Tiered filtering logic
- [ ] qc/report.py - HTML report generation
- [ ] Summary statistics output

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
