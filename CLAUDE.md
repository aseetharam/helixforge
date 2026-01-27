# CLAUDE.md - HelixForge Development Guide

## Project Vision and Goals

HelixForge is a modular toolkit designed to transform Helixer's raw gene predictions into publication-quality genome annotations. While Helixer provides excellent initial predictions using deep learning, real-world applications require refinement using empirical evidence, quality control, and expert validation.

### Core Objectives

1. **Evidence Integration**: Incorporate RNA-seq alignments, protein homology, and domain annotations to validate and refine gene models
2. **Confidence Scoring**: Provide transparent, interpretable confidence metrics for each predicted gene
3. **Isoform Reconstruction**: Use splice junction evidence to reconstruct alternative transcripts (v0.2)
4. **Quality Control**: Flag problematic predictions and generate comprehensive QC reports
5. **HPC-Ready**: Scale from single genomes to pan-genome datasets on cluster infrastructure
6. **Reproducibility**: Ensure all analyses are configurable, logged, and reproducible

## Current Development Phase

| Phase | Description | Status |
|-------|-------------|--------|
| 0 | Project scaffold | COMPLETE |
| 1 | I/O Foundation | COMPLETE |
| 2 | Confidence Scoring | COMPLETE |
| 3 | Splice Refinement | COMPLETE |
| 4 | Parallel Execution | COMPLETE |
| 5 | Homology Pipeline | COMPLETE |
| 6 | QC System | COMPLETE |
| 7 | Visualization polish | PENDING |
| 8 | CLI integration | PENDING |
| 9 | Testing and documentation | PENDING |
| 10 | Isoform support (v0.2) | PENDING |

**Current Status**: Core v0.1 functionality complete. Ready for real-world testing and refinement.

## Key Design Decisions

### HDF5 Requirement

HelixForge requires Helixer's HDF5 predictions file (not just GFF3) to access per-base confidence scores. Users must preserve this file during Helixer runs by running internal scripts separately:

```bash
# Step 1: Convert FASTA to H5
apptainer exec ${sif} fasta2h5.py \
    --species ${species} \
    --h5-output-path ${outdir}/${species}_input.h5 \
    --fasta-path ${genome}

# Step 2: Run predictions (preserves H5)
apptainer exec --nv ${sif} HybridModel.py \
    --load-model-path ${model_path} \
    --test-data ${outdir}/${species}_input.h5 \
    --prediction-output-path ${outdir}/${species}_predictions.h5 \
    --batch-size 8

# Step 3: Generate GFF3
apptainer exec ${sif} helixer_post_bin \
    ${outdir}/${species}_predictions.h5 \
    ${outdir}/${species}_helixer.gff3
```

### Parallel Execution: HyperShell over SLURM

We generate simple task files (one command per line) rather than SLURM-specific scripts because:

- SLURM configurations vary significantly across HPC sites
- HyperShell provides scheduler-agnostic parallelization
- Users configure their environment once; HelixForge generates portable task lists

```bash
# HelixForge generates task file
helixforge parallel tasks --chunk-plan chunks.json --command "..." --output tasks.txt

# User runs with their preferred tool
hs launch --parallelism 50 < tasks.txt  # HyperShell (recommended)
parallel -j 32 < tasks.txt               # GNU Parallel
```

### Chunk-aware Processing

Core commands (confidence, splice, homology) support region-based processing for parallelization:

| Flag | Description | Example |
|------|-------------|---------|
| `--region` | Process specific region (1-based) | `--region chr1:1000-2000` |
| `--scaffold` | Process entire scaffold | `--scaffold chr1` |
| `--chunk-id` | Identifier for logging/output | `--chunk-id chunk_042` |

### Coordinate Conventions

| Context | Convention | Example |
|---------|------------|---------|
| CLI input | 1-based inclusive | chr1:1000-2000 |
| Internal storage | 0-based half-open | (999, 2000) |
| GFF3 files | 1-based inclusive | 1000\t2000 |
| BED files | 0-based half-open | 999\t2000 |
| HDF5 arrays | 0-based | index 999 |

### QC Tier System

Genes are classified into four tiers based on aggregated evidence:

| Tier | Description | Typical criteria |
|------|-------------|------------------|
| high | Publication-ready | Confidence ≥0.85, homology support, RNA-seq support |
| medium | Usable with caveats | Confidence ≥0.70, no critical flags |
| low | Requires review | Confidence <0.70 or error-level flags |
| reject | Likely false positive | TE overlap, internal stops, chimeric |

### QC Flag System

Flags organized by category and severity:

| Category | Example flags |
|----------|---------------|
| confidence | LOW_CONF, VERY_LOW_CONF, WEAK_EXON, UNCERTAIN_BOUND |
| splice | SPLICE_CORR, NO_RNASEQ, NONCANON_SPLICE |
| homology | NO_HOMOL, PARTIAL_HOMOL, TE_OVERLAP, CHIMERIC, FRAGMENTED |
| structure | NO_START, NO_STOP, INTERNAL_STOP, SHORT_INTRON |

Severities: INFO < WARNING < ERROR < CRITICAL

### Evidence Handling (Future)

Under consideration for `--add-evidence` flag:

- Embed minimal evidence as GFF3 attributes (junction counts, best hit)
- Keep full evidence tracks as separate export option
- Avoid bloating primary GFF3 output

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                          CLI (cli.py)                           │
│              User-facing commands and orchestration             │
└─────────────────────────────────────────────────────────────────┘
                                │
        ┌───────────────────────┼───────────────────────┐
        ▼                       ▼                       ▼
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│   IO Module   │       │  Core Module  │       │   Parallel    │
│ hdf5,gff,fasta│       │ confidence,   │       │ chunker,exec, │
│     bam       │       │ splice,bounds │       │ taskgen,slurm │
└───────────────┘       └───────────────┘       └───────────────┘
        │                       │                       │
        │       ┌───────────────┼───────────────┐       │
        │       ▼               ▼               ▼       │
        │ ┌───────────┐   ┌───────────┐   ┌───────────┐ │
        │ │ Homology  │   │ Isoforms  │   │    QC     │ │
        │ │  search,  │   │ (v0.2)    │   │  flags,   │ │
        │ │ validate  │   │           │   │ aggregate,│ │
        │ │           │   │           │   │  report   │ │
        │ └───────────┘   └───────────┘   └───────────┘ │
        │       │               │               │       │
        └───────┴───────────────┴───────────────┴───────┘
                                │
                                ▼
                    ┌───────────────────┐
                    │   Viz Module      │
                    │ locus, genome     │
                    └───────────────────┘
```

### Module Responsibilities

#### `io/` - Input/Output Handlers

- **hdf5.py**: Chunked, parallel HDF5 reader for Helixer predictions
- **gff.py**: GFF3 parsing (GFF3Parser) and writing (GFF3Writer) with attribute preservation
- **fasta.py**: Indexed FASTA access via pyfaidx (GenomeAccessor)
- **bam.py**: RNA-seq junction extraction (JunctionExtractor), multi-sample aggregation, STAR SJ.out.tab support

#### `core/` - Core Refinement Logic

- **refine.py**: Main refinement pipeline (RefinePipeline) - combines all components
- **confidence.py**: Per-gene confidence scoring from HDF5 probabilities
- **splice.py**: Splice site refinement using RNA-seq junctions and PWM scoring
- **boundaries.py**: Start/stop codon adjustment
- **evidence.py**: RNA-seq evidence scoring (EvidenceScorer, AED calculation)
- **evidence_output.py**: Evidence report writers (TSV, GFF3 updates)

#### `homology/` - Homology-Based Validation

- **search.py**: Diamond/MMseqs2 interface, protein extraction
- **validate.py**: Coverage-based validation, chimera/fragment detection
- **domains.py**: Optional InterProScan integration

#### `qc/` - Quality Control

- **flags.py**: QC flag definitions (Flags registry, GeneQC)
- **aggregate.py**: Combine results from all modules (QCAggregator)
- **filters.py**: Gene filtering by criteria/presets (GeneFilter)
- **report.py**: HTML report generation (QCReportGenerator)

#### `parallel/` - Parallelization

- **chunker.py**: Genome partitioning (GenomeChunker, ChunkPlan)
- **executor.py**: Local multiprocessing (ParallelExecutor)
- **taskgen.py**: Task file generation for HyperShell/GNU Parallel (TaskGenerator)
- **slurm.py**: Minimal SLURM utilities (environment detection only)

#### `viz/` - Visualization

- **locus.py**: Single-gene plots with evidence tracks
- **genome.py**: Genome-wide summary plots

#### `utils/` - Utilities

- **regions.py**: Region parsing and coordinate handling
- **intervals.py**: Genomic interval operations
- **sequences.py**: Translation, codon tables
- **logging.py**: Structured logging

## Commands Reference

### Refine Pipeline (Recommended)

The `refine` command is the primary pipeline that combines splice correction, boundary adjustment, confidence scoring, and evidence scoring:

```bash
# Basic refinement (single BAM)
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r refine_report.tsv

# Multi-tissue with full reports
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam liver.bam,brain.bam,heart.bam \
    --min-tissues 2 \
    --min-reads 3 \
    -o refined.gff3 \
    -r refine_report.tsv \
    --splice-details splice_corrections.tsv

# Using STAR junction files
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --junctions-bed sample1_SJ.out.tab,sample2_SJ.out.tab \
    -o refined.gff3

# Parallel processing by region
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    --region chr1:1-10000000 \
    --chunk-id chunk_001 \
    -o chunk_001.gff3
```

**Output attributes:** `confidence`, `evidence_level`, `aed` (MAKER-compatible 0-1), `junction_support`, `splice_corrections`, `refine_flags`.

### Confidence Scoring (Standalone)

For confidence scoring only without refinement:

```bash
helixforge confidence \
    -p helixer_predictions.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    -o confidence.tsv \
    --bed confidence.bed \
    --low-conf-bed low_confidence_regions.bed \
    --threshold 0.7 \
    --region chr1:1-10000000 \
    --chunk-id chunk_001 \
    -j 8
```

### Evidence Scoring (Standalone)

For evidence scoring only without refinement (no HDF5 required):

```bash
# Basic evidence scoring
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3

# Multiple BAMs with detailed reports
helixforge evidence \
    -g predictions.gff3 \
    -b tissue1.bam,tissue2.bam,tissue3.bam \
    -o annotated.gff3 \
    --report evidence_report.tsv \
    --summary evidence_summary.tsv \
    --junction-details junction_details.tsv \
    --exon-details exon_details.tsv

# Junction-only mode (faster, no coverage analysis)
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3 \
    --no-coverage
```

**Output attributes:** `evidence_level` (full/partial/minimal/none), `aed` (0-1), `junction_support`, `exon_coverage`, `evidence_flags`.

### Splice Refinement (DEPRECATED)

> **Note:** The standalone `splice` command is deprecated. Use `helixforge refine` instead.

```bash
# DEPRECATED - use 'helixforge refine' instead
helixforge splice \
    --helixer-gff helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r splice_report.tsv
```

### Homology Validation

```bash
# Extract proteins from gene models
helixforge homology extract-proteins \
    --gff helixer.gff3 \
    --genome genome.fa \
    -o proteins.fa \
    --longest-isoform

# Format database
helixforge homology format-db \
    --proteins swissprot_viridiplantae.fa \
    -o swissprot.dmnd \
    --tool diamond

# Run search
helixforge homology search \
    --proteins proteins.fa \
    --database swissprot.dmnd \
    -o hits.tsv \
    --threads 16 \
    --evalue 1e-5

# Validate gene models
helixforge homology validate \
    --search-results hits.tsv \
    --gff helixer.gff3 \
    -o validation.tsv \
    --te-regions repeats.bed \
    --fragments-output fragments.tsv \
    --chimeras-output chimeras.tsv
```

### QC Aggregation and Reporting

```bash
# Aggregate from refine output (recommended)
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv validation.tsv \
    -o qc_results.tsv

# Or aggregate from separate module outputs
helixforge qc aggregate \
    --confidence-tsv confidence.tsv \
    --splice-tsv splice_report.tsv \
    --homology-tsv validation.tsv \
    -o qc_results.tsv

# Generate HTML report
helixforge qc report \
    --qc-tsv qc_results.tsv \
    -o qc_report.html \
    --title "Genome Annotation QC"

# Filter genes by preset
helixforge qc filter \
    --gff helixer.gff3 \
    --qc-results qc_results.tsv \
    -o filtered.gff3 \
    --preset publication_ready

# Create tiered outputs
helixforge qc tiered-output \
    --gff helixer.gff3 \
    --qc-results qc_results.tsv \
    -o tiered_output/ \
    --prefix genes

# List available flags
helixforge qc list-flags --category confidence --severity warning
```

Filter presets: `high_confidence`, `publication_ready`, `has_homology`, `has_rnaseq_support`, `needs_review`, `reject`, `not_rejected`

### Parallel Execution

```bash
# Create chunk plan
helixforge parallel plan \
    --genome genome.fa \
    --gff helixer.gff3 \
    --strategy adaptive \
    --target-chunks 100 \
    -o chunks.json \
    --summary

# Generate task file
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge confidence -p helixer.h5 -g helixer.gff3 --genome genome.fa --region {seqid}:{start}-{end} --chunk-id {chunk_id} -o outputs/{chunk_id}.tsv" \
    -o tasks.txt

# Execute with HyperShell
hs launch --parallelism 64 < tasks.txt

# Or with GNU Parallel
parallel -j 32 < tasks.txt

# Aggregate outputs
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern "*.tsv" \
    -o combined.tsv \
    --type merge_tsv
```

Chunking strategies: `scaffold`, `size`, `genes`, `adaptive`

## Coding Conventions

### Type Hints

All public functions MUST have complete type annotations:

```python
def calculate_confidence(
    gene: GeneModel,
    predictions: np.ndarray,
    thresholds: ConfidenceThresholds | None = None,
) -> GeneConfidence:
    """Calculate confidence score for a gene model."""
    ...
```

### Docstrings

Use Google-style docstrings:

```python
def extract_splice_junctions(
    bam_path: Path,
    region: GenomicRegion,
    min_support: int = 3,
) -> list[SpliceJunction]:
    """Extract splice junctions from RNA-seq alignments.

    Args:
        bam_path: Path to indexed BAM file.
        region: Genomic region to extract junctions from.
        min_support: Minimum number of reads supporting a junction.

    Returns:
        List of SpliceJunction objects with coordinates and support counts.

    Raises:
        FileNotFoundError: If BAM file or index doesn't exist.
        ValueError: If region is invalid.
    """
    ...
```

### HPC-Aware Design

- Accept configurable chunk sizes
- Use generators for large datasets
- Provide memory estimates where possible
- Support region-based processing for parallelization

```python
def process_genes(
    genes: Iterable[GeneModel],
    chunk_size: int = 1000,
) -> Iterator[ProcessedGene]:
    """Process genes in memory-efficient chunks."""
    ...
```

### Naming Conventions

- Classes: `PascalCase` (GeneModel, SpliceJunction)
- Functions/methods: `snake_case` (calculate_confidence)
- Constants: `UPPER_SNAKE_CASE` (DEFAULT_CHUNK_SIZE)
- Private: `_leading_underscore`

## Dependencies

### Core Dependencies

| Package | Min Version | Purpose |
|---------|-------------|---------|
| h5py | 3.8 | HDF5 I/O |
| pysam | 0.21 | BAM parsing |
| pyfaidx | 0.7 | FASTA indexing |
| gffutils | 0.12 | GFF3 parsing |
| numpy | 1.24 | Numerical ops |
| pandas | 2.0 | Data manipulation |
| scipy | 1.10 | Statistics |
| matplotlib | 3.7 | Plotting |
| plotly | 5.14 | Interactive plots |
| jinja2 | 3.1 | Report templates |
| click | 8.1 | CLI framework |
| rich | 13.0 | Terminal output |
| attrs | 23.0 | Data classes |
| ncls | 0.0.68 | Interval queries |

### External Tools

| Tool | Required for | Installation |
|------|--------------|--------------|
| Diamond | Homology search | `conda install -c bioconda diamond` |
| samtools | BAM indexing | `conda install -c bioconda samtools` |
| HyperShell | HPC parallelization | `pip install hypershell` |

### Python Version

- Minimum: Python 3.10
- Tested: 3.10, 3.11, 3.12

## Known Limitations

### Helixer H5 File Preservation

Helixer deletes intermediate H5 files by default. Users must run Helixer's internal scripts separately to preserve the predictions file. See "HDF5 Requirement" section above.

### Memory Requirements

- HDF5 reading: ~20 bytes/base for predictions
- Large chromosomes (>200 Mb): Use chunked processing
- Full maize genome (~2.3 Gb): Recommend 64+ GB RAM or parallel chunks

### Single Isoform per Gene

Helixer outputs only primary transcripts. Isoform reconstruction planned for v0.2.

## Testing

### Running Tests

```bash
# All tests
pytest tests/ -v

# With coverage
pytest --cov=helixforge --cov-report=html

# Specific module
pytest tests/test_confidence.py -v
```

### Testing Checklist

Before release, verify:

```bash
# 1. Unit tests pass
pytest tests/ -v

# 2. Single chromosome test
helixforge confidence \
    -p test.h5 -g test.gff3 --genome test.fa \
    --scaffold chr1 -o test_conf.tsv

# 3. Parallel workflow
helixforge parallel plan --genome test.fa --strategy scaffold -o chunks.json
helixforge parallel tasks --chunk-plan chunks.json --command "echo {chunk_id}" -o tasks.txt
head -1 tasks.txt | bash

# 4. QC report renders
helixforge qc report --qc-results test_qc.tsv -o test_report.html
# Open in browser, verify charts

# 5. CLI help is accurate
helixforge --help
helixforge confidence --help
helixforge splice --help
helixforge homology --help
helixforge qc --help
helixforge parallel --help
```

## Next Development Targets

### v0.1 Polish (Phases 7-9)

1. Visualization improvements
2. CLI usability and error messages
3. Documentation and tutorials
4. Test coverage >80%

### v0.2 Features (Phase 10)

1. Isoform reconstruction from junction evidence
2. Alternative transcript scoring
3. Representative transcript selection

**Note:** Multi-tissue RNA-seq integration is now available in v0.1 via `--min-tissues` flag.

### Future Considerations

- Long-read RNA-seq support (IsoSeq, ONT)
- Protein-to-genome alignment (miniprot)
- Comparative annotation across genomes
- Integration with Apollo/WebApollo for manual curation
