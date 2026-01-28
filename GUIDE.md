# CLAUDE.md - HelixForge Development Guide

## Project Vision and Goals

HelixForge is a modular toolkit designed to transform Helixer's raw gene predictions into publication-quality genome annotations. While Helixer provides excellent initial predictions using deep learning, real-world applications require refinement using empirical evidence, quality control, and expert validation.

### Core Objectives

1. **Evidence Integration**: Incorporate RNA-seq alignments, protein homology, and domain annotations to validate and refine gene models
2. **Confidence Scoring**: Provide transparent, interpretable confidence metrics for each predicted gene
3. **Isoform Reconstruction**: Use splice junction evidence to reconstruct alternative transcripts
4. **Quality Control**: Flag problematic predictions and generate comprehensive QC reports
5. **HPC-Ready**: Scale from single genomes to pan-genome datasets on cluster infrastructure
6. **Reproducibility**: Ensure all analyses are configurable, logged, and reproducible

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                          CLI (cli.py)                           │
│              User-facing commands and orchestration             │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────┐
│                       Config (config.py)                        │
│           Configuration management and validation               │
└─────────────────────────────────────────────────────────────────┘
                                │
        ┌───────────────────────┼───────────────────────┐
        ▼                       ▼                       ▼
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│   IO Module   │       │  Core Module  │       │   Parallel    │
│ hdf5,gff,fasta│       │ confidence,   │       │   chunker,    │
│     bam       │       │ splice,bounds │       │   executor    │
└───────────────┘       └───────────────┘       └───────────────┘
        │                       │
        │       ┌───────────────┼───────────────┐
        │       ▼               ▼               ▼
        │ ┌───────────┐   ┌───────────┐   ┌───────────┐
        │ │ Homology  │   │ Isoforms  │   │    QC     │
        │ │  search,  │   │ evidence, │   │  flags,   │
        │ │ validate  │   │reconstruct│   │  report   │
        │ └───────────┘   └───────────┘   └───────────┘
        │       │               │               │
        └───────┴───────────────┴───────────────┘
                                │
                                ▼
                    ┌───────────────────┐
                    │   Viz Module      │
                    │ locus, genome,    │
                    │   interactive     │
                    └───────────────────┘
```

### Module Responsibilities

#### `io/` - Input/Output Handlers
- **hdf5.py**: Read/write Helixer HDF5 predictions and intermediate results
- **gff.py**: Parse and write GFF3/GTF files with validation
- **fasta.py**: Efficient genome sequence access with pyfaidx
- **bam.py**: RNA-seq alignment parsing and coverage extraction

#### `core/` - Core Refinement Logic
- **confidence.py**: Multi-factor confidence scoring for gene models
- **splice.py**: Splice site analysis and canonical junction verification
- **boundaries.py**: Start/stop codon refinement and UTR extension
- **merge_split.py**: Gene model merging (fragments) and splitting (fusions)

#### `homology/` - Homology-Based Validation
- **search.py**: Interface to BLAST/DIAMOND for homology searches
- **validate.py**: Validate gene models against protein databases
- **domains.py**: InterProScan integration for domain annotation

#### `isoforms/` - Alternative Splicing
- **evidence.py**: Collect splice junction evidence from BAM files
- **reconstruct.py**: Build alternative transcript models
- **select.py**: Select representative transcripts per gene

#### `qc/` - Quality Control
- **flags.py**: Define and apply QC flags to gene models
- **report.py**: Generate HTML/PDF QC reports
- **filters.py**: Configurable filtering criteria

#### `viz/` - Visualization
- **locus.py**: Single-gene visualization with evidence tracks
- **genome.py**: Genome-wide statistics and plots
- **interactive.py**: Dash/Panel-based interactive browser

#### `parallel/` - Parallelization
- **chunker.py**: Divide genome into processable chunks
- **executor.py**: Local multiprocessing execution with memory monitoring
- **taskgen.py**: Task file generation for HyperShell/GNU Parallel
- **slurm.py**: Minimal SLURM utilities (environment detection, example scripts)

#### `utils/` - Utilities
- **intervals.py**: Genomic interval operations (overlap, merge, subtract)
- **sequences.py**: Sequence manipulation and translation
- **logging.py**: Structured logging configuration

## Current Development Phase

**Phase 0: Initialization** (Complete)
- [x] Project structure and configuration
- [x] Core data structures (GeneModel, Transcript, Exon)
- [x] Basic I/O for GFF3 and FASTA
- [x] CLI skeleton with help text

**Phase 1: Foundation** (Complete)
- [x] HDF5 reader for Helixer predictions
- [x] Confidence scoring framework
- [x] Basic QC flags
- [x] Unit test infrastructure

**Phase 2: Evidence Integration** (Complete)
- [x] BAM parsing for RNA-seq evidence
- [x] Splice junction extraction
- [x] Coverage-based confidence

**Phase 3: Splice Refinement** (Complete)
- [x] Splice site data structures (SpliceSiteType, SpliceSite, IntronModel, SpliceCorrection, GeneSpliceReport)
- [x] Position Weight Matrix (PWM) scoring for splice motifs
- [x] SpliceRefiner engine with RNA-seq junction matching
- [x] BoundaryAdjuster for start/stop codon optimization
- [x] Plant PWM data (Arabidopsis/rice/maize derived)
- [x] CLI `splice` subcommand
- [x] Splice visualization functions (locus and genome-wide)
- [x] Comprehensive test suite with synthetic fixtures
- [x] Documentation (docs/splice.md)

**Phase 4: Advanced Refinement** (Next)
- [ ] Gene merge/split detection
- [ ] Isoform reconstruction from junction evidence
- [ ] UTR extension using coverage data

## Coding Conventions

### Type Hints
All public functions MUST have complete type annotations:

```python
def calculate_confidence(
    gene: GeneModel,
    evidence: EvidenceCollection,
    weights: dict[str, float] | None = None,
) -> ConfidenceScore:
    """Calculate confidence score for a gene model."""
    ...
```

### Docstrings
Use Google-style docstrings for all public functions and classes:

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
        ValueError: If region is invalid or outside reference bounds.
    """
    ...
```

### HPC-Aware Design
- All processing functions should accept configurable chunk sizes
- Use generators for large datasets to minimize memory footprint
- Provide memory usage estimates where possible
- Support both streaming and batch modes

```python
def process_genes(
    genes: Iterable[GeneModel],
    chunk_size: int = 1000,
    max_memory_mb: int | None = None,
) -> Iterator[ProcessedGene]:
    """Process genes in memory-efficient chunks."""
    ...
```

### Parallel-Ready Functions
Design functions to work on genomic regions for easy parallelization:

```python
def refine_region(
    region: GenomicRegion,
    genome: FastaFile,
    predictions: HDF5File,
    evidence: EvidenceCollection,
) -> list[GeneModel]:
    """Refine all genes in a genomic region.

    This function is designed to be called in parallel across
    non-overlapping genomic regions.
    """
    ...
```

### Naming Conventions
- Classes: `PascalCase` (e.g., `GeneModel`, `SpliceJunction`)
- Functions/methods: `snake_case` (e.g., `calculate_confidence`)
- Constants: `UPPER_SNAKE_CASE` (e.g., `DEFAULT_CHUNK_SIZE`)
- Private attributes: `_leading_underscore`
- Module-level type aliases: `PascalCase` (e.g., `Strand = Literal['+', '-']`)

## Key Data Structures

### Core Types (to be implemented in `core/`)

```python
@attrs.define
class GenomicRegion:
    """A region on a chromosome."""
    seqid: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    strand: Strand = "+"

@attrs.define
class Exon:
    """A single exon within a transcript."""
    region: GenomicRegion
    exon_number: int
    phase: int | None = None

@attrs.define
class Transcript:
    """A transcript (mRNA) with exons and CDS."""
    transcript_id: str
    gene_id: str
    exons: list[Exon]
    cds_start: int | None = None
    cds_end: int | None = None
    attributes: dict[str, str] = attrs.Factory(dict)

@attrs.define
class GeneModel:
    """A gene with one or more transcripts."""
    gene_id: str
    region: GenomicRegion
    transcripts: list[Transcript]
    source: str = "helixer"
    gene_type: str = "protein_coding"
    confidence: ConfidenceScore | None = None
    qc_flags: list[QCFlag] = attrs.Factory(list)
    attributes: dict[str, str] = attrs.Factory(dict)

@attrs.define
class ConfidenceScore:
    """Multi-component confidence score."""
    overall: float  # 0.0 - 1.0
    components: dict[str, float]
    evidence_summary: dict[str, Any]
```

### Evidence Types (to be implemented in `isoforms/`)

```python
@attrs.define
class SpliceJunction:
    """A splice junction with evidence."""
    region: GenomicRegion
    donor: int
    acceptor: int
    read_count: int
    is_canonical: bool

@attrs.define
class CoverageProfile:
    """Per-base coverage for a region."""
    region: GenomicRegion
    values: np.ndarray
    strand_specific: bool = False
```

## Testing Requirements

### Test Organization
- Unit tests in `tests/test_<module>.py`
- Integration tests in `tests/integration/`
- Test data in `tests/data/`

### Coverage Targets
- Minimum 80% coverage for core modules
- 100% coverage for data structure classes
- All public functions must have at least one test

### Test Data Guidelines
- Use minimal synthetic examples for unit tests
- Include small real-world examples for integration tests
- Never commit large (>1MB) test files; use fixtures or downloads

### Running Tests
```bash
# All tests
pytest

# With coverage
pytest --cov=helixforge --cov-report=html

# Specific module
pytest tests/test_confidence.py -v
```

## Dependencies and Version Constraints

### Core Dependencies (Required)
| Package | Min Version | Purpose |
|---------|-------------|---------|
| h5py | 3.8 | HDF5 file I/O |
| pysam | 0.21 | BAM file parsing |
| pyfaidx | 0.7 | FASTA indexing |
| gffutils | 0.12 | GFF3/GTF parsing |
| numpy | 1.24 | Numerical operations |
| pandas | 2.0 | Data manipulation |
| scipy | 1.10 | Statistical functions |
| matplotlib | 3.7 | Static plotting |
| plotly | 5.14 | Interactive plots |
| jinja2 | 3.1 | Report templating |
| click | 8.1 | CLI framework |
| rich | 13.0 | Terminal output |
| attrs | 23.0 | Data classes |

### Python Version
- Minimum: Python 3.10 (for modern type hints, match statements)
- Tested: Python 3.10, 3.11, 3.12

## Next Development Targets

### Immediate (Phase 4: Advanced Refinement)
1. Implement gene merge detection for fragmented predictions
2. Implement gene split detection for fused predictions
3. Add merge_split.py module to core/
4. Integrate merge/split analysis into CLI

### Short-term (Isoform Reconstruction)
1. Build alternative transcript models from junction combinations
2. Implement transcript scoring based on evidence support
3. Select representative transcripts per gene
4. Add isoforms/ module structure

### Medium-term (UTR Extension)
1. Extend 5' UTR using RNA-seq coverage
2. Extend 3' UTR using coverage and poly-A signals
3. Validate UTR boundaries against known patterns

## Commands Reference

### Development Setup

```bash
# Clone and install
git clone https://github.com/aseetharam/helixforge.git
cd helixforge
pip install -e ".[dev]"

# Run pre-commit hooks
pre-commit install
pre-commit run --all-files

# Run tests
pytest -v

# Type checking
mypy src/helixforge
```

### CLI Commands

```bash
# Main refinement pipeline (recommended)
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r refine_report.tsv

# Multi-tissue refinement
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam liver.bam,brain.bam,heart.bam \
    --min-tissues 2 \
    --min-reads 3 \
    -o refined.gff3 \
    -r refine_report.tsv

# Standalone confidence scoring (no RNA-seq required)
helixforge confidence \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    -o confidence.tsv

# Standalone evidence scoring (no HDF5 required)
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3 \
    --report evidence_report.tsv

# Protein extraction and homology search
helixforge homology extract-proteins \
    --gff refined.gff3 \
    --genome genome.fa \
    -o proteins.fa

helixforge homology search \
    --query proteins.fa \
    --database uniprot.dmnd \
    -o hits.tsv \
    --threads 16

helixforge homology validate \
    --search-results hits.tsv \
    --gff refined.gff3 \
    -o validation.tsv

# QC aggregation and reporting
# Calculates combined AED from RNA-seq + homology + confidence evidence
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv validation.tsv \
    -o qc_results.tsv

# Custom AED weights (for species with limited RNA-seq data)
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv validation.tsv \
    -o qc_results.tsv \
    --aed-rnaseq-weight 0.3 \
    --aed-homology-weight 0.5 \
    --aed-confidence-weight 0.2

helixforge qc report \
    --qc-tsv qc_results.tsv \
    -o qc_report.html

helixforge qc tiered-output \
    --qc-results qc_results.tsv \
    --gff refined.gff3 \
    -o tiered/

# Parallel processing for large genomes
helixforge parallel plan \
    --genome genome.fa \
    --strategy adaptive \
    -o chunks.json

helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command 'helixforge refine -p pred.h5 -g pred.gff3 --genome genome.fa --rnaseq-bam rna.bam --region {seqid}:{start}-{end} --chunk-id {chunk_id} -o outputs/{chunk_id}.gff3' \
    -o tasks.txt

helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*.gff3' \
    -o combined.gff3 \
    --type merge_gff
```
