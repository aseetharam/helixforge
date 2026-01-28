# RNA-seq Evidence Scoring Guide

> **Note**: For most workflows, use the `helixforge refine` command which combines evidence scoring with splice correction, boundary adjustment, and confidence scoring. The standalone `evidence` command is useful when you only have RNA-seq data without the Helixer HDF5 predictions file.

This document explains the RNA-seq evidence scoring functionality in HelixForge, which evaluates how well predicted gene models match empirical RNA-seq data.

## Overview

The evidence scoring module (`helixforge.core.evidence`) evaluates gene predictions by:

1. **Junction Support**: Matching predicted splice sites to observed RNA-seq junctions
2. **Exon Coverage**: Measuring expression levels across predicted exons
3. **Boundary Agreement**: Checking coverage at start/stop codon positions
4. **Strand Consistency**: Verifying expression matches predicted strand

The output includes an Annotation Edit Distance (AED) score for each gene, where 0 = perfect evidence support and 1 = no support.

## Quick Start

```bash
# Basic usage
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3

# With detailed reports
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3 \
    --report evidence_report.tsv \
    --summary evidence_summary.tsv
```

## How Evidence Scoring Works

### Step 1: Junction Extraction

Splice junctions are extracted from RNA-seq BAM files using reads with CIGAR 'N' operations (introns). Each junction is characterized by:

- **Position**: Intron start and end coordinates
- **Support**: Number of reads spanning the junction
- **Strand**: Inferred from XS tag or read orientation

### Step 2: Junction Matching

For each predicted intron in a gene model:

1. Look for an exact match in extracted junctions
2. If not found, look for near-matches within `--boundary-tolerance` bp
3. Record read count and whether match is exact or shifted

A junction is considered "supported" if it has at least `--min-reads` reads.

### Step 3: Exon Coverage Analysis

For each predicted exon (unless `--no-coverage` is used):

1. Calculate per-base read coverage
2. Compute mean, median, min, and max coverage
3. Determine fraction of bases with coverage >= `--min-coverage`

An exon is considered "expressed" if median coverage >= `--min-coverage`.

### Step 4: Boundary Support

Start and stop codon positions are checked for coverage support:

- **Start supported**: First exon (for + strand) or last exon (for - strand) has minimum coverage
- **Stop supported**: Last exon (for + strand) or first exon (for - strand) has minimum coverage

### Step 5: AED Calculation

The Annotation Edit Distance (AED) is a weighted combination of support metrics:

```
AED = junction_weight * (1 - junction_ratio)
    + coverage_weight * (1 - coverage_ratio)
    + boundary_weight * (1 - boundary_ratio)
```

Default weights: junction=0.5, coverage=0.3, boundary=0.2

### Step 6: Evidence Level Classification

Genes are classified into four evidence levels:

| Level | Multi-exon Criteria | Single-exon Criteria |
|-------|---------------------|---------------------|
| FULL | 100% junction support AND ≥80% exon coverage | ≥80% exon coverage |
| PARTIAL | ≥80% junctions OR (≥50% junctions AND ≥50% coverage) | ≥50% exon coverage |
| MINIMAL | Any junction OR coverage support | Any coverage support |
| NONE | No evidence support | No coverage support |

## Using the CLI

### Basic Usage

```bash
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3
```

### Multiple BAM Files

```bash
# Comma-separated
helixforge evidence \
    -g predictions.gff3 \
    -b tissue1.bam,tissue2.bam,tissue3.bam \
    -o annotated.gff3

# Repeated flags
helixforge evidence \
    -g predictions.gff3 \
    -b tissue1.bam \
    -b tissue2.bam \
    -b tissue3.bam \
    -o annotated.gff3

# From file list
helixforge evidence \
    -g predictions.gff3 \
    --bam-list bam_files.txt \
    -o annotated.gff3
```

### With All Reports

```bash
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3 \
    --report evidence_report.tsv \
    --summary evidence_summary.tsv \
    --junction-details junction_details.tsv \
    --exon-details exon_details.tsv
```

### Junction-Only Mode

For faster analysis when only junction support matters:

```bash
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3 \
    --no-coverage
```

### Full Options

| Option | Default | Description |
|--------|---------|-------------|
| `-g/--gff` | Required | Input GFF3 file |
| `-b/--bam` | Required | BAM file(s), comma-separated or repeated |
| `--bam-list` | - | File with BAM paths (one per line) |
| `-o/--output` | Required | Output GFF3 with evidence attributes |
| `--report` | - | Detailed evidence report TSV |
| `--summary` | - | Summary TSV for QC aggregation |
| `--junction-details` | - | Per-junction details TSV |
| `--exon-details` | - | Per-exon details TSV |
| `--min-reads` | 3 | Minimum reads for junction support |
| `--min-coverage` | 5 | Minimum coverage for expression |
| `--boundary-tolerance` | 10 | Max bp shift for near-match junctions |
| `--no-coverage` | false | Skip exon coverage (junction-only) |

## Output Files

### Annotated GFF3

The output GFF3 includes evidence attributes on gene features:

```
chr1    HelixForge    gene    1000    5000    .    +    .    ID=gene1;evidence_level=full;aed=0.0500;junction_support=1.0000;exon_coverage=0.9500
```

| Attribute | Description |
|-----------|-------------|
| `evidence_level` | full/partial/minimal/none |
| `aed` | Annotation Edit Distance (0-1) |
| `junction_support` | Fraction of supported junctions |
| `exon_coverage` | Fraction of expressed exons |
| `evidence_flags` | Warning flags (if any) |

### Evidence Report TSV

Detailed per-gene metrics (columns):

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| transcript_id | Primary transcript ID |
| seqid | Scaffold/chromosome |
| strand | Strand (+ or -) |
| n_introns | Total introns |
| n_introns_supported | Introns with junction support |
| n_introns_exact | Introns with exact coordinate match |
| junction_support_ratio | Fraction supported (0-1) |
| n_exons | Total exons |
| n_exons_expressed | Exons with expression |
| mean_exon_coverage | Mean coverage across exons |
| exon_coverage_ratio | Fraction expressed (0-1) |
| start_supported | yes/no |
| stop_supported | yes/no |
| strand_consistent | yes/no |
| evidence_level | full/partial/minimal/none |
| aed | Annotation Edit Distance |
| flags | Warning flags |

### Junction Details TSV

Per-junction information:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| transcript_id | Transcript ID |
| intron_start | Intron start (0-based) |
| intron_end | Intron end (0-based, exclusive) |
| intron_length | Intron length in bp |
| observed | yes/no - has RNA-seq support |
| read_count | Supporting read count |
| unique_count | Uniquely mapped reads |
| exact_match | yes/no - exact coordinate match |
| shift_donor | Shift in donor site (bp) |
| shift_acceptor | Shift in acceptor site (bp) |

### Exon Details TSV

Per-exon coverage information:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| transcript_id | Transcript ID |
| exon_start | Exon start (0-based) |
| exon_end | Exon end (0-based, exclusive) |
| exon_length | Exon length in bp |
| mean_coverage | Mean read coverage |
| median_coverage | Median read coverage |
| min_coverage | Minimum coverage |
| max_coverage | Maximum coverage |
| fraction_covered | Fraction ≥ threshold |
| coverage_uniformity | Coefficient of variation |
| is_expressed | yes/no |

## Evidence Flags

| Flag | Description |
|------|-------------|
| NO_JUNCTION_SUPPORT | No introns have junction support |
| PARTIAL_JUNCTION_SUPPORT | Some but not all introns supported |
| NO_EXON_EXPRESSION | No exons show expression |
| LOW_EXON_COVERAGE | Less than half of exons expressed |
| START_UNSUPPORTED | No coverage at start codon |
| STOP_UNSUPPORTED | No coverage at stop codon |
| STRAND_INCONSISTENT | Expression pattern doesn't match strand |

## Recommended Parameters

### High-depth RNA-seq (>50M reads)

```bash
--min-reads 5 --min-coverage 10
```

Higher thresholds reduce false positives from background noise.

### Medium-depth RNA-seq (10-50M reads)

```bash
--min-reads 3 --min-coverage 5
```

Default settings work well for typical experiments.

### Low-depth RNA-seq (<10M reads)

```bash
--min-reads 2 --min-coverage 3
```

Lower thresholds to capture real but low-abundance genes.

### Conservative Analysis

```bash
--min-reads 10 --min-coverage 20 --boundary-tolerance 5
```

Only report high-confidence evidence.

## Integration with QC Pipeline

Evidence scores integrate with the HelixForge QC system. The recommended approach is to use the `refine` command which produces a comprehensive report:

```bash
# Recommended: Use the refine command (includes evidence scoring)
helixforge refine \
    -p helixer_predictions.h5 \
    -g predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r refine_report.tsv

# Aggregate with homology results
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv validation.tsv \
    -o qc_aggregated.tsv

# Generate report
helixforge qc report \
    --qc-tsv qc_aggregated.tsv \
    -o qc_report.html
```

### RNA-seq AED vs Combined AED

**Important:** The AED calculated by `refine` and `evidence` commands is based **only on RNA-seq evidence**. This can be misleading for genes with good homology support but no expression in sampled tissues.

The `qc aggregate` command calculates a **combined AED** that incorporates:
- RNA-seq evidence (junction support + coverage)
- Homology evidence (query coverage × identity)
- Helixer confidence (prediction probability)

This ensures genes are evaluated using all available evidence:

| Scenario | RNA-seq AED | Combined AED |
|----------|-------------|--------------|
| Good RNA-seq, good homology | 0.1 | ~0.10 |
| No RNA-seq, complete homology | 0.5 | ~0.26 |
| Good RNA-seq, no homology | 0.2 | ~0.28 |
| Neither | 0.5 | ~0.58 |

The combined AED weights are configurable via CLI options:
- `--aed-rnaseq-weight` (default: 0.4)
- `--aed-homology-weight` (default: 0.4)
- `--aed-confidence-weight` (default: 0.2)

If you only have RNA-seq data (no HDF5 predictions), use the standalone `evidence` command:

```bash
# Standalone evidence scoring (no HDF5 required)
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3 \
    --report evidence_report.tsv

# Aggregate with other results
helixforge qc aggregate \
    --evidence-tsv evidence_report.tsv \
    --homology-tsv validation.tsv \
    -o qc_aggregated.tsv
```

## API Reference

### EvidenceScorer Class

```python
from helixforge.core.evidence import EvidenceScorer, EvidenceScorerConfig
from helixforge.io.bam import JunctionExtractor
from helixforge.io.gff import GFF3Parser

# Configure scorer
config = EvidenceScorerConfig(
    min_junction_reads=3,
    min_exon_coverage=5,
    boundary_tolerance=10,
)

scorer = EvidenceScorer(config)

# Load data
genes = list(GFF3Parser("predictions.gff3").iter_genes())

with JunctionExtractor("rnaseq.bam") as extractor:
    junctions = extractor.extract_all()

    # Score all genes
    for score in scorer.score_genes(genes, junctions, extractor):
        print(f"{score.gene_id}: AED={score.aed:.3f}")
```

### EvidenceScore

```python
@attrs.define
class EvidenceScore:
    gene_id: str
    transcript_id: str
    seqid: str
    strand: str

    # Junction metrics
    n_introns: int
    n_introns_supported: int
    junction_support_ratio: float

    # Exon metrics
    n_exons: int
    n_exons_expressed: int
    exon_coverage_ratio: float

    # Boundary metrics
    start_supported: bool
    stop_supported: bool

    # Overall
    evidence_level: EvidenceLevel
    aed: float
    flags: list[str]

    # Detailed evidence
    junction_evidence: list[JunctionEvidence]
    exon_evidence: list[ExonEvidence]
```

### EvidenceLevel Enum

```python
class EvidenceLevel(Enum):
    FULL = "full"       # Fully supported
    PARTIAL = "partial" # Partially supported
    MINIMAL = "minimal" # Weak support
    NONE = "none"       # No support
```

## Troubleshooting

### Low Junction Support

If many genes show no junction support:

1. Verify BAM file uses same reference as GFF3
2. Check BAM is properly sorted and indexed
3. Try lowering `--min-reads` threshold
4. Check if RNA-seq is from the same species/tissue

### No Exon Coverage

If exon coverage is unexpectedly low:

1. Verify BAM alignments cover annotated regions
2. Check strand-specificity of RNA-seq library
3. Try lowering `--min-coverage` threshold
4. Consider using `--no-coverage` for junction-only analysis

### Many Near-Match Junctions

If most junctions are near-matches rather than exact:

1. Coordinate systems may be misaligned (check 0-based vs 1-based)
2. Try increasing `--boundary-tolerance`
3. Original predictions may need splice refinement first

### Memory Issues

For large genomes:

1. Process by chromosome using the splice command's `--region` option
2. Use `--no-coverage` to reduce memory footprint
3. Consider pre-extracting junctions to BED format
