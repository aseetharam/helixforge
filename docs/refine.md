# Gene Refinement Pipeline Guide

This document explains the main `refine` command in HelixForge, which combines splice site correction, boundary adjustment, confidence scoring, and evidence scoring into a single integrated pipeline.

## Overview

The `refine` command is the primary HelixForge pipeline for transforming Helixer predictions into publication-quality annotations. It requires:

1. **Helixer HDF5 predictions** - For confidence scoring from neural network outputs
2. **RNA-seq evidence** - For splice refinement and evidence scoring

The pipeline performs these steps automatically:

1. **Splice Site Correction**: Adjust splice sites to match RNA-seq junctions
2. **Boundary Adjustment**: Refine start/stop codon positions
3. **Confidence Scoring**: Calculate confidence from Helixer HDF5 predictions
4. **Evidence Scoring**: Evaluate RNA-seq support (AED calculation)
5. **Flag Aggregation**: Combine quality flags from all steps

## Quick Start

```bash
# Basic refinement with single BAM
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3

# Multi-tissue with full reports
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam liver.bam,brain.bam,heart.bam \
    --min-tissues 2 \
    -o refined.gff3 \
    -r refine_report.tsv
```

## Pipeline Steps

### Step 1: Splice Site Correction

The pipeline first corrects splice sites using RNA-seq junction evidence:

1. Extract splice junctions from BAM files
2. For each predicted intron, search for matching observed junctions
3. If a nearby junction exists (within `--max-shift` bp), adjust the splice site
4. Use position weight matrix (PWM) scoring to prefer canonical GT-AG sites
5. Track corrections for reporting

**Key parameters:**
- `--max-shift`: Maximum distance to shift splice sites (default: 15 bp)
- `--min-reads`: Minimum reads supporting a junction (default: 3)
- `--min-tissues`: For multi-BAM, minimum samples supporting a junction

### Step 2: Boundary Adjustment

If `--adjust-boundaries` is enabled (default), the pipeline refines start and stop codons:

1. Search within `--boundary-window` bp for better codon positions
2. Evaluate candidates using Kozak context for start codons
3. Ensure proper ATG/stop codon sequences
4. Preserve reading frame

**Key parameters:**
- `--adjust-boundaries/--no-adjust-boundaries`: Enable boundary adjustment (default: enabled)
- `--boundary-window`: Search window size (default: 30 bp)

### Step 3: Confidence Scoring

Confidence scores are calculated from Helixer HDF5 predictions:

1. Extract per-base class probabilities for each gene region
2. Calculate mean confidence across exons, introns, and boundaries
3. Identify weak regions with low confidence
4. Flag genes below confidence thresholds

**Key parameters:**
- `--confidence-threshold`: Threshold for LOW_CONF flag (default: 0.5)

### Step 4: Evidence Scoring

RNA-seq evidence is evaluated for each refined gene:

1. Check junction support from RNA-seq alignments
2. Calculate exon coverage profiles
3. Evaluate boundary support at start/stop codons
4. Compute Annotation Edit Distance (AED)

The AED score is MAKER-compatible: 0 = perfect support, 1 = no support.

> **Note:** This AED is based on **RNA-seq evidence only**. For a combined AED that incorporates homology and Helixer confidence, use `helixforge qc aggregate` which calculates `combined_aed` from all available evidence types. This is important for genes with good homology support but no expression in sampled tissues.

**Key parameters:**
- `--min-coverage`: Minimum coverage for exon expression (default: 5)
- `--boundary-tolerance`: Maximum shift for near-match junctions (default: 10 bp)
- `--no-coverage`: Skip exon coverage analysis (junction-only mode)

## Command Line Options

### Required Inputs

| Option | Description |
|--------|-------------|
| `-p/--helixer-h5` | Helixer HDF5 predictions file |
| `-g/--helixer-gff` | Helixer GFF3 predictions file |
| `--genome` | Reference genome FASTA file |
| `-o/--output` | Output refined GFF3 file |

### RNA-seq Evidence (at least one required)

| Option | Description |
|--------|-------------|
| `--rnaseq-bam` | BAM file(s), comma-separated or repeated |
| `--rnaseq-bam-list` | File containing BAM paths (one per line) |
| `--junctions-bed` | Pre-extracted junction files (BED or STAR SJ.out.tab) |
| `--junctions-list` | File containing junction file paths |

### Output Reports

| Option | Description |
|--------|-------------|
| `-r/--report` | Refine report TSV with all scores |
| `--splice-details` | Detailed splice corrections TSV |
| `--evidence-details` | Per-gene evidence details TSV |
| `--unsupported-bed` | BED of introns without RNA-seq support |

### Refinement Parameters

| Option | Default | Description |
|--------|---------|-------------|
| `--max-shift` | 15 | Maximum splice site correction distance |
| `--min-reads` | 3 | Minimum reads for junction support |
| `--min-tissues` | 1 | Minimum samples supporting a junction |
| `--adjust-boundaries/--no-adjust-boundaries` | enabled | Enable boundary adjustment |
| `--boundary-window` | 30 | Search window for boundaries |
| `--min-coverage` | 5 | Minimum coverage for exon expression |
| `--boundary-tolerance` | 10 | Tolerance for near-match junctions |
| `--confidence-threshold` | 0.5 | Threshold for LOW_CONF flag |
| `--no-coverage` | false | Skip exon coverage analysis |

### Parallel Processing

| Option | Description |
|--------|-------------|
| `-j/--workers` | Number of parallel workers |
| `--region` | Process specific region (seqid:start-end) |
| `--scaffold` | Process entire scaffold |
| `--chunk-id` | Chunk identifier for parallel logging |

## Output Files

### Refined GFF3

The output GFF3 includes these attributes on gene features:

```
chr1    HelixForge    gene    1000    5000    .    +    .    ID=gene1;confidence_score=0.9234;evidence_score=0.9500;aed=0.0500;junction_support=1.0000;mean_coverage=45.3;flags=SPLICE_CORRECTED
```

| Attribute | Description |
|-----------|-------------|
| `confidence_score` | HDF5-based confidence score (0-1) |
| `evidence_score` | RNA-seq evidence score (0-1), equals 1 - AED |
| `aed` | RNA-seq Annotation Edit Distance (0-1, MAKER-compatible) |
| `junction_support` | Fraction of supported junctions (e.g., "5/5" or ratio) |
| `mean_coverage` | Mean exon coverage from RNA-seq |
| `flags` | Comma-separated quality flags |

> **Note:** The `aed` attribute here is based only on RNA-seq evidence. The `qc aggregate` command calculates a `combined_aed` that also incorporates homology and Helixer confidence.
| `evidence_level` | Optional: full/partial/minimal/none |
| `splice_corrections` | Optional: number of corrected splice sites |
| `boundary_adjusted` | Optional: yes if boundaries were adjusted |

### Refine Report TSV

The `-r/--report` option generates a comprehensive TSV with columns:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| seqid | Scaffold/chromosome |
| start | Gene start position |
| end | Gene end position |
| strand | Strand (+ or -) |
| n_exons | Number of exons |
| confidence_score | HDF5 confidence (0-1) |
| evidence_score | Evidence score = 1 - AED |
| aed | Annotation Edit Distance |
| junction_support | Supported introns / total introns |
| mean_coverage | Mean exon coverage |
| splice_corrections | Number of splice corrections |
| boundary_adjusted | true/false |
| flags | Quality flags |

## Quality Flags

Genes may receive these quality flags:

### Confidence Flags
- `VERY_LOW_CONF` - Confidence < 0.5
- `LOW_CONF` - Confidence < 0.7
- `WEAK_EXON` - Worst exon score < 0.6

### Evidence Flags
- `NO_JUNCTION_SUPPORT` - No introns have RNA-seq support
- `PARTIAL_JUNCTION_SUPPORT` - Some introns unsupported
- `NO_EXON_EXPRESSION` - No exons show expression
- `LOW_EXON_COVERAGE` - Less than half of exons expressed

### Splice Flags
- `SPLICE_CORRECTED` - Splice sites were adjusted
- `NONCANONICAL_SPLICE` - Gene has non-GT-AG splice sites

### Boundary Flags
- `BOUNDARY_ADJUSTED` - Start/stop codons were adjusted

## Multi-Tissue Analysis

When using multiple BAM files (e.g., different tissues):

```bash
helixforge refine \
    -p helixer.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam liver.bam,brain.bam,heart.bam \
    --min-tissues 2 \
    --min-reads 3 \
    -o refined.gff3
```

With `--min-tissues 2`:
- Junctions require support from at least 2 different samples
- `--min-reads` is applied per-sample (each supporting sample needs ≥3 reads)
- Improves confidence by reducing tissue-specific artifacts

## Integration with QC Pipeline

The refine report integrates directly with the QC aggregation system:

```bash
# 1. Run refinement
helixforge refine \
    -p helixer.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r refine_report.tsv

# 2. Aggregate with homology results
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv validation.tsv \
    -o qc_results.tsv

# 3. Generate report
helixforge qc report \
    --qc-tsv qc_results.tsv \
    -o qc_report.html
```

## Recommended Workflows

### Standard Refinement

For typical RNA-seq with adequate depth (>20M reads):

```bash
helixforge refine \
    -p helixer.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r refine_report.tsv \
    --splice-details splice_details.tsv
```

### High-Stringency (Publication)

For high-confidence publication annotations:

```bash
helixforge refine \
    -p helixer.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam tissue1.bam,tissue2.bam,tissue3.bam \
    --min-tissues 2 \
    --min-reads 5 \
    --max-shift 10 \
    --min-coverage 10 \
    -o refined.gff3 \
    -r refine_report.tsv
```

### Low-Depth RNA-seq

For limited sequencing depth (<10M reads):

```bash
helixforge refine \
    -p helixer.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    --min-reads 2 \
    --min-coverage 3 \
    --max-shift 20 \
    -o refined.gff3
```

### Junction-Only (Fast Mode)

When exon coverage analysis is not needed:

```bash
helixforge refine \
    -p helixer.h5 \
    -g helixer.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    --no-coverage \
    -o refined.gff3
```

## Troubleshooting

### Low Splice Correction Rate

If few genes show splice corrections:

1. Verify BAM file uses same reference as GFF3
2. Check BAM is properly sorted and indexed
3. Try increasing `--max-shift` to 20-25
4. Lower `--min-reads` threshold

### Many Unsupported Junctions

If many introns lack RNA-seq support:

1. Verify RNA-seq is from the same species
2. Check tissue/condition matches expectations
3. Use `--unsupported-bed` to examine problematic regions
4. Consider adding more tissue samples

### Memory Issues

For large genomes:

1. Process by scaffold using `--scaffold`
2. Use `--no-coverage` to reduce memory
3. Use parallel processing with `--region` and `--chunk-id`

## API Reference

### RefinePipeline Class

```python
from helixforge.core.refine import RefinePipeline, RefineConfig
from helixforge.io.fasta import GenomeAccessor
from helixforge.io.hdf5 import HelixerHDF5Reader

# Configure
config = RefineConfig(
    max_shift=15,
    min_junction_reads=3,
    adjust_boundaries=True,
)

# Initialize
pipeline = RefinePipeline(
    genome=GenomeAccessor("genome.fa"),
    hdf5_reader=HelixerHDF5Reader("helixer.h5"),
    bam_files=[Path("rnaseq.bam")],
    config=config,
)

# Refine genes
for result in pipeline.refine_genes(genes):
    print(f"{result.gene.gene_id}: conf={result.confidence_score:.3f}, AED={result.aed:.3f}")
```

### RefinedGene Result

```python
@attrs.define
class RefinedGene:
    gene: GeneModel           # Refined gene model
    original_gene: GeneModel  # Original before refinement

    # Scores
    confidence: GeneConfidence | None
    evidence: EvidenceScore | None
    splice_report: GeneSpliceReport | None

    # Details
    splice_corrections: int   # Number of splice corrections
    boundary_adjusted: bool   # Whether boundaries changed
    flags: list[str]          # Quality flags

    # Convenience properties
    confidence_score: float | None  # HDF5 confidence
    evidence_score: float | None    # 1 - AED
    aed: float | None               # Annotation Edit Distance
```
