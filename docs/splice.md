# Splice Site Refinement Guide

This document explains the splice site refinement functionality in HelixForge, which uses RNA-seq junction evidence and Position Weight Matrix (PWM) scoring to refine Helixer gene predictions.

## Overview

The splice refinement module (`helixforge.core.splice`) corrects predicted splice sites by:

1. **RNA-seq Junction Matching**: Aligning predicted introns to splice junctions observed in RNA-seq data
2. **PWM Scoring**: Evaluating splice site sequences against Position Weight Matrices trained on canonical splice motifs
3. **Canonical Site Verification**: Ensuring splice sites conform to GT-AG, GC-AG, or AT-AC dinucleotide patterns

## How Splice Refinement Works

### Step 1: Intron Extraction

For each gene, the refiner extracts introns from the transcript models. An intron is defined as the genomic region between consecutive exons.

```
Exon 1          Exon 2          Exon 3
[======]--------[======]--------[======]
        Intron 1        Intron 2
        GT...AG         GT...AG
```

### Step 2: Junction Matching

Each intron is compared against RNA-seq splice junctions within a configurable search window (`max_shift`, default 15bp). Junctions are scored by:

- **Exact match**: Intron boundaries match junction exactly
- **Near match**: Intron can be shifted <=max_shift bp to match a junction
- **No match**: No supporting junction found

### Step 3: PWM Scoring

If no exact junction match is found, the refiner uses PWM scoring to identify optimal splice sites. The PWM evaluates:

**Donor site (5' splice site)**:
- 9-position window centered on the GT dinucleotide
- Positions -3 to +6 relative to the splice site

**Acceptor site (3' splice site)**:
- 16-position window ending at the AG dinucleotide
- Positions -14 to +2 relative to the splice site

### Step 4: Site Correction

Splice sites are corrected when:
1. An RNA-seq junction provides strong evidence (>= `min_junction_reads`)
2. PWM scoring identifies a better canonical site within the search window
3. The correction maintains reading frame consistency

## Using the CLI

### Basic Usage

```bash
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    --output-gff refined.gff3 \
    --report splice_report.tsv
```

### With Pre-computed Junctions

If you have junctions from multiple samples, provide them as a BED file:

```bash
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --junctions-bed all_junctions.bed \
    --output-gff refined.gff3
```

### Full Options

```bash
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    --output-gff refined.gff3 \
    --report splice_report.tsv \
    --max-shift 15 \          # Maximum bp to shift splice site
    --min-reads 3 \           # Minimum junction read support
    --adjust-boundaries \     # Also refine start/stop codons
    --workers 4               # Parallel processing threads
```

## PWM Scoring Interpretation

### Score Ranges

PWM scores are log-likelihood ratios comparing the observed sequence to a background model:

| Score Range | Interpretation |
|-------------|----------------|
| > 5.0 | Strong canonical signal |
| 2.0 - 5.0 | Good canonical signal |
| 0.0 - 2.0 | Weak signal |
| < 0.0 | Non-canonical or poor match |

### Comparing Sites

When evaluating alternative splice positions, higher PWM scores indicate better matches to canonical splice motifs. The refiner selects positions that:

1. Have RNA-seq junction support (preferred)
2. Maximize combined donor + acceptor PWM scores
3. Maintain canonical GT-AG, GC-AG, or AT-AC dinucleotides

## Splice Site Types

HelixForge recognizes three splice site types:

### GT-AG (Major Spliceosome)
- ~99% of introns in most eukaryotes
- Donor: GT (GU in RNA)
- Acceptor: AG
- Processed by U2-type spliceosome

### GC-AG (Minor Canonical)
- ~0.5-1% of introns
- Donor: GC
- Acceptor: AG
- Processed by U2-type spliceosome

### AT-AC (U12 Minor Spliceosome)
- ~0.1-0.5% of introns
- Donor: AT
- Acceptor: AC
- Processed by U12-type spliceosome
- Often found in conserved genes with regulatory functions

## Recommended Parameters

### High-depth RNA-seq (>50M reads)

```bash
--max-shift 15 --min-reads 5
```

With high coverage, require more reads for confidence. Junction detection is reliable.

### Medium-depth RNA-seq (10-50M reads)

```bash
--max-shift 15 --min-reads 3
```

Default settings work well for most applications.

### Low-depth RNA-seq (<10M reads)

```bash
--max-shift 10 --min-reads 2
```

Reduce shift distance to avoid false corrections. Accept lower read thresholds.

### No RNA-seq Data

```bash
--max-shift 10 --min-reads 0
```

Rely entirely on PWM scoring. Only correct to canonical sites with strong motif scores.

### Conservative Refinement

```bash
--max-shift 5 --min-reads 10
```

Make minimal changes. Only correct obvious errors with strong evidence.

## Output Files

### Refined GFF3

The refined GFF3 file contains corrected gene models. Modified features include:

- Corrected exon boundaries
- Updated CDS coordinates (maintaining phase)
- Attributes documenting corrections made

### Splice Report TSV

The report file contains per-gene statistics:

| Column | Description |
|--------|-------------|
| gene_id | Gene identifier |
| n_introns | Total introns in gene |
| n_supported | Introns with junction support |
| n_corrected | Introns that were corrected |
| n_canonical | Introns with canonical splice sites |
| n_noncanonical | Introns with non-canonical sites |
| support_ratio | Fraction of supported introns |
| flags | Warning flags (e.g., "noncanonical:2") |

## Visualization

### Locus-level Plot

```python
from helixforge.viz.locus import plot_splice_refinement

fig = plot_splice_refinement(
    original_gene=gene,
    refined_gene=refined_gene,
    report=report,
    junctions=junctions,
    genome=genome,
)
fig.savefig("gene_refinement.png")
```

### Genome-wide Summary

```python
from helixforge.viz.genome import plot_splice_summary

fig = plot_splice_summary(reports)
fig.savefig("splice_summary.png")
```

## Troubleshooting

### Many Non-canonical Sites

If the report shows many non-canonical splice sites:

1. Check genome/annotation coordinate systems match
2. Verify BAM file uses same reference
3. Consider if species has unusual splice patterns (e.g., some fungi)

### Low Junction Support

If few introns have junction support:

1. Verify RNA-seq is from same species/tissue
2. Check BAM file indexing
3. Consider lowering `--min-reads` threshold

### Over-correction

If too many corrections are being made:

1. Increase `--min-reads` threshold
2. Decrease `--max-shift` distance
3. Review original Helixer predictions for systematic issues

## API Reference

### SpliceRefiner Class

```python
from helixforge.core.splice import SpliceRefiner

refiner = SpliceRefiner(
    genome=genome_accessor,      # GenomeAccessor instance
    junctions=junction_dict,     # Dict mapping seqid -> list[SpliceJunction]
    donor_pwm=None,              # Optional custom donor PWM
    acceptor_pwm=None,           # Optional custom acceptor PWM
    max_shift=15,                # Maximum correction distance
    min_junction_reads=3,        # Minimum reads for junction support
    pwm_threshold=0.0,           # Minimum PWM score for correction
)

# Refine single gene
refined_gene, report = refiner.refine_gene(gene)

# Parallel processing
results = refiner.refine_genes_parallel(genes, n_workers=4)
for refined_gene, report in results:
    ...
```

### GeneSpliceReport

```python
@attrs.define
class GeneSpliceReport:
    gene_id: str
    n_introns: int
    n_supported: int
    n_corrected: int
    n_canonical: int
    n_noncanonical: int
    corrections: list[SpliceCorrection] = []
    unsupported_introns: list[int] = []
    flags: list[str] = []

    @property
    def support_ratio(self) -> float:
        """Fraction of introns with RNA-seq support."""

    @property
    def fully_supported(self) -> bool:
        """True if all introns have support."""
```

### PositionWeightMatrix

```python
from helixforge.core.splice import PositionWeightMatrix

# Load plant defaults
donor_pwm, acceptor_pwm = PositionWeightMatrix.load_plant_defaults()

# Score a sequence
score = donor_pwm.score("CAGGTAAGT")

# Score all positions in a longer sequence
scores = donor_pwm.score_all_positions(long_sequence)
best_position = scores.argmax()

# Build custom PWM from training sequences
custom_pwm = PositionWeightMatrix.from_sequences(
    sequences=["CAGGTAAGT", "AAGGTAAGT", ...],
    site_type="donor",
)
```
