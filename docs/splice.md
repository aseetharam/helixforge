# Splice Site Refinement Guide

> **DEPRECATED**: The standalone `splice` command has been deprecated. Use `helixforge refine` instead, which combines splice correction with boundary adjustment, confidence scoring, and evidence scoring in a single integrated pipeline.
>
> See [docs/refine.md](refine.md) and [docs/workflow.md](workflow.md) for the recommended workflow.

---

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

> **Note**: The examples below use the deprecated `splice` command. For new projects, use the `refine` command:
>
> ```bash
> # Recommended: Use the refine command
> helixforge refine \
>     -p predictions.h5 \
>     -g predictions.gff3 \
>     --genome genome.fa \
>     --rnaseq-bam rnaseq.bam \
>     -o refined.gff3 \
>     -r refine_report.tsv
> ```

### Basic Usage (DEPRECATED)

```bash
# DEPRECATED - use 'helixforge refine' instead
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    --output-gff refined.gff3 \
    --report splice_report.tsv
```

### With Multiple Tissues/Samples

When RNA-seq data from multiple tissues or developmental stages is available, provide BAM files to improve splice junction confidence and enable tissue-aware filtering.

**Multiple input formats are supported:**

```bash
# Comma-separated list
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam liver.bam,brain.bam,heart.bam,root.bam \
    --min-tissues 2 \
    --min-reads 3 \
    --output-gff refined.gff3

# Repeated flag (equivalent)
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam liver.bam \
    --rnaseq-bam brain.bam \
    --rnaseq-bam heart.bam \
    --min-tissues 2 \
    --output-gff refined.gff3

# From a file list (one path per line)
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam-list bam_files.txt \
    --min-tissues 2 \
    --output-gff refined.gff3

# Mixed (all can be combined)
helixforge splice \
    --rnaseq-bam liver.bam,brain.bam \
    --rnaseq-bam heart.bam \
    --rnaseq-bam-list additional_bams.txt \
    ...
```

**Benefits of multi-BAM mode:**
- **Tissue-specific junction evidence**: Track which tissues support each splice junction
- **Weighted confidence**: Junctions supported by multiple tissues are more reliable
- **Avoids dilution**: Rare tissue-specific isoforms aren't drowned out by high-expression tissues
- **Flexible thresholds**: Apply `--min-reads` per-sample when `--min-tissues` > 1

**How thresholds work in multi-BAM mode:**
- When `--min-tissues 1` (default): Junctions are aggregated and `--min-reads` applies to the total count across all samples
- When `--min-tissues > 1`: `--min-reads` applies per-sample; a junction must have >= `--min-reads` in at least `--min-tissues` samples

### With Pre-computed Junctions

HelixForge supports junction files in BED format or STAR's `SJ.out.tab` format. The format is auto-detected.

```bash
# Single junction file (BED or STAR SJ.out.tab)
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --junctions-bed junctions.bed \
    --output-gff refined.gff3

# Multiple junction files (comma-separated)
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --junctions-bed liver_SJ.out.tab,brain_SJ.out.tab,heart_SJ.out.tab \
    --min-tissues 2 \
    --output-gff refined.gff3

# From a file list
helixforge splice \
    --helixer-gff predictions.gff3 \
    --genome genome.fa \
    --junctions-list junction_files.txt \
    --min-tissues 2 \
    --output-gff refined.gff3
```

### Generating Junction Files

**From STAR aligner:**

STAR automatically generates `*_SJ.out.tab` files during alignment. These can be used directly:

```bash
# STAR alignment (generates SJ.out.tab)
STAR --genomeDir genome_index \
    --readFilesIn reads_R1.fq reads_R2.fq \
    --outFileNamePrefix sample_ \
    --outSAMtype BAM SortedByCoordinate

# Use the junction file directly
helixforge splice \
    --junctions-bed sample_SJ.out.tab \
    ...
```

**STAR SJ.out.tab format** (tab-separated, 9 columns):
1. Chromosome
2. First base of intron (1-based)
3. Last base of intron (1-based)
4. Strand (0=undefined, 1=+, 2=-)
5. Intron motif (0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC, 6=GT/AT)
6. Annotation status (0=unannotated, 1=annotated)
7. Uniquely mapping reads
8. Multi-mapping reads
9. Maximum overhang

**From regtools:**

```bash
# Extract junctions from BAM
regtools junctions extract -s 0 -a 8 -m 50 -M 500000 \
    aligned.bam -o junctions.bed
```

**From portcullis:**

```bash
# Full portcullis pipeline
portcullis full --threads 8 genome.fa aligned.bam -o portcullis_out

# Use the filtered junctions
helixforge splice \
    --junctions-bed portcullis_out/portcullis.filtered.pass.junctions.bed \
    ...
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
    --min-tissues 1 \         # Minimum samples supporting junction
    --adjust-boundaries \     # Also refine start/stop codons
    --workers 4               # Parallel processing threads
```

### Option Reference

| Option | Default | Description |
|--------|---------|-------------|
| `--rnaseq-bam` | - | BAM file(s). Comma-separated or repeat flag |
| `--rnaseq-bam-list` | - | File with BAM paths (one per line) |
| `--junctions-bed` | - | Junction file(s) in BED or STAR SJ.out.tab format |
| `--junctions-list` | - | File with junction file paths (one per line) |
| `--max-shift` | 15 | Maximum bp to shift splice sites |
| `--min-reads` | 3 | Minimum reads supporting a junction |
| `--min-tissues` | 1 | Minimum samples supporting a junction |
| `--adjust-boundaries` | false | Also refine start/stop codons |
| `--workers` | 1 | Parallel processing threads |
| `--report` | - | Output splice report TSV |
| `--corrections-detail` | - | Output detailed corrections TSV |
| `--unsupported-bed` | - | Output BED of unsupported introns |

**Notes:**
- When `--min-tissues > 1`, the `--min-reads` threshold applies per-sample rather than to the total count
- Junction file format (BED vs STAR SJ.out.tab) is auto-detected
- You can mix BAM and junction file inputs; junctions will be aggregated from all sources

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

### Multi-tissue Analysis

When you have RNA-seq from multiple tissues/conditions:

```bash
# High confidence: require support from at least 2 tissues
--min-tissues 2 --min-reads 3

# Very high confidence: require support from 3+ tissues
--min-tissues 3 --min-reads 2

# Relaxed: aggregate all samples (default behavior)
--min-tissues 1 --min-reads 3
```

**Rationale**: Junctions observed across multiple independent samples are more likely to be real biological splice sites rather than alignment artifacts. Tissue-specific isoforms (supported by only one tissue) are still retained when `--min-tissues 1` is used.

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
