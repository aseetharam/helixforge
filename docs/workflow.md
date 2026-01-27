# HelixForge Workflow Guide

This guide shows the recommended order of operations for refining Helixer gene predictions into publication-quality annotations.

## Quick Reference

```
                          ┌─────────────────────┐
                          │   Helixer Output    │
                          │  predictions.h5     │
                          │  predictions.gff3   │
                          └──────────┬──────────┘
                                     │
                                     │  + RNA-seq BAM(s)
                                     ▼
                          ┌─────────────────────┐
                          │       refine        │
                          │  (main pipeline)    │
                          │                     │
                          │ • Splice correction │
                          │ • Boundary adjust   │
                          │ • Confidence score  │
                          │ • Evidence score    │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │      validate       │
                          │ (protein homology)  │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │    qc aggregate     │
                          │    qc report        │
                          └──────────┬──────────┘
                                     │
                                     ▼
                          ┌─────────────────────┐
                          │   qc tiered-output  │
                          │   (final GFFs)      │
                          └─────────────────────┘
```

## Step-by-Step Workflow

### Prerequisites

You need the following files from Helixer:
- `*_predictions.h5` - HDF5 file with softmax probabilities
- `*_input.h5` - HDF5 file with coordinate mapping (usually auto-detected)
- `*.gff3` - Gene predictions in GFF3 format

Optional but recommended:
- `genome.fa` - Reference genome FASTA (indexed with `samtools faidx`)
- `rnaseq.bam` - RNA-seq alignments (indexed with `samtools index`)
- `proteins.fa` - Reference protein database (UniProt, SwissProt, etc.)

---

### Step 1: Refine Gene Predictions (Main Pipeline)

The `refine` command is the primary HelixForge pipeline. It combines splice correction, boundary adjustment, confidence scoring, and evidence scoring into a single integrated workflow.

**RNA-seq evidence is required** for this step.

```bash
# Basic refinement (single BAM)
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam rnaseq.bam \
    -o refined.gff3 \
    -r refine_report.tsv
```

```bash
# Multi-tissue refinement (higher confidence)
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam liver.bam,brain.bam,heart.bam \
    --min-tissues 2 \
    --min-reads 3 \
    -o refined.gff3 \
    -r refine_report.tsv \
    --splice-details splice_corrections.tsv

# Or from a file list
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    --rnaseq-bam-list tissue_bams.txt \
    --min-tissues 2 \
    -o refined.gff3 \
    -r refine_report.tsv
```

**Key options:**
- `-p/--helixer-h5` - Helixer HDF5 predictions (required)
- `-g/--helixer-gff` - Helixer GFF3 predictions (required)
- `--rnaseq-bam` - BAM file(s), comma-separated or repeated
- `--rnaseq-bam-list` - File containing BAM paths
- `--junctions-bed` - Pre-computed junction files (STAR SJ.out.tab)
- `--max-shift` - Maximum bp to shift splice sites (default: 15)
- `--min-reads` - Minimum reads for junction support (default: 3)
- `--min-tissues` - Minimum samples supporting a junction (default: 1)
- `--adjust-boundaries/--no-adjust-boundaries` - Refine start/stop codons (default: on)
- `-r/--report` - Output comprehensive refine report TSV

**Alternative: Using STAR junction files:**

```bash
helixforge refine \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    --genome genome.fa \
    --junctions-bed sample1_SJ.out.tab,sample2_SJ.out.tab \
    --min-tissues 2 \
    -o refined.gff3
```

**Output:**
- Refined GFF3 with corrected splice sites and quality attributes
- Report TSV with confidence, evidence, and splice correction metrics

---

### Alternative: Standalone Commands

If you don't have RNA-seq data or prefer separate steps, you can use:

#### Confidence Scoring Only

```bash
helixforge confidence \
    -p helixer_predictions.h5 \
    -g helixer_predictions.gff3 \
    -o confidence_scores.tsv \
    --bed confidence.bed
```

#### Evidence Scoring Only (no HDF5 required)

```bash
helixforge evidence \
    -g predictions.gff3 \
    -b rnaseq.bam \
    -o annotated.gff3 \
    --report evidence_report.tsv
```

---

### Step 2: Validate with Protein Homology

Search predicted proteins against a reference database to validate gene models.

```bash
# Extract proteins from refined predictions
helixforge homology extract-proteins \
    --gff refined.gff3 \
    --genome genome.fa \
    -o predicted_proteins.fa \
    --longest-isoform

# Search against database
helixforge homology search \
    --proteins predicted_proteins.fa \
    --database uniprot_plants.dmnd \
    -o homology_hits.tsv \
    --threads 16 \
    --evalue 1e-5

# Validate gene models
helixforge homology validate \
    --search-results homology_hits.tsv \
    --gff refined.gff3 \
    -o homology_validation.tsv \
    --te-regions repeats.bed
```

**Output:** TSV with homology validation metrics per gene.

---

### Step 3: Aggregate QC Results

Combine results from refine and homology into unified QC metrics.

```bash
# Using refine report (recommended)
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv homology_validation.tsv \
    -o qc_aggregated.tsv
```

**Output:** Combined TSV with all QC metrics and tier classifications per gene.

---

### Step 4: Generate QC Report

Create a comprehensive HTML report summarizing annotation quality.

```bash
helixforge qc report \
    --qc-tsv qc_aggregated.tsv \
    -o qc_report.html \
    --title "Genome Annotation QC"
```

**Output:** Interactive HTML report with statistics and visualizations.

---

### Step 5: Create Tiered Output

Generate filtered gene sets based on quality tiers.

```bash
helixforge qc tiered-output \
    --qc-results qc_aggregated.tsv \
    --gff refined.gff3 \
    -o tiered_output/
```

**Output directory contents:**

- `high_confidence.gff3` - Highest quality genes
- `medium_confidence.gff3` - Moderate quality genes
- `low_confidence.gff3` - Lower quality, needs review
- `rejected.gff3` - Problematic genes
- `tier_summary.tsv` - Summary statistics

---

### Step 6: Visualize Results (optional)

Inspect specific genes or regions interactively.

```bash
# Static plot for a gene
helixforge viz \
    -g refined.gff3 \
    --genome genome.fa \
    --bam rnaseq.bam \
    --gene GENE_ID \
    -o gene_plot.png

# Interactive browser
helixforge viz \
    -g refined.gff3 \
    --genome genome.fa \
    --bam rnaseq.bam \
    --interactive \
    --port 8050
```

---

## Parallel Processing for Large Genomes

For genomes larger than ~500 Mb, use parallel processing:

```bash
# 1. Create chunk plan
helixforge parallel plan \
    --genome genome.fa \
    --target-size 50000000 \
    -o chunks.json

# 2. Generate task file for confidence scoring
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command 'helixforge confidence -p predictions.h5 -g predictions.gff3 --region {seqid}:{start}-{end} --chunk-id {chunk_id} -o outputs/{chunk_id}.tsv' \
    --output tasks.txt \
    --output-dir outputs/

# 3. Execute with HyperShell or GNU Parallel
hs cluster tasks.txt --num-tasks 32
# or
parallel -j 32 < tasks.txt

# 4. Aggregate chunk results
helixforge parallel aggregate \
    --chunk-plan chunks.json \
    --input-dir outputs/ \
    --pattern '*.tsv' \
    -o confidence_scores.tsv
```

---

## Complete Example Pipeline

Here's a complete pipeline script:

```bash
#!/bin/bash
set -euo pipefail

# Input files
PREDICTIONS_H5="helixer_output/predictions.h5"
PREDICTIONS_GFF="helixer_output/predictions.gff3"
GENOME="genome.fa"
RNASEQ_BAM="rnaseq.bam"
PROTEINS="uniprot_plants.fa"

# Output directory
OUTDIR="helixforge_output"
mkdir -p "$OUTDIR"

echo "Step 1: Confidence scoring..."
helixforge confidence \
    -p "$PREDICTIONS_H5" \
    -g "$PREDICTIONS_GFF" \
    -o "$OUTDIR/confidence.tsv" \
    --bed "$OUTDIR/confidence.bed" \
    -j 4

echo "Step 2: Splice refinement..."
# For single BAM:
helixforge splice \
    --helixer-gff "$PREDICTIONS_GFF" \
    --genome "$GENOME" \
    --rnaseq-bam "$RNASEQ_BAM" \
    --output-gff "$OUTDIR/splice_refined.gff3" \
    --report "$OUTDIR/splice_report.tsv" \
    --adjust-boundaries \
    --workers 4

# For multiple BAMs (multi-tissue), replace above with:
# helixforge splice \
#     --helixer-gff "$PREDICTIONS_GFF" \
#     --genome "$GENOME" \
#     --rnaseq-bam liver.bam,brain.bam,heart.bam \
#     --min-tissues 2 --min-reads 3 \
#     --output-gff "$OUTDIR/splice_refined.gff3" \
#     --report "$OUTDIR/splice_report.tsv" \
#     --adjust-boundaries --workers 4

echo "Step 3: Protein homology validation..."
helixforge validate \
    -g "$OUTDIR/splice_refined.gff3" \
    -p "$PROTEINS" \
    -o "$OUTDIR/validated.gff3" \
    --tool diamond

echo "Step 4: QC aggregation..."
helixforge qc aggregate \
    --gff "$OUTDIR/validated.gff3" \
    --confidence "$OUTDIR/confidence.tsv" \
    --splice "$OUTDIR/splice_report.tsv" \
    -o "$OUTDIR/qc_aggregated.tsv"

echo "Step 5: Generate report..."
helixforge qc report \
    --aggregated "$OUTDIR/qc_aggregated.tsv" \
    -o "$OUTDIR/qc_report.html"

echo "Step 6: Create tiered output..."
helixforge qc tiered-output \
    --aggregated "$OUTDIR/qc_aggregated.tsv" \
    --gff "$OUTDIR/validated.gff3" \
    --output-dir "$OUTDIR/tiered/"

echo "Done! Results in $OUTDIR/"
```

---

## Minimal Workflow (No RNA-seq)

If you don't have RNA-seq data:

```bash
# 1. Confidence scoring only
helixforge confidence \
    -p predictions.h5 \
    -g predictions.gff3 \
    -o confidence.tsv

# 2. Validate with proteins
helixforge validate \
    -g predictions.gff3 \
    -p proteins.fa \
    -o validated.gff3

# 3. QC and filter
helixforge qc aggregate \
    --gff validated.gff3 \
    --confidence confidence.tsv \
    -o qc_aggregated.tsv

helixforge qc tiered-output \
    --aggregated qc_aggregated.tsv \
    --gff validated.gff3 \
    --output-dir tiered/
```

---

## Command Summary

| Step | Command | Purpose |
|------|---------|---------|
| 1 | `confidence` | Score genes using Helixer probabilities |
| 2 | `splice` | Refine splice sites with RNA-seq |
| 3 | `add-evidence` | Annotate with RNA-seq support |
| 4 | `validate` | Protein homology validation |
| 5 | `qc aggregate` | Combine all QC metrics |
| 6 | `qc report` | Generate HTML report |
| 7 | `qc tiered-output` | Create quality-filtered GFFs |
| 8 | `viz` | Visualize genes and evidence |

**For large genomes:** Use `parallel plan`, `parallel tasks`, and `parallel aggregate` to distribute work.
