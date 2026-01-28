# HelixForge Workflow Guide

This guide shows the recommended order of operations for refining Helixer gene predictions into publication-quality annotations.

## Available Commands

```
helixforge
├── refine        # Main pipeline (recommended) - combines all refinement steps
├── confidence    # Standalone confidence scoring from HDF5
├── evidence      # Standalone RNA-seq evidence scoring
├── homology      # Protein homology validation
│   ├── extract-proteins
│   ├── search
│   └── validate
├── validate      # Quick protein validation
├── qc            # Quality control and reporting
│   ├── aggregate
│   ├── report
│   ├── filter
│   ├── tiered-output
│   └── list-flags
├── parallel      # Large genome parallelization
│   ├── plan
│   ├── tasks
│   ├── aggregate
│   └── suggest
└── viz           # Visualization (optional)
```

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
                          │   homology validate │
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

---

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

The `refine` command is the **primary HelixForge pipeline**. It combines splice correction, boundary adjustment, confidence scoring, and evidence scoring into a single integrated workflow.

```
┌─────────────────────────────────────────────────────────────────┐
│                     refine Pipeline                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Input: HDF5 + GFF3 + BAM(s)                                    │
│                    │                                             │
│                    ▼                                             │
│  ┌─────────────────────────────┐                                │
│  │  1. Junction Extraction     │  Extract splice junctions      │
│  │     from RNA-seq BAMs       │  from all BAM files            │
│  └─────────────┬───────────────┘                                │
│                ▼                                                 │
│  ┌─────────────────────────────┐                                │
│  │  2. Splice Correction       │  Shift splice sites to match  │
│  │     (max_shift, min_reads)  │  RNA-seq evidence + PWM       │
│  └─────────────┬───────────────┘                                │
│                ▼                                                 │
│  ┌─────────────────────────────┐                                │
│  │  3. Boundary Adjustment     │  Find optimal start/stop      │
│  │     (Kozak context scoring) │  codons in search window      │
│  └─────────────┬───────────────┘                                │
│                ▼                                                 │
│  ┌─────────────────────────────┐                                │
│  │  4. Confidence Scoring      │  Score from HDF5 softmax      │
│  │     (from HDF5 predictions) │  probabilities                │
│  └─────────────┬───────────────┘                                │
│                ▼                                                 │
│  ┌─────────────────────────────┐                                │
│  │  5. Evidence Scoring        │  Calculate AED, junction      │
│  │     (AED, junction support) │  support ratio, coverage      │
│  └─────────────┬───────────────┘                                │
│                ▼                                                 │
│  Output: Refined GFF3 + Report TSV                              │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

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
| Option | Description | Default |
|--------|-------------|---------|
| `-p/--helixer-h5` | Helixer HDF5 predictions | Required |
| `-g/--helixer-gff` | Helixer GFF3 predictions | Required |
| `--genome` | Reference genome FASTA | Required |
| `--rnaseq-bam` | BAM file(s), comma-separated | - |
| `--rnaseq-bam-list` | File containing BAM paths | - |
| `--junctions-bed` | Pre-computed junction files | - |
| `--max-shift` | Maximum bp to shift splice sites | 15 |
| `--min-reads` | Minimum reads for junction support | 3 |
| `--min-tissues` | Minimum samples supporting junction | 1 |
| `--adjust-boundaries` | Refine start/stop codons | on |
| `-r/--report` | Output comprehensive report TSV | - |

**Output GFF3 attributes:**
- `confidence_score` - Overall confidence (0-1)
- `evidence_score` - RNA-seq evidence score (0-1)
- `aed` - Annotation Edit Distance (0-1, lower is better)
- `junction_support` - Fraction of junctions supported
- `mean_coverage` - Mean exon coverage
- `flags` - Quality flags

---

### Alternative: Standalone Commands

If you don't have RNA-seq data or prefer separate steps:

#### Confidence Scoring Only (requires HDF5)

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

```
┌─────────────────────────────────────────────────────────────────┐
│                   Homology Validation                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌─────────────────┐    ┌─────────────────┐                     │
│  │ extract-proteins│───▶│    proteins.fa  │                     │
│  │  (from GFF)     │    │                 │                     │
│  └─────────────────┘    └────────┬────────┘                     │
│                                  │                               │
│                                  ▼                               │
│  ┌─────────────────┐    ┌─────────────────┐                     │
│  │ Reference DB    │───▶│     search      │──▶ hits.tsv        │
│  │ (UniProt, etc)  │    │  (Diamond/MMseqs)│                    │
│  └─────────────────┘    └─────────────────┘                     │
│                                  │                               │
│                                  ▼                               │
│                         ┌─────────────────┐                     │
│                         │    validate     │──▶ validation.tsv   │
│                         │ (chimera/frag   │                     │
│                         │  detection)     │                     │
│                         └─────────────────┘                     │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

```bash
# Extract proteins from refined predictions
helixforge homology extract-proteins \
    --gff refined.gff3 \
    --genome genome.fa \
    -o predicted_proteins.fa \
    --longest-isoform

# Search against database
helixforge homology search \
    --query predicted_proteins.fa \
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

**Note:** The validation output includes ALL genes from the GFF, not just those with homology hits. Genes without hits are marked as `status: no_hit` with the `no_homology` flag. This ensures orphan genes and species-specific genes are retained in QC aggregation.

---

### Step 3: Aggregate QC Results

Combine results from refine and homology into unified QC metrics. This step calculates a **combined AED** that incorporates RNA-seq evidence, homology evidence, and Helixer confidence.

```bash
# Using refine report (recommended)
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv homology_validation.tsv \
    -o qc_aggregated.tsv

# Custom AED weights (for species with limited RNA-seq data)
helixforge qc aggregate \
    --refine-tsv refine_report.tsv \
    --homology-tsv homology_validation.tsv \
    -o qc_aggregated.tsv \
    --aed-rnaseq-weight 0.3 \
    --aed-homology-weight 0.5 \
    --aed-confidence-weight 0.2
```

**AED weight options:**
| Option | Default | Description |
|--------|---------|-------------|
| `--aed-rnaseq-weight` | 0.4 | Weight for RNA-seq evidence |
| `--aed-homology-weight` | 0.4 | Weight for homology evidence |
| `--aed-confidence-weight` | 0.2 | Weight for Helixer confidence |

**Why combined AED matters:**
The traditional AED (from `refine`) is based only on RNA-seq evidence. This can be misleading for genes with good homology but no expression in sampled tissues. The combined AED ensures genes are evaluated using all available evidence.

| Scenario | RNA-seq AED | Combined AED |
|----------|-------------|--------------|
| Good RNA-seq, good homology | 0.1 | ~0.10 |
| No RNA-seq, complete homology | 0.5 | ~0.26 |
| Good RNA-seq, no homology | 0.2 | ~0.28 |
| Neither | 0.5 | ~0.58 |

**Output columns:**
- `gene_id`, `tier`, `flag_count`, `flag_codes`, `max_severity`
- `confidence_score`, `splice_score`, `homology_score`
- `rnaseq_aed` - RNA-seq only AED (from refine report)
- `combined_aed` - Multi-evidence AED (RNA-seq + homology + confidence)

---

### Step 4: Generate QC Report

Create a comprehensive HTML report summarizing annotation quality.

```bash
helixforge qc report \
    --qc-tsv qc_aggregated.tsv \
    -o qc_report.html \
    --title "Genome Annotation QC"
```

**Report Contents:**

| Section | Description |
|---------|-------------|
| **Summary** | Tier counts, average scores, mean RNA-seq AED, mean Combined AED |
| **Distribution Charts** | Tier pie chart, confidence histogram, AED histogram (RNA-seq vs Combined) |
| **Flag Analysis** | Top 15 flags by gene count, severity distribution, category distribution |
| **QC Flags Reference** | Table of all possible flags with descriptions and severity levels |
| **Gene Details** | Searchable table (limited to 1000 genes for HTML performance) |

> **Note:** The gene table is limited to 1000 genes to keep the HTML file size manageable. Use the TSV output from `qc aggregate` for complete data on all genes.

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
- `high_confidence.gff3` - Highest quality genes (conf ≥0.85, no errors)
- `medium_confidence.gff3` - Moderate quality genes (conf ≥0.70)
- `low_confidence.gff3` - Lower quality, needs review
- `rejected.gff3` - Problematic genes (TE overlap, internal stops)
- `tier_summary.tsv` - Summary statistics

---

## Parallel Processing for Large Genomes

For genomes larger than ~500 Mb, use parallel processing:

```bash
# 1. Create chunk plan
helixforge parallel plan \
    --genome genome.fa \
    --target-size 50000000 \
    -o chunks.json

# 2. Generate task file
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command 'helixforge refine -p predictions.h5 -g predictions.gff3 --genome genome.fa --rnaseq-bam rna.bam --region {seqid}:{start}-{end} --chunk-id {chunk_id} -o outputs/{chunk_id}.gff -r outputs/{chunk_id}.tsv' \
    --output tasks.txt \
    --output-dir outputs/

# 3. Execute with HyperShell or GNU Parallel
hs cluster tasks.txt --num-tasks 32
# or
parallel -j 32 < tasks.txt

# 4. Aggregate chunk results
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*.gff' \
    -o refined_combined.gff3 \
    --type merge_gff
```

---

## Complete Example Pipeline

```bash
#!/bin/bash
set -euo pipefail

# Input files
PREDICTIONS_H5="helixer_output/predictions.h5"
PREDICTIONS_GFF="helixer_output/predictions.gff3"
GENOME="genome.fa"
RNASEQ_BAM="rnaseq.bam"
PROTEINS="uniprot_plants.dmnd"

# Output directory
OUTDIR="helixforge_output"
mkdir -p "$OUTDIR"

echo "Step 1: Refine predictions with RNA-seq evidence..."
helixforge refine \
    -p "$PREDICTIONS_H5" \
    -g "$PREDICTIONS_GFF" \
    --genome "$GENOME" \
    --rnaseq-bam "$RNASEQ_BAM" \
    -o "$OUTDIR/refined.gff3" \
    -r "$OUTDIR/refine_report.tsv" \
    --verbose

echo "Step 2: Extract and search proteins..."
helixforge homology extract-proteins \
    --gff "$OUTDIR/refined.gff3" \
    --genome "$GENOME" \
    -o "$OUTDIR/proteins.fa"

helixforge homology search \
    --query "$OUTDIR/proteins.fa" \
    --database "$PROTEINS" \
    -o "$OUTDIR/homology_hits.tsv" \
    --threads 8

echo "Step 3: Validate with homology..."
helixforge homology validate \
    --search-results "$OUTDIR/homology_hits.tsv" \
    --gff "$OUTDIR/refined.gff3" \
    -o "$OUTDIR/homology_validation.tsv"

echo "Step 4: Aggregate QC results..."
helixforge qc aggregate \
    --refine-tsv "$OUTDIR/refine_report.tsv" \
    --homology-tsv "$OUTDIR/homology_validation.tsv" \
    -o "$OUTDIR/qc_aggregated.tsv"

echo "Step 5: Generate QC report..."
helixforge qc report \
    --qc-tsv "$OUTDIR/qc_aggregated.tsv" \
    -o "$OUTDIR/qc_report.html"

echo "Step 6: Create tiered output..."
helixforge qc tiered-output \
    --qc-results "$OUTDIR/qc_aggregated.tsv" \
    --gff "$OUTDIR/refined.gff3" \
    -o "$OUTDIR/tiered/"

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
helixforge homology extract-proteins --gff predictions.gff3 --genome genome.fa -o proteins.fa
helixforge homology search --query proteins.fa --database uniprot.dmnd -o hits.tsv
helixforge homology validate --search-results hits.tsv --gff predictions.gff3 -o validation.tsv

# 3. QC and filter
helixforge qc aggregate \
    --confidence-tsv confidence.tsv \
    --homology-tsv validation.tsv \
    -o qc_aggregated.tsv

helixforge qc tiered-output \
    --qc-results qc_aggregated.tsv \
    --gff predictions.gff3 \
    -o tiered/
```

---

## Command Summary

| Command | Purpose | Requires |
|---------|---------|----------|
| `refine` | **Main pipeline** - splice correction, boundaries, confidence, evidence | HDF5 + BAM |
| `confidence` | Score genes using Helixer HDF5 probabilities | HDF5 |
| `evidence` | Score genes using RNA-seq evidence only | BAM |
| `homology extract-proteins` | Extract protein sequences from GFF | GFF + genome |
| `homology search` | Search proteins against database | proteins + DB |
| `homology validate` | Validate genes with homology results | hits TSV |
| `validate` | Quick protein validation (all-in-one) | GFF + proteins |
| `qc aggregate` | Combine all QC metrics | TSV files |
| `qc report` | Generate HTML report | aggregated TSV |
| `qc tiered-output` | Create quality-filtered GFFs | aggregated TSV + GFF |
| `qc filter` | Filter genes by criteria | aggregated TSV |
| `parallel plan` | Create genome chunk plan | genome |
| `parallel tasks` | Generate task file for parallel execution | chunk plan |
| `parallel aggregate` | Combine chunked outputs | chunk outputs |
| `viz` | Visualize genes and evidence | GFF + genome |

**For large genomes:** Use `parallel plan`, `parallel tasks`, and `parallel aggregate` to distribute work across HPC nodes.
