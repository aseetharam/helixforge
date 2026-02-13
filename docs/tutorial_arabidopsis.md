# Real-World Tutorial: Arabidopsis thaliana Annotation with HelixForge

This tutorial demonstrates a complete HelixForge workflow using *Arabidopsis thaliana* with 56 RNA-seq samples from diverse tissues. This serves as a reference for testing all HelixForge features and debugging edge cases.

## Overview

| Item | Details |
| ------ | --------- |
| **Species** | *Arabidopsis thaliana* (TAIR10/Araport11) |
| **Genome size** | ~135 Mb (5 chromosomes + organelles) |
| **RNA-seq** | 56 samples, paired-end, 2x76 bp (PRJEB32665) |
| **Tissues** | Roots, leaves, flowers, seeds, embryos, etc. |
| **Reference annotation** | Araport11 (for validation comparison) |

## Directory Structure

```text
arabidopsis_helixforge/
├── genome/
│   └── athaliana.fasta              # Reference genome (TAIR10)
├── helixer_output/
│   ├── Arabidopsis-thaliana_input.h5
│   ├── Arabidopsis-thaliana_predictions.h5
│   └── Arabidopsis-thaliana_helixer.gff3
├── rnaseq/
│   ├── fastq_files/                 # Raw FASTQ files
│   ├── star_index_araport11/        # STAR index
│   └── sample_list.txt              # Sample names
├── bam/                             # Aligned BAM files
├── databases/                       # Protein databases
├── outputs/                         # HelixForge chunk outputs
├── results/                         # Final aggregated results
└── scripts/                         # Helper scripts
```

---

## Part 1: Prerequisites

### 1.1 Software Requirements

```bash
# HelixForge (install from source)
pip install -e /path/to/helixforge

# External tools
module load biocontainers
module load star
module load samtools
module load diamond

# Or via conda
conda install -c bioconda star samtools diamond
```

### 1.2 Download Reference Data

```bash
# Create working directory
mkdir -p arabidopsis_helixforge && cd arabidopsis_helixforge
mkdir -p genome helixer_output rnaseq bam databases outputs results scripts

# Download TAIR10 genome
cd genome
wget https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
# need to gunzip and then rename
gunzip download\?filePath\=Genes%2FTAIR10_genome_release%2FTAIR10_chromosome_files%2FTAIR10_chr_all.fas.gz 
mv download\?filePath\=Genes%2FTAIR10_genome_release%2FTAIR10_chromosome_files%2FTAIR10_chr_all.fas athaliana.fasta
samtools faidx athaliana.fasta
cd ..

# Download Araport11 annotation (for comparison)
wget https://www.arabidopsis.org/api/download-files/download?filePath=Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.20250813.gff.gz
gunzip download?filePath=Genes%2FAraport11_genome_release%2FAraport11_GFF3_genes_transposons.20250813.gff
mv download\?filePath\=Genes%2FAraport11_genome_release%2FAraport11_GFF3_genes_transposons.20250813.gff reference_araport11.gff3
```

---

## Part 2: Helixer Gene Prediction

### 2.1 Pull Helixer Container

```bash
apptainer pull docker://gglyptodon/helixer-docker:helixer_v0.3.6_cuda_12.2.2-cudnn8_1
```

### 2.2 Run Helixer (Preserving H5 Files)

**Important:** Run Helixer steps separately to preserve the H5 predictions file required by HelixForge.

```bash
#!/bin/bash
# scripts/run_helixer.sh

# Configuration
GENOME="genome/athaliana.fasta"
SPECIES="Arabidopsis-thaliana"
SIF="helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8_1.sif"
OUTDIR="helixer_output"
MODEL_PATH="${HOME}/.local/share/Helixer/models/land_plant/land_plant_v0.3_a_0080.h5"

mkdir -p ${OUTDIR}

# Step 1: Convert FASTA to H5
echo "Step 1: Converting FASTA to H5..."
apptainer exec ${SIF} fasta2h5.py \
    --species ${SPECIES} \
    --h5-output-path ${OUTDIR}/${SPECIES}_input.h5 \
    --fasta-path ${GENOME}

# Step 2: Run neural network predictions (requires GPU)
echo "Step 2: Running predictions..."
apptainer exec --nv ${SIF} HybridModel.py \
    --load-model-path ${MODEL_PATH} \
    --test-data ${OUTDIR}/${SPECIES}_input.h5 \
    --prediction-output-path ${OUTDIR}/${SPECIES}_predictions.h5 \
    --batch-size 8

# Step 3: Generate GFF3 from predictions
echo "Step 3: Generating GFF3..."
apptainer exec ${SIF} helixer_post_bin ${OUTDIR}/${SPECIES}_input.h5 ${OUTDIR}/${SPECIES}_predictions.h5 100 0.1 0.8 100 ${OUTDIR}/${SPECIES}_helixer.gff3

echo "Helixer complete. Files:"
ls -lh ${OUTDIR}/
```

**Expected output files:**

- `Arabidopsis-thaliana_input.h5` (~123 MB) - Coordinate mapping
- `Arabidopsis-thaliana_predictions.h5` (~2.5 GB) - Neural network predictions
- `Arabidopsis-thaliana_helixer.gff3` - Gene predictions

---

## Part 3: RNA-seq Alignment with STAR

### 3.1 Build STAR Index

```bash
#!/bin/bash
# scripts/build_star_index.sh

module load biocontainers star

GENOME="genome/athaliana.fasta"
GFF="reference_araport11.gff3"
INDEX_DIR="rnaseq/star_index_araport11"

mkdir -p ${INDEX_DIR}

STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir ${INDEX_DIR} \
     --genomeFastaFiles ${GENOME} \
     --sjdbGTFfile ${GFF} \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbOverhang 75 \
     --genomeSAindexNbases 12
```

### 3.2 Create Sample List

```bash
# List all sample prefixes (without _R1/_R2 suffixes)
ls rnaseq/fastq_files/*_R1.fq.gz | xargs -n1 basename | sed 's/_R1.fq.gz//' > rnaseq/sample_list.txt

# Verify sample count (should be 56)
wc -l rnaseq/sample_list.txt
```

**Sample list contents (56 samples):**

```text
sepal_flwr15
stem_node1_adult
stem_internode2_adult
cauline_leaf_adult
rosette_lf7_distal_adult
...
pollen_flwr15
```

### 3.3 STAR Alignment Array Job

```bash
#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --array=1-56
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/star_%A_%a.log

# scripts/star_array.sh

set -euo pipefail

# Configuration
STAR_INDEX="rnaseq/star_index_araport11"
FASTQ_DIR="rnaseq/fastq_files"
OUT_DIR="bam"
SAMPLE_LIST="rnaseq/sample_list.txt"

# Get sample name for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No sample found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing sample: ${SAMPLE}"
echo "Start time: $(date)"

# Input files
R1="${FASTQ_DIR}/${SAMPLE}_R1.fq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.fq.gz"

# Verify inputs exist
for f in "${R1}" "${R2}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Input file not found: ${f}"
        exit 1
    fi
done

module --force purge
module load biocontainers star samtools

# Run STAR alignment
STAR --runMode alignReads \
     --runThreadN ${SLURM_CPUS_PER_TASK} \
     --genomeDir "${STAR_INDEX}" \
     --readFilesIn "${R1}" "${R2}" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS NM MD XS \
     --outSAMstrandField intronMotif \
     --twopassMode Basic \
     --alignIntronMin 20 \
     --alignIntronMax 10000 \
     --alignMatesGapMax 10000 \
     --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 6 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --sjdbOverhang 75 \
     --quantMode GeneCounts

# Index BAM
samtools index -@ ${SLURM_CPUS_PER_TASK} "${OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

# Cleanup STAR temp files
rm -rf "${OUT_DIR}/${SAMPLE}__STARgenome" "${OUT_DIR}/${SAMPLE}__STARpass1"

echo "Finished: ${SAMPLE} at $(date)"
```

**Submit the job:**

```bash
mkdir -p logs bam
sbatch scripts/star_array.sh
```

### 3.4 Verify Alignments

```bash
# Check all BAMs exist and are indexed
ls bam/*.bam | wc -l           # Should be 56
ls bam/*.bam.bai | wc -l       # Should be 56

# Create BAM file list for HelixForge
ls bam/*.bam > bam_files.txt
```

---

## Part 4: HelixForge Refinement

### 4.1 Create Chunk Plan

```bash
# Plan genome chunks for parallel processing
helixforge parallel plan \
    --genome genome/athaliana.fasta \
    --gff helixer_output/Arabidopsis-thaliana_helixer.gff3 \
    --strategy adaptive \
    --target-chunks 20 \
    -o chunks.json 

# Expected output: ~20 chunks covering all chromosomes
```

### 4.2 Generate Task File

```bash
# Create output directory
mkdir -p outputs

# Generate tasks for parallel execution
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge refine \
        --helixer-h5 helixer_output/Arabidopsis-thaliana_predictions.h5 \
        --helixer-gff helixer_output/Arabidopsis-thaliana_helixer.gff3 \
        --rnaseq-bam-list bam_files.txt \
        --genome genome/athaliana.fasta \
        --region {seqid}:{start}-{end} \
        --chunk-id {chunk_id} \
        --output outputs/{chunk_id}_refined.gff3 \
        --report outputs/{chunk_id}_refine_report.tsv \
        --splice-details outputs/{chunk_id}_splice_details.tsv \
        --max-shift 15 \
        --min-reads 5 \
        --min-tissues 2 \
        --adjust-boundaries \
        --min-coverage 5 \
        --verbose" \
    --output tasks.txt \
    --output-dir outputs/

# View generated tasks
head -3 tasks.txt
```

### 4.3 Execute Tasks

#### Option A: Using HyperShell (recommended for HPC)

```bash
# Install HyperShell if needed
pip install hypershell

# Run with 8 parallel tasks
hs cluster tasks.txt --num-tasks 8
```

#### Option B: Using GNU Parallel

```bash
parallel -j 8 < tasks.txt
```

#### Option C: SLURM Array Job

```bash
#!/bin/bash
#SBATCH --job-name=helixforge_refine
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/refine_%A_%a.log

# scripts/refine_array.sh

TASK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tasks.txt)
eval ${TASK}
```

### 4.4 Aggregate Chunk Results

```bash
# Merge GFF files
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*_refined.gff3' \
    -o results/refined_genes.gff3 \
    --type merge_gff

# Merge refine reports
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*_refine_report.tsv' \
    -o results/refine_report.tsv \
    --type merge_tsv

# Merge splice details
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*_splice_details.tsv' \
    -o results/splice_details.tsv \
    --type merge_tsv

# Verify outputs
wc -l results/*.tsv
grep -c "gene" results/refined_genes.gff3
```

---

## Part 5: Homology Validation

### 5.1 Setup Protein Database

```bash
# List available databases
helixforge homology list-databases

# Download Swiss-Prot plants subset
helixforge homology download-db \
    --database swissprot_plants \
    --output-dir databases/

# Format for Diamond
helixforge homology format-db \
    --input databases/swissprot_plants/uniprot_sprot_plants.fasta \
    --output databases/swissprot_plants.dmnd \
    --tool diamond \
    --threads 12
```

### 5.2 Extract and Search Proteins

```bash
# Extract protein sequences from refined genes
helixforge homology extract-proteins \
    --gff results/refined_genes.gff3 \
    --genome genome/athaliana.fasta \
    -o results/predicted_proteins.fa \
    --longest-isoform

# Search against Swiss-Prot
helixforge homology search \
    --query results/predicted_proteins.fa \
    --database databases/swissprot_plants.dmnd \
    -o results/homology_hits.tsv \
    --threads 16 \
    --evalue 1e-5
```

### 5.3 Validate Gene Models

```bash
# Run validation (includes genes without hits)
helixforge homology validate \
    --search-results results/homology_hits.tsv \
    --gff results/refined_genes.gff3 \
    -o results/homology_validation.tsv \
    --thresholds relaxed \
    --chimera-report results/chimeras.tsv \
    --fragment-report results/fragments.tsv \
    --output-gff results/genes_with_homology.gff3

# Check results
echo "Validation summary:"
head -1 results/homology_validation.tsv
wc -l results/homology_validation.tsv
wc -l results/chimeras.tsv
wc -l results/fragments.tsv
```

---

## Part 6: QC Aggregation and Reporting

### 6.1 Aggregate QC Results

```bash
# Combine all evidence into unified QC metrics
helixforge qc aggregate \
    --refine-tsv results/refine_report.tsv \
    --homology-tsv results/homology_validation.tsv \
    -o results/qc_aggregated.tsv \
    --aed-rnaseq-weight 0.4 \
    --aed-homology-weight 0.4 \
    --aed-confidence-weight 0.2

# Check output columns
head -1 results/qc_aggregated.tsv | tr '\t' '\n'
```

**Output columns explained:**

- `gene_id` - Gene identifier
- `tier` - Quality tier (high/medium/low/reject)
- `confidence_score` - Helixer prediction confidence
- `splice_score` - Junction support ratio
- `homology_score` - Homology coverage × identity
- `rnaseq_aed` - RNA-seq only AED
- `combined_aed` - Multi-evidence AED
- `flag_codes` - Quality flags

### 6.2 Generate QC Report

```bash
# Generate interactive HTML report
helixforge qc report \
    --qc-tsv results/qc_aggregated.tsv \
    -o results/qc_report.html \
    --title "Arabidopsis thaliana Annotation QC - HelixForge"

# Open in browser
firefox results/qc_report.html &
```

**Report sections:**

- Summary statistics with AED means
- Tier distribution pie chart
- Confidence and AED histograms
- Top 15 flags by occurrence
- Searchable gene table (limited to 1000)

### 6.3 Create Tiered Output

```bash
# Generate quality-filtered GFF files
helixforge qc tiered-output \
    --qc-results results/qc_aggregated.tsv \
    --gff results/refined_genes.gff3 \
    -o results/tiered/ \
    --prefix athaliana

# Check results
ls -la results/tiered/
wc -l results/tiered/*.gff3
```

**Output files:**

- `athaliana_high.gff3` - Publication-ready genes
- `athaliana_medium.gff3` - Usable with caveats
- `athaliana_low.gff3` - Requires manual review
- `athaliana_rejected.gff3` - Likely false positives
- `athaliana_summary.tsv` - Tier statistics

---

## Part 7: Quality Assessment

### 7.1 Compare with Reference Annotation

```bash
# Count genes per tier
echo "=== Gene Counts by Tier ==="
for tier in high medium low rejected; do
    count=$(grep -c "^[^#]" results/tiered/athaliana_${tier}.gff3 2>/dev/null | grep "gene" || echo 0)
    echo "${tier}: ${count}"
done

# Compare with Araport11 gene count
araport_genes=$(grep -c "gene" reference_araport11.gff3)
helixforge_genes=$(grep -c "gene" results/refined_genes.gff3)
echo ""
echo "Araport11 genes: ${araport_genes}"
echo "HelixForge genes: ${helixforge_genes}"
```

### 7.2 Analyze Flag Distribution

```bash
# Extract flag statistics from QC TSV
echo "=== Most Common Flags ==="
cut -f4 results/qc_aggregated.tsv | tail -n +2 | tr ',' '\n' | \
    sort | uniq -c | sort -rn | head -20
```

### 7.3 AED Distribution Analysis

```bash
# RNA-seq AED distribution
echo "=== RNA-seq AED Distribution ==="
cut -f9 results/qc_aggregated.tsv | tail -n +2 | \
    awk '{if($1!="") print int($1*10)/10}' | sort | uniq -c

# Combined AED distribution
echo ""
echo "=== Combined AED Distribution ==="
cut -f10 results/qc_aggregated.tsv | tail -n +2 | \
    awk '{if($1!="") print int($1*10)/10}' | sort | uniq -c
```

---

## Part 8: Testing Edge Cases

### 8.1 Single Chromosome Processing

```bash
# Test on chromosome 1 only
helixforge refine \
    --helixer-h5 helixer_output/Arabidopsis-thaliana_predictions.h5 \
    --helixer-gff helixer_output/Arabidopsis-thaliana_helixer.gff3 \
    --rnaseq-bam-list bam_files.txt \
    --genome genome/athaliana.fasta \
    --scaffold Chr1 \
    --output results/chr1_test.gff3 \
    --report results/chr1_test_report.tsv \
    --verbose
```

### 8.2 Single BAM Processing

```bash
# Test with single tissue
helixforge refine \
    --helixer-h5 helixer_output/Arabidopsis-thaliana_predictions.h5 \
    --helixer-gff helixer_output/Arabidopsis-thaliana_helixer.gff3 \
    --rnaseq-bam bam/flower_flwr15_Aligned.sortedByCoord.out.bam \
    --genome genome/athaliana.fasta \
    --scaffold Chr1 \
    --output results/chr1_single_bam.gff3 \
    --report results/chr1_single_bam_report.tsv
```

### 8.3 No RNA-seq (Confidence Only)

```bash
# Test confidence scoring without RNA-seq
helixforge confidence \
    --helixer-h5 helixer_output/Arabidopsis-thaliana_predictions.h5 \
    --helixer-gff helixer_output/Arabidopsis-thaliana_helixer.gff3 \
    --genome genome/athaliana.fasta \
    --scaffold Chr1 \
    -o results/chr1_confidence_only.tsv
```

### 8.4 Custom AED Weights

```bash
# Test with homology-prioritized weights (for species with limited RNA-seq)
helixforge qc aggregate \
    --refine-tsv results/refine_report.tsv \
    --homology-tsv results/homology_validation.tsv \
    -o results/qc_homology_priority.tsv \
    --aed-rnaseq-weight 0.2 \
    --aed-homology-weight 0.6 \
    --aed-confidence-weight 0.2

# Compare AED distributions
echo "Default weights:"
cut -f10 results/qc_aggregated.tsv | tail -n +2 | awk '{sum+=$1; n++} END {print "Mean:", sum/n}'

echo "Homology-priority weights:"
cut -f10 results/qc_homology_priority.tsv | tail -n +2 | awk '{sum+=$1; n++} END {print "Mean:", sum/n}'
```

### 8.5 Filter Presets

```bash
# Test different filter presets
for preset in high_confidence publication_ready has_homology needs_review; do
    helixforge qc filter \
        --gff results/refined_genes.gff3 \
        --qc-results results/qc_aggregated.tsv \
        -o results/filtered_${preset}.gff3 \
        --preset ${preset}

    count=$(grep -c "gene" results/filtered_${preset}.gff3)
    echo "${preset}: ${count} genes"
done
```

---

## Part 9: Troubleshooting

### 9.1 Memory Issues

```bash
# For large genomes, increase chunk count
helixforge parallel plan \
    --genome genome/athaliana.fasta \
    --strategy size \
    --target-size 10000000 \  # 10 Mb chunks
    -o chunks_small.json
```

### 9.2 Missing H5 File

If Helixer H5 file was not preserved:

```bash
# Re-run Helixer Step 2 only (requires original input.h5)
apptainer exec --nv ${SIF} HybridModel.py \
    --load-model-path ${MODEL_PATH} \
    --test-data helixer_output/Arabidopsis-thaliana_input.h5 \
    --prediction-output-path helixer_output/Arabidopsis-thaliana_predictions.h5 \
    --batch-size 8
```

### 9.3 BAM Index Missing

```bash
# Index all BAM files
for bam in bam/*.bam; do
    if [[ ! -f "${bam}.bai" ]]; then
        echo "Indexing: ${bam}"
        samtools index "${bam}"
    fi
done
```

### 9.4 Check Log Files

```bash
# Review logs for errors
grep -i "error\|warning\|fail" logs/*.log | head -50
```

---

## Part 10: Expected Results Summary

For Arabidopsis thaliana with 56 tissue RNA-seq samples:

| Metric | Expected Range |
| -------- | ---------------- |
| Total predicted genes | 25,000 - 30,000 |
| High confidence | 60-75% |
| Medium confidence | 15-25% |
| Low confidence | 5-10% |
| Rejected | 2-5% |
| Mean confidence score | 0.65 - 0.75 |
| Mean RNA-seq AED | 0.25 - 0.40 |
| Mean combined AED | 0.20 - 0.35 |
| Genes with homology | 85-95% |
| Chimeric genes | < 1% |
| Fragmented genes | 1-3% |

---

## Quick Reference Commands

```bash
# Full pipeline (copy-paste ready)

# 1. Refine with parallel processing
helixforge parallel plan --genome genome/athaliana.fasta --gff helixer_output/Arabidopsis-thaliana_helixer.gff3 --strategy adaptive --target-chunks 20 -o chunks.json
helixforge parallel tasks --chunk-plan chunks.json --command "helixforge refine --helixer-h5 helixer_output/Arabidopsis-thaliana_predictions.h5 --helixer-gff helixer_output/Arabidopsis-thaliana_helixer.gff3 --rnaseq-bam-list bam_files.txt --genome genome/athaliana.fasta --region {seqid}:{start}-{end} --chunk-id {chunk_id} --output outputs/{chunk_id}_refined.gff3 --report outputs/{chunk_id}_report.tsv --min-reads 5 --min-tissues 2" --output tasks.txt
hs cluster tasks.txt --num-tasks 8

# 2. Aggregate results
helixforge parallel aggregate --input-dir outputs/ --pattern '*_refined.gff3' -o results/refined_genes.gff3 --type merge_gff
helixforge parallel aggregate --input-dir outputs/ --pattern '*_report.tsv' -o results/refine_report.tsv --type merge_tsv

# 3. Homology validation
helixforge homology extract-proteins --gff results/refined_genes.gff3 --genome genome/athaliana.fasta -o results/proteins.fa --longest-isoform
helixforge homology search --query results/proteins.fa --database databases/swissprot_plants.dmnd -o results/hits.tsv --threads 16
helixforge homology validate --search-results results/hits.tsv --gff results/refined_genes.gff3 -o results/homology.tsv

# 4. QC and reporting
helixforge qc aggregate --refine-tsv results/refine_report.tsv --homology-tsv results/homology.tsv -o results/qc.tsv
helixforge qc report --qc-tsv results/qc.tsv -o results/report.html
helixforge qc tiered-output --qc-results results/qc.tsv --gff results/refined_genes.gff3 -o results/tiered/
```

---

## Contact and Support

For issues or questions:

- GitHub: [HelixForge Issues](https://github.com/aseetharam/helixforge/issues)
- Documentation: See `docs/workflow.md` for detailed command references
