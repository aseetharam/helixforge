# `helixforge confidence`

Inspect (read-only): score genes against the Helixer HDF5 confidence track.

A standalone, HDF5-only scorer: it reads the Helixer softmax predictions and
computes multi-factor confidence metrics per gene (class probabilities,
Shannon entropy, boundary sharpness, CDS coding consistency, per-exon
scores). It needs *only* the Helixer HDF5 + a GFF3, no Mikado, no RNA-seq,
no external toolchain. Use ``evidence`` instead to score against RNA-seq /
protein evidence, and ``reconcile`` to actually build/fix models; this
command never modifies a model, it only annotates confidence.

Genes are classed high (>=0.85), medium (>=0.70), or low (<0.70), but the
scores are genome-relative, prefer a cutoff from the printed distribution
(or ``--summary-tsv``) over the fixed class thresholds.

Coordinate mapping comes from (in order of preference): ``--input-h5``
(strand-aware), auto-detected ``*_input.h5`` next to the predictions, or
``--genome`` FASTA + ``.fai``. For parallel runs use ``--region`` /
``--scaffold`` to subset and ``--chunk-id`` for output organization.

Outputs:
- -o/--output TSV, one row per gene, columns:
  gene_id, seqid, start, end, strand, mean_prob, min_prob, median_prob,
  entropy, boundary_sharpness, coding_consistency, worst_exon_score,
  overall_score, confidence_class, flags, n_low_conf_regions, n_exons
- --bed: BED12 of genes coloured by confidence; --low-conf-bed: BED of
  low-confidence sub-regions
- --summary-tsv: long-format metric/value distribution (n, mean, median,
  std, min, max, p5/p25/p50/p75/p95; per-component and per-exon stats;
  n_high/n_medium/n_low)
- --distribution-plot / --plot-dir: genome-wide and per-gene plots

## Options

| Option | Description |
|---|---|
| `-p`, `--predictions` | Helixer HDF5 predictions file. (required) |
| `-g`, `--gff` | GFF3 file with gene predictions. (required) |
| `--genome` | Reference genome FASTA file. Optional if --input-h5 is provided or auto-detected. |
| `--input-h5` | Helixer input HDF5 file (contains coordinate mapping). If not provided, will auto-detect from predictions file location. Required for strand-aware prediction retrieval. |
| `-o`, `--output` | Output TSV file for confidence scores. (required) |
| `--bed` | Output BED file for genome browser visualization. |
| `--low-conf-bed` | Output BED file for low-confidence regions within genes. |
| `--plot-dir` | Directory for per-gene confidence plots. Use with --plot-threshold or --max-plots to limit output for large genomes. |
| `--plot-threshold` | Only plot genes with confidence below this threshold (e.g., 0.7). Recommended for large genomes. |
| `--max-plots` | Maximum number of gene plots to generate. Plots lowest-confidence genes first. |
| `--distribution-plot` | Output path for genome-wide distribution plot. |
| `--summary-tsv` | Write the on-screen distribution summary (metric/value, long format) to this TSV. |
| `--threshold` | Threshold for low-confidence regions. (default: `0.7`) |
| `-j`, `--threads` | Number of parallel threads. (default: `1`) |
| `--format` | Format for visualization output. (choices: html/png/pdf; default: `html`) |
| `--region` | Process only this region (format: seqid:start-end, 1-based inclusive). |
| `--chunk-id` | Chunk identifier for logging and output naming in parallel mode. |
| `--scaffold` | Process only this scaffold (simpler alternative to --region). |

## Examples

```bash
  # explicit input HDF5 (preferred, enables strand-aware coordinate mapping)
  helixforge confidence -p predictions.h5 -g genes.gff3 \
      --input-h5 input.h5 -o scores.tsv

  # auto-detect *_input.h5 next to the predictions file
  helixforge confidence -p predictions.h5 -g genes.gff3 -o scores.tsv

  # FASTA-index fallback + a genome-wide distribution plot and summary TSV
  helixforge confidence -p predictions.h5 -g genes.gff3 --genome genome.fa \
      -o scores.tsv --distribution-plot dist.html --summary-tsv summary.tsv

  # one scaffold only, low-confidence plots, for a parallel chunk
  helixforge confidence -p predictions.h5 -g genes.gff3 --input-h5 input.h5 \
      -o chr1.scores.tsv --scaffold chr1 --plot-dir plots --plot-threshold 0.7 \
      --chunk-id chr1
```
