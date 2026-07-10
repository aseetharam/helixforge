# `helixforge utils qc`

Genome-wide QC report (HTML/JSON/TSV).

## Options

| Option | Description |
|---|---|
| `--gff3` | Reconciled GFF3. (required) |
| `--genome` | Reference genome FASTA. |
| `--helixer-h5` | Helixer HDF5. |
| `--out` | Output report path. (required) |
| `--format` | Output format [default: html]. (choices: html/json/tsv; default: `html`) |

## Examples

```bash
  helixforge utils qc --gff3 helixforge.gff3 --out report.html

  helixforge utils qc --gff3 helixforge.gff3 --genome genome.fa \
      --helixer-h5 predictions.h5 --out report.json --format json
```
