# `helixforge stats`

Inspect (read-only): before/after table of a Helixer vs HelixForge GFF3.

Works directly off two GFF3 paths (no rerun, no Mikado): it summarises the
Helixer input and the reconciled output and tabulates the change. Pass
``--helixer-h5`` to add the mean-Helixer-support row. Use this to show what
reconciliation did; use ``benchmark`` for external accuracy metrics against a
reference / BUSCO / OMArk.

Output (--out Markdown), one row per metric, columns:
metric, helixer, helixforge, delta. Metric rows include gene and transcript
counts, coding genes, isoforms/gene (mean/max), mono/multi-exon genes, CDS
length (mean/median), % complete ORFs, mean AED, and mean Helixer support.

## Options

| Option | Description |
|---|---|
| `--helixer` | Helixer GFF3 (before). (required) |
| `--reconciled` | HelixForge GFF3 (after). (required) |
| `--genome` | Genome FASTA (optional). |
| `--helixer-h5` | Helixer HDF5 → populates the Helixer-support before/after column. |
| `--out` | Markdown summary output path. (default: `before_after.md`) |

## Examples

```bash
  # before/after on two GFF3s
  helixforge stats --helixer helixer.gff3 --reconciled helixforge.gff3 \
      --out before_after.md

  # add the Helixer-support column (needs the HDF5)
  helixforge stats --helixer helixer.gff3 --reconciled helixforge.gff3 \
      --helixer-h5 helixer.h5 --out before_after.md
```
