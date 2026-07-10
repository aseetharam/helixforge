# `helixforge utils extract-proteins`

Translate CDS from GFF3 + genome to protein FASTA.

## Options

| Option | Description |
|---|---|
| `--gff3` | Input GFF3. (required) |
| `--genome` | Reference genome FASTA. (required) |
| `--out` | Output protein FASTA. (required) |
| `--longest-only`, `--all-isoforms` | Keep only longest isoform per gene. (default: `True`) |
| `--transl-table` | NCBI translation table [default: 1]. (default: `1`) |
| `--min-length` | Minimum protein length in aa [default: 30]. (default: `30`) |

## Examples

```bash
  helixforge utils extract-proteins --gff3 genes.gff3 --genome genome.fa --out proteins.fa

  helixforge utils extract-proteins --gff3 genes.gff3 --genome genome.fa \
      --out proteins.fa --all-isoforms --min-length 50 --transl-table 11
```
