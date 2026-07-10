# `helixforge utils align`

Run miniprot protein-to-genome alignment.

## Options

| Option | Description |
|---|---|
| `--genome` | Reference genome FASTA. (required) |
| `--proteins` | Protein FASTA. (required) |
| `--out` | Output GFF. (required) |
| `--threads` | Number of threads [default: 4]. (default: `4`) |
| `--miniprot-bin` | miniprot binary [default: miniprot]. (default: `miniprot`) |

## Examples

```bash
  helixforge utils align --genome genome.fa --proteins uniprot.fa --out miniprot.gff

  helixforge utils align --genome genome.fa --proteins uniprot.fa \
      --out miniprot.gff --threads 16 --miniprot-bin /opt/bin/miniprot
```
