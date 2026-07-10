# `helixforge utils convert`

GFF3 <-> GTF format conversion.

## Options

| Option | Description |
|---|---|
| `--input` | Input file. (required) |
| `--out` | Output file. (required) |
| `--from-format` | Source format: gff3 \| gtf (inferred from extension). |
| `--to-format` | Target format: gff3 \| gtf (inferred from extension). |

## Examples

```bash
  helixforge utils convert --input genes.gff3 --out genes.gtf

  helixforge utils convert --input genes.gtf --out genes.gff3 \
      --from-format gtf --to-format gff3
```
