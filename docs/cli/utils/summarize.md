# `helixforge utils summarize`

Annotation statistics table.

## Options

| Option | Description |
|---|---|
| `--gff3` | Input GFF3. (required) |
| `--out` | Output file (default: stdout). |
| `--format` | Output format [default: tsv]. (choices: tsv/json/markdown; default: `tsv`) |

## Examples

```bash
  helixforge utils summarize --gff3 helixforge.gff3

  helixforge utils summarize --gff3 helixforge.gff3 --out stats.json --format json
```
