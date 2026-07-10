# `helixforge utils filter`

Tiered output filtering by confidence/evidence/biotype.

## Options

| Option | Description |
|---|---|
| `--gff3` | Input GFF3. (required) |
| `--out` | Output filtered GFF3. (required) |
| `--preset` | Filter preset [default: custom]. (choices: high_confidence/publication_ready/custom; default: `custom`) |
| `--min-confidence` | Minimum confidence score. |
| `--min-tier` | Minimum tier (inclusive). |
| `--max-tier` | Maximum tier (inclusive). |
| `--exclude-biotype` | Biotype(s) to exclude, repeatable. (repeatable) |

## Examples

```bash
  helixforge utils filter --gff3 helixforge.gff3 --out tier1.gff3 --max-tier 1

  helixforge utils filter --gff3 helixforge.gff3 --out pub.gff3 \
      --preset publication_ready --exclude-biotype transposable_element
```
