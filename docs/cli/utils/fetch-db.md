# `helixforge utils fetch-db`

Download, decompress, and Diamond-format reference protein databases.

## Options

| Option | Description |
|---|---|
| `--db` | Database name or local FASTA path. |
| `--cache-dir` | Cache directory [default: ~/.helixforge/databases]. |
| `--force` | Re-download even if cached. (default: `False`) |
| `--list` | List available databases and exit. (default: `False`) |
| `--format-tool` | Formatting tool [default: diamond]. (default: `diamond`) |

## Examples

```bash
  helixforge utils fetch-db --list

  helixforge utils fetch-db --db uniprot_sprot --cache-dir /scratch/dbs

  helixforge utils fetch-db --db /local/proteins.fa --format-tool diamond --force
```
