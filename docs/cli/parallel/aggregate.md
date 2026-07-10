# `helixforge parallel aggregate`

Collect per-chunk outputs by pattern → one annotation; verify + fail closed.

Step 4 of the chunked path. Merges the per-chunk outputs (by ``--input-dir`` +
``--pattern``, or explicit ``--chunk-output`` prefixes) into one genome-wide
annotation and folds the chunk id_maps into ``--master-id-map``. The merge is
verifying: it asserts globally-unique gene AND transcript ids and that every
Helixer locus is covered exactly once, raising on any violation.

Outputs (under the resolved prefix):
- <prefix>.gff3 (+ .tier1/2/3.gff3) and <prefix>.report.tsv — the merged
  annotation + per-gene report (same columns as ``reconcile``)
- the master id_map.json (created if absent)
- a printed summary: num_genes, num_loci, tier_counts, origin_counts

## Options

| Option | Description |
|---|---|
| `--input-dir` | Directory of per-chunk outputs (use with --pattern). |
| `--pattern` | Glob for the per-chunk GFF3s under --input-dir. (default: `*.gff3`) |
| `--chunk-output` | Explicit per-chunk output prefix (repeatable; alternative to --input-dir). (repeatable) |
| `-o`, `--output` | Merged GFF3 path, e.g. combined.gff3 ('.gff3' stripped for the prefix). |
| `--out-prefix` | Output prefix (alternative to -o/--output). |
| `--master-id-map` | Master id_map.json to fold chunk maps into (created if absent). |

## Examples

```bash
  helixforge parallel aggregate --input-dir chunks/ --pattern '*.gff3' \
      -o maize_helixforge.gff3 --master-id-map master_id_map.json
```
