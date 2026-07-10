# `helixforge parallel plan`

Partition the genome (v1 strategies) + reserve disjoint HFG ranges → plan.json.

Step 1 of the chunked path (plan → tasks → run → aggregate). Cuts the genome
only in inter-locus gaps wider than ``--min-boundary-gap`` so no gene — and no
Mikado merge — is ever split, then reserves each chunk a disjoint, contiguous
HFG number range (seeded from ``--id-map`` for run-stable ids). The plan is
consumed by ``parallel tasks``.

Output: plan.json with top-level {strategy, min_boundary_gap, flank,
total_loci} and a "chunks" list; each chunk has: chunk_id, regions,
num_loci, locus_ids, id_base, id_range, novel_base, novel_range,
est_resources.

## Options

| Option | Description |
|---|---|
| `--genome` | Genome FASTA or .fai index (for scaffold sizes). (required) |
| `--helixer`, `--gff` | Helixer GFF3 (master loci; drives gene-aware cuts + id ranges). (required) |
| `--strategy` | Chunking strategy (v1 vocabulary). (choices: scaffold/size/genes/adaptive; default: `scaffold`) |
| `--chunk-size` | Bases per chunk for 'size'; loci per chunk for 'genes'. |
| `--min-chunk-size` | Smallest chunk (bp) for the size-based strategies. (default: `100000`) |
| `--max-chunk-size` | Split a scaffold longer than this (bp) under 'scaffold'. |
| `--target-chunks` | Target chunk count for 'adaptive' (from `parallel suggest`). |
| `--min-boundary-gap` | Only cut in inter-locus gaps at least this wide (bp), so a gene is never split. Default max(--flank, 1000). |
| `--flank` | (default: `200`) |
| `--id-map` | Master id_map.json to seed stable HFG ranges from. |
| `-o`, `--output` | Output plan path. (default: `plan.json`) |

## Examples

```bash
  # one chunk per scaffold (default), stable ids seeded from a prior run
  helixforge parallel plan --genome maize.fa --gff helixer.gff3 \
      --id-map master_id_map.json -o plan.json

  # ~2000 loci per chunk
  helixforge parallel plan --genome maize.fa --gff helixer.gff3 \
      --strategy genes --chunk-size 2000 -o plan.json
```
