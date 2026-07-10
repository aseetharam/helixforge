# `helixforge parallel suggest`

Recommend chunk count + per-chunk resources (heuristic; prints trade-offs).

Step 0 of the chunked path. From the genome size/scaffold profile and your
node spec (``--cores-per-node``, ``--mem-per-node``, ``--max-array-size``,
``--walltime-cap``) it suggests a granularity and per-chunk resources to feed
into ``parallel plan`` / ``run --scatter``. Heuristic only, it reads sizes,
runs nothing, and modifies nothing.

Output (printed recommendation): target_chunks, target_loci_per_chunk,
min_boundary_gap, and per-chunk procs / threads / mem_gb / walltime_min,
with a rationale explaining the trade-offs.

## Options

| Option | Description |
|---|---|
| `--genome` | Genome FASTA or .fai index. (required) |
| `--cores-per-node` | (default: `16`) |
| `--mem-per-node` | RAM per node, GB. (default: `64`) |
| `--max-array-size` | (default: `1000`) |
| `--walltime-cap` | Walltime cap, hours. (default: `24`) |

## Examples

```bash
  # default node profile
  helixforge parallel suggest --genome maize.fa

  # tune to a specific node + array/walltime caps
  helixforge parallel suggest --genome maize.fa --cores-per-node 32 \
      --mem-per-node 128 --max-array-size 500 --walltime-cap 12
```
