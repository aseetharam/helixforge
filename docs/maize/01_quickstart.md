# 01 — Quickstart (the chunked happy path)

This is the condensed end-to-end path, assuming the inputs from
[`00_inputs_checklist.md`](00_inputs_checklist.md) are ready. Replace the
placeholder paths. The full explanation of each step — and the **required**
single-chromosome dry run you should do *first* — is in
[`02_end_to_end_chunked.md`](02_end_to_end_chunked.md). **Do the dry run before
you trust the full-genome path.**

The workflow is four transparent steps: **suggest → plan → tasks → (run the task
file) → aggregate**. `plan` partitions the genome; `tasks` expands a command
template into a plain task file you can read and run with *any* executor (GNU
parallel, a Slurm array, xargs, HyperShell); `aggregate` merges the per-chunk
outputs back into one verified annotation.

```bash
# (one-time) build a file-of-filenames for the per-sample StringTie GTFs.
ls stringtie.v2/*.gtf > stringtie.list

# 0) Preflight: tools resolve, inputs are well-formed.
helixforge doctor \
  --genome maize.fa --helixer helixer.gff3 --helixer-h5 helixer_predictions.h5 \
  --stringtie root.gtf --stringtie leaf.gtf

# 1) Recommend chunk count + resources for a ~2.1 Gb genome.
#    Its printed --target-chunks / --min-boundary-gap drop straight into step 2.
helixforge parallel suggest --genome maize.fa.fai

# 2) Partition the genome + reserve disjoint per-chunk HFG id ranges.
#    plan only needs genome + helixer (NO pipeline options here).
#    --strategy genes groups ~2000 master loci per chunk; scaffold/size/adaptive
#    are the other v1 strategies.
helixforge parallel plan \
  --genome maize.fa --gff helixer.gff3 \
  --strategy genes --chunk-size 2000 \
  --min-boundary-gap 50000 \
  --id-map master_id_map.json -o plan.json

# 3) Expand a command template into a task file. ALL pipeline inputs/flags attach
#    HERE — they are baked into the default per-chunk `reconcile` template, with
#    {region}/{id_start}/{output_dir} filled per chunk. --stringtie is repeatable /
#    comma-separated; for many samples pass the FOFN via --stringtie-list.
helixforge parallel tasks \
  --plan plan.json -o tasks.txt --output-dir chunks \
  --genome maize.fa --helixer helixer.gff3 --helixer-h5 helixer_predictions.h5 \
  --stringtie-list stringtie.list \
  --trace-primary

# 4) Run the task file with your executor of choice (one line per chunk):
parallel -j 16 < tasks.txt                  # local / GNU parallel
#   hs cluster tasks.txt --num-tasks 16       # HyperShell (SLURM/PBS/SGE/local)
#   sbatch --array=0-$(($(wc -l < tasks.txt)-1)) run_array.sbatch   # Slurm array

# 5) After all chunks finish: merge them → one genome-wide annotation,
#    verifying global HFG (gene + transcript) uniqueness + no-gene-lost.
helixforge parallel aggregate \
  --input-dir chunks --pattern '*.gff3' \
  -o maize_helixforge.gff3 --master-id-map master_id_map.json
```

`aggregate` prints the genome-wide `num_genes`, `num_loci`, and per-`tier` /
per-`origin` tallies, and writes `maize_helixforge.gff3` plus
`maize_helixforge.tier{1,2,3}.gff3` and `maize_helixforge.report.tsv`. See
[`02_end_to_end_chunked.md`](02_end_to_end_chunked.md) §5 for what the outputs
mean and the EDTA/TE post-filter note.

**Full control when you want it:** pass `--command '<template>'` to `tasks` to
replace the default per-chunk command wholesale (v1 model), and `--wrapper
run_chunk.sh --wrapper-setup 'module load helixforge'` to wrap every task in a
one-time environment setup. See [`02_end_to_end_chunked.md`](02_end_to_end_chunked.md).
