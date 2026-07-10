# 02 — End-to-end chunked run (full walkthrough)

This is the full maize path, in order. Each step shows the real command, what it
produces, and how to check it. Inputs are assumed ready
([`00_inputs_checklist.md`](00_inputs_checklist.md)). The complete flag list for
any command is in the per-command reference under [`../cli/`](../cli/README.md).

The chunked path is four transparent steps: **suggest** the granularity →
**plan** the partition → **tasks** expand a command template into a plain task
file → **run that task file** with any executor → **aggregate** the chunks back
into one annotation. Pipeline inputs/flags attach to `tasks` (they are baked into
the default per-chunk command); `plan` only needs the genome + Helixer.

The model is deliberately v1-style and executor-agnostic: `tasks` writes one
command per chunk to a file *you can read*, and you run it however you like — GNU
parallel, a Slurm array, `xargs`, or HyperShell. Nothing is hidden inside a
generated submit script.

---

## Step 1 — `helixforge doctor` (environment + input readiness)

```bash
helixforge doctor \
  --genome maize.fa \
  --helixer helixer.gff3 \
  --helixer-h5 helixer_predictions.h5 \
  --stringtie root.gtf --stringtie leaf.gtf
```

`doctor` resolves the external tools (Mikado/DIAMOND/etc.) with versions and
runs input-integrity checks (format, seqid concordance, BAM/CRAM index
presence). Fix anything it flags before launching a multi-day run.

**If any RNA-seq input is CRAM**, add the decode reference so doctor can confirm
an offline reference exists (never rely on a remote ENA fetch mid-run):

```bash
helixforge doctor --genome maize.fa --helixer helixer.gff3 --reference maize.fa
```

`doctor` also flags a **foreign `helixforge` console-script shim** (the
`mikado-env` gotcha) — see [`03_troubleshooting.md`](03_troubleshooting.md).

---

## Step 2 — `helixforge parallel suggest` (chunk count + resources)

```bash
helixforge parallel suggest --genome maize.fa.fai
```

Prints a recommended chunk count and per-chunk resources for the genome size.
Tune the recommendation to your cluster with:

```bash
helixforge parallel suggest \
  --genome maize.fa.fai \
  --cores-per-node 36 --mem-per-node 180 \
  --max-array-size 1000 --walltime-cap 48
```

It prints a ready-to-paste `parallel plan` line: the recommended
`--target-chunks` and `--min-boundary-gap` are **directly usable** as `plan`
flags (Step 4). For a ~2.1 Gb genome expect it to recommend many chunks; pick a
strategy/granularity in Step 4 that keeps per-chunk walltime under your cap.

---

## Step 3 — Single-chromosome dry run **FIRST** (required de-risking step)

**Do not skip this.** Run the *entire* plan → tasks → run → aggregate path on
**one maize chromosome forced into two chunks** before the full genome. Reason:
the genome-wide `aggregate` merge has been exercised in unit tests but **not yet
on a real multi-chunk genome**. Validating it small — does it merge, are HFG IDs
globally unique, are gene counts conserved — catches problems in minutes
instead of after a multi-day full run.

Subset the inputs to one chromosome (plain shell, not HelixForge):

```bash
# pick chr1; build a chr1-only genome + index, and filter Helixer to chr1
samtools faidx maize.fa chr1 > chr1.fa
samtools faidx chr1.fa
awk -F'\t' '$1=="chr1" || /^#/' helixer.gff3 > helixer.chr1.gff3
```

Plan it into **two** chunks so `aggregate` actually has a multi-chunk merge to
verify, then run the chunks **locally** (no Slurm needed for the dry run):

```bash
# plan: genome + helixer only, forced to 2 chunks (adaptive strategy)
helixforge parallel plan \
  --genome chr1.fa --gff helixer.chr1.gff3 \
  --strategy adaptive --target-chunks 2 \
  --min-boundary-gap 50000 \
  -o plan.chr1.json

# tasks: expand the per-chunk reconcile template into a task file
# (pipeline flags attach here, baked into the default command)
helixforge parallel tasks \
  --plan plan.chr1.json -o tasks.chr1.txt --output-dir chunks \
  --genome chr1.fa --helixer helixer.chr1.gff3 --helixer-h5 helixer_predictions.h5 \
  --stringtie root.gtf --stringtie leaf.gtf \
  --trace-primary

# run both chunks locally (one command per line) and wait
parallel -j 2 < tasks.chr1.txt

# aggregate the chunk outputs and confirm the safety checks pass
helixforge parallel aggregate \
  --input-dir chunks --pattern '*.gff3' \
  -o chr1_dryrun.gff3 --master-id-map master_id_map.chr1.json
```

**Pass criteria:** `aggregate` exits 0 (it raises on any duplicate HFG gene/
transcript id or missing/duplicated locus), prints `num_genes` / `num_loci` /
tier / origin tallies, and the gene count equals the sum across the two chunk
reports. If it *fails closed*, that is the safety check doing its job — read the
message and fix the cause (usually a boundary-split gene → raise
`--min-boundary-gap`) rather than forcing past it.

> `--input-dir chunks --pattern '*.gff3'` discovers the per-chunk outputs the
> `tasks` step wrote under `--output-dir chunks` (the default per-chunk command
> sets `--output-prefix chunks/<chunk_id>`). The `*.gff3` glob matches both the
> main and the tier GFF3s; aggregate collapses them to one prefix per chunk.

---

## Step 4 — Full genome

### 4a. Plan (genome + Helixer only — **no** pipeline options here)

```bash
helixforge parallel plan \
  --genome maize.fa --gff helixer.gff3 \
  --strategy genes --chunk-size 2000 \
  --min-boundary-gap 50000 \
  --flank 2000 \
  --id-map master_id_map.json \
  -o plan.json
```

- `--strategy` — v1's chunking vocabulary: `scaffold` (one chunk per scaffold,
  split long scaffolds with `--max-chunk-size`), `size` (`--chunk-size` bp
  windows), `genes` (`--chunk-size` master loci per chunk, shown here), or
  `adaptive` (`--target-chunks`, from `suggest`). Every strategy cuts **only in
  inter-locus gaps**, so a gene is never split.
- `--min-boundary-gap` — only cut in inter-locus gaps **at least this wide**, so
  a gene (or a Mikado merge) is never split across a chunk boundary. **Set it
  generously for maize** (wide intergenic gaps make this safe); too small and
  genes split at boundaries (the `aggregate` merge then fails closed).
- `--id-map master_id_map.json` — seed disjoint, **stable** HFG id ranges from a
  prior run's master map, so reruns keep their gene IDs. Omit on the first run.
- `-o plan.json` — the partition + per-chunk reserved id ranges (incl. the
  `id_start` each chunk's reconcile uses).

### 4b. Tasks → a task file (**all pipeline flags attach here**)

```bash
# With many tissues, collect the per-sample GTFs into a file-of-filenames.
# --stringtie is also repeatable / comma-separated; all three forms merge.
ls stringtie.v2/*.gtf > stringtie.list

helixforge parallel tasks \
  --plan plan.json \
  -o tasks.txt --output-dir chunks \
  --helixforge-bin helixforge \
  --genome maize.fa \
  --helixer helixer.gff3 \
  --helixer-h5 helixer_predictions.h5 \
  --stringtie-list stringtie.list \
  --miniprot proteins.miniprot.gff \
  --protein-db proteins.fasta \
  --trace-primary
```

`tasks` prints the **resolved command template** and writes one command per chunk
to `tasks.txt`. By default that command is `helixforge reconcile` built from the
inputs you attached, with `{region}`, `{id_start}` and `{output_dir}/{chunk_id}`
filled per chunk — so it works out of the box. `--trace-primary` elects the
canonical isoform per gene by ranked-choice voting over the per-sample
assemblies. See [`../cli/tasks.md`](../cli/parallel/tasks.md) for the full set.

**Default vs custom `--command`.** To control the per-chunk command exactly (v1
model), pass your own template; the placeholders are `{chunk_id} {region} {seqid}
{start} {end} {start_0} {end_0} {size} {id_start} {novel_start} {output_dir}`:

```bash
helixforge parallel tasks --plan plan.json -o tasks.txt --output-dir chunks \
  --genome maize.fa --helixer helixer.gff3 \
  --command 'helixforge reconcile --genome maize.fa --helixer helixer.gff3 \
             --helixer-h5 helixer_predictions.h5 --stringtie-list stringtie.list \
             --region {region} --id-base {id_start} --novel-base {novel_start} \
             --output-prefix {output_dir}/{chunk_id} \
             --id-map-path {output_dir}/{chunk_id}.id_map.json --trace-primary'
```

Add `--wrapper run_chunk.sh --wrapper-setup 'module load helixforge'` (repeatable)
to wrap every task in a one-time environment setup, and `--include-logging` to
redirect each task's stdout/stderr to `chunks/logs/<chunk_id>.log`.

### 4c. Run the task file

Pick whichever executor your site uses — `tasks.txt` is one command per chunk:

```bash
# Slurm array (one task per line):
sbatch --array=0-$(($(wc -l < tasks.txt)-1)) --wrap \
  "sed -n \"\$((SLURM_ARRAY_TASK_ID+1))p\" tasks.txt | bash"
# or GNU parallel on a big node:
parallel -j 32 < tasks.txt
# or HyperShell (scheduler-agnostic):
hs cluster tasks.txt --num-tasks 32
```

Each task runs reconcile on its chunk and writes that chunk's tiered GFF3 +
report + id_map under `chunks/<chunk_id>`.

### 4d. Aggregate (after **all** tasks finish)

```bash
helixforge parallel aggregate \
  --input-dir chunks --pattern '*.gff3' \
  -o maize_helixforge.gff3 \
  --master-id-map master_id_map.json
```

`aggregate` collects the per-chunk GFF3s by pattern, merges them into one
genome-wide annotation, and **verifies genome-wide HFG gene + transcript id
uniqueness and no-gene-lost across the whole genome**, failing closed on any
violation. (You can instead pass explicit `--chunk-output <prefix>` flags, one
per chunk, if you did not use a single output directory.)

---

## Step 5 — Outputs

`aggregate --out-prefix maize_helixforge` writes:

| File | Contents |
|---|---|
| `maize_helixforge.gff3` | **All** genes, all tiers, isoform-aware. |
| `maize_helixforge.tier1.gff3` | **Tier 1** — evidence-backed (RNA-seq/protein homology). |
| `maize_helixforge.tier2.gff3` | **Tier 2** — partial / weaker evidence. |
| `maize_helixforge.tier3.gff3` | **Tier 3** — Helixer-only backstop (rescued silent genes; `HELIXER_ONLY`). |
| `maize_helixforge.report.tsv` | Per-gene report: tier, origin, QC flags, AS events. |
| `maize_helixforge.id_map.json` | Merged master HFG id map (also folded into `--master-id-map`). |

**Tiers** rank evidence strength: Tier 1 is what you trust most (junction- and
homology-supported); Tier 3 is the completeness backstop carried over from
Helixer with no independent evidence. `aggregate` also prints the genome-wide
`num_genes`, `num_loci`, and per-`tier` / per-`origin` tallies to stdout — record
these as your run summary.

**Interpreting maize gene counts:** the count will be **inflated by TE-derived
models** (no TE module yet). Post-filter the output against your **EDTA TE
annotation** to flag/remove TE-overlapping models — see the honest note in
[`00_inputs_checklist.md`](00_inputs_checklist.md) §6 and the caveats in
[`04_caveats_known_limits.md`](04_caveats_known_limits.md).

Trouble? → [`03_troubleshooting.md`](03_troubleshooting.md).
