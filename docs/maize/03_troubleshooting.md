# 03 — Troubleshooting (the real gotchas)

Each entry is **symptom → cause → fix**. The most likely first-real-chunked-run
failure is the `aggregate` safety gate (last-but-one entry) — that is the check
working, not a bug.

---

### Mikado `serialise` crashes / pandas or SQLAlchemy errors; or `diamond` not found

**Cause.** Mikado pins pandas in its conda env, and `diamond` typically lives
**only** in that env. Cold runs (not resuming prior Mikado artifacts) need both.
Separately, a stale `helixforge` console-script shim inside `mikado-env` can
point at an old tree.

**Fix.**
- Run inside the `mikado` env (its pinned pandas), and pass `diamond` explicitly
  if it lives elsewhere:

  ```bash
  helixforge parallel tasks --plan plan.json -o tasks.txt \
    --genome maize.fa --helixer helixer.gff3 \
    --mikado-bin mikado --diamond-bin diamond
  ```
- If the `helixforge` command in `mikado-env` is a **foreign shim**, reinstall
  the package into that env (`pip install -e .`) or call the user shim directly
  (`~/.local/bin/helixforge`). `helixforge doctor` flags a foreign shim.

---

### CRAM RNA-seq input fails to decode / tries to fetch a remote reference

**Cause.** CRAM stores no sequence; it needs the reference. With no local
reference it falls back to a remote (ENA) fetch — slow and fragile mid-run.

**Fix.** Provide the reference offline: `--reference maize.fa` (or set
`REF_CACHE` / `REF_PATH`). Confirm with `helixforge doctor --reference maize.fa`
before launching. Never rely on the remote fetch.

---

### Run is pathologically slow, or fails writing the Mikado SQLite DB

**Cause.** `mikado serialise`'s SQLite DB insertion is single-threaded and dies
on a slow/over-quota network filesystem.

**Fix.** Put each chunk's outputs/work dir on **fast local/scratch** with free
space — point `tasks --output-dir` at scratch (the default per-chunk command
derives each chunk's work dir from its output prefix). Never put the Mikado DB on
network FS.

---

### Organellar (plastid/mito) models look wrong / CDS rejected

**Cause.** The default genetic code (table 1, standard nuclear) is wrong for
plastid/some mito scaffolds, so valid CDS fail the codon gate.

**Fix.** Set the genetic code per seqid on `tasks`:

```bash
helixforge parallel tasks --plan plan.json -o tasks.txt \
  --genome maize.fa --helixer helixer.gff3 \
  --transl-table-map organelle_codes.txt
```

with `chrPt=11` / `chrMt=1` lines (see `00_inputs_checklist.md` §5), or
`--transl-table 11` for a whole organellar run.

---

### `helixforge parallel aggregate` raises on a duplicate HFG id or missing/duplicated locus

**Cause — this is the safety check, not a defect.** `aggregate` fails closed
when genome-wide HFG uniqueness or no-gene-lost is violated. The usual root
cause on a first real chunked run is a **gene split across a chunk boundary**
(so two chunks each emit part of it / reuse an id).

**Fix.** **Report the message**, don't force past it. Then:
- Raise `--min-boundary-gap` on `parallel plan` so cuts fall only in wide
  inter-locus gaps (maize's wide intergenic spacing makes this safe), re-plan,
  re-run the affected chunks, and re-aggregate.
- Confirm you seeded the **same** `--id-map`/`--master-id-map` across the run so
  reserved id ranges stay disjoint.
- This is exactly what the single-chromosome dry run (`02` §3) is meant to
  surface in minutes.

---

### Inflated gene count / odd single-exon models in TE-rich regions

**Cause.** Expected — maize is ~85% TE and **there is no TE-overlap module
yet**, so Helixer-derived TE models pass through.

**Fix.** Post-filter the output against your **EDTA TE annotation** (e.g.
`bedtools intersect` of the tiered GFF3 vs the EDTA GFF3) to flag/remove
TE-overlapping models. Do not treat the inflated count as a tool error. See
`00_inputs_checklist.md` §6 and `04_caveats_known_limits.md`.

---

### 27 RNA-seq samples → an unwieldy command line of `--bam`/`--sj`/`--stringtie`

**Cause.** Every multi-file evidence flag is repeatable, but typing 27 of each
is error-prone.

**Fix.** Each flag also accepts a **comma-separated** value and a `-list`
**file-of-filenames** (FOFN) companion (`--bam-list`, `--sj-list`,
`--stringtie-list`); all three forms merge, paths are validated, and a FOFN of
bare names resolves against the list file's own directory. Build the lists once:

```bash
ls bamfiles.v2/*.bam       > bam.list
ls bamfiles.v2/*SJ.out.tab > sj.list
ls stringtie.v2/*.gtf      > stringtie.list

# standalone evidence scoring (the pipeline flags take --*-list the same way):
helixforge evidence --gff3 helixer.gff3 \
  --bam-list bam.list --sj-list sj.list --stringtie-list stringtie.list \
  --out evidence.tsv
```

`helixforge evidence` reports **two separate AED scores** (each in `[0, 1]`,
lower = better), never a fused one — a column is filled only when its evidence
was supplied (missing evidence is neutral, not a penalty):

- `rna_aed` (+ `rna_junction_ratio` / `rna_coverage_ratio` / `rna_boundary_ratio`)
  from `--bam`/`--sj`;
- `protein_aed` (+ `protein_struct_ratio` / `protein_cds_cov_ratio` /
  `protein_prot_cov_ratio`) from a reference proteome aligned with miniprot.

Add the protein axis with `--proteins` (give a **TE-filtered** proteome — TE
proteins align well and would hand TE models a deceptively low `protein_aed`),
and use `-j/--threads` to parallelize on one node (thread-parallel extraction
then process-parallel scoring; the output is identical to `-j 1`). A per-gene
table (best-supported transcript per gene) is written alongside `--out`:

```bash
helixforge evidence --gff3 helixer.gff3 \
  --bam-list bam.list --sj-list sj.list \
  --proteins swissprot.te_filtered.faa --genome genome.fa \
  -j 32 --out evidence.tsv --gene-out evidence.gene.tsv
```

> Re-runs can skip realigning by passing the saved miniprot alignment with
> `--miniprot-gff` instead of `--proteins`/`--genome`.

> `--stringtie` used to be the sample-list file itself; it is now repeatable
> individual GTFs, with the list-file moved to `--stringtie-list`. The old form
> still works for one release (with a deprecation warning) — switch to
> `--stringtie-list`.

---

Still stuck? Capture the artifacts listed in
[`04_caveats_known_limits.md`](04_caveats_known_limits.md) ("How to report
issues") and send them along.
