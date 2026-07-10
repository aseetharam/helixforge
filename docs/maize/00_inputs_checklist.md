# 00 — Inputs checklist (prepare these before HelixForge)

HelixForge does **not** align reads, call genes, or assemble transcripts. It
*reconciles* evidence you produce upstream. Prepare the following first; the
right column is the HelixForge flag each one plugs into (flags attach to
`helixforge parallel tasks` — see `02_end_to_end_chunked.md`).

| Input | Produced by (you run this) | HelixForge flag |
|---|---|---|
| Genome FASTA + `.fai` | `samtools faidx maize.fa` | `--genome` (FASTA), `suggest`/`plan` accept the `.fai` |
| Helixer GFF3 | Helixer (GPU) | `--helixer` (**required**) |
| Helixer prediction HDF5 | Helixer (GPU) | `--helixer-h5` (confidence prior) |
| StringTie GTF, **one per sample** | align reads → StringTie | `--stringtie` (**repeatable / comma**) or `--stringtie-list` (FOFN) |
| Protein GFF (optional) | miniprot | `--miniprot` |
| Protein FASTA (optional) | your proteome | `--protein-db` (DIAMOND homology) |
| Organellar genetic-code map (optional) | hand-written | `--transl-table` / `--transl-table-map` |

## 1. Genome FASTA + `.fai`

```bash
samtools faidx maize.fa        # writes maize.fa.fai
```

`parallel suggest` and `parallel plan` read scaffold sizes; either pass the
FASTA or the `.fai` index (the `.fai` is enough for `suggest` and is what the
walkthrough uses, since it is tiny).

## 2. Helixer run on maize (external, GPU — a prerequisite you run yourself)

Run Helixer with the land-plant model on the maize genome. It produces:

- the **Helixer GFF3** → `--helixer` (the gene set HelixForge anchors to and
  rescues silent genes from; **required**), and
- the **prediction HDF5** → `--helixer-h5` (the per-base confidence track that
  becomes the `helixer_support` / `helixer_locus_conf` external metric Mikado
  scores against).

HelixForge reads the native split-Helixer HDF5 pair directly; pass the
`*_predictions.h5` (it links its `*_input.h5` metadata automatically).

## 3. RNA-seq evidence → StringTie per sample

Align each RNA-seq sample (STAR/HISAT2), then assemble **per sample** with
StringTie, and pass each sample's GTF. `--stringtie` is repeatable, accepts a
comma-separated list, and has a `--stringtie-list` file-of-filenames (FOFN)
companion — all three forms merge, so with many samples a FOFN is easiest:

```text
# any of these (they merge):
--stringtie root.gtf --stringtie leaf.gtf --stringtie tassel.gtf
--stringtie root.gtf,leaf.gtf,tassel.gtf
--stringtie-list stringtie.list          # one GTF path per line
```

```bash
ls stringtie.v2/*.gtf > stringtie.list   # build the FOFN for 27 samples
```

Per-sample GTFs are **both** the isoform evidence and the **TRaCE voters** for
the canonical-isoform election. **More samples / tissues = better isoform
discovery and more confident canonical calls.** Do not merge all samples into
one GTF — that throws away the per-sample votes.

## 4. Optional protein evidence (miniprot)

A miniprot proteome→genome GFF (`--miniprot`) is the primary source of backstop
CDS: it rescues a coding model for Helixer-only genes that have no RNA-seq
support. Without it, those backstop genes stay CDS-less and are flagged
`HELIXER_ONLY`. A protein FASTA (`--protein-db`) drives DIAMOND homology scoring.

## 5. Organellar scaffolds (if your assembly includes them)

If plastid/mito scaffolds are in the genome, set the genetic code so CDS
validation does not reject them:

- `--transl-table 11` for a whole run that is bacterial/plastid, **or**
- `--transl-table-map codes.txt` for per-seqid overrides, one `seqid=table` per
  line, e.g.:

  ```text
  chrMt=1
  chrPt=11
  ```

  This overrides the default `--transl-table` only for the named seqids.

## 6. EDTA TE annotation — honest note (read this)

**The current pipeline does NOT consume the TE annotation. There is no
TE-overlap module yet.** Consequences for maize specifically:

- Maize is ~85% transposable elements. Helixer (like all *ab initio* callers)
  emits gene-like models over some TEs, so the **raw HelixForge gene set will
  contain TE-derived models** and the **gene count will be inflated**. This is
  expected behaviour, **not** a tool defect.
- Keep your **EDTA TE annotation to post-filter / interpret** the output: flag
  or remove output genes whose loci overlap annotated TEs (e.g. `bedtools
  intersect` between the tiered GFF3 and the EDTA GFF3). Many TE models land in
  Tier 3 (Helixer-only backstop), but do not assume tier alone separates them —
  use the EDTA overlap.
- Expect the inflated count until native TE flagging lands (see
  `04_caveats_known_limits.md`).

Next: [`01_quickstart.md`](01_quickstart.md).
