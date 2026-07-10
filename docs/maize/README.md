# HelixForge on maize — chunked end-to-end handoff

HelixForge refines [Helixer](https://github.com/weberlab-hhu/Helixer) gene
annotations into an **isoform-aware, tiered GFF3** by reconciling RNA-seq and
protein evidence through [Mikado](https://github.com/EI-CoreBioinformatics/mikado),
anchored to the Helixer gene set (stable gene IDs, no Helixer locus silently
lost). Isoforms come only from evidence reconciliation on a splice graph; the
canonical/primary isoform per gene is elected by ranked-choice voting (TRaCE)
over the per-sample assemblies. Output genes are tiered by evidence strength
(Tier 1 evidence-backed → Tier 3 Helixer-only backstop).

> **Status:** HelixForge is **validated on Arabidopsis (TAIR10)**. **Maize is the
> first large-genome (~2.1 Gb) test.** These docs cover the *only* path that
> scales to maize — **scatter-gather (chunking)** via `helixforge parallel` plus a
> Slurm array. Read `04_caveats_known_limits.md` before you start: some of the
> chunked path (the genome-wide `aggregate` merge) has not yet run on a real
> multi-chunk genome, which is exactly why the walkthrough makes you do a
> single-chromosome dry run first.

## Document index

| File | What it covers |
|---|---|
| [`00_inputs_checklist.md`](00_inputs_checklist.md) | Everything to prepare **before** HelixForge, and how each input plugs in (incl. the EDTA/TE honest note). |
| [`01_quickstart.md`](01_quickstart.md) | The condensed chunked happy path end to end (suggest → plan → tasks → sbatch → aggregate). |
| [`02_end_to_end_chunked.md`](02_end_to_end_chunked.md) | The full walkthrough: doctor, suggest, the **required single-chromosome dry run**, the full-genome run, and reading the outputs. |
| [`03_troubleshooting.md`](03_troubleshooting.md) | The real gotchas (symptom → cause → fix). |
| [`04_caveats_known_limits.md`](04_caveats_known_limits.md) | Honest caveats and how to report issues back. |

## Complete option reference

Every subcommand and every flag is listed, auto-generated from the live CLI, in
the per-command reference under [`../cli/`](../cli/README.md) (one page per
command). These maize docs show the one happy path and the gotchas; those pages
are the authoritative source for anything not shown here. (They are regenerated
by `python scripts/gen_cli_reference.py`; a CI drift gate keeps them in sync with
the code, so they never lie about a flag.)

## Install

HelixForge is installed from the GitHub repo plus the conda/Apptainer
environment (it is **not** a published PyPI package). Install the package
editable into the same environment that holds Mikado/TransDecoder/DIAMOND, or
use the Apptainer image. See `docs/EXTERNAL_TOOLS.md` for the pinned tool
versions and `03_troubleshooting.md` for the `mikado-env` console-script gotcha.
