# HelixForge Command-line reference

> Auto-generated from the `helixforge` click app by `scripts/gen_cli_reference.py`. Do not edit by hand, run the script (or let the pre-commit/CI drift gate regenerate it). Each command's description, options, and examples come straight from its docstring and epilog, so this reference cannot drift from the code.

HelixForge v3: isoform-aware refinement of Helixer annotations.

| Command | Description |
|---|---|
| [`helixforge confidence`](confidence.md) | Inspect (read-only): score genes against the Helixer HDF5 confidence track. |
| [`helixforge doctor`](doctor.md) | Preflight: validate inputs + resolve external tools before a run. |
| [`helixforge evidence`](evidence.md) | Inspect (read-only): score any GFF3 against RNA-seq + protein evidence. |
| [`helixforge parallel aggregate`](parallel/aggregate.md) | Collect per-chunk outputs by pattern → one annotation; verify + fail closed. |
| [`helixforge parallel example-sbatch`](parallel/example-sbatch.md) | Write a copy-paste SBATCH wrapper that runs a ``tasks`` file on one node. |
| [`helixforge parallel plan`](parallel/plan.md) | Partition the genome (v1 strategies) + reserve disjoint HFG ranges → plan.json. |
| [`helixforge parallel suggest`](parallel/suggest.md) | Recommend chunk count + per-chunk resources (heuristic; prints trade-offs). |
| [`helixforge parallel tasks`](parallel/tasks.md) | Expand a command template over the plan → an executor-agnostic task file. |
| [`helixforge reconcile`](reconcile.md) | Reconcile Helixer models with evidence → tiered GFF3 + per-gene report. |
| [`helixforge stats`](stats.md) | Inspect (read-only): before/after table of a Helixer vs HelixForge GFF3. |
| [`helixforge utils align`](utils/align.md) | Run miniprot protein-to-genome alignment. |
| [`helixforge utils convert`](utils/convert.md) | GFF3 <-> GTF format conversion. |
| [`helixforge utils extract-proteins`](utils/extract-proteins.md) | Translate CDS from GFF3 + genome to protein FASTA. |
| [`helixforge utils fetch-db`](utils/fetch-db.md) | Download, decompress, and Diamond-format reference protein databases. |
| [`helixforge utils filter`](utils/filter.md) | Tiered output filtering by confidence/evidence/biotype. |
| [`helixforge utils qc`](utils/qc.md) | Genome-wide QC report (HTML/JSON/TSV). |
| [`helixforge utils summarize`](utils/summarize.md) | Annotation statistics table. |
| [`helixforge viz`](viz.md) | Render per-locus plots, interactive pages, or browser tracks for a run. |
