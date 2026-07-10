# 04 — Caveats & known limits (read before you start)

You are the **first external user running HelixForge on maize** — the first
large genome it has seen. Honest, up front:

## 1. Validated only on Arabidopsis so far

HelixForge's accuracy/completeness has been validated on **Arabidopsis (TAIR10)**
against Araport11 (see `docs/BENCHMARK.md`). **Maize (~2.1 Gb) is the first
large-genome run.** Treat headline behaviour as expected-but-unverified at maize
scale, and lean on the dry run + the `aggregate` safety gate.

## 2. The genome-wide `aggregate` merge has not run on a real multi-chunk genome

The chunked **`aggregate`** positive path is exercised in unit tests but has
**not yet run on a real multi-chunk genome**. This is the single biggest reason
the walkthrough makes you do the **single-chromosome dry run first**
([`02_end_to_end_chunked.md`](02_end_to_end_chunked.md) §3): validate the merge
(does it merge, are IDs globally unique, are gene counts conserved) on two
chunks in minutes before committing a multi-day full-genome run. `aggregate`
fails closed on any violation — trust that, and report the message.

## 3. No TE-overlap handling yet

The pipeline **does not consume a TE annotation**; there is no TE-overlap
module. Maize is ~85% TE, so the gene set **will include TE-derived models** and
the gene count **will be inflated**. Interpret/post-filter with your **EDTA**
annotation (`00_inputs_checklist.md` §6, `03_troubleshooting.md`). This is
expected behaviour, not a defect.

## 4. Accuracy tradeoffs carry over from the Arabidopsis benchmark

From `docs/BENCHMARK.md` (Arabidopsis), and expected to carry over:

- **Precision is strong; reconciliation trades a little completeness for it.**
  compleasm Complete drops only slightly under reconciliation (99.09% → 97.38%)
  while gene-level F1 rises (67.0 → 74.7). Expect the same precision/recall
  shape on maize: a modest completeness cost for cleaner, evidence-backed models.
- **Isoform recall is bounded by RNA-seq depth.** HelixForge leads on
  **precision of the isoforms it predicts**, but full recovery of the true
  alternative-splicing catalog is bounded by how deep/broad your short-read
  RNA-seq is. **More tissues/samples → more of the catalog** (each per-sample
  StringTie GTF is also a TRaCE voter). Do not expect isoforms with no
  short-read evidence.

## 5. Install / packaging

HelixForge is **not** a published PyPI package. Install from the GitHub repo plus
the conda/Apptainer environment, into the same env as Mikado/TransDecoder/
DIAMOND (see `docs/EXTERNAL_TOOLS.md` for pinned versions, and
`03_troubleshooting.md` for the `mikado-env` console-script gotcha).

## How to report issues back

When something fails, capture and send:

1. The **exact failing command** (full argv) and its stderr tail.
2. `helixforge doctor` output (with the same inputs).
3. The **`aggregate` error message** verbatim if it failed closed.
4. The relevant `plan.json` and the per-chunk report/id_map JSON for the
   affected chunk(s).
5. Tool versions (the `doctor` preflight prints them) and which conda env /
   container you ran in.

These five items are usually enough to reproduce and diagnose without the full
genome.
