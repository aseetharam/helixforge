---
title: Known limitations
description: What the HelixForge benchmarks do and do not support, stated plainly.
---

Stated plainly, because they bound what the [benchmark numbers](/benchmarks/) support.

- **Novel content is junctions, not genes.** HelixForge calls coding loci that NAM lacks,
  and most of those are not RNA-supported. Treat them as candidates, not discoveries.
- **ORF completeness trails the reference.** Roughly **86%** of HelixForge ORFs are complete
  versus **93%** for NAM, driven by 5′-partial models. This is not parity.
- **Improvement is net, not universal.** Reconciliation adds far more exact transcript
  matches than it degrades, but a small number of models do get worse.
- **Unsupported novel junctions are ambiguous.** They are either unexpressed in the sampled
  tissues or genuine over-calls. The available evidence cannot separate the two.
- **Isoform precision is the open axis.** The extra isoforms HelixForge emits are largely
  plausible, but tightening isoform admission would raise chain-level precision.
- **Requires protein evidence.** Without a protein database (`--protein-db`) the Mikado
  chain does not run and the output is backstop-only, with no isoform discovery.
