# HelixForge

**Evidence-based refinement of deep-learning gene annotations.**

[Helixer](https://github.com/weberlab-hhu/Helixer) is a deep-learning *ab initio* gene
finder that works across species without retraining, but it emits exactly one transcript
per gene and carries no evidence support. HelixForge keeps that gene set as a stable
anchor, reconciles it against RNA-seq assemblies, splice junctions, and protein
alignments through [Mikado](https://github.com/EI-CoreBioinformatics/mikado), and writes a
tiered, isoform-aware, ID-stable annotation.

The design principle: the predictor is the backbone. Evidence refines a model only where
it demonstrably improves it, and nothing the predictor gets right is discarded.

---

## What it does

- **Reconciles, does not replace.** RNA-seq and protein evidence are scored against each
  Helixer model; a model changes only when the evidence outscores the original.
- **Adds the missing isoforms.** Helixer is structurally capped at one transcript per
  gene. HelixForge reconstructs alternatives where transcript evidence supports them.
- **Uses the predictor's own confidence.** Helixer per-base confidence is injected as an
  external metric, targeting correction where the predictor is least certain.
- **Rescues silent genes.** Predicted genes without RNA-seq support are preserved when
  they carry a valid ORF, so recall is not traded for precision.
- **Tiered output.** Every gene is labeled by evidence tier, so downstream users filter
  by confidence rather than trusting a flat annotation.
- **No species-specific training, no manual curation.** Scales across HPC via genome
  chunking.

---

## Installation

Requires Python 3.10 or newer.

```bash
pip install helixforge
```

From source:

```bash
git clone https://github.com/aseetharam/helixforge
cd helixforge
pip install -e ".[cli,stats]"
```

Optional extras:

| Extra | Adds |
|---|---|
| `cli` | command-line interface (click, rich) |
| `stats` | tidy stats tables |
| `viz` | plotting and locus figures |
| `bigwig`, `tracks` | bigWig coverage support |
| `hpc` | [HyperShell](https://hypershell.readthedocs.io) distributed execution backend |
| `dev` | test and type-check baseline |

### External tools

HelixForge drives several bioinformatics tools as subprocesses. They are not Python
dependencies and must be installed separately (conda or container): Mikado, TransDecoder,
DIAMOND, Portcullis, miniprot, StringTie, STAR, and samtools.

Verify the environment before a long run:

```bash
helixforge doctor
```

This resolves every external tool, checks work-directory writability, and validates the
scoring profile, so a misconfigured environment fails in seconds rather than hours.

---

## Quick start

```bash
helixforge reconcile \
  --genome genome.fa \
  --helixer helixer.gff3 \
  --helixer-h5 helixer_predictions.h5 \
  --stringtie-list stringtie.fofn \
  --bam-list bam.fofn \
  --star-sj-list sj.fofn \
  --miniprot proteins_vs_genome.gff3 \
  --protein-db uniprot_plants.fasta \
  --scoring-profile permissive \
  --threads 16 \
  --output-prefix sample
```

`--protein-db` gates the Mikado reconciliation chain and is required for isoform
discovery. `--miniprot` supplies backstop CDS only and is not a substitute for it.

### Whole genomes

Chunk the genome and distribute the work:

```bash
helixforge reconcile ... --scatter 64 --hpc slurm
helixforge reconcile ... --scatter 64 --hpc hypershell
```

Chunked runs reserve disjoint gene-ID ranges, stitch loci across chunk boundaries, and
verify genome-wide ID uniqueness on aggregation. Long runs are resumable: completed
stages are detected and reused, so a failure late in the pipeline does not repeat the
expensive Mikado stage.

---

## Commands

| Command | Purpose |
|---|---|
| `reconcile` | the main pipeline: prep, Mikado, reconcile, output |
| `doctor` | preflight: external tools, paths, config, scoring profiles |
| `evidence` | standalone RNA-seq and protein evidence scoring (AED) |
| `confidence` | per-locus confidence from the Helixer HDF5 predictions |
| `stats` | annotation summary statistics |
| `viz` | locus and genome-level figures |
| `utils` | format conversion, filtering, protein extraction, database fetch |
| `parallel` | chunk planning, task generation, aggregation |

Each command's full option list is available via `--help`.

---

## Results

Benchmarked on *Zea mays* B73, chromosomes 1 to 10, against the NAM v5 reference
annotation (`Zm00001eb.1`). Comparisons use gffcompare 0.12.10 and Mikado 2.3.4, with raw
Helixer as the before baseline. RNA-seq evidence is 21 B73 samples.

### Structural agreement with the reference

F1 against NAM, all isoforms:

| Level | Helixer | HelixForge |
|---|---|---|
| Exon | 54.9 | **61.4** |
| Intron | 62.7 | **72.1** |
| Intron chain | 24.6 | **42.5** |
| Locus | 37.5 | **55.5** |

Transcripts whose entire intron chain exactly matches a NAM transcript rise from 12,308
to 27,096. The gain holds when both annotations are reduced to one transcript per gene,
so it is not an artifact of emitting more transcripts, and it holds at the coding (CDS)
level, so it is not confined to UTRs. Improvement is uniform across all ten chromosomes
(intron-chain sensitivity standard deviation 0.6), with no regressing chromosome.

### Isoform reconstruction

| | Genes | Transcripts | Transcripts per gene |
|---|---|---|---|
| Helixer | 44,480 | 44,480 | 1.00 |
| HelixForge | 44,307 | 71,639 | 1.62 |
| NAM (reference) | 39,035 | 71,791 | 1.84 |

HelixForge recovers most of the reference's isoform multiplicity from a single-isoform
prior.

### RNA-seq validates the added splice junctions

Splice junctions absent from the reference, adjudicated against the RNA-seq evidence:

| | Novel junctions with RNA support | Reference RNA-supported junctions recovered |
|---|---|---|
| Helixer | 5.0% | 77.9% (32,267 missed) |
| HelixForge | **22.4%** | **93.6%** (9,353 missed) |

The junctions HelixForge adds are RNA-supported 4.5 times more often than the input
predictor's, and it simultaneously recovers more of the reference's real junctions rather
than trading one for the other. The result is robust to the read-depth threshold, and
98.5% of supported novel junctions carry canonical GT-AG motifs.

### Coding structure

Median CDS length is 1,026 nt for HelixForge and 1,023 nt for NAM; raw Helixer under-calls
at 855 nt. Reconciliation restores reference-like coding length, not only UTR extent.

---

## Known limitations

Stated plainly, because they bound what the numbers above support.

- **Novel content is junctions, not genes.** HelixForge calls coding loci that NAM lacks,
  and most of those are not RNA-supported. Treat them as candidates, not discoveries.
- **ORF completeness trails the reference.** Roughly 86% of HelixForge ORFs are complete
  versus 93% for NAM, driven by 5'-partial models. This is not parity.
- **Improvement is net, not universal.** Reconciliation adds far more exact transcript
  matches than it degrades, but a small number of models do get worse.
- **Unsupported novel junctions are ambiguous.** They are either unexpressed in the
  sampled tissues or genuine over-calls. The available evidence cannot separate the two.
- **Isoform precision is the open axis.** The extra isoforms HelixForge emits are largely
  plausible, but tightening isoform admission would raise chain-level precision.
- **Requires protein evidence.** Without a protein database the Mikado chain does not run
  and the output is backstop-only, with no isoform discovery.

---

## How it works

```
Helixer prediction  ─┐
  + per-base confidence
                     │
RNA-seq              ├──▶  Mikado-anchored reconciliation  ──▶  tiered, isoform-aware
  StringTie + junctions                                          evidence-backed GFF3
                     │     score before modifying
Protein homology     ┘     keep a change only if it outscores
  miniprot alignments      call isoforms where evidence supports
                           confidence-gated splice correction
                           rescue silent genes with valid ORFs
```

The pipeline runs in four stages: prep (evidence preparation), Mikado (reconciliation
engine), reconcile (model selection, backstop rescue, tiering), and output (validated
GFF3 plus reports). All CDS are stop-inclusive, following the Helixer and Ensembl
convention. A gene whose CDS cannot be made structurally coherent is flagged and emitted
without a CDS rather than repaired heuristically.

---

## Citing

If HelixForge is useful in your work, please cite it. Distributed execution on HPC uses
HyperShell:

> Lentner, G. and Gorenstein, L. (2022). HyperShell v2: Distributed Task Execution for
> HPC. *Practice and Experience in Advanced Research Computing (PEARC '22)*, ACM.
> doi:10.1145/3491418.3535138

---

## License

MIT.
