# `helixforge reconcile`

Reconcile Helixer models with evidence → tiered GFF3 + per-gene report.

The main HelixForge command. It runs PREP → MIKADO → RECONCILE → OUTPUT end
to end: prepares Mikado inputs from the Helixer models plus RNA-seq and
protein evidence, reconciles Mikado loci against the Helixer gene set
(rescuing silent genes, stabilising HFG ids, calling isoforms), applies the
structural codon gate, and writes the tiered annotation. Unlike ``evidence``
(which only *scores* an existing GFF3) this command *builds and fixes* models.

Execution mode:
- default (``--scatter off``): one in-process pass over the whole genome
  or a single ``--region``.
- ``--scatter auto|<N>``: partition the genome into gene-safe chunks, run
  them (``--hpc local`` now, or ``--hpc slurm`` to emit an array to submit),
  and aggregate into one verified annotation (+ a run manifest). Chunking
  changes memory/speed, not the result. For a fully scheduler-driven
  distributed run, drive ``parallel plan/tasks/aggregate`` directly.

Outputs (``--output-prefix`` is the stem):
- <prefix>.gff3                       full reconciled annotation (1-based GFF3)
- <prefix>.tier1/2/3.gff3             cumulative per-tier GFF3 subsets
- <prefix>.report.tsv                 per-gene report (one row per gene)
- <prefix>.id_map.json                Helixer-locus → HFG id map (stable reruns)
- <prefix>.report.json / .report.html / .run_stats.json   run summaries
- (chunked) <out-prefix>.manifest.json   plan + chunk prefixes + id_map

Per-gene report (<prefix>.report.tsv) columns:
gene_id, seqid, start, end, strand, tier, origin, biotype,
primary_transcript_id, num_isoforms, num_as_events, has_cds, protein_id,
max_tpm, junction_support, helixer_support, combined_score, aed, flags.

## Options

| Option | Description |
|---|---|
| `--genome` | Genome FASTA (required). (required) |
| `--helixer` | Helixer GFF3 (required). (required) |
| `--helixer-h5` | Helixer confidence HDF5: a combined file, the '*_input.h5' metadata half, or the '*_predictions.h5' softmax half (the partner half is auto-detected as the sibling, or named via --helixer-input-h5). |
| `--helixer-input-h5` | Helixer '*_input.h5' metadata half (coordinate mapping). Needed only when --helixer-h5 is a bare '*_predictions.h5' whose '*_input.h5' sibling is not co-located. Mirrors 'confidence --input-h5'. |
| `--stringtie` | StringTie GTF, one per sample. Comma-separated or repeated. (repeatable) |
| `--stringtie-list` | File of StringTie GTF paths, one per line (blank lines and # comments ignored). |
| `--bam` | RNA-seq BAM. Comma-separated or repeated. (repeatable) |
| `--bam-list` | File of BAM paths, one per line (blank lines and # comments ignored). |
| `--bigwig` | Coverage bigWig. Comma-separated or repeated. (repeatable) |
| `--bigwig-list` | File of bigWig paths, one per line (blank lines and # comments ignored). |
| `--star-sj` | STAR SJ.out.tab. Comma-separated or repeated. (repeatable) |
| `--star-sj-list` | File of SJ.out.tab paths, one per line (blank lines and # comments ignored). |
| `--miniprot` | miniprot GFF for backstop CDS projection. Not a substitute for --protein-db; does not enable Mikado. |
| `--protein-db` | Protein FASTA for DIAMOND homology. Required (with StringTie) to enable Mikado reconciliation and isoform discovery. |
| `--te-annotation` | EDTA TE GFF3. Flags TE-overlapping models (TE_OVERLAP) and reclassifies good-ORF genes above --te-overlap-threshold as transposable_element. Uses the EDTA Classification order; knob/satellite/low-complexity never gate. |
| `--te-overlap-threshold` | Model-fraction TE overlap at/above which a good-ORF gene is reclassified transposable_element (default 0.5). (default: `0.5`) |
| `--te-class` | EDTA Classification order to treat as a TE (repeatable). Default covers true TE orders (LTR, DNA, MITE, TIR, Helitron, LINE, SINE) and excludes knob/satellite/low-complexity. (repeatable) |
| `--transdecoder-bin-dir` |  |
| `--backstop-transdecoder`, `--no-backstop-transdecoder` | (default: `False`) |
| `--portcullis-bin` |  |
| `--mikado-bin` | (default: `mikado`) |
| `--diamond-bin` | (default: `diamond`) |
| `--mikado-configure`, `--no-mikado-configure` | (default: `True`) |
| `--functional-annotation`, `--no-functional-annotation` | Run the opt-in functional-annotation hook (InterProScan / eggNOG). (default: `False`) |
| `--functional-tool` | Annotation tool(s) to run for --functional-annotation. (choices: interproscan/eggnog/both; default: `interproscan`) |
| `--functional-db` | eggNOG data dir (required for --functional-tool eggnog/both). |
| `--interproscan-bin` | InterProScan binary (bare command or absolute path). (default: `interproscan.sh`) |
| `--eggnog-bin` | eggNOG-mapper binary (bare command or absolute path). (default: `emapper.py`) |
| `--transl-table` | NCBI genetic-code table id for CDS validation (default 1 = standard nuclear code). (default: `1`) |
| `--transl-table-map` | Per-seqid genetic-code overrides; lines 'seqid=table' (e.g. chrMt=1, chrPt=11). 'seqid<TAB/space>table' also accepted. Overrides --transl-table for the named seqids. |
| `--scoring-profile` | (choices: strict/permissive; default: `strict`) |
| `--helixer-reference`, `--no-helixer-reference` | Inject Helixer as is_reference=true. (default: `True`) |
| `--helixer-support-weight` | Scale external helixer_support (0 disables coupling). (default: `1.0`) |
| `--output-prefix` | (default: `helixforge`) |
| `--report-path` |  |
| `--id-map-path` |  |
| `--work-dir` |  |
| `--region` | seqid or seqid:start-end. |
| `--chunk-id` |  |
| `--id-base` | Lowest HFG number for this chunk's reserved range. (default: `1`) |
| `--novel-base` | Lowest HFG number for this chunk's novel sub-range. (default: `90000`) |
| `--procs` | (default: `1`) |
| `--threads` | (default: `4`) |
| `--min-tpm` | (default: `0.5`) |
| `--min-samples` | (default: `1`) |
| `--coverage-threshold` | (default: `2.0`) |
| `--near-zero-coverage` | (default: `0.1`) |
| `--as-report`, `--no-as-report` | (default: `True`) |
| `--only-confirmed-introns`, `--allow-unconfirmed-introns` | (default: `True`) |
| `--max-isoforms` | (default: `5`) |
| `--keep-retained-introns`, `--drop-retained-introns` | (default: `False`) |
| `--pad`, `--no-pad` | Unify isoform termini (the strongest consistency lever). (default: `True`) |
| `--chimera-split`, `--no-chimera-split` | (default: `True`) |
| `--flank` | (default: `200`) |
| `--reciprocal-overlap` | (default: `0.5`) |
| `--min-cds-overlap` | (default: `0.6`) |
| `--min-cdna-overlap` | (default: `0.6`) |
| `--admit-novel`, `--no-admit-novel` | (default: `False`) |
| `--novel-evidence-floor` |  |
| `--trace-primary`, `--no-trace-primary` | Elect the canonical/primary isoform per gene by ranked-choice voting (TRaCE) instead of highest combined_score. (default: `False`) |
| `--trace-max-aed` | TRaCE: max AED for a sample to vote for a candidate. (default: `0.5`) |
| `--trace-min-tpm` | TRaCE: min TPM for a sample's assembled transcript to vote. (default: `0.5`) |
| `--trace-min-overlap` | TRaCE: min proportion-overlap for a sample transcript to vote. (default: `0.5`) |
| `--trace-weight-domain` | TRaCE: domain-coverage voter weight. (default: `9.0`) |
| `--trace-weight-protein` | TRaCE: protein (CDS) length voter weight. (default: `6.0`) |
| `--trace-weight-cdna` | TRaCE: transcript (cDNA) length voter weight. (default: `3.0`) |
| `--trace-use-domain`, `--no-trace-use-domain` | TRaCE: include the domain-coverage voter when per-isoform domain coverage is available. (default: `True`) |
| `--short-cds-threshold` | (default: `300`) |
| `--short-exon-threshold` | (default: `10`) |
| `--long-intron-threshold` | (default: `100000`) |
| `--junction-tolerance` | (default: `0`) |
| `--junction-min-reads` | (default: `3`) |
| `--scatter` | Chunked execution: off \| auto \| <N chunks>. 'off' is a single in-process pass; 'auto'/<N> partition the genome into gene-safe chunks (same result, less memory). For a scheduler-driven distributed run drive `parallel plan/tasks/aggregate` yourself. (default: `off`) |
| `--hpc` | Chunked mode only: run chunks locally now, run them now via HyperShell (`hs cluster`), or emit a Slurm array to submit. (choices: local/hypershell/slurm; default: `local`) |
| `--workers` | Chunked --hpc local: process-pool size (also the default HyperShell --num-tasks concurrency). (default: `4`) |
| `--hs-bin` | Chunked --hpc hypershell: the HyperShell executable (install via `pip install helixforge[hpc]`). (default: `hs`) |
| `--target-loci-per-chunk` |  |
| `--min-boundary-gap` |  |
| `--out-prefix` | Chunked mode: prefix for the merged outputs + manifest (defaults to --output-prefix). |
| `--master-id-map` | Chunked mode: genome-wide master id_map.json (defaults to --id-map-path). |
| `--stitch`, `--no-stitch` | Chunked mode: report boundary-recoverable cross-chunk merges after aggregate. (default: `False`) |
| `--finalize-workers` | Worker processes for the per-gene finalize stage (intra-chunk). (default: `1`) |
| `--seqid-aliases` | Chunked mode: seqid alias map (JSON {alias: canonical} or 2-col TSV) for the integrity preflight. |
| `--skip-preflight` | Chunked mode: skip the input-integrity preflight (not recommended). (default: `False`) |
| `--resume`, `--no-resume` | Reuse completed stage outputs (PREP, Mikado) from the work dir instead of re-running them. A stage is skipped only if its outputs exist and validate (non-empty + content hash). --no-resume forces a fresh run even when valid outputs exist. (default: `True`) |
| `--mikado-loci` | Path to an externally-run Mikado loci GFF3 (e.g. mikado.loci.gff3). Skips the entire PREP + Mikado stage; companion files (*.metrics.tsv, *.scores.tsv) are auto-detected next to the GFF3 and used when present. The rest of the pipeline (backstop, reconcile, output) runs normally on the supplied loci. |

## Examples

```bash
  # RNA-seq (StringTie + BAM) + protein homology, strict scoring (single pass)
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \
      --helixer-h5 helixer.h5 --stringtie sampleA.gtf,sampleB.gtf \
      --bam sampleA.bam,sampleB.bam --protein-db proteins.fa \
      --output-prefix run1

  # Protein-only backstop CDS from a precomputed miniprot GFF, permissive scoring
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \
      --miniprot miniprot.gff --scoring-profile permissive \
      --output-prefix run1

  # Chunked whole genome: partition into gene-safe chunks, run locally, aggregate
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \
      --helixer-h5 helixer.h5 --stringtie-list stringtie.list \
      --scatter auto --workers 16 --output-prefix maize

  # Chunked on a cluster: emit a Slurm array (submit it, then `parallel aggregate`)
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \
      --scatter 64 --hpc slurm --output-prefix maize

  # One region only with a reserved id range (the shape distributed runs emit)
  helixforge reconcile --genome genome.fa --helixer helixer.gff3 \
      --helixer-h5 helixer.h5 --region chr1:1-2000000 \
      --id-base 1 --output-prefix chunk_chr1
```
