# `helixforge parallel tasks`

Expand a command template over the plan → an executor-agnostic task file.

Step 2 of the chunked path. Expands a command template over every chunk in
plan.json; the default template is a per-chunk ``helixforge reconcile`` built
from the inputs you attach, so it works out of the box. Run the resulting file
with any executor (GNU parallel, a Slurm array, xargs, HyperShell), then
``parallel aggregate`` the per-chunk outputs.

Template placeholders, substituted per chunk:
{chunk_id} {region} {seqid} {start} {end} {start_0} {end_0} {size}
{id_start} {novel_start} {output_dir}
({region}/{start_0}/{end_0} are internal 0-based; {seqid}/{start}/{end} are
1-based inclusive for the CLI.)

Output:
- -o/--output: the task file, one command per chunk (one line each)
- --wrapper: an optional reusable per-task wrapper script (with
  --wrapper-setup lines) so cluster setup runs once per task

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
| `--plan` | plan.json from `parallel plan`. (required) |
| `--command` | Command template with chunk placeholders ({chunk_id} {region} {seqid} {start} {end} {start_0} {end_0} {size} {id_start} {novel_start} {output_dir}). Default: a per-chunk `helixforge reconcile` built from the attached inputs. |
| `-o`, `--output` | Task file to write (one command per chunk). (default: `tasks.txt`) |
| `--output-dir` | Per-chunk output directory ({output_dir} placeholder). (default: `chunks`) |
| `--wrapper` | Also write a reusable per-task wrapper script (setup once). |
| `--wrapper-setup` | Setup line for the wrapper (module load / conda activate); repeatable. (repeatable) |
| `--include-logging` | Redirect each task's stdout/stderr to <output-dir>/logs/<chunk>.log. (default: `False`) |
| `--helixforge-bin` | helixforge binary used in the default command template. (default: `helixforge`) |

## Examples

```bash
  # default template (per-chunk reconcile from the attached inputs)
  helixforge parallel tasks --plan plan.json \
      --genome maize.fa --helixer helixer.gff3 --helixer-h5 preds.h5 \
      --stringtie-list stringtie.list -o tasks.txt
  parallel -j 16 < tasks.txt          # or: hs cluster tasks.txt --num-tasks 16

  # full override + cluster setup wrapper
  helixforge parallel tasks --plan plan.json --genome g.fa --helixer h.gff3 \
      --command 'helixforge reconcile --genome g.fa --helixer h.gff3 \
                 --region {region} --id-base {id_start} \
                 --output-prefix {output_dir}/{chunk_id}' \
      --wrapper run_chunk.sh --wrapper-setup 'module load helixforge' -o tasks.txt
```
