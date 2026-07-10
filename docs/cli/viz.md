# `helixforge viz`

Render per-locus plots, interactive pages, or browser tracks for a run.

Visualizes a reconciled run. Because the rich ``ReconciledGene`` set is not
serialised to disk, ``viz`` takes the same inputs as ``reconcile`` and
rebuilds the gene set via the pipeline before plotting — so it needs the full
pipeline option set, not a finished GFF3.

Output (written under --out-dir), by --mode:
- static:      one figure per locus (--fmt svg/pdf/png); --top-n keeps the N
               most AS-complex loci, --gene restricts to a single gene id
- interactive: an HTML index of interactive per-locus pages
- tracks:      a BED12 browser track (<output-prefix>.bed)

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
| `--out-dir` | Directory for the figures / tracks. (required) |
| `--mode` | Output kind. (choices: static/interactive/tracks; default: `static`) |
| `--gene` | Restrict to a single gene id. |
| `--top-n` | Static mode: only the N most AS-complex loci. |
| `--fmt` | Static figure format (svg/pdf/png). (default: `svg`) |

## Examples

```bash
  # static SVG of the 50 most AS-complex loci
  helixforge viz --genome genome.fa --helixer helixer.gff3 \
      --helixer-h5 helixer.h5 --stringtie-list stringtie.list \
      --out-dir figures --mode static --top-n 50

  # one gene, as a PDF
  helixforge viz --genome genome.fa --helixer helixer.gff3 \
      --out-dir figures --mode static --gene HFG_00042 --fmt pdf

  # BED12 browser track for the whole run
  helixforge viz --genome genome.fa --helixer helixer.gff3 \
      --out-dir tracks --mode tracks
```
