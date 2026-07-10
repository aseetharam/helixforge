# `helixforge evidence`

Inspect (read-only): score any GFF3 against RNA-seq + protein evidence.

A standalone, evidence-only scorer. It works on *any* GFF3 (not just
HelixForge output) and needs only RNA-seq (BAM / STAR SJ / StringTie) and/or
protein evidence — no Mikado, no Helixer HDF5, no external toolchain beyond
an optional miniprot run. This command only *scores*; it never modifies a
model (use ``reconcile`` for that), and it scores against evidence rather
than the Helixer track (use ``confidence`` for that).

Reports two **separate** AED scores (never fused): ``rna_aed`` (junctions +
coverage + boundary, from BAM/SJ) and ``protein_aed`` (intron structure + CDS
coverage + reference-protein coverage, from miniprot). Each is in [0, 1],
lower = better, and is populated only when its evidence was supplied —
missing evidence is neutral, not a penalty.

Outputs:
- --out: per-transcript TSV, one row per transcript, columns:
  gene_id, transcript_id, seqid, strand, start, end, num_exons, num_introns,
  supported, contradicted, novel_in_data, junction_support_fraction,
  intron_precision, intron_recall, intron_f1, mean_coverage, tpm,
  then the RNA-AED block (when BAM coverage given): rna_aed,
  rna_junction_ratio, rna_coverage_ratio, rna_boundary_ratio,
  and the protein-AED block (when proteins/miniprot given): protein_id,
  protein_aed, protein_struct_ratio, protein_cds_cov_ratio,
  protein_prot_cov_ratio
- --gene-out: per-gene TSV (best-supported transcript per gene), same columns;
  defaults to '<out>.gene.tsv'
- a printed summary table (n_transcripts, fraction_fully_supported, mean
  intron F1, mean rna_aed / protein_aed, …)

## Options

| Option | Description |
|---|---|
| `--gff3` | Annotation GFF3 to score (any GFF3, not just HelixForge output). (required) |
| `--bam` | RNA-seq BAM (sorted + indexed). Comma-separated or repeated. (repeatable) |
| `--bam-list` | File of BAM paths, one per line (blank lines and # comments ignored). |
| `--sj` | STAR SJ.out.tab. Comma-separated or repeated. (repeatable) |
| `--sj-list` | File of SJ.out.tab paths, one per line (blank lines and # comments ignored). |
| `--stringtie` | StringTie GTF, one per sample (TPM). Comma-separated or repeated. (repeatable) |
| `--stringtie-list` | File of StringTie GTF paths, one per line (blank lines and # comments ignored). |
| `--proteins` | Reference proteome FASTA for the protein-AED axis. Comma-separated or repeated. Aligned with miniprot (needs --genome) unless --miniprot-gff is supplied. NOTE: use a TE-filtered proteome — TE proteins align well and would give TE models a deceptively low protein_aed. (repeatable) |
| `--proteins-list` | File of proteome FASTA paths, one per line (blank lines and # comments ignored). |
| `--miniprot-gff` | Precomputed miniprot GFF3 for --proteins; skips realigning on re-runs (parsed directly, no miniprot needed). |
| `--genome` | Genome FASTA — the miniprot target when aligning --proteins from scratch. Not needed with --miniprot-gff. |
| `--reference` | Genome FASTA for CRAM decode (passed as reference_filename to pysam, so CRAM never triggers a remote ENA fetch). Ignored for BAM. If omitted, htslib honours a REF_CACHE / REF_PATH env cache. |
| `--region` | Restrict to seqid or seqid:start-end (1-based). |
| `--out` | Per-transcript TSV output path. (default: `evidence.tsv`) |
| `--gene-out` | Per-gene TSV (best-supported transcript per gene). Defaults to '<out>.gene.tsv'. |
| `--min-reads` | Minimum junction read support to qualify as evidence. (default: `3`) |
| `--min-mapq` | Minimum read MAPQ for BAM junction extraction. (default: `10`) |
| `--min-overhang` | Minimum spliced-read overhang (bp) on each side of a junction. (default: `8`) |
| `-j`, `--threads` | Parallelism dial: process-parallel evidence extraction (one BAM/SJ per worker — a single coverage pass per locus, not per exon — plus the miniprot run) then process-parallel per-gene scoring. Output is identical to -j 1. (default: `1`) |

## Examples

```bash
  # RNA-seq only (junctions + coverage) → rna_aed
  helixforge evidence --gff3 annotation.gff3 \
      --bam sampleA.bam,sampleB.bam --sj sampleA.SJ.out.tab \
      --out evidence.tsv

  # add protein evidence by aligning a TE-filtered proteome with miniprot
  helixforge evidence --gff3 annotation.gff3 --bam sampleA.bam \
      --proteins refprot.fa --genome genome.fa --out evidence.tsv

  # reuse a precomputed miniprot GFF (no realignment) + StringTie TPM
  helixforge evidence --gff3 annotation.gff3 --stringtie-list stringtie.list \
      --proteins refprot.fa --miniprot-gff miniprot.gff \
      --out evidence.tsv --gene-out evidence.gene.tsv
```
