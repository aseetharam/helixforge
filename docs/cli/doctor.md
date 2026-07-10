# `helixforge doctor`

Preflight: validate inputs + resolve external tools before a run.

The single input-validation entry point, run it before a real
``reconcile``. With no inputs it resolves and versions every external tool
the pipeline shells out to (mikado, diamond, transdecoder, portcullis, …)
and checks environment hygiene (console-script shim, CRAM reference). Given
``--genome`` + ``--helixer`` (and optional evidence) it also runs the
input-integrity preflight: seqid concordance across inputs, GFF3/GTF/FASTA/BAM
format sniffing, and BAM/CRAM index (CSI) presence. ``--check-config``
verifies an emitted Mikado config/scoring file against the detected Mikado
version.

Nothing is modified and nothing is annotated, this command only inspects.
It exits non-zero if a required tool is missing, an emitted config fails
schema verification, or the input-integrity gate fails.

Reports (printed to the console; advisory checks warn, gates fail):
- external-tool table: tool, resolved path, detected version, status
- environment hygiene: console-script shim, CRAM-reference readiness
- (optional) Mikado config/scoring verification verdict
- (optional) input-integrity report: seqid concordance, format, index

## Options

| Option | Description |
|---|---|
| `--mikado-bin` | Override the mikado binary. |
| `--diamond-bin` | Override the diamond binary. |
| `--check-config` | Verify an emitted Mikado configuration.yaml against the detected version. |
| `--check-scoring` | Scoring YAML to verify alongside --check-config. |
| `--genome` | Genome FASTA (integrity check). |
| `--helixer` | Helixer GFF3 (integrity check). |
| `--helixer-h5` | Helixer HDF5 (integrity check). |
| `--stringtie` | StringTie GTF (repeatable). (repeatable) |
| `--bam` | RNA-seq BAM (repeatable). (repeatable) |
| `--star-sj` | STAR SJ.out.tab (repeatable). (repeatable) |
| `--miniprot` | miniprot GFF (integrity check). |
| `--reference` | Genome FASTA for CRAM decode, doctor reports whether CRAM BAM inputs have an offline reference (--reference or REF_CACHE). |
| `--seqid-aliases` | Seqid alias map (JSON {alias: canonical} or 2-col TSV). |

## Examples

```bash
  # just the external-tool preflight (mikado, diamond, …)
  helixforge doctor

  # tools + input-integrity (seqid concordance / format sniff / BAM index)
  helixforge doctor --genome genome.fa --helixer helixer.gff3 \
      --helixer-h5 helixer.h5 --bam sampleA.bam --star-sj sampleA.SJ.out.tab

  # verify an emitted Mikado config/scoring against the detected mikado version
  helixforge doctor --check-config configuration.yaml --check-scoring scoring.yaml
```
