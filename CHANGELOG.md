# Changelog

All notable changes to HelixForge are recorded here. The format loosely follows
[Keep a Changelog](https://keepachangelog.com/); versions follow PEP 440.

## [4.0.0.dev0] - 2026-06-28

### Added
- `helixforge utils` subcommand group with 8 standalone utilities
- `utils fetch-db`: download and Diamond-format reference protein databases
- `utils extract-proteins`: translate CDS from GFF3 + genome to protein FASTA
- `utils align`: miniprot protein-to-genome alignment wrapper
- `utils qc`: genome-wide QC report generation (HTML/JSON/TSV)
- `utils validate-inputs`: input format sniffing and seqid concordance
- `utils convert`: GFF3 <-> GTF format conversion
- `utils filter`: tiered output filtering by confidence, evidence, biotype
- `utils summarize`: annotation statistics (gene counts, biotypes, CDS lengths)
- `utils/databases.py`: DatabaseManager with UniProt/UniRef download and
  Diamond formatting (ported from v1)
- `utils/filters.py`: GeneFilter with preset profiles (ported from v1)
- `stats/genome_report.py`: genome-wide QC report with embedded charts

### Changed
- Removed 153 CLAUDE.md inline references from source code
- Replaced stale GUIDE.md (described v1 architecture) with current summary
- Normalized docstring tone across codebase
- Version bumped from 3.0.0.dev0 to 4.0.0.dev0

## [Unreleased]

### Fixed
- **`reconcile` now assigns StringTie TPM to every gene.** Mikado consumes
  StringTie TPM only as an external scoring metric and never emits it back, so
  reconciled transcripts carried `tpm=None` and the report's
  `frac_genes_tpm_pass` was always 0 even with valid StringTie input. The
  pipeline now assigns each transcript the TPM of the best exonic-overlap
  StringTie structure — the same rule the `evidence` scorer uses, now shared via
  `io.stringtie.build_tpm_overlap_index` / `best_overlapping_tpm` (one
  implementation, no drift). Count-neutral (sets `tpm` only).
- **Biotype is now driven by the model's own ORF, not by the presence of
  evidence.** A Helixer-only (backstop) gene previously lost the CDS Helixer
  itself predicted (it was rebuilt CDS-less) and was mis-called
  lncRNA/ncRNA whenever it lacked expression and homology — so only
  miniprot-rescued genes were coding. The backstop transcript now carries the
  Helixer model's intrinsic CDS, so a gene with a complete ORF is
  `protein_coding` regardless of expression/homology (absence of corroboration
  is neutral, not a non-coding signal). Expression/homology adjust confidence
  *within* coding: a good-ORF gene with neither is still coding but flagged
  `PUTATIVE_CODING`; miniprot homology now *confirms* a backstop's intrinsic ORF
  (re-tier to 1, record the accession) instead of overwriting its coordinates.

### Added
- **Optional EDTA TE gating (`--te-annotation`).** EDTA is the only signal that
  calls a transposable element. With a TE GFF3, each model's overlap with TE
  features is flagged `TE_OVERLAP`, and a good-ORF gene whose model-fraction TE
  overlap is at/above `--te-overlap-threshold` (default 0.5) is reclassified
  `transposable_element` and demoted to Tier 4 — gating the intrinsic-ORF coding
  call on TE-encoded transposase/gag-pol ORFs. The EDTA `Classification` order is
  used (tunable via `--te-class`; default LTR/DNA/MITE/TIR/Helitron/LINE/SINE),
  so knob/satellite/centromere/rDNA/low-complexity repeats never gate. Flag
  emission and reclassification are distinct (inspect before it changes output).
  With no `--te-annotation` it is a no-op — behavior is identical to before.

### Changed
- **Multi-file evidence inputs now follow the v1 convention.** Every multi-file
  evidence flag (`--bam`, `--sj`/`--star-sj`, `--stringtie`) accepts three
  merging forms: the **repeatable** flag, a **comma-separated** value, and a new
  `-list` **file-of-filenames (FOFN)** companion (`--bam-list`, `--sj-list` /
  `--star-sj-list`, `--stringtie-list`). All forms merge order-preserving and
  de-duplicated; every path is validated to exist (a FOFN error names the list
  file + line number). FOFN entries given as bare names resolve against the list
  file's own directory. Implemented once in `helixforge.io.fofn.expand_file_args`
  and wired through `evidence` and every `pipeline_options` command (`reconcile`,
  `viz`, `benchmark ablation`, `parallel tasks`, `run`).

### Breaking
- **`--stringtie` is now repeatable individual GTFs, not a sample-list file.**
  The sample-list file moves to the new `--stringtie-list` companion, matching
  `--bam`/`--sj`. A one-release back-compat shim is kept: if `--stringtie`
  receives a single existing non-`.gtf`/`.gff` file whose first real line
  resolves to an existing GTF, it is treated as a list-file with a deprecation
  warning pointing at `--stringtie-list`. The standalone `evidence` scorer's
  `score_annotation(..., stringtie_list=<file>)` keyword is likewise renamed to
  `stringtie_gtfs=<list of GTF paths>`.
