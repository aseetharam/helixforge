#!/usr/bin/env python3
"""M3 validation — full HelixForge v3 pipeline (Phase 8) on a real genome.

Unlike ``scripts/m2_validate.py`` (which wired the building blocks by hand), this
driver calls the **single** Phase-8 entry point ``run_pipeline(PipelineConfig)``
on the real Arabidopsis set in ``helixforge_testing/``. It is the M3 milestone:
the orchestrator runs stages A→D end-to-end and produces the tiered GFF3s + TSV
report.

Reusing the M2 Mikado run (``--resume``)
----------------------------------------
The expensive MIKADO steps (prepare, TransDecoder, DIAMOND, serialise, pick)
were already computed at M2 and live under ``m2_out/mikado_run/`` (see
``docs/M2_VALIDATION.md``). ``--resume`` (the default) monkeypatches the
pipeline's ``run_*`` step functions to return those existing artifacts instead
of re-running hours of external tools, so M3 exercises the **new** Phase-8 code
on the real reconciled gene set: orchestration, Phase-7 finalisation (backstop
CDS + junction correction + structural codon gate), tiered GFF3 output, and the
report. ``--no-resume`` runs the whole chain afresh (needs the conda env + GPU-
free external tools; expect hours).

Environment (docs/M2_VALIDATION.md — read it):
  source ~/miniforge3/etc/profile.d/conda.sh && conda activate mikado-env
  export MIKADO_BIN=mikado DIAMOND_BIN=diamond \
         TRANSDECODER_BIN_DIR=$HOME/miniforge3/envs/mikado-env/bin
  /usr/bin/python3 -u scripts/m3_run.py            # resume from m2_out (default)
  /usr/bin/python3 -u scripts/m3_run.py --no-resume

Orchestration runs under system python3 (it has the helixforge deps); external
tools resolve from the active conda env via PATH / the env vars above.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Pin to THIS repo's package (an unrelated helixforge.v4 checkout may be the
# installed one on sys.path; we must not import that).
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import helixforge.reconcile.pipeline as pl
from helixforge.reconcile.pipeline import PipelineConfig, run_pipeline
from helixforge.stats.summary import summarize_genes

# ──────────────────────────────────────────────────────────────────────────
# Real-data paths (helixforge_testing/, gitignored)
# ──────────────────────────────────────────────────────────────────────────
DATA = REPO / "helixforge_testing"

GENOME      = DATA / "genome" / "athaliana.fasta"
HELIXER_GFF = DATA / "helixer_output" / "Arabidopsis-thaliana_helixer.gff3"
PROTEIN_FA  = DATA / "swissprot" / "uniprot_sprot.fasta"

HELIXER_INPUT_H5 = DATA / "helixer_output" / "Arabidopsis-thaliana_input.h5"
HELIXER_PRED_H5  = DATA / "helixer_output" / "Arabidopsis-thaliana_predictions.h5"

STRINGTIE_DIR = DATA / "stringtie"
SJ_DIR        = DATA / "SJ_out"
BIGWIG_DIR    = DATA / "bigwig"
BAM_DIR       = DATA / "bamfiles"

MIKADO_BIN          = os.environ.get("MIKADO_BIN", "mikado")
DIAMOND_BIN         = os.environ.get("DIAMOND_BIN", "diamond")
TRANSDECODER_BINDIR = os.environ.get("TRANSDECODER_BIN_DIR")


def install_resume(reuse_dir):
    """Monkeypatch the pipeline's Mikado step funcs to reuse M2 artifacts.

    Each replacement returns the existing path under ``reuse_dir`` rather than
    invoking the external tool. ``_build_external_scores`` is short-circuited too
    (its TSV feeds ``serialise``, which is skipped on resume), so M3 does not
    recompute Helixer support for ~194k prepared transcripts.
    """
    reuse_dir = Path(reuse_dir)

    def reuse(name, label):
        p = reuse_dir / name
        if not p.exists():
            sys.exit(f"ERROR: --resume needs {label} but it is missing: {p}\n"
                     f"Run M2 first (scripts/m2_validate.py) or use --no-resume.")
        print(f"   [resume] reuse {label}: {p.name}")
        return p

    def fake_prepare(cfg, out_dir, procs=1, mikado_bin="mikado"):
        return (reuse("mikado_prepared.gtf", "prepared GTF"),
                reuse("mikado_prepared.fasta", "prepared FASTA"))

    def fake_transdecoder(prepared_fasta, out_dir, **kw):
        return reuse("mikado_prepared.fasta.transdecoder.bed", "TransDecoder BED")

    def fake_diamond(prepared_fasta, db, out_dir, **kw):
        return reuse("mikado_diamond.xml", "DIAMOND XML")

    def fake_serialise(*a, **k):
        return reuse("mikado.db", "mikado.db")

    def fake_pick(cfg, scoring, prepared_gtf, out_dir, procs=1, mikado_bin="mikado"):
        return reuse("mikado.loci.gff3", "loci GFF3")

    # metrics/scores TSVs are read by parse_loci_gff3 from {work}/mikado_run;
    # on resume that directory is empty, so point the parse at reuse_dir by
    # symlinking the three files the parser needs into the work run dir.
    pl.run_prepare = fake_prepare
    pl.run_transdecoder = fake_transdecoder
    pl.run_diamond = fake_diamond
    pl.run_serialise = fake_serialise
    pl.run_pick = fake_pick
    pl._build_external_scores = lambda *a, **k: {}


def build_config(outdir, *, reuse_dir, id_map, resume=True, region=None,
                 procs=12, threads=12):
    """Build the M3 ``PipelineConfig`` (and apply resume wiring if requested).

    Shared by ``main()`` and the golden-output regression test
    (``tests/test_m3_golden.py``) so both drive the pipeline through an
    identical config. On ``resume`` this monkeypatches the Mikado step funcs
    and symlinks the reused loci/metrics/scores into the work dir.
    """
    for p in (GENOME, HELIXER_GFF, HELIXER_INPUT_H5, HELIXER_PRED_H5, PROTEIN_FA):
        if not p.exists():
            sys.exit(f"ERROR: missing input: {p}")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    work_dir = outdir / "work"

    if resume:
        print(f">> RESUME mode: reusing Mikado artifacts from {reuse_dir}")
        install_resume(reuse_dir)
        # parse_loci_gff3 reads loci.gff3 + metrics + scores from {work}/mikado_run;
        # symlink the three reused files there so the parser finds them.
        run_dir = work_dir / "mikado_run"
        run_dir.mkdir(parents=True, exist_ok=True)
        for name in ("mikado.loci.gff3", "mikado.loci.metrics.tsv",
                     "mikado.loci.scores.tsv"):
            link = run_dir / name
            src = Path(reuse_dir) / name
            if not link.exists() and src.exists():
                link.symlink_to(src.resolve())

    # Coverage fallback for the ~SILENT loci with no StringTie overlap: prefer
    # bigWig (fast) if pyBigWig is installed, else the BAMs (pysam, slower).
    import importlib.util
    if importlib.util.find_spec("pyBigWig") is not None:
        bigwig_paths = [str(p) for p in sorted(BIGWIG_DIR.glob("*.bw"))]
        bam_paths = []
        print(f"   coverage fallback: {len(bigwig_paths)} bigWigs")
    else:
        bigwig_paths = []
        bam_paths = [str(p) for p in sorted(BAM_DIR.glob("*.bam"))]
        print(f"   coverage fallback: {len(bam_paths)} BAMs "
              "(install pyBigWig for the fast path)")

    return PipelineConfig(
        genome_fasta=str(GENOME),
        helixer_gff3=str(HELIXER_GFF),
        # Phase 13: HDF5ConfidenceReader reads the native Helixer split layout
        # directly (auto-detects the *_predictions.h5 sibling of this *_input.h5),
        # so no hand-built combined view is needed any more.
        helixer_h5=str(HELIXER_INPUT_H5),
        stringtie_list=[str(p) for p in sorted(STRINGTIE_DIR.glob("*.gtf"))],
        star_sj_paths=[str(p) for p in sorted(SJ_DIR.glob("*_SJ.out.tab"))],
        bam_paths=bam_paths,
        bigwig_paths=bigwig_paths,
        protein_db=str(PROTEIN_FA),
        transdecoder_bin_dir=TRANSDECODER_BINDIR,
        mikado_bin=MIKADO_BIN,
        diamond_bin=DIAMOND_BIN,
        scoring_profile="strict",
        output_prefix=str(outdir / "athaliana"),
        report_path=str(outdir / "athaliana.report.tsv"),
        id_map_path=str(id_map),
        work_dir=str(work_dir),
        region=region,
        procs=procs,
        threads=threads,
        # On resume the templated config is never consumed (pick is skipped);
        # avoids a `mikado configure` subprocess. A fresh run should set this True.
        use_mikado_configure=not resume,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default=str(REPO / "m3_out"),
                    help="output directory (GFF3s, report, work/) [m3_out]")
    ap.add_argument("--reuse-mikado-run", default=str(REPO / "m2_out" / "mikado_run"),
                    help="dir holding the M2 Mikado artifacts to reuse on --resume")
    ap.add_argument("--id-map", default=str(REPO / "m2_out" / "id_map.json"),
                    help="persisted Helixer→HFG id map (reused for stable IDs)")
    ap.add_argument("--resume", dest="resume", action="store_true", default=True,
                    help="reuse the M2 Mikado artifacts (default)")
    ap.add_argument("--no-resume", dest="resume", action="store_false",
                    help="run the full Mikado chain afresh (needs external tools)")
    ap.add_argument("--region", default=None,
                    help="restrict to one scaffold, e.g. 'chr1' (smoke test)")
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--threads", type=int, default=12)
    # Phase 33b: re-elect each multi-isoform gene's canonical transcript via the
    # TRaCE ranked-choice election (RNA-seq + length + optional domain voters).
    # Default OFF so _build_config / the golden gate stay on combined_score
    # primaries; this is an opt-in passthrough to PipelineConfig.trace_primary.
    ap.add_argument("--trace-primary", dest="trace_primary", action="store_true",
                    default=False,
                    help="elect canonical isoforms with TRaCE (Phase 33b)")
    args = ap.parse_args()

    config = build_config(
        args.outdir,
        reuse_dir=args.reuse_mikado_run,
        id_map=args.id_map,
        resume=args.resume,
        region=args.region,
        procs=args.procs,
        threads=args.threads,
    )
    config.trace_primary = args.trace_primary
    if config.trace_primary:
        print(">> TRaCE primary election: ON (--trace-primary)")

    print(">> Running pipeline (stages A→D)...")
    genes = run_pipeline(config)

    # ── M3 summary ────────────────────────────────────────────────────────────
    # summarize_genes() is the SAME helper the golden-output regression test
    # asserts against (CLAUDE.md §15, gate clause 3); keep counting here in one
    # place so the script and the test never drift.
    s = summarize_genes(genes)
    n = s["genes"]

    print("\n" + "=" * 60)
    print("M3 PIPELINE SUMMARY")
    print("=" * 60)
    print(f"Reconciled genes:       {n}")
    print(f"Total transcripts:      {s['total_tx']}")
    print(f"Multi-isoform genes:    {s['multi']} ({100*s['multi']/max(n,1):.1f}%)")
    print(f"Coding genes (any CDS): {s['coding']}")
    print(f"Tier:                   {s['tier']}")
    print(f"Origin:                 {s['origin']}")
    print(f"AS events:              {s['as_kinds']}")
    print(f"Flags:                  {s['flags']}")
    print(f"\nOutputs in {args.outdir}:")
    for suffix in (".gff3", ".tier1.gff3", ".tier2.gff3", ".tier3.gff3"):
        f = Path(f"{config.output_prefix}{suffix}")
        if f.exists():
            print(f"   {f.name}  ({f.stat().st_size/1e6:.1f} MB)")
    rep = Path(config.report_path)
    if rep.exists():
        print(f"   {rep.name}  ({rep.stat().st_size/1e6:.1f} MB)")
    print(f"   id_map: {config.id_map_path}")


if __name__ == "__main__":
    main()
