#!/usr/bin/env python3
"""M2 validation — end-to-end HelixForge v3 (Phases 0-6) on a real genome.

Runs the full PREP -> MIKADO -> RECONCILE chain on the Arabidopsis test set in
``helixforge_testing/`` and prints the M2 sanity summary (CLAUDE.md / DESIGN.md
"What to validate at M2"). There is no CLI yet (Phase 11) and no pipeline
orchestrator (Phase 8), so this script wires the building blocks directly.

Data flow:
  1. (Phase 1-3) load Helixer loci + HDF5 confidence; parse StringTie; classify.
  2. (Phase 4)   emit Mikado inputs: GTFs, junctions (from STAR SJ.out.tab),
                 list.txt, scoring profile, configuration.yaml.
  3. (Phase 5)   mikado prepare  ->  [compute external scores on PREPARED tids]
                 ->  TransDecoder + DIAMOND  ->  serialise  ->  pick.
  4. (Phase 5-6) parse mikado.loci.gff3; reconcile onto the Helixer gene set.

External metrics MUST be keyed to the *prepared* transcript ids (the output of
``mikado prepare``), which is why the Mikado run is broken into steps here rather
than calling ``run_mikado`` in one shot.

Requirements on PATH (activate your conda env / container first):
  mikado >= 2.3.2, TransDecoder ~5.7, diamond ~2.1  (portcullis optional)
Override any binary via env vars: MIKADO_BIN, DIAMOND_BIN, TRANSDECODER_BIN_DIR.

Usage:
  cd /home/arnstrm/svn/helixforge.v3
  python scripts/m2_validate.py                 # full genome
  python scripts/m2_validate.py --region Chr1   # one scaffold (fast smoke test)
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from collections import Counter
from pathlib import Path

# Pin to THIS repo's package (an unrelated helixforge.v4 checkout may be the
# installed one on sys.path; we must not import that).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from helixforge.io.hdf5 import HDF5ConfidenceReader
from helixforge.io.stringtie import StringTieParser
from helixforge.io.bam import parse_star_sj_tab
from helixforge.reconcile.models import Exon, SpliceJunction
from helixforge.reconcile.locus import load_helixer_loci, filter_loci_by_region
from helixforge.reconcile.classify import classify_loci
from helixforge.mikado.emit_gtf import helixer_gff3_to_gtf, stringtie_to_labelled_gtf
from helixforge.mikado.emit_junctions import junctions_to_portcullis_tab
from helixforge.mikado.emit_external import (
    helixer_support,
    helixer_locus_conf,
    normalize_tpm,
    write_external_scores_tsv,
)
from helixforge.mikado.config import (
    write_input_list,
    write_configuration,
    install_scoring_profile,
)
from helixforge.mikado.run import (
    run_prepare,
    run_transdecoder,
    run_diamond,
    run_serialise,
    run_pick,
)
from helixforge.mikado.parse import parse_loci_gff3
from helixforge.reconcile.mikado_integrate import reconcile

# ──────────────────────────────────────────────────────────────────────────
# Paths — everything is under helixforge_testing/
# ──────────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "helixforge_testing"

GENOME      = DATA / "genome" / "athaliana.fasta"
HELIXER_GFF = DATA / "helixer_output" / "Arabidopsis-thaliana_helixer.gff3"
PROTEIN_FA  = DATA / "swissprot" / "uniprot_sprot.fasta"

# Helixer writes two HDF5s: metadata (data/seqids, data/start_ends, data/X) in
# *_input.h5 and softmax in *_predictions.h5. The Phase 1 HDF5ConfidenceReader
# expects /predictions + /seqids + /start_ends in ONE file, so we build a tiny
# combined view (metadata copied, predictions external-linked → zero copy).
HELIXER_INPUT_H5 = DATA / "helixer_output" / "Arabidopsis-thaliana_input.h5"
HELIXER_PRED_H5  = DATA / "helixer_output" / "Arabidopsis-thaliana_predictions.h5"
HELIXER_H5       = DATA / "helixer_output" / "Arabidopsis-thaliana_combined.h5"

STRINGTIE_DIR = DATA / "stringtie"
SJ_DIR        = DATA / "SJ_out"
BIGWIG_DIR    = DATA / "bigwig"
BAM_DIR       = DATA / "bamfiles"

# External tool binaries (override via env for conda/container layouts)
MIKADO_BIN          = os.environ.get("MIKADO_BIN", "mikado")
DIAMOND_BIN         = os.environ.get("DIAMOND_BIN", "diamond")
TRANSDECODER_BINDIR = os.environ.get("TRANSDECODER_BIN_DIR")  # None -> PATH

SCORING_PROFILE = "strict"

# AS / isoform-quality knobs (CLAUDE.md §7) patched into the generated config.
PICK_AS_KNOBS = {
    "report": True,
    "only_confirmed_introns": True,
    "min_cds_overlap": 0.6,
    "min_cdna_overlap": 0.6,
    "max_isoforms": 5,
    "keep_retained_introns": False,
    "pad": True,            # strongest lever for *consistent* isoform termini
}


def _strip(name, *suffixes):
    for s in suffixes:
        if name.endswith(s):
            return name[: -len(s)]
    return name


def sample_id_from(path, kind):
    n = path.name
    if kind == "stringtie":
        return _strip(n, "_Aligned.sortedByCoord.out.gtf")
    if kind == "sj":
        return _strip(n, "_SJ.out.tab")
    if kind == "bigwig":
        return _strip(n, "_Signal.Unique.str1.out.bw")
    return n


# ──────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────
def ensure_combined_h5():
    """Build the combined Helixer HDF5 (metadata + external-linked predictions).

    The Phase 1 reader needs /predictions, /seqids, /start_ends in one file;
    real Helixer output splits them across *_input.h5 (``data/seqids``,
    ``data/start_ends``) and *_predictions.h5 (``/predictions``). We copy the two
    small metadata datasets and **external-link** predictions (no array copy).
    Rebuilds if missing or stale.
    """
    import h5py

    if HELIXER_H5.exists() and (
        HELIXER_H5.stat().st_mtime >= HELIXER_INPUT_H5.stat().st_mtime
        and HELIXER_H5.stat().st_mtime >= HELIXER_PRED_H5.stat().st_mtime
    ):
        return HELIXER_H5
    if not HELIXER_INPUT_H5.exists() or not HELIXER_PRED_H5.exists():
        sys.exit(f"ERROR: need both {HELIXER_INPUT_H5.name} and "
                 f"{HELIXER_PRED_H5.name} to build the combined Helixer HDF5.")
    with h5py.File(HELIXER_INPUT_H5, "r") as inp, h5py.File(HELIXER_H5, "w") as comb:
        inp.copy("data/seqids", comb, name="seqids")
        inp.copy("data/start_ends", comb, name="start_ends")
        comb["predictions"] = h5py.ExternalLink(
            str(HELIXER_PRED_H5.resolve()), "predictions")
    print(f"   built combined HDF5 view: {HELIXER_H5.name}")
    return HELIXER_H5


def merge_junctions(sj_paths):
    """Union STAR SJ.out.tab junctions across samples (sum reads, count samples)."""
    agg = {}  # (seqid, donor, acceptor, strand) -> [reads, samples]
    for p in sj_paths:
        for j in parse_star_sj_tab(str(p)):
            key = (j.seqid, j.donor, j.acceptor, j.strand)
            if key in agg:
                agg[key][0] += j.read_count
                agg[key][1] += 1
            else:
                agg[key] = [j.read_count, 1]
    out = [
        SpliceJunction(seqid=s, donor=d, acceptor=a, strand=st,
                       read_count=reads, samples=samples)
        for (s, d, a, st), (reads, samples) in agg.items()
    ]
    out.sort(key=lambda j: (j.seqid, j.donor, j.acceptor))
    return out


def _structure_key(seqid, strand, exons):
    """Replicates StringTieParser._structure_key so prepared tx can look up TPM."""
    exon_tuples = tuple((e.start, e.end) for e in exons)
    return f"{seqid}:{strand}:{exon_tuples}"


def parse_prepared_gtf(gtf_path):
    """Minimal GTF reader for mikado_prepared.gtf.

    Returns ``{tid: (seqid, strand, [Exon...])}`` in internal 0-based half-open
    coordinates. Unlike StringTieParser this keeps every transcript (no TPM
    filter) — every prepared transcript needs an external-scores row.
    """
    records = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            seqid, _src, _feat, start_s, end_s, _sc, strand, _fr, attrs = cols[:9]
            tid = None
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("transcript_id"):
                    tid = field.split('"')[1] if '"' in field else field.split()[-1]
                    break
            if tid is None:
                continue
            rec = records.setdefault(tid, {"seqid": seqid, "strand": strand, "exons": []})
            rec["exons"].append(Exon(int(start_s) - 1, int(end_s)))  # 1-based -> internal
    out = {}
    for tid, rec in records.items():
        exons = sorted(rec["exons"], key=lambda e: e.start)
        out[tid] = (rec["seqid"], rec["strand"], exons)
    return out


def build_external_scores(prepared_gtf, h5_reader, struct_tpm, global_max_tpm):
    """One external-scores row per prepared transcript (all three metrics)."""
    prepared = parse_prepared_gtf(prepared_gtf)
    rows = {}
    for tid, (seqid, strand, exons) in prepared.items():
        hs = helixer_support(exons, seqid, strand, h5_reader)
        start, end = exons[0].start, exons[-1].end
        hlc = helixer_locus_conf(seqid, start, end, h5_reader)
        raw_tpm = struct_tpm.get(_structure_key(seqid, strand, exons), 0.0)
        tpm = normalize_tpm(raw_tpm, global_max_tpm)
        rows[tid] = {"helixer_support": hs, "helixer_locus_conf": hlc, "tpm": tpm}
    return rows


def patch_config_as_knobs(config_path):
    """Inject pick.alternative_splicing knobs into a mikado-configure YAML.

    ``mikado configure`` writes a complete valid config but with default AS
    settings; the §7 knobs (esp. ``pad``) must be set afterward. Uses PyYAML if
    available (Mikado's own env ships it); warns and skips otherwise.
    """
    try:
        import yaml
    except ImportError:
        print("  ! PyYAML unavailable — skipping AS-knob patch; using mikado "
              "configure defaults. Set pick.alternative_splicing manually for "
              "the strict/permissive benchmark.", file=sys.stderr)
        return
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    cfg.setdefault("pick", {}).setdefault("alternative_splicing", {}).update(PICK_AS_KNOBS)
    with open(config_path, "w") as fh:
        yaml.safe_dump(cfg, fh, default_flow_style=False, sort_keys=False)
    print(f"  patched AS knobs into {Path(config_path).name}: {PICK_AS_KNOBS}")


def check_tools(need_diamond=True, need_transdecoder=True):
    """Verify the external binaries for the steps that will actually run.

    In --resume mode the caller passes need_* = False for any step whose output
    already exists (e.g. diamond when the XML is present), so those tools need
    not be installed in the active env.
    """
    missing = []
    if shutil.which(MIKADO_BIN) is None:           # serialise + pick always run
        missing.append(MIKADO_BIN)
    if need_diamond and shutil.which(DIAMOND_BIN) is None:
        missing.append(DIAMOND_BIN)
    if need_transdecoder:
        td = "TransDecoder.LongOrfs"
        td_path = (Path(TRANSDECODER_BINDIR) / td) if TRANSDECODER_BINDIR else td
        if shutil.which(str(td_path)) is None:
            missing.append(str(td_path))
    if missing:
        print("ERROR: external tools not found on PATH: " + ", ".join(missing),
              file=sys.stderr)
        print("Activate your Mikado conda env / container, or set MIKADO_BIN / "
              "DIAMOND_BIN / TRANSDECODER_BIN_DIR.", file=sys.stderr)
        sys.exit(1)


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--region", default=None,
                    help="restrict to one scaffold, e.g. 'Chr1' or 'Chr1:0-2000000' "
                         "(internal 0-based coords). Subsets Helixer loci only; "
                         "Mikado still runs genome-wide unless you also subset inputs.")
    ap.add_argument("--outdir", default=str(REPO / "m2_out"))
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-mikado", action="store_true",
                    help="stop after emitting Mikado inputs (PREP only)")
    ap.add_argument("--resume", action="store_true",
                    help="skip prepare/external/TransDecoder/DIAMOND steps whose "
                         "outputs already exist under --outdir (resume a run)")
    args = ap.parse_args()

    for p in (GENOME, HELIXER_GFF, HELIXER_INPUT_H5, HELIXER_PRED_H5, PROTEIN_FA):
        if not p.exists():
            sys.exit(f"ERROR: missing input: {p}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mik_in = outdir / "mikado_inputs"
    mik_in.mkdir(exist_ok=True)
    mik_run = outdir / "mikado_run"
    mik_run.mkdir(exist_ok=True)

    if not args.skip_mikado:
        # In --resume mode, a step's tool isn't needed if its output exists.
        td_done = (mik_run / "mikado_prepared.fasta.transdecoder.bed").exists()
        xml_done = (mik_run / "mikado_diamond.xml").exists()
        check_tools(
            need_diamond=not (args.resume and xml_done),
            need_transdecoder=not (args.resume and td_done),
        )

    # ── Phase 1-3: Helixer loci + StringTie + classification ────────────────
    print(">> Building combined Helixer HDF5 view...")
    ensure_combined_h5()
    print(">> Loading Helixer loci (+ HDF5 confidence)...")
    helixer_loci = load_helixer_loci(str(HELIXER_GFF), h5_path=str(HELIXER_H5))
    if args.region:
        if ":" in args.region:
            seqid, span = args.region.split(":", 1)
            lo, hi = span.split("-")
            helixer_loci = filter_loci_by_region(helixer_loci, seqid, int(lo), int(hi))
        else:
            helixer_loci = [g for g in helixer_loci if g.seqid == args.region]
    print(f"   {len(helixer_loci)} Helixer loci")

    print(">> Parsing StringTie assemblies...")
    st_parser = StringTieParser()
    per_sample = {}            # sample_id -> [StringTieTranscript]
    all_st = []
    for gtf in sorted(STRINGTIE_DIR.glob("*.gtf")):
        sid = sample_id_from(gtf, "stringtie")
        tx = st_parser.parse_gtf(str(gtf), sid)
        per_sample[sid] = tx
        all_st.extend(tx)
    st_agg = st_parser.aggregate_across_samples(all_st)
    print(f"   {len(all_st)} transcripts / {len(per_sample)} samples / "
          f"{len(st_agg)} unique structures")

    # structure -> max TPM, for the prepared-transcript external TPM metric
    struct_tpm = {k: v["max_tpm"] for k, v in st_agg.items()}
    global_max_tpm = max(struct_tpm.values(), default=0.0)

    print(">> Classifying loci (EXPRESSED / LOW / SILENT)...")
    # Coverage fallback for loci with no StringTie overlap: prefer bigWig (fast)
    # if pyBigWig is installed, else fall back to the BAMs (pysam, slower).
    import importlib.util
    if importlib.util.find_spec("pyBigWig") is not None:
        cov_kwargs = {"bigwig_paths": [str(p) for p in sorted(BIGWIG_DIR.glob("*.bw"))]}
        print(f"   coverage fallback: {len(cov_kwargs['bigwig_paths'])} bigWigs")
    else:
        cov_kwargs = {"bam_paths": [str(p) for p in sorted(BAM_DIR.glob("*.bam"))]}
        print(f"   coverage fallback: {len(cov_kwargs['bam_paths'])} BAMs "
              "(install pyBigWig for a faster bigWig fallback)")
    classifications = classify_loci(
        helixer_loci,
        stringtie_aggregated=st_agg,
        stringtie_transcripts=all_st,
        **cov_kwargs,
    )
    status_counts = Counter(c.status for c in classifications)
    print(f"   {dict(status_counts)}")

    # ── Phase 4: emit Mikado inputs ─────────────────────────────────────────
    print(">> Emitting Mikado inputs...")
    helixer_gtf = helixer_gff3_to_gtf(str(HELIXER_GFF), str(mik_in / "helixer.gtf"))

    st_entries = []
    for sid, tx in per_sample.items():
        if not tx:
            continue
        out_gtf = mik_in / f"stringtie_{sid}.gtf"
        stringtie_to_labelled_gtf(tx, str(out_gtf), sid)
        st_entries.append((str(out_gtf), sid))

    print("   merging STAR junctions...")
    junctions = merge_junctions(sorted(SJ_DIR.glob("*_SJ.out.tab")))
    junctions_bed = junctions_to_portcullis_tab(junctions, str(mik_in / "junctions.bed"))
    print(f"   {len(junctions)} junctions")

    list_entries = [{"file": str(helixer_gtf), "label": "helixer",
                     "is_reference": True, "score": 0}]
    for gtf_path, sid in st_entries:
        list_entries.append({"file": gtf_path, "label": sid, "score": 0})
    list_path = write_input_list(list_entries, mik_in / "list.txt")

    scoring_path = install_scoring_profile(SCORING_PROFILE, mik_in)

    if args.skip_mikado:
        print(">> --skip-mikado: stopping after PREP. Inputs are ready in "
              f"{mik_in}")
        return

    config_path = write_configuration(
        genome_fa=str(GENOME),
        list_path=str(list_path),
        scoring_path=str(scoring_path),
        junctions_path=str(junctions_bed),
        out_path=str(mik_in / "configuration.yaml"),
        use_subprocess=True,            # real `mikado configure` -> complete config
        mikado_bin=MIKADO_BIN,
    )
    patch_config_as_knobs(config_path)
    print(f"   inputs in {mik_in}")

    # ── Phase 5: prepare -> external -> transdecoder/diamond -> serialise -> pick
    # --resume skips any step whose output already exists (steps are expensive).
    prepared_gtf = mik_run / "mikado_prepared.gtf"
    prepared_fasta = mik_run / "mikado_prepared.fasta"
    if args.resume and prepared_gtf.exists() and prepared_fasta.exists():
        print(">> mikado prepare... [resume: reusing existing prepared files]")
    else:
        print(">> mikado prepare...")
        prepared_gtf, prepared_fasta = run_prepare(
            str(config_path), mik_run, procs=args.procs, mikado_bin=MIKADO_BIN)

    external_tsv = mik_in / "external_scores.tsv"
    if args.resume and external_tsv.exists():
        print(">> external scores... [resume: reusing existing TSV]")
    else:
        print(">> computing external scores on prepared transcripts...")
        with HDF5ConfidenceReader(str(HELIXER_H5)) as h5:
            rows = build_external_scores(prepared_gtf, h5, struct_tpm, global_max_tpm)
        external_tsv = write_external_scores_tsv(rows, external_tsv)
        print(f"   {len(rows)} external rows -> {external_tsv}")

    orfs_bed = mik_run / f"{prepared_fasta.name}.transdecoder.bed"
    if args.resume and orfs_bed.exists():
        print(">> TransDecoder... [resume: reusing existing ORFs bed]")
    else:
        print(">> TransDecoder...")
        orfs_bed = run_transdecoder(prepared_fasta, mik_run,
                                    transdecoder_bin_dir=TRANSDECODER_BINDIR)

    blast_out = mik_run / "mikado_diamond.xml"
    if args.resume and blast_out.exists():
        print(">> DIAMOND blastx... [resume: reusing existing XML]")
    else:
        print(">> DIAMOND blastx...")
        blast_out = run_diamond(prepared_fasta, str(PROTEIN_FA), mik_run,
                                threads=args.threads, diamond_bin=DIAMOND_BIN)
    mikado_db = mik_run / "mikado.db"
    if args.resume and mikado_db.exists():
        print(">> mikado serialise... [resume: reusing existing mikado.db]")
    else:
        print(">> mikado serialise...")
        run_serialise(str(config_path), prepared_fasta, orfs_bed, blast_out,
                      str(PROTEIN_FA), str(junctions_bed), str(external_tsv),
                      str(GENOME), mik_run, mikado_bin=MIKADO_BIN)

    loci_gff3 = mik_run / "mikado.loci.gff3"
    if args.resume and loci_gff3.exists():
        print(">> mikado pick... [resume: reusing existing loci GFF3]")
    else:
        print(">> mikado pick...")
        loci_gff3 = run_pick(str(config_path), str(scoring_path), prepared_gtf,
                             mik_run, procs=args.procs, mikado_bin=MIKADO_BIN)

    # ── Phase 5-6: parse + reconcile ────────────────────────────────────────
    print(">> parsing mikado.loci.gff3...")
    mikado_loci = parse_loci_gff3(
        str(loci_gff3),
        metrics_tsv=str(mik_run / "mikado.loci.metrics.tsv"),
        scores_tsv=str(mik_run / "mikado.loci.scores.tsv"),
    )
    print(f"   {len(mikado_loci)} Mikado loci")

    # ID-stability: reuse a persisted id_map if present (rerun-stable HFGs)
    id_map_path = outdir / "id_map.json"
    id_map = json.loads(id_map_path.read_text()) if id_map_path.exists() else None

    print(">> reconciling onto the Helixer gene set...")
    genes, id_map, admissions = reconcile(helixer_loci, classifications, mikado_loci,
                                          id_map=id_map)
    id_map_path.write_text(json.dumps(id_map, indent=0, sort_keys=True))

    # ── M2 summary ──────────────────────────────────────────────────────────
    tier = Counter(g.tier for g in genes)
    origin = Counter(g.origin for g in genes)
    flags = Counter(f.name for g in genes for f in g.flags)
    as_kinds = Counter(e.kind for g in genes for e in g.as_events)
    total_tx = sum(len(g.transcripts) for g in genes)
    multi = sum(1 for g in genes if len(g.transcripts) > 1)
    retained = sum(1 for a in admissions if a.admitted)
    dropped = len(admissions) - retained

    # Correct "no Helixer locus lost" check: every Helixer locus must have been
    # assigned an HFG (it lives in id_map). Raw gene COUNT legitimately drops
    # below the input count because merges collapse >=2 Helixer loci into one
    # gene (recorded in merged_from), so a count comparison is the wrong test.
    covered = sum(1 for h in helixer_loci if h.gene_id in id_map)
    none_lost = covered == len(helixer_loci)

    print("\n" + "=" * 60)
    print("M2 VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Helixer loci in:        {len(helixer_loci)}")
    print(f"Mikado loci:            {len(mikado_loci)}")
    print(f"Reconciled genes:       {len(genes)}  (merges collapse loci; see merged_from)")
    print(f"  Helixer loci covered: {covered}/{len(helixer_loci)}  none lost: {none_lost}")
    print(f"Total transcripts:      {total_tx}")
    print(f"Multi-isoform genes:    {multi} ({100*multi/max(len(genes),1):.1f}%)")
    print(f"Tier:                   {dict(sorted(tier.items()))}")
    print(f"Origin:                 {dict(sorted(origin.items()))}")
    print(f"AS events:              {dict(sorted(as_kinds.items()))}")
    print(f"Flags:                  {dict(sorted(flags.items()))}")
    print(f"Isoform admissions:     {retained} retained / {dropped} dropped")
    print(f"id_map persisted:       {id_map_path}")
    print("\nManual M2 checks: spot-check LOCUS_SPLIT / LOCUS_MERGE genes in IGV "
          "against RNA-seq; rerun this script to confirm HFG IDs are stable.")


if __name__ == "__main__":
    main()
