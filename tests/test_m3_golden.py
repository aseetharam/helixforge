"""Golden-output regression gate (CLAUDE.md §15, regression-gate clause 3).

Re-runs the validated M3 resume pipeline on the real *A. thaliana* set and
asserts its gene / tier / origin / AS-event / flag counts are **identical** to
the recorded snapshot in ``tests/golden/m3_counts.json``. Any phase touching
``reconcile/``, ``mikado/``, ``io/``, or ``pipeline.py`` must keep these counts
unchanged; a 3 bp or ID change that moves them is a regression unless intended.

This test is **opt-in and heavy**: it needs the gitignored real data
(``helixforge_testing/`` + the M2 Mikado artifacts under ``m2_out/``) and runs
the full reconciliation, so it is skipped unless ``HELIXFORGE_GOLDEN=1`` is set
*and* the data is present. The default ``pytest tests/ -q`` suite stays fast and
green everywhere. Run the gate with:

    HELIXFORGE_GOLDEN=1 pytest tests/test_m3_golden.py -v

If a deliberate, documented change is meant to move the counts, regenerate the
snapshot via ``/usr/bin/python3 scripts/m3_run.py --outdir m3_out`` and update
``tests/golden/m3_counts.json`` in the same commit.
"""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
GOLDEN = REPO / "tests" / "golden" / "m3_counts.json"
M3_SCRIPT = REPO / "scripts" / "m3_run.py"

# Real-data prerequisites (all gitignored). Mirror the paths m3_run.py uses.
DATA = REPO / "helixforge_testing"
REQUIRED = [
    DATA / "genome" / "athaliana.fasta",
    DATA / "helixer_output" / "Arabidopsis-thaliana_helixer.gff3",
    DATA / "helixer_output" / "Arabidopsis-thaliana_input.h5",
    DATA / "helixer_output" / "Arabidopsis-thaliana_predictions.h5",
    DATA / "swissprot" / "uniprot_sprot.fasta",
    REPO / "m2_out" / "mikado_run" / "mikado.loci.gff3",
    REPO / "m2_out" / "id_map.json",
]

_OPT_IN = os.environ.get("HELIXFORGE_GOLDEN") == "1"
_DATA_PRESENT = all(p.exists() for p in REQUIRED)

pytestmark = [
    pytest.mark.skipif(not _OPT_IN, reason="set HELIXFORGE_GOLDEN=1 to run the M3 golden gate"),
    pytest.mark.skipif(not _DATA_PRESENT, reason="real M3 data (helixforge_testing/ + m2_out/) not present"),
]


def _load_m3_module():
    """Import scripts/m3_run.py as a module (it is not on the package path)."""
    spec = importlib.util.spec_from_file_location("m3_run", M3_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def m3_summary(tmp_path_factory):
    """Run the M3 resume pipeline once and return (summary, id_map_loci)."""
    from helixforge.reconcile.pipeline import run_pipeline
    from helixforge.stats.summary import summarize_genes

    m3 = _load_m3_module()
    outdir = tmp_path_factory.mktemp("m3_golden")
    config = m3.build_config(
        outdir,
        reuse_dir=REPO / "m2_out" / "mikado_run",
        id_map=REPO / "m2_out" / "id_map.json",
        resume=True,
    )
    genes = run_pipeline(config)
    summary = summarize_genes(genes)
    id_map = json.loads((REPO / "m2_out" / "id_map.json").read_text())
    return summary, len(id_map), genes


@pytest.fixture(scope="module")
def golden():
    return json.loads(GOLDEN.read_text())


def test_m3_summary_matches_golden(m3_summary, golden):
    summary, _, _ = m3_summary
    # JSON keys are strings; the live tier dict has int keys — normalize both.
    live = dict(summary)
    live["tier"] = {str(k): v for k, v in live["tier"].items()}
    assert live == golden["summary"]


def test_m3_gene_and_tier_counts(m3_summary):
    """Spot-check the headline numbers with concrete literals (CLAUDE.md §12)."""
    summary, _, _ = m3_summary
    # TPM/biotype/TE fix re-baseline: biotype is driven by the model's own ORF.
    # The 836 helixer_backstop genes carry their intrinsic Helixer CDS, turning
    # most protein_coding (Tier 2, no homology). The Mikado-compare CDS-coherence
    # fix then rejects the intrinsic CDS of the 85 backstops whose junction
    # correction stranded an internal CDS boundary off its splice site (a CDS
    # Mikado's finalizer cannot reconcile): coding 26999->26914, tier 2
    # 6368->6283, tier 3 187->272, biotype protein_coding 26999->26914
    # (lncRNA 145->204, ncRNA_undetermined 42->68). genes/origin/tier1 unchanged.
    assert summary["genes"] == 27186
    assert summary["tier"] == {1: 20631, 2: 6283, 3: 272}
    assert summary["origin"]["helixer_backstop"] == 836
    # A coding gene with a complete ORF is protein_coding regardless of evidence;
    # only genuinely ORF-less loci (incl. the 85 rejected incoherent backstops)
    # stay non-coding; TAIR10 has 0 pseudogenes.
    assert summary["coding"] == 26914
    assert summary["biotype"] == {
        "lncRNA": 204, "ncRNA_undetermined": 68, "protein_coding": 26914,
    }


def test_m3_no_gene_lost(m3_summary, golden):
    """No Helixer locus lost: id_map coverage matches the recorded count."""
    _, id_map_loci, _ = m3_summary
    assert id_map_loci == 27198
    assert id_map_loci == golden["id_map_loci"]


@pytest.fixture(scope="module")
def m3_trace_runs(tmp_path_factory):
    """Run the M3 pipeline twice with ``trace_primary=True`` (Phase 33b golden variant).

    Returns ``(summary, [genes_run1, genes_run2])``. The TRaCE election reorders
    isoforms (changing ``.N`` ids / primary) but must leave the gene / tier /
    origin / biotype / AS-event-set counts identical to the default golden run,
    keep global id-uniqueness, and be reproducible across invocations.
    """
    import dataclasses

    from helixforge.reconcile.pipeline import run_pipeline
    from helixforge.stats.summary import summarize_genes

    m3 = _load_m3_module()
    runs = []
    for i in (1, 2):
        outdir = tmp_path_factory.mktemp(f"m3_trace_{i}")
        config = m3.build_config(
            outdir,
            reuse_dir=REPO / "m2_out" / "mikado_run",
            id_map=REPO / "m2_out" / "id_map.json",
            resume=True,
        )
        config = dataclasses.replace(config, trace_primary=True)
        runs.append(run_pipeline(config))
    return summarize_genes(runs[0]), runs


def test_m3_trace_on_summary_unchanged(m3_trace_runs, golden):
    """TRaCE-on keeps the gene/tier/origin/AS-event-set counts identical."""
    summary, _ = m3_trace_runs
    live = dict(summary)
    live["tier"] = {str(k): v for k, v in live["tier"].items()}
    assert live == golden["summary"]


def test_m3_trace_on_ids_unique(m3_trace_runs):
    """TRaCE-on still emits globally-unique gene + transcript IDs (CLAUDE.md §11)."""
    _, runs = m3_trace_runs
    genes = runs[0]
    gene_ids = [g.gene_id for g in genes]
    tx_ids = [t.transcript_id for g in genes for t in g.transcripts]
    assert len(set(gene_ids)) == len(gene_ids), "duplicate gene IDs under TRaCE"
    assert len(set(tx_ids)) == len(tx_ids), "duplicate transcript IDs under TRaCE"


def test_m3_trace_on_is_reproducible(m3_trace_runs):
    """Two TRaCE-on invocations elect the same primary + ordering for every gene."""
    _, runs = m3_trace_runs

    def fingerprint(genes):
        return {
            g.gene_id: (
                g.primary_transcript_id,
                tuple((t.transcript_id, t.trace_rank) for t in g.transcripts),
            )
            for g in genes
        }

    assert fingerprint(runs[0]) == fingerprint(runs[1])


def test_m3_unique_gene_and_transcript_ids(m3_summary):
    """Every emitted gene/transcript ID is globally unique (CLAUDE.md §11).

    Regression guard for the MERGE_REJECTED duplicate-HFG bug that crashed
    `scorer:stats` on TAIR10 (assessment-v4.md §2.3): two rejected merges each
    released two backstop loci under one shared HFG, so the GFF3 carried
    duplicate gene + transcript IDs. The in-memory golden counts still matched
    (they count objects, not id strings), so only an explicit id-uniqueness
    assertion catches it.
    """
    _, _, genes = m3_summary
    gene_ids = [g.gene_id for g in genes]
    tx_ids = [t.transcript_id for g in genes for t in g.transcripts]
    dup_genes = {x for x in gene_ids if gene_ids.count(x) > 1}
    assert not dup_genes, f"duplicate gene IDs emitted: {sorted(dup_genes)}"
    assert len(set(tx_ids)) == len(tx_ids), "duplicate transcript IDs emitted"
