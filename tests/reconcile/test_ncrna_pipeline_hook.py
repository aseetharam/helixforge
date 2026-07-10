"""Defensive-fallback tests for ``_maybe_add_structured_ncrna`` (M-REALIZE-5).

The structured-ncRNA hook is a "completeness extra" (CLAUDE.md §15): a missing
or failing tRNAscan-SE / Infernal must **warn and continue**, never sink an
otherwise-finished run — mirroring the sibling ``_maybe_annotate_function``
hook. These tests pin that contract (both strands; concrete coordinates).
"""

from __future__ import annotations

import pytest

from helixforge.prep._subprocess import ToolError
from helixforge.reconcile import pipeline as pl
from helixforge.reconcile.models import (
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)


def _gene(gid: str, strand: str, seqid: str = "chr1", start: int = 1000) -> ReconciledGene:
    end = start + 500
    exons = [Exon(start, start + 200), Exon(start + 300, end)]
    tx = TranscriptCandidate(
        transcript_id=f"{gid}.1", locus_id=gid, source="helixer", seqid=seqid,
        start=start, end=end, strand=strand, exons=exons, cds=None,
        is_primary=True,
    )
    return ReconciledGene(
        gene_id=gid, seqid=seqid, start=start, end=end, strand=strand, tier=3,
        transcripts=[tx], primary_transcript_id=f"{gid}.1",
        classification=LocusClassification(locus_id=gid, status="SILENT", max_tpm=0.0),
        origin="helixer_backstop",
    )


def _config(tmp_path) -> pl.PipelineConfig:
    return pl.PipelineConfig(
        genome_fasta="genome.fa",
        helixer_gff3="helixer.gff3",
        ncrna_scan=True,
        ncrna_tool="trnascan",
        work_dir=str(tmp_path),
    )


def test_disabled_hook_is_identity(tmp_path, monkeypatch):
    # ncrna_scan False → returns genes unchanged and never imports/calls the tool.
    def boom(*a, **k):
        raise AssertionError("scan must not run when ncrna_scan is False")

    monkeypatch.setattr("helixforge.prep.ncrna.scan_structured_ncrna", boom)
    cfg = _config(tmp_path)
    cfg.ncrna_scan = False
    genes = [_gene("HFG_00001", "+"), _gene("HFG_00002", "-")]
    assert pl._maybe_add_structured_ncrna(genes, cfg) == genes


@pytest.mark.parametrize("exc", [
    FileNotFoundError("tool not found: 'tRNAscan-SE'"),
    ToolError("tRNAscan-SE exited 1"),
])
def test_missing_or_failing_tool_warns_and_returns_genes(tmp_path, monkeypatch, caplog, exc):
    # The defect: an absent/failing ncRNA tool must NOT propagate out of the hook.
    def fail(*a, **k):
        raise exc

    monkeypatch.setattr("helixforge.prep.ncrna.scan_structured_ncrna", fail)
    genes = [_gene("HFG_00001", "+"), _gene("HFG_00002", "-")]
    with caplog.at_level("WARNING"):
        out = pl._maybe_add_structured_ncrna(genes, _config(tmp_path))
    # genes returned unchanged (both strands preserved) — pipeline can complete.
    assert out == genes
    assert {g.strand for g in out} == {"+", "-"}
    assert any("structured-ncRNA hook skipped" in r.message for r in caplog.records)


def test_success_path_merges_and_sorts(tmp_path, monkeypatch):
    # When the tool yields loci they are merged and the set stays (seqid, start, end)-sorted.
    ncrna_gene = _gene("HFG_95000", "+", start=2000)  # falls between the input genes
    incoming = [_gene("HFG_00001", "+", start=100), _gene("HFG_00009", "-", start=5000)]

    monkeypatch.setattr(
        "helixforge.prep.ncrna.scan_structured_ncrna",
        lambda *a, **k: [ncrna_gene],
    )
    out = pl._maybe_add_structured_ncrna(incoming, _config(tmp_path))
    assert [g.start for g in out] == [100, 2000, 5000]
