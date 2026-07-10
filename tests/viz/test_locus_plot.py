"""Phase 10 D1 — static per-locus figures.

No pixel assertions (CLAUDE.md / phase spec): we check that a non-empty Figure is
returned, files are produced, evidence is optional, and the lazy-import error
message is clear when matplotlib is absent. Both strands exercised. Concrete
literal coordinates only (CLAUDE.md §12).
"""

import builtins

import pytest

from helixforge.reconcile.models import (
    ASEvent,
    CDSSegment,
    Exon,
    HelixerLocus,
    LocusClassification,
    ReconciledGene,
    SpliceJunction,
    TranscriptCandidate,
)
from helixforge.viz.locus_plot import plot_loci, plot_locus


def make_tx(tid, strand="+", exon_bounds=((1000, 1200), (1300, 1500), (1600, 1800)),
            cds=((1050, 1200, 0), (1300, 1500, 0), (1600, 1649, 1)), is_primary=True):
    return TranscriptCandidate(
        transcript_id=tid,
        locus_id="HFG_00001",
        source="mikado",
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        exons=[Exon(s, e) for s, e in exon_bounds],
        cds=[CDSSegment(*c) for c in cds] if cds else None,
        tpm=10.0,
        junction_support_fraction=1.0,
        confidence=0.8,
        is_primary=is_primary,
    )


def make_gene(strand="+", txs=None, as_events=None):
    txs = txs or [make_tx("HFG_00001.1", strand=strand)]
    return ReconciledGene(
        gene_id="HFG_00001",
        seqid="chr1",
        start=min(t.start for t in txs),
        end=max(t.end for t in txs),
        strand=strand,
        tier=1,
        transcripts=txs,
        primary_transcript_id=txs[0].transcript_id,
        classification=LocusClassification(locus_id="HFG_00001", status="EXPRESSED"),
        origin="mikado_1to1",
        as_events=as_events or [ASEvent("ES", "chr1", 1300, 1500, "+", 20)],
    )


class FakeH5:
    """Minimal HDF5ConfidenceReader stand-in returning a flat (L,4) array."""

    seqids = ["chr1"]

    def get_per_base_predictions(self, seqid, start, end):
        import numpy as np

        a = np.zeros((end - start, 4))
        a[:, 2] = 0.7  # CDS channel
        return a


class FakeCoverage:
    def region_coverage_array(self, seqid, start, end):
        return [5] * (end - start)


def test_plot_locus_returns_figure():
    fig = plot_locus(make_gene())
    assert fig is not None
    assert len(fig.axes) >= 1


def test_plot_locus_minus_strand():
    fig = plot_locus(make_gene(strand="-",
                               txs=[make_tx("HFG_00001.1", strand="-")]))
    assert len(fig.axes) >= 1


def test_plot_locus_multi_isoform():
    g = make_gene(txs=[
        make_tx("HFG_00001.1"),
        make_tx("HFG_00001.2", exon_bounds=((1000, 1200), (1600, 1800)),
                cds=((1050, 1200, 0), (1600, 1648, 0)), is_primary=False),
    ])
    fig = plot_locus(g)
    assert len(fig.axes) >= 1


def test_plot_locus_with_all_evidence():
    junctions = [
        SpliceJunction("chr1", 1200, 1300, "+", 30),
        SpliceJunction("chr1", 1500, 1600, "+", 12),
    ]
    fig = plot_locus(make_gene(), junctions=junctions, coverage=FakeCoverage(),
                     h5_reader=FakeH5())
    assert fig is not None


def test_plot_locus_coverage_as_list():
    fig = plot_locus(make_gene(), coverage=[1, 2, 3, 4, 5])
    assert fig is not None


def test_plot_locus_handles_missing_evidence():
    # No junctions/coverage/h5 — must still return a figure.
    fig = plot_locus(make_gene(), junctions=None, coverage=None, h5_reader=None)
    assert fig is not None


def test_plot_locus_before_after_panel():
    before = HelixerLocus(
        gene_id="HELIXER_1", seqid="chr1", start=1000, end=1800, strand="+",
        exons=[Exon(1000, 1500), Exon(1600, 1800)],
        cds=[CDSSegment(1050, 1500, 0), CDSSegment(1600, 1649, 1)],
    )
    fig = plot_locus(make_gene(), before=before)
    # Before/after panel → at least two axes (Helixer + reconciled).
    assert len(fig.axes) >= 2


def test_plot_locus_before_reconciled_gene():
    before = make_gene()  # ReconciledGene as the "before"
    fig = plot_locus(make_gene(), before=before)
    assert len(fig.axes) >= 2


def test_plot_locus_saves_file(tmp_path):
    out = tmp_path / "locus.svg"
    plot_locus(make_gene(), out_path=out)
    assert out.exists() and out.stat().st_size > 0


def test_plot_loci_batch(tmp_path):
    genes = [make_gene()]
    paths = plot_loci(genes, tmp_path, fmt="svg")
    assert len(paths) == 1
    assert paths[0].exists()


def test_plot_loci_top_n_by_complexity(tmp_path):
    simple = make_gene(as_events=[])
    simple = ReconciledGene(
        gene_id="HFG_00002", seqid="chr1", start=1000, end=1800, strand="+",
        tier=1, transcripts=[make_tx("HFG_00002.1")],
        primary_transcript_id="HFG_00002.1",
        classification=LocusClassification(locus_id="HFG_00002", status="EXPRESSED"),
        origin="mikado_1to1", as_events=[],
    )
    complex_gene = make_gene(as_events=[
        ASEvent("ES", "chr1", 1300, 1500, "+"),
        ASEvent("IR", "chr1", 1200, 1300, "+"),
    ])
    paths = plot_loci([simple, complex_gene], tmp_path, top_n_by_as_complexity=1)
    assert len(paths) == 1
    # The AS-complex gene (HFG_00001) is the one kept.
    assert paths[0].name == "HFG_00001.svg"


def test_plot_locus_lazy_import_error(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "matplotlib" or name.startswith("matplotlib"):
            raise ImportError("no matplotlib")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match=r"helixforge\[viz\]"):
        plot_locus(make_gene())
