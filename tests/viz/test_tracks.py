"""Phase 10 D3 — genome-browser track export.

BED12 column correctness (both strands, thick=CDS), junction BED scores, the
bigWig writer (pyBigWig mocked — it is not a core dep), and session/config files
referencing the right paths. Concrete literal coordinates only (CLAUDE.md §12).
"""

import builtins
import json
import sys
import types

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    SpliceJunction,
    TranscriptCandidate,
)
from helixforge.viz.tracks import (
    write_bed12,
    write_bigbed,
    write_confidence_bigwig,
    write_igv_session,
    write_jbrowse_config,
    write_junction_bed,
)


def make_tx(tid, strand="+", exon_bounds=((1000, 1200), (1300, 1500), (1600, 1800)),
            cds=((1050, 1200, 0), (1300, 1500, 0), (1600, 1649, 1))):
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
    )


def make_gene(gid="HFG_00001", strand="+", tier=1, txs=None):
    txs = txs or [make_tx("HFG_00001.1", strand=strand)]
    return ReconciledGene(
        gene_id=gid,
        seqid="chr1",
        start=min(t.start for t in txs),
        end=max(t.end for t in txs),
        strand=strand,
        tier=tier,
        transcripts=txs,
        primary_transcript_id=txs[0].transcript_id,
        classification=LocusClassification(locus_id=gid, status="EXPRESSED"),
        origin="mikado_1to1",
    )


# ---------------------------------------------------------------------------
# BED12
# ---------------------------------------------------------------------------


def test_bed12_columns_plus(tmp_path):
    out = write_bed12([make_gene()], tmp_path / "models.bed")
    cols = out.read_text().splitlines()[0].split("\t")
    assert len(cols) == 12
    assert cols[0] == "chr1"
    assert cols[1] == "1000"          # chromStart (0-based)
    assert cols[2] == "1800"          # chromEnd
    assert cols[3] == "HFG_00001.1"   # name
    assert cols[5] == "+"             # strand
    assert cols[6] == "1050"          # thickStart = CDS start
    assert cols[7] == "1649"          # thickEnd = CDS end
    assert cols[9] == "3"             # blockCount
    assert cols[10] == "200,200,200,"  # blockSizes
    assert cols[11] == "0,300,600,"    # blockStarts (relative to chromStart)


def test_bed12_minus_strand(tmp_path):
    out = write_bed12([make_gene(strand="-")], tmp_path / "models.bed")
    cols = out.read_text().splitlines()[0].split("\t")
    assert cols[5] == "-"
    # Block order is always genomic low→high regardless of strand.
    assert cols[10] == "200,200,200,"
    assert cols[11] == "0,300,600,"


def test_bed12_noncoding_thick_collapses(tmp_path):
    tx = make_tx("HFG_00001.1", exon_bounds=((1000, 1800),), cds=None)
    out = write_bed12([make_gene(txs=[tx])], tmp_path / "models.bed")
    cols = out.read_text().splitlines()[0].split("\t")
    # No CDS → thickStart == thickEnd == chromStart.
    assert cols[6] == "1000"
    assert cols[7] == "1000"


def test_bed12_sorted_by_position(tmp_path):
    late = make_gene("HFG_00009", txs=[make_tx("HFG_00009.1",
                     exon_bounds=((9000, 9200),), cds=None)])
    early = make_gene("HFG_00001")
    out = write_bed12([late, early], tmp_path / "models.bed")
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[1] == "1000"
    assert lines[-1].split("\t")[1] == "9000"


def test_bed12_score_by_tier(tmp_path):
    out = write_bed12([make_gene(tier=3)], tmp_path / "models.bed")
    cols = out.read_text().splitlines()[0].split("\t")
    assert cols[4] == "500"  # tier 3 → score 500


# ---------------------------------------------------------------------------
# junction BED
# ---------------------------------------------------------------------------


def test_junction_bed_scores_and_coords(tmp_path):
    junctions = [
        SpliceJunction("chr1", 1200, 1300, "+", 42),
        SpliceJunction("chr1", 1500, 1600, "-", 5000),  # capped at 1000
    ]
    out = write_junction_bed(junctions, tmp_path / "junc.bed")
    rows = [ln.split("\t") for ln in out.read_text().splitlines()]
    assert rows[0][1] == "1200" and rows[0][2] == "1300"
    assert rows[0][4] == "42"            # read count score
    assert rows[1][4] == "1000"          # capped
    assert rows[1][5] == "-"
    # anchor-style two-block intron representation
    assert rows[0][9] == "2"
    assert rows[0][10] == "1,1,"
    assert rows[0][11] == "0,99,"        # acceptor-donor-1 = 99


def test_junction_bed_sorted(tmp_path):
    junctions = [
        SpliceJunction("chr1", 1500, 1600, "+", 10),
        SpliceJunction("chr1", 1200, 1300, "+", 10),
    ]
    out = write_junction_bed(junctions, tmp_path / "junc.bed")
    rows = [ln.split("\t") for ln in out.read_text().splitlines()]
    assert rows[0][1] == "1200"
    assert rows[1][1] == "1500"


# ---------------------------------------------------------------------------
# bigWig (pyBigWig mocked)
# ---------------------------------------------------------------------------


class _FakeBigWig:
    def __init__(self, path):
        self.path = path
        self.header = None
        self.entries = []

    def addHeader(self, header):
        self.header = header

    def addEntries(self, chrom, start, values=None, span=1, step=1):
        self.entries.append((chrom, start, list(values), span, step))

    def close(self):
        # Write a stub so the file exists on disk.
        with open(self.path, "w") as fh:
            fh.write("bigwig")


class FakeH5:
    seqids = ["chr1"]

    def get_per_base_predictions(self, seqid, start, end):
        import numpy as np

        a = np.zeros((end - start, 4))
        a[:, 1] = 0.3  # UTR
        a[:, 2] = 0.6  # CDS
        return a


@pytest.fixture
def fake_pybigwig(monkeypatch):
    created = {}

    def _open(path, mode):
        bw = _FakeBigWig(path)
        created["bw"] = bw
        return bw

    mod = types.ModuleType("pyBigWig")
    mod.open = _open
    monkeypatch.setitem(sys.modules, "pyBigWig", mod)
    return created


def test_confidence_bigwig_genic(tmp_path, fake_pybigwig):
    out = write_confidence_bigwig(FakeH5(), {"chr1": 50}, tmp_path / "conf.bw")
    assert out.exists()
    bw = fake_pybigwig["bw"]
    assert bw.header == [("chr1", 50)]
    chrom, start, values, span, step = bw.entries[0]
    assert chrom == "chr1"
    # genic = max(CDS=0.6, UTR=0.3) = 0.6
    assert values[0] == pytest.approx(0.6)
    assert len(values) == 50


def test_confidence_bigwig_channel(tmp_path, fake_pybigwig):
    write_confidence_bigwig(FakeH5(), {"chr1": 10}, tmp_path / "utr.bw", channel="utr")
    values = fake_pybigwig["bw"].entries[0][2]
    assert values[0] == pytest.approx(0.3)


def test_confidence_bigwig_bad_channel(tmp_path, fake_pybigwig):
    with pytest.raises(ValueError, match="channel"):
        write_confidence_bigwig(FakeH5(), {"chr1": 10}, tmp_path / "x.bw", channel="bogus")


def test_confidence_bigwig_import_error(tmp_path, monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pyBigWig":
            raise ImportError("no pyBigWig")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match=r"helixforge\[tracks\]"):
        write_confidence_bigwig(FakeH5(), {"chr1": 10}, tmp_path / "x.bw")


# ---------------------------------------------------------------------------
# bigBed
# ---------------------------------------------------------------------------


def test_bigbed_missing_tool_raises(tmp_path):
    with pytest.raises(RuntimeError, match="bedToBigBed"):
        write_bigbed([make_gene()], tmp_path / "models.bb", {"chr1": 9999},
                     bedtobigbed_bin="definitely-not-a-real-binary-xyz")


def test_bigbed_invokes_tool(tmp_path, monkeypatch):
    calls = {}

    def fake_run(cmd, check, capture_output):
        calls["cmd"] = cmd
        open(cmd[3], "w").write("bb")  # the output path

        class R:
            returncode = 0
        return R()

    monkeypatch.setattr("helixforge.viz.tracks.subprocess.run", fake_run)
    out = write_bigbed([make_gene()], tmp_path / "models.bb", {"chr1": 9999})
    assert out.exists()
    assert calls["cmd"][0] == "bedToBigBed"
    # the intermediate BED + sizes files were produced as inputs.
    assert (tmp_path / "models.bed").exists()
    assert (tmp_path / "models.sizes").exists()


# ---------------------------------------------------------------------------
# session / config
# ---------------------------------------------------------------------------


def test_igv_session_references_tracks(tmp_path):
    bed = tmp_path / "models.bed"
    bed.write_text("x")
    out = write_igv_session([bed], "/genome/at.fa", tmp_path / "session.xml")
    text = out.read_text()
    assert "models.bed" in text
    assert "/genome/at.fa" in text


def test_jbrowse_config_track_types(tmp_path):
    tracks = {"models": tmp_path / "models.bed", "conf": tmp_path / "conf.bw"}
    out = write_jbrowse_config(tracks, tmp_path / "config.json",
                               genome_fa="/genome/at.fa")
    cfg = json.loads(out.read_text())
    types_by_id = {t["trackId"]: t["adapter"]["type"] for t in cfg["tracks"]}
    assert types_by_id["models"] == "BedAdapter"
    assert types_by_id["conf"] == "BigWigAdapter"
    assert cfg["assemblies"][0]["sequence"]["adapter"]["fastaLocation"]["uri"] == "/genome/at.fa"
