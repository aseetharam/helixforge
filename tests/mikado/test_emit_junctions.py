"""Tests for mikado/emit_junctions.py (Phase 4). Floor: 10.

Subprocess (portcullis) is mocked — no real Portcullis.
"""

import os

import pytest

from helixforge.mikado import emit_junctions
from helixforge.mikado.emit_junctions import junctions_to_portcullis_tab, run_portcullis
from helixforge.reconcile.models import SpliceJunction


def _parse_bed(path):
    return [line.rstrip("\n").split("\t") for line in open(path) if line.strip()]


def _intron_from_bed(row):
    """Recover (donor, acceptor) from a BED12 junction row."""
    chrom_start = int(row[1])
    block_sizes = [int(x) for x in row[10].split(",")]
    block_starts = [int(x) for x in row[11].split(",")]
    donor = chrom_start + block_sizes[0]
    acceptor = chrom_start + block_starts[1]
    return donor, acceptor


def test_junctions_tab_returns_path(tmp_path):
    juncs = [SpliceJunction("chr1", 150, 250, "+", read_count=5)]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    assert out.exists()


def test_junctions_tab_row_count(tmp_path):
    juncs = [
        SpliceJunction("chr1", 150, 250, "+", read_count=5),
        SpliceJunction("chr1", 450, 550, "-", read_count=3),
    ]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    assert len(_parse_bed(out)) == 2


def test_junctions_tab_plus_recovers_intron(tmp_path):
    juncs = [SpliceJunction("chr1", 150, 250, "+", read_count=5)]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    row = _parse_bed(out)[0]
    assert _intron_from_bed(row) == (150, 250)


def test_junctions_tab_minus_recovers_intron(tmp_path):
    juncs = [SpliceJunction("chr1", 450, 550, "-", read_count=3)]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    row = _parse_bed(out)[0]
    assert _intron_from_bed(row) == (450, 550)


def test_junctions_tab_strand_column(tmp_path):
    juncs = [
        SpliceJunction("chr1", 150, 250, "+", read_count=5),
        SpliceJunction("chr1", 450, 550, "-", read_count=3),
    ]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    rows = _parse_bed(out)
    assert rows[0][5] == "+"
    assert rows[1][5] == "-"


def test_junctions_tab_score_is_read_count(tmp_path):
    juncs = [SpliceJunction("chr1", 150, 250, "+", read_count=7)]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    assert _parse_bed(out)[0][4] == "7"


def test_junctions_tab_score_capped_1000(tmp_path):
    juncs = [SpliceJunction("chr1", 150, 250, "+", read_count=5000)]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    assert _parse_bed(out)[0][4] == "1000"


def test_junctions_tab_blockcount_two(tmp_path):
    juncs = [SpliceJunction("chr1", 150, 250, "+", read_count=5)]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    assert _parse_bed(out)[0][9] == "2"


def test_junctions_tab_preserves_samples_in_name(tmp_path):
    juncs = [SpliceJunction("chr1", 150, 250, "+", read_count=5, samples=3)]
    out = junctions_to_portcullis_tab(juncs, tmp_path / "j.bed")
    assert "_s3" in _parse_bed(out)[0][3]


def test_junctions_tab_empty(tmp_path):
    out = junctions_to_portcullis_tab([], tmp_path / "j.bed")
    assert out.read_text() == ""


# --- run_portcullis (mocked) ---

def test_run_portcullis_missing_binary_raises(monkeypatch):
    monkeypatch.setattr(emit_junctions.shutil, "which", lambda b: None)
    with pytest.raises(FileNotFoundError):
        run_portcullis("genome.fa", ["a.bam"], "out")


def test_run_portcullis_assembles_argv(monkeypatch, tmp_path):
    # Phase 23 D4: run_portcullis now routes through the shared run_tool, so the
    # argv (and stderr-log streaming) match every other external tool.
    monkeypatch.setattr(emit_junctions.shutil, "which", lambda b: "/usr/bin/portcullis")
    captured = {}

    def fake_run_tool(argv, **kwargs):
        captured["argv"] = [str(a) for a in argv]
        captured["kwargs"] = kwargs

    monkeypatch.setattr(emit_junctions, "run_tool", fake_run_tool)
    out = run_portcullis("genome.fa", ["a.bam", "b.bam"], tmp_path, threads=8)
    argv = captured["argv"]
    assert argv[0] == "portcullis"
    assert argv[1] == "full"
    assert "-t" in argv and argv[argv.index("-t") + 1] == "8"
    assert any(a.endswith("genome.fa") for a in argv)
    assert all(os.path.isabs(a) for a in argv[-2:])
    assert argv[-2].endswith("a.bam") and argv[-1].endswith("b.bam")
    # full stderr is streamed to a per-step log under the output dir
    assert str(captured["kwargs"]["stderr_log"]).endswith("portcullis.stderr.log")
    assert str(out).endswith("portcullis_filtered.pass.junctions.tab")


def test_run_portcullis_custom_binary(monkeypatch, tmp_path):
    monkeypatch.setattr(emit_junctions.shutil, "which", lambda b: "/opt/portcullis")
    captured = {}
    monkeypatch.setattr(emit_junctions, "run_tool",
                        lambda argv, **k: captured.update(argv=[str(a) for a in argv]))
    run_portcullis("g.fa", ["a.bam"], tmp_path, portcullis_bin="portcullis2")
    assert captured["argv"][0] == "portcullis2"
