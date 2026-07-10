"""Tests for prep skip-if-exists idempotency (Phase 22 D3). Floor: 3.

A re-run of a long maize-scale prep should skip already-completed work: a wrapper
whose output exists, is non-empty, and is newer than its inputs is skipped unless
``force=True``. Tools are fully mocked (no real binary runs).
"""

import os
import subprocess

import pytest

from helixforge.prep import assemble as asm
from helixforge.prep import protein_align as pa
from helixforge.prep import samtools as st
from helixforge.prep._subprocess import output_is_fresh


@pytest.fixture
def calls(monkeypatch):
    recorded = []

    def fake_run_tool(argv, **kw):
        recorded.append([str(a) for a in argv])
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    # Patch the shared name in each wrapper module.
    for mod in (asm, pa, st):
        monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    return recorded


def _older(path, ref):
    """Make ``path`` older than ``ref`` by backdating its mtime."""
    ref_mtime = os.stat(ref).st_mtime
    os.utime(path, (ref_mtime - 100, ref_mtime - 100))


# --- output_is_fresh helper -------------------------------------------------

def test_output_is_fresh_missing_or_empty(tmp_path):
    missing = tmp_path / "nope.gtf"
    assert output_is_fresh(missing, []) is False
    empty = tmp_path / "empty.gtf"
    empty.write_text("")
    assert output_is_fresh(empty, []) is False


def test_output_is_fresh_newer_and_stale(tmp_path):
    inp = tmp_path / "in.bam"
    inp.write_text("reads")
    out = tmp_path / "out.gtf"
    out.write_text("assembled")
    # out written after in ⇒ fresh.
    assert output_is_fresh(out, [inp]) is True
    # backdate out before in ⇒ stale.
    _older(out, inp)
    assert output_is_fresh(out, [inp]) is False


# --- existing output is skipped ---------------------------------------------

def test_existing_output_skips_run(calls, tmp_path):
    bam = tmp_path / "s.bam"
    bam.write_text("reads")
    out_gtf = tmp_path / "s.gtf"
    out_gtf.write_text("already-assembled")  # newer than bam
    result = asm.run_stringtie(str(bam), out_gtf, "sampA")
    assert calls == []                       # tool NOT invoked
    assert str(result) == str(out_gtf)


def test_force_reruns_even_if_output_exists(calls, tmp_path):
    bam = tmp_path / "s.bam"
    bam.write_text("reads")
    out_gtf = tmp_path / "s.gtf"
    out_gtf.write_text("already-assembled")
    asm.run_stringtie(str(bam), out_gtf, "sampA", force=True)
    assert len(calls) == 1                    # force overrides the skip
    assert calls[0][0] == "stringtie"


def test_stale_output_reruns(calls, tmp_path):
    bam = tmp_path / "s.bam"
    bam.write_text("reads")
    out_gtf = tmp_path / "s.gtf"
    out_gtf.write_text("old-assembly")
    _older(out_gtf, bam)                       # out older than its input bam
    asm.run_stringtie(str(bam), out_gtf, "sampA")
    assert len(calls) == 1                     # stale ⇒ re-run


def test_miniprot_skips_when_fresh(calls, tmp_path):
    genome = tmp_path / "g.fa"
    genome.write_text(">c\nACGT\n")
    prot = tmp_path / "p.fa"
    prot.write_text(">p\nMK\n")
    out = tmp_path / "mp.gff"
    out.write_text("##gff-version 3\n")        # fresh result already present
    pa.run_miniprot(str(genome), str(prot), out)
    assert calls == []
    # force re-runs.
    pa.run_miniprot(str(genome), str(prot), out, force=True)
    assert len(calls) == 1


def test_samtools_faidx_skips_when_fresh(calls, tmp_path):
    fasta = tmp_path / "g.fa"
    fasta.write_text(">c\nACGT\n")
    fai = tmp_path / "g.fa.fai"
    fai.write_text("c\t4\t3\t4\t5\n")          # fresh .fai present
    st.faidx(str(fasta))
    assert calls == []
    st.faidx(str(fasta), force=True)
    assert len(calls) == 1
