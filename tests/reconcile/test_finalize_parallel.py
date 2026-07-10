"""Phase 24 D2 — parallel per-gene finalize + batched TransDecoder. Floor: 4.

Parallel finalize must equal the serial path **exactly** (counts / IDs / flags /
order), independent of worker scheduling. Both strands are represented (g1 + and
g2 -). The batched TransDecoder must produce the same CDS the per-gene subprocess
did. Synthetic fixtures + concrete literal coordinates (CLAUDE.md §12); the
ThreadPoolExecutor seam observes the per-worker genome re-open in-process.
"""

from concurrent.futures import ThreadPoolExecutor

import pytest

from helixforge.reconcile import cds as cds_mod
from helixforge.reconcile import pipeline as pl
from helixforge.reconcile.classify import classify_loci
from helixforge.reconcile.mikado_integrate import IdAllocator, reconcile
from helixforge.reconcile.pipeline import PipelineConfig
from helixforge.reconcile.runstats import RunStats

# g1 chr1 + 101-200 (single exon, len 100); g2 chr1 - 301-600 (two exons,
# intron 400-500); g3 chr2 + 101-250 (single exon, second scaffold).
_HELIXER_GFF = """\
##gff-version 3
chr1\tHelixer\tgene\t101\t200\t.\t+\t.\tID=g1
chr1\tHelixer\tmRNA\t101\t200\t.\t+\t.\tID=g1.m;Parent=g1
chr1\tHelixer\texon\t101\t200\t.\t+\t.\tID=g1.e1;Parent=g1.m
chr1\tHelixer\tgene\t301\t600\t.\t-\t.\tID=g2
chr1\tHelixer\tmRNA\t301\t600\t.\t-\t.\tID=g2.m;Parent=g2
chr1\tHelixer\texon\t301\t400\t.\t-\t.\tID=g2.e1;Parent=g2.m
chr1\tHelixer\texon\t501\t600\t.\t-\t.\tID=g2.e2;Parent=g2.m
chr2\tHelixer\tgene\t101\t250\t.\t+\t.\tID=g3
chr2\tHelixer\tmRNA\t101\t250\t.\t+\t.\tID=g3.m;Parent=g3
chr2\tHelixer\texon\t101\t250\t.\t+\t.\tID=g3.e1;Parent=g3.m
"""


@pytest.fixture
def genome_fasta(tmp_path):
    fa = tmp_path / "genome.fasta"
    fa.write_text(">chr1\n" + "ACGT" * 200 + "\n>chr2\n" + "ACGT" * 100 + "\n")
    return str(fa)


@pytest.fixture
def helixer_gff3(tmp_path):
    p = tmp_path / "helixer.gff3"
    p.write_text(_HELIXER_GFF)
    return str(p)


def _config(genome_fasta, helixer_gff3, tmp_path, **kw):
    return PipelineConfig(
        genome_fasta=genome_fasta, helixer_gff3=helixer_gff3,
        output_prefix=str(tmp_path / "out"), work_dir=str(tmp_path / "work"), **kw,
    )


def _genes_and_inputs(config):
    """Reconcile (basic mode → all helixer_backstop) → (genes, open inputs)."""
    inputs = pl._load_inputs(config)
    loci = inputs["helixer_loci"]
    cls = classify_loci(loci)
    genes, _id_map, _adm = reconcile(loci, cls, [], id_map=None,
                                     allocator=IdAllocator())
    return genes, inputs


def _summary(genes):
    return [
        (g.gene_id, g.seqid, g.start, g.end, g.strand, g.tier, g.origin,
         tuple(f.name for f in g.flags))
        for g in genes
    ]


# ---------------------------------------------------------------------------
# Parallel == serial
# ---------------------------------------------------------------------------

def test_parallel_finalize_equals_serial(genome_fasta, helixer_gff3, tmp_path):
    serial_cfg = _config(genome_fasta, helixer_gff3, tmp_path, finalize_workers=1)
    par_cfg = _config(genome_fasta, helixer_gff3, tmp_path, finalize_workers=3)

    genes, inputs = _genes_and_inputs(serial_cfg)
    try:
        serial = pl._finalize_genes(genes, serial_cfg, inputs)
        parallel = pl._finalize_genes(
            genes, par_cfg, inputs, _executor_cls=ThreadPoolExecutor
        )
    finally:
        inputs["genome"].close()

    # Identical genes, IDs, flags, AND order — worker scheduling changes nothing.
    assert _summary(parallel) == _summary(serial)
    # Both strands present (the minus-strand g2 must survive identically).
    assert any(s[4] == "-" for s in _summary(parallel))


def test_parallel_finalize_three_genes_two_scaffolds(genome_fasta, helixer_gff3,
                                                      tmp_path):
    par_cfg = _config(genome_fasta, helixer_gff3, tmp_path, finalize_workers=2)
    genes, inputs = _genes_and_inputs(par_cfg)
    try:
        out = pl._finalize_genes(genes, par_cfg, inputs,
                                 _executor_cls=ThreadPoolExecutor)
    finally:
        inputs["genome"].close()
    ids = [g.gene_id for g in out]
    assert ids == sorted(ids)               # genomic order preserved
    assert len(out) == 3                     # no gene lost
    assert {g.seqid for g in out} == {"chr1", "chr2"}


def test_parallel_finalize_merges_worker_stats(genome_fasta, helixer_gff3, tmp_path):
    serial_cfg = _config(genome_fasta, helixer_gff3, tmp_path, finalize_workers=1)
    par_cfg = _config(genome_fasta, helixer_gff3, tmp_path, finalize_workers=3)

    genes, inputs = _genes_and_inputs(serial_cfg)
    s_stats, p_stats = RunStats(), RunStats()
    try:
        pl._finalize_genes(genes, serial_cfg, inputs, stats=s_stats)
        pl._finalize_genes(genes, par_cfg, inputs, stats=p_stats,
                           _executor_cls=ThreadPoolExecutor)
    finally:
        inputs["genome"].close()
    # No miniprot/TransDecoder evidence → all three counted as "rescued none",
    # and the parallel path sums the per-worker counters back identically.
    assert p_stats.backstop_rescued_none == s_stats.backstop_rescued_none == 3


# ---------------------------------------------------------------------------
# Batched TransDecoder == per-gene
# ---------------------------------------------------------------------------

def _fake_transdecoder_factory():
    """A run_transdecoder mock: ORF = [0, (len//3)*3) for every FASTA record."""
    import os

    def fake(prepared_fasta, out_dir, **kw):
        records = {}
        name = None
        with open(prepared_fasta) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    records[name] = ""
                elif name is not None:
                    records[name] += line
        bed = os.path.join(out_dir, "out.transdecoder.bed")
        with open(bed, "w") as fh:
            for nm, seq in records.items():
                thick_end = (len(seq) // 3) * 3
                fh.write(
                    f"{nm}\t0\t{len(seq)}\t{nm}.p1\t0\t+\t0\t{thick_end}\t"
                    f"0\t1\t{len(seq)}\t0\n"
                )
        return bed

    return fake


def test_batched_transdecoder_equals_per_gene(genome_fasta, helixer_gff3, tmp_path,
                                              monkeypatch):
    cfg = _config(genome_fasta, helixer_gff3, tmp_path)
    genes, inputs = _genes_and_inputs(cfg)
    genome = inputs["genome"]
    g1 = next(g for g in genes if g.seqid == "chr1" and g.strand == "+")

    monkeypatch.setattr(cds_mod, "run_transdecoder", _fake_transdecoder_factory())
    try:
        # Per-gene path (the old behavior, transdecoder_bin_dir passed through).
        per_gene = cds_mod.assign_backstop_cds(
            g1, [], genome=genome, transdecoder_bin_dir="X"
        )
        # Batched path (one multi-FASTA invocation).
        batched = cds_mod.batch_backstop_transdecoder([g1], genome, "X")[0]
    finally:
        genome.close()

    assert per_gene.transcripts[0].cds is not None
    pg = [(s.start, s.end, s.phase) for s in per_gene.transcripts[0].cds]
    bt = [(s.start, s.end, s.phase) for s in batched.transcripts[0].cds]
    assert pg == bt                          # identical CDS from both paths
    assert batched.tier == per_gene.tier


def test_batched_transdecoder_skips_cdsbearing_and_mikado(genome_fasta, helixer_gff3,
                                                          tmp_path, monkeypatch):
    cfg = _config(genome_fasta, helixer_gff3, tmp_path)
    genes, inputs = _genes_and_inputs(cfg)
    genome = inputs["genome"]
    monkeypatch.setattr(cds_mod, "run_transdecoder", _fake_transdecoder_factory())
    try:
        stats = RunStats()
        out = cds_mod.batch_backstop_transdecoder(genes, genome, "X", stats=stats)
    finally:
        genome.close()
    # Every gene is a CDS-less backstop here → all three get a TransDecoder CDS,
    # and the rescue moves them from the "none" bucket into "transdecoder".
    assert all(g.transcripts[0].cds is not None for g in out)
    assert stats.backstop_rescued_transdecoder == 3
    assert stats.backstop_rescued_none == -3   # net of the move (no prior bumps)
