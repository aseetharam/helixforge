"""Synthetic fixtures for the parallel (scatter-gather) tests (Phase 16).

Concrete literal coordinates only (CLAUDE.md §12); both strands represented. No
real biological data — a tiny FASTA + Helixer GFF3 with **hand-placed inter-locus
gaps** so the partition's cut points are predictable.
"""

import pytest

# Helixer-style genes as (gene_id, seqid, start1, end1, strand) — 1-based GFF3.
# chr1 carries 5 plus-strand loci with two large (>=1000 bp) gaps at the right
# places; chr2 is a single minus-strand locus (whole-scaffold chunk).
_GENES = [
    ("G1", "chr1", 1, 500, "+"),       # internal (0, 500)
    ("G2", "chr1", 601, 1000, "+"),    # internal (600, 1000); gap 100 to G1
    ("G3", "chr1", 6001, 6500, "+"),   # internal (6000, 6500); gap 5000 to G2
    ("G4", "chr1", 6601, 7000, "+"),   # internal (6600, 7000); gap 100 to G3
    ("G5", "chr1", 13001, 13500, "+"), # internal (13000, 13500); gap 6000 to G4
    ("G6", "chr2", 1001, 1500, "-"),   # internal (1000, 1500); minus strand
]
_SCAFFOLD_LEN = {"chr1": 20000, "chr2": 5000}


def _write_fasta(path, lengths):
    lines = []
    for sid, length in lengths.items():
        lines.append(f">{sid}")
        lines.append("A" * length)
    path.write_text("\n".join(lines) + "\n")


def _write_gff(path, genes):
    lines = ["##gff-version 3"]
    for gid, seqid, s, e, strand in genes:
        mrna = f"{gid}.t1"
        lines.append(f"{seqid}\tTest\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={gid}")
        lines.append(f"{seqid}\tTest\tmRNA\t{s}\t{e}\t.\t{strand}\t.\tID={mrna};Parent={gid}")
        lines.append(f"{seqid}\tTest\texon\t{s}\t{e}\t.\t{strand}\t.\tID={mrna}.exon1;Parent={mrna}")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture
def synthetic_genome(tmp_path):
    """Write the synthetic FASTA + Helixer GFF3; return their paths + metadata."""
    fasta = tmp_path / "genome.fasta"
    gff3 = tmp_path / "helixer.gff3"
    _write_fasta(fasta, _SCAFFOLD_LEN)
    _write_gff(gff3, _GENES)
    return {
        "fasta": str(fasta),
        "gff3": str(gff3),
        "scaffold_len": dict(_SCAFFOLD_LEN),
        "genes": list(_GENES),
        "total_loci": len(_GENES),
        "tmp_path": tmp_path,
    }
