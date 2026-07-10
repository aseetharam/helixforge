"""Phase 9 D3 — export writers (GTF / FASTA / AGAT-clean GFF3 / per-gene JSON).

Concrete literal coordinates only (CLAUDE.md §12); FASTA translation is checked
on **both strands**. Synthetic fixtures + a tiny ``MockGenome``.
"""

import json

import pytest

from helixforge.export.writers import (
    build_gene_records,
    write_agat_clean_gff3,
    write_cdna_fasta,
    write_cds_fasta,
    write_gtf,
    write_per_gene_json,
    write_protein_fasta,
)
from helixforge.io.gff import GFF3Parser
from helixforge.reconcile.models import (
    ASEvent,
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.utils.sequences import reverse_complement


# chrP: plus ORF, CDS [3,18) = ATG AAA CCC GGG TTT → protein MKPGF.
_CHRP = "CCC" + "ATGAAACCCGGGTTT" + "TAA" + "G" * 176
# chrM: minus ORF, reading '-' over CDS [3,18) yields ATG AAA CCC GGG TTT.
_CHRM = "TTA" + reverse_complement("ATGAAACCCGGGTTT") + "G" * 182
# chrE: two exons for cDNA splicing. [10,13)=AAA, [20,23)=GGG.
_CHRE = list("C" * 60)
_CHRE[10:13] = "AAA"
_CHRE[20:23] = "GGG"
_CHRE = "".join(_CHRE)


class MockGenome:
    def __init__(self, seqs):
        self.seqs = dict(seqs)

    def get_sequence(self, seqid, start, end, strand="+"):
        seq = self.seqs[seqid][start:end]
        return reverse_complement(seq) if strand == "-" else seq


@pytest.fixture
def genome():
    return MockGenome({"chrP": _CHRP, "chrM": _CHRM, "chrE": _CHRE})


def make_tx(tid, seqid="chr1", strand="+", exon_bounds=((1000, 1200),), cds=None,
            cds_partial=False, protein_id=None):
    return TranscriptCandidate(
        transcript_id=tid,
        locus_id=tid.rsplit(".", 1)[0],
        source="mikado",
        seqid=seqid,
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        exons=[Exon(s, e) for s, e in exon_bounds],
        cds=[CDSSegment(*c) for c in cds] if cds else None,
        cds_partial=cds_partial,
        tpm=10.0,
        junction_support_fraction=1.0,
        confidence=0.8,
        protein_id=protein_id,
        combined_score=15.0,
        is_primary=True,
    )


def make_gene(gid, txs, seqid="chr1", strand="+", origin="mikado_1to1", tier=1,
              as_events=None, flags=None):
    if not isinstance(txs, list):
        txs = [txs]
    return ReconciledGene(
        gene_id=gid,
        seqid=seqid,
        start=min(t.start for t in txs),
        end=max(t.end for t in txs),
        strand=strand,
        tier=tier,
        transcripts=txs,
        primary_transcript_id=txs[0].transcript_id,
        classification=LocusClassification(locus_id=gid, status="EXPRESSED", max_tpm=10.0),
        origin=origin,
        as_events=as_events or [],
        flags=flags or [],
    )


def _read_fasta(path):
    records, name, seq = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(seq)
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name is not None:
        records[name] = "".join(seq)
    return records


# ---------------------------------------------------------------------------
# protein / cds FASTA — both strands
# ---------------------------------------------------------------------------


def test_protein_fasta_plus(tmp_path, genome):
    gene = make_gene("HFG_00001",
                     make_tx("HFG_00001.1", seqid="chrP", exon_bounds=((0, 60),),
                             cds=((3, 18, 0),)), seqid="chrP")
    out = tmp_path / "prot.fa"
    write_protein_fasta([gene], genome, out)
    rec = _read_fasta(out)
    assert rec["HFG_00001.1"] == "MKPGF"


def test_protein_fasta_minus(tmp_path, genome):
    gene = make_gene("HFG_00002",
                     make_tx("HFG_00002.1", seqid="chrM", strand="-",
                             exon_bounds=((0, 60),), cds=((3, 18, 0),)),
                     seqid="chrM", strand="-")
    out = tmp_path / "prot.fa"
    write_protein_fasta([gene], genome, out)
    rec = _read_fasta(out)
    assert rec["HFG_00002.1"] == "MKPGF"


def test_cds_fasta_plus(tmp_path, genome):
    gene = make_gene("HFG_00001",
                     make_tx("HFG_00001.1", seqid="chrP", exon_bounds=((0, 60),),
                             cds=((3, 18, 0),)), seqid="chrP")
    out = tmp_path / "cds.fa"
    write_cds_fasta([gene], genome, out)
    rec = _read_fasta(out)
    assert rec["HFG_00001.1"] == "ATGAAACCCGGGTTT"


def test_cds_fasta_minus(tmp_path, genome):
    gene = make_gene("HFG_00002",
                     make_tx("HFG_00002.1", seqid="chrM", strand="-",
                             exon_bounds=((0, 60),), cds=((3, 18, 0),)),
                     seqid="chrM", strand="-")
    out = tmp_path / "cds.fa"
    write_cds_fasta([gene], genome, out)
    rec = _read_fasta(out)
    assert rec["HFG_00002.1"] == "ATGAAACCCGGGTTT"


def test_cdna_fasta_plus(tmp_path, genome):
    gene = make_gene("HFG_00003",
                     make_tx("HFG_00003.1", seqid="chrE", exon_bounds=((10, 13), (20, 23))),
                     seqid="chrE")
    out = tmp_path / "cdna.fa"
    write_cdna_fasta([gene], genome, out)
    rec = _read_fasta(out)
    assert rec["HFG_00003.1"] == "AAAGGG"


def test_cdna_fasta_minus(tmp_path, genome):
    gene = make_gene("HFG_00004",
                     make_tx("HFG_00004.1", seqid="chrE", strand="-",
                             exon_bounds=((10, 13), (20, 23))),
                     seqid="chrE", strand="-")
    out = tmp_path / "cdna.fa"
    write_cdna_fasta([gene], genome, out)
    rec = _read_fasta(out)
    # coding order reversed, each piece RC'd: RC(GGG)+RC(AAA) = CCC+TTT.
    assert rec["HFG_00004.1"] == "CCCTTT"


def test_protein_fasta_skips_cds_less(tmp_path, genome):
    gene = make_gene("HFG_00005", make_tx("HFG_00005.1", exon_bounds=((1000, 1300),)))
    out = tmp_path / "prot.fa"
    write_protein_fasta([gene], genome, out)
    assert _read_fasta(out) == {}


# ---------------------------------------------------------------------------
# GTF
# ---------------------------------------------------------------------------


def test_gtf_coordinates_are_1based(tmp_path):
    gene = make_gene("HFG_00001",
                     make_tx("HFG_00001.1", exon_bounds=((1000, 1200),),
                             cds=((1050, 1200, 0),)))
    out = tmp_path / "out.gtf"
    write_gtf([gene], out)
    lines = [ln.split("\t") for ln in out.read_text().splitlines()]
    exon = next(c for c in lines if c[2] == "exon")
    assert exon[3] == "1001"   # internal 1000 → GTF 1001
    assert exon[4] == "1200"
    assert 'transcript_id "HFG_00001.1"' in exon[8]


def test_gtf_cds_phase(tmp_path):
    gene = make_gene("HFG_00001",
                     make_tx("HFG_00001.1", exon_bounds=((1000, 1200),),
                             cds=((1050, 1199, 1),), cds_partial=True))
    out = tmp_path / "out.gtf"
    write_gtf([gene], out)
    lines = [ln.split("\t") for ln in out.read_text().splitlines()]
    cds = next(c for c in lines if c[2] == "CDS")
    assert cds[7] == "1"


def test_gtf_minus_strand(tmp_path):
    gene = make_gene("HFG_00002",
                     make_tx("HFG_00002.1", strand="-", exon_bounds=((2000, 2300),)),
                     strand="-")
    out = tmp_path / "out.gtf"
    write_gtf([gene], out)
    lines = [ln.split("\t") for ln in out.read_text().splitlines()]
    assert all(c[6] == "-" for c in lines)


def test_gtf_emits_biotype(tmp_path):
    # Phase 30 D2: GTF carries Ensembl gene_biotype / transcript_biotype.
    import attrs

    tx = make_tx("HFG_00003.1", exon_bounds=((1000, 1200), (1400, 1600)))
    tx = attrs.evolve(tx, biotype="lncRNA", cds=None)
    gene = attrs.evolve(make_gene("HFG_00003", tx, tier=3), biotype="lncRNA")
    out = tmp_path / "biotype.gtf"
    write_gtf([gene], out)
    text = out.read_text()
    assert 'gene_biotype "lncRNA"' in text
    assert 'transcript_biotype "lncRNA"' in text


# ---------------------------------------------------------------------------
# AGAT-clean GFF3
# ---------------------------------------------------------------------------


def test_agat_clean_gff3_reparses_with_biotype(tmp_path):
    # Phase 30 D2: AGAT-clean output still reparses with biotype attributes set
    # (a CDS-less lncRNA gene + a coding gene), both strands.
    import attrs

    lnc_tx = attrs.evolve(
        make_tx("HFG_00002.1", seqid="chr2", strand="-",
                exon_bounds=((500, 700), (800, 1000))),
        biotype="lncRNA", cds=None,
    )
    g1 = attrs.evolve(
        make_gene("HFG_00002", lnc_tx, seqid="chr2", strand="-", tier=3,
                  origin="helixer_backstop"),
        biotype="lncRNA",
    )
    coding_tx = attrs.evolve(
        make_tx("HFG_00001.1", exon_bounds=((1000, 1200),), cds=((1050, 1200, 0),)),
        biotype="protein_coding",
    )
    g2 = attrs.evolve(make_gene("HFG_00001", coding_tx), biotype="protein_coding")
    out = tmp_path / "agat.gff3"
    write_agat_clean_gff3([g1, g2], out)
    text = out.read_text()
    assert "gene_biotype=lncRNA" in text
    assert "gene_biotype=protein_coding" in text
    parsed = GFF3Parser(str(out)).parse_helixer_genes()
    assert {g.gene_id for g in parsed} == {"HFG_00001", "HFG_00002"}


def test_per_gene_json_has_biotype(tmp_path):
    import attrs

    tx = attrs.evolve(make_tx("HFG_00004.1", exon_bounds=((1000, 1300),)),
                      biotype="ncRNA_undetermined", cds=None)
    gene = attrs.evolve(make_gene("HFG_00004", tx, tier=3), biotype="ncRNA_undetermined")
    rec = build_gene_records([gene])[0]
    assert rec["biotype"] == "ncRNA_undetermined"
    assert rec["transcripts"][0]["biotype"] == "ncRNA_undetermined"


def test_agat_clean_gff3_reparses(tmp_path):
    g1 = make_gene("HFG_00002",
                   make_tx("HFG_00002.1", seqid="chr2", exon_bounds=((500, 800),),
                           cds=((500, 800, 0),)), seqid="chr2")
    g2 = make_gene("HFG_00001",
                   make_tx("HFG_00001.1", exon_bounds=((1000, 1200),),
                           cds=((1050, 1200, 0),)))
    out = tmp_path / "agat.gff3"
    write_agat_clean_gff3([g1, g2], out)
    parsed = GFF3Parser(str(out)).parse_helixer_genes()
    ids = [g.gene_id for g in parsed]
    assert set(ids) == {"HFG_00001", "HFG_00002"}


def test_agat_clean_gff3_sorted_by_position(tmp_path):
    g_late = make_gene("HFG_00009",
                       make_tx("HFG_00009.1", exon_bounds=((9000, 9200),)))
    g_early = make_gene("HFG_00001",
                        make_tx("HFG_00001.1", exon_bounds=((1000, 1200),)))
    out = tmp_path / "agat.gff3"
    write_agat_clean_gff3([g_late, g_early], out)
    gene_lines = [ln for ln in out.read_text().splitlines()
                  if "\tgene\t" in ln]
    assert "HFG_00001" in gene_lines[0]
    assert "HFG_00009" in gene_lines[1]


# ---------------------------------------------------------------------------
# per-gene JSON
# ---------------------------------------------------------------------------


def test_per_gene_json_roundtrip(tmp_path):
    gene = make_gene(
        "HFG_00001",
        [
            make_tx("HFG_00001.1", exon_bounds=((1000, 1200), (1300, 1500)),
                    cds=((1050, 1200, 0), (1300, 1450, 0)), protein_id="P1"),
            make_tx("HFG_00001.2", exon_bounds=((1000, 1500),)),
        ],
        as_events=[ASEvent("ES", "chr1", 1300, 1500, "+", support_read_count=12)],
    )
    out = tmp_path / "genes.json"
    write_per_gene_json([gene], out)
    loaded = json.load(open(out))
    assert loaded == build_gene_records([gene])


def test_per_gene_json_structure(tmp_path):
    gene = make_gene(
        "HFG_00001",
        make_tx("HFG_00001.1", exon_bounds=((1000, 1200), (1300, 1500)),
                cds=((1050, 1200, 0), (1300, 1450, 0))),
        as_events=[ASEvent("IR", "chr1", 1200, 1300, "+")],
    )
    out = tmp_path / "genes.json"
    write_per_gene_json([gene], out)
    rec = json.load(open(out))[0]
    assert rec["gene_id"] == "HFG_00001"
    assert rec["num_isoforms"] == 1
    assert rec["tier"] == 1
    assert rec["origin"] == "mikado_1to1"
    assert rec["transcripts"][0]["exons"][0] == {"start": 1000, "end": 1200}
    assert rec["transcripts"][0]["cds"][0] == {"start": 1050, "end": 1200, "phase": 0}
    assert rec["as_events"][0]["kind"] == "IR"


def test_per_gene_json_coords_are_internal(tmp_path):
    gene = make_gene("HFG_00001", make_tx("HFG_00001.1", exon_bounds=((1000, 1200),)))
    out = tmp_path / "genes.json"
    write_per_gene_json([gene], out)
    rec = json.load(open(out))[0]
    # internal 0-based half-open preserved (no +1 conversion in JSON).
    assert rec["start"] == 1000
    assert rec["end"] == 1200
