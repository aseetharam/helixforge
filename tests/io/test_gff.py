"""Tests for GFF3Parser + GFF3Writer (Phase 1). Floor: 23.

Round-trip coordinate-fidelity tests are critical (write ReconciledGene →
parse → coordinates match exactly).
"""

import pytest

from helixforge.io.gff import GFF3Parser, GFF3Writer
from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def test_parse_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        GFF3Parser("/no/such/file.gff3")


def test_parse_helixer_gene_count(helixer_gff_path):
    loci = GFF3Parser(helixer_gff_path).parse_helixer_genes()
    assert len(loci) == 3


def test_parse_helixer_sorted_by_seqid_start(helixer_gff_path):
    loci = GFF3Parser(helixer_gff_path).parse_helixer_genes()
    assert [g.gene_id for g in loci] == ["gene1", "gene2", "gene3"]


def test_parse_helixer_gene1_coords_converted(helixer_gff_path):
    # GFF 11-40 (1-based inclusive) -> internal (10, 40)
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[0]
    assert g.start == 10
    assert g.end == 40


def test_parse_helixer_gene1_strand_plus(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[0]
    assert g.strand == "+"


def test_parse_helixer_gene1_exons(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[0]
    assert [(e.start, e.end) for e in g.exons] == [(10, 20), (30, 40)]


def test_parse_helixer_gene1_cds(helixer_gff_path):
    # CDS 12-20 -> (11,20); CDS 31-39 -> (30,39)
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[0]
    assert [(c.start, c.end) for c in g.cds] == [(11, 20), (30, 39)]


def test_parse_helixer_cds_phase(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[0]
    assert all(c.phase == 0 for c in g.cds)


def test_parse_helixer_gene2_minus_strand(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[1]
    assert g.gene_id == "gene2"
    assert g.strand == "-"


def test_parse_helixer_gene2_exons_minus(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[1]
    assert [(e.start, e.end) for e in g.exons] == [(10, 25), (30, 40)]


def test_parse_helixer_gene2_cds_minus(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[1]
    assert [(c.start, c.end) for c in g.cds] == [(10, 25), (30, 39)]


def test_parse_helixer_single_exon_gene(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[2]
    assert g.gene_id == "gene3"
    assert [(e.start, e.end) for e in g.exons] == [(49, 58)]


def test_parse_helixer_cds_optional_is_none(helixer_gff_path):
    g = GFF3Parser(helixer_gff_path).parse_helixer_genes()[2]
    assert g.cds is None


def test_parse_genes_generic_returns_dicts(helixer_gff_path):
    genes = GFF3Parser(helixer_gff_path).parse_genes_generic()
    assert len(genes) == 3
    g1 = genes[0]
    assert g1["gene_id"] == "gene1"
    assert g1["start"] == 10
    assert g1["end"] == 40
    assert g1["strand"] == "+"


def test_parse_genes_generic_transcript_exons(helixer_gff_path):
    genes = GFF3Parser(helixer_gff_path).parse_genes_generic()
    tx = genes[0]["transcripts"][0]
    assert [(e.start, e.end) for e in tx["exons"]] == [(10, 20), (30, 40)]


def test_get_features_in_region_overlap(helixer_gff_path):
    parser = GFF3Parser(helixer_gff_path)
    feats = parser.get_features_in_region("chr1", 10, 20, featuretype="gene")
    assert [f.id for f in feats] == ["gene1"]


def test_get_features_in_region_no_overlap(helixer_gff_path):
    parser = GFF3Parser(helixer_gff_path)
    feats = parser.get_features_in_region("chr1", 45, 60, featuretype="gene")
    assert feats == []


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------

def test_write_produces_gff_version_header(tmp_path, sample_reconciled_gene):
    out = tmp_path / "out.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    assert out.read_text().splitlines()[0] == "##gff-version 3"


def test_write_uses_gene_separator(tmp_path, sample_reconciled_gene):
    out = tmp_path / "out.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene, sample_reconciled_gene])
    assert out.read_text().count("###") == 2


def test_write_exon_and_cds_ids(tmp_path, sample_reconciled_gene):
    out = tmp_path / "out.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    text = out.read_text()
    assert "ID=HFG_00001.1.exon1" in text
    assert "ID=HFG_00001.1.CDS1" in text


def test_write_mrna_attributes_present(tmp_path, sample_reconciled_gene):
    out = tmp_path / "out.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    text = out.read_text()
    assert "structure_source=mikado" in text
    assert "tier=1" in text
    assert "origin=mikado_1to1" in text


def test_write_gene_and_transcript_biotype(tmp_path, sample_reconciled_gene):
    # Phase 30 D2: gene_biotype / transcript_biotype emitted (Ensembl convention).
    import attrs

    gene = attrs.evolve(
        sample_reconciled_gene,
        biotype="protein_coding",
        transcripts=[attrs.evolve(t, biotype="protein_coding")
                     for t in sample_reconciled_gene.transcripts],
    )
    out = tmp_path / "biotype.gff3"
    GFF3Writer(out).write_genes([gene])
    text = out.read_text()
    assert "gene_biotype=protein_coding" in text
    assert "transcript_biotype=protein_coding" in text


def test_write_lncrna_biotype_both_strands(tmp_path):
    # A CDS-less lncRNA gene on each strand carries gene/transcript biotype=lncRNA.
    import attrs

    for strand in ("+", "-"):
        tx = TranscriptCandidate(
            "HFG_00002.1", "HFG_00002", "stringtie", "chr1", 1000, 1600, strand,
            [Exon(1000, 1200), Exon(1400, 1600)], is_primary=True, biotype="lncRNA",
        )
        gene = ReconciledGene(
            "HFG_00002", "chr1", 1000, 1600, strand, 3, [tx], "HFG_00002.1",
            LocusClassification("HFG_00002", "EXPRESSED"), "helixer_backstop",
            biotype="lncRNA",
        )
        out = tmp_path / f"lnc_{strand}.gff3"
        GFF3Writer(out).write_genes([gene])
        text = out.read_text()
        assert "gene_biotype=lncRNA" in text
        assert "transcript_biotype=lncRNA" in text


def test_write_no_biotype_when_unset(tmp_path, sample_reconciled_gene):
    # biotype=None (not yet classified) -> attribute omitted (no empty value).
    out = tmp_path / "nobiotype.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    text = out.read_text()
    assert "gene_biotype=" not in text
    assert "transcript_biotype=" not in text


def test_write_url_encodes_special_chars(tmp_path):
    # protein_id with '=' must be escaped to %3D in the attribute value.
    tx = TranscriptCandidate(
        "HFG_00009.1", "HFG_00009", "miniprot", "chr1", 1000, 1200, "+",
        [Exon(1000, 1200)], protein_id="sp|P1=2",
    )
    cls = LocusClassification("HFG_00009", "LOW")
    gene = ReconciledGene(
        "HFG_00009", "chr1", 1000, 1200, "+", 3, [tx], "HFG_00009.1", cls, "novel"
    )
    out = tmp_path / "enc.gff3"
    GFF3Writer(out).write_genes([gene])
    text = out.read_text()
    assert "sp|P1%3D2" in text
    assert "sp|P1=2" not in text


# --------------------------------------------------------------------------
# Round-trip coordinate fidelity (CRITICAL)
# --------------------------------------------------------------------------

def test_roundtrip_coordinate_fidelity_plus(tmp_path, sample_reconciled_gene):
    out = tmp_path / "rt_plus.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    loci = GFF3Parser(str(out)).parse_helixer_genes()
    assert len(loci) == 1
    g = loci[0]
    assert (g.start, g.end) == (1000, 1800)
    assert [(e.start, e.end) for e in g.exons] == [(1000, 1200), (1300, 1500), (1600, 1800)]
    assert [(c.start, c.end) for c in g.cds] == [(1050, 1200), (1300, 1500), (1600, 1649)]


def test_roundtrip_preserves_strand_plus(tmp_path, sample_reconciled_gene):
    out = tmp_path / "rt_plus2.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    g = GFF3Parser(str(out)).parse_helixer_genes()[0]
    assert g.strand == "+"
    assert g.gene_id == "HFG_00001"


def test_roundtrip_preserves_cds_phase(tmp_path, sample_reconciled_gene):
    # sample_cds_segments phases are [0, 0, 1]
    out = tmp_path / "rt_phase.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    g = GFF3Parser(str(out)).parse_helixer_genes()[0]
    assert [c.phase for c in g.cds] == [0, 0, 1]


def test_roundtrip_coordinate_fidelity_minus(tmp_path):
    exons = [Exon(2000, 2200), Exon(2300, 2500), Exon(2600, 2800)]
    cds = [CDSSegment(2050, 2200, 0), CDSSegment(2300, 2500, 0), CDSSegment(2600, 2649, 1)]
    tx = TranscriptCandidate(
        "HFG_00002.1", "HFG_00002", "mikado", "chr1", 2000, 2800, "-", exons, cds=cds
    )
    cls = LocusClassification("HFG_00002", "EXPRESSED")
    gene = ReconciledGene(
        "HFG_00002", "chr1", 2000, 2800, "-", 1, [tx], "HFG_00002.1", cls, "mikado_1to1"
    )
    out = tmp_path / "rt_minus.gff3"
    GFF3Writer(out).write_genes([gene])
    g = GFF3Parser(str(out)).parse_helixer_genes()[0]
    assert g.strand == "-"
    assert (g.start, g.end) == (2000, 2800)
    assert [(e.start, e.end) for e in g.exons] == [(2000, 2200), (2300, 2500), (2600, 2800)]
    assert [(c.start, c.end) for c in g.cds] == [(2050, 2200), (2300, 2500), (2600, 2649)]


def test_roundtrip_single_exon_no_cds(tmp_path):
    tx = TranscriptCandidate(
        "HFG_00003.1", "HFG_00003", "helixer", "chr1", 5000, 5300, "+", [Exon(5000, 5300)]
    )
    cls = LocusClassification("HFG_00003", "SILENT")
    gene = ReconciledGene(
        "HFG_00003", "chr1", 5000, 5300, "+", 4, [tx], "HFG_00003.1", cls, "helixer_backstop"
    )
    out = tmp_path / "rt_single.gff3"
    GFF3Writer(out).write_genes([gene])
    g = GFF3Parser(str(out)).parse_helixer_genes()[0]
    assert [(e.start, e.end) for e in g.exons] == [(5000, 5300)]
    assert g.cds is None


def test_roundtrip_output_reparses_cleanly(tmp_path, sample_reconciled_gene):
    out = tmp_path / "rt_clean.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    # generic parser must also read it without error
    genes = GFF3Parser(str(out)).parse_genes_generic()
    assert genes[0]["gene_id"] == "HFG_00001"


# --------------------------------------------------------------------------
# Mikado-finalizer safety: a written multi-exon coding transcript must emit
# per-exon CDS segments whose derived introns coincide with the exon introns.
# Mikado's transcript finalizer asserts ``len(cds_introns) > 0``; a single CDS
# spanning multiple exons, or an internal CDS boundary not on a splice site,
# yields zero matching CDS introns and crashes it. Both strands (CLAUDE.md §12).
# --------------------------------------------------------------------------

def _mikado_cds_introns_match(g):
    """Replicate Mikado's coherence rule: every CDS-derived intron is an exon
    intron. Returns (num_cds_features, all_match) for a reparsed locus."""
    cds = sorted((c.start, c.end) for c in (g.cds or []))
    ex = sorted((e.start, e.end) for e in g.exons)
    exon_introns = {(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)}
    cds_introns = [(cds[i][1], cds[i + 1][0]) for i in range(len(cds) - 1)]
    return len(cds), all(ci in exon_introns for ci in cds_introns)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_written_coding_cds_is_finalizer_safe(tmp_path, strand):
    # 3-exon coding transcript; CDS split per exon at the splice sites, flush at
    # every internal boundary (5' UTR in exon1, 3' UTR in exon3). 150+200+49=399.
    exons = [Exon(1000, 1200), Exon(1300, 1500), Exon(1600, 1800)]
    cds = [CDSSegment(1050, 1200, 0), CDSSegment(1300, 1500, 0), CDSSegment(1600, 1649, 1)]
    tx = TranscriptCandidate(
        "HFG_00009.1", "HFG_00009", "helixer", "chr1", 1000, 1800, strand, exons, cds=cds
    )
    cls = LocusClassification("HFG_00009", "EXPRESSED")
    gene = ReconciledGene(
        "HFG_00009", "chr1", 1000, 1800, strand, 2, [tx], "HFG_00009.1", cls, "mikado_1to1"
    )
    out = tmp_path / f"safe_{strand}.gff3"
    GFF3Writer(out).write_genes([gene])
    g = GFF3Parser(str(out)).parse_helixer_genes()[0]
    n_cds, all_match = _mikado_cds_introns_match(g)
    assert n_cds == 3  # per-exon CDS segments, not one spanning feature
    assert all_match  # every CDS intron coincides with an exon intron
    assert [c.phase for c in g.cds] == [0, 0, 1]  # phase preserved through I/O


# --------------------------------------------------------------------------
# Phase 19 D3: gffutils DB reuse via dbfn / keep_db (build once, reuse on disk).
# Default dbfn=None preserves the historical in-memory behavior exactly.
# --------------------------------------------------------------------------

def _loci_tuples(loci):
    return [
        (
            g.gene_id, g.seqid, g.start, g.end, g.strand,
            [(e.start, e.end) for e in g.exons],
            [(c.start, c.end, c.phase) for c in (g.cds or [])],
        )
        for g in loci
    ]


def test_persistent_dbfn_matches_memory(helixer_gff_path, tmp_path):
    in_memory = GFF3Parser(helixer_gff_path).parse_helixer_genes()
    dbfn = tmp_path / "helixer.gffutils.db"
    # first build: file absent + keep_db=True -> create_db(force=False) on disk
    on_disk = GFF3Parser(helixer_gff_path, dbfn=str(dbfn), keep_db=True).parse_helixer_genes()
    assert dbfn.exists()
    assert _loci_tuples(on_disk) == _loci_tuples(in_memory)


def test_dbfn_reuse_path_returns_same_loci(helixer_gff_path, tmp_path):
    dbfn = tmp_path / "helixer.gffutils.db"
    first = GFF3Parser(helixer_gff_path, dbfn=str(dbfn), keep_db=True).parse_helixer_genes()
    mtime = dbfn.stat().st_mtime_ns
    # second construction reuses the on-disk DB (FeatureDB open; no rebuild)
    reused = GFF3Parser(helixer_gff_path, dbfn=str(dbfn), keep_db=True).parse_helixer_genes()
    assert dbfn.stat().st_mtime_ns == mtime  # not rebuilt
    assert _loci_tuples(reused) == _loci_tuples(first)


# --------------------------------------------------------------------------
# Phase 31 D1: functional attributes (Ontology_term / Dbxref) round-trip
# --------------------------------------------------------------------------

class _FuncRec:
    """Minimal stand-in for prep.function.FunctionalRecord (go_terms/dbxrefs)."""

    def __init__(self, go_terms=(), dbxrefs=()):
        self.go_terms = tuple(go_terms)
        self.dbxrefs = tuple(dbxrefs)


def test_write_emits_ontology_term_and_dbxref(tmp_path, sample_reconciled_gene):
    out = tmp_path / "func.gff3"
    functional = {
        "HFG_00001.1": _FuncRec(
            go_terms=("GO:0004672", "GO:0005524"),
            dbxrefs=("InterPro:IPR000719", "Pfam:PF00069"),
        )
    }
    GFF3Writer(out).write_genes([sample_reconciled_gene], functional=functional)
    text = out.read_text()
    assert "Ontology_term=GO:0004672,GO:0005524" in text
    assert "Dbxref=InterPro:IPR000719,Pfam:PF00069" in text
    # AGAT-clean: the generic parser still reparses the hierarchy without error.
    genes = GFF3Parser(str(out)).parse_genes_generic()
    assert genes[0]["gene_id"] == "HFG_00001"


def test_write_omits_functional_attrs_when_absent(tmp_path, sample_reconciled_gene):
    # Default functional=None keeps the output free of the reserved attributes.
    out = tmp_path / "nofunc.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    text = out.read_text()
    assert "Ontology_term=" not in text
    assert "Dbxref=" not in text


# --------------------------------------------------------------------------
# Phase 32 D1: GFF3 spec completeness — ##sequence-region + provenance + ##FASTA
# --------------------------------------------------------------------------

def test_write_sequence_region_per_contig(tmp_path, sample_reconciled_gene):
    # Each contig gets a 1-based-inclusive ##sequence-region with its real length.
    out = tmp_path / "seqreg.gff3"
    GFF3Writer(out, sequence_regions={"chr1": 30427671, "chr2": 19698289}).write_genes(
        [sample_reconciled_gene]
    )
    lines = out.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert "##sequence-region chr1 1 30427671" in lines
    assert "##sequence-region chr2 1 19698289" in lines
    # sorted, deterministic order
    assert lines.index("##sequence-region chr1 1 30427671") < lines.index(
        "##sequence-region chr2 1 19698289"
    )


def test_write_sequence_region_minus_strand_gene(tmp_path):
    # A minus-strand gene still gets the contig directive; coordinates unchanged.
    tx = TranscriptCandidate(
        "HFG_00007.1", "HFG_00007", "mikado", "chr3", 2000, 2800, "-",
        [Exon(2000, 2800)], cds=[CDSSegment(2050, 2650, 0)],
    )
    gene = ReconciledGene(
        "HFG_00007", "chr3", 2000, 2800, "-", 1, [tx], "HFG_00007.1",
        LocusClassification("HFG_00007", "EXPRESSED"), "mikado_1to1",
    )
    out = tmp_path / "seqreg_minus.gff3"
    GFF3Writer(out, sequence_regions={"chr3": 23459830}).write_genes([gene])
    text = out.read_text()
    assert "##sequence-region chr3 1 23459830" in text
    g = GFF3Parser(str(out)).parse_helixer_genes()[0]
    assert g.strand == "-"
    assert (g.start, g.end) == (2000, 2800)


def test_write_provenance_preamble_present_and_parseable(tmp_path, sample_reconciled_gene):
    from helixforge.provenance import build_provenance

    prov = build_provenance(
        params={"scoring_profile": "strict", "pad": True},
        tool_bins={"mikado": "definitely-not-a-real-binary-xyz"},
        input_files={},
    )
    out = tmp_path / "prov.gff3"
    GFF3Writer(out, provenance=prov).write_genes([sample_reconciled_gene])
    lines = out.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert any(l.startswith("#!helixforge-version ") for l in lines)
    assert any(l.startswith("#!param-hash ") for l in lines)
    # an unresolvable tool is recorded as 'unresolved', never a crash
    assert "#!tool mikado=unresolved" in lines
    # the feature hierarchy still parses (comments ignored by gffutils)
    g = GFF3Parser(str(out)).parse_helixer_genes()[0]
    assert g.gene_id == "HFG_00001"


def test_write_default_header_byte_identical(tmp_path, sample_reconciled_gene):
    # No metadata args -> exactly the historical single-line header (count-neutral).
    out = tmp_path / "plain.gff3"
    GFF3Writer(out).write_genes([sample_reconciled_gene])
    lines = out.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert not any(l.startswith("##sequence-region") for l in lines)
    assert not any(l.startswith("#!") for l in lines)
    assert "##FASTA" not in out.read_text()


def test_write_embed_fasta_appends_sequence(tmp_path, sample_reconciled_gene):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGTAC\n")
    out = tmp_path / "embed.gff3"
    GFF3Writer(out, embed_fasta=str(fasta)).write_genes([sample_reconciled_gene])
    text = out.read_text()
    assert "##FASTA" in text
    assert ">chr1" in text
    assert "ACGTACGTAC" in text
    # ##FASTA must come after the feature block
    assert text.index("##FASTA") > text.index("HFG_00001")


def test_write_biotype_retained_with_sequence_region(tmp_path):
    # D1 must not disturb the Phase-30 biotype attributes.
    import attrs

    tx = TranscriptCandidate(
        "HFG_00008.1", "HFG_00008", "stringtie", "chr1", 1000, 1600, "+",
        [Exon(1000, 1200), Exon(1400, 1600)], is_primary=True, biotype="lncRNA",
    )
    gene = ReconciledGene(
        "HFG_00008", "chr1", 1000, 1600, "+", 3, [tx], "HFG_00008.1",
        LocusClassification("HFG_00008", "EXPRESSED"), "helixer_backstop",
        biotype="lncRNA",
    )
    out = tmp_path / "biotype_seqreg.gff3"
    GFF3Writer(out, sequence_regions={"chr1": 30427671}).write_genes([gene])
    text = out.read_text()
    assert "gene_biotype=lncRNA" in text
    assert "transcript_biotype=lncRNA" in text
    assert "##sequence-region chr1 1 30427671" in text
