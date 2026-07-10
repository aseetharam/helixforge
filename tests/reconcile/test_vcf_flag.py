"""Tests for the experimental VCF / haplotype scaffolding (Phase 32 D4, §1.8).

Concrete literal coordinates; both strands. The flag path is EXPERIMENTAL and
off by default — these tests pin (a) decomposed-only parsing, (b) high-impact CDS
overlap → VARIANT_IMPACTED, (c) the disabled / no-overlap identity.
"""

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.reconcile.vcf import (
    Variant,
    flag_variant_impacted_genes,
    load_vcf,
    parse_vcf_line,
    variant_is_high_impact,
)

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _gene(strand, *, gene_id="HFG_00001", seqid="chr1", cds_start=1050, cds_end=1200):
    tx = TranscriptCandidate(
        f"{gene_id}.1", gene_id, "mikado", seqid, 1000, 1500, strand,
        [Exon(1000, 1500)], cds=[CDSSegment(cds_start, cds_end, 0)], is_primary=True,
    )
    return ReconciledGene(
        gene_id, seqid, 1000, 1500, strand, 1, [tx], f"{gene_id}.1",
        LocusClassification(gene_id, "EXPRESSED"), "mikado_1to1",
    )


# --------------------------------------------------------------------------
# Decomposition requirement
# --------------------------------------------------------------------------

def test_multiallelic_rejected_with_clear_message(tmp_path):
    vcf = tmp_path / "multi.vcf"
    vcf.write_text(_VCF_HEADER + "chr1\t1100\t.\tA\tG,T\t.\t.\tIMPACT=HIGH\n")
    with pytest.raises(ValueError, match="bcftools norm -m -"):
        load_vcf(str(vcf))


def test_decomposed_biallelic_parses(tmp_path):
    vcf = tmp_path / "ok.vcf"
    vcf.write_text(
        _VCF_HEADER
        + "chr1\t1100\t.\tA\tG\t.\t.\tIMPACT=HIGH\n"
        + "chr1\t1300\t.\tC\tT\t.\t.\tANN=T|missense_variant|MODERATE|GENE\n"
    )
    variants = load_vcf(str(vcf))
    assert len(variants) == 2
    # 1-based POS 1100 -> internal start 1099, SNV end 1100
    assert (variants[0].seqid, variants[0].start, variants[0].end) == ("chr1", 1099, 1100)
    assert variants[0].high_impact is True
    assert variants[1].high_impact is False


def test_parse_line_skips_headers_and_blanks():
    assert parse_vcf_line("##fileformat=VCFv4.2") is None
    assert parse_vcf_line("#CHROM\tPOS") is None
    assert parse_vcf_line("") is None


def test_indel_ref_span(tmp_path):
    # REF 'ACG' at POS 1100 spans internal [1099, 1102).
    v = parse_vcf_line("chr1\t1100\t.\tACG\tA\t.\t.\tIMPACT=HIGH")
    assert (v.start, v.end) == (1099, 1102)


# --------------------------------------------------------------------------
# High-impact classification
# --------------------------------------------------------------------------

def test_high_impact_from_explicit_impact():
    assert variant_is_high_impact({"IMPACT": "HIGH"}) is True
    assert variant_is_high_impact({"IMPACT": "MODERATE"}) is False


def test_high_impact_from_ann_consequence():
    assert variant_is_high_impact({"ANN": "G|stop_gained|HIGH|GENE1|..."}) is True
    assert variant_is_high_impact({"CSQ": "T|splice_donor_variant|..."}) is True
    assert variant_is_high_impact({"ANN": "G|synonymous_variant|LOW|GENE1"}) is False


def test_no_annotation_is_not_high_impact():
    assert variant_is_high_impact({}) is False


# --------------------------------------------------------------------------
# Flagging genes (both strands)
# --------------------------------------------------------------------------

def test_high_impact_over_cds_flags_gene_plus():
    gene = _gene("+")  # CDS [1050, 1200)
    var = Variant("chr1", 1100, 1101, "A", "G", high_impact=True)
    out = flag_variant_impacted_genes([gene], [var], enabled=True)
    assert any(f.name == "VARIANT_IMPACTED" for f in out[0].flags)
    # structure preserved — flag only, no re-typing
    assert out[0].tier == 1
    assert out[0].biotype is None
    assert out[0].transcripts[0].cds == gene.transcripts[0].cds


def test_high_impact_over_cds_flags_gene_minus():
    gene = _gene("-")
    var = Variant("chr1", 1100, 1101, "A", "G", high_impact=True)
    out = flag_variant_impacted_genes([gene], [var], enabled=True)
    assert any(f.name == "VARIANT_IMPACTED" for f in out[0].flags)


def test_variant_outside_cds_not_flagged():
    gene = _gene("+")  # CDS [1050, 1200)
    # variant in the 3' UTR at 1300 -> no CDS overlap
    var = Variant("chr1", 1300, 1301, "A", "G", high_impact=True)
    out = flag_variant_impacted_genes([gene], [var], enabled=True)
    assert not any(f.name == "VARIANT_IMPACTED" for f in out[0].flags)


def test_moderate_impact_variant_not_flagged():
    gene = _gene("+")
    var = Variant("chr1", 1100, 1101, "A", "G", high_impact=False)
    out = flag_variant_impacted_genes([gene], [var], enabled=True)
    assert not any(f.name == "VARIANT_IMPACTED" for f in out[0].flags)


def test_disabled_is_identity():
    gene = _gene("+")
    var = Variant("chr1", 1100, 1101, "A", "G", high_impact=True)
    out = flag_variant_impacted_genes([gene], [var], enabled=False)
    assert out[0] is gene  # untouched object — pure identity when off
