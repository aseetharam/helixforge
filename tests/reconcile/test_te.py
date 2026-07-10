"""Tests for reconcile/te.py — optional EDTA TE gating.

EDTA is the only TE signal. A good-ORF gene overlapping a configured TE *class*
above threshold is reclassified transposable_element; a gene overlapping a
knob/satellite is never reclassified (the Classification filter); below-threshold
overlap flags but does not gate; with no annotation behavior is unchanged. Both
strands, concrete literal coordinates.
"""

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.reconcile.te import (
    DEFAULT_TE_CLASSES,
    gate_te,
    parse_edta_te_intervals,
    te_overlap_fraction,
)


# --- EDTA GFF3 fixtures (1-based inclusive, like the real file) -------------

_TE_GFF3 = """\
##gff-version 3
chr1\tEDTA\tknob\t76\t255\t1445\t-\t.\tID=TE_homo_0;Name=knob180;Classification=knob/knob180
chr1\tEDTA\tLTR_retrotransposon\t1001\t2000\t999\t+\t.\tID=TE_0;Name=g;Classification=LTR/Gypsy
chr1\tEDTA\tCACTA_TIR_transposon\t3001\t4000\t999\t-\t.\tID=TE_1;Name=d;Classification=DNA/DTC
chr1\tEDTA\trepeat_region\t5001\t6000\t999\t+\t.\tID=R_0;Classification=Low_complexity
"""


def _write(tmp_path, text):
    p = tmp_path / "te.gff3"
    p.write_text(text)
    return str(p)


def _gene(strand, *, exons, cds=None, biotype="protein_coding", tier=2, seqid="chr1"):
    start = min(s for s, _ in exons)
    end = max(e for _, e in exons)
    tid = "HFG_00001.1"
    cds_segs = [CDSSegment(*c) for c in cds] if cds else None
    t = TranscriptCandidate(
        tid, "HFG_00001", "helixer", seqid, start, end, strand,
        [Exon(s, e) for s, e in exons], cds=cds_segs, is_primary=True,
    )
    return ReconciledGene(
        "HFG_00001", seqid, start, end, strand, tier, [t], tid,
        LocusClassification("HFG_00001", "SILENT"), "helixer_backstop",
        biotype=biotype,
    )


# --- parsing: only TE classes, coordinate conversion ------------------------

def test_parse_keeps_te_classes_drops_knob_and_lowcomplexity(tmp_path):
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3))
    ivals = index["chr1"].query_with_data(0, 10000)
    spans = sorted((s, e) for s, e, *_ in ivals)
    # LTR (1001-2000) + DNA/DTC (3001-4000) kept; knob + Low_complexity dropped.
    # 1-based inclusive -> 0-based half-open: 1001..2000 -> (1000, 2000).
    assert spans == [(1000, 2000), (3000, 4000)]


def test_parse_custom_classes(tmp_path):
    # Restrict to LTR only -> the DNA/DTC feature is dropped too.
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3), te_classes=["ltr"])
    spans = sorted((s, e) for s, e, *_ in index["chr1"].query_with_data(0, 10000))
    assert spans == [(1000, 2000)]


def test_default_classes_exclude_satellites():
    assert "knob" not in DEFAULT_TE_CLASSES
    assert "ltr" in DEFAULT_TE_CLASSES and "dna" in DEFAULT_TE_CLASSES


# --- overlap fraction (both strands) ----------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_overlap_fraction_full_and_partial(tmp_path, strand):
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3))
    idx = index["chr1"]
    # Fully inside the LTR interval (1000,2000): fraction 1.0.
    assert te_overlap_fraction([Exon(1200, 1400)], idx) == pytest.approx(1.0)
    # 100 of 400 exonic bases overlap the LTR -> 0.25.
    frac = te_overlap_fraction([Exon(1900, 2300)], idx)
    assert frac == pytest.approx(100 / 400)
    # No overlap -> 0.0.
    assert te_overlap_fraction([Exon(8000, 8500)], idx) == 0.0


# --- gating (both strands) --------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_good_orf_over_te_class_is_reclassified(tmp_path, strand):
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3))
    gene = _gene(strand, exons=[(1000, 1300)], cds=[(1000, 1300, 0)])  # in LTR
    out, flagged, reclass = gate_te([gene], index, threshold=0.5)
    assert flagged == 1 and reclass == 1
    g = out[0]
    assert g.biotype == "transposable_element"
    assert g.transcripts[0].biotype == "transposable_element"
    assert g.tier == 4
    assert any(f.name == "TE_OVERLAP" for f in g.flags)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_good_orf_over_knob_stays_coding(tmp_path, strand):
    # The knob feature is not a TE class -> no overlap -> stays protein_coding.
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3))
    gene = _gene(strand, exons=[(100, 250)], cds=[(100, 250, 0)])  # over the knob
    out, flagged, reclass = gate_te([gene], index, threshold=0.5)
    assert flagged == 0 and reclass == 0
    assert out[0].biotype == "protein_coding"
    assert not any(f.name == "TE_OVERLAP" for f in out[0].flags)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_below_threshold_flags_but_does_not_reclassify(tmp_path, strand):
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3))
    # 100/400 = 0.25 overlap, below the 0.5 threshold.
    gene = _gene(strand, exons=[(1900, 2300)], cds=[(1900, 2200, 0)])
    out, flagged, reclass = gate_te([gene], index, threshold=0.5)
    assert flagged == 1 and reclass == 0
    g = out[0]
    assert g.biotype == "protein_coding"  # gating distinct from flagging
    assert any(f.name == "TE_OVERLAP" for f in g.flags)


def test_noncoding_gene_flagged_not_reclassified(tmp_path):
    # Only protein_coding genes are gated; a non-coding overlap is flag-only.
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3))
    gene = _gene("+", exons=[(1000, 1300)], cds=None, biotype="lncRNA")
    out, flagged, reclass = gate_te([gene], index, threshold=0.5)
    assert flagged == 1 and reclass == 0
    assert out[0].biotype == "lncRNA"
    assert any(f.name == "TE_OVERLAP" for f in out[0].flags)


def test_no_te_on_other_seqid_is_noop(tmp_path):
    index = parse_edta_te_intervals(_write(tmp_path, _TE_GFF3))
    gene = _gene("+", exons=[(1000, 1300)], cds=[(1000, 1300, 0)], seqid="chrX")
    out, flagged, reclass = gate_te([gene], index, threshold=0.5)
    assert flagged == 0 and reclass == 0
    assert out[0] is gene  # untouched


# --- a slice of the real maize EDTA sample file -----------------------------

def test_parse_real_sample_slice():
    """The committed maize EDTA head is all knob/knob180 — none should survive."""
    import pathlib

    sample = pathlib.Path(
        "/home/arnstrm/svn/helixforge.v3/helixforge_testing/"
        "Zm-B73-REFERENCE-NAM-5.0.TE.gff3"
    )
    if not sample.exists():
        pytest.skip("maize EDTA sample not present")
    # Read just the header rows (all knob) into a temp file.
    lines = []
    with sample.open() as fh:
        for line in fh:
            lines.append(line)
            if len(lines) >= 40:
                break
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".gff3", delete=False) as tf:
        tf.writelines(lines)
        path = tf.name
    index = parse_edta_te_intervals(path)
    # The first 40 rows are knob/knob180 (a satellite) -> nothing kept.
    assert all(idx.is_empty for idx in index.values()) or not index
