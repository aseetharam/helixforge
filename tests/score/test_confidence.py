"""Tests for the standalone HDF5-only ``confidence`` scorer (Phase 14).

Metrics are checked to **match** ``mikado.emit_external`` on the same structures
(the scorer reuses that code, so it must agree bit for bit), every value is
asserted inside [0, 1], and both strands plus single/multi-exon are covered.
Coordinates in fixtures are literal (CLAUDE.md §12).
"""

import pytest

from helixforge.io.hdf5 import HDF5ConfidenceReader
from helixforge.mikado.emit_external import (
    helixer_locus_conf as ee_locus_conf,
    helixer_support as ee_support,
)
from helixforge.reconcile.models import Exon, TranscriptCandidate
from helixforge.score.confidence import (
    TSV_COLUMNS,
    score_annotation,
    score_transcript_confidence,
    summarize_confidence,
    write_confidence_tsv,
)

# Literal internal (0-based half-open) structures matching the GFF3 fixture.
PLUS_EXONS = [Exon(10, 20), Exon(40, 50)]      # chr1 +, intron [20,40)
MINUS_EXONS = [Exon(20, 40), Exon(60, 64)]     # chr2 -, intron [40,60)
SINGLE_EXON = [Exon(10, 20)]                    # chr1 +


@pytest.fixture
def reader(conf_hdf5_path):
    r = HDF5ConfidenceReader(conf_hdf5_path)
    yield r
    r.close()


def _tx(seqid, strand, exons):
    return {
        "seqid": seqid, "strand": strand, "exons": exons,
        "start": min(e.start for e in exons), "end": max(e.end for e in exons),
    }


# --- equality with emit_external ------------------------------------------

def test_support_matches_emit_external_plus(reader):
    got = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader)
    expected = ee_support(PLUS_EXONS, "chr1", "+", reader)
    assert got["helixer_support"] == expected


def test_support_matches_emit_external_minus(reader):
    got = score_transcript_confidence(_tx("chr2", "-", MINUS_EXONS), reader)
    expected = ee_support(MINUS_EXONS, "chr2", "-", reader)
    assert got["helixer_support"] == expected


def test_locus_conf_matches_emit_external_plus(reader):
    got = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader)
    expected = ee_locus_conf("chr1", 10, 50, reader)
    assert got["helixer_locus_conf"] == expected


def test_locus_conf_matches_emit_external_minus(reader):
    got = score_transcript_confidence(_tx("chr2", "-", MINUS_EXONS), reader)
    expected = ee_locus_conf("chr2", 20, 64, reader)
    assert got["helixer_locus_conf"] == expected


# --- bounds ---------------------------------------------------------------

@pytest.mark.parametrize("seqid,strand,exons", [
    ("chr1", "+", PLUS_EXONS),
    ("chr2", "-", MINUS_EXONS),
    ("chr1", "+", SINGLE_EXON),
])
def test_all_values_in_unit_interval(reader, seqid, strand, exons):
    m = score_transcript_confidence(_tx(seqid, strand, exons), reader)
    for key in ("helixer_support", "helixer_locus_conf",
                "exon_confidence", "intron_coincidence_fraction"):
        assert 0.0 <= m[key] <= 1.0, key


# --- structure counts -----------------------------------------------------

def test_multi_exon_counts(reader):
    m = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader)
    assert m["num_exons"] == 2
    assert m["num_introns"] == 1


def test_single_exon_has_no_introns(reader):
    m = score_transcript_confidence(_tx("chr1", "+", SINGLE_EXON), reader)
    assert m["num_exons"] == 1
    assert m["num_introns"] == 0
    assert m["intron_coincidence_fraction"] == 0.0


def test_single_exon_support_equals_exon_confidence(reader):
    m = score_transcript_confidence(_tx("chr1", "+", SINGLE_EXON), reader)
    # With no introns helixer_support collapses to the exon component.
    assert m["helixer_support"] == m["exon_confidence"]


def test_single_exon_lands_on_cds_channel(reader):
    # exon [10,20) sits entirely on the literal CDS channel (max(CDS,UTR)=0.90).
    m = score_transcript_confidence(_tx("chr1", "+", SINGLE_EXON), reader)
    assert m["exon_confidence"] == pytest.approx(0.90, abs=1e-4)


def test_supported_intron_coincidence_is_one(reader):
    # The plus intron [20,40) sits on the high intron channel (0.85 >= 0.5).
    m = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader)
    assert m["intron_coincidence_fraction"] == 1.0


# --- strand-agnostic value ------------------------------------------------

def test_support_is_strand_agnostic(reader):
    plus = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader)
    # Same structure scored as if minus must yield the identical value.
    minus = score_transcript_confidence(_tx("chr1", "-", PLUS_EXONS), reader)
    assert plus["helixer_support"] == minus["helixer_support"]


# --- adapter accepts both inputs ------------------------------------------

def test_transcript_candidate_accepted(reader):
    cand = TranscriptCandidate(
        transcript_id="t1", locus_id="L1", source="helixer",
        seqid="chr1", start=10, end=50, strand="+", exons=PLUS_EXONS,
    )
    from_obj = score_transcript_confidence(cand, reader)
    from_dict = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader)
    assert from_obj == from_dict


def test_generic_dict_accepted(reader):
    from helixforge.io.gff import GFF3Parser
    # parse_genes_generic transcript dicts carry exons but not seqid/strand;
    # the scorer's adapter takes those from the enriched dict.
    m = score_transcript_confidence(_tx("chr2", "-", MINUS_EXONS), reader)
    assert m["num_exons"] == 2


def test_no_exons_raises(reader):
    with pytest.raises(ValueError):
        score_transcript_confidence({"seqid": "chr1", "strand": "+", "exons": []}, reader)


# --- exon_weight / intron_cutoff parameters -------------------------------

def test_exon_weight_zero_uses_only_intron_fraction(reader):
    m = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader, exon_weight=0.0)
    assert m["helixer_support"] == pytest.approx(m["intron_coincidence_fraction"])


def test_exon_weight_one_uses_only_exon_confidence(reader):
    m = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader, exon_weight=1.0)
    assert m["helixer_support"] == pytest.approx(m["exon_confidence"])


def test_high_intron_cutoff_drops_coincidence(reader):
    # Intron channel peaks at 0.85; a cutoff above that means no intron counts.
    m = score_transcript_confidence(_tx("chr1", "+", PLUS_EXONS), reader, intron_high_cutoff=0.99)
    assert m["intron_coincidence_fraction"] == 0.0


# --- score_annotation -----------------------------------------------------

def test_score_annotation_all_rows(conf_gff3_path, conf_hdf5_path):
    df = score_annotation(conf_gff3_path, conf_hdf5_path)
    assert list(df.columns) == TSV_COLUMNS
    assert len(df) == 3
    assert set(df["transcript_id"]) == {"g1.t1", "g2.t1", "g3.t1"}


def test_score_annotation_bounds(conf_gff3_path, conf_hdf5_path):
    df = score_annotation(conf_gff3_path, conf_hdf5_path)
    for col in ("helixer_support", "helixer_locus_conf",
                "exon_confidence", "intron_coincidence_fraction"):
        assert (df[col] >= 0.0).all() and (df[col] <= 1.0).all()


def test_score_annotation_matches_transcript_scorer(conf_gff3_path, conf_hdf5_path, reader):
    df = score_annotation(conf_gff3_path, conf_hdf5_path)
    row = df[df["transcript_id"] == "g2.t1"].iloc[0]
    direct = score_transcript_confidence(_tx("chr2", "-", MINUS_EXONS), reader)
    assert row["helixer_support"] == pytest.approx(direct["helixer_support"])
    assert row["helixer_locus_conf"] == pytest.approx(direct["helixer_locus_conf"])


def test_region_filter_by_seqid(conf_gff3_path, conf_hdf5_path):
    df = score_annotation(conf_gff3_path, conf_hdf5_path, region="chr1")
    assert set(df["transcript_id"]) == {"g1.t1", "g3.t1"}


def test_region_filter_by_span_excludes_within_seqid(conf_gff3_path, conf_hdf5_path):
    # chr1:41-50 (internal [40,50)) overlaps g1 (10-50) but not g3 (10-20).
    df = score_annotation(conf_gff3_path, conf_hdf5_path, region="chr1:41-50")
    assert set(df["transcript_id"]) == {"g1.t1"}


def test_region_filter_other_seqid(conf_gff3_path, conf_hdf5_path):
    df = score_annotation(conf_gff3_path, conf_hdf5_path, region="chr2:21-40")
    assert set(df["transcript_id"]) == {"g2.t1"}


def test_score_annotation_split_matches_combined(conf_gff3_path, conf_hdf5_path, conf_split_paths):
    inp, pred = conf_split_paths
    combined = score_annotation(conf_gff3_path, conf_hdf5_path)
    split = score_annotation(conf_gff3_path, inp, predictions_h5=pred)
    assert split["helixer_support"].tolist() == pytest.approx(combined["helixer_support"].tolist())


def test_score_annotation_closes_reader_on_error(conf_gff3_path, conf_hdf5_path, monkeypatch):
    import helixforge.score.confidence as mod
    opened = {}
    real_init = HDF5ConfidenceReader.__init__

    def spy_init(self, *a, **k):
        real_init(self, *a, **k)
        opened["reader"] = self

    monkeypatch.setattr(HDF5ConfidenceReader, "__init__", spy_init)
    monkeypatch.setattr(mod, "GFF3Parser", lambda p: (_ for _ in ()).throw(RuntimeError("boom")))
    with pytest.raises(RuntimeError):
        score_annotation(conf_gff3_path, conf_hdf5_path)
    # finally-clause must have closed the reader (handles released).
    assert opened["reader"]._meta_h5 is None


# --- summary + TSV --------------------------------------------------------

def test_summarize_confidence_on_known_set(conf_gff3_path, conf_hdf5_path):
    df = score_annotation(conf_gff3_path, conf_hdf5_path)
    summary = summarize_confidence(df)
    assert summary["n_transcripts"] == 3
    s = df["helixer_support"]
    assert summary["mean_helixer_support"] == pytest.approx(float(s.mean()))
    assert summary["min_helixer_support"] == pytest.approx(float(s.min()))
    assert summary["max_helixer_support"] == pytest.approx(float(s.max()))
    assert 0.0 <= summary["mean_helixer_support"] <= 1.0


def test_summarize_low_confidence_count(conf_gff3_path, conf_hdf5_path):
    df = score_annotation(conf_gff3_path, conf_hdf5_path)
    cutoff = 0.95  # all three transcripts fall below this
    summary = summarize_confidence(df, low_confidence_cutoff=cutoff)
    assert summary["low_confidence_cutoff"] == cutoff
    assert summary["n_below_cutoff"] == int((df["helixer_support"] < cutoff).sum())


def test_summarize_empty_dataframe():
    import pandas as pd
    df = pd.DataFrame(columns=TSV_COLUMNS)
    summary = summarize_confidence(df)
    assert summary["n_transcripts"] == 0
    assert summary["mean_helixer_support"] is None
    assert summary["n_below_cutoff"] == 0


def test_write_confidence_tsv_roundtrip(conf_gff3_path, conf_hdf5_path, tmp_path):
    import pandas as pd
    df = score_annotation(conf_gff3_path, conf_hdf5_path)
    out = tmp_path / "conf.tsv"
    returned = write_confidence_tsv(df, out)
    assert str(returned) == str(out)
    back = pd.read_csv(out, sep="\t")
    assert list(back.columns) == TSV_COLUMNS
    assert len(back) == 3
