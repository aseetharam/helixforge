"""CRAM offline-reference threading (Phase 33 D1). Floor: 8.

A CRAM decodes only against its genome FASTA — without it pysam attempts a
remote ENA MD5 fetch that fails offline. These tests build a tiny CRAM from a
BAM against a *local* reference (no network), and assert:

- ``reference_filename`` is threaded through every BAM/CRAM opener the evidence
  scorer uses (``_bam_targets`` / ``collect_junctions`` / ``score_annotation``)
  and through ``overall_mapping_rate``;
- the CRAM pileup path passes ``multiple_iterators=False`` (htslib does not
  implement multiple iterators for CRAM) while the BAM path is left unchanged;
- a CRAM and its source BAM yield identical junctions/coverage — on **both
  strands** — and no remote fetch is attempted when a reference is given.

Concrete literal coordinates throughout; the fixture carries a ``+`` junction at
(150, 250) and a ``-`` junction at (350, 450).
"""

import warnings

import pysam
import pytest

from helixforge.io.bam import (
    CoverageCalculator,
    JunctionExtractor,
    _is_cram,
    _open_alignment,
    _pileup_kwargs,
    overall_mapping_rate,
)

_BAM_HDR_600 = {"HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "chr1", "LN": 600}]}
_QUERY_OPS_T = {0, 1, 4, 7, 8}


def _build_simple_bam(path, reads):
    reads = sorted(reads, key=lambda r: r["start"])
    with pysam.AlignmentFile(path, "wb", header=_BAM_HDR_600) as af:
        for spec in reads:
            a = pysam.AlignedSegment(af.header)
            a.query_name = spec["name"]
            a.reference_id = 0
            a.reference_start = spec["start"]
            a.mapping_quality = spec["mapq"]
            a.cigartuples = spec["cigar"]
            qlen = sum(l for op, l in spec["cigar"] if op in _QUERY_OPS_T)
            a.query_sequence = "A" * qlen
            a.flag = spec.get("flag", 0)
            if spec.get("xs"):
                a.set_tag("XS", spec["xs"], "A")
            af.write(a)
    pysam.index(str(path))
    return str(path)


def _ref_600(tmp_path):
    p = tmp_path / "ref600.fa"
    seq = ("ACGT" * 150)  # exactly 600 bp
    with open(p, "w") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return str(p)


def _xs_reads_600():
    reads = []
    for i in range(5):  # donor 150 acceptor 250 (+, 5 reads)
        reads.append(dict(name=f"p{i}", start=100, cigar=[(0, 50), (3, 100), (0, 50)],
                          mapq=60, xs="+"))
    for i in range(3):  # donor 350 acceptor 450 (-, 3 reads)
        reads.append(dict(name=f"m{i}", start=300, cigar=[(0, 50), (3, 100), (0, 50)],
                          mapq=60, xs="-"))
    return reads


def _bam_to_cram(bam_path, cram_path, reference_filename):
    with pysam.AlignmentFile(bam_path, "rb") as src:
        with pysam.AlignmentFile(
            cram_path, "wc", template=src, reference_filename=reference_filename
        ) as dst:
            for read in src.fetch(until_eof=True):
                dst.write(read)
    pysam.index(str(cram_path))
    return str(cram_path)


@pytest.fixture
def bam_cram_pair(tmp_path):
    ref = _ref_600(tmp_path)
    bam = _build_simple_bam(tmp_path / "pair.bam", _xs_reads_600())
    cram = _bam_to_cram(bam, str(tmp_path / "pair.cram"), ref)
    return bam, cram, ref


# --- _pileup_kwargs branches ------------------------------------------------

def test_pileup_kwargs_cram_disables_multiple_iterators(bam_cram_pair):
    _bam, cram, ref = bam_cram_pair
    with _open_alignment(cram, reference_filename=ref) as handle:
        assert _pileup_kwargs(handle) == {"multiple_iterators": False}


def test_pileup_kwargs_bam_unchanged(bam_cram_pair):
    bam, _cram, _ref = bam_cram_pair
    with _open_alignment(bam) as handle:
        # BAM keeps pysam's default (no key injected) → behaviour byte-identical.
        assert _pileup_kwargs(handle) == {}


def test_cram_coverage_emits_no_multiple_iterators_warning(bam_cram_pair):
    bam, cram, ref = bam_cram_pair
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any pysam UserWarning becomes a failure
        with CoverageCalculator.from_bam(cram, reference_filename=ref) as cov:
            cm = cov.mean_coverage("chr1", 100, 150)
    with CoverageCalculator.from_bam(bam) as bam_cov:
        bm = bam_cov.mean_coverage("chr1", 100, 150)
    assert cm == pytest.approx(bm)
    assert cm > 0


# --- evidence-scorer opener threading --------------------------------------

def test_bam_targets_decodes_cram_with_reference(bam_cram_pair):
    from helixforge.score.evidence import _bam_targets

    _bam, cram, ref = bam_cram_pair
    assert _is_cram(cram)
    targets = _bam_targets([cram], None, None, None, reference_filename=ref)
    assert targets == [("chr1", 0, 600)]


def test_collect_junctions_cram_equals_bam_both_strands(bam_cram_pair):
    from helixforge.score.evidence import collect_junctions

    bam, cram, ref = bam_cram_pair
    bam_j = {(j.seqid, j.donor, j.acceptor, j.strand): j.read_count
             for j in collect_junctions(bam_paths=[bam])}
    cram_j = {(j.seqid, j.donor, j.acceptor, j.strand): j.read_count
              for j in collect_junctions(bam_paths=[cram], reference_filename=ref)}
    assert cram_j == bam_j
    assert ("chr1", 150, 250, "+") in cram_j  # plus-strand junction
    assert ("chr1", 350, 450, "-") in cram_j  # minus-strand junction


def test_score_annotation_threads_reference_to_openers(bam_cram_pair, monkeypatch):
    import helixforge.io.bam as bam_mod

    bam, cram, ref = bam_cram_pair
    seen: list[str | None] = []
    real_open = bam_mod._open_alignment

    def spy(path, *, reference_filename=None):
        seen.append(reference_filename)
        return real_open(path, reference_filename=reference_filename)

    monkeypatch.setattr(bam_mod, "_open_alignment", spy)

    # A 1-transcript GFF3 over chr1 so score_annotation has something to score.
    import textwrap
    gff = str(bam.replace("pair.bam", "ann.gff3"))
    with open(gff, "w") as fh:
        fh.write(textwrap.dedent("""\
            chr1\t.\tgene\t101\t500\t.\t+\t.\tID=g1
            chr1\t.\tmRNA\t101\t500\t.\t+\t.\tID=t1;Parent=g1
            chr1\t.\texon\t101\t150\t.\t+\t.\tID=e1;Parent=t1
            chr1\t.\texon\t251\t500\t.\t+\t.\tID=e2;Parent=t1
            """))

    from helixforge.score.evidence import score_annotation

    df = score_annotation(gff, bam_paths=[cram], reference_filename=ref)
    assert len(df) == 1
    # Every CRAM open saw the reference → no remote ENA fetch was possible.
    assert seen  # at least one opener fired
    assert all(r == ref for r in seen)


def test_remote_fetch_not_attempted_when_reference_given(bam_cram_pair, monkeypatch):
    import helixforge.io.bam as bam_mod

    _bam, cram, ref = bam_cram_pair
    captured: list[dict] = []
    real_af = pysam.AlignmentFile

    def fake_af(path, mode="r", **kwargs):
        if "c" in mode:  # CRAM open must carry the local reference
            captured.append({"path": str(path), "ref": kwargs.get("reference_filename")})
        return real_af(path, mode, **kwargs)

    monkeypatch.setattr(bam_mod.pysam, "AlignmentFile", fake_af)
    with JunctionExtractor(cram, reference_filename=ref) as je:
        juncs = je.extract_junctions("chr1", 0, 600)
    assert juncs
    assert captured and all(c["ref"] == ref for c in captured)


def test_overall_mapping_rate_threads_reference(bam_cram_pair, monkeypatch):
    import helixforge.io.bam as bam_mod

    _bam, cram, ref = bam_cram_pair
    seen: list[str | None] = []
    real_stats = bam_mod.bam_mapping_stats

    def spy(path, *, reference_filename=None):
        seen.append(reference_filename)
        return real_stats(path, reference_filename=reference_filename)

    monkeypatch.setattr(bam_mod, "bam_mapping_stats", spy)
    out = overall_mapping_rate([cram], reference_filename=ref)
    assert out is not None
    assert seen == [ref]
