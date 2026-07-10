"""Tests for BAM/STAR evidence I/O (Phase 2). Floor: 24.

The synthetic BAM fixture is validated first, then junctions (both strands,
filters, multisample), coverage, STAR SJ conversion, and the bigWig ImportError.
"""

import builtins

import pysam
import pytest

from helixforge.io.bam import (
    CoverageCalculator,
    JunctionExtractor,
    bam_mapping_stats,
    overall_mapping_rate,
    parse_star_sj_tab,
)


# --------------------------------------------------------------------------
# Validate the BAM fixture itself first
# --------------------------------------------------------------------------

def test_bam_fixture_is_indexed(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as af:
        assert af.has_index()


def test_bam_fixture_read_count(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as af:
        # 5 plus + 3 minus + lowmapq + shortoh + secondary on chr1, noxs on chr2 = 12
        assert af.count(until_eof=True) == 12


def test_bam_fixture_chr1_has_reads(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as af:
        assert af.count("chr1", 100, 300) > 0


# --------------------------------------------------------------------------
# JunctionExtractor construction
# --------------------------------------------------------------------------

def test_junction_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        JunctionExtractor("/no/such.bam")


def test_junction_missing_index_raises(bam_no_index):
    with pytest.raises(ValueError):
        JunctionExtractor(bam_no_index)


# --------------------------------------------------------------------------
# Junction extraction
# --------------------------------------------------------------------------

def test_extract_junctions_default_count(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000)
        # plus (150,250) and minus (450,550); lowmapq/shortoh/secondary excluded
        assert len(juncs) == 2


def test_extract_junctions_plus_coordinates(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000)
        plus = next(j for j in juncs if j.strand == "+")
        assert plus.donor == 150
        assert plus.acceptor == 250


def test_extract_junctions_plus_read_count(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000)
        plus = next(j for j in juncs if j.strand == "+")
        assert plus.read_count == 5


def test_extract_junctions_minus_strand(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000)
        minus = next(j for j in juncs if j.strand == "-")
        assert minus.donor == 450
        assert minus.acceptor == 550
        assert minus.read_count == 3


def test_extract_junctions_intron_length(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000)
        assert all(j.intron_length == 100 for j in juncs)


def test_extract_junctions_samples_default_one(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000)
        assert all(j.samples == 1 for j in juncs)


def test_mapq_filter_excludes_low_mapq(bam_path):
    with JunctionExtractor(bam_path) as je:
        # default min_mapq=10 excludes the mapq=5 read at donor 650
        assert all(j.donor != 650 for j in je.extract_junctions("chr1", 0, 1000))


def test_mapq_filter_relaxed_includes_low_mapq(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000, min_mapq=1)
        assert any(j.donor == 650 for j in juncs)


def test_overhang_filter_excludes_short(bam_path):
    with JunctionExtractor(bam_path) as je:
        # default min_overhang=8 excludes the 5M-flanked junction at donor 805
        assert all(j.donor != 805 for j in je.extract_junctions("chr1", 0, 1000))


def test_overhang_filter_relaxed_includes_short(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000, min_overhang=3)
        assert any(j.donor == 805 for j in juncs)


def test_secondary_read_excluded(bam_path):
    # secondary read starts at 110 -> junction donor 160; must not appear
    with JunctionExtractor(bam_path) as je:
        assert all(j.donor != 160 for j in je.extract_junctions("chr1", 0, 1000))


def test_unknown_strand_excluded(bam_path):
    # the chr2 read has no XS tag -> excluded
    with JunctionExtractor(bam_path) as je:
        assert je.extract_junctions("chr2", 0, 1000) == []


def test_extract_junctions_sorted(bam_path):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions("chr1", 0, 1000)
        donors = [j.donor for j in juncs]
        assert donors == sorted(donors)


# --------------------------------------------------------------------------
# Multisample merge
# --------------------------------------------------------------------------

def test_multisample_merge_read_count_summed(bam_path, bam_path2):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions_multisample([bam_path, bam_path2], "chr1", 0, 1000)
        plus = next(j for j in juncs if j.donor == 150)
        # sample1: 5, sample2: 2 -> 7
        assert plus.read_count == 7


def test_multisample_merge_samples_counted(bam_path, bam_path2):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions_multisample([bam_path, bam_path2], "chr1", 0, 1000)
        plus = next(j for j in juncs if j.donor == 150)
        assert plus.samples == 2


def test_multisample_unique_junction_single_sample(bam_path, bam_path2):
    with JunctionExtractor(bam_path) as je:
        juncs = je.extract_junctions_multisample([bam_path, bam_path2], "chr1", 0, 1000)
        # sample2-only junction donor 350
        uniq = next(j for j in juncs if j.donor == 350)
        assert uniq.read_count == 4
        assert uniq.samples == 1


# --------------------------------------------------------------------------
# Coverage
# --------------------------------------------------------------------------

def test_coverage_known_depth(bam_path):
    # 5 plus reads each 50M at [100,150) -> depth 5
    with CoverageCalculator.from_bam(bam_path) as cov:
        assert cov.mean_coverage("chr1", 100, 150) == pytest.approx(5.0)


def test_coverage_zero_when_no_reads(bam_path):
    # [970,1000) is past every aligned block (max M end is 955)
    with CoverageCalculator.from_bam(bam_path) as cov:
        assert cov.mean_coverage("chr1", 970, 1000) == 0.0


def test_coverage_zero_in_intron_gap(bam_path):
    # [150,250) is spanned by N (no aligned bases) -> 0
    with CoverageCalculator.from_bam(bam_path) as cov:
        assert cov.mean_coverage("chr1", 150, 250) == 0.0


def test_coverage_array_length(bam_path):
    with CoverageCalculator.from_bam(bam_path) as cov:
        arr = cov.region_coverage_array("chr1", 100, 150)
        assert len(arr) == 50
        assert all(d == 5 for d in arr)


def test_coverage_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        CoverageCalculator.from_bam("/no/such.bam")


def test_coverage_multisample_mean(bam_path, bam_path2):
    # sample1 depth 5 at [100,150); sample2 depth 2 -> mean of means = 3.5
    mean = CoverageCalculator.mean_coverage_multisample(
        [bam_path, bam_path2], "chr1", 100, 150
    )
    assert mean == pytest.approx(3.5)


# --------------------------------------------------------------------------
# region_coverage_arrays — batched, one pileup per merged span, byte-identical
# to per-region region_coverage_array. Both strands of the fixture are covered:
# chr1 plus reads cover [100,150)+[250,300); minus reads cover [400,450)+[550,600).
# --------------------------------------------------------------------------

def test_region_coverage_arrays_matches_per_region(bam_path):
    regions = [
        ("chr1", 100, 150),  # plus exon, depth 5
        ("chr1", 250, 300),  # plus exon, depth 5
        ("chr1", 400, 450),  # minus exon, depth 3
        ("chr1", 550, 600),  # minus exon, depth 3
        ("chr1", 150, 250),  # intron gap, depth 0
    ]
    with CoverageCalculator.from_bam(bam_path) as cov:
        batched = cov.region_coverage_arrays(regions)
        for seqid, s, e in regions:
            per_region = cov.region_coverage_array(seqid, s, e)
            assert list(batched[(seqid, s, e)]) == per_region


def test_region_coverage_arrays_known_depths(bam_path):
    with CoverageCalculator.from_bam(bam_path) as cov:
        out = cov.region_coverage_arrays(
            [("chr1", 100, 150), ("chr1", 400, 450), ("chr1", 150, 250)]
        )
    assert list(out[("chr1", 100, 150)]) == [5.0] * 50  # plus, depth 5
    assert list(out[("chr1", 400, 450)]) == [3.0] * 50  # minus, depth 3
    assert list(out[("chr1", 150, 250)]) == [0.0] * 100  # intron gap


def test_region_coverage_arrays_one_pileup_per_merged_span(bam_path, monkeypatch):
    # [100,150) and [120,200) overlap -> one merged span -> ONE pileup; the distant
    # [400,450) is a second span. So 3 regions collapse to 2 pileup calls.
    calls = []
    orig = CoverageCalculator.region_coverage_array

    def spy(self, seqid, start, end, **kw):
        calls.append((seqid, start, end))
        return orig(self, seqid, start, end, **kw)

    monkeypatch.setattr(CoverageCalculator, "region_coverage_array", spy)
    with CoverageCalculator.from_bam(bam_path) as cov:
        cov.region_coverage_arrays(
            [("chr1", 100, 150), ("chr1", 120, 200), ("chr1", 400, 450)]
        )
    assert calls == [("chr1", 100, 200), ("chr1", 400, 450)]


def test_region_coverage_arrays_rejects_empty_region(bam_path):
    with CoverageCalculator.from_bam(bam_path) as cov:
        with pytest.raises(ValueError):
            cov.region_coverage_arrays([("chr1", 150, 150)])


# --------------------------------------------------------------------------
# STAR SJ.out.tab
# --------------------------------------------------------------------------

def test_star_sj_count(star_sj_path):
    juncs = parse_star_sj_tab(star_sj_path)
    # (150,250,+) and (450,550,-) pass; others filtered
    assert len(juncs) == 2


def test_star_sj_plus_conversion(star_sj_path):
    juncs = parse_star_sj_tab(star_sj_path)
    plus = next(j for j in juncs if j.strand == "+")
    assert plus.donor == 150
    assert plus.acceptor == 250
    assert plus.read_count == 8


def test_star_sj_minus_conversion(star_sj_path):
    juncs = parse_star_sj_tab(star_sj_path)
    minus = next(j for j in juncs if j.strand == "-")
    assert minus.donor == 450
    assert minus.acceptor == 550


def test_star_sj_excludes_undefined_strand(star_sj_path):
    juncs = parse_star_sj_tab(star_sj_path)
    assert all(j.donor != 700 for j in juncs)


def test_star_sj_min_unique_reads_filter(star_sj_path):
    juncs = parse_star_sj_tab(star_sj_path)
    # the unique=2 junction (900,1000) is filtered
    assert all(j.read_count >= 3 for j in juncs)


def test_star_sj_overhang_filter(star_sj_path):
    juncs = parse_star_sj_tab(star_sj_path)
    # the overhang=5 junction (donor 99) is filtered
    assert all(j.donor != 99 for j in juncs)


def test_star_sj_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        parse_star_sj_tab("/no/such/SJ.out.tab")


# --------------------------------------------------------------------------
# bigWig optionality
# --------------------------------------------------------------------------

def test_bigwig_import_error(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pyBigWig":
            raise ImportError("no pyBigWig")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="helixforge\\[bigwig\\]"):
        CoverageCalculator.from_bigwig("/whatever.bw")


# ==========================================================================
# Phase 27 D4 — N-masked coverage denominators (assessment §1.7).
# A gene spanning an assembly gap (reference N-run, depth 0) must not have its
# mean coverage diluted by the gap; the N positions are excluded from the
# denominator. With no mask the behaviour is byte-identical to before.
# ==========================================================================

class _ArrayCoverage(CoverageCalculator):
    """CoverageCalculator over a fixed per-base depth array (no BAM needed)."""

    def __init__(self, depth):
        self._handle = None
        self.source_type = "array"
        self._depth = list(depth)

    def region_coverage_array(self, seqid, start, end, **kwargs):
        # Phase 28 D3: mean_coverage now forwards min_mapq/max_depth; the fake
        # ignores them (a fixed array has no per-read MAPQ notion).
        return list(self._depth)


class _MaskGenome:
    """Genome accessor whose plus-strand sequence is a fixed string."""

    def __init__(self, seq):
        self._seq = seq

    def get_sequence(self, seqid, start, end, strand="+"):
        return self._seq[start:end]


def test_mean_coverage_no_mask_unchanged():
    # Baseline: half the region is an assembly gap (depth 0). Without a mask the
    # mean is diluted: (10*5 + 0*5) / 10 = 5.0.
    cov = _ArrayCoverage([10] * 5 + [0] * 5)
    assert cov.mean_coverage("chr1", 0, 10) == pytest.approx(5.0)


def test_mean_coverage_genome_nmask_excludes_gap():
    # The last 5 bases are reference 'N' (the gap). Excluding them, the mean is
    # taken over the 5 covered bases only: 50/5 = 10.0 (not 5.0).
    cov = _ArrayCoverage([10] * 5 + [0] * 5)
    genome = _MaskGenome("ACGTA" + "NNNNN")
    assert cov.mean_coverage("chr1", 0, 10, genome=genome) == pytest.approx(10.0)


def test_mean_coverage_explicit_nmask_excludes_gap():
    cov = _ArrayCoverage([10] * 5 + [0] * 5)
    mask = [False] * 5 + [True] * 5
    assert cov.mean_coverage("chr1", 0, 10, n_mask=mask) == pytest.approx(10.0)


def test_mean_coverage_no_n_case_identical_with_genome():
    # No reference N anywhere: passing a genome must not change the result.
    cov = _ArrayCoverage([4] * 10)
    genome = _MaskGenome("ACGTACGTAC")
    assert cov.mean_coverage("chr1", 0, 10) == pytest.approx(4.0)
    assert cov.mean_coverage("chr1", 0, 10, genome=genome) == pytest.approx(4.0)


def test_mean_coverage_whole_region_is_gap_returns_zero():
    # Entire region is reference N -> no mappable bases -> 0.0 (not a div-by-zero).
    cov = _ArrayCoverage([0] * 6)
    genome = _MaskGenome("NNNNNN")
    assert cov.mean_coverage("chr1", 0, 6, genome=genome) == 0.0


# ==========================================================================
# Phase 28 — splice canonicity, strand-from-motif, read-filtering policy
# ==========================================================================

from helixforge.io.bam import (  # noqa: E402
    classify_splice_motif,
    infer_strand_from_motif,
)
from helixforge.io.fasta import GenomeAccessor  # noqa: E402

# --- D1: classify_splice_motif (both strands) ---


def test_classify_motif_gt_ag_plus():
    # Genomic-forward GT..AG on the + strand is the canonical major motif.
    assert classify_splice_motif("GT", "AG", "+") == "GT-AG"


def test_classify_motif_gt_ag_minus():
    # A minus-strand GT-AG intron reads as CT..AC on the forward genome strand
    # (revcomp of AG..GT). The classifier must orient by strand.
    assert classify_splice_motif("CT", "AC", "-") == "GT-AG"
    # ...and the same forward bases on the + strand are NOT GT-AG.
    assert classify_splice_motif("CT", "AC", "+") == "non-canonical"


def test_classify_motif_gc_ag_both_strands():
    assert classify_splice_motif("GC", "AG", "+") == "GC-AG"
    # minus GC-AG: forward = revcomp(AG)=CT donor, revcomp(GC)=GC acceptor
    assert classify_splice_motif("CT", "GC", "-") == "GC-AG"


def test_classify_motif_at_ac_both_strands():
    assert classify_splice_motif("AT", "AC", "+") == "AT-AC"
    # minus AT-AC: forward = revcomp(AC)=GT donor, revcomp(AT)=AT acceptor
    assert classify_splice_motif("GT", "AT", "-") == "AT-AC"


def test_classify_motif_non_canonical_and_n_tolerant():
    assert classify_splice_motif("AA", "TT", "+") == "non-canonical"
    assert classify_splice_motif("AA", "TT", "-") == "non-canonical"
    # An N/IUPAC dinucleotide never raises — it is simply non-canonical.
    assert classify_splice_motif("NN", "AG", "+") == "non-canonical"


# --- D1: STAR motif column preserved + genome cross-check ---


def _write_one_seq_fasta(path, seqid, seq):
    with open(path, "w") as fh:
        fh.write(f">{seqid}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return str(path)


def _motif_genome(tmp_path):
    """A chr1 of length 600 with canonical/non-canonical dinucleotides embedded.

    Junction A (+ GT-AG):   donor 150 -> 'GT', acceptor 250 -> 'AG'
    Junction B (- GT-AG):   donor 350 -> 'CT', acceptor 450 -> 'AC'
    Junction C (non-canon): donor 490 -> 'AA', acceptor 540 -> 'TT'
    """
    seq = list("A" * 600)

    def put(pos, dn):
        seq[pos], seq[pos + 1] = dn[0], dn[1]

    put(150, "GT"); put(248, "AG")   # A donor [150:152], acceptor [248:250]
    put(350, "CT"); put(448, "AC")   # B donor [350:352], acceptor [448:450]
    put(490, "AA"); put(538, "TT")   # C donor [490:492], acceptor [538:540]
    p = _write_one_seq_fasta(tmp_path / "motif.fa", "chr1", "".join(seq))
    return GenomeAccessor(p)


def test_star_motif_preserved_as_canonical(star_sj_path):
    # The STAR fixture's first kept row has motif code 1 (GT/AG), the second
    # motif code 2 (CT/AC = minus GT-AG). Both map to the GT-AG canonical class.
    juncs = parse_star_sj_tab(star_sj_path)
    by_pos = {(j.donor, j.acceptor): j for j in juncs}
    assert by_pos[(150, 250)].canonical == "GT-AG"
    assert by_pos[(450, 550)].canonical == "GT-AG"
    # The multimap column (col 7) is carried too: row 1 has 2 multi-mappers.
    assert by_pos[(150, 250)].multimap_reads == 2
    assert by_pos[(450, 550)].multimap_reads == 0


def test_star_motif_genome_cross_check(tmp_path):
    # A STAR row claims motif 1 (GT-AG) at donor 150/acceptor 250 on '+'; the
    # genome agrees (GT..AG) -> canonical stays GT-AG. A second row claims motif
    # 1 at a position the genome shows as non-canonical (AA..TT) -> genome wins.
    sj = tmp_path / "SJ.out.tab"
    sj.write_text(
        "chr1\t151\t250\t1\t1\t0\t9\t0\t30\n"   # genome GT..AG -> GT-AG
        "chr1\t491\t540\t1\t1\t0\t9\t0\t30\n"   # genome AA..TT -> non-canonical
    )
    genome = _motif_genome(tmp_path)
    juncs = parse_star_sj_tab(str(sj), genome=genome)
    by_pos = {(j.donor, j.acceptor): j for j in juncs}
    assert by_pos[(150, 250)].canonical == "GT-AG"
    # STAR said GT-AG but the genome is the ground truth -> non-canonical.
    assert by_pos[(490, 540)].canonical == "non-canonical"


# --- D2: strand inference from motif when XS absent ---

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


# 50M 100N 50M from start 100 -> donor 150, acceptor 250 (junction A, + GT-AG)
# 50M 100N 50M from start 300 -> donor 350, acceptor 450 (junction B, - GT-AG)
# 30M 50N 30M  from start 460 -> donor 490, acceptor 540 (junction C, non-canon)
def _xsless_reads():
    reads = []
    for i in range(4):
        reads.append(dict(name=f"a{i}", start=100, cigar=[(0, 50), (3, 100), (0, 50)], mapq=60))
    for i in range(4):
        reads.append(dict(name=f"b{i}", start=300, cigar=[(0, 50), (3, 100), (0, 50)], mapq=60))
    for i in range(4):
        reads.append(dict(name=f"c{i}", start=460, cigar=[(0, 30), (3, 50), (0, 30)], mapq=60))
    return reads


def test_xsless_read_strand_from_motif(tmp_path):
    # No read carries XS; with the default xs_then_motif policy + a genome, the
    # canonical motif orients them: junction A -> '+', junction B -> '-'.
    bam = _build_simple_bam(tmp_path / "xsless.bam", _xsless_reads())
    genome = _motif_genome(tmp_path)
    with JunctionExtractor(bam) as je:
        juncs = je.extract_junctions("chr1", 0, 600, genome=genome)
    by_pos = {(j.donor, j.acceptor): j for j in juncs}
    assert by_pos[(150, 250)].strand == "+"
    assert by_pos[(150, 250)].canonical == "GT-AG"
    assert by_pos[(350, 450)].strand == "-"
    assert by_pos[(350, 450)].canonical == "GT-AG"


def test_xsless_non_canonical_dropped(tmp_path):
    # Junction C is non-canonical and the reads have no XS -> cannot be oriented
    # -> dropped (assessment §1.4: only discard when motif non-canonical AND no XS).
    bam = _build_simple_bam(tmp_path / "xsless.bam", _xsless_reads())
    genome = _motif_genome(tmp_path)
    with JunctionExtractor(bam) as je:
        juncs = je.extract_junctions("chr1", 0, 600, genome=genome)
    assert all((j.donor, j.acceptor) != (490, 540) for j in juncs)


def test_strand_source_policy_respected(tmp_path):
    # An XS='+' read over junction B (motif is minus-canonical). xs_then_motif
    # trusts XS ('+'); strand_source='motif' ignores XS and uses the motif ('-').
    reads = [dict(name="x", start=300, cigar=[(0, 50), (3, 100), (0, 50)], mapq=60, xs="+")]
    bam = _build_simple_bam(tmp_path / "policy.bam", reads)
    genome = _motif_genome(tmp_path)
    with JunctionExtractor(bam) as je:
        xs_then = je.extract_junctions("chr1", 0, 600, genome=genome,
                                       strand_source="xs_then_motif")
        motif_only = je.extract_junctions("chr1", 0, 600, genome=genome,
                                          strand_source="motif")
        xs_only = je.extract_junctions("chr1", 0, 600, genome=genome,
                                       strand_source="xs")
    assert xs_then[0].strand == "+"      # XS wins
    assert motif_only[0].strand == "-"   # motif wins
    assert xs_only[0].strand == "+"      # XS only


def test_xs_only_drops_xsless_even_with_genome(tmp_path):
    # strand_source='xs' is the legacy strict policy: XS-less reads dropped even
    # when a genome could orient them.
    bam = _build_simple_bam(tmp_path / "xsless.bam", _xsless_reads())
    genome = _motif_genome(tmp_path)
    with JunctionExtractor(bam) as je:
        juncs = je.extract_junctions("chr1", 0, 600, genome=genome, strand_source="xs")
    assert juncs == []


def test_infer_strand_from_motif_helper():
    assert infer_strand_from_motif("GT", "AG") == "+"
    assert infer_strand_from_motif("CT", "AC") == "-"
    assert infer_strand_from_motif("AA", "CC") is None


# --- D3: read-filtering policy (MAPQ / secondary / dup / max_depth) ---

def _cov_reads():
    """20 primary + 1 secondary + 1 duplicate, all 30M from pos 100; plus 5
    low-MAPQ (mapq 3) reads at pos 100 to exercise the unique-vs-all split."""
    reads = []
    for i in range(20):
        reads.append(dict(name=f"u{i}", start=100, cigar=[(0, 30)], mapq=60))
    for i in range(5):
        reads.append(dict(name=f"lq{i}", start=100, cigar=[(0, 30)], mapq=3))
    reads.append(dict(name="sec", start=100, cigar=[(0, 30)], mapq=60, flag=256))
    reads.append(dict(name="dup", start=100, cigar=[(0, 30)], mapq=60, flag=1024))
    return reads


def test_coverage_excludes_secondary_and_duplicate(tmp_path):
    # 20 unique-MAPQ + 5 low-MAPQ primaries = 25 counted; secondary + duplicate
    # are excluded by stepper='all'. (min_mapq=0 default counts all primaries.)
    bam = _build_simple_bam(tmp_path / "cov.bam", _cov_reads())
    with CoverageCalculator.from_bam(bam) as cov:
        arr = cov.region_coverage_array("chr1", 110, 120)
    assert arr[0] == 25


def test_coverage_min_mapq_honored(tmp_path):
    # With min_mapq=UNIQUE_MIN_MAPQ (10) the 5 mapq-3 reads drop -> 20 remain.
    from helixforge.constants import UNIQUE_MIN_MAPQ
    bam = _build_simple_bam(tmp_path / "cov.bam", _cov_reads())
    with CoverageCalculator.from_bam(bam) as cov:
        arr = cov.region_coverage_array("chr1", 110, 120, min_mapq=UNIQUE_MIN_MAPQ)
    assert arr[0] == 20


def test_coverage_max_depth_not_silently_capping(tmp_path):
    # 25 primary reads stack at pos 100. A tiny max_depth caps the pileup; the
    # large default does not -> all 25 counted (no silent htslib ~8000 cap).
    bam = _build_simple_bam(tmp_path / "cov.bam", _cov_reads())
    with CoverageCalculator.from_bam(bam) as cov:
        capped = cov.region_coverage_array("chr1", 110, 120, max_depth=5)
        full = cov.region_coverage_array("chr1", 110, 120)
    assert capped[0] <= 5
    assert full[0] == 25


def test_coverage_unique_vs_all_modes_differ(tmp_path):
    # region_coverage_both reports all (25) and unique-only (20, MAPQ>=10) in one
    # pass; the two differ exactly by the 5 multi-mapper-grade (MAPQ 3) reads.
    bam = _build_simple_bam(tmp_path / "cov.bam", _cov_reads())
    with CoverageCalculator.from_bam(bam) as cov:
        all_arr, uniq_arr = cov.region_coverage_both("chr1", 110, 120)
        all_mean, uniq_mean = cov.mean_coverage_both("chr1", 110, 120)
    assert all_arr[0] == 25
    assert uniq_arr[0] == 20
    assert all_mean == pytest.approx(25.0)
    assert uniq_mean == pytest.approx(20.0)


# --------------------------------------------------------------------------
# Mapping-rate stats (Phase 31 D4) — the synthetic fixture is fully mapped
# --------------------------------------------------------------------------

def test_bam_mapping_stats_fully_mapped(bam_path):
    s = bam_mapping_stats(bam_path)
    assert s["unmapped"] == 0
    assert s["mapped"] > 0
    assert s["mapping_rate"] == pytest.approx(1.0)


def test_bam_mapping_stats_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        bam_mapping_stats("/no/such/file.bam")


def test_overall_mapping_rate_aggregates(bam_path, bam_path2):
    agg = overall_mapping_rate([bam_path, bam_path2])
    assert agg is not None and agg["num_bams"] == 2
    assert agg["mapping_rate"] == pytest.approx(1.0)


def test_overall_mapping_rate_empty_is_none():
    assert overall_mapping_rate([]) is None


# --------------------------------------------------------------------------
# Phase 32 D2: CRAM input (detect by extension/magic; reference_filename) — a
# CRAM and its source BAM give identical junctions/coverage on a paired fixture.
# --------------------------------------------------------------------------

def _ref_600(tmp_path):
    """A 600 bp chr1 reference FASTA matching _BAM_HDR_600 (for CRAM decode)."""
    p = tmp_path / "ref600.fa"
    seq = ("ACGT" * 150)  # exactly 600 bp
    with open(p, "w") as fh:
        fh.write(">chr1\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return str(p)


def _xs_reads_600():
    # 50M 100N 50M from 100 -> donor 150 acceptor 250 (+, 5 reads)
    # 50M 100N 50M from 300 -> donor 350 acceptor 450 (-, 3 reads)
    reads = []
    for i in range(5):
        reads.append(dict(name=f"p{i}", start=100, cigar=[(0, 50), (3, 100), (0, 50)],
                          mapq=60, xs="+"))
    for i in range(3):
        reads.append(dict(name=f"m{i}", start=300, cigar=[(0, 50), (3, 100), (0, 50)],
                          mapq=60, xs="-"))
    return reads


def _bam_to_cram(bam_path, cram_path, reference_filename):
    """Transcode a BAM to an indexed CRAM against ``reference_filename``."""
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
    """(bam, cram, reference) carrying the same alignments (Phase 32 D2)."""
    ref = _ref_600(tmp_path)
    bam = _build_simple_bam(tmp_path / "pair.bam", _xs_reads_600())
    cram = _bam_to_cram(bam, str(tmp_path / "pair.cram"), ref)
    return bam, cram, ref


def test_cram_detected_and_opens(bam_cram_pair):
    from helixforge.io.bam import _is_cram

    _bam, cram, ref = bam_cram_pair
    assert _is_cram(cram)
    with JunctionExtractor(cram, reference_filename=ref) as je:
        juncs = je.extract_junctions("chr1", 0, 600)
    assert juncs  # decoded fine


def test_cram_junctions_equal_bam(bam_cram_pair):
    bam, cram, ref = bam_cram_pair
    with JunctionExtractor(bam) as je:
        bam_juncs = {(j.donor, j.acceptor, j.strand): j.read_count
                     for j in je.extract_junctions("chr1", 0, 600)}
    with JunctionExtractor(cram, reference_filename=ref) as je:
        cram_juncs = {(j.donor, j.acceptor, j.strand): j.read_count
                      for j in je.extract_junctions("chr1", 0, 600)}
    assert cram_juncs == bam_juncs
    # both strands present in the fixture
    assert (150, 250, "+") in cram_juncs
    assert (350, 450, "-") in cram_juncs


@pytest.mark.filterwarnings(
    "ignore:multiple_iterators not implemented for CRAM:UserWarning"
)
def test_cram_coverage_equals_bam(bam_cram_pair):
    bam, cram, ref = bam_cram_pair
    with CoverageCalculator.from_bam(bam) as bam_cov:
        bm = bam_cov.mean_coverage("chr1", 100, 150)
    with CoverageCalculator.from_bam(cram, reference_filename=ref) as cram_cov:
        cm = cram_cov.mean_coverage("chr1", 100, 150)
    assert cm == pytest.approx(bm)
    assert cm > 0


def test_cram_missing_crai_errors_clearly(tmp_path):
    ref = _ref_600(tmp_path)
    bam = _build_simple_bam(tmp_path / "noidx.bam", _xs_reads_600())
    cram = str(tmp_path / "noidx.cram")
    with pysam.AlignmentFile(bam, "rb") as src:
        with pysam.AlignmentFile(
            cram, "wc", template=src, reference_filename=ref
        ) as dst:
            for read in src.fetch(until_eof=True):
                dst.write(read)
    # deliberately NOT indexed -> clear index-missing error
    with pytest.raises(ValueError, match="index"):
        JunctionExtractor(cram, reference_filename=ref)
