"""Seqid + length concordance tests (biology.v4 Phase 26 D1). Floor: ≥6.

The #1 silent multi-tool genomics failure is a seqid naming mismatch (`Chr1`
vs `1` vs `chr1`); here it becomes a loud error. Disjoint sets → error; partial
overlap → warning naming the missing seqids; shared-contig length mismatch →
error; an alias map resolves a deliberate translation. Synthetic fixtures only.
"""

import pysam
import pytest

from helixforge.prep.preflight import check_reference_concordance


def _genome(tmp_path, seqs):
    p = tmp_path / "genome.fa"
    with open(p, "w") as fh:
        for sid, seq in seqs.items():
            fh.write(f">{sid}\n{seq}\n")
    return str(p)


def _gff3(tmp_path, name, seqids):
    p = tmp_path / name
    lines = ["##gff-version 3"]
    for i, sid in enumerate(seqids):
        # one plus-strand gene + one minus-strand gene per seqid (both strands)
        lines.append(f"{sid}\tHelixer\tgene\t100\t200\t.\t+\t.\tID=g{i}p")
        lines.append(f"{sid}\tHelixer\tgene\t300\t400\t.\t-\t.\tID=g{i}m")
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def _bam(tmp_path, name, sq):
    p = str(tmp_path / name)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    with pysam.AlignmentFile(p, "wb", header=header) as af:
        a = pysam.AlignedSegment(af.header)
        a.query_name = "r0"
        a.reference_id = 0
        a.reference_start = 10
        a.mapping_quality = 60
        a.cigartuples = [(0, 20)]
        a.query_sequence = "A" * 20
        af.write(a)
    pysam.index(p)
    return p


def _cfg(**kw):
    c = type("C", (), {})()
    c.genome_fasta = kw.get("genome_fasta")
    c.helixer_gff3 = kw.get("helixer_gff3")
    c.helixer_h5 = kw.get("helixer_h5")
    c.stringtie_list = kw.get("stringtie_list", [])
    c.bam_paths = kw.get("bam_paths", [])
    c.star_sj_paths = kw.get("star_sj_paths", [])
    return c


def test_all_concordant_passes(tmp_path):
    g = _genome(tmp_path, {"chr1": "ACGT" * 16, "chr2": "ACGT" * 16})
    gff = _gff3(tmp_path, "h.gff3", ["chr1", "chr2"])
    report = check_reference_concordance(_cfg(genome_fasta=g, helixer_gff3=gff))
    assert report.ok()
    assert not report.warnings


def test_disjoint_seqids_error(tmp_path):
    # genome names chr1/chr2; Helixer names 1/2 → the classic Chr1-vs-1 mismatch.
    g = _genome(tmp_path, {"chr1": "ACGT" * 16, "chr2": "ACGT" * 16})
    gff = _gff3(tmp_path, "h.gff3", ["1", "2"])
    report = check_reference_concordance(_cfg(genome_fasta=g, helixer_gff3=gff))
    assert not report.ok()
    assert any("NO seqid" in e for e in report.errors)


def test_partial_overlap_warns_with_names(tmp_path):
    g = _genome(tmp_path, {"chr1": "ACGT" * 16, "chr2": "ACGT" * 16})
    # Helixer has chr1 (shared) + chrX (absent from genome) → warning naming chrX.
    gff = _gff3(tmp_path, "h.gff3", ["chr1", "chrX"])
    report = check_reference_concordance(_cfg(genome_fasta=g, helixer_gff3=gff))
    assert report.ok()  # partial overlap is a warning, not an error
    assert report.warnings
    assert any("chrX" in w for w in report.warnings)


def test_length_mismatch_error(tmp_path):
    # genome chr1 length 64; a BAM declares chr1 length 999 → same name, diff
    # assembly. (chr1 is shared so it is not a disjoint error.)
    g = _genome(tmp_path, {"chr1": "ACGT" * 16})
    bam = _bam(tmp_path, "s.bam", [{"SN": "chr1", "LN": 999}])
    report = check_reference_concordance(_cfg(genome_fasta=g, bam_paths=[bam]))
    assert not report.ok()
    assert any("length disagrees" in e for e in report.errors)


def test_length_match_ok(tmp_path):
    g = _genome(tmp_path, {"chr1": "ACGT" * 16})  # 64 bp
    bam = _bam(tmp_path, "s.bam", [{"SN": "chr1", "LN": 64}])
    report = check_reference_concordance(_cfg(genome_fasta=g, bam_paths=[bam]))
    assert report.ok()


def test_alias_map_resolves_mismatch(tmp_path):
    g = _genome(tmp_path, {"chr1": "ACGT" * 16, "chr2": "ACGT" * 16})
    gff = _gff3(tmp_path, "h.gff3", ["1", "2"])
    alias = {"1": "chr1", "2": "chr2"}
    report = check_reference_concordance(
        _cfg(genome_fasta=g, helixer_gff3=gff), alias_map=alias)
    assert report.ok()
    assert not report.warnings


def test_alias_map_partial_still_warns(tmp_path):
    # alias resolves "1"→chr1 but "Z" has no alias and is absent → warning.
    g = _genome(tmp_path, {"chr1": "ACGT" * 16})
    gff = _gff3(tmp_path, "h.gff3", ["1", "Z"])
    report = check_reference_concordance(
        _cfg(genome_fasta=g, helixer_gff3=gff), alias_map={"1": "chr1"})
    assert report.ok()
    assert any("Z" in w for w in report.warnings)
