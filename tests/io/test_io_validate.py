"""Format-sniffer tests (biology.v4 Phase 26 D2). Floor for this file: ≥5.

Each sniffer accepts a good file and raises a line-located :class:`FormatError`
on a malformed one (truncated SJ, bad strand code, non-``>`` FASTA, missing
``@SQ``). Synthetic fixtures only (CLAUDE.md §12).
"""

import pytest

from helixforge.io.validate import (
    FormatError,
    sniff_fasta,
    sniff_gff3,
    sniff_gtf,
    sniff_star_sj,
)

# --- FASTA --------------------------------------------------------------


def test_sniff_fasta_good(tmp_path):
    p = tmp_path / "g.fa"
    p.write_text(">chr1\nACGTACGTAC\n>chr2\nTTTTGGGGCC\n")
    assert sniff_fasta(p) is None


def test_sniff_fasta_non_header_start(tmp_path):
    p = tmp_path / "bad.fa"
    p.write_text("ACGTACGT\n>chr1\nACGT\n")
    with pytest.raises(FormatError) as ei:
        sniff_fasta(p)
    assert ei.value.line == 1
    assert ">" in str(ei.value)


def test_sniff_fasta_empty_record(tmp_path):
    p = tmp_path / "empty_rec.fa"
    p.write_text(">chr1\n>chr2\nACGT\n")
    with pytest.raises(FormatError) as ei:
        sniff_fasta(p)
    assert ei.value.line == 1  # the empty header
    assert "empty FASTA record" in str(ei.value)


def test_sniff_fasta_bare_cr(tmp_path):
    p = tmp_path / "cr.fa"
    p.write_bytes(b">chr1\rACGTACGT\r")  # bare CR line endings
    with pytest.raises(FormatError) as ei:
        sniff_fasta(p)
    assert "carriage return" in str(ei.value)


# --- GFF3 / GTF ---------------------------------------------------------


def test_sniff_gff3_good(tmp_path):
    p = tmp_path / "a.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tID=g1\n"
        "chr2\tsrc\tCDS\t1050\t1200\t.\t-\t0\tID=c1\n"
    )
    assert sniff_gff3(p) is None


def test_sniff_gff3_missing_version(tmp_path):
    p = tmp_path / "nover.gff3"
    p.write_text("chr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tID=g1\n")
    with pytest.raises(FormatError) as ei:
        sniff_gff3(p)
    assert "gff-version" in str(ei.value)


def test_sniff_gff3_bad_strand(tmp_path):
    p = tmp_path / "badstrand.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t1000\t2000\t.\tX\t.\tID=g1\n"
    )
    with pytest.raises(FormatError) as ei:
        sniff_gff3(p)
    assert ei.value.line == 2
    assert ei.value.value == "X"


def test_sniff_gff3_too_few_fields(tmp_path):
    p = tmp_path / "short.gff3"
    p.write_text("##gff-version 3\nchr1\tsrc\tgene\t1000\t2000\n")
    with pytest.raises(FormatError) as ei:
        sniff_gff3(p)
    assert ei.value.line == 2
    assert "tab fields" in str(ei.value)


def test_sniff_gtf_good_no_version_required(tmp_path):
    p = tmp_path / "a.gtf"
    # Minus-strand transcript — both strands exercised across this file.
    p.write_text(
        'chr3\tStringTie\texon\t3000\t3200\t.\t-\t.\ttranscript_id "t1";\n'
        'chr3\tStringTie\texon\t3300\t3500\t.\t-\t.\ttranscript_id "t1";\n'
    )
    assert sniff_gtf(p) is None


def test_sniff_gtf_bad_phase(tmp_path):
    p = tmp_path / "badphase.gtf"
    p.write_text('chr1\tsrc\tCDS\t1000\t2000\t.\t+\t5\ttranscript_id "t1";\n')
    with pytest.raises(FormatError) as ei:
        sniff_gtf(p)
    assert ei.value.value == "5"


# --- STAR SJ.out.tab ----------------------------------------------------


def _sj_line(chrom="chr1", start=1001, end=1200, strand=1, motif=1,
             annot=0, uniq=10, multi=0, overhang=30):
    return f"{chrom}\t{start}\t{end}\t{strand}\t{motif}\t{annot}\t{uniq}\t{multi}\t{overhang}\n"


def test_sniff_star_sj_good(tmp_path):
    p = tmp_path / "SJ.out.tab"
    p.write_text(_sj_line(strand=1) + _sj_line(chrom="chr2", strand=2))
    assert sniff_star_sj(p) is None


def test_sniff_star_sj_truncated(tmp_path):
    p = tmp_path / "trunc.tab"
    p.write_text("chr1\t1001\t1200\t1\t1\t0\n")  # only 6 columns
    with pytest.raises(FormatError) as ei:
        sniff_star_sj(p)
    assert ei.value.line == 1
    assert "columns" in str(ei.value)


def test_sniff_star_sj_bad_strand_code(tmp_path):
    p = tmp_path / "badstrand.tab"
    p.write_text(_sj_line(strand=3))  # strand code must be 0/1/2
    with pytest.raises(FormatError) as ei:
        sniff_star_sj(p)
    assert ei.value.value == "3"
    assert "strand code" in str(ei.value)


def test_sniff_star_sj_non_numeric(tmp_path):
    p = tmp_path / "header.tab"
    # A header line accidentally left in — non-numeric coordinate.
    p.write_text("chrom\tstart\tend\tstrand\tmotif\tannot\tuniq\tmulti\toverhang\n")
    with pytest.raises(FormatError) as ei:
        sniff_star_sj(p)
    assert "non-numeric" in str(ei.value)


# --- BAM ----------------------------------------------------------------


def _write_bam(path, header):
    import pysam

    with pysam.AlignmentFile(str(path), "wb", header=header) as af:
        if header.get("SQ"):
            a = pysam.AlignedSegment(af.header)
            a.query_name = "r0"
            a.reference_id = 0
            a.reference_start = 10
            a.mapping_quality = 60
            a.cigartuples = [(0, 20)]
            a.query_sequence = "A" * 20
            af.write(a)


def test_sniff_bam_good(tmp_path):
    from helixforge.io.validate import sniff_bam

    p = tmp_path / "good.bam"
    _write_bam(p, {"HD": {"VN": "1.6", "SO": "coordinate"},
                   "SQ": [{"SN": "chr1", "LN": 1000}]})
    assert sniff_bam(p) is None


def test_sniff_bam_missing_sq(tmp_path):
    from helixforge.io.validate import sniff_bam

    p = tmp_path / "nosq.bam"
    _write_bam(p, {"HD": {"VN": "1.6", "SO": "coordinate"}})  # no @SQ
    with pytest.raises(FormatError) as ei:
        sniff_bam(p)
    assert "@SQ" in str(ei.value)


def test_sniff_bam_not_coordinate_sorted(tmp_path):
    from helixforge.io.validate import sniff_bam

    p = tmp_path / "qsort.bam"
    _write_bam(p, {"HD": {"VN": "1.6", "SO": "queryname"},
                   "SQ": [{"SN": "chr1", "LN": 1000}]})
    with pytest.raises(FormatError) as ei:
        sniff_bam(p)
    assert "coordinate" in str(ei.value)
