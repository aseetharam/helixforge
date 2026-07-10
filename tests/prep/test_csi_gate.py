"""CSI-gate tests (biology.v4 Phase 26 D3). Floor for this file: ≥2.

A reference contig > 2^29 bp needs a ``.csi`` BAM index; ``.bai`` cannot address
those coordinates and pysam fetch silently returns nothing (assessment §3.4).
To avoid materialising a 512 Mb fixture, the 2^29 threshold is monkeypatched to
a tiny value so a normal small contig counts as "large" — the logic under test
is identical. Both the open-path enforcement (`JunctionExtractor` /
`CoverageCalculator`) and the preflight check are covered.
"""

import pysam
import pytest

from helixforge.io import bam as bam_mod
from helixforge.io.bam import CoverageCalculator, JunctionExtractor
from helixforge.io.validate import FormatError
from helixforge.prep import preflight as pf

_HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000}],  # plus & minus reads both at small coords
}


def _build_indexed_bam(path):
    with pysam.AlignmentFile(path, "wb", header=_HEADER) as af:
        for i, (start, strand) in enumerate([(100, "+"), (400, "-")]):
            a = pysam.AlignedSegment(af.header)
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, 50)]
            a.query_sequence = "A" * 50
            a.set_tag("XS", strand, "A")
            af.write(a)
    pysam.index(path)  # creates a .bai


@pytest.fixture
def bai_only_bam(tmp_path):
    p = str(tmp_path / "big.bam")
    _build_indexed_bam(p)
    return p


# --- open-path enforcement (io/bam.py) ----------------------------------


def test_junction_extractor_requires_csi_for_large_contig(monkeypatch, bai_only_bam):
    monkeypatch.setattr(bam_mod, "CSI_REQUIRED_THRESHOLD", 100)
    with pytest.raises(FormatError) as ei:
        JunctionExtractor(bai_only_bam)
    assert ".csi" in str(ei.value)


def test_junction_extractor_accepts_csi(monkeypatch, bai_only_bam):
    monkeypatch.setattr(bam_mod, "CSI_REQUIRED_THRESHOLD", 100)
    pysam.index("-c", bai_only_bam)  # build a real .csi sibling
    with JunctionExtractor(bai_only_bam) as je:
        # opens cleanly; junctions list is empty (no spliced reads) but no raise
        assert je.extract_junctions("chr1", 0, 1000) == []


def test_coverage_calculator_requires_csi(monkeypatch, bai_only_bam):
    monkeypatch.setattr(bam_mod, "CSI_REQUIRED_THRESHOLD", 100)
    with pytest.raises(FormatError):
        CoverageCalculator.from_bam(bai_only_bam)


def test_small_contig_with_bai_is_fine(bai_only_bam):
    # Default 2^29 threshold: a 1000 bp contig is small → .bai is sufficient.
    with JunctionExtractor(bai_only_bam) as je:
        assert je.extract_junctions("chr1", 0, 1000) == []


# --- preflight check (prep/preflight.py) --------------------------------


def _genome(tmp_path):
    p = tmp_path / "g.fa"
    p.write_text(">chr1\n" + "ACGT" * 16 + "\n")  # 64 bp
    return str(p)


def test_check_csi_requirement_errors_with_only_bai(monkeypatch, tmp_path, bai_only_bam):
    monkeypatch.setattr(pf, "CSI_REQUIRED_THRESHOLD", 10)
    cfg = type("C", (), {})()
    cfg.genome_fasta = _genome(tmp_path)
    cfg.bam_paths = [bai_only_bam]
    errors = pf.check_csi_requirement(cfg)
    assert len(errors) == 1
    assert ".csi" in errors[0]


def test_check_csi_requirement_ok_with_csi(monkeypatch, tmp_path, bai_only_bam):
    monkeypatch.setattr(pf, "CSI_REQUIRED_THRESHOLD", 10)
    open(bai_only_bam + ".csi", "wb").close()
    cfg = type("C", (), {})()
    cfg.genome_fasta = _genome(tmp_path)
    cfg.bam_paths = [bai_only_bam]
    assert pf.check_csi_requirement(cfg) == []


def test_check_csi_requirement_no_large_contig_ok(tmp_path, bai_only_bam):
    # Default threshold (2^29): no large contig → no CSI required.
    cfg = type("C", (), {})()
    cfg.genome_fasta = _genome(tmp_path)
    cfg.bam_paths = [bai_only_bam]
    assert pf.check_csi_requirement(cfg) == []
