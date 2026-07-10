"""Tests for region parsing + coordinate conversion (Phase 0). Floor: 15."""

import pytest

from helixforge.utils.regions import (
    format_region,
    gff3_to_internal,
    internal_to_gff3,
    parse_region,
)


def test_parse_region_basic():
    assert parse_region("chr1:1000-2000") == ("chr1", 1000, 2000)


def test_parse_region_seqid_only():
    assert parse_region("chr1") == ("chr1", None, None)


def test_parse_region_seqid_with_colon_in_name():
    # rpartition keeps everything before the last colon as seqid
    assert parse_region("scaffold:1:100-200") == ("scaffold:1", 100, 200)


def test_parse_region_single_base():
    assert parse_region("chr1:500-500") == ("chr1", 500, 500)


def test_parse_region_rejects_inverted():
    with pytest.raises(ValueError):
        parse_region("chr1:2000-1000")


def test_parse_region_rejects_non_integer():
    with pytest.raises(ValueError):
        parse_region("chr1:abc-200")


def test_parse_region_rejects_missing_dash():
    with pytest.raises(ValueError):
        parse_region("chr1:1000")


def test_parse_region_rejects_empty():
    with pytest.raises(ValueError):
        parse_region("")


def test_format_region_basic():
    assert format_region("chr1", 1000, 2000) == "chr1:1000-2000"


def test_format_region_seqid_only():
    assert format_region("chr1") == "chr1"


def test_format_region_rejects_partial():
    with pytest.raises(ValueError):
        format_region("chr1", 1000)


def test_format_region_rejects_inverted():
    with pytest.raises(ValueError):
        format_region("chr1", 2000, 1000)


def test_parse_format_roundtrip():
    assert format_region(*parse_region("chr5:123-456")) == "chr5:123-456"


def test_gff3_to_internal():
    assert gff3_to_internal(1000, 2000) == (999, 2000)


def test_gff3_to_internal_first_base():
    assert gff3_to_internal(1, 3) == (0, 3)


def test_gff3_to_internal_rejects_zero_start():
    with pytest.raises(ValueError):
        gff3_to_internal(0, 100)


def test_gff3_to_internal_rejects_end_before_start():
    with pytest.raises(ValueError):
        gff3_to_internal(100, 50)


def test_internal_to_gff3():
    assert internal_to_gff3(999, 2000) == (1000, 2000)


def test_internal_to_gff3_first_base():
    assert internal_to_gff3(0, 3) == (1, 3)


def test_internal_to_gff3_rejects_negative():
    with pytest.raises(ValueError):
        internal_to_gff3(-1, 100)


def test_internal_to_gff3_rejects_zero_width():
    with pytest.raises(ValueError):
        internal_to_gff3(100, 100)


def test_conversion_roundtrip():
    assert internal_to_gff3(*gff3_to_internal(1000, 2000)) == (1000, 2000)
