"""Tests for the v1-style multi-file evidence expansion (``io/fofn.py``).

Covers the three merging forms (repeatable + comma + ``-list`` FOFN), FOFN
blank/comment skipping, relative-path resolution against the list dir, missing
paths raising with a named source, cross-form de-duplication, and order.
"""

from __future__ import annotations

import os

import pytest

from helixforge.io.fofn import expand_file_args, expand_stringtie_args


def _touch(path) -> str:
    path.write_text("x")
    return str(path)


# ---------------------------------------------------------------------------
# expand_file_args — the three forms
# ---------------------------------------------------------------------------


def test_comma_split_of_single_value(tmp_path):
    a = _touch(tmp_path / "a.bam")
    b = _touch(tmp_path / "b.bam")
    assert expand_file_args((f"{a},{b}",), None, label="BAM file") == [a, b]


def test_repeated_values(tmp_path):
    a = _touch(tmp_path / "a.bam")
    b = _touch(tmp_path / "b.bam")
    assert expand_file_args((a, b), None, label="BAM file") == [a, b]


def test_merge_repeated_comma_and_list(tmp_path):
    a = _touch(tmp_path / "a.bam")
    b = _touch(tmp_path / "b.bam")
    c = _touch(tmp_path / "c.bam")
    d = _touch(tmp_path / "d.bam")
    fofn = tmp_path / "more.list"
    fofn.write_text(f"{c}\n{d}\n")
    # repeated `a`, comma `b`, then FOFN `c`,`d` — all merge, in order.
    out = expand_file_args((a, b), str(fofn), label="BAM file")
    assert out == [a, b, c, d]


def test_blank_and_comment_lines_skipped_in_list(tmp_path):
    a = _touch(tmp_path / "a.bam")
    b = _touch(tmp_path / "b.bam")
    fofn = tmp_path / "x.list"
    fofn.write_text(f"# header comment\n\n{a}\n   \n# {b} commented out\n{b}\n")
    assert expand_file_args((), str(fofn), label="BAM file") == [a, b]


def test_relative_paths_resolved_against_list_dir(tmp_path):
    sub = tmp_path / "data"
    sub.mkdir()
    _touch(sub / "a.bam")
    _touch(sub / "b.bam")
    fofn = sub / "bams.list"
    # bare names, no directory — must resolve against the list file's own dir.
    fofn.write_text("a.bam\nb.bam\n")
    out = expand_file_args((), str(fofn), label="BAM file")
    assert out == [str(sub / "a.bam"), str(sub / "b.bam")]


def test_missing_direct_path_raises_with_label(tmp_path):
    with pytest.raises(FileNotFoundError) as exc:
        expand_file_args((str(tmp_path / "nope.bam"),), None, label="BAM file")
    assert "BAM file not found" in str(exc.value)


def test_missing_list_path_raises_with_source_and_line(tmp_path):
    a = _touch(tmp_path / "a.bam")
    fofn = tmp_path / "x.list"
    fofn.write_text(f"{a}\n/no/such/file.bam\n")
    with pytest.raises(FileNotFoundError) as exc:
        expand_file_args((), str(fofn), label="BAM file")
    msg = str(exc.value)
    assert "BAM file not found" in msg
    assert "x.list" in msg
    assert "line 2" in msg


def test_dedup_across_forms_by_absolute_path(tmp_path):
    a = _touch(tmp_path / "a.bam")
    b = _touch(tmp_path / "b.bam")
    fofn = tmp_path / "x.list"
    fofn.write_text(f"{a}\n{b}\n")  # `a` repeats the direct flag
    out = expand_file_args((a,), str(fofn), label="BAM file")
    assert out == [a, b]  # `a` appears once


def test_order_preserved_direct_then_list(tmp_path):
    a = _touch(tmp_path / "a.bam")
    b = _touch(tmp_path / "b.bam")
    c = _touch(tmp_path / "c.bam")
    fofn = tmp_path / "x.list"
    fofn.write_text(f"{c}\n")
    out = expand_file_args((b, a), str(fofn), label="BAM file")
    assert out == [b, a, c]  # direct order kept, then FOFN entries


def test_empty_inputs_return_empty():
    assert expand_file_args((), None, label="BAM file") == []


# ---------------------------------------------------------------------------
# expand_stringtie_args — the --stringtie back-compat shim
# ---------------------------------------------------------------------------


def test_stringtie_repeatable_gtfs_no_shim(tmp_path):
    a = _touch(tmp_path / "a.gtf")
    b = _touch(tmp_path / "b.gtf")
    seen: list[str] = []
    out = expand_stringtie_args((a, b), None, warn=seen.append)
    assert out == [a, b]
    assert seen == []  # not a list-file → no deprecation warning


def test_stringtie_list_companion(tmp_path):
    a = _touch(tmp_path / "a.gtf")
    b = _touch(tmp_path / "b.gtf")
    fofn = tmp_path / "samples.list"
    fofn.write_text(f"{a}\n{b}\n")
    seen: list[str] = []
    out = expand_stringtie_args((), str(fofn), warn=seen.append)
    assert out == [a, b]
    assert seen == []  # the proper companion flag — never deprecated


def test_stringtie_legacy_listfile_shim_warns(tmp_path):
    a = _touch(tmp_path / "a.gtf")
    b = _touch(tmp_path / "b.gtf")
    legacy = tmp_path / "samples.txt"  # not a .gtf/.gff → looks like a FOFN
    legacy.write_text(f"# samples\n{a}\n{b}\n")
    seen: list[str] = []
    out = expand_stringtie_args((str(legacy),), None, warn=seen.append)
    assert out == [a, b]
    assert len(seen) == 1
    assert "--stringtie-list" in seen[0]


def test_stringtie_legacy_shim_emits_deprecationwarning_by_default(tmp_path):
    a = _touch(tmp_path / "a.gtf")
    legacy = tmp_path / "samples.txt"
    legacy.write_text(f"{a}\n")
    with pytest.warns(DeprecationWarning):
        out = expand_stringtie_args((str(legacy),), None)
    assert out == [a]


def test_stringtie_single_gtf_is_not_shimmed(tmp_path):
    a = _touch(tmp_path / "a.gtf")
    seen: list[str] = []
    out = expand_stringtie_args((a,), None, warn=seen.append)
    assert out == [a]
    assert seen == []  # a lone real .gtf is a GTF, never a list-file


def test_stringtie_gz_extension_not_shimmed(tmp_path):
    a = _touch(tmp_path / "a.gtf")
    gz = _touch(tmp_path / "sample.gtf.gz")  # .gz stripped → .gtf → real GTF
    seen: list[str] = []
    out = expand_stringtie_args((gz,), None, warn=seen.append)
    assert out == [gz]
    assert seen == []
    assert os.path.exists(a)  # fixture sanity
