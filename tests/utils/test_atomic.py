"""Tests for utils/atomic.py (Phase 22 D1). Floor: 3.

The atomic-write guarantee: a crash mid-write never corrupts the destination —
a reader sees either the complete old file or the complete new one.
"""

import os
import stat

import pytest

from helixforge.utils.atomic import atomic_write


def test_temp_then_rename_leaves_no_partial_file(tmp_path):
    target = tmp_path / "out.txt"
    with atomic_write(target) as fh:
        fh.write("hello world\n")
    assert target.read_text() == "hello world\n"
    # No leftover temp files beside the destination.
    leftovers = [p.name for p in tmp_path.iterdir() if p.name != "out.txt"]
    assert leftovers == []


def test_crash_mid_write_keeps_old_file(tmp_path):
    target = tmp_path / "id_map.json"
    target.write_text('{"old": 1}')

    class Boom(RuntimeError):
        pass

    with pytest.raises(Boom):
        with atomic_write(target) as fh:
            fh.write('{"new": 2, "partial":')  # half-written
            raise Boom("kill -9 mid-write")

    # The prior file is intact; the partial write never reached the target.
    assert target.read_text() == '{"old": 1}'
    leftovers = [p.name for p in tmp_path.iterdir() if p.name != "id_map.json"]
    assert leftovers == []


def test_overwrite_preserves_existing_permissions(tmp_path):
    target = tmp_path / "report.tsv"
    target.write_text("v1\n")
    os.chmod(target, 0o640)
    with atomic_write(target) as fh:
        fh.write("v2\n")
    assert target.read_text() == "v2\n"
    assert stat.S_IMODE(target.stat().st_mode) == 0o640


def test_new_file_is_not_private_only(tmp_path):
    # A fresh file must respect the umask default, not mkstemp's 0o600.
    target = tmp_path / "fresh.txt"
    with atomic_write(target) as fh:
        fh.write("x")
    mode = stat.S_IMODE(target.stat().st_mode)
    cur = os.umask(0)
    os.umask(cur)
    assert mode == (0o666 & ~cur)


def test_binary_mode_writes_bytes(tmp_path):
    target = tmp_path / "blob.bin"
    payload = b"\x00\x01\x02\xff"
    with atomic_write(target, mode="wb") as fh:
        fh.write(payload)
    assert target.read_bytes() == payload


def test_creates_missing_parent_directory(tmp_path):
    target = tmp_path / "nested" / "deep" / "out.txt"
    with atomic_write(target) as fh:
        fh.write("ok\n")
    assert target.read_text() == "ok\n"
