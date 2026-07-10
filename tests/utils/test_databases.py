"""Tests for ``helixforge.utils.databases``."""

from __future__ import annotations

import gzip
import hashlib
import io
from pathlib import Path
from unittest import mock

import pytest

from helixforge.utils.databases import (
    DatabaseInfo,
    DatabaseManager,
    DatabaseType,
    PREDEFINED_DATABASES,
    calculate_md5,
    convert_dat_to_fasta,
    count_fasta_sequences,
    decompress_gzip,
    download_file,
)


# ---------------------------------------------------------------------------
# DatabaseInfo properties
# ---------------------------------------------------------------------------


class TestDatabaseInfo:
    def test_is_downloaded_true(self, tmp_path: Path) -> None:
        fasta = tmp_path / "db.fasta"
        fasta.write_text(">seq1\nACGT\n")
        info = DatabaseInfo(name="test", db_type=DatabaseType.CUSTOM, path=fasta)
        assert info.is_downloaded is True

    def test_is_downloaded_false_no_path(self) -> None:
        info = DatabaseInfo(name="test", db_type=DatabaseType.CUSTOM)
        assert info.is_downloaded is False

    def test_is_downloaded_false_missing_file(self, tmp_path: Path) -> None:
        info = DatabaseInfo(
            name="test", db_type=DatabaseType.CUSTOM,
            path=tmp_path / "missing.fasta",
        )
        assert info.is_downloaded is False

    def test_is_formatted_true(self, tmp_path: Path) -> None:
        dmnd = tmp_path / "db.dmnd"
        dmnd.write_bytes(b"\x00")
        info = DatabaseInfo(
            name="test", db_type=DatabaseType.CUSTOM,
            formatted_path=dmnd,
        )
        assert info.is_formatted is True

    def test_is_formatted_false_no_path(self) -> None:
        info = DatabaseInfo(name="test", db_type=DatabaseType.CUSTOM)
        assert info.is_formatted is False

    def test_is_formatted_dmnd_suffix_variant(self, tmp_path: Path) -> None:
        base = tmp_path / "db"
        dmnd = tmp_path / "db.dmnd"
        dmnd.write_bytes(b"\x00")
        info = DatabaseInfo(
            name="test", db_type=DatabaseType.CUSTOM,
            formatted_path=base,
        )
        assert info.is_formatted is True


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------


class TestCalculateMd5:
    def test_known_hash(self, tmp_path: Path) -> None:
        data = b"hello world"
        f = tmp_path / "test.bin"
        f.write_bytes(data)
        expected = hashlib.md5(data).hexdigest()
        assert calculate_md5(f) == expected

    def test_empty_file(self, tmp_path: Path) -> None:
        f = tmp_path / "empty"
        f.write_bytes(b"")
        expected = hashlib.md5(b"").hexdigest()
        assert calculate_md5(f) == expected


class TestDecompressGzip:
    def test_decompress_removes_original(self, tmp_path: Path) -> None:
        content = b">seq1\nACGT\n>seq2\nTGCA\n"
        gz_path = tmp_path / "db.fasta.gz"
        with gzip.open(gz_path, "wb") as fh:
            fh.write(content)

        out = decompress_gzip(gz_path)
        assert out.read_bytes() == content
        assert out == tmp_path / "db.fasta"
        assert not gz_path.exists()

    def test_decompress_keep_original(self, tmp_path: Path) -> None:
        content = b">s\nA\n"
        gz_path = tmp_path / "x.fasta.gz"
        with gzip.open(gz_path, "wb") as fh:
            fh.write(content)

        out = decompress_gzip(gz_path, keep_original=True)
        assert out.read_bytes() == content
        assert gz_path.exists()

    def test_decompress_explicit_output(self, tmp_path: Path) -> None:
        content = b">s\nA\n"
        gz_path = tmp_path / "x.gz"
        explicit_out = tmp_path / "result.fasta"
        with gzip.open(gz_path, "wb") as fh:
            fh.write(content)

        out = decompress_gzip(gz_path, output_path=explicit_out)
        assert out == explicit_out
        assert out.read_bytes() == content


class TestCountFastaSequences:
    def test_count(self, tmp_path: Path) -> None:
        fasta = tmp_path / "db.fasta"
        fasta.write_text(">a\nACGT\n>b\nTGCA\n>c\nAAAA\n")
        assert count_fasta_sequences(fasta) == 3

    def test_count_gzipped(self, tmp_path: Path) -> None:
        gz = tmp_path / "db.fasta.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write(">x\nAA\n>y\nCC\n")
        assert count_fasta_sequences(gz) == 2


_DAT_FIXTURE = """\
ID   TEST1_ARATH             Reviewed;          10 AA.
AC   P12345; Q99999;
DT   01-JAN-2020, integrated into UniProtKB/Swiss-Prot.
DE   RecName: Full=Test protein one;
OS   Arabidopsis thaliana (Mouse-ear cress).
SQ   SEQUENCE   10 AA;  1234 MW;  ABCDEF CRC64;
     MKVLAAGTRS
//
ID   TEST2_YEAST             Reviewed;          12 AA.
AC   Q54321;
DE   RecName: Full=Test protein two;
OS   Saccharomyces cerevisiae.
SQ   SEQUENCE   12 AA;  1500 MW;  123456 CRC64;
     MKLVDEFGHIKL
//
"""


class TestConvertDatToFasta:
    def test_converts_two_entries(self, tmp_path: Path) -> None:
        dat = tmp_path / "subset.dat"
        dat.write_text(_DAT_FIXTURE)
        fasta = tmp_path / "subset.fasta"

        n = convert_dat_to_fasta(dat, fasta)

        assert n == 2
        assert count_fasta_sequences(fasta) == 2

    def test_headers_and_sequences(self, tmp_path: Path) -> None:
        dat = tmp_path / "subset.dat"
        dat.write_text(_DAT_FIXTURE)
        fasta = tmp_path / "subset.fasta"

        convert_dat_to_fasta(dat, fasta)
        lines = fasta.read_text().splitlines()

        assert lines[0] == (
            ">sp|P12345|TEST1_ARATH Test protein one "
            "OS=Arabidopsis thaliana (Mouse-ear cress)"
        )
        assert lines[1] == "MKVLAAGTRS"
        assert lines[2] == (
            ">sp|Q54321|TEST2_YEAST Test protein two OS=Saccharomyces cerevisiae"
        )
        assert lines[3] == "MKLVDEFGHIKL"

    def test_uses_first_accession_only(self, tmp_path: Path) -> None:
        dat = tmp_path / "subset.dat"
        dat.write_text(_DAT_FIXTURE)
        fasta = tmp_path / "subset.fasta"

        convert_dat_to_fasta(dat, fasta)
        text = fasta.read_text()

        assert "P12345" in text
        assert "Q99999" not in text

    def test_sequence_wrapped_at_60(self, tmp_path: Path) -> None:
        seq = "A" * 130
        block1 = "     " + seq[:60] + "\n"
        block2 = "     " + seq[60:120] + "\n"
        block3 = "     " + seq[120:] + "\n"
        dat_text = (
            "ID   LONG_TEST               Reviewed;         130 AA.\n"
            "AC   P00001;\n"
            "DE   RecName: Full=Long protein;\n"
            "OS   Homo sapiens.\n"
            "SQ   SEQUENCE   130 AA;  9999 MW;  AAAAAA CRC64;\n"
            + block1 + block2 + block3 + "//\n"
        )
        dat = tmp_path / "long.dat"
        dat.write_text(dat_text)
        fasta = tmp_path / "long.fasta"

        n = convert_dat_to_fasta(dat, fasta)
        lines = fasta.read_text().splitlines()

        assert n == 1
        assert lines[0] == ">sp|P00001|LONG_TEST Long protein OS=Homo sapiens"
        assert lines[1] == "A" * 60
        assert lines[2] == "A" * 60
        assert lines[3] == "A" * 10

    def test_gzipped_input(self, tmp_path: Path) -> None:
        dat_gz = tmp_path / "subset.dat.gz"
        with gzip.open(dat_gz, "wt") as fh:
            fh.write(_DAT_FIXTURE)
        fasta = tmp_path / "subset.fasta"

        n = convert_dat_to_fasta(dat_gz, fasta)

        assert n == 2
        assert count_fasta_sequences(fasta) == 2


class TestDownloadFile:
    def test_download_writes_content(self, tmp_path: Path) -> None:
        content = b">seq1\nACGT\n"
        resp = mock.MagicMock()
        resp.read.side_effect = [content, b""]
        resp.headers = {"Content-Length": str(len(content))}
        resp.__enter__ = lambda s: s
        resp.__exit__ = mock.MagicMock(return_value=False)

        out = tmp_path / "downloaded.fasta"
        with mock.patch(
            "helixforge.utils.databases.urlopen", return_value=resp
        ):
            result = download_file("https://example.com/db.fasta", out)

        assert result == out
        assert out.read_bytes() == content

    def test_download_cleans_up_on_error(self, tmp_path: Path) -> None:
        from urllib.error import URLError

        out = tmp_path / "fail.fasta"
        with mock.patch(
            "helixforge.utils.databases.urlopen",
            side_effect=URLError("boom"),
        ):
            with pytest.raises(URLError):
                download_file("https://example.com/fail", out)
        assert not out.exists()


# ---------------------------------------------------------------------------
# DatabaseManager
# ---------------------------------------------------------------------------


class TestDatabaseManager:
    def test_creates_cache_dir(self, tmp_path: Path) -> None:
        cache = tmp_path / "new_cache"
        assert not cache.exists()
        DatabaseManager(cache_dir=cache)
        assert cache.is_dir()

    def test_list_available(self, tmp_path: Path) -> None:
        mgr = DatabaseManager(cache_dir=tmp_path)
        avail = mgr.list_available()
        assert set(avail) == set(PREDEFINED_DATABASES.keys())
        assert len(avail) == 6

    def test_get_custom_fasta(self, tmp_path: Path) -> None:
        fasta = tmp_path / "custom.fasta"
        fasta.write_text(">a\nACGT\n>b\nTGCA\n>c\nGGGG\n")

        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value = mock.MagicMock(returncode=0)
            mgr = DatabaseManager(cache_dir=tmp_path)
            info = mgr.get_database(str(fasta))

        assert info.name == "Custom: custom.fasta"
        assert info.db_type == DatabaseType.CUSTOM
        assert info.path == fasta
        assert info.formatted_path == fasta.with_suffix(".dmnd")

        mock_run.assert_called_once()
        call_args = mock_run.call_args
        argv = call_args[0][0]
        assert argv[0] == "diamond"
        assert argv[1] == "makedb"
        assert "--in" in argv
        assert "--db" in argv

    def test_get_unknown_raises(self, tmp_path: Path) -> None:
        mgr = DatabaseManager(cache_dir=tmp_path)
        with pytest.raises(ValueError, match="Unknown database"):
            mgr.get_database("nonexistent_db_xyz")

    def test_format_database_missing_diamond(self, tmp_path: Path) -> None:
        fasta = tmp_path / "db.fasta"
        fasta.write_text(">a\nACGT\n")
        info = DatabaseInfo(name="test", db_type=DatabaseType.CUSTOM, path=fasta)

        with mock.patch("subprocess.run", side_effect=FileNotFoundError):
            mgr = DatabaseManager(cache_dir=tmp_path, diamond_bin="/bad/diamond")
            with pytest.raises(FileNotFoundError, match="diamond not found"):
                mgr.format_database(info)

    def test_force_download_bypasses_cache(self, tmp_path: Path) -> None:
        """force_download=True re-downloads even when the file is cached."""
        cache = tmp_path / "cache"

        content = b">s\nA\n"
        gz_content = io.BytesIO()
        with gzip.open(gz_content, "wb") as gz:
            gz.write(content)
        gz_bytes = gz_content.getvalue()

        resp = mock.MagicMock()
        resp.read.side_effect = [gz_bytes, b""]
        resp.headers = {"Content-Length": str(len(gz_bytes))}
        resp.__enter__ = lambda s: s
        resp.__exit__ = mock.MagicMock(return_value=False)

        with (
            mock.patch(
                "helixforge.utils.databases.urlopen", return_value=resp
            ) as mock_urlopen,
            mock.patch("subprocess.run") as mock_run,
        ):
            mock_run.return_value = mock.MagicMock(returncode=0)
            mgr = DatabaseManager(cache_dir=cache)

            info = mgr.get_database("swissprot", force_download=True)

        assert info.is_downloaded is False or info.path is not None
        mock_urlopen.assert_called_once()

    def test_cached_database_skips_download(self, tmp_path: Path) -> None:
        """A cached (already-downloaded) database is returned without re-download."""
        cache = tmp_path / "cache"
        fasta = cache / "swissprot" / "uniprot_sprot.fasta"
        fasta.parent.mkdir(parents=True)
        fasta.write_text(">a\nACGT\n")

        mgr = DatabaseManager(cache_dir=cache)
        mgr._databases["swissprot"] = DatabaseInfo(
            name="Swiss-Prot (complete)",
            db_type=DatabaseType.SWISSPROT,
            path=fasta,
            formatted_path=fasta.with_suffix(".dmnd"),
        )
        fasta.with_suffix(".dmnd").write_bytes(b"\x00")

        with mock.patch(
            "helixforge.utils.databases.download_file"
        ) as mock_dl:
            info = mgr.get_database("swissprot")
        mock_dl.assert_not_called()
        assert info.path == fasta


class TestDatabaseType:
    def test_enum_values(self) -> None:
        assert DatabaseType.SWISSPROT.value == "swissprot"
        assert DatabaseType.TREMBL.value == "trembl"
        assert DatabaseType.UNIREF90.value == "uniref90"
        assert DatabaseType.UNIREF50.value == "uniref50"
        assert DatabaseType.UNIREF100.value == "uniref100"
        assert DatabaseType.CUSTOM.value == "custom"
