"""Tests for GenomeAccessor (Phase 1). Floor: 17.

Both-strand extraction with literal expected strings; 0-based behavior; bounds.
"""

import pytest

from helixforge.io.fasta import GenomeAccessor


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        GenomeAccessor("/no/such/genome.fa")


def test_get_seqids_sorted(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_seqids() == ["chr1", "chr2"]


def test_get_length(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_length("chr1") == 64
        assert g.get_length("chr2") == 64


def test_get_length_unknown_seqid(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        with pytest.raises(KeyError):
            g.get_length("chrX")


def test_get_scaffold_lengths(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_scaffold_lengths() == {"chr1": 64, "chr2": 64}


def test_contains_true(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert "chr1" in g


def test_contains_false(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert "chrX" not in g


def test_get_sequence_plus_first_block(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr1", 0, 4, "+") == "AAAA"


def test_get_sequence_plus_second_block(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr1", 4, 8, "+") == "CCCC"


def test_get_sequence_plus_third_fourth_block(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr1", 8, 12, "+") == "GGGG"
        assert g.get_sequence("chr1", 12, 16, "+") == "TTTT"


def test_get_sequence_zero_based_first_base(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr1", 0, 1, "+") == "A"


def test_get_sequence_minus_single_block(fasta_path):
    # RC of "AAAA" is "TTTT"
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr1", 0, 4, "-") == "TTTT"


def test_get_sequence_minus_two_blocks(fasta_path):
    # RC of "AAAACCCC" is "GGGGTTTT"
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr1", 0, 8, "-") == "GGGGTTTT"


def test_get_sequence_uppercases_lowercase_fasta(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr2", 0, 4, "+") == "ACGT"


def test_get_sequence_minus_lowercase(fasta_path):
    # chr2 "acgt" -> upper "ACGT" -> RC "ACGT"
    with GenomeAccessor(fasta_path) as g:
        assert g.get_sequence("chr2", 0, 4, "-") == "ACGT"


def test_get_sequence_out_of_bounds_raises_index_error(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        with pytest.raises(IndexError):
            g.get_sequence("chr1", 0, 65, "+")


def test_get_sequence_negative_start_raises_value_error(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        with pytest.raises(ValueError):
            g.get_sequence("chr1", -1, 10, "+")


def test_get_sequence_end_le_start_raises_value_error(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        with pytest.raises(ValueError):
            g.get_sequence("chr1", 20, 20, "+")


def test_get_sequence_bad_strand_raises_value_error(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        with pytest.raises(ValueError):
            g.get_sequence("chr1", 0, 4, ".")


def test_get_sequence_unknown_seqid_raises_key_error(fasta_path):
    with GenomeAccessor(fasta_path) as g:
        with pytest.raises(KeyError):
            g.get_sequence("chrX", 0, 4, "+")


def test_context_manager_closes(fasta_path):
    g = GenomeAccessor(fasta_path)
    g.close()
    # idempotent double-close should not raise
    g.close()
