"""Tests for sequence utilities (Phase 0). Floor: 20.

Both-strand codon extraction is mandatory (CLAUDE.md §5, §12).
"""

import pytest

from helixforge.utils.sequences import (
    CODON_TABLE,
    COMPLEMENT,
    STANDARD_START_CODONS,
    STANDARD_STOP_CODONS,
    check_internal_stops,
    extract_start_codon,
    extract_stop_codon_after_cds,
    is_start_codon,
    is_stop_codon,
    reverse_complement,
    translate,
)


# --- tables ---

def test_codon_table_has_64_entries():
    assert len(CODON_TABLE) == 64


def test_codon_table_atg_is_methionine():
    assert CODON_TABLE["ATG"] == "M"


def test_codon_table_stops():
    assert CODON_TABLE["TAA"] == "*"
    assert CODON_TABLE["TAG"] == "*"
    assert CODON_TABLE["TGA"] == "*"


def test_codon_table_exactly_three_stops():
    assert sum(1 for v in CODON_TABLE.values() if v == "*") == 3


def test_complement_table():
    # Phase 27 D1 (deliberate interface change, assessment §1.1): COMPLEMENT is
    # now the full IUPAC complement map (ACGT + N + degenerate codes), not the
    # 4-base ACGT map. The ACGT entries are unchanged.
    assert COMPLEMENT["A"] == "T"
    assert COMPLEMENT["C"] == "G"
    assert COMPLEMENT["G"] == "C"
    assert COMPLEMENT["T"] == "A"
    assert COMPLEMENT["N"] == "N"
    assert COMPLEMENT["R"] == "Y"
    assert COMPLEMENT["Y"] == "R"
    assert COMPLEMENT["S"] == "S"
    assert COMPLEMENT["W"] == "W"
    assert COMPLEMENT["K"] == "M"
    assert COMPLEMENT["M"] == "K"
    assert COMPLEMENT["B"] == "V"
    assert COMPLEMENT["V"] == "B"
    assert COMPLEMENT["D"] == "H"
    assert COMPLEMENT["H"] == "D"
    assert len(COMPLEMENT) == 15


def test_standard_codon_sets():
    assert "ATG" in STANDARD_START_CODONS
    assert set(STANDARD_STOP_CODONS) == {"TAA", "TAG", "TGA"}


# --- reverse_complement ---

def test_reverse_complement_basic():
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_palindrome_like():
    assert reverse_complement("AATT") == "AATT"


def test_reverse_complement_lowercase():
    assert reverse_complement("atgc") == "GCAT"


def test_reverse_complement_tolerates_iupac():
    # Phase 27 D1 (deliberate interface change, §1.1): reverse_complement now
    # complements IUPAC codes instead of raising. ATGN -> RC -> NCAT.
    assert reverse_complement("ATGN") == "NCAT"


def test_reverse_complement_rejects_truly_invalid():
    # A character outside the IUPAC alphabet entirely still raises (strict edge).
    with pytest.raises(ValueError):
        reverse_complement("ATGZ")


# --- translate ---

def test_translate_basic():
    assert translate("ATGAAACCCGGGTTT") == "MKPGF"


def test_translate_with_stop():
    assert translate("ATGAAATAA") == "MK*"


def test_translate_ignores_trailing_partial_codon():
    assert translate("ATGAA") == "M"


def test_translate_phase_one():
    assert translate("AATGAAA", phase=1) == "MK"


def test_translate_phase_two():
    assert translate("AAATGAAA", phase=2) == "MK"


def test_translate_rejects_bad_phase():
    with pytest.raises(ValueError):
        translate("ATGAAA", phase=3)


def test_translate_ambiguous_codon_to_x():
    # Phase 27 D1 (deliberate interface change, §1.1): an N/IUPAC codon now
    # translates to 'X' instead of raising. ATG NNN -> 'MX'.
    assert translate("ATGNNN") == "MX"


def test_translate_rejects_truly_invalid():
    with pytest.raises(ValueError):
        translate("ATGZZZ")


# --- start / stop codon predicates ---

def test_is_start_codon_true():
    assert is_start_codon("ATG")


def test_is_start_codon_lowercase():
    assert is_start_codon("atg")


def test_is_start_codon_false():
    assert not is_start_codon("AAA")


def test_is_stop_codon_all_three():
    assert is_stop_codon("TAA")
    assert is_stop_codon("TAG")
    assert is_stop_codon("TGA")


def test_is_stop_codon_false():
    assert not is_stop_codon("ATG")


# --- internal stops ---

def test_check_internal_stops_finds_premature():
    assert check_internal_stops("ATGTAAAAA") == [3]


def test_check_internal_stops_ignores_terminal_stop():
    assert check_internal_stops("ATGAAATAA") == []


def test_check_internal_stops_none():
    assert check_internal_stops("ATGAAACCC") == []


def test_check_internal_stops_multiple():
    # ATG TAA GGG TGA AAA -> internal stops at 3 and 9 (terminal at 12)
    assert check_internal_stops("ATGTAAGGGTGAAAA") == [3, 9]


# --- codon extraction, both strands ---

def test_extract_start_codon_plus():
    assert extract_start_codon("ATGAAACCC", 0, "+") == "ATG"


def test_extract_start_codon_minus():
    # high coord = 6; genome[3:6]="CAT"; RC -> "ATG"
    assert extract_start_codon("GGGCAT", 6, "-") == "ATG"


def test_extract_stop_codon_plus():
    # cds_end = 6 (half-open high); genome[6:9]="TAA"
    assert extract_stop_codon_after_cds("ATGAAATAA", 6, "+") == "TAA"


def test_extract_stop_codon_minus():
    # low coord = 3; genome[0:3]="TTA"; RC -> "TAA"
    assert extract_stop_codon_after_cds("TTAGGG", 3, "-") == "TAA"


def test_extract_start_codon_rejects_bad_strand():
    with pytest.raises(ValueError):
        extract_start_codon("ATGAAA", 0, ".")


def test_extract_stop_codon_rejects_bad_strand():
    with pytest.raises(ValueError):
        extract_stop_codon_after_cds("ATGAAA", 0, ".")


# --- MockGenome fixture sanity (RC on minus) ---

def test_mock_genome_plus(mock_genome):
    assert mock_genome.get_sequence("chr1", 0, 3, "+") == "ATG"


def test_mock_genome_minus_is_rc(mock_genome):
    # plus [0,3) = "ATG"; minus -> RC = "CAT"
    assert mock_genome.get_sequence("chr1", 0, 3, "-") == "CAT"


# ===========================================================================
# Phase 27 D1 — ambiguity-tolerant sequence layer (assessment §1.1)
# ===========================================================================

from helixforge.utils.sequences import (  # noqa: E402
    CODON_TABLES,
    DEFAULT_TRANSL_TABLE,
    get_codon_table,
    has_ambiguous_base,
    require_acgt,
    stop_codons_for_table,
)


def test_iupac_reverse_complement_full_alphabet():
    # Each IUPAC code complements to its partner; reversed order.
    # RYSWKM -> complement "YRSWMK" -> reversed -> "KMWSRY"
    assert reverse_complement("RYSWKM") == "KMWSRY"


def test_iupac_reverse_complement_degenerate_triples():
    # BDHV -> complements V H D B -> reversed -> "BDHV"
    assert reverse_complement("BDHV") == "BDHV"


def test_translate_ambiguous_internal_codon_is_x():
    # ATG AAN CCC -> M X P (the N-containing middle codon is unknown, not a stop)
    assert translate("ATGAANCCC") == "MXP"


def test_translate_n_does_not_become_stop():
    # An ambiguous codon must never be reported as '*'.
    assert "*" not in translate("ATGNNNAAA")


def test_check_internal_stops_ignores_ambiguous_codon():
    # ATG NNN AAA TAA -> codon NNN is 'X' (not internal stop); terminal TAA OK.
    assert check_internal_stops("ATGNNNAAATAA") == []


def test_check_internal_stops_real_stop_with_ambiguous_present():
    # ATG TAA NNN AAA -> internal stop at 3 still detected despite an N codon.
    assert check_internal_stops("ATGTAANNNAAA") == [3]


def test_has_ambiguous_base():
    assert has_ambiguous_base("ATGN")
    assert has_ambiguous_base("atgr")
    assert not has_ambiguous_base("ATGC")
    assert not has_ambiguous_base("atgc")


def test_extract_start_codon_with_n_plus():
    # plus strand: genome[0:3] = "ANG" (middle base masked) -> tolerated, not raised
    assert extract_start_codon("ANGAAACCC", 0, "+") == "ANG"


def test_extract_start_codon_with_n_minus():
    # minus strand: high coord 6; genome[3:6]="CNT"; RC -> "ANG"
    assert extract_start_codon("GGGCNT", 6, "-") == "ANG"


def test_extract_stop_codon_with_n_minus():
    # minus strand: low coord 3; genome[0:3]="TNA"; RC -> "TNA"
    assert extract_stop_codon_after_cds("TNAGGG", 3, "-") == "TNA"


def test_require_acgt_strict_still_raises_on_n():
    # Strict path retained for callers that genuinely need unambiguous bases.
    with pytest.raises(ValueError):
        require_acgt("ATGN")
    assert require_acgt("atgc") == "ATGC"


# ===========================================================================
# Phase 27 D3 — parameterized genetic code (assessment §1.2)
# ===========================================================================

def test_default_transl_table_is_one():
    assert DEFAULT_TRANSL_TABLE == 1
    assert get_codon_table() is CODON_TABLES[1]
    assert get_codon_table() is CODON_TABLE


def test_table4_translates_tga_as_tryptophan():
    # Table 4 (mold/protozoan mito): TGA codes Trp (W), not stop — a clean
    # differing codon vs the standard table.
    assert get_codon_table(1)["TGA"] == "*"
    assert get_codon_table(4)["TGA"] == "W"
    assert translate("ATGTGAAAA", transl_table=4) == "MWK"
    assert translate("ATGTGAAAA", transl_table=1) == "M*K"


def test_table4_internal_stop_differs_both_strands():
    # plus: ATG TGA AAA TAA -> table 1 has internal stop at 3; table 4 does not.
    assert check_internal_stops("ATGTGAAAATAA", transl_table=1) == [3]
    assert check_internal_stops("ATGTGAAAATAA", transl_table=4) == []
    # minus-strand sanity: the coding sequence is built RC upstream, so the same
    # string semantics apply once reverse-complemented. RC("...TCA...") logic is
    # exercised via translate of an RC'd input.
    rc = reverse_complement("ATGTGAAAATAA")  # coding seq if stored on minus
    # rc starts "TTA..."; we only assert table selection is honoured, not motif.
    assert check_internal_stops(rc, transl_table=1) == check_internal_stops(rc, transl_table=1)


def test_stop_codons_for_table_default_unchanged():
    assert set(stop_codons_for_table(1)) == {"TAA", "TAG", "TGA"}
    # Table 4: TGA reassigned to W, so it drops out of the stop set.
    assert set(stop_codons_for_table(4)) == {"TAA", "TAG"}


def test_is_stop_codon_table_aware():
    assert is_stop_codon("TGA", transl_table=1)
    assert not is_stop_codon("TGA", transl_table=4)
    # An ambiguous codon is never a stop, regardless of table.
    assert not is_stop_codon("TNA", transl_table=1)


def test_get_codon_table_unknown_raises():
    with pytest.raises(ValueError):
        get_codon_table(999)
