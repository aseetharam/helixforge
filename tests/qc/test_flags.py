"""Tests for the QC flag registry (Phase 0). Floor: 12."""

import pytest

from helixforge.qc import flags as flagmod
from helixforge.qc.flags import (
    ALL_FLAGS,
    CDS_DISAGREE,
    HELIXER_ONLY,
    LOCUS_SPLIT,
    NO_EXPRESSION,
    NO_START,
    as_event_flag,
    dedup_flags,
    get_flag,
    tier_flag,
)
from helixforge.reconcile.models import QCFlag


def test_all_flags_are_qcflags():
    assert all(isinstance(f, QCFlag) for f in ALL_FLAGS.values())


def test_all_flags_keyed_by_name():
    for name, flag in ALL_FLAGS.items():
        assert name == flag.name


def test_no_duplicate_names():
    names = [f.name for f in ALL_FLAGS.values()]
    assert len(names) == len(set(names))


def test_registry_contains_required_locus_flags():
    for name in ("LOCUS_SPLIT", "LOCUS_MERGE", "MERGE_REJECTED", "NOVEL_LOCUS", "CDS_DISAGREE"):
        assert name in ALL_FLAGS


def test_registry_contains_structure_flags():
    for name in ("NO_START", "NO_STOP", "INTERNAL_STOP", "SHORT_CDS", "SHORT_EXON", "LONG_INTRON"):
        assert name in ALL_FLAGS


def test_registry_contains_splice_flags():
    for name in ("ALL_JUNCTIONS_SUPPORTED", "PARTIAL_JUNCTION_SUPPORT", "NO_JUNCTION_SUPPORT"):
        assert name in ALL_FLAGS


def test_registry_count_matches_constants():
    # 3 evidence + 2 confidence + 10 structure + 2 homology + 4 splice + 8 locus = 29
    # (Phase 18: +BACKSTOP_RESCUED structure, +FROM_TANGLED_LOCUS locus;
    #  Phase 27: +AMBIGUOUS_CODON structure;
    #  Phase 28: +NON_CANONICAL_SPLICE splice;
    #  Phase 29: +PSEUDOGENE_CANDIDATE locus;
    #  Phase 31: +DOMAIN_COMPLETE homology;
    #  Phase 32: +VARIANT_IMPACTED structure;
    #  TPM/biotype/TE fix: +PUTATIVE_CODING evidence, +TE_OVERLAP locus)
    assert len(ALL_FLAGS) == 29


def test_constant_categories_valid():
    assert HELIXER_ONLY.category == "evidence"
    assert NO_START.category == "structure"
    assert LOCUS_SPLIT.category == "locus"
    assert CDS_DISAGREE.category == "locus"


def test_get_flag_returns_constant():
    assert get_flag("LOCUS_SPLIT") is LOCUS_SPLIT


def test_get_flag_missing_raises_keyerror():
    with pytest.raises(KeyError):
        get_flag("NOT_A_FLAG")


def test_tier_flag_factory():
    f = tier_flag(2)
    assert f.name == "TIER_2"
    assert isinstance(f, QCFlag)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_tier_flag_all_valid(n):
    assert tier_flag(n).name == f"TIER_{n}"


def test_tier_flag_rejects_bad_tier():
    with pytest.raises(ValueError):
        tier_flag(5)


def test_as_event_flag_factory():
    f = as_event_flag("ES")
    assert f.name == "AS_ES"
    assert f.category == "splice"


def test_as_event_flag_rejects_bad_kind():
    with pytest.raises(ValueError):
        as_event_flag("RI")


def test_all_constants_pass_validation():
    # constructing the module already validated them; re-confirm severity domain
    for f in ALL_FLAGS.values():
        assert f.severity in ("INFO", "WARNING", "ERROR", "CRITICAL")


def test_registry_is_copy_not_internal():
    # mutating ALL_FLAGS must not corrupt the private registry
    snapshot = dict(ALL_FLAGS)
    ALL_FLAGS["BOGUS"] = None
    assert "BOGUS" not in flagmod._REGISTRY
    ALL_FLAGS.clear()
    ALL_FLAGS.update(snapshot)


# ---------------------------------------------------------------------------
# Phase 19 D1: dedup_flags unified here (was duplicated 4x in reconcile/*).
# ---------------------------------------------------------------------------

def test_dedup_flags_order_preserving_first_wins():
    # Two distinct QCFlag objects sharing a name: dedup is by name, first-wins,
    # so the *first* instance is the one retained (identity check).
    first = QCFlag("DUP", "structure", "INFO", "first")
    second = QCFlag("DUP", "structure", "INFO", "second")
    out = dedup_flags([NO_START, first, NO_EXPRESSION, second, NO_START])
    # first occurrence of each name is kept, in input order; later dups dropped
    assert [f.name for f in out] == ["NO_START", "DUP", "NO_EXPRESSION"]
    assert out[1] is first  # first-wins (not `second`)


def test_dedup_flags_idempotent_and_empty():
    once = dedup_flags([HELIXER_ONLY, NO_START, HELIXER_ONLY])
    twice = dedup_flags(once)
    assert [f.name for f in twice] == [f.name for f in once] == ["HELIXER_ONLY", "NO_START"]
    assert dedup_flags([]) == []
