"""Phase 19 D2: the central constants module + its re-exports.

Every constant collected into ``helixforge.constants`` is re-exported from its
original module under its original public name. These tests pin both the
canonical value (against the literal it had before centralization) and the
re-export identity, so a value drift or a broken re-export fails loudly. Phase 19
is count-neutral — a moved constant with a new value is a bug.
"""

from __future__ import annotations

from helixforge import constants


def test_canonical_values_unchanged():
    # Literal values captured from the original sites (assesment §1.6).
    assert constants.NOVEL_ID_BASE == 90000
    assert constants.CROSS_CHECK_OVERLAP == 0.8
    assert constants.MIN_EXON_BP == 3
    assert constants.MIN_INTRON_BP == 20
    assert constants.PLAN_DEFAULT_MIN_BOUNDARY_GAP == 1000
    assert constants.SUGGEST_DEFAULT_MIN_BOUNDARY_GAP == 2000
    assert constants.DEFAULT_EXON_WEIGHT == 0.7
    assert constants.DEFAULT_INTRON_HIGH_CUTOFF == 0.5
    assert constants.STDERR_TAIL == 2000
    assert constants.TARGET_CHUNK_BP == 40_000_000
    assert constants.MEM_GB_PER_MB_SEQUENCE == 0.05
    assert constants.MIN_CHUNK_MEM_GB == 4
    assert constants.WALLTIME_MIN_PER_MB == 1.5
    assert constants.MIN_WALLTIME_MIN == 60


def test_reexports_equal_original_literals():
    # Each module keeps its original public name, sourced from constants.
    from helixforge.mikado.emit_external import (
        DEFAULT_EXON_WEIGHT,
        DEFAULT_INTRON_HIGH_CUTOFF,
    )
    from helixforge.parallel.plan import DEFAULT_MIN_BOUNDARY_GAP as PLAN_GAP
    from helixforge.parallel.suggest import (
        DEFAULT_MIN_BOUNDARY_GAP as SUGGEST_GAP,
        MEM_GB_PER_MB_SEQUENCE,
        MIN_CHUNK_MEM_GB,
        MIN_WALLTIME_MIN,
        TARGET_CHUNK_BP,
        WALLTIME_MIN_PER_MB,
    )
    from helixforge.prep._subprocess import _STDERR_TAIL
    from helixforge.reconcile.cds import _CROSS_CHECK_OVERLAP
    from helixforge.reconcile.fallback import MIN_EXON_BP, MIN_INTRON_BP
    from helixforge.reconcile.mikado_integrate import NOVEL_ID_BASE

    assert NOVEL_ID_BASE == constants.NOVEL_ID_BASE == 90000
    assert _CROSS_CHECK_OVERLAP == constants.CROSS_CHECK_OVERLAP == 0.8
    assert MIN_EXON_BP == constants.MIN_EXON_BP == 3
    assert MIN_INTRON_BP == constants.MIN_INTRON_BP == 20
    # Two distinct boundary-gap defaults: plan's (1000) and suggest's (2000).
    assert PLAN_GAP == constants.PLAN_DEFAULT_MIN_BOUNDARY_GAP == 1000
    assert SUGGEST_GAP == constants.SUGGEST_DEFAULT_MIN_BOUNDARY_GAP == 2000
    assert DEFAULT_EXON_WEIGHT == constants.DEFAULT_EXON_WEIGHT == 0.7
    assert DEFAULT_INTRON_HIGH_CUTOFF == constants.DEFAULT_INTRON_HIGH_CUTOFF == 0.5
    assert _STDERR_TAIL == constants.STDERR_TAIL == 2000
    assert TARGET_CHUNK_BP == constants.TARGET_CHUNK_BP == 40_000_000
    assert MEM_GB_PER_MB_SEQUENCE == constants.MEM_GB_PER_MB_SEQUENCE == 0.05
    assert MIN_CHUNK_MEM_GB == constants.MIN_CHUNK_MEM_GB == 4
    assert WALLTIME_MIN_PER_MB == constants.WALLTIME_MIN_PER_MB == 1.5
    assert MIN_WALLTIME_MIN == constants.MIN_WALLTIME_MIN == 60
