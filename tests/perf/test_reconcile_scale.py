"""Phase 21 D4 — M8 performance / scale regression (opt-in, ``slow``).

Proves the code.v3 hot-path rewrites (Phase 20 O(1) ``IdAllocator``; Phase 21
``JunctionIndex``) keep a large reconcile **sub-quadratic** and inside a generous
wall-clock / peak-memory envelope. This is the milestone-M8 evidence; the numbers
captured here are recorded in ``docs/CODE_V3_NOTES.md``.

It is **opt-in** — marked ``slow`` and deselected from the default
``pytest tests/`` run by ``addopts = -m 'not slow'`` (pyproject). Run it with::

    pytest tests/perf -q -m slow

The envelopes are deliberately loose (CI machines vary wildly); the meaningful
assertions are the *shape* checks: doubling the locus count must not quadruple
the time (rules out an O(N²) regression), and matching 50 K introns against a
100 K-junction set must finish in seconds (impossible with the old linear scan,
which would be ~5×10⁹ comparisons).
"""

from __future__ import annotations

import time
import tracemalloc

import pytest

from helixforge.reconcile.fallback import JunctionIndex, find_matching_junction
from helixforge.reconcile.mikado_integrate import IdAllocator, reconcile
from helixforge.reconcile.models import (
    Exon,
    HelixerLocus,
    Interval,
    LocusClassification,
    SpliceJunction,
)

# Phase 25 D3: the heavy 50 K-locus / 100 K-junction envelopes stay opt-in
# (``@pytest.mark.slow``, deselected by ``addopts = -m 'not slow'``), but the new
# ``test_reconcile_ci_budget`` below is intentionally *un*marked so it runs in the
# default ``pytest tests/`` and fails loudly if a change reintroduces O(N²) into
# the reconcile hot path. (Previously the whole module was module-level ``slow``.)

_SPACING = 1000  # disjoint loci, one per 1 kb slot


def _build_loci(n, strand="+"):
    """``n`` disjoint two-exon Helixer loci on one scaffold + SILENT classifications."""
    loci, cls = [], []
    for i in range(n):
        base = i * _SPACING
        gid = f"g{i:06d}"
        exons = [Exon(base, base + 200), Exon(base + 400, base + 600)]
        loci.append(HelixerLocus(gid, "chr1", base, base + 600, strand, exons=exons))
        cls.append(LocusClassification(gid, "SILENT"))
    return loci, cls


def _time_reconcile(n):
    loci, cls = _build_loci(n)
    t0 = time.perf_counter()
    genes, id_map, _ = reconcile(loci, cls, [], allocator=IdAllocator())
    elapsed = time.perf_counter() - t0
    assert len(genes) == n            # every locus carried through as a backstop gene
    assert len(id_map) == n           # one stable HFG per Helixer locus
    return elapsed


def test_reconcile_ci_budget():
    """Phase 25 D3 — CI-runnable perf-regression guard (NOT ``slow``).

    Runs in the default suite. Small enough to be fast (a few k loci ≈ tens of ms
    on the Phase 20 O(1) allocator) yet large enough that a reintroduced O(N²)
    allocator/overlap would blow the generous wall-clock budget and the
    doubling-ratio shape check. Envelopes are deliberately loose (CI machines vary
    by ~10×); the meaningful signal is *shape*: 2× loci must not ~4× the time.
    """
    _time_reconcile(1_000)               # warm import / allocation paths

    t_small = _time_reconcile(4_000)
    t_big = _time_reconcile(8_000)

    # Absolute budget: hugely generous (real run is ~27k loci; 8k is a fraction).
    # An O(N²) allocator at 8k loci would take many seconds; linear is ~10 ms.
    assert t_big < 5.0, f"8k-locus reconcile took {t_big:.3f}s (CI budget 5s)"

    # Shape check: only assert the ratio when t_small is large enough to time
    # reliably (sub-millisecond timings are dominated by noise on a busy CI box).
    if t_small > 0.02:
        ratio = t_big / t_small
        assert ratio < 3.0, f"reconcile scaled ~{ratio:.1f}× for 2× loci (quadratic?)"


@pytest.mark.slow
def test_reconcile_50k_is_subquadratic_and_bounded():
    # Warm up import/allocation paths so the ratio reflects the algorithm, not setup.
    _time_reconcile(2_000)

    # Timing is measured WITHOUT tracemalloc active — its per-allocation tracing
    # would dominate this allocation-heavy path and corrupt the ratio.
    t_small = _time_reconcile(25_000)
    t_big = _time_reconcile(50_000)
    ratio = t_big / max(t_small, 1e-6)

    # Peak memory is measured in a separate, untimed run.
    loci, cls = _build_loci(50_000)
    tracemalloc.start()
    reconcile(loci, cls, [], allocator=IdAllocator())
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_mb = peak / 1e6

    print(
        f"\n[M8] reconcile 25k={t_small:.2f}s 50k={t_big:.2f}s "
        f"ratio={ratio:.2f} peak={peak_mb:.0f}MB"
    )

    # Absolute envelope (very generous; the real Arabidopsis run is ~27k loci).
    assert t_big < 60.0, f"50k reconcile took {t_big:.2f}s (envelope 60s)"
    assert peak_mb < 2000.0, f"50k reconcile peaked at {peak_mb:.0f}MB (envelope 2GB)"

    # Shape check: doubling N must not quadruple the time. Quadratic ⇒ ~4×; linear
    # ⇒ ~2×. Only assert the ratio when t_small is large enough to time reliably.
    if t_small > 0.05:
        assert ratio < 3.0, f"reconcile scaled ~{ratio:.1f}× for 2× loci (quadratic?)"


@pytest.mark.slow
def test_junction_index_match_50k_introns_against_100k_junctions():
    # 100k junctions: exact entries at every 1 kb slot + decoys nearby.
    n = 50_000
    junctions = []
    for i in range(n):
        base = i * _SPACING
        junctions.append(SpliceJunction("chr1", base + 200, base + 400, "+", read_count=9))
        junctions.append(SpliceJunction("chr1", base + 210, base + 400, "+", read_count=4))

    t0 = time.perf_counter()
    index = JunctionIndex.from_junctions(junctions)
    build_s = time.perf_counter() - t0
    assert len(index) == 2 * n

    introns = [Interval(i * _SPACING + 200, i * _SPACING + 400) for i in range(n)]
    t0 = time.perf_counter()
    matched = 0
    for intron in introns:
        j = find_matching_junction(intron, index, "chr1", "+")
        if j is not None and j.read_count == 9:
            matched += 1
    query_s = time.perf_counter() - t0

    print(f"\n[M8] junction index build={build_s:.2f}s query50k={query_s:.2f}s")
    assert matched == n  # every intron found its exact, higher-read-count junction
    # A linear scan would be n×2n ≈ 5×10⁹ comparisons; the index must be seconds.
    assert query_s < 10.0, f"indexed 50k-intron match took {query_s:.2f}s"
