"""Phase 7 D2 — backstop intron junction correction (transaction pattern).

Concrete literal coordinates only (CLAUDE.md §12); both strands mandatory.
The transaction-revert tests are the most important in the repo (CLAUDE.md §6,
§12): any correction that would create a 0-width exon, a <3 bp exon, a <20 bp
intron, an inverted/overlapping exon, or an out-of-span exon MUST revert to
``None`` and leave the original transcript untouched.
"""

import pytest

from helixforge.reconcile.fallback import (
    CONTRADICTED,
    JunctionIndex,
    NO_EVIDENCE,
    SUPPORTED,
    apply_intron_corrections,
    classify_introns,
    find_contradicting_junction,
    find_matching_junction,
    refine_backstop_gene,
    _verify_bounds,
)
from helixforge.reconcile.models import (
    Exon,
    Interval,
    LocusClassification,
    ReconciledGene,
    SpliceJunction,
    TranscriptCandidate,
)

# Base structure: 3 exons, introns at (1200,1300) and (1500,1600).
BASE_EXONS = [(1000, 1200), (1300, 1500), (1600, 1800)]


def exons(bounds):
    return [Exon(s, e) for s, e in bounds]


def make_tx(exon_bounds=BASE_EXONS, strand="+", source="helixer", cds=None):
    from helixforge.reconcile.models import CDSSegment

    return TranscriptCandidate(
        transcript_id="HFG_00001.1",
        locus_id="HFG_00001",
        source=source,
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand=strand,
        exons=exons(exon_bounds),
        cds=[CDSSegment(*c) for c in cds] if cds else None,
    )


def make_gene(transcript, origin="helixer_backstop"):
    return ReconciledGene(
        gene_id="HFG_00001",
        seqid="chr1",
        start=transcript.start,
        end=transcript.end,
        strand=transcript.strand,
        tier=4,
        transcripts=[transcript],
        primary_transcript_id="HFG_00001.1",
        classification=LocusClassification("HFG_00001", "SILENT"),
        origin=origin,
    )


def jn(donor, acceptor, reads=5, strand="+", seqid="chr1"):
    return SpliceJunction(seqid, donor, acceptor, strand, read_count=reads)


def exon_bounds(tx):
    return [(e.start, e.end) for e in tx.exons]


def flag_names(gene):
    return {f.name for f in gene.flags}


# ---------------------------------------------------------------------------
# find_matching_junction
# ---------------------------------------------------------------------------

def test_matching_exact():
    intron = Interval(1200, 1300)
    j = find_matching_junction(intron, [jn(1200, 1300, 5)], "chr1", "+")
    assert j is not None and j.donor == 1200


def test_matching_below_min_reads():
    intron = Interval(1200, 1300)
    assert find_matching_junction(intron, [jn(1200, 1300, 2)], "chr1", "+", min_reads=3) is None


def test_matching_within_tolerance():
    intron = Interval(1200, 1300)
    assert find_matching_junction(intron, [jn(1201, 1299, 5)], "chr1", "+", tolerance=1) is not None


def test_matching_outside_tolerance():
    intron = Interval(1200, 1300)
    assert find_matching_junction(intron, [jn(1201, 1299, 5)], "chr1", "+", tolerance=0) is None


def test_matching_wrong_strand():
    intron = Interval(1200, 1300)
    assert find_matching_junction(intron, [jn(1200, 1300, 5, strand="-")], "chr1", "+") is None


def test_matching_wrong_seqid():
    intron = Interval(1200, 1300)
    assert find_matching_junction(intron, [jn(1200, 1300, 5, seqid="chr2")], "chr1", "+") is None


def test_matching_best_by_read_count():
    intron = Interval(1200, 1300)
    j = find_matching_junction(intron, [jn(1200, 1300, 5), jn(1200, 1300, 8)], "chr1", "+")
    assert j.read_count == 8


# ---------------------------------------------------------------------------
# find_contradicting_junction
# ---------------------------------------------------------------------------

def test_contradicting_shares_donor():
    intron = Interval(1200, 1300)
    j = find_contradicting_junction(intron, [jn(1200, 1320, 10)], "chr1", "+")
    assert j is not None and j.acceptor == 1320


def test_contradicting_shares_acceptor():
    intron = Interval(1200, 1300)
    j = find_contradicting_junction(intron, [jn(1180, 1300, 10)], "chr1", "+")
    assert j is not None and j.donor == 1180


def test_contradicting_full_match_not_contradiction():
    intron = Interval(1200, 1300)
    assert find_contradicting_junction(intron, [jn(1200, 1300, 10)], "chr1", "+") is None


def test_contradicting_no_shared_boundary_ignored():
    intron = Interval(1200, 1300)
    assert find_contradicting_junction(intron, [jn(1210, 1320, 10)], "chr1", "+") is None


def test_contradicting_best_by_read_count():
    intron = Interval(1200, 1300)
    j = find_contradicting_junction(
        intron, [jn(1200, 1320, 5), jn(1200, 1340, 9)], "chr1", "+"
    )
    assert j.read_count == 9


# ---------------------------------------------------------------------------
# classify_introns
# ---------------------------------------------------------------------------

def test_classify_all_supported():
    tx = make_tx()
    js = [jn(1200, 1300, 5), jn(1500, 1600, 5)]
    result = classify_introns(tx, js)
    assert [r[0] for r in result] == [SUPPORTED, SUPPORTED]


def test_classify_mixed():
    tx = make_tx()
    js = [jn(1200, 1300, 5), jn(1500, 1620, 8)]  # second contradicts
    result = classify_introns(tx, js)
    assert [r[0] for r in result] == [SUPPORTED, CONTRADICTED]


def test_classify_no_evidence():
    tx = make_tx()
    result = classify_introns(tx, [])
    assert [r[0] for r in result] == [NO_EVIDENCE, NO_EVIDENCE]


def test_classify_minus_strand():
    tx = make_tx(strand="-")
    js = [jn(1200, 1300, 5, strand="-"), jn(1500, 1600, 5, strand="-")]
    result = classify_introns(tx, js)
    assert [r[0] for r in result] == [SUPPORTED, SUPPORTED]


# ---------------------------------------------------------------------------
# apply_intron_corrections — success
# ---------------------------------------------------------------------------

def test_apply_single_correction_success():
    tx = make_tx()
    out = apply_intron_corrections(tx, [(0, 1200, 1320)])
    assert out is not None
    assert exon_bounds(out) == [(1000, 1200), (1320, 1500), (1600, 1800)]


def test_apply_multiple_corrections_success():
    tx = make_tx()
    out = apply_intron_corrections(tx, [(0, 1200, 1320), (1, 1480, 1600)])
    assert out is not None
    assert exon_bounds(out) == [(1000, 1200), (1320, 1480), (1600, 1800)]


def test_apply_empty_corrections_returns_same():
    tx = make_tx()
    assert apply_intron_corrections(tx, []) is tx


def test_apply_minus_strand_success():
    tx = make_tx(strand="-")
    out = apply_intron_corrections(tx, [(1, 1520, 1600)])
    assert out is not None
    assert exon_bounds(out) == [(1000, 1200), (1300, 1520), (1600, 1800)]


# ---------------------------------------------------------------------------
# apply_intron_corrections — TRANSACTION REVERTS (most important)
# ---------------------------------------------------------------------------

def test_revert_zero_width_exon():
    tx = make_tx()
    assert apply_intron_corrections(tx, [(0, 1000, 1300)]) is None  # exon0 -> (1000,1000)


def test_revert_short_exon_below_3bp():
    tx = make_tx()
    assert apply_intron_corrections(tx, [(0, 1002, 1300)]) is None  # exon0 -> (1000,1002)


def test_revert_short_intron_below_20bp():
    tx = make_tx()
    assert apply_intron_corrections(tx, [(0, 1290, 1300)]) is None  # intron0 = 10 bp


def test_revert_inverted_overlap():
    tx = make_tx()
    assert apply_intron_corrections(tx, [(0, 1350, 1300)]) is None  # exon0 end > exon1 start


def test_revert_inverted_exon():
    tx = make_tx()
    # exon1 becomes (1240,1500) but exon0 end pushed to 1250 -> overlap/inverted
    assert apply_intron_corrections(tx, [(0, 1250, 1240)]) is None


def test_revert_nonexistent_intron_index():
    tx = make_tx()
    assert apply_intron_corrections(tx, [(5, 1200, 1300)]) is None


def test_revert_atomic_one_bad_correction():
    tx = make_tx()
    # first correction valid, second creates a <20 bp intron -> WHOLE thing reverts
    out = apply_intron_corrections(tx, [(0, 1200, 1320), (1, 1480, 1495)])
    assert out is None


def test_revert_breaks_cds_within_exon():
    # CDS sits in exon1 (1300,1500); moving exon1 start past the CDS breaks
    # containment, so model construction raises -> transaction reverts.
    tx = make_tx(cds=[(1350, 1449)])
    assert apply_intron_corrections(tx, [(0, 1200, 1400)]) is None


def test_revert_leaves_original_untouched():
    tx = make_tx()
    before = exon_bounds(tx)
    apply_intron_corrections(tx, [(0, 1000, 1300)])  # reverts
    assert exon_bounds(tx) == before  # original transcript object unchanged


# ---------------------------------------------------------------------------
# _verify_bounds — out-of-span guard
# ---------------------------------------------------------------------------

def test_verify_bounds_ok():
    assert _verify_bounds([(1000, 1200), (1300, 1500)], 1000, 1500) is True


def test_verify_bounds_below_span():
    assert _verify_bounds([(990, 1200), (1300, 1500)], 1000, 1500) is False


def test_verify_bounds_above_span():
    assert _verify_bounds([(1000, 1200), (1300, 1600)], 1000, 1500) is False


# ---------------------------------------------------------------------------
# refine_backstop_gene
# ---------------------------------------------------------------------------

def test_refine_all_supported():
    gene = make_gene(make_tx())
    js = [jn(1200, 1300, 5), jn(1500, 1600, 5)]
    out = refine_backstop_gene(gene, js)
    assert "ALL_JUNCTIONS_SUPPORTED" in flag_names(out)
    assert exon_bounds(out.transcripts[0]) == BASE_EXONS
    assert out.transcripts[0].junction_support_fraction == 1.0


def test_refine_no_evidence():
    gene = make_gene(make_tx())
    out = refine_backstop_gene(gene, [])
    assert "NO_JUNCTION_SUPPORT" in flag_names(out)
    assert out.transcripts[0].junction_support_fraction == 0.0


def test_refine_partial_support():
    gene = make_gene(make_tx())
    js = [jn(1200, 1300, 5)]  # only intron0 supported
    out = refine_backstop_gene(gene, js)
    assert "PARTIAL_JUNCTION_SUPPORT" in flag_names(out)
    assert out.transcripts[0].junction_support_fraction == 0.5


def test_refine_corrects_contradicted_intron():
    gene = make_gene(make_tx())
    js = [jn(1200, 1320, 10), jn(1500, 1600, 5)]  # intron0 contradicted, intron1 ok
    out = refine_backstop_gene(gene, js)
    assert exon_bounds(out.transcripts[0]) == [(1000, 1200), (1320, 1500), (1600, 1800)]
    assert "ALL_JUNCTIONS_SUPPORTED" in flag_names(out)
    assert out.transcripts[0].junction_support_fraction == 1.0


def test_refine_correction_reverts_keeps_structure():
    gene = make_gene(make_tx())
    # contradicting junction would create a <20 bp intron (1200..1215) ->
    # correction reverts; structure stays, that intron counts as unsupported.
    js = [jn(1200, 1215, 10), jn(1500, 1600, 5)]
    out = refine_backstop_gene(gene, js)
    assert exon_bounds(out.transcripts[0]) == BASE_EXONS  # unchanged
    assert "PARTIAL_JUNCTION_SUPPORT" in flag_names(out)


def test_refine_mikado_origin_untouched():
    gene = make_gene(make_tx(source="mikado"), origin="mikado_1to1")
    out = refine_backstop_gene(gene, [jn(1200, 1320, 10)])
    assert out is gene


def test_refine_single_exon_untouched():
    gene = make_gene(make_tx(exon_bounds=[(1000, 1200)]))
    out = refine_backstop_gene(gene, [jn(1050, 1150, 10)])
    assert out is gene  # no introns -> never insert one


def test_refine_minus_strand_corrects():
    gene = make_gene(make_tx(strand="-"))
    js = [jn(1200, 1300, 5, strand="-"), jn(1500, 1620, 9, strand="-")]
    out = refine_backstop_gene(gene, js)
    assert exon_bounds(out.transcripts[0]) == [(1000, 1200), (1300, 1500), (1620, 1800)]
    assert "ALL_JUNCTIONS_SUPPORTED" in flag_names(out)


# ==========================================================================
# Phase 28 D1 — canonicity is a hard gate on backstop junction correction
# (assessment §1.3 / CLAUDE.md §6: the one place coordinates are mutated)
# ==========================================================================

def jn_c(donor, acceptor, canonical, reads=10, strand="+", seqid="chr1"):
    return SpliceJunction(
        seqid, donor, acceptor, strand, read_count=reads, canonical=canonical
    )


def test_refine_rejects_non_canonical_correction():
    # The same contradicting junction that test_refine_corrects_contradicted_intron
    # accepts (when canonical is unknown) is now REJECTED because its motif is
    # non-canonical: structure is left untouched, the intron counts unsupported,
    # and NON_CANONICAL_SPLICE is flagged.
    gene = make_gene(make_tx())
    js = [jn_c(1200, 1320, "non-canonical", 10), jn(1500, 1600, 5)]
    out = refine_backstop_gene(gene, js)
    assert exon_bounds(out.transcripts[0]) == BASE_EXONS  # NOT corrected
    assert "NON_CANONICAL_SPLICE" in flag_names(out)
    assert "PARTIAL_JUNCTION_SUPPORT" in flag_names(out)  # intron0 unsupported


def test_refine_accepts_canonical_correction():
    # A canonical (GT-AG) correction target IS applied — the gate only blocks
    # non-canonical motifs. canonical=None (legacy) is also accepted (covered by
    # test_refine_corrects_contradicted_intron), keeping the old path identical.
    gene = make_gene(make_tx())
    js = [jn_c(1200, 1320, "GT-AG", 10), jn(1500, 1600, 5)]
    out = refine_backstop_gene(gene, js)
    assert exon_bounds(out.transcripts[0]) == [(1000, 1200), (1320, 1500), (1600, 1800)]
    assert "ALL_JUNCTIONS_SUPPORTED" in flag_names(out)
    assert "NON_CANONICAL_SPLICE" not in flag_names(out)


def test_refine_minus_strand_rejects_non_canonical():
    # Both strands (CLAUDE.md §12): a non-canonical correction on the minus
    # strand is rejected just the same.
    gene = make_gene(make_tx(strand="-"))
    js = [jn(1200, 1300, 5, strand="-"),
          jn_c(1500, 1620, "non-canonical", 9, strand="-")]
    out = refine_backstop_gene(gene, js)
    assert exon_bounds(out.transcripts[0]) == BASE_EXONS  # intron1 not corrected
    assert "NON_CANONICAL_SPLICE" in flag_names(out)
    # intron0 matched exactly -> partial support overall
    assert "PARTIAL_JUNCTION_SUPPORT" in flag_names(out)


# ==========================================================================
# Phase 21 D2 — JunctionIndex equals the linear scan (both strands, tol 0 & >0)
# ==========================================================================

# A mixed junction set across strands/scaffolds with exact + near-miss entries,
# read-count ties, and contradicting (one-boundary) junctions.
_JSET = [
    jn(1200, 1300, 5),                 # exact, +
    jn(1200, 1300, 8),                 # exact tie-break (higher reads wins)
    jn(1201, 1299, 7),                 # near (tolerance window), +
    jn(1200, 1500, 9),                 # shares donor only -> contradiction, +
    jn(1260, 1300, 6),                 # shares acceptor only -> contradiction, +
    jn(1200, 1300, 4, strand="-"),     # exact on the minus strand
    jn(1199, 1300, 6, strand="-"),     # near donor, minus
    jn(1200, 1300, 5, seqid="chr2"),   # other scaffold
    jn(1700, 1850, 2),                 # below default min_reads
]


def _probe_introns():
    return [Interval(1200, 1300), Interval(1450, 1500), Interval(1700, 1850)]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("tolerance", [0, 2])
def test_index_matching_equals_linear_scan(strand, tolerance):
    index = JunctionIndex.from_junctions(_JSET)
    for intron in _probe_introns():
        for min_reads in (1, 3, 5):
            linear = find_matching_junction(
                intron, _JSET, "chr1", strand, tolerance=tolerance, min_reads=min_reads
            )
            indexed = find_matching_junction(
                intron, index, "chr1", strand, tolerance=tolerance, min_reads=min_reads
            )
            assert linear is indexed  # same object, incl. the read-count tie-break


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("tolerance", [0, 2])
def test_index_contradiction_equals_linear_scan(strand, tolerance):
    index = JunctionIndex.from_junctions(_JSET)
    for intron in _probe_introns():
        for min_reads in (1, 3):
            linear = find_contradicting_junction(
                intron, _JSET, "chr1", strand, tolerance=tolerance, min_reads=min_reads
            )
            indexed = find_contradicting_junction(
                intron, index, "chr1", strand, tolerance=tolerance, min_reads=min_reads
            )
            assert linear is indexed


def test_index_empty_and_unknown_bucket():
    index = JunctionIndex.from_junctions([])
    assert len(index) == 0
    assert find_matching_junction(Interval(1200, 1300), index, "chrX", "+") is None
    assert find_contradicting_junction(Interval(1200, 1300), index, "chrX", "+") is None


def test_refine_backstop_index_matches_list_plus_strand():
    # The whole backstop refine must be identical whether fed a list or an index.
    js = [jn(1200, 1320, 7), jn(1500, 1600, 9)]  # one contradiction, one exact
    via_list = refine_backstop_gene(make_gene(make_tx()), js)
    via_index = refine_backstop_gene(make_gene(make_tx()), JunctionIndex.from_junctions(js))
    assert exon_bounds(via_index.transcripts[0]) == exon_bounds(via_list.transcripts[0])
    assert flag_names(via_index) == flag_names(via_list)
    tx_l, tx_i = via_list.transcripts[0], via_index.transcripts[0]
    assert tx_i.junction_support_fraction == tx_l.junction_support_fraction


def test_refine_backstop_index_matches_list_minus_strand():
    js = [jn(1200, 1300, 5, strand="-"), jn(1500, 1620, 9, strand="-")]
    via_list = refine_backstop_gene(make_gene(make_tx(strand="-")), js)
    via_index = refine_backstop_gene(
        make_gene(make_tx(strand="-")), JunctionIndex.from_junctions(js)
    )
    assert exon_bounds(via_index.transcripts[0]) == exon_bounds(via_list.transcripts[0])
    assert flag_names(via_index) == flag_names(via_list)


def test_classify_introns_index_equals_list():
    tx = make_tx()  # introns (1200,1300) and (1500,1500)? see BASE_EXONS
    js = [jn(1200, 1300, 6), jn(1500, 1620, 9)]
    via_list = classify_introns(tx, js)
    via_index = classify_introns(tx, JunctionIndex.from_junctions(js))
    assert [s for s, _ in via_list] == [s for s, _ in via_index]
    assert [j for _, j in via_list] == [j for _, j in via_index]


# ---------------------------------------------------------------------------
# Intrinsic CDS must not block a junction correction (strip + re-project)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_intrinsic_cds_does_not_block_correction(strand):
    # A backstop gene now carries its Helixer CDS *before* junction correction.
    # A correction that shrinks exon1's end (1200 -> 1150) would put the old CDS
    # (1000-1180) outside the corrected exon and, pre-fix, revert the whole fix on
    # CDS-containment. The structure must still be corrected (CDS is stripped for
    # the fix, then re-projected onto the corrected exons).
    tx = make_tx(cds=[(1000, 1180, 0), (1300, 1420, 0)], strand=strand)
    gene = make_gene(tx)
    js = [jn(1150, 1300, 5, strand), jn(1500, 1600, 5, strand)]  # contradict i0, support i1
    out = refine_backstop_gene(gene, js)
    new = out.transcripts[0]
    # correction applied despite the CDS
    assert exon_bounds(new)[0] == (1000, 1150)
    # CDS re-projected onto the corrected exons (1000-1150 clip + 1300-1420 = 270)
    assert new.cds is not None
    assert new.total_cds_length == 270
    assert "ALL_JUNCTIONS_SUPPORTED" in flag_names(out)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_cds_preserved_when_correction_misses_it(strand):
    # A correction far from the CDS leaves the CDS untouched (fast re-attach path).
    tx = make_tx(cds=[(1050, 1200, 0), (1300, 1450, 0)], strand=strand)  # 150+150=300
    gene = make_gene(tx)
    js = [jn(1200, 1300, 5, strand), jn(1500, 1620, 5, strand)]  # support i0, correct i1
    out = refine_backstop_gene(gene, js)
    new = out.transcripts[0]
    assert exon_bounds(new)[2] == (1620, 1800)  # exon3 start corrected
    assert [(c.start, c.end) for c in new.cds] == [(1050, 1200), (1300, 1450)]


# ---------------------------------------------------------------------------
# A correction that GROWS an exon outward past a flush CDS boundary strands that
# boundary inside the grown exon. The stale CDS is contained (so it survives the
# model's within-exon check) but no longer coherent — its CDS-derived intron no
# longer matches any exon intron, which crashes Mikado's finalizer
# (assert len(cds_introns) > 0). Such a CDS must be rejected, not emitted.
# ---------------------------------------------------------------------------

def _cds_introns_match_exon_introns(tx):
    """Replicates Mikado's coherence requirement: every CDS-derived intron is an
    exon intron. Used to assert emitted coding transcripts are finalizer-safe."""
    cds = sorted((c.start, c.end) for c in (tx.cds or []))
    if len(cds) < 2:
        return True
    ex = sorted((e.start, e.end) for e in tx.exons)
    exon_introns = {(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)}
    return all(
        (cds[i][1], cds[i + 1][0]) in exon_introns for i in range(len(cds) - 1)
    )


@pytest.mark.parametrize("strand", ["+", "-"])
def test_outward_grow_strands_flush_cds_is_dropped(strand):
    # intron0 acceptor 1300 -> 1250 (exon2 grows outward: 1300..1500 -> 1250..1500,
    # intron shrinks (1200,1300) -> (1200,1250)). The CDS segment in exon2 starts
    # flush at the OLD boundary 1300, so after correction it is 50 bp inside the
    # grown exon: CDS intron (1200,1300) no longer matches exon intron (1200,1250).
    # Clipping cannot extend it to 1250, so the incoherent CDS is dropped.
    tx = make_tx(cds=[(1050, 1200, 0), (1300, 1450, 0)], strand=strand)  # 150+150=300
    gene = make_gene(tx)
    js = [jn(1200, 1250, 10, strand), jn(1500, 1600, 5, strand)]  # contradict i0, support i1
    out = refine_backstop_gene(gene, js)
    new = out.transcripts[0]
    # structure corrected (evidence honored)…
    assert exon_bounds(new)[1] == (1250, 1500)
    assert "ALL_JUNCTIONS_SUPPORTED" in flag_names(out)
    # …but the stranded intrinsic CDS is rejected rather than emitted incoherently.
    assert new.cds is None


@pytest.mark.parametrize("strand", ["+", "-"])
def test_inward_shrink_reclips_cds_stays_coherent(strand):
    # The mirror case (already handled): intron0 acceptor 1300 -> 1303 shrinks
    # exon2 inward by 3 bp; the CDS segment (1300-1450) now overhangs into the new
    # intron and is re-clipped flush to the new boundary 1303, staying coherent and
    # mod-3 (150 + 147 = 297) — so it is kept, not dropped.
    tx = make_tx(cds=[(1050, 1200, 0), (1300, 1450, 0)], strand=strand)
    gene = make_gene(tx)
    js = [jn(1200, 1303, 10, strand), jn(1500, 1600, 5, strand)]  # contradict i0, support i1
    out = refine_backstop_gene(gene, js)
    new = out.transcripts[0]
    assert exon_bounds(new)[1] == (1303, 1500)
    assert new.cds is not None
    assert [(c.start, c.end) for c in new.cds] == [(1050, 1200), (1303, 1450)]
    assert _cds_introns_match_exon_introns(new)  # finalizer-safe (CDS intron == exon intron)
