"""Poster step 2b — matched-locus + junction-level isoform accuracy.

Pure-Python metric (no ``mikado`` binary): tested directly on synthetic gene
dicts shaped like ``GFF3Parser.parse_genes_generic()`` output. Concrete literal
coordinates; both strands. See ``prompts/helixforge-poster/02b_*.md``.
"""

from __future__ import annotations

from typing import Any

from helixforge.bench.isoform import (
    intron_chain,
    intron_set,
    match_genes_1to1,
    matched_locus_isoform_accuracy,
)
from helixforge.reconcile.models import Exon


# ---------------------------------------------------------------------------
# Builders — gene dicts in the parse_genes_generic shape
# ---------------------------------------------------------------------------


def _gene(
    gene_id: str,
    seqid: str,
    strand: str,
    transcripts: list[tuple[str, list[tuple[int, int]]]],
) -> dict[str, Any]:
    """Build a gene dict; ``transcripts`` = list of ``(transcript_id, exon_pairs)``."""
    txs: list[dict[str, Any]] = []
    starts: list[int] = []
    ends: list[int] = []
    for tid, exon_pairs in transcripts:
        exons = [Exon(s, e) for s, e in exon_pairs]
        txs.append({"transcript_id": tid, "exons": exons, "cds": None})
        starts += [s for s, _ in exon_pairs]
        ends += [e for _, e in exon_pairs]
    return {
        "gene_id": gene_id,
        "seqid": seqid,
        "start": min(starts),
        "end": max(ends),
        "strand": strand,
        "transcripts": txs,
    }


# ---------------------------------------------------------------------------
# D1 — intron_chain / intron_set
# ---------------------------------------------------------------------------


def test_intron_chain_known_coords():
    exons = [Exon(100, 200), Exon(300, 400), Exon(500, 600)]
    assert intron_chain(exons) == ((200, 300), (400, 500))


def test_intron_chain_terminus_independent_plus():
    # Same internal introns, different terminal-exon lengths -> equal chain.
    short_termini = [Exon(100, 200), Exon(300, 400), Exon(500, 600)]
    long_termini = [Exon(10, 200), Exon(300, 400), Exon(500, 999)]
    assert intron_chain(short_termini) == intron_chain(long_termini) == ((200, 300), (400, 500))


def test_intron_chain_terminus_independent_minus():
    # Strand-agnostic by construction: coords identical -> equal chain on minus too.
    a = [Exon(1000, 1100), Exon(1300, 1400), Exon(1600, 1700)]
    b = [Exon(950, 1100), Exon(1300, 1400), Exon(1600, 1800)]
    assert intron_chain(a) == intron_chain(b) == ((1100, 1300), (1400, 1600))


def test_intron_chain_single_exon_is_empty():
    assert intron_chain([Exon(100, 200)]) == ()
    assert intron_chain([]) == ()


def test_intron_chain_unsorted_input():
    exons = [Exon(500, 600), Exon(100, 200), Exon(300, 400)]
    assert intron_chain(exons) == ((200, 300), (400, 500))


def test_intron_set_union_over_transcripts():
    g = _gene(
        "G", "chr1", "+",
        [
            ("G.1", [(100, 200), (300, 400), (500, 600)]),  # introns 200-300, 400-500
            ("G.2", [(100, 200), (500, 600)]),              # intron 200-500
        ],
    )
    assert intron_set(g["transcripts"]) == {(200, 300), (400, 500), (200, 500)}


# ---------------------------------------------------------------------------
# D1 — match_genes_1to1
# ---------------------------------------------------------------------------


def test_match_picks_best_reciprocal_overlap_above_threshold():
    # pred strongly overlaps RG1, only a tiny (<10%) sliver of RG2 -> matches RG1.
    pred = _gene("PG", "chr1", "+", [("PG.1", [(100, 200), (300, 400), (500, 600)])])
    rg1 = _gene("RG1", "chr1", "+", [("RG1.1", [(100, 200), (300, 400), (500, 600)])])
    rg2 = _gene("RG2", "chr1", "+", [("RG2.1", [(590, 600)])])  # 10 bp inside pred
    matches = match_genes_1to1([pred], [rg1, rg2])
    assert len(matches) == 1
    assert matches[0][1]["gene_id"] == "RG1"


def test_match_excludes_split_merge_and_counts():
    # One pred spanning two ref genes (both above threshold) -> dropped + counted.
    pred = _gene(
        "PGX", "chr1", "+",
        [("PGX.1", [(100, 200), (300, 400), (1000, 1100), (1300, 1400)])],
    )
    rga = _gene("RGa", "chr1", "+", [("RGa.1", [(100, 200), (300, 400)])])
    rgb = _gene("RGb", "chr1", "+", [("RGb.1", [(1000, 1100), (1300, 1400)])])
    matches = match_genes_1to1([pred], [rga, rgb])
    assert matches == []
    res = matched_locus_isoform_accuracy([pred], [rga, rgb])
    assert res["split_merge_excluded_n"] == 1
    assert res["matched_genes_n"] == 0


def test_match_strand_mismatch_no_match():
    pred = _gene("PG", "chr1", "+", [("PG.1", [(100, 200), (300, 400)])])
    ref = _gene("RG", "chr1", "-", [("RG.1", [(100, 200), (300, 400)])])
    assert match_genes_1to1([pred], [ref]) == []


def test_match_deterministic_order_by_overlap_then_id():
    # Two clean pairs; PG_b has the larger overlap so it sorts first.
    pg_a = _gene("PG_a", "chr1", "+", [("PG_a.1", [(100, 200)])])
    rg_a = _gene("RG_a", "chr1", "+", [("RG_a.1", [(100, 200)])])
    pg_b = _gene("PG_b", "chr1", "+", [("PG_b.1", [(5000, 5300)])])
    rg_b = _gene("RG_b", "chr1", "+", [("RG_b.1", [(5000, 5300)])])
    matches = match_genes_1to1([pg_a, pg_b], [rg_a, rg_b])
    assert [p["gene_id"] for p, _ in matches] == ["PG_b", "PG_a"]


# ---------------------------------------------------------------------------
# D2 — matched-locus alternative-isoform precision / recall
# ---------------------------------------------------------------------------


def _matched_multiiso_pair_plus():
    """A 1:1-matched multi-iso plus-strand pair with one terminus-shifted alt."""
    pred = _gene(
        "PG1", "chr1", "+",
        [
            ("PG1.1", [(100, 200), (300, 400), (500, 600)]),  # primary, chain (200,300),(400,500)
            ("PG1.2", [(100, 200), (500, 600)]),              # alt, chain (200,500)  (exon skip)
        ],
    )
    ref = _gene(
        "RG1", "chr1", "+",
        [
            ("RG1.1", [(100, 200), (300, 400), (500, 600)]),  # chain (200,300),(400,500)
            ("RG1.2", [(150, 200), (500, 650)]),              # chain (200,500) — diff termini, same intron
        ],
    )
    return pred, ref


def test_matched_locus_alt_precision_hit():
    pred, ref = _matched_multiiso_pair_plus()
    res = matched_locus_isoform_accuracy([pred], [ref])
    # One predicted alternative (PG1.2); its chain (200,500) matches RG1.2 -> 1/1.
    assert res["pred_alt_isoforms_n"] == 1
    assert res["pred_alt_isoforms_hit_n"] == 1
    assert res["alt_isoform_precision"] == 1.0
    assert res["matched_multiiso_genes_n"] == 1


def test_matched_locus_alt_precision_miss():
    # Alternative chain present in neither ref isoform -> precision 0 over denom 1.
    pred = _gene(
        "PG1", "chr1", "+",
        [
            ("PG1.1", [(100, 200), (300, 400), (500, 600)]),
            ("PG1.2", [(100, 200), (320, 400), (500, 600)]),  # chain (200,320),(400,500) — bogus
        ],
    )
    ref = _gene(
        "RG1", "chr1", "+",
        [
            ("RG1.1", [(100, 200), (300, 400), (500, 600)]),
            ("RG1.2", [(100, 200), (500, 600)]),
        ],
    )
    res = matched_locus_isoform_accuracy([pred], [ref])
    assert res["pred_alt_isoforms_n"] == 1
    assert res["pred_alt_isoforms_hit_n"] == 0
    assert res["alt_isoform_precision"] == 0.0


def test_matched_locus_precision_denominator_is_predicted_alternatives_only():
    # Three predicted isoforms -> denominator counts the two NON-primary ones only.
    pred = _gene(
        "PG1", "chr1", "+",
        [
            ("PG1.1", [(100, 200), (300, 400), (500, 600)]),  # primary (excluded from denom)
            ("PG1.2", [(100, 200), (500, 600)]),
            ("PG1.3", [(100, 200), (300, 400), (700, 800)]),
        ],
    )
    ref = _gene(
        "RG1", "chr1", "+",
        [
            ("RG1.1", [(100, 200), (300, 400), (500, 600)]),
            ("RG1.2", [(100, 200), (500, 600)]),
        ],
    )
    res = matched_locus_isoform_accuracy([pred], [ref])
    assert res["pred_alt_isoforms_n"] == 2  # PG1.2, PG1.3 — not the .1 primary


def test_matched_locus_recall_fair_denominator_excludes_global_catalog():
    # The reference catalog has many extra multi-iso loci that are NOT matched
    # (different seqid). The recall denominator must count only ref isoforms in
    # the one matched locus (2), never the global catalog.
    pred, ref_matched = _matched_multiiso_pair_plus()
    extra = [
        _gene(f"RX{i}", "chr9", "+",
              [(f"RX{i}.1", [(10, 20), (40, 50)]), (f"RX{i}.2", [(10, 20), (60, 70)])])
        for i in range(50)
    ]
    res = matched_locus_isoform_accuracy([pred], [ref_matched, *extra])
    assert res["ref_isoforms_n"] == 2          # only RG1.1 + RG1.2
    assert res["ref_isoforms_hit_n"] == 2      # both reproduced by a predicted chain
    assert res["isoform_recall"] == 1.0


def test_matched_locus_single_iso_pred_excluded_from_d2():
    # pred has only 1 transcript -> not a multi-iso locus -> excluded from D2,
    # but it is still a 1:1 match (counts in matched_genes_n / junction D3).
    pred = _gene("PG1", "chr1", "+", [("PG1.1", [(100, 200), (300, 400), (500, 600)])])
    ref = _gene(
        "RG1", "chr1", "+",
        [("RG1.1", [(100, 200), (300, 400), (500, 600)]),
         ("RG1.2", [(100, 200), (500, 600)])],
    )
    res = matched_locus_isoform_accuracy([pred], [ref])
    assert res["matched_genes_n"] == 1
    assert res["matched_multiiso_genes_n"] == 0
    assert res["pred_alt_isoforms_n"] == 0


# ---------------------------------------------------------------------------
# D2b — matched-locus PRIMARY (canonical) precision; the TRaCE-on metric
# ---------------------------------------------------------------------------

# A reference-correct chain and a bogus chain; which one is elected primary
# (``.1``) is exactly what TRaCE changes.
_CORRECT = [(100, 200), (300, 400), (500, 600)]   # chain (200,300),(400,500) == RG1.1
_BOGUS = [(100, 200), (320, 400), (500, 600)]     # chain (200,320),(400,500) — in neither ref
_REF_PLUS = _gene(
    "RG1", "chr1", "+",
    [("RG1.1", [(100, 200), (300, 400), (500, 600)]),
     ("RG1.2", [(100, 200), (500, 600)])],
)


def test_primary_precision_hit_when_correct_is_canonical():
    # Correct chain elected primary (.1): primary precision 1/1; the bogus alt misses.
    pred = _gene("PG1", "chr1", "+", [("PG1.1", _CORRECT), ("PG1.2", _BOGUS)])
    res = matched_locus_isoform_accuracy([pred], [_REF_PLUS])
    assert res["pred_primary_n"] == 1
    assert res["pred_primary_hit_n"] == 1
    assert res["primary_precision"] == 1.0
    assert res["alt_isoform_precision"] == 0.0     # only the bogus alt remains


def test_primary_precision_miss_when_bogus_is_canonical():
    # SAME two transcripts, but the bogus chain is elected primary (.1): the
    # primary now misses (0/1) and the correct chain scores as the alt (1/1).
    # This is the TRaCE swap: re-electing the canonical moves a hit between buckets.
    pred = _gene("PG1", "chr1", "+", [("PG1.1", _BOGUS), ("PG1.2", _CORRECT)])
    res = matched_locus_isoform_accuracy([pred], [_REF_PLUS])
    assert res["pred_primary_hit_n"] == 0
    assert res["primary_precision"] == 0.0
    assert res["alt_isoform_precision"] == 1.0
    # Conservation: primary_hit + alt_hit is invariant under re-election (==1 here).
    assert res["pred_primary_hit_n"] + res["pred_alt_isoforms_hit_n"] == 1


def test_primary_precision_minus_strand_swap():
    # Minus strand; correct chain (1100,1300),(1400,1600) == RG2.1.
    correct = [(1000, 1100), (1300, 1400), (1600, 1700)]
    bogus = [(1000, 1100), (1320, 1400), (1600, 1700)]   # chain (1100,1320),(1400,1600)
    ref = _gene(
        "RG2", "chr2", "-",
        [("RG2.1", [(1000, 1100), (1300, 1400), (1600, 1700)]),
         ("RG2.2", [(1000, 1100), (1600, 1700)])],
    )
    canonical_correct = _gene("PG2", "chr2", "-", [("PG2.1", correct), ("PG2.2", bogus)])
    canonical_bogus = _gene("PG2", "chr2", "-", [("PG2.1", bogus), ("PG2.2", correct)])
    assert matched_locus_isoform_accuracy([canonical_correct], [ref])["primary_precision"] == 1.0
    assert matched_locus_isoform_accuracy([canonical_bogus], [ref])["primary_precision"] == 0.0


# ---------------------------------------------------------------------------
# D3 — junction-level (intron-set) precision / recall
# ---------------------------------------------------------------------------


def test_junction_level_pr_known_overlap():
    # pred introns {200-300, 400-500, 200-500}; ref introns {200-300, 400-500, 700-800}.
    # intersection = {200-300, 400-500} = 2.  precision 2/3, recall 2/3.
    pred = _gene(
        "PG1", "chr1", "+",
        [("PG1.1", [(100, 200), (300, 400), (500, 600)]),
         ("PG1.2", [(100, 200), (500, 600)])],
    )
    ref = _gene(
        "RG1", "chr1", "+",
        [("RG1.1", [(100, 200), (300, 400), (500, 600)]),
         ("RG1.2", [(100, 200), (300, 400), (700, 800)])],
    )
    res = matched_locus_isoform_accuracy([pred], [ref])
    assert res["pred_introns_n"] == 3
    assert res["ref_introns_n"] == 3
    assert res["pred_introns_hit_n"] == 2
    assert res["ref_introns_hit_n"] == 2
    assert abs(res["junction_precision"] - 2 / 3) < 1e-9
    assert abs(res["junction_recall"] - 2 / 3) < 1e-9
    assert abs(res["junction_f1"] - 2 / 3) < 1e-9


def test_junction_tolerance_admits_off_by_one():
    # Donor shifted by 1 bp: a miss at tol=0, a hit at tol=1.
    pred = _gene("PG1", "chr1", "+", [("PG1.1", [(100, 200), (301, 400)])])  # intron 200-301
    ref = _gene("RG1", "chr1", "+", [("RG1.1", [(100, 200), (300, 400)])])   # intron 200-300
    strict = matched_locus_isoform_accuracy([pred], [ref], junction_tolerance=0)
    assert strict["pred_introns_hit_n"] == 0
    assert strict["junction_precision"] == 0.0
    loose = matched_locus_isoform_accuracy([pred], [ref], junction_tolerance=1)
    assert loose["pred_introns_hit_n"] == 1
    assert loose["junction_precision"] == 1.0


def test_junction_tolerance_admits_off_by_one_in_chain_match():
    # Same 1 bp shift inside a multi-iso locus: alt chain matches only with tol>=1.
    pred = _gene(
        "PG1", "chr1", "+",
        [("PG1.1", [(100, 200), (300, 400), (500, 600)]),
         ("PG1.2", [(100, 200), (501, 600)])],  # alt intron 200-501
    )
    ref = _gene(
        "RG1", "chr1", "+",
        [("RG1.1", [(100, 200), (300, 400), (500, 600)]),
         ("RG1.2", [(100, 200), (500, 600)])],  # ref intron 200-500
    )
    assert matched_locus_isoform_accuracy([pred], [ref])["alt_isoform_precision"] == 0.0
    assert matched_locus_isoform_accuracy(
        [pred], [ref], junction_tolerance=1)["alt_isoform_precision"] == 1.0


# ---------------------------------------------------------------------------
# Worked both-strand fixtures — every metric reproduced by hand
# ---------------------------------------------------------------------------


def test_worked_fixture_plus_strand_all_metrics():
    pred, ref = _matched_multiiso_pair_plus()
    res = matched_locus_isoform_accuracy([pred], [ref])
    assert res["matched_genes_n"] == 1
    assert res["matched_multiiso_genes_n"] == 1
    assert res["split_merge_excluded_n"] == 0
    # D2: 1 predicted alt, hit; recall 2/2.
    assert res["alt_isoform_precision"] == 1.0
    assert res["isoform_recall"] == 1.0
    assert res["isoform_f1"] == 1.0
    # D2b: the primary PG1.1 chain == RG1.1 -> primary precision 1/1.
    assert res["pred_primary_n"] == 1
    assert res["primary_precision"] == 1.0
    # D3: union introns {200-300, 400-500, 200-500} identical both sides.
    assert res["pred_introns_n"] == 3
    assert res["ref_introns_n"] == 3
    assert res["junction_precision"] == 1.0
    assert res["junction_recall"] == 1.0


def test_worked_fixture_minus_strand_all_metrics():
    # Minus strand, terminus-shifted alternative; chain matching is strand-agnostic.
    pred = _gene(
        "PG2", "chr2", "-",
        [
            ("PG2.1", [(1000, 1100), (1300, 1400), (1600, 1700)]),  # chain (1100,1300),(1400,1600)
            ("PG2.2", [(1000, 1100), (1600, 1700)]),                # alt chain (1100,1600)
        ],
    )
    ref = _gene(
        "RG2", "chr2", "-",
        [
            ("RG2.1", [(1000, 1100), (1300, 1400), (1600, 1700)]),
            ("RG2.2", [(1050, 1100), (1600, 1750)]),                # chain (1100,1600), diff termini
        ],
    )
    res = matched_locus_isoform_accuracy([pred], [ref])
    assert res["matched_genes_n"] == 1
    assert res["alt_isoform_precision"] == 1.0   # PG2.2 (1100,1600) == RG2.2 (1100,1600)
    assert res["primary_precision"] == 1.0       # PG2.1 chain == RG2.1
    assert res["isoform_recall"] == 1.0
    assert res["pred_introns_n"] == 3            # {1100-1300, 1400-1600, 1100-1600}
    assert res["ref_introns_n"] == 3
    assert res["junction_precision"] == 1.0
    assert res["junction_recall"] == 1.0


def test_no_matches_yields_nan_metrics():
    # Disjoint loci -> no matches -> precision/recall undefined (nan), counts zero.
    import math

    pred = _gene("PG", "chr1", "+", [("PG.1", [(100, 200), (300, 400)])])
    ref = _gene("RG", "chr3", "+", [("RG.1", [(100, 200), (300, 400)])])
    res = matched_locus_isoform_accuracy([pred], [ref])
    assert res["matched_genes_n"] == 0
    assert math.isnan(res["alt_isoform_precision"])
    assert math.isnan(res["junction_recall"])
