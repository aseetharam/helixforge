"""Isoform-level accuracy metric."""

from __future__ import annotations

import bisect
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from helixforge.bench.wrappers import run_mikado_compare, run_omark
from helixforge.io.gff import GFF3Parser
from helixforge.reconcile.as_events import overlap_bases, reciprocal_overlap
from helixforge.utils.logging import get_logger
from helixforge.utils.regions import internal_to_gff3

_log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Parsing helpers — count isoforms, identify multi-isoform loci
# ---------------------------------------------------------------------------


def _parse_genes(gff3: str | Path) -> list[dict[str, Any]]:
    """Parse a GFF3 into generic gene dicts (internal coords). See GFF3Parser."""
    return GFF3Parser(str(gff3)).parse_genes_generic()


def _transcript_counts(genes: list[dict[str, Any]]) -> list[int]:
    """Per-gene transcript (isoform) counts, in gene order."""
    return [len(g["transcripts"]) for g in genes]


def _multiiso_genes(
    genes: list[dict[str, Any]], min_isoforms: int
) -> list[dict[str, Any]]:
    """Genes with ``min_isoforms`` or more transcripts (the AS-bearing loci)."""
    return [g for g in genes if len(g["transcripts"]) >= min_isoforms]


# ---------------------------------------------------------------------------
# Subset GFF3 writer (I/O boundary — internal → 1-based inclusive GFF3)
# ---------------------------------------------------------------------------


def _write_genes_gff3(genes: list[dict[str, Any]], out_path: str | Path) -> int:
    """Write ``genes`` (generic gene dicts) as a minimal valid GFF3.

    Re-emits gene → mRNA → exon, CDS with proper ``ID``/``Parent`` links so
    ``mikado compare`` sees the same structure. Coordinates are converted from
    internal 0-based half-open to 1-based inclusive (the only place this module
    touches coordinates). Returns the number of genes written.
    """
    out_path = Path(out_path)
    lines: list[str] = ["##gff-version 3"]
    for g in sorted(genes, key=lambda x: (x["seqid"], x["start"], x["end"])):
        gid = g["gene_id"]
        g_start, g_end = internal_to_gff3(g["start"], g["end"])
        lines.append(
            f"{g['seqid']}\thf\tgene\t{g_start}\t{g_end}\t.\t{g['strand']}\t.\tID={gid}"
        )
        for tx in g["transcripts"]:
            tid = tx["transcript_id"]
            exons = sorted(tx["exons"], key=lambda e: (e.start, e.end))
            # mRNA span = union of its exons (fall back to gene span if none).
            if exons:
                t_lo, t_hi = internal_to_gff3(exons[0].start, exons[-1].end)
            else:
                t_lo, t_hi = g_start, g_end
            lines.append(
                f"{g['seqid']}\thf\tmRNA\t{t_lo}\t{t_hi}\t.\t{g['strand']}\t.\t"
                f"ID={tid};Parent={gid}"
            )
            for i, e in enumerate(exons, start=1):
                e_lo, e_hi = internal_to_gff3(e.start, e.end)
                lines.append(
                    f"{g['seqid']}\thf\texon\t{e_lo}\t{e_hi}\t.\t{g['strand']}\t.\t"
                    f"ID={tid}.exon{i};Parent={tid}"
                )
            for i, c in enumerate(
                sorted(tx["cds"] or [], key=lambda s: (s.start, s.end)), start=1
            ):
                c_lo, c_hi = internal_to_gff3(c.start, c.end)
                lines.append(
                    f"{g['seqid']}\thf\tCDS\t{c_lo}\t{c_hi}\t.\t{g['strand']}\t{c.phase}\t"
                    f"ID={tid}.CDS{i};Parent={tid}"
                )
    out_path.write_text("\n".join(lines) + "\n")
    return len(genes)


# ---------------------------------------------------------------------------
# Level selection from a parsed mikado-compare .stats dict
# ---------------------------------------------------------------------------

# mikado emits several "Transcript level (...)" rows; the >=80% base-F1 row is the
# conventional comparison metric. Fall back through the others if a version omits it.
_TRANSCRIPT_KEYS = (
    "transcript_80_base_f1",
    "transcript_95_base_f1",
    "transcript_stringent",
)
_EMPTY = {"sn": float("nan"), "pr": float("nan"), "f1": float("nan")}


def _pick_transcript_level(levels: dict[str, dict[str, float]]) -> dict[str, float]:
    """The canonical transcript-level Sn/Pr/F1 (>=80% base F1 if present)."""
    for key in _TRANSCRIPT_KEYS:
        if key in levels:
            return levels[key]
    return dict(_EMPTY)


def _select(levels: dict[str, dict[str, float]]) -> dict[str, Any]:
    """Pull the isoform-relevant levels out of a full parsed .stats dict."""
    return {
        "intron_chain": levels.get("intron_chain", dict(_EMPTY)),
        "transcript": _pick_transcript_level(levels),
        "levels": levels,
    }


# ---------------------------------------------------------------------------
# Isoform-count summary
# ---------------------------------------------------------------------------


def _summary(counts: list[int], min_isoforms: int) -> dict[str, float]:
    """mean/median/multi-isoform-gene count for a list of per-gene isoform counts."""
    if not counts:
        return {
            "mean": float("nan"),
            "median": float("nan"),
            "multiiso_genes": 0,
            "total_genes": 0,
        }
    return {
        "mean": float(statistics.mean(counts)),
        "median": float(statistics.median(counts)),
        "multiiso_genes": sum(1 for c in counts if c >= min_isoforms),
        "total_genes": len(counts),
    }


# ---------------------------------------------------------------------------
# Matched-locus + junction-level isoform accuracy (poster step 2b)
#
# The poster-defensible isoform numbers. Two problems with the step-2 subset F1
# (25.84): it conflates "loci HelixForge did not split into isoforms" with
# "isoforms it missed inside loci it did handle", and its recall denominator is
# the whole Araport11 multi-isoform catalog (long-read-defined). The fix is to
# condition on **1:1-matched genes** and match on **internal splice structure
# only** (intron chains/junctions), which is terminus-independent by construction
# — Araport11 termini come from long-read 5'/3' data we do not have, so requiring
# identical termini would tank precision for a reason unrelated to splicing.
# Pure Python: no ``mikado`` binary, fully unit-testable on synthetic GFF3s.
# ---------------------------------------------------------------------------


def intron_chain(exons: list[Any]) -> tuple[tuple[int, int], ...]:
    """Ordered internal introns (gaps between consecutive exons), internal coords.

    The intron chain is the ordered list of ``(donor, acceptor)`` internal
    introns. It is **terminus-independent by construction**: two isoforms with
    identical internal introns but different terminal-exon lengths produce the
    same chain. A single-exon transcript has no introns ⇒ empty tuple.

    ``exons`` is a list of objects/tuples carrying ``.start``/``.end`` (Exon) or
    indexable as ``(start, end)``; 0-based half-open, so the gap between exon
    ``[s1, e1)`` and the next exon ``[s2, e2)`` is the intron ``[e1, s2)``.
    """
    es = sorted(
        (e.start, e.end) if hasattr(e, "start") else (e[0], e[1]) for e in exons
    )
    return tuple((es[i][1], es[i + 1][0]) for i in range(len(es) - 1))


def intron_set(transcripts: list[dict[str, Any]]) -> set[tuple[int, int]]:
    """Union of all introns across a gene's transcripts (for junction-level P/R)."""
    out: set[tuple[int, int]] = set()
    for tx in transcripts:
        out.update(intron_chain(tx["exons"]))
    return out


def _gene_footprint(gene: dict[str, Any]) -> list[tuple[int, int]]:
    """Disjoint genomic footprint (merged exon union) of all the gene's exons.

    ``reciprocal_overlap``/``overlap_bases`` assume disjoint interval lists, but
    a gene's per-transcript exons overlap each other — so merge the exon union
    into disjoint intervals first. Falls back to the gene span if it has no exons.
    """
    ivs = sorted(
        (e.start, e.end) if hasattr(e, "start") else (e[0], e[1])
        for tx in gene["transcripts"]
        for e in tx["exons"]
    )
    merged: list[tuple[int, int]] = []
    for s, e in ivs:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged or [(gene["start"], gene["end"])]


def _primary_tid(gene: dict[str, Any]) -> str | None:
    """The predicted primary transcript id (``.1`` convention), else first.

    Only the **prediction** side has a known primary; ``parse_genes_generic`` does
    not surface the ``primary=true`` attribute, so we use the stable ``.N`` id
    convention (primary = ``.1``). Used solely to split predicted *alternative*
    (non-primary) isoforms for the precision numerator — never on the reference
    (we do not invent a reference canonical transcript).
    """
    txs = gene["transcripts"]
    for t in txs:
        if str(t["transcript_id"]).endswith(".1"):
            return str(t["transcript_id"])
    return str(txs[0]["transcript_id"]) if txs else None


def _match_genes(
    pred_genes: list[dict[str, Any]],
    ref_genes: list[dict[str, Any]],
    min_reciprocal: float,
) -> tuple[list[tuple[dict[str, Any], dict[str, Any]]], int]:
    """1:1 reciprocal-overlap gene correspondence; ``(matches, n_split_merge)``.

    Same-seqid, same-strand only. Edges are gene pairs whose exon-union footprints
    reach ``reciprocal_overlap >= min_reciprocal``. A gene touched by more than one
    above-threshold edge is part of a split/merge (a pred over >1 ref, or a ref
    over >1 pred) and is **excluded** — split/merge is a gene-boundary concern, not
    an isoform-quality one. Remaining edges are strictly 1:1 (each side degree 1).
    ``n_split_merge`` is the count of prediction genes that had a qualifying
    overlap but were dropped this way. Deterministic: matches are ordered by
    descending overlap, tie-broken by ``gene_id``.
    """
    pred_fp = {id(p): _gene_footprint(p) for p in pred_genes}

    # Bucket refs by (seqid, strand), sorted by start, with a bisect window.
    raw: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in ref_genes:
        raw[(r["seqid"], r["strand"])].append(r)
    buckets: dict[
        tuple[str, str],
        tuple[list[int], list[dict[str, Any]], list[list[tuple[int, int]]], int],
    ] = {}
    for key, lst in raw.items():
        lst.sort(key=lambda g: g["start"])
        fps = [_gene_footprint(r) for r in lst]
        starts = [r["start"] for r in lst]
        max_span = max((r["end"] - r["start"] for r in lst), default=0)
        buckets[key] = (starts, lst, fps, max_span)

    edges: list[tuple[int, dict[str, Any], dict[str, Any]]] = []
    for p in pred_genes:
        bucket = buckets.get((p["seqid"], p["strand"]))
        if bucket is None:
            continue
        starts, genes_sorted, fps, max_span = bucket
        lo = bisect.bisect_left(starts, p["start"] - max_span)
        hi = bisect.bisect_right(starts, p["end"])
        p_fp = pred_fp[id(p)]
        for idx in range(lo, hi):
            r = genes_sorted[idx]
            if r["end"] <= p["start"] or r["start"] >= p["end"]:
                continue  # spans do not overlap
            if reciprocal_overlap(p_fp, fps[idx]) >= min_reciprocal:
                edges.append((overlap_bases(p_fp, fps[idx]), p, r))

    pred_deg = Counter(id(p) for _, p, _ in edges)
    ref_deg = Counter(id(r) for _, _, r in edges)
    clean = [e for e in edges if pred_deg[id(e[1])] == 1 and ref_deg[id(e[2])] == 1]
    clean.sort(key=lambda t: (-t[0], str(t[1]["gene_id"]), str(t[2]["gene_id"])))
    matches = [(p, r) for _, p, r in clean]

    preds_with_edge = {id(p) for _, p, _ in edges}
    matched_pred = {id(p) for p, _ in matches}
    n_split_merge = len(preds_with_edge) - len(matched_pred)
    return matches, n_split_merge


def match_genes_1to1(
    pred_genes: list[dict[str, Any]],
    ref_genes: list[dict[str, Any]],
    *,
    min_reciprocal: float = 0.10,
) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    """Greedy 1:1 reciprocal-overlap gene correspondence (split/merge excluded).

    See :func:`_match_genes`. Returns only the clean 1:1 ``(pred_gene, ref_gene)``
    pairs; split/merge loci are dropped (their count is reported by the
    matched-locus metric as ``split_merge_excluded_n``).
    """
    matches, _ = _match_genes(pred_genes, ref_genes, min_reciprocal=min_reciprocal)
    return matches


def _chains_match(
    a: tuple[tuple[int, int], ...],
    b: tuple[tuple[int, int], ...],
    tol: int,
) -> bool:
    """Two intron chains match (exact when ``tol==0``; ±``tol`` bp per splice site)."""
    if len(a) != len(b):
        return False
    if tol == 0:
        return a == b
    return all(
        abs(ai[0] - bi[0]) <= tol and abs(ai[1] - bi[1]) <= tol for ai, bi in zip(a, b)
    )


def _chain_matches_any(
    chain: tuple[tuple[int, int], ...],
    others: list[tuple[tuple[int, int], ...]],
    tol: int,
) -> bool:
    return any(_chains_match(chain, o, tol) for o in others)


def _junction_overlap(
    pred_introns: set[tuple[int, int]],
    ref_introns: set[tuple[int, int]],
    tol: int,
) -> tuple[int, int]:
    """``(pred_matched, ref_matched)`` intron counts (exact set-∩ when ``tol==0``)."""
    if tol == 0:
        inter = len(pred_introns & ref_introns)
        return inter, inter
    pm = sum(
        1
        for p in pred_introns
        if any(abs(p[0] - r[0]) <= tol and abs(p[1] - r[1]) <= tol for r in ref_introns)
    )
    rm = sum(
        1
        for r in ref_introns
        if any(
            abs(p[0] - r[0]) <= tol and abs(p[1] - r[1]) <= tol for p in pred_introns
        )
    )
    return pm, rm


def _ratio(num: int, den: int) -> float:
    return float(num) / den if den else float("nan")


def _harmonic_f1(precision: float, recall: float) -> float:
    if math.isnan(precision) or math.isnan(recall) or (precision + recall) == 0:
        return float("nan")
    return 2.0 * precision * recall / (precision + recall)


def matched_locus_isoform_accuracy(
    pred_genes: list[dict[str, Any]],
    ref_genes: list[dict[str, Any]],
    *,
    min_isoforms: int = 2,
    min_reciprocal: float = 0.10,
    junction_tolerance: int = 0,
) -> dict[str, Any]:
    """Matched-locus alternative-isoform + junction-level accuracy (poster step 2b).

    Conditions on 1:1-matched, same-strand genes (see :func:`match_genes_1to1`) and
    matches **internal intron chains/junctions only** (terminus-independent).

    * **Matched-locus alt-isoform precision** (the headline-quality number): over
      matched loci where *both* sides have ``>= min_isoforms`` transcripts, the
      fraction of predicted **alternative** (non-primary) isoforms whose intron
      chain exactly matches some reference isoform's chain. "When HelixForge
      proposes an alternative isoform, how often is it a real Araport11 chain."
    * **Matched-locus primary precision** (``primary_precision``): over the same
      multi-iso matched loci, the fraction of predicted **primary/canonical** (the
      ``.1``) transcripts whose intron chain matches some reference chain — "is the
      *elected canonical* transcript real?". This is the metric that responds to
      TRaCE primary re-election (``--trace-primary``): TRaCE promotes the
      evidence-best transcript to ``.1``, so it *raises* this precision while
      mechanically lowering the alt-precision (the promoted, often-correct chain
      leaves the non-primary pool). Denominator = one primary per multi-iso locus.
    * **Matched-locus isoform recall** (the *fair* recall): over the same matched
      multi-iso loci, the fraction of *all* reference isoforms whose chain is
      matched by some predicted isoform. Denominator = reference isoforms in
      shared genes, **not** the whole 10,737-locus catalog.
    * **Junction-level precision/recall** over *all* 1:1-matched genes (not only
      multi-iso): introns pooled per matched gene, ``|pred∩ref| / |pred|`` and
      ``/ |ref|``. Granular and terminus-independent — credits individual correct
      splice events; one wrong intron does not zero a 10-intron transcript.

    ``junction_tolerance`` (default 0 = exact) allows ±N bp on each splice site,
    for sensitivity checks. Returns a dict of metrics and their raw denominators
    (every reported value traces to a count).
    """
    matches, n_split_merge = _match_genes(
        pred_genes, ref_genes, min_reciprocal=min_reciprocal
    )

    # D2 — matched-locus alternative-isoform precision/recall (multi-iso both sides).
    pred_alt_total = pred_alt_hit = 0
    ref_iso_total = ref_iso_hit = 0
    pred_primary_total = pred_primary_hit = 0
    n_multiiso_matched = 0
    for p, r in matches:
        if len(p["transcripts"]) < min_isoforms or len(r["transcripts"]) < min_isoforms:
            continue
        n_multiiso_matched += 1
        ref_chains = [intron_chain(t["exons"]) for t in r["transcripts"]]
        pred_chains = [intron_chain(t["exons"]) for t in p["transcripts"]]
        primary = _primary_tid(p)
        for t in p["transcripts"]:
            if str(t["transcript_id"]) == primary:
                # D2b — primary/canonical precision: is the *elected* canonical
                # transcript a real Araport11 chain? This is the metric TRaCE moves
                # (it changes which transcript is primary); the alt-precision below
                # scores only the remaining non-primary isoforms.
                pred_primary_total += 1
                if _chain_matches_any(
                    intron_chain(t["exons"]), ref_chains, junction_tolerance
                ):
                    pred_primary_hit += 1
                continue
            pred_alt_total += 1
            if _chain_matches_any(
                intron_chain(t["exons"]), ref_chains, junction_tolerance
            ):
                pred_alt_hit += 1
        for rch in ref_chains:
            ref_iso_total += 1
            if _chain_matches_any(rch, pred_chains, junction_tolerance):
                ref_iso_hit += 1

    alt_precision = _ratio(pred_alt_hit, pred_alt_total)
    iso_recall = _ratio(ref_iso_hit, ref_iso_total)
    primary_precision = _ratio(pred_primary_hit, pred_primary_total)

    # D3 — junction-level (intron-set) precision/recall over ALL matched genes.
    pred_intr_total = pred_intr_hit = 0
    ref_intr_total = ref_intr_hit = 0
    for p, r in matches:
        ps = intron_set(p["transcripts"])
        rs = intron_set(r["transcripts"])
        pm, rm = _junction_overlap(ps, rs, junction_tolerance)
        pred_intr_total += len(ps)
        pred_intr_hit += pm
        ref_intr_total += len(rs)
        ref_intr_hit += rm

    junction_precision = _ratio(pred_intr_hit, pred_intr_total)
    junction_recall = _ratio(ref_intr_hit, ref_intr_total)

    return {
        "matched_genes_n": len(matches),
        "matched_multiiso_genes_n": n_multiiso_matched,
        "split_merge_excluded_n": n_split_merge,
        "alt_isoform_precision": alt_precision,
        "primary_precision": primary_precision,
        "isoform_recall": iso_recall,
        "isoform_f1": _harmonic_f1(alt_precision, iso_recall),
        "junction_precision": junction_precision,
        "junction_recall": junction_recall,
        "junction_f1": _harmonic_f1(junction_precision, junction_recall),
        # Raw denominators / numerators — every reported value is traceable.
        "pred_alt_isoforms_n": pred_alt_total,
        "pred_alt_isoforms_hit_n": pred_alt_hit,
        "pred_primary_n": pred_primary_total,
        "pred_primary_hit_n": pred_primary_hit,
        "ref_isoforms_n": ref_iso_total,
        "ref_isoforms_hit_n": ref_iso_hit,
        "pred_introns_n": pred_intr_total,
        "pred_introns_hit_n": pred_intr_hit,
        "ref_introns_n": ref_intr_total,
        "ref_introns_hit_n": ref_intr_hit,
    }


# ---------------------------------------------------------------------------
# The metric
# ---------------------------------------------------------------------------


def isoform_accuracy(
    reconciled_gff3: str | Path,
    reference_gff3: str | Path,
    out_prefix: str | Path,
    *,
    mikado_bin: str = "mikado",
    min_isoforms: int = 2,
    helixer_gff3: str | Path | None = None,
    proteins_fa: str | Path | None = None,
    omadb: str | Path | None = None,
    omark_bin: str = "omark",
) -> dict[str, Any]:
    """Isoform-level accuracy of ``reconciled_gff3`` against ``reference_gff3``.

    Restricts both annotations to **multi-isoform loci** (genes with
    ``min_isoforms`` or more transcripts) and runs ``mikado compare`` on that
    subset, returning intron-chain and transcript-level Sn/Pr/F1 for the subset
    *and* for the whole set. Also returns the isoform-count summary (mean/median
    isoforms per gene for the Helixer input vs the HelixForge output, and the
    number of genes that gained isoforms) and, when an OMA DB is supplied,
    OMArk's spurious/fragmented-isoform numbers.

    Parameters
    ----------
    reconciled_gff3, reference_gff3
        HelixForge output and the reference (Araport11) GFF3.
    out_prefix
        Prefix for the ``mikado compare`` outputs; ``<prefix>_multiiso.stats`` and
        ``<prefix>_whole.stats`` are written, and the multi-isoform subset GFF3s
        next to them.
    min_isoforms
        Threshold for "multi-isoform" (default 2 — any AS at all).
    helixer_gff3
        Optional Helixer input GFF3 for the real before-counts. When omitted the
        Helixer side is taken as exactly 1 transcript/gene (its documented design,
        Helixer emits one transcript per gene), so ``counts['helixer_*']`` is
        reported as 1.0.
    proteins_fa, omadb, omark_bin
        When both ``proteins_fa`` and ``omadb`` are given, OMArk is run on the
        protein set and its numbers attached under ``result['omark']``.

    Returns
    -------
    dict with keys ``whole``, ``multiiso``, ``counts`` and (optional) ``omark``.
    """
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    pred_genes = _parse_genes(reconciled_gff3)
    ref_genes = _parse_genes(reference_gff3)

    # Matched-locus + junction-level isoform accuracy (poster step 2b). Pure
    # Python (no mikado); computed before the mikado comparisons so it is the
    # poster-defensible isoform number even on the multi-isoform-subset path.
    matched = matched_locus_isoform_accuracy(
        pred_genes,
        ref_genes,
        min_isoforms=min_isoforms,
    )

    # Whole-set comparison.
    whole_levels = run_mikado_compare(
        reference_gff3,
        reconciled_gff3,
        str(out_prefix) + "_whole",
        mikado_bin=mikado_bin,
    )

    # Multi-isoform subset comparison.
    pred_multi = _multiiso_genes(pred_genes, min_isoforms)
    ref_multi = _multiiso_genes(ref_genes, min_isoforms)
    pred_sub = str(out_prefix) + "_pred_multiiso.gff3"
    ref_sub = str(out_prefix) + "_ref_multiiso.gff3"
    _write_genes_gff3(pred_multi, pred_sub)
    _write_genes_gff3(ref_multi, ref_sub)
    _log.info(
        "isoform_accuracy: %d/%d prediction and %d/%d reference multi-isoform loci "
        "(min_isoforms=%d)",
        len(pred_multi),
        len(pred_genes),
        len(ref_multi),
        len(ref_genes),
        min_isoforms,
    )
    multi_levels = run_mikado_compare(
        ref_sub,
        pred_sub,
        str(out_prefix) + "_multiiso",
        mikado_bin=mikado_bin,
    )

    # Isoform-count summary.
    hf_counts = _transcript_counts(pred_genes)
    ref_counts = _transcript_counts(ref_genes)
    hf_sum = _summary(hf_counts, min_isoforms)
    ref_sum = _summary(ref_counts, min_isoforms)
    if helixer_gff3 is not None:
        hx_sum = _summary(_transcript_counts(_parse_genes(helixer_gff3)), min_isoforms)
    else:
        # Helixer emits exactly one transcript per gene by design.
        hx_sum = {
            "mean": 1.0,
            "median": 1.0,
            "multiiso_genes": 0,
            "total_genes": hf_sum["total_genes"],
        }

    counts = {
        "helixer_mean": hx_sum["mean"],
        "helixer_median": hx_sum["median"],
        "helixforge_mean": hf_sum["mean"],
        "helixforge_median": hf_sum["median"],
        "reference_mean": ref_sum["mean"],
        "reference_median": ref_sum["median"],
        # Genes that gained isoforms relative to Helixer's one-per-gene baseline:
        # any HelixForge output gene with >= min_isoforms transcripts.
        "genes_gained_isoforms": hf_sum["multiiso_genes"],
        "helixforge_multiiso_genes": hf_sum["multiiso_genes"],
        "reference_multiiso_genes": ref_sum["multiiso_genes"],
    }

    result: dict[str, Any] = {
        "whole": _select(whole_levels),
        "multiiso": {
            **_select(multi_levels),
            "n_pred_genes": len(pred_multi),
            "n_ref_genes": len(ref_multi),
        },
        "counts": counts,
        "matched": matched,
    }

    if proteins_fa is not None and omadb is not None:
        result["omark"] = run_omark(
            proteins_fa,
            omadb,
            str(out_prefix) + "_omark",
            omark_bin=omark_bin,
        )

    return result


# ---------------------------------------------------------------------------
# Flatten → benchmark.tsv rows
# ---------------------------------------------------------------------------


def isoform_result_to_rows(result: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten an :func:`isoform_accuracy` result into ``benchmark.tsv`` rows.

    Each row is ``{tool: "isoform", metric, value, status: "ok"}`` with a clear
    ``isoform_`` metric prefix; the ``_multiiso`` suffix marks the hard subset.
    Every value traces back to a ``mikado compare`` ``.stats`` (the whole/subset
    runs) or to the parsed isoform counts.
    """
    rows: list[dict[str, Any]] = []

    def _emit(metric: str, value: Any) -> None:
        rows.append(
            {"tool": "isoform", "metric": metric, "value": float(value), "status": "ok"}
        )

    for level, suffix in (("intron_chain", ""), ("transcript", "")):
        for axis in ("sn", "pr", "f1"):
            _emit(f"isoform_{level}_{axis}", result["whole"][level][axis])
            _emit(f"isoform_{level}_{axis}_multiiso", result["multiiso"][level][axis])

    _emit("isoform_multiiso_pred_genes", result["multiiso"]["n_pred_genes"])
    _emit("isoform_multiiso_ref_genes", result["multiiso"]["n_ref_genes"])

    # Matched-locus + junction-level rows (poster step 2b). Headline = the
    # alt-isoform precision-on-predicted; junction P/R/F1 are the granular numbers.
    matched = result.get("matched")
    if matched is not None:
        _emit("isoform_matched_locus_alt_precision", matched["alt_isoform_precision"])
        _emit("isoform_matched_locus_primary_precision", matched["primary_precision"])
        _emit("isoform_matched_locus_recall", matched["isoform_recall"])
        _emit("isoform_matched_locus_f1", matched["isoform_f1"])
        _emit("isoform_junction_precision", matched["junction_precision"])
        _emit("isoform_junction_recall", matched["junction_recall"])
        _emit("isoform_junction_f1", matched["junction_f1"])
        _emit("isoform_matched_genes_n", matched["matched_genes_n"])
        _emit("isoform_matched_multiiso_genes_n", matched["matched_multiiso_genes_n"])
        _emit("isoform_split_merge_excluded_n", matched["split_merge_excluded_n"])
        # Raw denominators / numerators (traceability).
        _emit("isoform_matched_pred_alt_n", matched["pred_alt_isoforms_n"])
        _emit("isoform_matched_pred_alt_hit_n", matched["pred_alt_isoforms_hit_n"])
        _emit("isoform_matched_pred_primary_n", matched["pred_primary_n"])
        _emit("isoform_matched_pred_primary_hit_n", matched["pred_primary_hit_n"])
        _emit("isoform_matched_ref_iso_n", matched["ref_isoforms_n"])
        _emit("isoform_matched_ref_iso_hit_n", matched["ref_isoforms_hit_n"])
        _emit("isoform_junction_pred_n", matched["pred_introns_n"])
        _emit("isoform_junction_pred_hit_n", matched["pred_introns_hit_n"])
        _emit("isoform_junction_ref_n", matched["ref_introns_n"])
        _emit("isoform_junction_ref_hit_n", matched["ref_introns_hit_n"])

    for key, value in result["counts"].items():
        _emit(f"isoform_count_{key}", value)

    if "omark" in result:
        for key, value in result["omark"].items():
            _emit(f"isoform_omark_{key}", value)

    return rows
