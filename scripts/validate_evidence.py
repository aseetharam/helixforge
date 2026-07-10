#!/usr/bin/env python3
"""Invariant checker for ``helixforge evidence`` output TSVs.

Reads an evidence TSV (the per-transcript table from ``helixforge evidence`` or
its ``*.gene.tsv`` rollup) and asserts the numeric invariants the scorer is
supposed to guarantee, reporting *every* violation with the offending
``transcript_id`` rather than stopping at the first. It also audits the
suspected-duplicate columns and the summary-table mean denominators, and prints
the value distributions that make the discrete/degenerate cases visible.

The AED weights and the single-exon renormalisation are imported from
``helixforge.score.evidence`` (never hardcoded), so this checker tracks the code.

Usage::

    python scripts/validate_evidence.py path/to/evidence.tsv
    python scripts/validate_evidence.py evidence.tsv --tol 1e-6

Exit code is the number of violation *categories* that fired (0 = clean), so it
is usable as a CI gate.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

# Pull the live weights + the renormalising AED helper from the scorer so this
# validator can never drift from the formulas it is checking.
from helixforge.score.evidence import (
    PROT_W_CDS,
    PROT_W_PROT,
    PROT_W_STRUCT,
    RNA_W_BOUNDARY,
    RNA_W_COVERAGE,
    RNA_W_JUNCTION,
    protein_aed_from_ratios,
    rna_aed_from_ratios,
)

# Columns that must lie in [0, 1] when present.
_RATIO_COLUMNS = (
    "junction_support_fraction",
    "intron_precision",
    "intron_recall",
    "intron_f1",
    "rna_aed",
    "rna_junction_ratio",
    "rna_coverage_ratio",
    "rna_boundary_ratio",
    "protein_aed",
    "protein_struct_ratio",
    "protein_cds_cov_ratio",
    "protein_prot_cov_ratio",
)

# Intron-chain metrics that are blank *iff* the model has no introns.
_STRICT_INTRON_COLUMNS = ("junction_support_fraction", "intron_precision")
# recall/f1 may *also* be blank when introns exist but no junction overlaps the
# locus (no qualifying evidence), so they get a relaxed rule (see _check_blanks).
_RELAXED_INTRON_COLUMNS = ("intron_recall", "intron_f1")


def _is_blank(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, float):
        return math.isnan(v)
    return str(v).strip() == ""


def _f(v: Any) -> float:
    return float(v)


class Report:
    """Accumulates violations per rule and the informational findings."""

    def __init__(self) -> None:
        self.violations: dict[str, list[str]] = {}
        self.findings: list[str] = []

    def violate(self, rule: str, tid: str, detail: str = "") -> None:
        msg = tid if not detail else f"{tid}: {detail}"
        self.violations.setdefault(rule, []).append(msg)

    def n_categories(self) -> int:
        return sum(1 for v in self.violations.values() if v)


# ---------------------------------------------------------------------------
# Rule 1 — AED formulas reconstruct from their components
# ---------------------------------------------------------------------------


def _check_aed_reconstruction(rows: list[dict[str, Any]], tol: float, rep: Report) -> None:
    for r in rows:
        tid = str(r["transcript_id"])

        if not _is_blank(r.get("rna_aed")):
            jr = r.get("rna_junction_ratio")
            cr = r.get("rna_coverage_ratio")
            br = r.get("rna_boundary_ratio")
            if any(_is_blank(x) for x in (jr, cr, br)):
                rep.violate(
                    "rna_aed_components_missing",
                    tid,
                    "rna_aed present but a component ratio is blank",
                )
            else:
                expect = rna_aed_from_ratios(
                    _f(jr), _f(cr), _f(br),
                    weights=(RNA_W_JUNCTION, RNA_W_COVERAGE, RNA_W_BOUNDARY),
                )
                if abs(expect - _f(r["rna_aed"])) > tol:
                    rep.violate(
                        "rna_aed_reconstruction", tid,
                        f"stored={_f(r['rna_aed']):.6f} expected={expect:.6f}",
                    )

        if not _is_blank(r.get("protein_aed")):
            cc = r.get("protein_cds_cov_ratio")
            pc = r.get("protein_prot_cov_ratio")
            sr = r.get("protein_struct_ratio")
            if any(_is_blank(x) for x in (cc, pc)):
                rep.violate(
                    "protein_aed_components_missing", tid,
                    "protein_aed present but cds/prot ratio is blank",
                )
            else:
                # struct_ratio blank -> single-exon model: AED renormalises over
                # the two present components (matches the scorer).
                struct = None if _is_blank(sr) else _f(sr)
                expect = protein_aed_from_ratios(
                    struct, _f(cc), _f(pc),
                    weights=(PROT_W_STRUCT, PROT_W_CDS, PROT_W_PROT),
                )
                if abs(expect - _f(r["protein_aed"])) > tol:
                    rep.violate(
                        "protein_aed_reconstruction", tid,
                        f"stored={_f(r['protein_aed']):.6f} expected={expect:.6f}",
                    )


# ---------------------------------------------------------------------------
# Rule 2 — every ratio in [0, 1]
# ---------------------------------------------------------------------------


def _check_ratio_ranges(rows: list[dict[str, Any]], rep: Report) -> None:
    for r in rows:
        tid = str(r["transcript_id"])
        for col in _RATIO_COLUMNS:
            v = r.get(col)
            if _is_blank(v):
                continue
            fv = _f(v)
            if fv < 0.0 or fv > 1.0:
                rep.violate("ratio_out_of_range", tid, f"{col}={fv}")


# ---------------------------------------------------------------------------
# Rule 3 — suspected duplicate columns
# ---------------------------------------------------------------------------


def _check_duplicate_columns(rows: list[dict[str, Any]], tol: float, rep: Report) -> None:
    """precision vs junction_support_fraction vs rna_junction_ratio.

    precision and junction_support_fraction are the *same* formula
    (supported/num_introns); rna_junction_ratio equals them for multi-exon models
    and defaults to 1.0 for single-exon. We report whether they coincide wherever
    all are defined and list any genuine divergence.
    """
    pj_compared = pj_diverge = 0
    pjr_compared = pjr_diverge = 0
    for r in rows:
        tid = str(r["transcript_id"])
        prec, jsf, jr = (
            r.get("intron_precision"),
            r.get("junction_support_fraction"),
            r.get("rna_junction_ratio"),
        )
        if not _is_blank(prec) and not _is_blank(jsf):
            pj_compared += 1
            if abs(_f(prec) - _f(jsf)) > tol:
                pj_diverge += 1
                rep.violate(
                    "duplicate_precision_jsf_diverge", tid,
                    f"precision={_f(prec)} jsf={_f(jsf)}",
                )
        if not _is_blank(prec) and not _is_blank(jr):
            pjr_compared += 1
            if abs(_f(prec) - _f(jr)) > tol:
                pjr_diverge += 1

    if pj_compared:
        verdict = "ALWAYS equal" if pj_diverge == 0 else f"{pj_diverge} divergent"
        rep.findings.append(
            f"intron_precision vs junction_support_fraction: {verdict} "
            f"over {pj_compared} rows where both defined "
            f"(definitionally identical: both = supported/num_introns)."
        )
    if pjr_compared:
        verdict = "ALWAYS equal" if pjr_diverge == 0 else f"{pjr_diverge} divergent"
        rep.findings.append(
            f"intron_precision vs rna_junction_ratio: {verdict} over "
            f"{pjr_compared} multi-exon rows (single-exon rna_junction_ratio "
            f"defaults to 1.0 where precision is blank — divergence by design)."
        )


# ---------------------------------------------------------------------------
# Rule 4 — blank-cell discipline
# ---------------------------------------------------------------------------


def _check_blanks(rows: list[dict[str, Any]], rep: Report) -> None:
    relaxed_blank_with_introns = 0
    for r in rows:
        tid = str(r["transcript_id"])
        ni = r.get("num_introns")
        has_introns = (not _is_blank(ni)) and int(ni) > 0
        supported = r.get("supported")
        sup = int(supported) if not _is_blank(supported) else 0

        # Strict intron metrics: blank iff no introns.
        for col in _STRICT_INTRON_COLUMNS:
            blank = _is_blank(r.get(col))
            if has_introns and blank:
                rep.violate("intron_metric_blank_with_introns", tid, col)
            if not has_introns and not blank:
                rep.violate("intron_metric_populated_no_introns", tid, col)

        # Relaxed intron metrics: must be blank with no introns; may be blank with
        # introns only when there is no junction evidence (supported == 0). If
        # supported > 0 a qualifying junction exists, so recall/f1 must be defined.
        for col in _RELAXED_INTRON_COLUMNS:
            blank = _is_blank(r.get(col))
            if not has_introns and not blank:
                rep.violate("intron_metric_populated_no_introns", tid, col)
            if has_introns and blank:
                if sup > 0:
                    rep.violate("relaxed_intron_blank_with_support", tid, col)
                else:
                    relaxed_blank_with_introns += 1

        # Protein columns: blank iff no protein_id, EXCEPT protein_struct_ratio is
        # also (legitimately) blank for a single-exon coding model.
        no_protein = _is_blank(r.get("protein_id"))
        for col in ("protein_aed", "protein_cds_cov_ratio", "protein_prot_cov_ratio"):
            blank = _is_blank(r.get(col))
            if no_protein and not blank:
                rep.violate("protein_metric_populated_no_hit", tid, col)
            if not no_protein and blank:
                rep.violate("protein_metric_blank_with_hit", tid, col)
        struct_blank = _is_blank(r.get("protein_struct_ratio"))
        if no_protein and not struct_blank:
            rep.violate("protein_metric_populated_no_hit", tid, "protein_struct_ratio")
        if not no_protein and not has_introns and not struct_blank:
            rep.violate(
                "struct_ratio_populated_no_introns", tid,
                "protein_struct_ratio should be blank for a no-intron model",
            )

    if relaxed_blank_with_introns:
        rep.findings.append(
            f"{relaxed_blank_with_introns} multi-exon rows have blank "
            f"intron_recall/intron_f1 with supported==0 — legitimate: no "
            f"qualifying junction overlaps the locus (recall denominator is 0)."
        )


# ---------------------------------------------------------------------------
# Rule 5 — count consistency
# ---------------------------------------------------------------------------


def _check_counts(rows: list[dict[str, Any]], tol: float, rep: Report) -> None:
    for r in rows:
        tid = str(r["transcript_id"])
        ni = r.get("num_introns")
        if _is_blank(ni):
            continue
        ni = int(ni)
        sup = int(r["supported"]) if not _is_blank(r.get("supported")) else 0
        con = int(r["contradicted"]) if not _is_blank(r.get("contradicted")) else 0
        nov = int(r["novel_in_data"]) if not _is_blank(r.get("novel_in_data")) else 0
        ne = int(r["num_exons"]) if not _is_blank(r.get("num_exons")) else None

        if sup + con + nov != ni:
            rep.violate(
                "count_partition", tid,
                f"supported+contradicted+novel={sup + con + nov} != num_introns={ni}",
            )
        if sup > ni:
            rep.violate("supported_exceeds_introns", tid, f"{sup} > {ni}")
        if ne is not None and ne != ni + 1:
            rep.violate("exon_intron_count", tid, f"num_exons={ne} num_introns={ni}")

        jsf = r.get("junction_support_fraction")
        if ni > 0 and not _is_blank(jsf):
            if abs(_f(jsf) - sup / ni) > tol:
                rep.violate(
                    "jsf_formula", tid, f"jsf={_f(jsf)} != supported/num={sup / ni}",
                )


# ---------------------------------------------------------------------------
# Rule 6 — intron_f1 == harmonic_mean(precision, recall)
# ---------------------------------------------------------------------------


def _check_f1(rows: list[dict[str, Any]], tol: float, rep: Report) -> None:
    for r in rows:
        tid = str(r["transcript_id"])
        p, rc, f1 = (
            r.get("intron_precision"),
            r.get("intron_recall"),
            r.get("intron_f1"),
        )
        if _is_blank(p) or _is_blank(rc):
            if not _is_blank(f1):
                rep.violate("f1_present_without_components", tid)
            continue
        p, rc = _f(p), _f(rc)
        expect = 0.0 if p + rc == 0 else 2 * p * rc / (p + rc)
        if _is_blank(f1):
            rep.violate("f1_missing", tid)
        elif abs(_f(f1) - expect) > tol:
            rep.violate("f1_reconstruction", tid, f"f1={_f(f1)} expected={expect}")


# ---------------------------------------------------------------------------
# Distributions + summary-mean denominator audit
# ---------------------------------------------------------------------------


def _distribution(rows: list[dict[str, Any]], col: str) -> Counter:
    c: Counter = Counter()
    for r in rows:
        v = r.get(col)
        c["<blank>" if _is_blank(v) else f"{_f(v):g}"] += 1
    return c


def _summary_denominators(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Audit the means a summary table would report: their denominators differ."""
    rna = [_f(r["rna_aed"]) for r in rows if not _is_blank(r.get("rna_aed"))]
    prot = [_f(r["protein_aed"]) for r in rows if not _is_blank(r.get("protein_aed"))]
    return {
        "n_all": len(rows),
        "n_rna_aed": len(rna),
        "mean_rna_aed": (sum(rna) / len(rna)) if rna else None,
        "n_protein_aed": len(prot),
        "mean_protein_aed": (sum(prot) / len(prot)) if prot else None,
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _read_tsv(path: str | Path) -> list[dict[str, Any]]:
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    return df.to_dict("records")


def validate(rows: list[dict[str, Any]], tol: float = 1e-6) -> Report:
    """Run every rule over already-parsed rows; return the populated Report."""
    rep = Report()
    _check_aed_reconstruction(rows, tol, rep)
    _check_ratio_ranges(rows, rep)
    _check_duplicate_columns(rows, tol, rep)
    _check_blanks(rows, rep)
    _check_counts(rows, tol, rep)
    _check_f1(rows, tol, rep)
    return rep


def _print_report(rows: list[dict[str, Any]], rep: Report, max_examples: int) -> None:
    print(f"rows checked: {len(rows)}")
    print()
    print("== violations per rule ==")
    if not rep.violations or rep.n_categories() == 0:
        print("  (none)")
    for rule, items in rep.violations.items():
        if not items:
            continue
        ex = ", ".join(items[:max_examples])
        more = "" if len(items) <= max_examples else f"  (+{len(items) - max_examples} more)"
        print(f"  {rule}: {len(items)}")
        print(f"      e.g. {ex}{more}")
    print()
    print("== findings (decide on; not failures) ==")
    for f in rep.findings:
        print(f"  - {f}")
    print()

    print("== distributions (discrete/degenerate cases) ==")
    for col in ("rna_boundary_ratio", "protein_struct_ratio"):
        dist = _distribution(rows, col)
        rendered = ", ".join(f"{k}:{v}" for k, v in sorted(dist.items()))
        print(f"  {col}: {rendered}")
    print()

    sd = _summary_denominators(rows)
    print("== summary-mean denominator audit ==")
    rna_m = "—" if sd["mean_rna_aed"] is None else f"{sd['mean_rna_aed']:.4f}"
    prot_m = "—" if sd["mean_protein_aed"] is None else f"{sd['mean_protein_aed']:.4f}"
    print(f"  mean_rna_aed     = {rna_m}  (n={sd['n_rna_aed']} of {sd['n_all']}, all scored)")
    print(f"  mean_protein_aed = {prot_m}  (n={sd['n_protein_aed']} with protein hit)")
    if sd["n_rna_aed"] != sd["n_protein_aed"]:
        print(
            "  NOTE: the two means use DIFFERENT denominators "
            "(RNA over all transcripts, protein over the hit subset) — "
            "do not compare them directly."
        )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=Path, help="evidence TSV to validate")
    ap.add_argument("--tol", type=float, default=1e-6, help="float tolerance")
    ap.add_argument(
        "--max-examples", type=int, default=10,
        help="offending transcript_ids to print per rule",
    )
    args = ap.parse_args(argv)

    rows = _read_tsv(args.tsv)
    rep = validate(rows, tol=args.tol)
    _print_report(rows, rep, args.max_examples)

    n = rep.n_categories()
    print()
    print(f"RESULT: {'PASS' if n == 0 else 'FAIL'} ({n} violation categories)")
    return n


if __name__ == "__main__":
    sys.exit(main())
