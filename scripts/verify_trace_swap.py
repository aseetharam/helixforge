#!/usr/bin/env python3
"""Verify a TRaCE-on GFF3 before swapping it in as the benchmark GFF3.

Poster step 2c (``prompts/helixforge-poster/02c_regenerate_gff3_trace_on.md``).
Loads the OLD (no-TRaCE, combined_score primary) and NEW (TRaCE-on) HelixForge
GFF3s and asserts:

  1. Structure unchanged  — identical gene count, tier counts, origin tallies,
     biotype tallies, and the per-gene SET of transcript structures (exon+CDS
     coordinate tuples). TRaCE only reorders/renumbers a gene's transcripts, so
     the multiset of transcript structures per gene must be byte-for-byte the
     same; only the ``.N`` labels and the ``primary`` flag move. (This is the
     GFF3-level analogue of ``summarize_genes()``, which needs in-memory
     ``ReconciledGene`` objects.) AS-event tallies are compared from the
     report.tsv ``num_as_events`` column when both reports are present.
  2. IDs unique  — every gene and transcript id in the NEW GFF3 is globally
     unique.
  3. TRaCE took effect  — count genes whose PRIMARY (``primary=true``) transcript
     structure differs old vs new. Must be > 0 (else the flag silently no-oped).

Exit status is nonzero (and the offending assertion is printed) if any check
fails, so a caller can gate the swap on ``verify_trace_swap.py ... && cp ...``.

Usage:
    python scripts/verify_trace_swap.py OLD.gff3 NEW.gff3 \\
        [--old-report OLD.report.tsv --new-report NEW.report.tsv]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path


def _attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in field.rstrip(";").split(";"):
        if not kv:
            continue
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def parse_gff3(path: Path) -> dict:
    """Parse a HelixForge GFF3 into a structural view.

    Returns a dict with:
      genes:        {gene_id: {tier, origin, biotype}}
      gene_struct:  {gene_id: frozenset of transcript-structure keys}
      primary:      {gene_id: primary transcript-structure key}
      tx_ids:       list of all transcript ids (for uniqueness)
      gene_ids:     list of all gene ids (for uniqueness)

    A transcript-structure key is (strand, exon-tuple, cds-tuple) — the geometry
    that must be invariant under a pure reorder; the ``.N`` id is excluded on
    purpose so a relabelled-but-identical transcript compares equal.
    """
    genes: dict[str, dict] = {}
    gene_ids: list[str] = []
    tx_ids: list[str] = []
    # transcript id -> [strand, primary?, exons[], cds[], gene_id]
    tx: dict[str, dict] = {}

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype, start, end, strand, attr = f[2], f[3], f[4], f[6], f[8]
            a = _attrs(attr)
            if ftype == "gene":
                gid = a["ID"]
                gene_ids.append(gid)
                genes[gid] = {
                    "tier": a.get("tier", ""),
                    "origin": a.get("origin", ""),
                    "biotype": a.get("gene_biotype", "none"),
                }
            elif ftype in ("mRNA", "transcript"):
                tid = a["ID"]
                tx_ids.append(tid)
                tx[tid] = {
                    "strand": strand,
                    "primary": a.get("primary", "") == "true",
                    "gene_id": a["Parent"],
                    "exons": [],
                    "cds": [],
                }
            elif ftype == "exon":
                tx[a["Parent"]]["exons"].append((int(start), int(end)))
            elif ftype == "CDS":
                tx[a["Parent"]]["cds"].append((int(start), int(end)))

    gene_struct: dict[str, list] = {g: [] for g in genes}
    primary: dict[str, tuple] = {}
    for tid, t in tx.items():
        key = (
            t["strand"],
            tuple(sorted(t["exons"])),
            tuple(sorted(t["cds"])),
        )
        gene_struct[t["gene_id"]].append(key)
        if t["primary"]:
            primary[t["gene_id"]] = key

    return {
        "genes": genes,
        # sorted tuple = multiset of structures, order-independent
        "gene_struct": {g: tuple(sorted(v)) for g, v in gene_struct.items()},
        "primary": primary,
        "tx_ids": tx_ids,
        "gene_ids": gene_ids,
    }


def num_as_events(report: Path) -> dict[str, int]:
    rows: dict[str, int] = {}
    with open(report) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi = header.index("gene_id")
        ai = header.index("num_as_events")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rows[c[gi]] = int(c[ai])
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("old", help="old (no-TRaCE) GFF3")
    ap.add_argument("new", help="new (TRaCE-on) GFF3")
    ap.add_argument("--old-report", default=None)
    ap.add_argument("--new-report", default=None)
    args = ap.parse_args()

    old = parse_gff3(Path(args.old))
    new = parse_gff3(Path(args.new))

    failures: list[str] = []

    def check(cond: bool, msg: str) -> None:
        status = "OK  " if cond else "FAIL"
        print(f"  [{status}] {msg}")
        if not cond:
            failures.append(msg)

    print("== 1. Structure unchanged ==")
    check(len(old["genes"]) == len(new["genes"]),
          f"gene count: old={len(old['genes'])} new={len(new['genes'])}")
    check(set(old["genes"]) == set(new["genes"]),
          "gene id set identical")

    for fld in ("tier", "origin", "biotype"):
        oc = Counter(g[fld] for g in old["genes"].values())
        nc = Counter(g[fld] for g in new["genes"].values())
        check(oc == nc, f"{fld} tally: old={dict(sorted(oc.items()))} "
                        f"new={dict(sorted(nc.items()))}")

    # per-gene multiset of transcript structures must be identical
    shared = set(old["genes"]) & set(new["genes"])
    diff_struct = [g for g in shared
                   if old["gene_struct"][g] != new["gene_struct"][g]]
    check(not diff_struct,
          f"per-gene transcript-structure set identical "
          f"({len(diff_struct)} genes differ)")
    if diff_struct:
        for g in diff_struct[:5]:
            print(f"        e.g. {g}")

    if args.old_report and args.new_report:
        oas = num_as_events(Path(args.old_report))
        nas = num_as_events(Path(args.new_report))
        check(oas == nas, "report num_as_events per gene identical")
    else:
        print("  [skip] AS-event tally (reports not supplied)")

    print("== 2. IDs unique (new GFF3) ==")
    gdup = [k for k, v in Counter(new["gene_ids"]).items() if v > 1]
    tdup = [k for k, v in Counter(new["tx_ids"]).items() if v > 1]
    check(not gdup, f"gene ids unique ({len(gdup)} dups)")
    check(not tdup, f"transcript ids unique ({len(tdup)} dups)")

    print("== 3. TRaCE took effect ==")
    reordered = [g for g in shared
                 if g in old["primary"] and g in new["primary"]
                 and old["primary"][g] != new["primary"][g]]
    multi = sum(1 for g in new["genes"]
                if len(new["gene_struct"][g]) > 1)
    pct = 100 * len(reordered) / max(multi, 1)
    print(f"        {len(reordered)} genes elected a new canonical transcript "
          f"({len(reordered)}/{multi} multi-isoform genes = {pct:.1f}%)")
    check(len(reordered) > 0,
          "TRaCE changed at least one primary (flag reached run_pipeline)")

    print()
    if failures:
        print(f"VERIFY FAILED — {len(failures)} assertion(s) failed; DO NOT SWAP.")
        return 1
    print("VERIFY PASSED — safe to back up and swap.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
