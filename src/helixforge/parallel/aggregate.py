"""Gather per-chunk outputs into one genome-wide annotation."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

from helixforge.utils.logging import get_logger
from helixforge.utils.regions import gff3_to_internal

_log = get_logger(__name__)

_TIERS = (1, 2, 3)


@dataclass
class AggregateResult:
    """Outcome of :func:`aggregate`, merged paths + genome-level tallies."""

    gff3_path: Path
    tier_paths: dict[int, Path]
    report_path: Path
    master_id_map_path: Path | None
    num_genes: int
    num_loci: int
    tier_counts: dict[str, int] = field(default_factory=dict)
    origin_counts: dict[str, int] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# GFF3 block handling
# ---------------------------------------------------------------------------


def _gene_id_from_attrs(attr_field: str) -> str | None:
    for part in attr_field.split(";"):
        part = part.strip()
        if part.startswith("ID="):
            return part[3:]
    return None


def _read_gene_blocks(path: Path | str) -> list[tuple[tuple[str, int], str, str]]:
    """Parse a HelixForge GFF3 into ``[(sort_key, gene_id, block_text)]``.

    A block is the lines of one gene record (gene→mRNA→exon/CDS) up to the ``###``
    separator the writer emits. ``sort_key = (seqid, internal_start)`` is taken
    from the block's ``gene`` line; ``gene_id`` from its ``ID=``.
    """
    blocks: list[tuple[tuple[str, int], str, str]] = []
    current: list[str] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.rstrip("\n")
            if stripped == "###":  # record separator, check BEFORE "##" skip
                if current:
                    blocks.append(_finalize_block(current, path))
                    current = []
                continue
            if line.startswith("#"):  # directives (##gff-version,
                # ##sequence-region, ##FASTA) + comments (#! provenance). The
                # ### record separator is already handled above.
                continue
            if not line.strip():
                continue
            current.append(line if line.endswith("\n") else line + "\n")
    if current:
        blocks.append(_finalize_block(current, path))
    return blocks


def _finalize_block(
    lines: list[str], path: Path | str
) -> tuple[tuple[str, int], str, str]:
    seqid: str | None = None
    start: int | None = None
    gene_id: str | None = None
    for line in lines:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        if cols[2] == "gene":
            seqid = cols[0]
            start = gff3_to_internal(int(cols[3]), int(cols[4]))[0]
            gene_id = _gene_id_from_attrs(cols[8])
            break
    if seqid is None:  # no gene line, fall back to the first feature
        cols = lines[0].rstrip("\n").split("\t")
        seqid = cols[0]
        start = gff3_to_internal(int(cols[3]), int(cols[4]))[0]
        gene_id = _gene_id_from_attrs(cols[8]) or f"{path}:{lines[0][:20]}"
    # After the fallback branch above, all three are always set:
    # either the loop set seqid/start/gene_id, or the fallback branch did.
    return ((seqid, cast(int, start)), cast(str, gene_id), "".join(lines))


def _write_merged_gff3(
    blocks: list[tuple[tuple[str, int], str, str]],
    out_path: Path | str,
) -> None:
    """Write blocks (already sorted) as one GFF3 with the writer's framing."""
    with open(out_path, "w") as fh:
        fh.write("##gff-version 3\n")
        for _key, _gid, text in blocks:
            fh.write(text)
            fh.write("###\n")


def _merge_gff3_files(paths: list[str], out_path: Path | str) -> list[str]:
    """Concatenate + re-sort the gene blocks of ``paths`` → ``out_path``.

    Returns the list of gene ids in genomic order. Raises on a duplicate gene id
    (global HFG-uniqueness violation).
    """
    blocks: list[tuple[tuple[str, int], str, str]] = []
    for path in paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"chunk GFF3 not found: {path}")
        blocks.extend(_read_gene_blocks(path))
    blocks.sort(key=lambda b: b[0])

    seen: dict[str, Path | str] = {}
    seen_tx: set[str] = set()
    for _key, gene_id, text in blocks:
        if gene_id in seen:
            raise ValueError(
                f"duplicate gene id across chunks: {gene_id!r} appears in "
                f"{seen[gene_id]} and again while merging, HFG ranges overlapped"
            )
        seen[gene_id] = out_path
        # Global transcript-id uniqueness: two chunks must never emit the same
        # isoform id (would mean reserved HFG ranges overlapped).
        for tid in _transcript_ids(text):
            if tid in seen_tx:
                raise ValueError(
                    f"duplicate transcript id across chunks: {tid!r}, "
                    "HFG ranges overlapped"
                )
            seen_tx.add(tid)
    _write_merged_gff3(blocks, out_path)
    return [gid for (_k, gid, _t) in blocks]


def _transcript_ids(block_text: str) -> list[str]:
    """Extract the mRNA ``ID=`` values from one gene block's lines."""
    tids: list[str] = []
    for line in block_text.splitlines():
        cols = line.split("\t")
        if len(cols) >= 9 and cols[2] == "mRNA":
            tid = _gene_id_from_attrs(cols[8])
            if tid is not None:
                tids.append(tid)
    return tids


def prefixes_from_pattern(input_dir: str | Path, pattern: str) -> list[str]:
    """Resolve a v1-style ``--input-dir`` + ``--pattern`` glob to chunk prefixes.

    Matches ``<input_dir>/<pattern>`` (expected to hit the per-chunk GFF3s), then
    derives each chunk's output **prefix** by stripping the ``.gff3`` suffix and an
    optional ``.tier{1,2,3}`` segment, so ``*.gff3`` matching both
    ``chunk_0000.gff3`` and ``chunk_0000.tier1.gff3`` collapses to the single
    prefix ``chunk_0000``. Returns sorted, de-duplicated prefixes (full paths).
    """
    matches = sorted(Path(input_dir).glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"no files match {pattern!r} in {input_dir}: point --input-dir at the "
            "per-chunk outputs and --pattern at their .gff3 files"
        )
    prefixes: list[str] = []
    seen: set[str] = set()
    for path in matches:
        name = path.name
        if not name.endswith(".gff3"):
            continue
        stem = name[: -len(".gff3")]
        for tier in ("tier1", "tier2", "tier3"):
            if stem.endswith("." + tier):
                stem = stem[: -(len(tier) + 1)]
                break
        prefix = str(path.parent / stem)
        if prefix not in seen:
            seen.add(prefix)
            prefixes.append(prefix)
    if not prefixes:
        raise FileNotFoundError(
            f"pattern {pattern!r} matched files in {input_dir} but none ended in "
            ".gff3, aggregate needs the per-chunk .gff3 outputs"
        )
    return sorted(prefixes)


# ---------------------------------------------------------------------------
# Report + id_map handling
# ---------------------------------------------------------------------------


def _merge_reports(
    paths: list[str],
    out_path: Path | str,
) -> tuple[dict[str, int], dict[str, int], int]:
    """Union per-chunk ``report.tsv`` rows, re-sorted by (seqid, start).

    Returns ``(tier_counts, origin_counts, num_rows)`` summed over all rows.
    """
    header: str | None = None
    rows: list[list[str]] = []
    for path in paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"chunk report not found: {path}")
        with open(path) as fh:
            file_header = fh.readline().rstrip("\n")
            if header is None:
                header = file_header
            elif file_header != header:
                raise ValueError(f"report header mismatch in {path}")
            for line in fh:
                line = line.rstrip("\n")
                if line:
                    rows.append(line.split("\t"))

    cols = header.split("\t") if header else []
    seqid_i = cols.index("seqid") if "seqid" in cols else 1
    start_i = cols.index("start") if "start" in cols else 2
    tier_i = cols.index("tier") if "tier" in cols else None
    origin_i = cols.index("origin") if "origin" in cols else None

    rows.sort(key=lambda r: (r[seqid_i], int(r[start_i])))
    tier_counts = Counter(r[tier_i] for r in rows) if tier_i is not None else Counter()
    origin_counts = (
        Counter(r[origin_i] for r in rows) if origin_i is not None else Counter()
    )

    with open(out_path, "w") as fh:
        if header is not None:
            fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return (
        dict(sorted(tier_counts.items())),
        dict(sorted(origin_counts.items())),
        len(rows),
    )


def _union_id_maps(paths: list[str]) -> dict[str, str]:
    """Load + union chunk id_maps; verify disjoint loci + globally-unique HFGs.

    Returns the merged ``{helixer_locus_id: HFG}``. Raises if a locus is owned by
    two chunks (gene split across chunks) or an HFG number is reused in two
    chunks (range overlap).
    """
    merged: dict[str, str] = {}
    key_owner: dict[str, str] = {}  # locus id -> chunk path that first claimed it
    value_owner: dict[str, str] = {}  # HFG -> chunk path that first allocated it
    for path in paths:
        if not Path(path).exists():
            raise FileNotFoundError(f"chunk id_map not found: {path}")
        chunk_map = json.loads(Path(path).read_text())
        for locus_id, hfg in chunk_map.items():
            if locus_id in key_owner and key_owner[locus_id] != str(path):
                raise ValueError(
                    f"Helixer locus {locus_id!r} covered by >1 chunk "
                    f"({key_owner[locus_id]} and {path}), partition split a gene"
                )
            if hfg in value_owner and value_owner[hfg] != str(path):
                raise ValueError(
                    f"HFG {hfg!r} allocated in two chunks "
                    f"({value_owner[hfg]} and {path}), reserved ranges overlapped"
                )
            key_owner[locus_id] = str(path)
            value_owner[hfg] = str(path)
            merged[locus_id] = hfg
    return merged


def _fold_into_master(
    merged: dict[str, str],
    master_id_map_path: str | Path | None,
) -> dict[str, str]:
    """Union ``merged`` into the persisted master, asserting no remapping."""
    master: dict[str, str] = {}
    if master_id_map_path and Path(master_id_map_path).exists():
        master = json.loads(Path(master_id_map_path).read_text())
    for locus_id, hfg in merged.items():
        if locus_id in master and master[locus_id] != hfg:
            raise ValueError(
                f"id_map conflict for {locus_id!r}: master has {master[locus_id]!r} "
                f"but chunk assigned {hfg!r}: id stability violated"
            )
        master[locus_id] = hfg
    if master_id_map_path:
        Path(master_id_map_path).write_text(
            json.dumps(master, indent=0, sort_keys=True)
        )
    return master


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def aggregate(
    chunk_outputs: list[Any],
    out_prefix: str | Path,
    master_id_map_path: str | Path | None = None,
) -> AggregateResult:
    """Merge per-chunk outputs into one genome-wide annotation + verify it.

    ``chunk_outputs`` is the list of per-chunk ``output_prefix`` values (as
    produced by :func:`~helixforge.parallel.tasks.run_local` or the drivers).
    For each prefix the files ``<prefix>.gff3``, ``<prefix>.tier{1,2,3}.gff3``,
    ``<prefix>.report.tsv`` and ``<prefix>.id_map.json`` are expected. Writes the
    merged ``<out_prefix>.*`` set, folds chunk id_maps into ``master_id_map_path``
    (if given), and returns an :class:`AggregateResult`. Raises on any HFG
    uniqueness / locus-coverage / id-stability violation.
    """
    prefixes = [str(p) for p in chunk_outputs]
    if not prefixes:
        raise ValueError("no chunk outputs to aggregate")

    # --- merge the full GFF3 (and check global gene-id uniqueness) ---
    gff3_path = Path(f"{out_prefix}.gff3")
    gene_ids = _merge_gff3_files([f"{p}.gff3" for p in prefixes], gff3_path)

    # --- merge tier GFF3s ---
    tier_paths: dict[int, Path] = {}
    for n in _TIERS:
        tp = Path(f"{out_prefix}.tier{n}.gff3")
        _merge_gff3_files([f"{p}.tier{n}.gff3" for p in prefixes], tp)
        tier_paths[n] = tp

    # --- union reports ---
    report_path = Path(f"{out_prefix}.report.tsv")
    tier_counts, origin_counts, num_rows = _merge_reports(
        [f"{p}.report.tsv" for p in prefixes], report_path
    )

    # --- union + verify id_maps, fold into master ---
    merged_map = _union_id_maps([f"{p}.id_map.json" for p in prefixes])
    _fold_into_master(merged_map, master_id_map_path)

    if num_rows != len(gene_ids):
        raise ValueError(
            f"report rows ({num_rows}) != merged genes ({len(gene_ids)}), "
            "a chunk's GFF3 and report disagree"
        )

    result = AggregateResult(
        gff3_path=gff3_path,
        tier_paths=tier_paths,
        report_path=report_path,
        master_id_map_path=Path(master_id_map_path) if master_id_map_path else None,
        num_genes=len(gene_ids),
        num_loci=len(merged_map),
        tier_counts=tier_counts,
        origin_counts=origin_counts,
    )
    _log.info(
        "aggregated %d chunks → %d genes, %d loci; tier=%s origin=%s",
        len(prefixes),
        result.num_genes,
        result.num_loci,
        result.tier_counts,
        result.origin_counts,
    )
    return result
