"""Functional-annotation hook: InterProScan / eggNOG-mapper."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import attrs

from helixforge.prep._subprocess import output_is_fresh, run_tool
from helixforge.qc.flags import DOMAIN_COMPLETE, dedup_flags
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene

_log = get_logger(__name__)

# Functional-annotation tools this hook understands.
FUNCTION_TOOLS = ("interproscan", "eggnog", "both")

# Biotypes whose proteins are worth functionally annotating: a real coding gene,
# or a homology-backed disabled homolog (pseudogene) whose residual domains are
# still informative. lncRNA / ncRNA / structured-ncRNA have no ORF to annotate.
ANNOTATABLE_BIOTYPES = ("protein_coding", "pseudogene")


@attrs.frozen
class FunctionalRecord:
    """Parsed functional annotation for one transcript (isoform).

    ``go_terms`` / ``dbxrefs`` are de-duplicated, sorted tuples so the emitted
    GFF3 is deterministic (diffable runs). ``domain_complete`` is
    the per-isoform ORF-credibility signal: True when at least one
    recognized domain match lies fully inside the translated protein.
    """

    transcript_id: str
    go_terms: tuple[str, ...] = ()
    dbxrefs: tuple[str, ...] = ()
    domain_complete: bool = False


def _norm(
    record_id: str, go: Iterable[str], xref: Iterable[str], domain_complete: bool
) -> FunctionalRecord:
    """Build a :class:`FunctionalRecord` with sorted/de-duplicated term lists."""
    return FunctionalRecord(
        transcript_id=record_id,
        go_terms=tuple(sorted({g for g in go if g})),
        dbxrefs=tuple(sorted({x for x in xref if x})),
        domain_complete=bool(domain_complete),
    )


# ---------------------------------------------------------------------------
# Parsers — tolerant, take a path, read tool output
# ---------------------------------------------------------------------------


def parse_interproscan_tsv(path: str | Path) -> dict[str, FunctionalRecord]:
    """Parse an InterProScan TSV → ``{transcript_id: FunctionalRecord}``.

    InterProScan's ``--formats TSV`` is a headerless, tab-separated file with one
    row per signature match. The columns used here (0-based):

    ``0`` protein accession (the FASTA header = our transcript id), ``2`` sequence
    length, ``3`` analysis (``Pfam`` / ``SUPERFAMILY`` / …), ``4`` signature
    accession, ``6`` match start, ``7`` match stop (both **protein** coordinates,
    1-based inclusive), ``11`` InterPro accession (optional), ``13`` GO terms
    (optional, ``|``-separated). Missing optional fields are ``-`` or absent.

    ``domain_complete`` is set when any recognized domain's match span lies fully
    within the protein (``1 <= start <= stop <= seq_len``) — the §2.2 signal that
    the ORF brackets a complete domain rather than a truncated fragment.
    """
    by_tx: dict[str, dict[str, Any]] = {}
    for raw in Path(path).read_text().splitlines():
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        tid = cols[0]
        analysis = cols[3]
        sig_acc = cols[4]
        try:
            seq_len = int(cols[2])
            start = int(cols[6])
            stop = int(cols[7])
        except ValueError:
            seq_len = start = stop = 0
        interpro_acc = cols[11] if len(cols) > 11 else ""
        go_field = cols[13] if len(cols) > 13 else ""

        acc = by_tx.setdefault(tid, {"go": set(), "xref": set(), "complete": False})
        if analysis.lower() == "pfam" and sig_acc and sig_acc != "-":
            acc["xref"].add(f"Pfam:{sig_acc}")
        if interpro_acc and interpro_acc.upper().startswith("IPR"):
            acc["xref"].add(f"InterPro:{interpro_acc}")
        for g in go_field.replace(",", "|").split("|"):
            g = g.strip()
            if g.startswith("GO:"):
                acc["go"].add(g)
        # A domain match fully bracketed by the ORF (recognized signature only).
        if sig_acc and sig_acc != "-" and seq_len > 0 and 1 <= start <= stop <= seq_len:
            acc["complete"] = True

    return {
        tid: _norm(tid, d["go"], d["xref"], d["complete"]) for tid, d in by_tx.items()
    }


def parse_eggnog_annotations(path: str | Path) -> dict[str, FunctionalRecord]:
    """Parse an eggNOG-mapper ``.emapper.annotations`` TSV → records.

    The file is tab-separated with ``#``-prefixed comment lines; the column header
    line begins ``#query``. The ``GOs`` column → ``go_terms``; the ``PFAMs``
    column → ``Pfam:<name>`` ``dbxrefs``; the ``seed_ortholog`` (a UniProt/RefSeq
    accession) → a ``UniProt:`` Dbxref when present. eggNOG carries no per-domain
    coordinates, so ``domain_complete`` is True when the row reports at least one
    Pfam family (a recognized domain) — a coarser signal than InterProScan's
    coordinate-aware one, documented as such.
    """
    header: list[str] | None = None
    by_tx: dict[str, dict[str, Any]] = {}
    for raw in Path(path).read_text().splitlines():
        line = raw.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            if line.lower().startswith("#query"):
                header = line.lstrip("#").split("\t")
            continue
        if header is None:
            continue
        cols = line.split("\t")
        row = {header[i]: cols[i] for i in range(min(len(header), len(cols)))}
        tid = row.get("query", cols[0]).strip()
        if not tid:
            continue

        def _vals(field: str) -> list[str]:
            v = (row.get(field, "") or "").strip()
            return (
                [] if v in ("", "-") else [t.strip() for t in v.split(",") if t.strip()]
            )

        gos = [g for g in _vals("GOs") if g.startswith("GO:")]
        pfams = _vals("PFAMs")
        xref = {f"Pfam:{p}" for p in pfams}
        seed = (row.get("seed_ortholog", "") or "").strip()
        if seed and seed != "-":
            xref.add(f"UniProt:{seed}")

        acc = by_tx.setdefault(tid, {"go": set(), "xref": set(), "complete": False})
        acc["go"].update(gos)
        acc["xref"].update(xref)
        if pfams:
            acc["complete"] = True

    return {
        tid: _norm(tid, d["go"], d["xref"], d["complete"]) for tid, d in by_tx.items()
    }


# ---------------------------------------------------------------------------
# Run wrappers — build argv, invoke, return the output path
# ---------------------------------------------------------------------------


def run_interproscan(
    proteins_fa: str | Path,
    out_path: str | Path,
    *,
    interproscan_bin: str = "interproscan.sh",
    applications: Iterable[str] | None = None,
    threads: int = 1,
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
) -> Path:
    """``interproscan.sh -i <proteins> -f TSV -o <out> [-appl ...] -cpu <n>``.

    Writes a TSV (``-f TSV``). ``applications`` restricts the member databases
    (e.g. ``["Pfam"]``). Skips when ``out_path`` is already fresh unless ``force``.
    Returns the output path.
    """
    out_path = Path(out_path)
    if not force and output_is_fresh(out_path, [proteins_fa]):
        _log.info("skip InterProScan: %s up to date", out_path.name)
        return out_path
    argv: list[str | Path | int] = [
        interproscan_bin,
        "-i",
        proteins_fa,
        "-f",
        "TSV",
        "-o",
        out_path,
        "-cpu",
        threads,
        # proteins are amino-acid sequences; disable nucleotide ORF prediction.
        "-dp",
    ]
    if applications:
        argv += ["-appl", ",".join(applications)]
    if extra_args:
        argv += list(extra_args)
    run_tool(argv)
    return out_path


def run_eggnog_mapper(
    proteins_fa: str | Path,
    out_dir: str | Path,
    *,
    out_prefix: str = "helixforge",
    eggnog_bin: str = "emapper.py",
    db: str | Path | None = None,
    threads: int = 1,
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
) -> Path:
    """``emapper.py -i <proteins> --itype proteins -o <prefix> --output_dir <dir>``.

    Returns the path to the ``<prefix>.emapper.annotations`` table. ``db`` points
    at the eggNOG data directory (``--data_dir``) when given. Skips when the
    annotations file is already fresh unless ``force``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    annotations = out_dir / f"{out_prefix}.emapper.annotations"
    if not force and output_is_fresh(annotations, [proteins_fa]):
        _log.info("skip eggNOG-mapper: %s up to date", annotations.name)
        return annotations
    argv: list[str | Path | int] = [
        eggnog_bin,
        "-i",
        proteins_fa,
        "--itype",
        "proteins",
        "-o",
        out_prefix,
        "--output_dir",
        out_dir,
        "--cpu",
        threads,
    ]
    if db is not None:
        argv += ["--data_dir", db]
    if extra_args:
        argv += list(extra_args)
    run_tool(argv)
    return annotations


# ---------------------------------------------------------------------------
# annotate_function — the opt-in orchestrator
# ---------------------------------------------------------------------------


def annotate_function(
    genes: list[ReconciledGene],
    proteins_fa: str | Path,
    *,
    out_dir: str | Path,
    tool: str = "interproscan",
    db: str | Path | None = None,
    enabled: bool = True,
    threads: int = 1,
    interproscan_bin: str = "interproscan.sh",
    eggnog_bin: str = "emapper.py",
    applications: Iterable[str] | None = None,
    extra_args: Iterable[str | Path] | None = None,
) -> dict[str, FunctionalRecord]:
    """Run the functional-annotation hook → ``{transcript_id: FunctionalRecord}``.

    **Off by default at the call site** (``enabled=False`` short-circuits to ``{}``
    and runs no subprocess) so the default/golden pipeline path is untouched.
    ``proteins_fa`` is the protein FASTA written upstream
    (:func:`export.writers.write_protein_fasta`) — its headers are transcript ids,
    so the returned records key on transcript id. ``tool`` ∈ :data:`FUNCTION_TOOLS`
    (``interproscan`` / ``eggnog`` / ``both``); with ``both`` the two records for a
    transcript are merged (union of terms; ``domain_complete`` OR-ed).

    ``genes`` is accepted for symmetry / future per-biotype scoping but is not
    mutated here — this function only reads the proteins and returns metrics;
    structural attachment happens in the GFF3 writer and
    :func:`apply_domain_credibility` (count-neutral).
    """
    if not enabled:
        return {}
    if tool not in FUNCTION_TOOLS:
        raise ValueError(f"tool must be one of {FUNCTION_TOOLS}, got {tool!r}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records: dict[str, FunctionalRecord] = {}
    if tool in ("interproscan", "both"):
        tsv = run_interproscan(
            proteins_fa,
            out_dir / "interproscan.tsv",
            interproscan_bin=interproscan_bin,
            applications=applications,
            threads=threads,
            extra_args=extra_args,
        )
        _merge_records(records, parse_interproscan_tsv(tsv))
    if tool in ("eggnog", "both"):
        ann = run_eggnog_mapper(
            proteins_fa,
            out_dir / "eggnog",
            eggnog_bin=eggnog_bin,
            db=db,
            threads=threads,
            extra_args=extra_args,
        )
        _merge_records(records, parse_eggnog_annotations(ann))

    _log.info(
        "functional annotation (%s): %d transcripts annotated (%d domain-complete)",
        tool,
        len(records),
        sum(1 for r in records.values() if r.domain_complete),
    )
    return records


def _merge_records(
    into: dict[str, FunctionalRecord], extra: dict[str, FunctionalRecord]
) -> None:
    """Union ``extra`` into ``into`` (terms unioned, ``domain_complete`` OR-ed)."""
    for tid, rec in extra.items():
        cur = into.get(tid)
        if cur is None:
            into[tid] = rec
            continue
        into[tid] = _norm(
            tid,
            set(cur.go_terms) | set(rec.go_terms),
            set(cur.dbxrefs) | set(rec.dbxrefs),
            cur.domain_complete or rec.domain_complete,
        )


# ---------------------------------------------------------------------------
# D2 — domain-completeness as a soft ORF-credibility signal
# ---------------------------------------------------------------------------


def apply_domain_credibility(
    genes: list[ReconciledGene],
    functional: dict[str, FunctionalRecord],
) -> list[ReconciledGene]:
    """Attach the :data:`~qc.flags.DOMAIN_COMPLETE` INFO flag to credible ORFs (D2).

    A gene whose primary transcript spans a **complete recognized domain** is more
    credible than a bare ORF (§2.2). This records that as an INFO flag on the gene
    — a soft, orthogonal credibility signal surfaced in QC reporting.

    **Count-neutral**: it adds an INFO flag only; it never re-tiers, re-types,
    or re-calls a gene, and an empty
    ``functional`` map (the default pipeline path) returns ``genes`` unchanged.
    The *intended future use* is to fold ``domain_complete`` into Tier-1
    eligibility and pseudogene discrimination alongside ``has_homology`` (a CDS
    bracketing a complete Pfam domain promoted over a homology-only ORF); that
    scoring change is deliberately deferred so this phase stays count-neutral.
    """
    if not functional:
        return genes
    out: list[ReconciledGene] = []
    for g in genes:
        rec = functional.get(g.primary_transcript_id)
        if rec is not None and rec.domain_complete:
            out.append(attrs.evolve(g, flags=dedup_flags([*g.flags, DOMAIN_COMPLETE])))
        else:
            out.append(g)
    return out
