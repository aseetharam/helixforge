"""Helixer-input vs HelixForge-output summary + delta table."""

from __future__ import annotations

import math
import os
from collections.abc import Callable
from statistics import mean, median
from typing import TYPE_CHECKING, Any

from helixforge.io.gff import GFF3Parser
from helixforge.stats.evidence_concordance import compute_aed, cds_completeness
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    import pandas as pd

    from helixforge.io.hdf5 import HDF5ConfidenceReader
    from helixforge.reconcile.models import ReconciledGene, SpliceJunction

_log = get_logger(__name__)

NAN = float("nan")

# A normalised gene dict produced by ``_normalize``.
_GeneDict = dict[str, Any]


# ---------------------------------------------------------------------------
# Input normalisation — ReconciledGene list OR Helixer GFF3 path → common shape
# ---------------------------------------------------------------------------


def _is_gff3_input(genes_or_gff3: Any) -> bool:
    """True when the argument is a path-like GFF3 input rather than a gene list."""
    return isinstance(genes_or_gff3, (str, bytes)) or hasattr(
        genes_or_gff3, "__fspath__"
    )


def _normalize(genes_or_gff3: Any) -> list[_GeneDict]:
    """Return a list of normalised gene dicts.

    Each: ``{gene_id, strand, transcripts}`` where each transcript is
    ``{transcript_id, exons, cds, cds_partial, tpm, junction_support, flags,
       protein_id, is_primary}``. ``ReconciledGene`` carries the rich fields;
    GFF3-parsed genes fill structural fields only (evidence fields default to
    ``None``/empty so before/after stays comparable).
    """
    if _is_gff3_input(genes_or_gff3):
        parsed = GFF3Parser(str(genes_or_gff3)).parse_genes_generic()
        out: list[_GeneDict] = []
        for g in parsed:
            txs = [
                {
                    "transcript_id": t["transcript_id"],
                    "exons": t["exons"],
                    "cds": t["cds"],
                    "cds_partial": False,
                    "tpm": None,
                    "junction_support": None,
                    "flags": set(),
                    "protein_id": None,
                    "is_primary": i == 0,
                }
                for i, t in enumerate(g["transcripts"])
            ]
            out.append(
                {
                    "gene_id": g["gene_id"],
                    "seqid": g["seqid"],
                    "strand": g["strand"],
                    "transcripts": txs,
                }
            )
        return out

    out = []
    for g in genes_or_gff3:
        flag_names = {f.name for f in g.flags}
        txs = [
            {
                "transcript_id": t.transcript_id,
                "exons": t.exons,
                "cds": t.cds,
                "cds_partial": t.cds_partial,
                "tpm": t.tpm,
                "junction_support": t.junction_support_fraction,
                "flags": flag_names,
                "protein_id": t.protein_id,
                "is_primary": t.transcript_id == g.primary_transcript_id,
            }
            for t in g.transcripts
        ]
        out.append(
            {
                "gene_id": g.gene_id,
                "seqid": g.seqid,
                "strand": g.strand,
                "transcripts": txs,
            }
        )
    return out


def _representative(gene: _GeneDict) -> dict[str, Any] | None:
    """The primary transcript dict (falls back to the first; ``None`` if empty).

    A well-formed gene always has >= 1 transcript, but a malformed GFF3 (e.g. a
    duplicate gene ID whose child mRNA was de-collided away from its renamed
    parent) can parse into a transcript-less gene. Returning ``None`` lets
    callers skip it instead of crashing with ``IndexError``.
    """
    txs: list[dict[str, Any]] = gene["transcripts"]
    if not txs:
        return None
    for t in txs:
        if t["is_primary"]:
            return t
    return txs[0]


def _cds_length(cds: Any) -> int:
    return sum(seg.end - seg.start for seg in cds) if cds else 0


def _introns_from_exons(exons: Any) -> list[tuple[int, int]]:
    """Intron ``(start, end)`` gaps between sorted exon-like intervals."""
    ex = sorted((e.start, e.end) for e in exons)
    return [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]


def _junction_support_from_set(
    seqid: str,
    strand: str,
    exons: Any,
    junction_set: set[tuple[str, int, int, str]],
) -> float | None:
    """Fraction of a transcript's introns confirmed by ``junction_set``.

    ``junction_set`` is a set of ``(seqid, donor, acceptor, strand)`` tuples
    (already read-count-filtered). ``None`` when the transcript is single-exon
    (no introns) so AED leaves the junction signal unused there.
    """
    introns = _introns_from_exons(exons)
    if not introns:
        return None
    supported = sum(1 for s, e in introns if (seqid, s, e, strand) in junction_set)
    return supported / len(introns)


def _helixer_support_for(
    exons: Any,
    seqid: str,
    strand: str,
    h5_reader: Any,
    intron_high_cutoff: float,
) -> float | None:
    """Per-transcript Helixer support from HDF5; ``None`` if the seqid is absent."""
    from helixforge.mikado.emit_external import helixer_support

    try:
        return helixer_support(
            exons,
            seqid,
            strand,
            h5_reader,
            intron_high_cutoff=intron_high_cutoff,
        )
    except (KeyError, IndexError, ValueError):
        return None


# ---------------------------------------------------------------------------
# D2.1 — annotation summary
# ---------------------------------------------------------------------------


def annotation_summary(
    genes_or_gff3: Any,
    genome: Any = None,
    h5_reader: Any = None,
    junction_set: set[tuple[str, int, int, str]] | None = None,
    protein_lengths: dict[str, int] | None = None,
    intron_high_cutoff: float = 0.5,
) -> dict[str, Any]:
    """Structural + evidence summary of an annotation set (gene list or Helixer GFF3).

    Returns a dict of scalar metrics plus an ``isoforms_per_gene`` distribution.

    - ``%complete ORFs`` uses codon checks when ``genome`` is given *and* the
      input is a ``ReconciledGene`` list; otherwise a structural proxy (has CDS,
      not partial, mod-3, no broken-ORF flag).
    - ``mean_aed`` blends junction support + expression + protein coverage
      (:func:`compute_aed`). For ``ReconciledGene`` input the model's stored
      fields drive it (protein coverage via ``cds_completeness`` when ``genome``
      + ``protein_lengths`` are present). For a Helixer GFF3 the only available
      signal is junction support, computed against ``junction_set`` when given —
      so the before/after AED delta is real; with no signal it stays ``NaN``.
    - ``mean_helixer_support`` is the per-gene (primary-transcript) Helixer
      support from ``h5_reader`` (the SAME HDF5 drives before and after), or
      ``NaN`` when no HDF5 is supplied.
    """
    is_recon = not _is_gff3_input(genes_or_gff3)
    genes = _normalize(genes_or_gff3)
    # Defensively drop transcript-less genes (malformed GFF3, e.g. a duplicate
    # gene ID de-collided from its child mRNA) so the summary degrades with a
    # warning instead of aborting the whole `stats` run.
    empty = [g for g in genes if not g["transcripts"]]
    if empty:
        shown = ", ".join(g["gene_id"] for g in empty[:5])
        _log.warning(
            "annotation_summary: skipping %d transcript-less gene(s) "
            "(likely duplicate/malformed IDs in the input GFF3): %s%s",
            len(empty),
            shown,
            " ..." if len(empty) > 5 else "",
        )
        genes = [g for g in genes if g["transcripts"]]
    n_genes = len(genes)

    iso_counts = [len(g["transcripts"]) for g in genes]
    n_tx = sum(iso_counts)
    distribution: dict[int, int] = {}
    for c in iso_counts:
        distribution[c] = distribution.get(c, 0) + 1

    mono = multi = 0
    cds_lengths: list[int] = []
    complete = 0
    coding = 0
    aeds: list[float] = []
    supports: list[float] = []

    # Codon-accurate completeness + protein-aware AED only for ReconciledGene +
    # genome. One cds_completeness pass per gene feeds both (keyed by tid).
    cc_by_tx: dict[str, dict[str, Any]] | None = None
    if is_recon and genome is not None:
        cc_by_tx = {}
        for g in genes_or_gff3:
            cc_by_tx.update(cds_completeness(g, genome, protein_lengths)["isoforms"])

    for g in genes:
        rep = _representative(g)
        if (
            rep is None
        ):  # transcript-less genes already filtered above; belt-and-suspenders
            continue
        is_multi = any(len(t["exons"]) > 1 for t in g["transcripts"])
        if is_multi:
            multi += 1
        else:
            mono += 1

        cds = rep["cds"]
        if cds:
            coding += 1
            length = _cds_length(cds)
            cds_lengths.append(length)
            if _is_complete(rep, cc_by_tx, length):
                complete += 1

        a = _representative_aed(rep, g, cc_by_tx, junction_set)
        if a is not None:
            aeds.append(a)

        if h5_reader is not None:
            hs = _helixer_support_for(
                rep["exons"], g["seqid"], g["strand"], h5_reader, intron_high_cutoff
            )
            if hs is not None:
                supports.append(hs)

    return {
        "gene_count": n_genes,
        "transcript_count": n_tx,
        "coding_gene_count": coding,
        "isoforms_per_gene_mean": mean(iso_counts) if iso_counts else NAN,
        "isoforms_per_gene_median": median(iso_counts) if iso_counts else NAN,
        "isoforms_per_gene_max": max(iso_counts) if iso_counts else 0,
        "isoforms_per_gene_distribution": distribution,
        "mono_exon_genes": mono,
        "multi_exon_genes": multi,
        "mean_cds_length": mean(cds_lengths) if cds_lengths else NAN,
        "median_cds_length": median(cds_lengths) if cds_lengths else NAN,
        "pct_complete_orfs": (100.0 * complete / coding) if coding else NAN,
        "mean_aed": mean(aeds) if aeds else NAN,
        "mean_helixer_support": mean(supports) if supports else NAN,
    }


def _representative_aed(
    rep: dict[str, Any],
    gene: _GeneDict,
    cc_by_tx: dict[str, dict[str, Any]] | None,
    junction_set: set[tuple[str, int, int, str]] | None,
) -> float | None:
    """AED for a gene's representative transcript (``None`` if no evidence signal).

    ReconciledGene + genome: reuse the protein-aware AED from ``cds_completeness``.
    Otherwise blend the model's junction support (or, for a GFF3, the junction set)
    with the expression presence signal.
    """
    if cc_by_tx is not None:
        cc = cc_by_tx.get(rep["transcript_id"])
        if cc is not None and cc.get("aed") is not None:
            return float(cc["aed"])
    junction = rep["junction_support"]
    if junction is None and junction_set is not None:
        junction = _junction_support_from_set(
            gene["seqid"], gene["strand"], rep["exons"], junction_set
        )
    expression: float | None = (
        None if rep["tpm"] is None else (1.0 if rep["tpm"] > 0 else 0.0)
    )
    return compute_aed(junction, expression, None)


_BROKEN_ORF_FLAGS = frozenset({"NO_START", "NO_STOP", "INTERNAL_STOP"})


def _is_complete(
    rep: dict[str, Any],
    codon: dict[str, dict[str, Any]] | None,
    length: int,
) -> bool:
    """Decide ORF completeness for the representative transcript dict."""
    if codon is not None:
        cc = codon.get(rep["transcript_id"])
        if cc is not None:
            return bool(cc["has_start"] and cc["has_stop"] and cc["internal_stop_free"])
    # Structural proxy: not partial, mod-3, no broken-ORF flag.
    if rep["cds_partial"]:
        return False
    if length % 3 != 0:
        return False
    return not (_BROKEN_ORF_FLAGS & rep["flags"])


# ---------------------------------------------------------------------------
# Completeness metrics — Phase 11 routes these through the real bench tools
# ---------------------------------------------------------------------------
#
# These run a real external tool ONLY when enough inputs are present: a
# ``ReconciledGene`` list (we can translate proteins), an open ``genome``, and a
# ``lineage``/``database``. Otherwise — including every GFF3-path input and the
# default no-lineage call from ``before_after_table`` — they return ``NaN`` (the
# Phase-9 behaviour, preserved). Any tool failure also degrades to ``NaN`` so the
# stats table is never sunk by a missing benchmarking binary.


def _proteome_metric(
    genes_or_gff3: Any,
    genome: Any,
    runner: Callable[[str, str], dict[str, Any]],
    pick: Callable[[dict[str, Any]], float | None],
) -> float:
    """Write proteins to a temp dir, run ``runner``, return ``pick(result)`` or NaN."""
    if genome is None or _is_gff3_input(genes_or_gff3):
        return NAN
    import tempfile

    from helixforge.bench.wrappers import BenchmarkError
    from helixforge.export.writers import write_protein_fasta

    try:
        with tempfile.TemporaryDirectory() as td:
            proteins = os.path.join(td, "proteins.fa")
            write_protein_fasta(genes_or_gff3, genome, proteins)
            result = runner(proteins, os.path.join(td, "out"))
        value = pick(result)
        return NAN if value is None else value
    except (BenchmarkError, OSError, ValueError) as exc:
        # Narrow the catch to the *expected* degradations —
        # a missing/failed external tool (BenchmarkError), an I/O error (OSError),
        # or an unparseable stats file (ValueError). A genuine code bug (TypeError,
        # AttributeError, …) must NOT masquerade as a NaN metric: it propagates.
        _log.warning("completeness metric degraded to NaN: %s", exc)
        return NAN


def busco_completeness(
    genes_or_gff3: Any, genome: Any = None, lineage: str | None = None
) -> float:
    """BUSCO complete % (C) on the protein set; NaN without a lineage/genome."""
    if lineage is None:
        return NAN
    from helixforge.bench.wrappers import run_busco

    return _proteome_metric(
        genes_or_gff3,
        genome,
        lambda fa, out: run_busco(fa, lineage, out),
        lambda res: res.get("complete"),
    )


def compleasm_score(
    genes_or_gff3: Any, genome: Any = None, lineage: str | None = None
) -> float:
    """compleasm complete % (S+D) on the protein set; NaN without a lineage/genome."""
    if lineage is None:
        return NAN
    from helixforge.bench.wrappers import run_compleasm

    return _proteome_metric(
        genes_or_gff3,
        genome,
        lambda fa, out: run_compleasm(fa, lineage, out),
        lambda res: res.get("complete"),
    )


def omark_score(
    genes_or_gff3: Any, genome: Any = None, database: str | None = None
) -> float:
    """OMArk proteome consistency %; NaN without an OMA database/genome."""
    if database is None:
        return NAN
    from helixforge.bench.wrappers import run_omark

    return _proteome_metric(
        genes_or_gff3,
        genome,
        lambda fa, out: run_omark(fa, database, out),
        lambda res: res.get("consistent"),
    )


# ---------------------------------------------------------------------------
# D2.2 — before/after delta table + Markdown summary
# ---------------------------------------------------------------------------

_DELTA_METRICS: tuple[tuple[str, str], ...] = (
    ("gene_count", "Genes"),
    ("transcript_count", "Transcripts"),
    ("coding_gene_count", "Coding genes"),
    ("isoforms_per_gene_mean", "Isoforms/gene (mean)"),
    ("isoforms_per_gene_max", "Isoforms/gene (max)"),
    ("mono_exon_genes", "Mono-exon genes"),
    ("multi_exon_genes", "Multi-exon genes"),
    ("mean_cds_length", "Mean CDS length"),
    ("median_cds_length", "Median CDS length"),
    ("pct_complete_orfs", "% complete ORFs"),
    ("mean_aed", "Mean AED"),
    ("mean_helixer_support", "Mean Helixer support"),
)


def _build_junction_set(
    junctions: list[SpliceJunction] | None,
    min_reads: int,
) -> set[tuple[str, int, int, str]] | None:
    """``{(seqid, donor, acceptor, strand)}`` from a ``SpliceJunction`` list."""
    if not junctions:
        return None
    return {
        (j.seqid, j.donor, j.acceptor, j.strand)
        for j in junctions
        if j.read_count >= min_reads
    }


def before_after_table(
    helixer_gff3: Any,
    reconciled_genes: list[ReconciledGene] | str,
    genome: Any = None,
    h5_path: str | None = None,
    junctions: list[SpliceJunction] | None = None,
    protein_lengths: dict[str, int] | None = None,
    min_reads: int = 3,
) -> pd.DataFrame:
    """Delta table: one row per metric (Helixer in | HelixForge out | Δ).

    Optional drivers populate the previously-stubbed columns:

    - ``h5_path`` — Helixer confidence HDF5 → the ``Mean Helixer support`` row is
      computed for **both** the Helixer input and the reconciled output from the
      same HDF5, so its Δ is meaningful (the ``no_helixer_support`` ablation aside,
      this is the headline coupling metric). Cleanly ``NaN`` when absent.
    - ``junctions`` — a ``SpliceJunction`` list → a real junction signal for the
      Helixer set's AED (the reconciled set carries its own).
    - ``protein_lengths`` + ``genome`` — protein coverage feeds the reconciled
      AED via ``cds_completeness``.

    BUSCO/compleasm/OMArk rows stay ``NaN`` unless a lineage/database is wired
    through the bench tools. Returns a ``pandas.DataFrame`` with columns
    ``metric, helixer, helixforge, delta``.
    """
    import pandas as pd

    from helixforge.io.hdf5 import HDF5ConfidenceReader

    junction_set = _build_junction_set(junctions, min_reads)
    h5_reader: HDF5ConfidenceReader | None = (
        HDF5ConfidenceReader(h5_path) if h5_path else None
    )
    try:
        kw: dict[str, Any] = dict(
            genome=genome,
            h5_reader=h5_reader,
            junction_set=junction_set,
            protein_lengths=protein_lengths,
        )
        before = annotation_summary(helixer_gff3, **kw)
        after = annotation_summary(reconciled_genes, **kw)
    finally:
        if h5_reader is not None:
            h5_reader.close()

    rows: list[dict[str, Any]] = []
    for key, label in _DELTA_METRICS:
        b, a = before[key], after[key]
        delta = a - b if _numeric(a) and _numeric(b) else NAN
        rows.append({"metric": label, "helixer": b, "helixforge": a, "delta": delta})

    bench_fns: list[tuple[str, Callable[..., float]]] = [
        ("BUSCO completeness %", busco_completeness),
        ("compleasm completeness %", compleasm_score),
        ("OMArk consistency %", omark_score),
    ]
    for label, fn in bench_fns:
        b, a = fn(helixer_gff3, genome), fn(reconciled_genes, genome)
        rows.append({"metric": label, "helixer": b, "helixforge": a, "delta": NAN})

    return pd.DataFrame(rows, columns=["metric", "helixer", "helixforge", "delta"])


def _numeric(v: Any) -> bool:
    return isinstance(v, (int, float)) and not (isinstance(v, float) and math.isnan(v))


def _fmt(v: Any) -> str:
    if isinstance(v, float):
        if math.isnan(v):
            return "—"
        return f"{v:.2f}" if v != int(v) else str(int(v))
    return str(v)


def write_summary(df: pd.DataFrame, md_path: str) -> str:
    """Render the delta table to a one-page Markdown table (manuscript figure)."""
    lines = [
        "# HelixForge v3 — before/after summary",
        "",
        "| Metric | Helixer | HelixForge | Δ |",
        "| --- | ---: | ---: | ---: |",
    ]
    for _, r in df.iterrows():
        lines.append(
            f"| {r['metric']} | {_fmt(r['helixer'])} | "
            f"{_fmt(r['helixforge'])} | {_fmt(r['delta'])} |"
        )
    lines.append("")
    with open(md_path, "w") as fh:
        fh.write("\n".join(lines))
    return str(md_path)
