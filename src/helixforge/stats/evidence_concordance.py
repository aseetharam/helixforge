"""Per-gene / per-isoform evidence-vs-model concordance stats."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from helixforge.reconcile.validate import (
    check_internal_stops as _check_internal_stops,
)
from helixforge.reconcile.validate import (
    extract_cds_sequence,
)
from helixforge.utils.sequences import is_start_codon, is_stop_codon, translate

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene, SpliceJunction

# ---------------------------------------------------------------------------
# Small interval helpers (shared with as_events conventions)
# ---------------------------------------------------------------------------


def _bounds(item: Any) -> tuple[int, int]:
    """(start, end) from an Exon/CDSSegment/Interval/SpliceJunction or a tuple."""
    if hasattr(item, "start"):
        return item.start, item.end
    if hasattr(item, "donor"):
        return item.donor, item.acceptor
    return item[0], item[1]


def _intron_bounds(tx: Any) -> list[tuple[int, int]]:
    """List of ``(start, end)`` intron tuples for a transcript-like object.

    Prefers the model ``.introns`` property; otherwise derives gaps between
    sorted exon-like intervals (``.exons`` or ``.cds_segments``).
    """
    if hasattr(tx, "introns") and not callable(tx.introns):
        return [(_bounds(i)[0], _bounds(i)[1]) for i in tx.introns]
    exons = _exon_like(tx)
    return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)]


def _exon_like(item: Any) -> list[tuple[int, int]]:
    """Sorted exon-like ``(start, end)`` intervals (``.exons`` else ``.cds_segments``)."""
    feats = None
    if hasattr(item, "exons") and item.exons:
        feats = item.exons
    elif hasattr(item, "cds_segments") and item.cds_segments:
        feats = item.cds_segments
    elif hasattr(item, "cds") and item.cds:
        feats = item.cds
    if not feats:
        return []
    return sorted(_bounds(f) for f in feats)


def _f1(precision: float | None, recall: float | None) -> float | None:
    if precision is None or recall is None:
        return None
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


# ---------------------------------------------------------------------------
# AED: composite Annotation Edit Distance against the evidence
# ---------------------------------------------------------------------------


def compute_aed(
    junction_support: float | None,
    expression: float | None,
    protein_coverage: float | None,
) -> float | None:
    """Composite AED in ``[0, 1]`` (0 = perfect agreement with the evidence).

    Each argument is an independent evidence-agreement signal already scaled to
    ``[0, 1]`` (or ``None`` when that evidence is unavailable for the isoform):

    - ``junction_support``, fraction of introns confirmed by splice junctions,
    - ``expression``, RNA-seq expression presence/level signal,
    - ``protein_coverage``, fraction of a homologous protein covered by the ORF.

    AED is ``1 - mean(available signals)``. With **no** signals available the
    distance is undefined → ``None``. This is the deliberately simple composite
    the design asks for; the per-source detail lives in
    :func:`evidence_agreement_matrix`.
    """
    signals = [
        min(1.0, max(0.0, float(v)))
        for v in (junction_support, expression, protein_coverage)
        if v is not None
    ]
    if not signals:
        return None
    return round(1.0 - sum(signals) / len(signals), 6)


# ---------------------------------------------------------------------------
# D1.1, intron concordance vs the filtered splice-junction set
# ---------------------------------------------------------------------------


def _relevant_junctions(
    gene: Any,
    junctions: list[SpliceJunction],
    min_reads: int,
) -> list[tuple[int, int]]:
    """Junctions on the gene's seqid+strand, inside its span, with enough reads."""
    out: list[tuple[int, int]] = []
    for j in junctions:
        if j.seqid != gene.seqid or j.strand != gene.strand:
            continue
        if j.read_count < min_reads:
            continue
        if j.donor >= gene.start and j.acceptor <= gene.end:
            out.append((j.donor, j.acceptor))
    return out


def intron_concordance(
    gene: Any,
    junctions: list[SpliceJunction],
    min_reads: int = 3,
) -> dict[str, Any]:
    """Per-isoform intron support vs a ``SpliceJunction`` set.

    For each isoform every intron is classified as:

    - **supported**, an exact junction (donor==intron.start, acceptor==intron.end,
      same strand, ``read_count >= min_reads``) exists;
    - **contradicted**, not supported, but a qualifying junction shares exactly
      one boundary (the evidence places a *different* splice site here);
    - **novel**, neither supported nor contradicted (no junction touches it).

    Treating the qualifying junction set as ground truth gives precision/recall/F1
    of the isoform's intron chain:

    - ``precision = supported / num_introns`` (model introns that are confirmed),
    - ``recall = matched_junctions / relevant_junctions`` (evidence reproduced),
    - ``f1``, their harmonic mean.

    Returns ``{"gene_id", "isoforms": {tid: {...}}}``. Mono-exon isoforms have no
    introns → counts are ``0`` and precision/recall/f1 are ``None``.
    """
    qualifying = set(_relevant_junctions(gene, junctions, min_reads))
    donors = {d for d, _ in qualifying}
    acceptors = {a for _, a in qualifying}

    per: dict[str, dict[str, Any]] = {}
    for tx in gene.transcripts:
        introns = _intron_bounds(tx)
        supported = contradicted = novel = 0
        matched: set[tuple[int, int]] = set()
        for s, e in introns:
            if (s, e) in qualifying:
                supported += 1
                matched.add((s, e))
            elif s in donors or e in acceptors:
                contradicted += 1
            else:
                novel += 1

        num = len(introns)
        precision: float | None = supported / num if num else None
        recall: float | None = (
            len(matched & qualifying) / len(qualifying) if qualifying else None
        )
        per[tx.transcript_id] = {
            "num_introns": num,
            "supported": supported,
            "contradicted": contradicted,
            "novel": novel,
            "precision": precision,
            "recall": recall,
            "f1": _f1(precision, recall),
        }
    return {"gene_id": gene.gene_id, "isoforms": per}


# ---------------------------------------------------------------------------
# D1.2, CDS completeness / ORF quality
# ---------------------------------------------------------------------------


def _start_codon_seq(tx: Any, genome: Any) -> str | None:
    """The 3-nt start codon in coding direction (``None`` if unreadable)."""
    seq = extract_cds_sequence(tx, genome)
    if seq is None or len(seq) < 3:
        return None
    return seq[:3]


def _stop_codon_seq(tx: Any, genome: Any) -> str | None:
    """The 3-nt codon immediately 3' of the CDS in coding direction.

    CDS excludes the stop codon, so the stop lives just past the 3'
    end. Window selection mirrors ``validate.check_stop_codon``; codon *identity*
    is decided by ``utils.sequences`` (no duplicated codon logic). ``None`` if the
    window runs off the contig.
    """
    if not tx.cds:
        return None
    first, last = tx.cds[0], tx.cds[-1]
    try:
        if tx.strand == "+":
            if last.end < 0:
                return None
            codon = genome.get_sequence(tx.seqid, last.end, last.end + 3, "+")
        else:
            if first.start - 3 < 0:
                return None
            codon = genome.get_sequence(tx.seqid, first.start - 3, first.start, "-")
    except (KeyError, IndexError, ValueError):
        return None
    return codon if len(codon) == 3 else None


def cds_completeness(
    gene: Any,
    genome: Any,
    protein_lengths: dict[str, int] | None = None,
) -> dict[str, Any]:
    """Per-isoform ORF completeness/quality.

    Per coding isoform returns ``has_start`` (ATG), ``has_stop``,
    ``internal_stop_free``, ``cds_length``, ``protein_coverage`` (translated
    length / ``protein_lengths[protein_id]`` when supplied, else ``None``),
    ``cds_partial``, and a composite ``aed`` (see :func:`compute_aed`). Isoforms
    with no CDS get ``has_cds=False`` and ``None`` codon fields.

    ``genome`` may be ``None`` (sequence checks then report ``None``).
    """
    protein_lengths = protein_lengths or {}
    per: dict[str, dict[str, Any]] = {}
    for tx in gene.transcripts:
        if not tx.cds:
            per[tx.transcript_id] = {
                "has_cds": False,
                "has_start": None,
                "has_stop": None,
                "internal_stop_free": None,
                "cds_length": 0,
                "cds_partial": tx.cds_partial,
                "protein_coverage": None,
                "aed": _isoform_aed(tx, None),
            }
            continue

        has_start: bool | None = None
        has_stop: bool | None = None
        internal_free: bool | None = None
        protein_cov: float | None = None
        if genome is not None:
            start_codon = _start_codon_seq(tx, genome)
            has_start = is_start_codon(start_codon) if start_codon else False
            stop_codon = _stop_codon_seq(tx, genome)
            has_stop = is_stop_codon(stop_codon) if stop_codon else False
            internal_free = _check_internal_stops(tx, genome) is None

            if tx.protein_id and tx.protein_id in protein_lengths:
                cds_seq = extract_cds_sequence(tx, genome)
                if cds_seq is not None:
                    first = tx.cds[0] if tx.strand == "+" else tx.cds[-1]
                    aa = translate(cds_seq, first.phase).rstrip("*")
                    ref_len = protein_lengths[tx.protein_id]
                    if ref_len:
                        protein_cov = min(1.0, len(aa) / ref_len)

        per[tx.transcript_id] = {
            "has_cds": True,
            "has_start": has_start,
            "has_stop": has_stop,
            "internal_stop_free": internal_free,
            "cds_length": tx.total_cds_length,
            "cds_partial": tx.cds_partial,
            "protein_coverage": protein_cov,
            "aed": _isoform_aed(tx, protein_cov),
        }
    return {"gene_id": gene.gene_id, "isoforms": per}


def _expression_signal(tx: Any) -> float | None:
    """Presence signal in [0,1] from TPM (``None`` when TPM is unknown)."""
    if tx.tpm is None:
        return None
    return 1.0 if tx.tpm > 0 else 0.0


def _isoform_aed(tx: Any, protein_coverage: float | None) -> float | None:
    return compute_aed(
        tx.junction_support_fraction, _expression_signal(tx), protein_coverage
    )


# ---------------------------------------------------------------------------
# D1.3, evidence agreement matrix
# ---------------------------------------------------------------------------


def _overlap_bases(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    total = 0
    for sa, ea in a:
        for sb, eb in b:
            ov = min(ea, eb) - max(sa, sb)
            if ov > 0:
                total += ov
    return total


def _coverage_fraction(
    model: list[tuple[int, int]], evidence: list[tuple[int, int]]
) -> float | None:
    """Fraction of the model intervals' bases covered by the evidence intervals."""
    span = sum(e - s for s, e in model)
    if span == 0:
        return None
    return _overlap_bases(model, evidence) / span


def evidence_agreement_matrix(
    gene: Any,
    sources: dict[str, list[Any]],
) -> dict[str, dict[str, dict[str, float | None]]]:
    """Exon/intron agreement of each isoform against each evidence source.

    ``sources`` maps a source name (``"helixer"``, a StringTie sample id,
    ``"miniprot"``) to a list of feature-bearing objects (transcripts / loci /
    alignments). For each isoform × source we report:

    - ``exon_overlap``, fraction of the isoform's exonic bases covered by any of
      the source's exon-like intervals (``.exons`` else ``.cds_segments``),
    - ``intron_match``, fraction of the isoform's introns that exactly match an
      intron implied by the source.

    Shows which evidence drove each model. Returns ``{tid: {source: {...}}}``.
    """
    # Pre-compute the merged exon/intron intervals per source.
    src_exons: dict[str, list[tuple[int, int]]] = {}
    src_introns: dict[str, set[tuple[int, int]]] = {}
    for name, items in sources.items():
        exons: list[tuple[int, int]] = []
        introns: set[tuple[int, int]] = set()
        for it in items:
            ex = _exon_like(it)
            exons.extend(ex)
            introns.update(_intron_bounds(it))
        src_exons[name] = exons
        src_introns[name] = introns

    out: dict[str, dict[str, dict[str, float | None]]] = {}
    for tx in gene.transcripts:
        tx_exons = _exon_like(tx)
        tx_introns = _intron_bounds(tx)
        row: dict[str, dict[str, float | None]] = {}
        for name in sources:
            exon_ov = _coverage_fraction(tx_exons, src_exons[name])
            intron_match: float | None
            if tx_introns:
                matched = sum(1 for i in tx_introns if i in src_introns[name])
                intron_match = matched / len(tx_introns)
            else:
                intron_match = None
            row[name] = {"exon_overlap": exon_ov, "intron_match": intron_match}
        out[tx.transcript_id] = row
    return out


# ---------------------------------------------------------------------------
# D1.4, tidy concordance table
# ---------------------------------------------------------------------------

CONCORDANCE_COLUMNS = (
    "gene_id",
    "transcript_id",
    "is_primary",
    "tier",
    "origin",
    "strand",
    "num_exons",
    "num_introns",
    "supported_introns",
    "contradicted_introns",
    "novel_introns",
    "intron_precision",
    "intron_recall",
    "intron_f1",
    "has_cds",
    "has_start",
    "has_stop",
    "internal_stop_free",
    "cds_length",
    "cds_partial",
    "protein_id",
    "protein_coverage",
    "tpm",
    "junction_support",
    "helixer_support",
    "combined_score",
    "aed",
    "has_homology",
    "flags",
)


def concordance_table(
    genes: list[ReconciledGene],
    junctions: list[SpliceJunction],
    genome: Any,
    sources: dict[str, list[Any]] | None = None,
    protein_lengths: dict[str, int] | None = None,
) -> Any:
    """One tidy row per isoform combining intron + CDS concordance.

    ``junctions`` is the splice-junction set, ``genome`` a ``GenomeAccessor``-like
    object (or ``None``). When ``sources`` is given, per-source agreement columns
    (``agree_<source>_exon`` / ``agree_<source>_intron``) are appended. Returns a
    ``pandas.DataFrame``.
    """
    import pandas as pd

    rows: list[dict[str, Any]] = []
    sources_nn: dict[str, list[Any]] = sources or {}
    source_names = list(sources_nn)
    for gene in genes:
        intron = intron_concordance(gene, junctions)["isoforms"]
        cds = cds_completeness(gene, genome, protein_lengths=protein_lengths)[
            "isoforms"
        ]
        agree = (
            evidence_agreement_matrix(gene, {n: sources_nn[n] for n in source_names})
            if source_names
            else {}
        )
        flag_names = ",".join(f.name for f in gene.flags)
        for tx in gene.transcripts:
            ic = intron[tx.transcript_id]
            cc = cds[tx.transcript_id]
            row: dict[str, Any] = {
                "gene_id": gene.gene_id,
                "transcript_id": tx.transcript_id,
                "is_primary": tx.transcript_id == gene.primary_transcript_id,
                "tier": gene.tier,
                "origin": gene.origin,
                "strand": gene.strand,
                "num_exons": tx.num_exons,
                "num_introns": ic["num_introns"],
                "supported_introns": ic["supported"],
                "contradicted_introns": ic["contradicted"],
                "novel_introns": ic["novel"],
                "intron_precision": ic["precision"],
                "intron_recall": ic["recall"],
                "intron_f1": ic["f1"],
                "has_cds": cc["has_cds"],
                "has_start": cc["has_start"],
                "has_stop": cc["has_stop"],
                "internal_stop_free": cc["internal_stop_free"],
                "cds_length": cc["cds_length"],
                "cds_partial": cc["cds_partial"],
                "protein_id": tx.protein_id,
                "protein_coverage": cc["protein_coverage"],
                "tpm": tx.tpm,
                "junction_support": tx.junction_support_fraction,
                "helixer_support": tx.confidence,
                "combined_score": tx.combined_score,
                "aed": cc["aed"],
                "has_homology": tx.has_homology,
                "flags": flag_names,
            }
            for name in source_names:
                a = agree[tx.transcript_id][name]
                row[f"agree_{name}_exon"] = a["exon_overlap"]
                row[f"agree_{name}_intron"] = a["intron_match"]
            rows.append(row)

    columns = list(CONCORDANCE_COLUMNS)
    for name in source_names:
        columns.extend([f"agree_{name}_exon", f"agree_{name}_intron"])
    return pd.DataFrame(rows, columns=columns)


def write_concordance_tsv(df: Any, path: str) -> str:
    """Write the concordance table to a TSV (tab-separated, no index)."""
    df.to_csv(path, sep="\t", index=False)
    return str(path)
