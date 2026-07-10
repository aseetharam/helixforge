"""StringTie GTF parsing."""

from __future__ import annotations

import os
from typing import Any, Iterable, Sequence

from helixforge.reconcile.models import Exon, StringTieTranscript

# One aggregated StringTie structure for overlap-based TPM lookup:
# ``(strand, start, end, exon_tuples, max_tpm)``.
TpmRecord = tuple[str, int, int, tuple[tuple[int, int], ...], float]


def _parse_attributes(attr_field: str) -> dict[str, str]:
    """Parse a GTF column-9 attribute string into a dict.

    Attributes look like ``key "value"; key2 "value2";`` — values may or may not
    be quoted. Keys are returned verbatim (caller lowercases as needed).
    """
    attrs: dict[str, str] = {}
    for chunk in attr_field.strip().split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        key, _, value = chunk.partition(" ")
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _get_ci(attrs: dict[str, str], name: str) -> str | None:
    """Case-insensitive attribute lookup (e.g. TPM vs tpm)."""
    name_low = name.lower()
    for key, value in attrs.items():
        if key.lower() == name_low:
            return value
    return None


class StringTieParser:
    """Parser for StringTie GTF assemblies."""

    def parse_gtf(
        self, gtf_path: str | os.PathLike[str], sample_id: str
    ) -> list[StringTieTranscript]:
        """Parse one StringTie GTF → ``list[StringTieTranscript]`` (TPM==0 dropped).

        Exons are grouped by ``transcript_id``; the transcript span is the
        min/max of its exons. TPM and coverage come from the transcript record
        (case-insensitive attribute names). Coordinates are converted to
        internal 0-based half-open.
        """
        if not os.path.exists(gtf_path):
            raise FileNotFoundError(f"StringTie GTF not found: {gtf_path}")

        # transcript_id -> accumulated record
        records: dict[str, dict[str, object]] = {}
        with open(gtf_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                (
                    seqid,
                    _src,
                    feature,
                    start_s,
                    end_s,
                    _score,
                    strand,
                    _frame,
                    attr_field,
                ) = cols[:9]
                attrs = _parse_attributes(attr_field)
                tid = attrs.get("transcript_id")
                if tid is None:
                    continue
                rec = records.setdefault(
                    tid,
                    {
                        "gene_id": attrs.get("gene_id", tid),
                        "seqid": seqid,
                        "strand": strand,
                        "exons": [],
                        "tpm": None,
                        "coverage": None,
                    },
                )
                if feature == "transcript":
                    tpm = _get_ci(attrs, "TPM")
                    cov = _get_ci(attrs, "cov")
                    if tpm is not None:
                        rec["tpm"] = float(tpm)
                    if cov is not None:
                        rec["coverage"] = float(cov)
                elif feature == "exon":
                    # 1-based inclusive -> 0-based half-open
                    exons_list = rec["exons"]
                    assert isinstance(exons_list, list)
                    exons_list.append(Exon(int(start_s) - 1, int(end_s)))
                    # capture TPM/cov if only present on exon lines
                    if rec["tpm"] is None:
                        tpm = _get_ci(attrs, "TPM")
                        if tpm is not None:
                            rec["tpm"] = float(tpm)
                    if rec["coverage"] is None:
                        cov = _get_ci(attrs, "cov")
                        if cov is not None:
                            rec["coverage"] = float(cov)

        transcripts: list[StringTieTranscript] = []
        for tid, rec in records.items():
            exons_list = rec["exons"]
            assert isinstance(exons_list, list)
            if not exons_list:
                continue
            # StringTie emits unstranded ('.') transcripts (single-exon, no
            # junction); the data model requires +/- and Mikado cannot use them
            # as stranded evidence — skip them.
            rec_strand = rec["strand"]
            assert isinstance(rec_strand, str)
            if rec_strand not in ("+", "-"):
                continue
            tpm_val = rec["tpm"] if rec["tpm"] is not None else 0.0
            assert isinstance(tpm_val, float)
            if tpm_val == 0:
                continue  # filter TPM == 0
            exons: list[Exon] = sorted(exons_list, key=lambda e: e.start)
            start = exons[0].start
            end = exons[-1].end
            rec_seqid = rec["seqid"]
            assert isinstance(rec_seqid, str)
            rec_cov = rec["coverage"]
            cov_val: float | None = (
                float(rec_cov) if isinstance(rec_cov, (int, float)) else None
            )
            transcripts.append(
                StringTieTranscript(
                    transcript_id=tid,
                    gene_id=str(rec["gene_id"]),
                    seqid=rec_seqid,
                    start=start,
                    end=end,
                    strand=rec_strand,
                    exons=exons,
                    tpm=tpm_val,
                    sample_id=sample_id,
                    coverage=cov_val,
                )
            )
        transcripts.sort(key=lambda t: (t.seqid, t.start))
        return transcripts

    def parse_sample_list(
        self, list_path: str | os.PathLike[str]
    ) -> list[StringTieTranscript]:
        """Parse a file of GTF paths (one per line) → all transcripts.

        Blank lines and ``#`` comments are skipped. ``sample_id`` is the GTF
        filename stem.
        """
        if not os.path.exists(list_path):
            raise FileNotFoundError(f"sample list not found: {list_path}")
        out: list[StringTieTranscript] = []
        with open(list_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                sample_id = os.path.splitext(os.path.basename(line))[0]
                out.extend(self.parse_gtf(line, sample_id))
        return out

    def parse_paths(
        self, paths: Iterable[str | os.PathLike[str]]
    ) -> list[StringTieTranscript]:
        """Parse a list of individual GTF paths → all transcripts.

        The per-sample counterpart to :meth:`parse_sample_list` (which reads the
        same paths from a FOFN): ``sample_id`` is each GTF's filename stem. Every
        path must exist (a clear error beats a silent skip).
        """
        out: list[StringTieTranscript] = []
        for path in paths:
            path_str = os.fspath(path)
            if not os.path.exists(path_str):
                raise FileNotFoundError(f"StringTie GTF not found: {path_str}")
            sample_id = os.path.splitext(os.path.basename(path_str))[0]
            out.extend(self.parse_gtf(path_str, sample_id))
        return out

    @staticmethod
    def _structure_key(t: StringTieTranscript) -> str:
        exon_tuples = tuple((e.start, e.end) for e in t.exons)
        return f"{t.seqid}:{t.strand}:{exon_tuples}"

    def aggregate_across_samples(
        self,
        transcripts: list[StringTieTranscript],
        min_tpm: float = 0.5,
        min_samples: int = 1,
    ) -> dict[str, dict[str, object]]:
        """Group transcripts by identical structure across samples.

        Returns ``{structure_key: {representative, max_tpm, mean_tpm,
        num_samples, sample_tpms}}``. Only structures with ``max_tpm >= min_tpm``
        and ``num_samples >= min_samples`` are kept. This feeds classification
        and the TPM external metric — **not** isoform selection (Mikado does the
        real cross-sample reconciliation later).
        """
        groups: dict[str, dict[str, object]] = {}
        for t in transcripts:
            key = self._structure_key(t)
            g = groups.setdefault(key, {"representative": t, "sample_tpms": {}})
            # keep highest-TPM transcript as representative
            rep = g["representative"]
            assert isinstance(rep, StringTieTranscript)
            if t.tpm > rep.tpm:
                g["representative"] = t
            # per-sample max TPM (a sample may assemble the structure more than once)
            sample_tpms = g["sample_tpms"]
            assert isinstance(sample_tpms, dict)
            prev = sample_tpms.get(t.sample_id)
            if prev is None or t.tpm > prev:
                sample_tpms[t.sample_id] = t.tpm

        result: dict[str, dict[str, object]] = {}
        for key, g in groups.items():
            sample_tpms = g["sample_tpms"]
            assert isinstance(sample_tpms, dict)
            tpms: list[float] = list(sample_tpms.values())
            max_tpm = max(tpms)
            num_samples = len(tpms)
            if max_tpm < min_tpm or num_samples < min_samples:
                continue
            result[key] = {
                "representative": g["representative"],
                "max_tpm": max_tpm,
                "mean_tpm": sum(tpms) / num_samples,
                "num_samples": num_samples,
                "sample_tpms": dict(sample_tpms),
            }
        return result


# ---------------------------------------------------------------------------
# Overlap-based per-model TPM lookup (shared by `evidence` and `reconcile`)
# ---------------------------------------------------------------------------
#
# A Helixer/reconciled model almost never matches a StringTie assembly
# exon-for-exon, so matching by exact structure hash leaves ``tpm`` empty for
# nearly every model. Both the standalone ``evidence`` scorer and the
# ``reconcile`` pipeline instead assign a model its TPM from the same-strand
# StringTie structure with the greatest exonic base overlap. These two functions
# are the single source of that logic so the two commands cannot drift apart.


def tpm_overlap_index_from_transcripts(
    transcripts: list[StringTieTranscript],
    *,
    min_tpm: float = 0.0,
    min_samples: int = 1,
) -> dict[str, list[TpmRecord]]:
    """Aggregate already-parsed StringTie transcripts into an overlap-lookup index.

    Returns ``{seqid: [(strand, start, end, exon_tuples, max_tpm), ...]}``.
    Structures are aggregated across samples (max TPM per identical structure)
    exactly as classification does; ``min_tpm=0.0`` keeps every structure so the
    threshold comparison happens later against the configured ``min_tpm``. Used by
    the reconcile pipeline, which has already parsed the GTFs once.
    """
    parser = StringTieParser()
    aggregated = parser.aggregate_across_samples(
        transcripts, min_tpm=min_tpm, min_samples=min_samples
    )
    index: dict[str, list[TpmRecord]] = {}
    for agg in aggregated.values():
        rep = agg["representative"]
        assert isinstance(rep, StringTieTranscript)
        exon_tuples = tuple((e.start, e.end) for e in rep.exons)
        index.setdefault(rep.seqid, []).append(
            (rep.strand, rep.start, rep.end, exon_tuples, float(agg["max_tpm"]))  # type: ignore[arg-type]
        )
    return index


def build_tpm_overlap_index(
    gtf_paths: Sequence[str | os.PathLike[str]],
    *,
    min_tpm: float = 0.0,
    min_samples: int = 1,
) -> dict[str, list[TpmRecord]]:
    """Parse StringTie GTFs and build the per-model TPM overlap-lookup index.

    Convenience wrapper over :func:`tpm_overlap_index_from_transcripts` for
    callers (the ``evidence`` scorer) that start from paths.
    """
    transcripts = StringTieParser().parse_paths(gtf_paths)
    return tpm_overlap_index_from_transcripts(
        transcripts, min_tpm=min_tpm, min_samples=min_samples
    )


def _exonic_overlap(
    model: list[tuple[int, int]], other: Sequence[tuple[int, int]]
) -> int:
    """Total base overlap between two lists of ``[start, end)`` intervals."""
    total = 0
    for s1, e1 in model:
        for s2, e2 in other:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi > lo:
                total += hi - lo
    return total


def best_overlapping_tpm(
    index: dict[str, list[TpmRecord]],
    seqid: str,
    strand: str,
    exons: Sequence[Any],
) -> float | None:
    """TPM of the same-strand StringTie structure with the greatest exonic overlap.

    The model's exons are intersected base-for-base with each StringTie structure
    on the same ``seqid``/``strand``; the one with the largest total exonic
    overlap wins (ties broken by higher TPM). ``None`` when nothing on the strand
    shares an exonic base.
    """
    model = sorted((e.start, e.end) for e in exons)
    if not model:
        return None
    m_lo, m_hi = model[0][0], model[-1][1]
    best_overlap = 0
    best_tpm: float | None = None
    for st_strand, st_start, st_end, st_exons, tpm in index.get(seqid, ()):
        if st_strand != strand or st_start >= m_hi or st_end <= m_lo:
            continue  # wrong strand or spans disjoint — cheap reject
        overlap = _exonic_overlap(model, st_exons)
        if overlap > best_overlap or (
            overlap == best_overlap
            and overlap > 0
            and (best_tpm is None or tpm > best_tpm)
        ):
            best_overlap, best_tpm = overlap, tpm
    return best_tpm if best_overlap > 0 else None
