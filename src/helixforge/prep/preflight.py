"""Input-integrity preflight."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from helixforge.io.bam import CSI_REQUIRED_THRESHOLD, _has_csi_index
from helixforge.io.fasta import GenomeAccessor
from helixforge.io.hdf5 import HDF5ConfidenceReader
from helixforge.io.validate import (
    FormatError,
    sniff_bam,
    sniff_fasta,
    sniff_gff3,
    sniff_gtf,
    sniff_star_sj,
)

# ---------------------------------------------------------------------------
# seqid collection
# ---------------------------------------------------------------------------


@dataclass
class SourceSeqids:
    """The seqids (and any known contig lengths) for one input source."""

    label: str
    seqids: set[str]
    lengths: dict[str, int] = field(default_factory=dict)


def _apply_alias(name: str, alias_map: dict[str, str] | None) -> str:
    return alias_map.get(name, name) if alias_map else name


def _gxf_seqids(path: str) -> set[str]:
    """Column-0 seqids of a GFF3/GTF (cheap scan; no gffutils DB)."""
    seqids: set[str] = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            seqid = line.split("\t", 1)[0]
            if seqid:
                seqids.add(seqid)
    return seqids


def _star_sj_seqids(path: str) -> set[str]:
    seqids: set[str] = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            seqid = line.split("\t", 1)[0]
            if seqid:
                seqids.add(seqid)
    return seqids


def _bam_seqids_lengths(path: str) -> tuple[set[str], dict[str, int]]:
    import pysam

    af = pysam.AlignmentFile(str(path), "rb", check_sq=False)
    try:
        refs = list(af.references or ())
        lengths = list(af.lengths or ())
    finally:
        af.close()
    return set(refs), {r: int(ln) for r, ln in zip(refs, lengths)}


def _collect_sources(
    config: object, *, alias_map: dict[str, str] | None
) -> list[SourceSeqids]:
    """Gather a :class:`SourceSeqids` per present input on ``config``.

    Each source's seqids are passed through ``alias_map`` (``{alias: canonical}``)
    so a deliberate translation is honored before any comparison. A source that
    cannot be opened is simply omitted here (the format sniffers report the
    malformed file separately).
    """
    sources: list[SourceSeqids] = []

    def alias_set(names: set[str]) -> set[str]:
        return {_apply_alias(n, alias_map) for n in names}

    def alias_lengths(lengths: dict[str, int]) -> dict[str, int]:
        return {_apply_alias(k, alias_map): v for k, v in lengths.items()}

    genome = getattr(config, "genome_fasta", None)
    if genome and Path(genome).exists():
        try:
            with GenomeAccessor(genome) as ga:
                lengths = ga.get_scaffold_lengths()
            sources.append(
                SourceSeqids(
                    "genome_fasta", alias_set(set(lengths)), alias_lengths(lengths)
                )
            )
        except (OSError, ValueError):
            pass

    helixer_gff3 = getattr(config, "helixer_gff3", None)
    if helixer_gff3 and Path(helixer_gff3).exists():
        sources.append(
            SourceSeqids("helixer_gff3", alias_set(_gxf_seqids(str(helixer_gff3))))
        )

    helixer_h5 = getattr(config, "helixer_h5", None)
    if helixer_h5 and Path(helixer_h5).exists():
        try:
            with HDF5ConfidenceReader(helixer_h5) as reader:
                sources.append(
                    SourceSeqids("helixer_h5", alias_set(set(reader.seqids)))
                )
        except (OSError, ValueError, FileNotFoundError):
            pass

    for gtf in getattr(config, "stringtie_list", None) or []:
        if Path(gtf).exists():
            sources.append(
                SourceSeqids(
                    f"stringtie:{Path(gtf).name}", alias_set(_gxf_seqids(str(gtf)))
                )
            )

    for bam in getattr(config, "bam_paths", None) or []:
        if Path(bam).exists():
            try:
                names, lengths = _bam_seqids_lengths(str(bam))
            except (OSError, ValueError):
                continue
            sources.append(
                SourceSeqids(
                    f"bam:{Path(bam).name}", alias_set(names), alias_lengths(lengths)
                )
            )

    for sj in getattr(config, "star_sj_paths", None) or []:
        if Path(sj).exists():
            sources.append(
                SourceSeqids(
                    f"star_sj:{Path(sj).name}", alias_set(_star_sj_seqids(str(sj)))
                )
            )

    return sources


# ---------------------------------------------------------------------------
# D1 — concordance
# ---------------------------------------------------------------------------


@dataclass
class ConcordanceReport:
    """Result of :func:`check_reference_concordance`."""

    sources: list[SourceSeqids] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def ok(self) -> bool:
        return not self.errors

    def render(self) -> str:
        lines = ["seqid concordance:"]
        for s in self.sources:
            lines.append(f"  {s.label}: {len(s.seqids)} seqids")
        for w in self.warnings:
            lines.append(f"  WARNING: {w}")
        for e in self.errors:
            lines.append(f"  ERROR: {e}")
        lines.append(f"concordance OK: {self.ok()}")
        return "\n".join(lines)


def _sample(names: set[str], n: int = 5) -> str:
    ordered = sorted(names)
    head = ", ".join(ordered[:n])
    return head + ("…" if len(ordered) > n else "")


def check_reference_concordance(
    config: object, *, alias_map: dict[str, str] | None = None
) -> ConcordanceReport:
    """Intersect seqids across every configured input.

    The genome FASTA seqids are the reference substrate. For each evidence
    source: an **entirely disjoint** seqid set (no shared name with the genome)
    is an **error** (a naming mismatch — the classic ``Chr1`` vs ``1`` silent
    failure); a **partial** overlap (some source seqids absent from the genome)
    is a **warning** listing the missing names. Where two sources both carry a
    contig length (FASTA ``.fai``, BAM ``@SQ:LN``), a disagreement on a shared
    contig is an **error** (same name, different assembly). An optional
    ``alias_map`` (``{alias: canonical}``) is applied to every seqid first.
    """
    report = ConcordanceReport()
    report.sources = _collect_sources(config, alias_map=alias_map)
    if not report.sources:
        return report

    ref = next((s for s in report.sources if s.label == "genome_fasta"), None)
    ref_seqids = ref.seqids if ref else set()

    if ref is not None:
        for src in report.sources:
            if src is ref or not src.seqids:
                continue
            shared = src.seqids & ref_seqids
            if ref_seqids and not shared:
                report.errors.append(
                    f"{src.label} shares NO seqid with the genome FASTA — likely "
                    f"a naming mismatch. {src.label}: [{_sample(src.seqids)}] vs "
                    f"genome: [{_sample(ref_seqids)}]"
                )
                continue
            missing = src.seqids - ref_seqids
            if missing:
                report.warnings.append(
                    f"{src.label} has {len(missing)} seqid(s) absent from the "
                    f"genome FASTA: [{_sample(missing)}]"
                )

    # Length concordance across any sources that carry lengths.
    length_owners = [s for s in report.sources if s.lengths]
    seen: dict[str, tuple[str, int]] = {}
    for src in length_owners:
        for seqid, ln in src.lengths.items():
            if seqid in seen:
                prev_label, prev_len = seen[seqid]
                if prev_len != ln:
                    report.errors.append(
                        f"contig {seqid!r} length disagrees: {prev_label}={prev_len} "
                        f"vs {src.label}={ln} (same name, different assembly?)"
                    )
            else:
                seen[seqid] = (src.label, ln)

    return report


# ---------------------------------------------------------------------------
# D3 — CSI gate
# ---------------------------------------------------------------------------


def check_csi_requirement(config: object) -> list[str]:
    """Return error strings for BAMs missing a ``.csi`` index on a large genome.

    If any reference contig exceeds :data:`CSI_REQUIRED_THRESHOLD` (2^29 bp), a
    ``.csi`` index is mandatory for every BAM (``.bai`` cannot address those
    coordinates and a fetch silently returns nothing — biological-assessment
    §3.4). Returns ``[]`` when no large contig is present or every BAM has CSI.
    """
    errors: list[str] = []
    genome = getattr(config, "genome_fasta", None)
    if not genome or not Path(genome).exists():
        return errors
    try:
        with GenomeAccessor(genome) as ga:
            lengths = ga.get_scaffold_lengths()
    except (OSError, ValueError):
        return errors
    big = [name for name, ln in lengths.items() if ln > CSI_REQUIRED_THRESHOLD]
    if not big:
        return errors
    for bam in getattr(config, "bam_paths", None) or []:
        if Path(bam).exists() and not _has_csi_index(str(bam)):
            errors.append(
                f"{Path(bam).name}: genome has contig(s) > 2^29 bp "
                f"([{_sample(set(big))}]) but this BAM has no .csi index "
                "(.bai cannot address those coordinates). Re-index with "
                "`samtools index -c`."
            )
    return errors


# ---------------------------------------------------------------------------
# D3b — work-dir writability
# ---------------------------------------------------------------------------


def check_work_dir(config: object) -> list[str]:
    """Verify the work dir and every stage output dir are creatable and writable.

    Creates the directories up front so a permissions / disk-full failure is
    caught in seconds, not after hours of Mikado + TransDecoder.  Returns a list
    of error strings (empty = OK).
    """
    import tempfile

    errors: list[str] = []
    work_dir = getattr(config, "work_dir", None)
    if not work_dir:
        return errors

    work = Path(work_dir).resolve()
    stage_dirs = [
        work,
        work / "mikado_inputs",
        work / "mikado_run",
    ]
    for d in stage_dirs:
        try:
            d.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            errors.append(f"cannot create stage dir {d}: {exc}")
            continue
        try:
            probe = d / ".helixforge_write_test"
            probe.write_text("ok")
            probe.unlink()
        except OSError as exc:
            errors.append(f"stage dir {d} not writable: {exc}")
    return errors


# ---------------------------------------------------------------------------
# D3c — tool-chain readiness for the configured pipeline
# ---------------------------------------------------------------------------


def check_tool_chain(config: object) -> list[str]:
    """Verify that every external tool the configured pipeline will invoke resolves.

    Returns a list of error strings for missing tools.  Checks only the tools
    the run will actually need (e.g. TransDecoder / DIAMOND are skipped in basic
    mode where no protein DB is configured).
    """
    import shutil

    errors: list[str] = []
    stringtie_list = getattr(config, "stringtie_list", None) or []
    protein_db = getattr(config, "protein_db", None)
    needs_mikado = bool(stringtie_list) and bool(protein_db)

    def _check(label: str, bin_name: str) -> None:
        if not shutil.which(bin_name):
            errors.append(
                f"{label} not found on PATH: {bin_name!r}. Install it or "
                "pass an absolute path via the corresponding *_bin argument."
            )

    if needs_mikado:
        _check("mikado", getattr(config, "mikado_bin", "mikado"))
        _check("diamond", getattr(config, "diamond_bin", "diamond"))
        td_dir = getattr(config, "transdecoder_bin_dir", None)
        td_bin = (
            str(Path(td_dir) / "TransDecoder.LongOrfs")
            if td_dir
            else "TransDecoder.LongOrfs"
        )
        _check("TransDecoder", td_bin)

    portcullis_bin = getattr(config, "portcullis_bin", None)
    bam_paths = getattr(config, "bam_paths", None) or []
    if portcullis_bin and bam_paths:
        _check("portcullis", portcullis_bin)

    return errors


# ---------------------------------------------------------------------------
# D4 / orchestration — the single reusable gate
# ---------------------------------------------------------------------------


@dataclass
class PreflightReport:
    """Combined input-integrity result: format + concordance + CSI + work-dir + tools."""

    format_errors: list[str] = field(default_factory=list)
    concordance: ConcordanceReport | None = None
    csi_errors: list[str] = field(default_factory=list)
    workdir_errors: list[str] = field(default_factory=list)
    tool_errors: list[str] = field(default_factory=list)

    def errors(self) -> list[str]:
        errs = list(self.format_errors)
        if self.concordance is not None:
            errs.extend(self.concordance.errors)
        errs.extend(self.csi_errors)
        errs.extend(self.workdir_errors)
        errs.extend(self.tool_errors)
        return errs

    def warnings(self) -> list[str]:
        return list(self.concordance.warnings) if self.concordance is not None else []

    def ok(self) -> bool:
        return not self.errors()

    def render(self) -> str:
        lines = ["input-integrity preflight:"]
        if self.format_errors:
            lines.append("  format:")
            lines.extend(f"    ERROR: {e}" for e in self.format_errors)
        else:
            lines.append("  format: all sniffed files OK")
        if self.concordance is not None:
            lines.append("  " + self.concordance.render().replace("\n", "\n  "))
        if self.csi_errors:
            lines.append("  CSI:")
            lines.extend(f"    ERROR: {e}" for e in self.csi_errors)
        if self.workdir_errors:
            lines.append("  work-dir:")
            lines.extend(f"    ERROR: {e}" for e in self.workdir_errors)
        if self.tool_errors:
            lines.append("  tool-chain:")
            lines.extend(f"    ERROR: {e}" for e in self.tool_errors)
        lines.append(f"preflight OK: {self.ok()}")
        return "\n".join(lines)


def _sniff_all(config: object) -> list[str]:
    """Run the per-format sniffers over every present input; collect errors."""
    errors: list[str] = []

    def run(sniffer: Callable[[str], None], path: str | None) -> None:
        if path and Path(path).exists():
            try:
                sniffer(path)
            except FormatError as exc:
                errors.append(str(exc))

    run(sniff_fasta, getattr(config, "genome_fasta", None))
    run(sniff_gff3, getattr(config, "helixer_gff3", None))
    miniprot = getattr(config, "miniprot_gff", None)
    if miniprot:
        run(sniff_gff3, miniprot)
    for gtf in getattr(config, "stringtie_list", None) or []:
        run(sniff_gtf, gtf)
    for bam in getattr(config, "bam_paths", None) or []:
        run(sniff_bam, bam)
    for sj in getattr(config, "star_sj_paths", None) or []:
        run(sniff_star_sj, sj)
    return errors


def run_preflight(
    config: object,
    *,
    alias_map: dict[str, str] | None = None,
    check_workdir: bool = True,
    check_tools: bool = True,
) -> PreflightReport:
    """The single input-integrity gate (format + concordance + CSI + work-dir + tools).

    Reused verbatim by ``helixforge doctor`` and ``helixforge run`` — never
    duplicated. Non-mutating: it only reads and reports, so a passing preflight
    leaves a run's gene/tier/origin/AS counts identical (count-neutral).

    ``check_workdir`` and ``check_tools`` gate the work-dir writability and
    tool-chain resolution checks (on by default; ``doctor`` without inputs
    skips them via the ``SimpleNamespace`` that has no ``work_dir``).
    """
    return PreflightReport(
        format_errors=_sniff_all(config),
        concordance=check_reference_concordance(config, alias_map=alias_map),
        csi_errors=check_csi_requirement(config),
        workdir_errors=check_work_dir(config) if check_workdir else [],
        tool_errors=check_tool_chain(config) if check_tools else [],
    )


class PreflightError(RuntimeError):
    """Raised by :func:`require_preflight` when the integrity gate fails."""


def require_preflight(
    config: object, *, alias_map: dict[str, str] | None = None
) -> PreflightReport:
    """Run :func:`run_preflight` and raise :class:`PreflightError` if not OK.

    The required form for ``helixforge run`` — a malformed/inconsistent input
    stops the run loudly instead of silently degrading the annotation.
    """
    report = run_preflight(config, alias_map=alias_map)
    if not report.ok():
        raise PreflightError("input-integrity preflight failed:\n" + report.render())
    return report
