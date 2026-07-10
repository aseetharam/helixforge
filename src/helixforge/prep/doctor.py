"""Preflight tool check: ``helixforge doctor``."""

from __future__ import annotations

import os
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from helixforge.mikado.run import mikado_version
from helixforge.prep._subprocess import tool_version

# --- the pinned matrix (single source of truth: docs/EXTERNAL_TOOLS.md) ---


@dataclass(frozen=True)
class ToolSpec:
    """One row of the pinned external-tool matrix."""

    key: str  # config key / override key, e.g. "mikado"
    default_bin: str  # bare command resolved on PATH by default
    stage: str  # "reconcile" | "prep" | "benchmark"
    pinned: str | None  # pinned version string, or None when not yet pinned
    required: bool  # a missing required tool ⇒ report not OK
    version_arg: str = "--version"


# Ordered by stage so the rendered table reads reconcile → prep → benchmark.
TOOL_MATRIX = (
    # reconcile core: required (the pipeline cannot run without these).
    # Pins are the versions the TAIR10 milestone runs actually used
    # (docs/EXTERNAL_TOOLS.md).
    ToolSpec("mikado", "mikado", "reconcile", "2.3.4", True),
    ToolSpec("diamond", "diamond", "reconcile", "2.1.16", True),
    ToolSpec("transdecoder", "TransDecoder.LongOrfs", "reconcile", "5.5.0", True),
    ToolSpec("portcullis", "portcullis", "reconcile", None, False),
    # prep: optional (only needed when staging evidence from raw reads)
    ToolSpec("star", "STAR", "prep", "2.7.11b", False),
    ToolSpec("hisat2", "hisat2", "prep", "2.2.1", False),
    ToolSpec("samtools", "samtools", "prep", "1.21", False),
    ToolSpec("stringtie", "stringtie", "prep", "3.0.3", False),
    ToolSpec("miniprot", "miniprot", "prep", "0.18", False),
    ToolSpec("helixer", "Helixer.py", "prep", "0.3.5", False),
    # ncRNA: optional (the opt-in structured-ncRNA hook; prep/ncrna.py)
    ToolSpec("trnascan", "tRNAscan-SE", "ncrna", "2.0.12", False),
    ToolSpec("infernal", "cmscan", "ncrna", "1.1.5", False, version_arg="-h"),
    # benchmark: optional (pins filled Phase 18; see docs/EXTERNAL_TOOLS.md)
    ToolSpec("mikado_compare", "mikado", "benchmark", "2.3.4", False),
    ToolSpec("gffcompare", "gffcompare", "benchmark", "0.12.6", False),
    ToolSpec("compleasm", "compleasm", "benchmark", "0.2.6", False),
    ToolSpec("busco", "busco", "benchmark", "5.7.1", False),
    ToolSpec("omark", "omark", "benchmark", "0.3.0", False),
    ToolSpec("agat", "agat_sp_statistics.pl", "benchmark", "1.4.1", False),
)

# Status tags.
FOUND = "found"
MISSING = "missing"
MISMATCH = "mismatch"


def _version_tuple(raw: str | None) -> tuple[int, ...] | None:
    """Pull a ``(major, minor[, patch])`` int tuple from a version string."""
    if not raw:
        return None
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", raw)
    if not m:
        return None
    return tuple(int(g) for g in m.groups() if g is not None)


def _matches_pin(pinned: str | None, detected_raw: str | None) -> bool:
    """Compare a detected version string against a pinned spec.

    Compares only as many components as the pin specifies (so pin ``2.1`` accepts
    detected ``2.1.11``). Returns ``True`` when they agree, ``False`` on a clear
    mismatch, and ``True`` when either side is unparseable (we do not flag a
    mismatch we cannot substantiate).
    """
    pin = _version_tuple(pinned)
    det = _version_tuple(detected_raw)
    if pin is None or det is None:
        return True
    n = min(len(pin), len(det))
    return pin[:n] == det[:n]


@dataclass
class ToolStatus:
    """The resolved state of one tool."""

    key: str
    bin: str
    stage: str
    required: bool
    pinned: str | None
    detected: str | None
    status: str


@dataclass
class DoctorReport:
    """The full preflight result; renders a found/missing/mismatch table."""

    statuses: list[ToolStatus] = field(default_factory=list)

    def ok(self) -> bool:
        """True iff no *required* tool is missing."""
        return not any(s.required and s.status == MISSING for s in self.statuses)

    def missing_required(self) -> list[ToolStatus]:
        return [s for s in self.statuses if s.required and s.status == MISSING]

    def drift_warnings(self) -> list[str]:
        """One warning per version-mismatch.

        A resolved tool whose version disagrees with the pinned matrix is a
        silent scientific-correctness risk, most acutely Mikado, whose
        scoring/config schema has changed across minor versions. The mismatch
        is a loud warning (not a hard failure: ``ok()`` ignores it), so a
        deliberate version bump is not blocked, only surfaced.
        """
        warnings: list[str] = []
        for s in self.statuses:
            if s.status != MISMATCH:
                continue
            extra = (
                ", Mikado's scoring/config schema changes across minors"
                "; re-verify the emitted config"
                if s.key in ("mikado", "mikado_compare")
                else ""
            )
            warnings.append(
                f"{s.key}: detected {s.detected!r} differs from pinned "
                f"{s.pinned!r}{extra}"
            )
        return warnings

    def render(self) -> str:
        """Render the report as a fixed-width text table."""
        headers = ("TOOL", "STAGE", "BIN", "REQ", "PINNED", "DETECTED", "STATUS")
        rows = [
            (
                s.key,
                s.stage,
                s.bin,
                "yes" if s.required else "no",
                s.pinned or "-",
                (s.detected or "-"),
                s.status,
            )
            for s in self.statuses
        ]
        widths = [
            max(len(headers[i]), *(len(r[i]) for r in rows))
            if rows
            else len(headers[i])
            for i in range(len(headers))
        ]

        def fmt(cols: tuple[str, ...]) -> str:
            return "  ".join(c.ljust(widths[i]) for i, c in enumerate(cols))

        lines = [fmt(headers), fmt(tuple("-" * w for w in widths))]
        lines.extend(fmt(r) for r in rows)
        for w in self.drift_warnings():
            lines.append(f"WARNING: version drift, {w}")
        n_missing = sum(1 for s in self.statuses if s.status == MISSING)
        n_mismatch = sum(1 for s in self.statuses if s.status == MISMATCH)
        n_found = sum(1 for s in self.statuses if s.status == FOUND)
        lines.append("")
        lines.append(
            f"{n_found} found, {n_missing} missing, {n_mismatch} version-mismatch; "
            f"required OK: {self.ok()}"
        )
        return "\n".join(lines)


def _resolve_bins(config_or_bins: object) -> dict[str, str]:
    """Map a tool key → bin name from a dict, an object with ``*_bin`` attrs, or None.

    Accepts: ``None`` (all defaults); a ``{key: bin}`` dict; or any object whose
    ``<key>_bin`` attribute overrides the default for that key (e.g. a
    ``PipelineConfig``). Unknown keys are ignored.
    """
    overrides: dict[str, str] = {}
    if config_or_bins is None:
        return overrides
    if isinstance(config_or_bins, dict):
        return dict(config_or_bins)
    for spec in TOOL_MATRIX:
        attr = f"{spec.key}_bin"
        val = getattr(config_or_bins, attr, None)
        if val:
            overrides[spec.key] = val
    return overrides


# ---------------------------------------------------------------------------
# Emitted-config schema verification
# ---------------------------------------------------------------------------
#
# mikado/config.py emits configuration.yaml + scoring YAML with the column
# order/schema *assumed*: and Mikado's scoring schema has changed across minor
# versions. A silent drift produces a config Mikado accepts but
# mis-interprets: a scientifically wrong run that still "passes". This check
# detects the Mikado version and asserts the emitted files carry the top-level
# sections that version expects (a heuristic, Mikado has no stable
# ``configure --check`` subcommand), warning loudly on drift.

# The pinned Mikado (single source: TOOL_MATRIX / docs/EXTERNAL_TOOLS.md).
_PINNED_MIKADO = next(s.pinned for s in TOOL_MATRIX if s.key == "mikado")

# Required top-level sections for the Mikado 2.x configuration + scoring schema.
_CONFIG_REQUIRED = ("reference", "prepare", "serialise", "pick")
_SCORING_REQUIRED = ("requirements", "scoring")


@dataclass
class ConfigVerification:
    """Result of :func:`verify_emitted_config`."""

    ok: bool
    version: tuple[int, ...] | None
    config_missing: list[str] = field(default_factory=list)
    scoring_missing: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def render(self) -> str:
        lines = [
            f"mikado version: "
            f"{'.'.join(map(str, self.version)) if self.version else 'undetected'}",
        ]
        if self.config_missing:
            lines.append(
                f"configuration.yaml missing sections: {', '.join(self.config_missing)}"
            )
        if self.scoring_missing:
            lines.append(
                f"scoring file missing sections: {', '.join(self.scoring_missing)}"
            )
        for w in self.warnings:
            lines.append(f"WARNING: {w}")
        lines.append(f"emitted-config schema OK: {self.ok}")
        return "\n".join(lines)


def _has_top_level_key(text: str, key: str) -> bool:
    """True if ``key:`` appears as a top-level (column-0) YAML mapping key."""
    return re.search(rf"(?m)^{re.escape(key)}\s*:", text) is not None


def verify_emitted_config(
    config_path: str | Path,
    scoring_path: str | Path | None = None,
    *,
    mikado_bin: str = "mikado",
    version: tuple[int, ...] | None = None,
) -> ConfigVerification:
    """Verify an emitted Mikado ``configuration.yaml`` (+ scoring) against the
    detected Mikado version's expected schema.

    Detects the Mikado version (``version`` overrides for testing), then asserts
    the configuration carries every required top-level section
    (:data:`_CONFIG_REQUIRED`) and, if ``scoring_path`` is given, that the scoring
    file carries :data:`_SCORING_REQUIRED`. A missing section ⇒ not OK. The
    detected version is compared to the pinned ``2.3.x``; a mismatch or an
    undetectable version is a loud **warning** (the schema may have drifted) but
    does not by itself fail the check. Never raises on a missing tool/file,
    those become warnings/missing-section findings.
    """
    if version is None:
        version = mikado_version(mikado_bin)

    warnings: list[str] = []
    if version is None:
        warnings.append(
            "could not detect Mikado version; assuming the pinned "
            f"{_PINNED_MIKADO} schema, verify manually before a real run"
        )
    elif not _matches_pin(_PINNED_MIKADO, ".".join(map(str, version))):
        warnings.append(
            f"detected Mikado {'.'.join(map(str, version))} differs from the "
            f"pinned {_PINNED_MIKADO}; the scoring/config schema may have drifted"
        )

    cfg_text = Path(config_path).read_text() if Path(config_path).exists() else ""
    if not cfg_text:
        warnings.append(f"configuration file not found: {config_path}")
    config_missing = [
        k for k in _CONFIG_REQUIRED if not _has_top_level_key(cfg_text, k)
    ]

    scoring_missing: list[str] = []
    if scoring_path is not None:
        sc_text = Path(scoring_path).read_text() if Path(scoring_path).exists() else ""
        if not sc_text:
            warnings.append(f"scoring file not found: {scoring_path}")
        scoring_missing = [
            k for k in _SCORING_REQUIRED if not _has_top_level_key(sc_text, k)
        ]

    ok = not config_missing and not scoring_missing and cfg_text != ""
    return ConfigVerification(
        ok=ok,
        version=version,
        config_missing=config_missing,
        scoring_missing=scoring_missing,
        warnings=warnings,
    )


def check_tools(
    config_or_bins: object = None,
    *,
    matrix: tuple[ToolSpec, ...] = TOOL_MATRIX,
) -> DoctorReport:
    """Probe every tool in ``matrix`` and return a :class:`DoctorReport`.

    ``config_or_bins`` supplies per-tool bin overrides (see :func:`_resolve_bins`).
    Each tool's version is probed via ``tool_version`` (never raises): a tool that
    does not resolve is ``missing``; one whose version disagrees with its pin is a
    ``mismatch``; otherwise ``found``.
    """
    overrides = _resolve_bins(config_or_bins)
    statuses: list[ToolStatus] = []
    for spec in matrix:
        bin_name = overrides.get(spec.key, spec.default_bin)
        detected = tool_version(bin_name, spec.version_arg)
        if detected is None:
            status = MISSING
        elif not _matches_pin(spec.pinned, detected):
            status = MISMATCH
        else:
            status = FOUND
        statuses.append(
            ToolStatus(
                key=spec.key,
                bin=bin_name,
                stage=spec.stage,
                required=spec.required,
                pinned=spec.pinned,
                detected=detected,
                status=status,
            )
        )
    return DoctorReport(statuses=statuses)


# --- environment hygiene checks ------------------------------


@dataclass
class ShimStatus:
    """Whether the imported ``helixforge`` matches the current repo checkout.

    The conda ``mikado-env`` console-script shim can point at a stale source tree
    so ``import helixforge`` resolves to foreign code (audit P2 / finding #3). We
    flag that: when the current directory *is* a helixforge checkout
    (``src/helixforge`` present) but the imported package lives elsewhere, the
    shim is foreign and the user is silently running the wrong code.
    """

    applicable: bool  # the cwd looks like a helixforge repo checkout
    foreign: bool  # imported package is NOT this checkout's src/helixforge
    imported: str  # resolved dir of the imported helixforge package
    expected: str  # the checkout's src/helixforge (when applicable)
    console_script: str | None  # resolved `helixforge` on PATH, if any

    def render(self) -> str:
        if not self.applicable:
            return (
                "shim check: n/a (run from a helixforge checkout to verify the "
                f"`helixforge` console script, resolved to {self.imported})"
            )
        if self.foreign:
            return (
                "shim check: FOREIGN, `import helixforge` resolves to "
                f"{self.imported}, not this checkout's {self.expected}. The "
                "`helixforge` console script points at a stale tree; "
                "`pip install -e .` into this env or call ~/.local/bin/helixforge."
            )
        return f"shim check: OK, helixforge imports from {self.imported}"


def check_shim(
    repo_root: str | os.PathLike[str] | None = None,
    *,
    module_file: str | os.PathLike[str] | None = None,
    which: object = None,
) -> ShimStatus:
    """Flag a mis-pointed ``helixforge`` console-script shim.

    ``module_file`` (default the live ``helixforge.__file__``) and ``which``
    (default :func:`shutil.which`) are injectable for testing. ``repo_root``
    defaults to the current working directory; the check only applies when that
    directory is a checkout (``src/helixforge`` exists).
    """
    import helixforge

    mod = Path(module_file or helixforge.__file__).resolve()
    pkg_dir = mod.parent
    root = Path(repo_root).resolve() if repo_root is not None else Path.cwd().resolve()
    expected = (root / "src" / "helixforge").resolve()
    applicable = expected.exists()
    foreign = applicable and pkg_dir != expected
    which_fn = which if which is not None else shutil.which
    console_script = which_fn("helixforge")  # type: ignore[operator]
    return ShimStatus(
        applicable=applicable,
        foreign=foreign,
        imported=str(pkg_dir),
        expected=str(expected),
        console_script=console_script,
    )


@dataclass
class CramRefStatus:
    """CRAM-reference readiness for the alignment inputs.

    A CRAM decodes only against its genome FASTA; with neither ``--reference``
    nor an htslib ``REF_CACHE`` / ``REF_PATH`` cache, pysam falls back to a
    remote ENA MD5 fetch that fails offline. ``ok`` is True when there are no
    CRAM inputs, or a reference/cache is available for them.
    """

    cram_paths: list[str]
    reference_set: bool
    ref_cache_set: bool

    @property
    def has_cram(self) -> bool:
        return bool(self.cram_paths)

    @property
    def ok(self) -> bool:
        return (not self.has_cram) or self.reference_set or self.ref_cache_set

    def render(self) -> str:
        if not self.has_cram:
            return "cram-reference: n/a (no CRAM inputs)"
        names = ", ".join(Path(p).name for p in self.cram_paths)
        if self.ok:
            via = "--reference" if self.reference_set else "REF_CACHE/REF_PATH"
            return f"cram-reference: OK ({names} → decode via {via})"
        return (
            f"cram-reference: MISSING, CRAM input(s) {names} but no --reference "
            "and no REF_CACHE/REF_PATH; pysam would attempt a remote ENA fetch "
            "(fails offline). Pass --reference <genome.fa> or seed a REF_CACHE."
        )


def check_cram_reference(
    bam_paths: Sequence[str | os.PathLike[str]] | None,
    *,
    reference: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
) -> CramRefStatus:
    """Report whether CRAM inputs have a usable offline reference."""
    from helixforge.io.bam import _is_cram

    env = env if env is not None else os.environ
    crams: list[str] = []
    for p in bam_paths or ():
        try:
            if _is_cram(p):
                crams.append(str(p))
        except OSError:  # unreadable path → defer to the opener's own error
            continue
    ref_cache = bool(env.get("REF_CACHE") or env.get("REF_PATH"))
    return CramRefStatus(
        cram_paths=crams,
        reference_set=reference is not None,
        ref_cache_set=ref_cache,
    )


# --- scoring-template availability check ---


@dataclass
class ScoringProfileStatus:
    """Whether the shipped scoring templates are installed and resolvable."""

    profiles: dict[str, bool]  # profile name → resolvable

    @property
    def ok(self) -> bool:
        return bool(self.profiles) and all(self.profiles.values())

    def render(self) -> str:
        if not self.profiles:
            return "scoring profiles: NONE FOUND, package data may not be installed"
        lines = ["scoring profiles:"]
        for name, available in sorted(self.profiles.items()):
            tag = "OK" if available else "MISSING"
            lines.append(f"  {name}: {tag}")
        return "\n".join(lines)


def check_scoring_profiles() -> ScoringProfileStatus:
    """Probe every shipped scoring profile for resolvability."""
    from helixforge.mikado.config import (
        list_scoring_profiles,
        validate_scoring_profile,
    )

    profiles: dict[str, bool] = {}
    for name in list_scoring_profiles():
        try:
            validate_scoring_profile(name)
            profiles[name] = True
        except (FileNotFoundError, ValueError):
            profiles[name] = False
    return ScoringProfileStatus(profiles=profiles)
