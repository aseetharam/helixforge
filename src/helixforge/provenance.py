"""Run-provenance preamble for self-describing GFF3 output."""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Mapping
from typing import Any

import attrs

from helixforge.prep._subprocess import tool_version

# External tools whose versions are worth recording in a HelixForge run's
# provenance. Keys are the provenance labels; the pipeline
# resolves each to the configured binary before probing.
PROVENANCE_TOOLS = ("mikado", "transdecoder", "diamond", "miniprot", "star")


def helixforge_version() -> str:
    """Return the installed HelixForge version (best-effort, never raises)."""
    try:
        import importlib.metadata as _md

        return _md.version("helixforge")
    except Exception:  # pragma: no cover - fallback for an uninstalled tree
        try:
            import helixforge

            return str(getattr(helixforge, "__version__", "unknown"))
        except Exception:
            return "unknown"


def file_md5(path: str | os.PathLike[str], *, chunk_size: int = 1 << 20) -> str | None:
    """Streaming MD5 hex digest of ``path`` (``None`` if the file is absent).

    Reads in ``chunk_size`` blocks so a multi-Gb genome FASTA never loads into
    memory. Returns ``None`` for a missing/unreadable file rather than raising,
    provenance is best-effort and must never sink a finished run.
    """
    if not path or not os.path.exists(path):
        return None
    h = hashlib.md5()
    try:
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(chunk_size), b""):
                h.update(block)
    except OSError:
        return None
    return h.hexdigest()


def param_hash(params: Mapping[str, Any] | None) -> str:
    """Deterministic short hash of the resolved run parameters.

    Hashes a canonical (sorted-key) JSON encoding of ``params`` so the same
    configuration always yields the same digest; non-JSON values fall back to
    ``str``. Returns the first 16 hex chars of the SHA-256 (``"none"`` for an
    empty/absent mapping).
    """
    if not params:
        return "none"
    encoded = json.dumps(params, sort_keys=True, default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def collect_tool_versions(
    tool_bins: Mapping[str, str | os.PathLike[str] | None],
) -> dict[str, str | None]:
    """Probe ``{label: binary}`` → ``{label: version-string-or-None}`` (best-effort).

    Each binary is resolved via :func:`prep._subprocess.tool_version` (runs
    ``<bin> --version``); a missing tool or nonzero exit yields ``None`` rather
    than raising, so provenance never blocks output. A ``None`` binary is skipped.
    """
    out: dict[str, str | None] = {}
    for label, binary in tool_bins.items():
        if binary is None:
            continue
        out[label] = tool_version(str(binary))
    return out


@attrs.frozen
class Provenance:
    """A resolved run-provenance record, renderable as a GFF3 ``#!`` preamble."""

    helixforge_version: str
    param_hash: str
    tool_versions: dict[str, str | None] = attrs.field(factory=dict)
    input_md5s: dict[str, str | None] = attrs.field(factory=dict)

    def to_gff3_lines(self) -> list[str]:
        """Render the provenance as a list of GFF3 ``#!`` comment lines.

        Deterministic ordering (version, param hash, then sorted tools, then
        sorted inputs) so the same run reproduces the same preamble byte-for-byte.
        """
        lines = [
            f"#!helixforge-version {self.helixforge_version}",
            f"#!param-hash {self.param_hash}",
        ]
        for label in sorted(self.tool_versions):
            ver = self.tool_versions[label]
            lines.append(f"#!tool {label}={ver if ver else 'unresolved'}")
        for label in sorted(self.input_md5s):
            md5 = self.input_md5s[label]
            lines.append(f"#!input {label} md5={md5 if md5 else 'NA'}")
        return lines


def build_provenance(
    *,
    params: Mapping[str, Any] | None = None,
    tool_bins: Mapping[str, str | os.PathLike[str] | None] | None = None,
    input_files: Mapping[str, str | os.PathLike[str] | None] | None = None,
) -> Provenance:
    """Assemble a :class:`Provenance` from resolved params, tool bins, and inputs.

    ``params`` is the resolved configuration (hashed); ``tool_bins`` maps a tool
    label to its configured binary (probed for a real version); ``input_files``
    maps an input label to a path (MD5-summed). All probes are best-effort,
    every value degrades to ``None``/``unresolved`` rather than raising, so a
    provenance build never fails a finished run.
    """
    versions = collect_tool_versions(tool_bins or {})
    md5s: dict[str, str | None] = {}
    for label, path in (input_files or {}).items():
        if path is None:
            continue
        md5s[label] = file_md5(path)
    return Provenance(
        helixforge_version=helixforge_version(),
        param_hash=param_hash(params),
        tool_versions=versions,
        input_md5s=md5s,
    )
