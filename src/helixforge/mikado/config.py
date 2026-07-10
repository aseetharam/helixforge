"""Mikado config emitters: list.txt, configuration.yaml, scoring."""

from __future__ import annotations

import importlib.resources as _resources
import shutil
import subprocess
from pathlib import Path
from typing import Any

# Default AS class codes admitted as alternative isoforms (validate per version).
DEFAULT_VALID_CCODES = ("j", "J", "C", "c", "=", "n", "h", "g", "G")

_SCORING_PKG = _resources.files("helixforge").joinpath("data").joinpath("scoring")


def _scoring_dir() -> Any:
    """Return the importlib.resources Traversable for the scoring package data."""
    return _SCORING_PKG


def list_scoring_profiles() -> list[str]:
    """Return the names of all shipped scoring profiles (e.g. ``['strict', 'permissive']``)."""
    prefix = "helixforge.plant."
    profiles: set[str] = set()
    try:
        for entry in _scoring_dir().iterdir():
            name = entry.name
            if name.startswith(prefix) and (
                name.endswith(".yaml") or name.endswith(".yaml.template")
            ):
                core = name[len(prefix) :]
                core = core.removesuffix(".template").removesuffix(".yaml")
                if core:
                    profiles.add(core)
    except (OSError, AttributeError):
        pass
    return sorted(profiles)


def _read_scoring_template(profile_name: str) -> tuple[str, str]:
    """Resolve a scoring profile's text + filename (``helixforge.plant.<p>.yaml``).

    Resolves via ``importlib.resources`` (works in an installed wheel and a
    source checkout with ``pip install -e .``). Returns ``(text, base_filename)``
    or raises ``FileNotFoundError``.
    """
    valid = list_scoring_profiles()
    if profile_name not in valid:
        raise FileNotFoundError(
            f"unknown scoring profile {profile_name!r} "
            f"(valid profiles: {', '.join(valid) or 'none found'})"
        )
    base = f"helixforge.plant.{profile_name}.yaml"
    for name in (base, base + ".template"):
        try:
            res = _scoring_dir().joinpath(name)
            if res.is_file():
                return res.read_text(), base
        except (FileNotFoundError, ModuleNotFoundError, AttributeError, OSError):
            pass
    raise FileNotFoundError(
        f"scoring profile {profile_name!r} files not installed "
        f"(expected helixforge/data/scoring/{base}[.template] in the package)"
    )


def _yaml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    return str(value)


def _emit_yaml(obj: Any, indent: int = 0) -> list[str]:
    """Minimal YAML emitter for dict/list/scalar trees (deterministic order)."""
    pad = "  " * indent
    lines: list[str] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(value, dict):
                lines.append(f"{pad}{key}:")
                lines.extend(_emit_yaml(value, indent + 1))
            elif isinstance(value, (list, tuple)):
                lines.append(f"{pad}{key}:")
                for item in value:
                    lines.append(f"{pad}- {_yaml_scalar(item)}")
            else:
                lines.append(f"{pad}{key}: {_yaml_scalar(value)}")
    else:
        lines.append(f"{pad}{_yaml_scalar(obj)}")
    return lines


def write_input_list(entries: list[dict[str, Any]], out_path: str | Path) -> Path:
    """Write Mikado ``list.txt``; return ``Path``.

    Each entry is a dict with keys ``file``, ``label`` and optional
    ``strand_specific`` (default True), ``score`` (0), ``is_reference`` (False),
    ``exclude_redundant`` (False). The **Helixer** entry must set
    ``is_reference=True``.

    Assumed tab column order (confirm against ``mikado configure --help`` for the
    pinned version): file, label, strand_specific, score, is_reference,
    exclude_redundant.
    """
    out_path = Path(out_path)
    with out_path.open("w") as fh:
        for e in entries:
            cols = [
                str(e["file"]),
                str(e["label"]),
                "True" if e.get("strand_specific", True) else "False",
                str(e.get("score", 0)),
                "True" if e.get("is_reference", False) else "False",
                "True" if e.get("exclude_redundant", False) else "False",
            ]
            fh.write("\t".join(cols) + "\n")
    return out_path


def write_configuration(
    genome_fa: str | Path,
    list_path: str | Path,
    scoring_path: str | Path,
    junctions_path: str | Path,
    out_path: str | Path,
    as_report: bool = True,
    only_confirmed_introns: bool = True,
    valid_ccodes: tuple[str, ...] = DEFAULT_VALID_CCODES,
    min_cds_overlap: float = 0.6,
    min_cdna_overlap: float = 0.6,
    max_isoforms: int = 5,
    keep_retained_introns: bool = False,
    pad: bool = True,
    chimera_split: bool = True,
    flank: int = 200,
    use_subprocess: bool = False,
    mikado_bin: str = "mikado",
    **pick_overrides: Any,
) -> Path:
    """Write (or generate via ``mikado configure``) ``configuration.yaml``.

    By default a configuration YAML is templated directly (testable without
    Mikado). With ``use_subprocess=True`` ``mikado configure`` is invoked instead
    (preferred for a real run; mocked in tests). The ``pick.alternative_splicing``
    knobs are the isoform-quality controls; each is a parameter and
    any may be overridden via ``pick_overrides`` (e.g. ``max_isoforms=10``).
    """
    out_path = Path(out_path)

    if use_subprocess:
        if shutil.which(mikado_bin) is None:
            raise FileNotFoundError(
                f"mikado binary not found: {mikado_bin!r}. Install Mikado or use "
                "the template path (use_subprocess=False)."
            )
        out_path = out_path.resolve()
        argv = [
            mikado_bin,
            "configure",
            "--list",
            str(Path(list_path).resolve()),
            "--reference",
            str(Path(genome_fa).resolve()),
            "--scoring",
            str(Path(scoring_path).resolve()),
            "--junctions",
            str(Path(junctions_path).resolve()),
            str(out_path),
        ]
        subprocess.run(argv, check=True, cwd=out_path.parent)
        return out_path

    alt_splicing: dict[str, Any] = {
        "report": as_report,
        "only_confirmed_introns": only_confirmed_introns,
        "valid_ccodes": list(valid_ccodes),
        "min_cds_overlap": min_cds_overlap,
        "min_cdna_overlap": min_cdna_overlap,
        "max_isoforms": max_isoforms,
        "keep_retained_introns": keep_retained_introns,
        "pad": pad,
    }
    # apply overrides for any alternative_splicing knob
    for key, value in pick_overrides.items():
        alt_splicing[key] = value

    config: dict[str, Any] = {
        "reference": {"genome": str(genome_fa)},
        "prepare": {"files": {"list": str(list_path)}},
        "serialise": {"files": {"junctions": [str(junctions_path)]}},
        "pick": {
            "scoring_file": str(scoring_path),
            "alternative_splicing": alt_splicing,
            "chimera_split": {"execute": bool(chimera_split)},
            "clustering": {"flank": flank},
        },
    }
    out_path.write_text("\n".join(_emit_yaml(config)) + "\n")
    return out_path


def validate_scoring_profile(profile_name: str) -> None:
    """Resolve + validate a scoring profile eagerly (fail-fast).

    Call at pipeline startup, before any expensive stage. Raises
    ``FileNotFoundError`` (unknown / not installed) or ``ValueError``
    (malformed template).
    """
    text, base = _read_scoring_template(profile_name)
    for section in ("requirements:", "scoring:"):
        if section not in text:
            raise ValueError(
                f"scoring profile {base} missing section {section!r} "
                "(schema must be validated against the pinned Mikado)"
            )


def install_scoring_profile(
    profile_name: str = "strict", dest_dir: str | Path | None = None
) -> Path:
    """Copy + lightly validate a scoring profile; return ``Path``.

    Resolves ``helixforge.plant.<profile>.yaml`` (or ``...yaml.template``) from
    the shipped package data via ``importlib.resources``. Validation checks the
    file has the expected top-level sections — the full schema must still be
    validated against the **pinned** Mikado version (see
    :func:`helixforge.prep.doctor.verify_emitted_config`).
    """
    text, base = _read_scoring_template(profile_name)
    for section in ("requirements:", "scoring:"):
        if section not in text:
            raise ValueError(
                f"scoring profile {base} missing section {section!r} "
                "(schema must be validated against the pinned Mikado)"
            )

    dest_dir = Path(dest_dir) if dest_dir is not None else Path.cwd()
    dest = dest_dir / base
    dest.write_text(text)
    return dest
