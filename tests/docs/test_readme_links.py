"""README front-door integrity: links resolve, citation present, version agrees,
and no claim of an unpublished package.
"""

from __future__ import annotations

import re
from pathlib import Path

import helixforge

ROOT = Path(__file__).resolve().parent.parent.parent
README = ROOT / "README.md"

# Markdown links/images: ![alt](target) and [text](target).
_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


def _relative_targets(text: str) -> list[str]:
    targets = []
    for m in _LINK.finditer(text):
        tgt = m.group(1).strip()
        if tgt.startswith(("http://", "https://", "mailto:", "#")):
            continue
        if tgt.startswith("<") and tgt.endswith(">"):  # placeholder like <repo-url>
            continue
        targets.append(tgt.split("#", 1)[0])  # drop any anchor
    return targets


def test_internal_links_resolve():
    targets = _relative_targets(README.read_text())
    assert targets, "expected at least one relative link in README"
    missing = [t for t in targets if not (ROOT / t).exists()]
    assert not missing, f"README links to nonexistent paths: {missing}"


def test_schematic_image_present():
    """The method schematic the README embeds exists on disk."""
    assert "presentation/method_schematic.svg" in README.read_text()
    assert (ROOT / "presentation" / "method_schematic.svg").exists()


def test_citation_block_present():
    text = README.read_text()
    assert "## Citation" in text, "README must carry a Citation section"
    assert "ISMB 2026" in text, "citation must reference the accepted ISMB abstract"


def test_version_matches_package():
    """The version string in README agrees with helixforge.__version__."""
    version = helixforge.__version__
    assert version == "4.0.0.dev0"  # single source of truth (pyproject mirrors it)
    assert version in README.read_text(), (
        f"README must state the package version {version}"
    )


def test_no_unpublished_pip_install():
    """README must not claim an unpublished PyPI package."""
    assert "pip install helixforge" not in README.read_text()
