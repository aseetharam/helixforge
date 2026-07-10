"""Phase 25 D4 — static-typing quality gate (M10).

Asserts ``src/helixforge`` is clean under ``mypy --strict`` (the milestone-M10
type gate) and that nobody silenced mypy wholesale: every ``# type: ignore`` in
the package must be *targeted* (an error code and/or a reason), never a blanket
``# type: ignore`` (What-NOT-to-do #3 in the phase prompt).

The strict check self-skips when mypy is not importable so a minimal test env
stays green; CI installs the ``dev`` extra (which pins mypy), so the gate runs
there. mypy reads the strict ``[tool.mypy]`` block in ``pyproject.toml``.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src" / "helixforge"


def _have_mypy() -> bool:
    try:
        import mypy  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(not _have_mypy(), reason="mypy not installed (dev extra)")
def test_mypy_strict_clean() -> None:
    """``mypy --strict src/helixforge`` exits 0 (reads the pyproject strict config)."""
    proc = subprocess.run(
        [sys.executable, "-m", "mypy", "--strict", "src/helixforge"],
        cwd=_REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, (
        "mypy --strict src/helixforge is not clean:\n"
        + proc.stdout[-4000:]
        + proc.stderr[-1000:]
    )


def test_no_blanket_type_ignores() -> None:
    """Every ``# type: ignore`` carries an error code or an inline reason.

    A bare ``# type: ignore`` (optionally with trailing whitespace / end-of-line)
    silences *all* errors on a line — exactly the blanket suppression the phase
    forbids. A targeted ignore looks like ``# type: ignore[union-attr]`` and/or
    is followed by a ``# reason`` comment.
    """
    bare = re.compile(r"#\s*type:\s*ignore\s*(?:#.*)?$")
    offenders = []
    for path in _SRC.rglob("*.py"):
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.rstrip()
            if "type: ignore" not in stripped:
                continue
            # Targeted forms carry a [code] right after "ignore".
            if re.search(r"#\s*type:\s*ignore\[", stripped):
                continue
            if bare.search(stripped):
                offenders.append(f"{path.relative_to(_REPO_ROOT)}:{lineno}: {stripped.strip()}")
    assert not offenders, "blanket '# type: ignore' (use # type: ignore[code]):\n" + "\n".join(offenders)
