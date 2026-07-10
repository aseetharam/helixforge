"""Independence invariant for the standalone ``confidence`` scorer (§C2).

``score/confidence.py`` must depend on **only** its one evidence type (the
Helixer HDF5) plus pure I/O — never on Mikado, the reconcile pipeline, RNA-seq
I/O, or the sibling ``score/evidence`` module. This is asserted two ways:
statically (no forbidden import in the source AST) and dynamically (importing the
module in a clean subprocess pulls none of the forbidden modules), plus a runtime
check that scoring works with only an HDF5 + a GFF3 present.
"""

import ast
import subprocess
import sys
import textwrap
from pathlib import Path

import helixforge.score.confidence as confidence

FORBIDDEN = [
    "helixforge.mikado.run",
    "helixforge.reconcile.pipeline",
    "helixforge.io.bam",
    "helixforge.io.stringtie",
    "helixforge.score.evidence",
]


def _imported_modules(source_path):
    tree = ast.parse(Path(source_path).read_text())
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                names.add(node.module)
    return names


def test_source_has_no_forbidden_module_level_import():
    imported = _imported_modules(confidence.__file__)
    for mod in FORBIDDEN:
        assert mod not in imported, f"{mod} must not be imported by score/confidence"


def test_source_imports_only_allowed_helixforge_modules():
    imported = _imported_modules(confidence.__file__)
    hf = {m for m in imported if m.startswith("helixforge.")}
    allowed = {
        "helixforge.io.gff",
        "helixforge.io.hdf5",
        "helixforge.mikado.emit_external",
        "helixforge.utils.regions",
        "helixforge.viz.tracks",  # optional, lazy — only for the bigWig track output
    }
    assert hf <= allowed, f"unexpected helixforge imports: {hf - allowed}"


def test_clean_import_pulls_no_forbidden_module():
    code = textwrap.dedent(
        """
        import sys
        import helixforge.score.confidence  # noqa: F401
        forbidden = %r
        leaked = [m for m in forbidden if m in sys.modules]
        print(",".join(leaked))
        """
    ) % FORBIDDEN
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, check=True,
    )
    leaked = result.stdout.strip()
    assert leaked == "", f"forbidden modules imported on clean import: {leaked}"


def test_runs_with_only_hdf5_and_gff3_present(conf_gff3_path, conf_hdf5_path):
    # No BAM / StringTie / SJ / Mikado artifact anywhere — scoring must succeed.
    df = confidence.score_annotation(conf_gff3_path, conf_hdf5_path)
    assert len(df) == 3
    assert (df["helixer_support"] >= 0.0).all()


def test_no_mikado_or_toolchain_attribute_leak():
    # The module namespace must not expose pipeline/run/bam handles.
    for forbidden_attr in ("run_pipeline", "PipelineConfig", "run_mikado", "BAMReader"):
        assert not hasattr(confidence, forbidden_attr)
