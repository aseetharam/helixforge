"""Independence invariant for the standalone ``evidence`` scorer (§C2).

``score/evidence.py`` must depend on **only** its RNA-seq evidence I/O plus the
shared concordance/region helpers — never on Mikado, the reconcile pipeline, the
Helixer HDF5 reader, or the sibling ``score/confidence`` module. Asserted both
statically (no forbidden import in the source AST) and dynamically (a clean
subprocess import pulls none of them), plus a runtime check that scoring works
with only BAM/SJ + a GFF3 present.
"""

import ast
import subprocess
import sys
import textwrap
from pathlib import Path

import helixforge.score.evidence as evidence

FORBIDDEN = [
    "helixforge.mikado.run",
    "helixforge.reconcile.pipeline",
    "helixforge.io.hdf5",
    "helixforge.score.confidence",
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
    imported = _imported_modules(evidence.__file__)
    for mod in FORBIDDEN:
        assert mod not in imported, f"{mod} must not be imported by score/evidence"


def test_source_imports_only_allowed_helixforge_modules():
    imported = _imported_modules(evidence.__file__)
    hf = {m for m in imported if m.startswith("helixforge.")}
    allowed = {
        "helixforge.io.bam",
        "helixforge.io.stringtie",
        "helixforge.io.gff",
        "helixforge.stats.evidence_concordance",
        "helixforge.utils.regions",
        # Protein-AED axis (poster phase): the miniprot parser + alignment
        # wrapper are evidence I/O, the protein analogue of io.bam — NOT the
        # forbidden pipeline modules. reconcile.models is a TYPE_CHECKING-only
        # import for the MiniprotAlignment type. The C2 invariant (no Mikado /
        # reconcile pipeline / HDF5 / sibling scorer; runs with one evidence
        # type) is still enforced by FORBIDDEN below.
        "helixforge.io.miniprot",
        "helixforge.prep.protein_align",
        "helixforge.reconcile.models",
    }
    assert hf <= allowed, f"unexpected helixforge imports: {hf - allowed}"


def test_clean_import_pulls_no_forbidden_module():
    code = textwrap.dedent(
        """
        import sys
        import helixforge.score.evidence  # noqa: F401
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


def test_runs_with_only_bam_and_gff3_present(ev_gff3_path, ev_bam_path):
    # No HDF5 / Mikado / StringTie anywhere — scoring must succeed.
    df = evidence.score_annotation(ev_gff3_path, bam_paths=[ev_bam_path])
    assert len(df) == 5
    assert (df["num_introns"] >= 0).all()


def test_runs_with_only_star_sj_and_gff3_present(ev_gff3_path, ev_star_sj_path):
    df = evidence.score_annotation(ev_gff3_path, star_sj_paths=[ev_star_sj_path])
    assert len(df) == 5


def test_no_mikado_or_hdf5_attribute_leak():
    for forbidden_attr in ("run_pipeline", "PipelineConfig", "HDF5ConfidenceReader"):
        assert not hasattr(evidence, forbidden_attr)
