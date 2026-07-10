"""HyperShell execution backend (subprocess-first).

HyperShell (https://github.com/hypershell/hypershell) is a distributed task
runner. This backend drives it as an **external CLI**, the same subprocess
contract every other external tool follows, so the core install needs no
``hypershell`` Python dependency (it is an optional ``helixforge[hpc]`` extra
that simply puts the ``hs`` command on PATH). Two entry points:

- :func:`run_hypershell`, expand the plan into the same executor-agnostic task
  file ``parallel tasks`` produces, then run it with ``hs cluster``. The commands
  are byte-identical to every other executor, so there is no new coordinate risk.
- :func:`write_hypershell_plan`, emit a HyperShell-native chunk plan JSON
  (targets the unreleased ``hs cluster --chunk-plan`` feature, hypershell#37) with
  per-chunk fields pre-resolved to CLI-ready strings, so HyperShell's ``{key}``
  substitution can never mis-handle coordinates. Dependency-free.

When HyperShell's ``hypershell.cluster`` Python API and the ``--chunk-plan``
ingestion land, the same public surface can switch to them internally.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Callable

from helixforge.parallel.taskgen import _chunk_variables, generate_task_file
from helixforge.parallel.tasks import default_reconcile_template
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.parallel.plan import Plan
    from helixforge.reconcile.pipeline import PipelineConfig

_log = get_logger(__name__)

# subprocess.run-compatible callable; injected in tests so no real `hs` is needed.
Runner = Callable[..., "subprocess.CompletedProcess[str]"]


def _hs_version(hs_bin: str, runner: Runner) -> str:
    """Best-effort ``hs --version`` for the run log (never raises)."""
    try:
        proc = runner(
            [hs_bin, "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
        return (proc.stdout or proc.stderr or "").strip() or "unknown"
    except (OSError, ValueError):
        return "unknown"


def run_hypershell(
    plan: Plan,
    base_config: PipelineConfig,
    out_prefix: str,
    *,
    num_tasks: int,
    hs_bin: str = "hs",
    helixforge_bin: str = "helixforge",
    output_dir: str | Path | None = None,
    tasks_path: str | Path | None = None,
    extra_hs_args: tuple[str, ...] | list[str] = (),
    _runner: Runner = subprocess.run,
) -> list[Path]:
    """Run every chunk of ``plan`` through ``hs cluster``; return chunk prefixes.

    Expands ``plan`` into an executor-agnostic task file (identical to
    ``parallel tasks``) and invokes ``hs cluster <tasks> --num-tasks <n>``.
    ``num_tasks`` bounds how many chunk commands run concurrently. Returns the
    per-chunk output prefixes (``<output_dir>/<chunk_id>``) for :func:`aggregate`.
    Raises ``RuntimeError`` on a nonzero ``hs`` exit (argv is included).
    """
    out_dir = Path(output_dir) if output_dir else Path(f"{out_prefix}.chunks")
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks_file = Path(tasks_path) if tasks_path else Path(f"{out_prefix}.tasks.txt")

    template = default_reconcile_template(base_config, helixforge_bin=helixforge_bin)
    task_file = generate_task_file(plan, template, tasks_file, output_dir=out_dir)

    _log.info("HyperShell %s → running %d tasks via `%s cluster`",
              _hs_version(hs_bin, _runner), task_file.n_tasks, hs_bin)

    argv = [
        hs_bin, "cluster", str(tasks_file),
        "--num-tasks", str(num_tasks),
        *extra_hs_args,
    ]
    try:
        proc = _runner(argv, text=True, check=False)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"`{hs_bin}` not found, install HyperShell (`pip install "
            f"helixforge[hpc]`) or point --hs-bin at it. argv: {argv}"
        ) from exc
    if proc.returncode != 0:
        tail = (getattr(proc, "stderr", None) or "").strip().splitlines()[-20:]
        raise RuntimeError(
            f"`hs cluster` exited {proc.returncode}. argv: {argv}"
            + ("\n" + "\n".join(tail) if tail else "")
        )

    return [out_dir / c.chunk_id for c in plan.chunks]


def write_hypershell_plan(
    plan: Plan,
    base_config: PipelineConfig,
    path: str | Path,
    *,
    helixforge_bin: str = "helixforge",
    output_dir: str = "chunks",
) -> tuple[Path, str]:
    """Emit a HyperShell-native chunk plan JSON + its command template.

    Targets ``hs cluster --chunk-plan <json> --template <tmpl>`` (hypershell#37,
    unreleased). Each chunk carries pre-resolved, CLI-ready fields (``chunk_id``,
    ``region``, ``id_start``, ``novel_start``, ``output_dir``) so HyperShell's
    ``{key}`` substitution produces the exact command the tasks.txt path would,
    HelixForge owns the coordinate conversion, not HyperShell. Single-region
    chunks only (packed multi-region chunks must use the tasks.txt path).

    Returns ``(json_path, template)``.
    """
    chunks_json: list[dict[str, str]] = []
    for c in plan.chunks:
        if len(c.regions) != 1:
            raise NotImplementedError(
                f"chunk {c.chunk_id} spans {len(c.regions)} regions; the HyperShell "
                "native plan supports single-region chunks only: partition with "
                "pack_small_scaffolds=False, or use `run_hypershell` (tasks.txt)."
            )
        v = _chunk_variables(c, Path(output_dir))
        chunks_json.append(
            {
                "chunk_id": v["chunk_id"],
                "region": v["region"],
                "id_start": v["id_start"],
                "novel_start": v["novel_start"],
                "output_dir": v["output_dir"],
            }
        )

    template = default_reconcile_template(base_config, helixforge_bin=helixforge_bin)
    payload = {
        "metadata": {
            "strategy": plan.strategy,
            "total_loci": plan.total_loci,
            "template": template,
        },
        "chunks": chunks_json,
    }
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    _log.info("wrote HyperShell-native plan (%d chunks) → %s", len(chunks_json), out)
    return out, template
