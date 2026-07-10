"""Phase 23 D4 — error-boundary hardening. Floor: 3.

The benchmark table must distinguish "tool absent" from "tool failed"; the
narrowed before_after catch must degrade an *expected* tool failure to NaN while
letting a genuine code bug propagate; and run_tool must stream the full stderr to
a per-step log while keeping only the tail in the exception.
"""

import math

import pytest

from helixforge.bench import wrappers
from helixforge.bench.wrappers import BenchmarkError, benchmark_all
from helixforge.prep._subprocess import run_tool


# --- benchmark_all status column: absent vs failed --------------------------

def test_benchmark_status_distinguishes_absent_from_failed(tmp_path, monkeypatch):
    monkeypatch.setattr(wrappers, "run_agat_stats",
                        lambda *a, **k: {"number_of_gene": 30000.0})

    def absent(*a, **k):
        raise BenchmarkError("compleasm not found", kind="tool-absent")

    def failed(*a, **k):
        raise BenchmarkError("busco exit 1", kind="failed")

    monkeypatch.setattr(wrappers, "run_compleasm", absent)
    monkeypatch.setattr(wrappers, "run_busco", failed)

    df = benchmark_all("pred.gff3", "proteins.fa", tmp_path / "bench",
                       lineage="brassicales_odb10")
    assert "status" in df.columns
    assert df[df.tool == "compleasm"]["status"].tolist() == ["tool-absent"]
    assert df[df.tool == "busco"]["status"].tolist() == ["failed"]
    # agat succeeded: its metric rows are status "ok".
    assert set(df[df.tool == "agat"]["status"]) == {"ok"}
    # the failed/absent tools contribute no metric rows, only a status row.
    assert df[(df.tool == "compleasm") & (df.metric != "status")].empty


# --- before_after narrowed except: degrade vs propagate ---------------------

def test_completeness_degrades_to_nan_on_tool_failure(monkeypatch):
    from helixforge.stats import before_after as ba

    monkeypatch.setattr(ba, "_is_gff3_input", lambda x: False)

    def raise_benchmark(fa, out):
        raise BenchmarkError("busco missing", kind="tool-absent")

    # genome present + a non-gff3 input → the runner is invoked and fails.
    monkeypatch.setattr("helixforge.export.writers.write_protein_fasta",
                        lambda *a, **k: None)
    val = ba._proteome_metric(["gene"], genome=object(),
                              runner=raise_benchmark, pick=lambda r: r)
    assert math.isnan(val)


def test_completeness_propagates_real_bug(monkeypatch):
    from helixforge.stats import before_after as ba

    monkeypatch.setattr(ba, "_is_gff3_input", lambda x: False)
    monkeypatch.setattr("helixforge.export.writers.write_protein_fasta",
                        lambda *a, **k: None)

    def real_bug(fa, out):
        raise TypeError("a genuine code bug, not a NaN")

    with pytest.raises(TypeError):
        ba._proteome_metric(["gene"], genome=object(),
                            runner=real_bug, pick=lambda r: r)


# --- run_tool streams full stderr to a per-step log -------------------------

def test_run_tool_streams_full_stderr_to_log(tmp_path, monkeypatch):
    import subprocess

    long_stderr = "HEAD-cause\n" + ("x" * 5000) + "\nTAIL-line"
    log_path = tmp_path / "step.stderr.log"

    def fake_run(argv, **kwargs):
        raise subprocess.CalledProcessError(2, argv, stderr=long_stderr)

    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(RuntimeError) as ei:
        run_tool(["faketool", "--go"], stderr_log=log_path)

    # full stderr (head + tail) is on disk; the exception keeps only the tail.
    written = log_path.read_text()
    assert "HEAD-cause" in written and "TAIL-line" in written
    assert "HEAD-cause" not in str(ei.value)   # head dropped from the message
    assert "TAIL-line" in str(ei.value)        # tail retained
