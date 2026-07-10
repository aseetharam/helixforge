"""Phase 11 D1 — benchmarking wrappers.

External tools are mocked (no real subprocess); parsers are exercised on tiny
synthetic stats outputs. argv is asserted to be a list with the right flags.
"""

import pytest

from helixforge.bench import wrappers
from helixforge.bench.wrappers import (
    BenchmarkError,
    benchmark_all,
    parse_agat_stats,
    parse_busco_summary,
    parse_compleasm_summary,
    parse_gffcompare_stats,
    parse_mikado_compare_stats,
    parse_omark_summary,
    run_agat_stats,
    run_busco,
    run_compleasm,
    run_gffcompare,
    run_mikado_compare,
    run_omark,
)

# ---------------------------------------------------------------------------
# Synthetic tool outputs
# ---------------------------------------------------------------------------

_MIKADO_STATS = """\
Command line: mikado compare -r ref.gff3 -p pred.gff3 -o cmp
    7 reference RNAs in 7 reference genes
    7 predicted RNAs in 7 predicted genes
--------------------------------- |   Sn |   Pr |   F1 |
                      Base level: 95.50  90.20  92.80
          Exon level (stringent): 88.00  84.00  85.95
                    Intron level: 99.00  97.00  97.99
              Intron chain level: 80.00  78.00  78.99
                Transcript level: 70.00  72.00  70.99
                      Gene level: 65.00  68.00  66.46
"""

_GFFCOMPARE_STATS = """\
# gffcompare v0.12.6
#= Summary for dataset: pred.gtf
        Base level:    95.5     |    90.2    |
        Exon level:    88.0     |    84.0    |
   Intron level:   99.0     |    97.0    |
Intron chain level:   80.0     |    78.0    |
  Transcript level:   70.0     |    72.0    |
       Locus level:    65.0     |    68.0    |
"""

_COMPLEASM_SUMMARY = """\
## lineage: brassicales_odb10
S:97.65%, 249
D:0.78%, 2
F:0.39%, 1
I:0.00%, 0
M:1.18%, 3
N:255
"""

_BUSCO_SUMMARY = """\
# BUSCO version is: 5.4.7
	C:97.8%[S:96.5%,D:1.3%],F:0.7%,M:1.5%,n:255
	249	Complete BUSCOs (C)
"""

_OMARK_SUMMARY = """\
The selected closest lineage is Brassicales.
Single:90.00% (13042)
Duplicated:5.00% (725)
Missing:5.00% (725)
Consistent:85.00% (29750)
Inconsistent:10.00% (3500)
Contaminants:0.00% (0)
Unknown:5.00% (1750)
"""

_AGAT_STATS = """\
--------------------------------------------------------------
Number of gene                              30000
Number of mrna                              35000
Number of cds                               34000
Number of exon                              150000
mean cds length (bp)                        1234.5
mean gene length (bp)                       2500
"""


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def test_parse_mikado_compare_stats(tmp_path):
    p = tmp_path / "cmp.stats"
    p.write_text(_MIKADO_STATS)
    out = parse_mikado_compare_stats(p)
    assert out["base"] == {"sn": 95.5, "pr": 90.2, "f1": 92.8}
    assert out["transcript"]["f1"] == 70.99
    assert out["intron_chain"]["sn"] == 80.0
    assert "gene" in out


def test_parse_gffcompare_stats(tmp_path):
    p = tmp_path / "cmp.stats"
    p.write_text(_GFFCOMPARE_STATS)
    out = parse_gffcompare_stats(p)
    assert out["transcript"] == {"sn": 70.0, "pr": 72.0}
    assert out["base"]["sn"] == 95.5
    assert out["locus"]["pr"] == 68.0


def test_parse_compleasm_summary(tmp_path):
    p = tmp_path / "summary.txt"
    p.write_text(_COMPLEASM_SUMMARY)
    out = parse_compleasm_summary(p)
    assert out["single"] == 97.65
    assert out["duplicated"] == 0.78
    assert out["missing"] == 1.18
    assert out["n"] == 255
    assert out["complete"] == pytest.approx(98.43)


def test_parse_busco_summary(tmp_path):
    p = tmp_path / "short_summary.txt"
    p.write_text(_BUSCO_SUMMARY)
    out = parse_busco_summary(p)
    assert out["complete"] == 97.8
    assert out["single"] == 96.5
    assert out["duplicated"] == 1.3
    assert out["missing"] == 1.5
    assert out["n"] == 255


def test_parse_omark_summary(tmp_path):
    p = tmp_path / "proteins.sum"
    p.write_text(_OMARK_SUMMARY)
    out = parse_omark_summary(p)
    assert out["single"] == 90.0
    assert out["consistent"] == 85.0
    assert out["inconsistent"] == 10.0
    assert out["complete"] == pytest.approx(95.0)


def test_parse_agat_stats(tmp_path):
    p = tmp_path / "agat.txt"
    p.write_text(_AGAT_STATS)
    out = parse_agat_stats(p)
    assert out["number_of_gene"] == 30000.0
    assert out["number_of_mrna"] == 35000.0
    assert out["mean_cds_length_bp"] == 1234.5


# ---------------------------------------------------------------------------
# Run wrappers — argv assembly (mocked subprocess) + parse round-trip
# ---------------------------------------------------------------------------


@pytest.fixture
def captured_argv(monkeypatch):
    calls = []

    def fake_run_tool(argv, step, cwd=None):
        calls.append((step, list(argv)))
        return None

    monkeypatch.setattr(wrappers, "_run_tool", fake_run_tool)
    return calls


def test_run_mikado_compare_argv_and_parse(tmp_path, captured_argv):
    (tmp_path / "cmp.stats").write_text(_MIKADO_STATS)
    out = run_mikado_compare("ref.gff3", "pred.gff3", str(tmp_path / "cmp"))
    step, argv = captured_argv[0]
    assert step == "mikado.compare"
    assert argv[:2] == ["mikado", "compare"]
    assert "-r" in argv and "ref.gff3" in argv
    assert "-p" in argv and "pred.gff3" in argv
    assert out["base"]["sn"] == 95.5


def test_run_gffcompare_argv_and_parse(tmp_path, captured_argv):
    (tmp_path / "gc.stats").write_text(_GFFCOMPARE_STATS)
    out = run_gffcompare("ref.gff3", "pred.gtf", str(tmp_path / "gc"))
    step, argv = captured_argv[0]
    assert step == "gffcompare"
    assert argv[0] == "gffcompare"
    assert argv[-1] == "pred.gtf"
    assert out["transcript"]["sn"] == 70.0


def test_run_compleasm_argv_and_parse(tmp_path, captured_argv):
    out_dir = tmp_path / "cl"
    out_dir.mkdir()
    (out_dir / "summary.txt").write_text(_COMPLEASM_SUMMARY)
    out = run_compleasm("proteins.fa", "brassicales_odb10", out_dir, threads=8)
    step, argv = captured_argv[0]
    assert step == "compleasm"
    assert argv[:2] == ["compleasm", "protein"]
    assert 8 in argv  # threads passed through (stringified later by _run_tool)
    assert out["complete"] == pytest.approx(98.43)


def test_run_busco_argv_and_parse(tmp_path, captured_argv):
    out_dir = tmp_path / "busco_run"
    out_dir.mkdir()
    (out_dir / "short_summary.specific.txt").write_text(_BUSCO_SUMMARY)
    out = run_busco("proteins.fa", "brassicales_odb10", out_dir)
    step, argv = captured_argv[0]
    assert step == "busco"
    assert "-m" in argv and "proteins" in argv
    assert out["complete"] == 97.8


def test_run_omark_argv_and_parse(tmp_path, captured_argv):
    out_dir = tmp_path / "omark_run"
    out_dir.mkdir()
    (out_dir / "proteins.sum").write_text(_OMARK_SUMMARY)
    out = run_omark("proteins.fa", "LUCA.h5", out_dir)
    step, argv = captured_argv[0]
    assert step == "omark"
    assert "-d" in argv and "LUCA.h5" in argv
    assert out["consistent"] == 85.0


def test_run_agat_stats_argv_and_parse(tmp_path, captured_argv):
    out_path = tmp_path / "agat.txt"
    out_path.write_text(_AGAT_STATS)
    out = run_agat_stats("annotation.gff3", out_path)
    step, argv = captured_argv[0]
    assert step == "agat.statistics"
    assert argv[0] == "agat_sp_statistics.pl"
    assert "--gff" in argv and "annotation.gff3" in argv
    assert out["number_of_gene"] == 30000.0


# ---------------------------------------------------------------------------
# _run_tool error handling
# ---------------------------------------------------------------------------


def test_run_tool_missing_binary_raises(monkeypatch):
    def boom(*a, **k):
        raise FileNotFoundError("no such tool")

    monkeypatch.setattr(wrappers.subprocess, "run", boom)
    with pytest.raises(BenchmarkError, match="not found"):
        wrappers._run_tool(["nope", "--help"], "test.step")


def test_run_tool_nonzero_exit_raises(monkeypatch):
    import subprocess as sp

    def boom(*a, **k):
        raise sp.CalledProcessError(2, a[0], stderr="kaboom")

    monkeypatch.setattr(wrappers.subprocess, "run", boom)
    with pytest.raises(BenchmarkError, match="exit 2"):
        wrappers._run_tool(["tool"], "test.step")


# ---------------------------------------------------------------------------
# benchmark_all — merges applicable tools into one table; skips failures
# ---------------------------------------------------------------------------


def test_benchmark_all_merges_tools(tmp_path, monkeypatch):
    monkeypatch.setattr(wrappers, "run_agat_stats",
                        lambda *a, **k: {"number_of_gene": 30000.0})
    monkeypatch.setattr(wrappers, "run_mikado_compare",
                        lambda *a, **k: {"transcript": {"f1": 70.99}})
    monkeypatch.setattr(wrappers, "run_gffcompare",
                        lambda *a, **k: {"transcript": {"sn": 70.0, "pr": 72.0}})
    monkeypatch.setattr(wrappers, "run_compleasm",
                        lambda *a, **k: {"complete": 98.43})
    monkeypatch.setattr(wrappers, "run_busco",
                        lambda *a, **k: {"complete": 97.8})
    monkeypatch.setattr(wrappers, "run_omark",
                        lambda *a, **k: {"consistent": 85.0})

    df = benchmark_all(
        "pred.gff3", "proteins.fa", tmp_path / "bench",
        reference_gff3="ref.gff3", lineage="brassicales_odb10", omadb="LUCA.h5",
    )
    tools = set(df["tool"])
    assert tools == {"agat", "mikado_compare", "gffcompare", "compleasm", "busco", "omark"}
    # nested dict flattened to namespaced metric
    row = df[(df.tool == "mikado_compare") & (df.metric == "transcript_f1")]
    assert float(row["value"].iloc[0]) == 70.99


def test_benchmark_all_skips_failing_tool(tmp_path, monkeypatch):
    monkeypatch.setattr(wrappers, "run_agat_stats",
                        lambda *a, **k: {"number_of_gene": 30000.0})

    def fail(*a, **k):
        raise BenchmarkError("compleasm not installed")

    monkeypatch.setattr(wrappers, "run_compleasm", fail)
    monkeypatch.setattr(wrappers, "run_busco", lambda *a, **k: {"complete": 97.8})

    df = benchmark_all("pred.gff3", "proteins.fa", tmp_path / "bench",
                       lineage="brassicales_odb10")
    # Phase 23 D4: a failing tool no longer silently vanishes — it gets a single
    # status-metric row, while its real metrics are still absent.
    metric_tools = set(df[df.metric != "status"]["tool"])
    assert "compleasm" not in metric_tools    # failed → no metric rows
    assert {"agat", "busco"} <= metric_tools   # the rest still computed
    assert df[df.tool == "compleasm"]["status"].tolist() == ["failed"]
    assert set(df[df.tool == "agat"]["status"]) == {"ok"}


def test_benchmark_all_only_agat_without_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(wrappers, "run_agat_stats",
                        lambda *a, **k: {"number_of_gene": 30000.0})
    df = benchmark_all("pred.gff3", None, tmp_path / "bench")
    assert set(df["tool"]) == {"agat"}
