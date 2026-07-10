"""Tests for the standalone ``helixforge confidence`` distribution summary.

Covers :func:`summarize_confidence_scores`, :func:`write_confidence_summary_tsv`
and the histogram + empirical-CDF distribution plot (no pie chart).

The ``overall_score`` weights (mean_prob 0.3, min_prob 0.2, boundary 0.2,
coding 0.15, worst_exon 0.15) sum to 1.0, so a gene whose every component equals
``v`` has ``overall_score == v`` exactly — letting us assert percentiles against
concrete literals rather than computed values.
"""

from pathlib import Path

import pytest

from helixforge.core.confidence import (
    GeneConfidence,
    summarize_confidence_scores,
    write_confidence_summary_tsv,
)
from helixforge.core.plots import plot_confidence_distribution


def _gene(gid: str, v: float, strand: str = "+", cls: str = "medium") -> GeneConfidence:
    """A gene whose overall_score == ``v`` exactly (all components = v)."""
    return GeneConfidence(
        gene_id=gid,
        seqid="chr1",
        start=100,
        end=400,
        strand=strand,
        mean_prob=v,
        min_prob=v,
        median_prob=v,
        entropy=v,
        boundary_sharpness=v,
        coding_consistency=v,
        worst_exon_score=v,
        exon_scores=[v, v],
        confidence_class=cls,
    )


@pytest.fixture
def scores() -> list[GeneConfidence]:
    # overall scores: 0.2, 0.4, 0.6, 0.8 (both strands represented)
    return [
        _gene("g1", 0.2, strand="+", cls="low"),
        _gene("g2", 0.4, strand="-", cls="low"),
        _gene("g3", 0.6, strand="+", cls="low"),
        _gene("g4", 0.8, strand="-", cls="medium"),
    ]


class TestSummaryStats:
    def test_percentiles_and_basics(self, scores) -> None:
        d = dict(summarize_confidence_scores(scores))
        assert d["n"] == 4
        assert d["mean"] == pytest.approx(0.5)
        assert d["median"] == pytest.approx(0.5)
        assert d["min"] == pytest.approx(0.2)
        assert d["max"] == pytest.approx(0.8)
        assert d["std"] == pytest.approx(0.2236068)
        # numpy linear-interpolation percentiles over [0.2,0.4,0.6,0.8]
        assert d["p5"] == pytest.approx(0.23)
        assert d["p25"] == pytest.approx(0.35)
        assert d["p50"] == pytest.approx(0.5)
        assert d["p75"] == pytest.approx(0.65)
        assert d["p95"] == pytest.approx(0.77)

    def test_component_subscores_present(self, scores) -> None:
        d = dict(summarize_confidence_scores(scores))
        # Each component reports mean/median/std; means equal the overall mean
        # here because every component == overall in the fixture.
        for label in ("entropy", "boundary_sharpness", "cds_consistency", "per_exon"):
            assert d[f"{label}_mean"] == pytest.approx(0.5)
            assert f"{label}_median" in d
            assert f"{label}_std" in d

    def test_class_counts_preserved(self, scores) -> None:
        d = dict(summarize_confidence_scores(scores))
        assert d["n_high"] == 0
        assert d["n_medium"] == 1
        assert d["n_low"] == 3

    def test_ordered_distribution_first(self, scores) -> None:
        rows = summarize_confidence_scores(scores)
        keys = [k for k, _ in rows]
        assert keys[0] == "n"
        # distribution stats come before class counts
        assert keys.index("mean") < keys.index("n_high")
        assert keys.index("p95") < keys.index("entropy_mean")

    def test_empty(self) -> None:
        d = dict(summarize_confidence_scores([]))
        assert d["n"] == 0
        for key in ("mean", "median", "std", "min", "max", "p5", "p95"):
            assert d[key] is None
        assert d["n_high"] == d["n_medium"] == d["n_low"] == 0

    def test_write_tsv_long_format(self, scores, tmp_path: Path) -> None:
        out = tmp_path / "stats.tsv"
        write_confidence_summary_tsv(summarize_confidence_scores(scores), out)
        lines = out.read_text().splitlines()
        assert lines[0] == "metric\tvalue"
        body = dict(line.split("\t") for line in lines[1:])
        assert body["n"] == "4"
        assert body["mean"] == "0.500000"
        assert body["p25"] == "0.350000"
        assert body["n_low"] == "3"

    def test_write_tsv_empty_blank_values(self, tmp_path: Path) -> None:
        out = tmp_path / "empty.tsv"
        write_confidence_summary_tsv(summarize_confidence_scores([]), out)
        body = dict(line.split("\t") for line in out.read_text().splitlines()[1:])
        assert body["n"] == "0"
        assert body["mean"] == ""  # None -> empty


class TestDistributionPlot:
    def test_matplotlib_has_two_panels_no_pie(self, scores, tmp_path: Path) -> None:
        import matplotlib
        matplotlib.use("Agg")

        out = tmp_path / "dist.png"
        fig = plot_confidence_distribution(scores, out, format="png")
        assert out.exists()
        assert len(fig.axes) == 2
        titles = [ax.get_title() for ax in fig.axes]
        assert any("CDF" in t for t in titles)
        # no pie: a pie axis would carry wedge patches, the CDF/hist axes don't
        import matplotlib.patches as mpatches
        wedges = [
            p for ax in fig.axes for p in ax.patches
            if isinstance(p, mpatches.Wedge)
        ]
        assert wedges == []

    def test_plotly_histogram_and_cdf_no_pie(self, scores, tmp_path: Path) -> None:
        pytest.importorskip("plotly")
        out = tmp_path / "dist.html"
        fig = plot_confidence_distribution(scores, out, format="html")
        assert out.exists()
        types = {t.type for t in fig.data}
        assert "pie" not in types
        assert "histogram" in types
        assert "scatter" in types  # the empirical CDF
