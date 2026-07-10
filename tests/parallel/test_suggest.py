"""Tests for parallel/suggest.py — granularity/resource heuristics (Phase 16 D4)."""

from helixforge.parallel.suggest import suggest_plan


def _write_fai(path, lengths):
    path.write_text("".join(f"{s}\t{n}\t0\t{n}\t{n + 1}\n" for s, n in lengths))


def test_suggest_chunk_count_by_size(tmp_path):
    fai = tmp_path / "g.fasta.fai"
    # 200 Mb in one scaffold → ~5 chunks at 40 Mb/chunk.
    _write_fai(fai, [("chr1", 200_000_000)])
    s = suggest_plan(str(fai), hpc_profile={"max_array_size": 1000})
    assert s.target_chunks == 5
    assert s.stats.genome_size == 200_000_000


def test_suggest_clamps_to_array_cap(tmp_path):
    fai = tmp_path / "g.fasta.fai"
    _write_fai(fai, [("chr1", 2_000_000_000)])  # would want 50 chunks
    s = suggest_plan(str(fai), hpc_profile={"max_array_size": 10})
    assert s.target_chunks == 10
    assert any("clamped" in note for note in s.rationale)


def test_suggest_resources_and_render(tmp_path):
    fai = tmp_path / "g.fasta.fai"
    _write_fai(fai, [("chr1", 120_000_000), ("chr2", 30_000_000)])
    s = suggest_plan(str(fai), hpc_profile={"cores_per_node": 16, "walltime_cap_hours": 12})
    assert s.threads == 8
    assert s.mem_gb >= 4
    assert s.walltime.count(":") == 2  # HH:MM:SS
    text = s.render()
    assert "HEURISTIC" in text
    assert "trade-offs" in text
    assert "N50" in text


def test_suggest_n50(tmp_path):
    fai = tmp_path / "g.fasta.fai"
    _write_fai(fai, [("a", 100), ("b", 100), ("c", 100), ("d", 100)])
    s = suggest_plan(str(fai))
    assert s.stats.scaffold_count == 4
    assert s.stats.n50 == 100
