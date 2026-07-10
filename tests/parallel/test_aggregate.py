"""Tests for parallel/aggregate.py — gather + verify (Phase 16 D3).

Floor: 10. Hand-built per-chunk output files (GFF3 + tier files + report +
id_map); aggregate must merge them, verify global HFG uniqueness + once-only
locus coverage, and raise loudly on any violation. Both strands present.
"""

import json

import pytest

from helixforge.parallel.aggregate import aggregate, prefixes_from_pattern

_HEADER = "gene_id\tseqid\tstart\tend\tstrand\ttier\torigin\tflags"


def _write_chunk(prefix, genes):
    """Write one chunk's {gff3, tier1/2/3.gff3, report.tsv, id_map.json}.

    ``genes``: list of (gene_id, seqid, start_internal, end_internal, strand,
    tier, origin, locus_id). Coordinates are internal 0-based half-open.
    """
    prefix = str(prefix)

    def _gff(path, subset):
        lines = ["##gff-version 3"]
        for gid, seqid, s, e, strand, tier, origin, _lid in subset:
            lines.append(
                f"{seqid}\tHelixForge\tgene\t{s + 1}\t{e}\t.\t{strand}\t.\t"
                f"ID={gid};tier={tier};origin={origin}"
            )
            lines.append(
                f"{seqid}\tHelixForge\tmRNA\t{s + 1}\t{e}\t.\t{strand}\t.\t"
                f"ID={gid}.1;Parent={gid}"
            )
            lines.append("###")
        path_obj = open(path, "w")
        path_obj.write("\n".join(lines) + "\n")
        path_obj.close()

    _gff(f"{prefix}.gff3", genes)
    for n in (1, 2, 3):
        _gff(f"{prefix}.tier{n}.gff3", [g for g in genes if g[5] <= n])

    with open(f"{prefix}.report.tsv", "w") as fh:
        fh.write(_HEADER + "\n")
        for gid, seqid, s, e, strand, tier, origin, _lid in genes:
            fh.write(f"{gid}\t{seqid}\t{s}\t{e}\t{strand}\t{tier}\t{origin}\t\n")

    id_map = {lid: gid for gid, _seqid, _s, _e, _st, _t, _o, lid in genes}
    with open(f"{prefix}.id_map.json", "w") as fh:
        json.dump(id_map, fh)


@pytest.fixture
def two_chunks(tmp_path):
    a = tmp_path / "run.chunk_0000"
    b = tmp_path / "run.chunk_0001"
    _write_chunk(a, [
        ("HFG_00001", "chr1", 0, 500, "+", 1, "mikado_1to1", "G1"),
        ("HFG_00002", "chr1", 6000, 6500, "+", 3, "helixer_backstop", "G2"),
    ])
    _write_chunk(b, [
        ("HFG_00003", "chr2", 1000, 1500, "-", 1, "mikado_1to1", "G6"),
    ])
    return tmp_path, [str(a), str(b)]


# --------------------------------------------------------------------------
# Happy path
# --------------------------------------------------------------------------

def test_aggregate_gene_count(two_chunks):
    tmp_path, prefixes = two_chunks
    result = aggregate(prefixes, str(tmp_path / "merged"))
    assert result.num_genes == 3


def test_aggregate_merged_gff3_sorted(two_chunks):
    tmp_path, prefixes = two_chunks
    result = aggregate(prefixes, str(tmp_path / "merged"))
    text = result.gff3_path.read_text()
    gene_lines = [ln for ln in text.splitlines() if "\tgene\t" in ln]
    order = [ln.split("ID=")[1].split(";")[0] for ln in gene_lines]
    assert order == ["HFG_00001", "HFG_00002", "HFG_00003"]
    assert text.startswith("##gff-version 3")


def test_aggregate_report_union(two_chunks):
    tmp_path, prefixes = two_chunks
    result = aggregate(prefixes, str(tmp_path / "merged"))
    rows = result.report_path.read_text().strip().splitlines()
    assert rows[0] == _HEADER
    assert len(rows) == 1 + 3  # header + 3 genes


def test_aggregate_tier_counts(two_chunks):
    tmp_path, prefixes = two_chunks
    result = aggregate(prefixes, str(tmp_path / "merged"))
    assert result.tier_counts == {"1": 2, "3": 1}


def test_aggregate_origin_counts(two_chunks):
    tmp_path, prefixes = two_chunks
    result = aggregate(prefixes, str(tmp_path / "merged"))
    assert result.origin_counts == {"helixer_backstop": 1, "mikado_1to1": 2}


def test_aggregate_tier1_file_filtered(two_chunks):
    tmp_path, prefixes = two_chunks
    result = aggregate(prefixes, str(tmp_path / "merged"))
    tier1 = result.tier_paths[1].read_text()
    assert "HFG_00001" in tier1
    assert "HFG_00003" in tier1
    assert "HFG_00002" not in tier1  # tier 3 excluded from tier1 file


def test_aggregate_num_loci(two_chunks):
    tmp_path, prefixes = two_chunks
    result = aggregate(prefixes, str(tmp_path / "merged"))
    assert result.num_loci == 3


def test_aggregate_master_id_map_union(two_chunks):
    tmp_path, prefixes = two_chunks
    master = tmp_path / "master.id_map.json"
    aggregate(prefixes, str(tmp_path / "merged"), master_id_map_path=str(master))
    loaded = json.loads(master.read_text())
    assert loaded == {"G1": "HFG_00001", "G2": "HFG_00002", "G6": "HFG_00003"}


# --------------------------------------------------------------------------
# Rerun stability
# --------------------------------------------------------------------------

def test_aggregate_rerun_stable(two_chunks):
    tmp_path, prefixes = two_chunks
    master = tmp_path / "master.id_map.json"
    r1 = aggregate(prefixes, str(tmp_path / "merged"), master_id_map_path=str(master))
    first = json.loads(master.read_text())
    # Same loci ⇒ same HFGs on a second full cycle (master fed back in).
    r2 = aggregate(prefixes, str(tmp_path / "merged2"), master_id_map_path=str(master))
    assert json.loads(master.read_text()) == first
    assert r1.num_genes == r2.num_genes == 3


# --------------------------------------------------------------------------
# Failure modes — must raise, never silently merge
# --------------------------------------------------------------------------

def test_aggregate_raises_duplicate_hfg(tmp_path):
    a = tmp_path / "c0"
    b = tmp_path / "c1"
    _write_chunk(a, [("HFG_00001", "chr1", 0, 500, "+", 1, "mikado_1to1", "G1")])
    _write_chunk(b, [("HFG_00001", "chr2", 0, 500, "+", 1, "mikado_1to1", "G2")])
    with pytest.raises(ValueError, match="duplicate gene id"):
        aggregate([str(a), str(b)], str(tmp_path / "merged"))


def test_aggregate_raises_locus_covered_twice(tmp_path):
    a = tmp_path / "c0"
    b = tmp_path / "c1"
    _write_chunk(a, [("HFG_00001", "chr1", 0, 500, "+", 1, "mikado_1to1", "G1")])
    # distinct gene id (passes gff merge) but the SAME locus id in its id_map.
    _write_chunk(b, [("HFG_00009", "chr2", 0, 500, "+", 1, "mikado_1to1", "G1")])
    with pytest.raises(ValueError, match="covered by >1 chunk"):
        aggregate([str(a), str(b)], str(tmp_path / "merged"))


def test_aggregate_raises_master_conflict(two_chunks):
    tmp_path, prefixes = two_chunks
    master = tmp_path / "master.id_map.json"
    master.write_text(json.dumps({"G1": "HFG_09999"}))  # G1 already mapped elsewhere
    with pytest.raises(ValueError, match="id_map conflict"):
        aggregate(prefixes, str(tmp_path / "merged"), master_id_map_path=str(master))


def test_aggregate_empty_raises(tmp_path):
    with pytest.raises(ValueError, match="no chunk outputs"):
        aggregate([], str(tmp_path / "merged"))


def test_aggregate_raises_duplicate_transcript(tmp_path):
    # Distinct gene ids but a colliding transcript id (mRNA ID) → fail closed.
    a = tmp_path / "c0"
    b = tmp_path / "c1"
    _write_chunk(a, [("HFG_00001", "chr1", 0, 500, "+", 1, "mikado_1to1", "G1")])
    _write_chunk(b, [("HFG_00001", "chr2", 0, 500, "-", 1, "mikado_1to1", "G2")])
    # _write_chunk derives mRNA id as "<gene>.1"; same gene id → same transcript id.
    with pytest.raises(ValueError, match="duplicate (gene|transcript) id"):
        aggregate([str(a), str(b)], str(tmp_path / "merged"))


# --------------------------------------------------------------------------
# v1-style pattern discovery (--input-dir / --pattern)
# --------------------------------------------------------------------------

def test_prefixes_from_pattern_dedups_tier_files(two_chunks):
    tmp_path, _prefixes = two_chunks
    # '*.gff3' matches main + tier gff3s; prefixes collapse to one per chunk.
    found = prefixes_from_pattern(tmp_path, "*.gff3")
    assert found == [str(tmp_path / "run.chunk_0000"), str(tmp_path / "run.chunk_0001")]


def test_pattern_then_verifying_merge(two_chunks):
    tmp_path, _prefixes = two_chunks
    prefixes = prefixes_from_pattern(tmp_path, "*.gff3")
    result = aggregate(prefixes, str(tmp_path / "merged"))
    assert result.num_genes == 3
    assert result.num_loci == 3


def test_prefixes_from_pattern_no_match_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="no files match"):
        prefixes_from_pattern(tmp_path, "*.gff3")
