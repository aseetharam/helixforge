"""Phase 10 D2 — self-contained interactive HTML + index.

We assert the HTML is written, is self-contained (embeds its data, no external
CDN), and contains the expected gene id / isoform count; the index links the
per-gene pages. No pixel/DOM-render assertions. Concrete literal coordinates
only (CLAUDE.md §12).
"""

import json

from helixforge.reconcile.models import (
    ASEvent,
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.viz.interactive import interactive_gene, interactive_index


def make_tx(tid, exon_bounds=((1000, 1200), (1300, 1500)),
            cds=((1050, 1200, 0), (1300, 1450, 0)), is_primary=True, tpm=10.0):
    return TranscriptCandidate(
        transcript_id=tid,
        locus_id=tid.rsplit(".", 1)[0],
        source="mikado",
        seqid="chr1",
        start=exon_bounds[0][0],
        end=exon_bounds[-1][1],
        strand="+",
        exons=[Exon(s, e) for s, e in exon_bounds],
        cds=[CDSSegment(*c) for c in cds] if cds else None,
        tpm=tpm,
        junction_support_fraction=1.0,
        confidence=0.8,
        combined_score=15.0,
        is_primary=is_primary,
    )


def make_gene(gid="HFG_00001", txs=None, tier=1, origin="mikado_1to1"):
    txs = txs or [make_tx(f"{gid}.1")]
    return ReconciledGene(
        gene_id=gid,
        seqid="chr1",
        start=min(t.start for t in txs),
        end=max(t.end for t in txs),
        strand="+",
        tier=tier,
        transcripts=txs,
        primary_transcript_id=txs[0].transcript_id,
        classification=LocusClassification(locus_id=gid, status="EXPRESSED"),
        origin=origin,
        as_events=[ASEvent("ES", "chr1", 1300, 1500, "+", 20)],
    )


def test_interactive_gene_writes_file(tmp_path):
    out = interactive_gene(make_gene(), tmp_path / "g.html")
    assert out.exists()
    assert out.stat().st_size > 0


def test_interactive_gene_contains_gene_id(tmp_path):
    out = interactive_gene(make_gene("HFG_00042"), tmp_path / "g.html")
    assert "HFG_00042" in out.read_text()


def test_interactive_gene_contains_isoform_count(tmp_path):
    gene = make_gene(txs=[make_tx("HFG_00001.1"),
                          make_tx("HFG_00001.2", is_primary=False)])
    text = interactive_gene(gene, tmp_path / "g.html").read_text()
    assert "isoforms: <b>2</b>" in text


def test_interactive_gene_is_self_contained(tmp_path):
    text = interactive_gene(make_gene(), tmp_path / "g.html").read_text()
    # Embedded JSON data, no external script/style CDN references.
    assert 'id="gene-data"' in text
    assert "http://" not in text.replace('"http://www.w3.org/2000/svg"', "")
    assert "https://" not in text


def test_interactive_gene_embeds_valid_json(tmp_path):
    text = interactive_gene(make_gene("HFG_00007"), tmp_path / "g.html").read_text()
    start = text.index('type="application/json">') + len('type="application/json">')
    end = text.index("</script>", start)
    record = json.loads(text[start:end])
    assert record["gene_id"] == "HFG_00007"
    assert record["transcripts"][0]["exons"][0] == {"start": 1000, "end": 1200}


def test_interactive_gene_shows_as_events(tmp_path):
    text = interactive_gene(make_gene(), tmp_path / "g.html").read_text()
    assert '"kind": "ES"' in text


def test_interactive_index_links_pages(tmp_path):
    genes = [make_gene("HFG_00001"), make_gene("HFG_00002", tier=2)]
    index = interactive_index(genes, tmp_path)
    assert index.name == "index.html"
    text = index.read_text()
    assert 'href="HFG_00001.html"' in text
    assert 'href="HFG_00002.html"' in text
    assert (tmp_path / "HFG_00001.html").exists()
    assert (tmp_path / "HFG_00002.html").exists()


def test_interactive_index_has_sort_columns(tmp_path):
    index = interactive_index([make_gene()], tmp_path)
    text = index.read_text()
    # Sortable by tier / isoforms / AED.
    assert "tier" in text
    assert "isoforms" in text
    assert "AED" in text
    assert "sortBy(" in text
