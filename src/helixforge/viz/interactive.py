"""Self-contained interactive HTML per gene + an index."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

from helixforge.export.writers import build_gene_records
from helixforge.stats.evidence_concordance import compute_aed

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene

# ---------------------------------------------------------------------------
# Per-gene page
# ---------------------------------------------------------------------------

_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HelixForge — {gene_id}</title>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 1.5rem; color: #222; }}
 h1 {{ font-size: 1.2rem; }}
 .meta span {{ display: inline-block; margin-right: 1rem; color: #555; }}
 #plot {{ border: 1px solid #ddd; width: 100%; height: 360px; }}
 .exon {{ fill: #cccccc; stroke: #888; }}
 .cds {{ fill: #2166ac; stroke: #1a4f86; }}
 .backbone {{ stroke: #555; stroke-width: 1; }}
 .as {{ fill: #b2182b; font-size: 10px; }}
 .lbl {{ font-size: 10px; fill: #333; }}
 #tip {{ position: fixed; pointer-events: none; background: #000a; color: #fff;
         padding: 4px 6px; border-radius: 3px; font-size: 11px; display: none;
         max-width: 320px; }}
 a {{ color: #2166ac; }}
</style>
</head>
<body>
<p><a href="index.html">&larr; index</a></p>
<h1>{gene_id}</h1>
<div class="meta">
 <span>seqid: <b>{seqid}</b></span>
 <span>strand: <b>{strand}</b></span>
 <span>tier: <b>{tier}</b></span>
 <span>origin: <b>{origin}</b></span>
 <span>isoforms: <b>{num_isoforms}</b></span>
</div>
<svg id="plot" preserveAspectRatio="xMinYMin meet"></svg>
<div id="tip"></div>
<script id="gene-data" type="application/json">{gene_json}</script>
<script>
const GENE = JSON.parse(document.getElementById("gene-data").textContent);
const svg = document.getElementById("plot");
const tip = document.getElementById("tip");
const NS = "http://www.w3.org/2000/svg";
const W = svg.clientWidth || 900, H = 360, PADX = 60, PADY = 30;
svg.setAttribute("viewBox", `0 0 ${{W}} ${{H}}`);

function rect(x, y, w, h, cls, title) {{
  const r = document.createElementNS(NS, "rect");
  r.setAttribute("x", x); r.setAttribute("y", y);
  r.setAttribute("width", Math.max(1, w)); r.setAttribute("height", h);
  r.setAttribute("class", cls);
  if (title) {{
    r.addEventListener("mousemove", e => {{
      tip.style.display = "block"; tip.style.left = (e.clientX + 12) + "px";
      tip.style.top = (e.clientY + 12) + "px"; tip.innerHTML = title;
    }});
    r.addEventListener("mouseleave", () => tip.style.display = "none");
  }}
  svg.appendChild(r);
}}
function line(x1, y1, x2, y2, cls) {{
  const l = document.createElementNS(NS, "line");
  l.setAttribute("x1", x1); l.setAttribute("y1", y1);
  l.setAttribute("x2", x2); l.setAttribute("y2", y2);
  l.setAttribute("class", cls); svg.appendChild(l);
}}
function text(x, y, s, cls) {{
  const t = document.createElementNS(NS, "text");
  t.setAttribute("x", x); t.setAttribute("y", y); t.setAttribute("class", cls);
  t.textContent = s; svg.appendChild(t);
}}

const txs = GENE.transcripts;
const lo = GENE.start, hi = GENE.end, span = (hi - lo) || 1;
const sx = x => PADX + (x - lo) / span * (W - 2 * PADX);
const laneH = Math.min(60, (H - 2 * PADY) / Math.max(1, txs.length));

txs.forEach((t, i) => {{
  const y = PADY + i * laneH + laneH / 2;
  line(sx(t.start), y, sx(t.end), y, "backbone");
  text(PADX, y - laneH / 2 + 10,
       t.transcript_id + (t.is_primary ? " *" : ""), "lbl");
  t.exons.forEach(e => rect(sx(e.start), y - 6, sx(e.end) - sx(e.start), 12,
      "exon", null));
  (t.cds || []).forEach(c => {{
    const tt = `<b>${{t.transcript_id}}</b><br>tpm: ${{t.tpm}}<br>`
      + `junction support: ${{t.junction_support}}<br>`
      + `helixer support: ${{t.helixer_support}}<br>`
      + `combined score: ${{t.combined_score}}<br>`
      + `homology: ${{t.has_homology}}`;
    rect(sx(c.start), y - 10, sx(c.end) - sx(c.start), 20, "cds", tt);
  }});
}});
(GENE.as_events || []).forEach(ev => {{
  const x = sx((ev.start + ev.end) / 2);
  text(x, H - 8, ev.kind, "as");
}});
</script>
</body>
</html>
"""


def interactive_gene(
    gene: ReconciledGene,
    out_path: str | Path,
    **evidence: Any,
) -> Path:
    """Write one self-contained interactive HTML file for ``gene``; return Path.

    ``**evidence`` is accepted for API symmetry with the static plot (future
    in-page tracks); the data layer is the Phase-9 per-gene JSON record.
    """
    record = build_gene_records([gene])[0]
    html = _PAGE_TEMPLATE.format(
        gene_id=record["gene_id"],
        seqid=record["seqid"],
        strand=record["strand"],
        tier=record["tier"],
        origin=record["origin"],
        num_isoforms=record["num_isoforms"],
        gene_json=json.dumps(record),
    )
    out_path = Path(out_path)
    out_path.write_text(html)
    return out_path


# ---------------------------------------------------------------------------
# Index page
# ---------------------------------------------------------------------------

_INDEX_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HelixForge — gene index</title>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 1.5rem; color: #222; }}
 table {{ border-collapse: collapse; width: 100%; }}
 th, td {{ border: 1px solid #ddd; padding: 4px 8px; text-align: left; }}
 th {{ cursor: pointer; background: #f4f4f4; }}
 a {{ color: #2166ac; }}
</style>
</head>
<body>
<h1>HelixForge gene index ({n} genes)</h1>
<table id="t">
<thead><tr>
 <th onclick="sortBy(0,false)">gene</th>
 <th onclick="sortBy(1,true)">tier</th>
 <th onclick="sortBy(2,true)">isoforms</th>
 <th onclick="sortBy(3,true)">AED</th>
 <th onclick="sortBy(4,false)">origin</th>
</tr></thead>
<tbody>
{rows}
</tbody>
</table>
<script>
function sortBy(col, numeric) {{
  const tb = document.querySelector("#t tbody");
  const rows = Array.from(tb.rows);
  const dir = tb.getAttribute("data-dir") === "asc" ? -1 : 1;
  tb.setAttribute("data-dir", dir === 1 ? "asc" : "desc");
  rows.sort((a, b) => {{
    let x = a.cells[col].textContent, y = b.cells[col].textContent;
    if (numeric) {{ x = parseFloat(x) || 0; y = parseFloat(y) || 0; return (x - y) * dir; }}
    return x.localeCompare(y) * dir;
  }});
  rows.forEach(r => tb.appendChild(r));
}}
</script>
</body>
</html>
"""

_INDEX_ROW = (
    '<tr><td><a href="{href}">{gene_id}</a></td><td>{tier}</td>'
    "<td>{num_isoforms}</td><td>{aed}</td><td>{origin}</td></tr>"
)


def _gene_aed(gene: ReconciledGene) -> float | None:
    """Composite AED of the primary transcript (reuses the Phase-9 helper)."""
    primary = next(
        (t for t in gene.transcripts if t.transcript_id == gene.primary_transcript_id),
        gene.transcripts[0],
    )
    expr = None if primary.tpm is None else (1.0 if primary.tpm > 0 else 0.0)
    return compute_aed(primary.junction_support_fraction, expr, None)


def interactive_index(
    genes: list[ReconciledGene],
    out_dir: str | Path,
) -> Path:
    """Write per-gene pages + a sortable ``index.html``; return the index Path.

    The index is sortable by tier / isoform count / AED (client-side, no server).
    Each row links to that gene's self-contained page in the same directory.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for gene in genes:
        href = f"{gene.gene_id}.html"
        interactive_gene(gene, out_dir / href)
        aed = _gene_aed(gene)
        rows.append(
            _INDEX_ROW.format(
                href=href,
                gene_id=gene.gene_id,
                tier=gene.tier,
                num_isoforms=len(gene.transcripts),
                aed="" if aed is None else f"{aed:.3f}",
                origin=gene.origin,
            )
        )

    index_path = out_dir / "index.html"
    index_path.write_text(_INDEX_TEMPLATE.format(n=len(genes), rows="\n".join(rows)))
    return index_path
