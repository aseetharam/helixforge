"""Genome-browser track export."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Any

from helixforge.utils.atomic import atomic_write

if TYPE_CHECKING:
    from helixforge.io.hdf5 import HDF5ConfidenceReader
    from helixforge.reconcile.models import ReconciledGene, SpliceJunction

_TRACKS_HINT = "bigWig export requires pyBigWig: pip install 'helixforge[tracks]'"

# itemRgb per tier (BED column 9), highest confidence greenest.
_TIER_RGB = {1: "27,120,55", 2: "127,191,123", 3: "217,179,101", 4: "179,88,6"}
_TIER_SCORE = {1: 1000, 2: 750, 3: 500, 4: 250}


# ---------------------------------------------------------------------------
# BED12 isoform models
# ---------------------------------------------------------------------------


def _bed12_line(gene: ReconciledGene, tx: Any) -> str:
    """One BED12 record for a transcript (thick = CDS, exons as blocks)."""
    chrom_start = tx.start
    chrom_end = tx.end
    if tx.cds:
        thick_start = min(c.start for c in tx.cds)
        thick_end = max(c.end for c in tx.cds)
    else:
        thick_start = thick_end = chrom_start  # non-coding convention

    exons = sorted(tx.exons, key=lambda e: e.start)
    block_sizes = ",".join(str(e.end - e.start) for e in exons) + ","
    block_starts = ",".join(str(e.start - chrom_start) for e in exons) + ","

    return "\t".join(
        str(v)
        for v in (
            tx.seqid,
            chrom_start,
            chrom_end,
            tx.transcript_id,
            _TIER_SCORE.get(gene.tier, 0),
            tx.strand,
            thick_start,
            thick_end,
            _TIER_RGB.get(gene.tier, "0,0,0"),
            len(exons),
            block_sizes,
            block_starts,
        )
    )


def write_bed12(genes: list[ReconciledGene], path: str | Path) -> Path:
    """Write a BED12 of every isoform (sorted by genomic position); return Path.

    Standards-valid for IGV/JBrowse2/UCSC: ``thickStart``/``thickEnd`` mark the
    CDS, exons are blocks, both strands handled (strand is column 6 only — block
    order is always genomic low→high).
    """
    rows = []
    for gene in genes:
        for tx in gene.transcripts:
            rows.append((tx.seqid, tx.start, _bed12_line(gene, tx)))
    rows.sort(key=lambda r: (r[0], r[1]))
    path = Path(path)
    with atomic_write(path) as fh:
        fh.write("".join(line + "\n" for _, _, line in rows))
    return path


def write_bigbed(
    genes: list[ReconciledGene],
    path: str | Path,
    chrom_sizes: dict[str, int],
    bedtobigbed_bin: str = "bedToBigBed",
) -> Path:
    """Convert the BED12 to bigBed via UCSC ``bedToBigBed``; return Path.

    ``chrom_sizes`` is ``{seqid: length}`` (written to a temp ``.sizes`` file).
    Raises ``RuntimeError`` with a clear hint if ``bedToBigBed`` is not on PATH —
    bigBed has no pure-Python writer.
    """
    path = Path(path)
    bed_path = path.with_suffix(".bed")
    sizes_path = path.with_suffix(".sizes")
    write_bed12(genes, bed_path)
    # bedToBigBed requires a coordinate-sorted BED (write_bed12 already sorts).
    sizes_path.write_text(
        "".join(f"{s}\t{n}\n" for s, n in sorted(chrom_sizes.items()))
    )
    try:
        subprocess.run(
            [bedtobigbed_bin, str(bed_path), str(sizes_path), str(path)],
            check=True,
            capture_output=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"{bedtobigbed_bin!r} not found — install UCSC bedToBigBed to write bigBed"
        ) from exc
    except subprocess.CalledProcessError as exc:  # pragma: no cover - tool-specific
        raise RuntimeError(
            f"bedToBigBed failed: {exc.stderr.decode(errors='replace')}"
        ) from exc
    return path


# ---------------------------------------------------------------------------
# Junction BED (Portcullis / regtools anchor style)
# ---------------------------------------------------------------------------


def write_junction_bed(junctions: list[SpliceJunction], path: str | Path) -> Path:
    """Write splice junctions as a BED12 (read count = score); return Path.

    Anchor style (two 1-bp blocks at the donor and acceptor) so IGV/JBrowse draw
    the intron as an arc. ``chromStart=donor``, ``chromEnd=acceptor`` (both
    0-based half-open, the internal convention). Score is the read count capped
    at the BED max of 1000.
    """
    path = Path(path)
    rows = sorted(junctions, key=lambda j: (j.seqid, j.donor, j.acceptor))
    lines = []
    for i, j in enumerate(rows, start=1):
        span = j.acceptor - j.donor
        line = "\t".join(
            str(v)
            for v in (
                j.seqid,
                j.donor,
                j.acceptor,
                f"JUNC{i:08d}",
                min(1000, j.read_count),
                j.strand,
                j.donor,
                j.donor,
                "0,0,0",
                2,
                "1,1,",
                f"0,{span - 1},",
            )
        )
        lines.append(line)
    with atomic_write(path) as fh:
        fh.write("".join(line + "\n" for line in lines))
    return path


# ---------------------------------------------------------------------------
# Confidence bigWig
# ---------------------------------------------------------------------------

_CHANNELS = {"intergenic": 0, "utr": 1, "cds": 2, "intron": 3}


def write_confidence_bigwig(
    h5_reader: HDF5ConfidenceReader,
    chrom_sizes: dict[str, int],
    path: str | Path,
    channel: str = "genic",
) -> Path:
    """Write Helixer per-base confidence as a bigWig; return Path.

    ``channel='genic'`` writes ``max(CDS, UTR)`` (the genic-confidence track used
    for ``helixer_locus_conf``); otherwise one of ``intergenic/utr/cds/intron``.
    Only seqids present in **both** ``chrom_sizes`` and the HDF5 are written.
    """
    import numpy as np

    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover - exercised via monkeypatch
        raise ImportError(_TRACKS_HINT) from exc

    if channel != "genic" and channel not in _CHANNELS:
        raise ValueError(
            f"channel must be 'genic' or one of {sorted(_CHANNELS)}, got {channel!r}"
        )

    path = Path(path)
    seqids = [s for s in h5_reader.seqids if s in chrom_sizes]
    bw = pyBigWig.open(str(path), "w")
    try:
        bw.addHeader([(s, int(chrom_sizes[s])) for s in sorted(seqids)])
        for s in sorted(seqids):
            length = int(chrom_sizes[s])
            preds = h5_reader.get_per_base_predictions(s, 0, length)
            if channel == "genic":
                values = np.maximum(
                    preds[:, _CHANNELS["cds"]], preds[:, _CHANNELS["utr"]]
                )
            else:
                values = preds[:, _CHANNELS[channel]]
            bw.addEntries(s, 0, values=[float(v) for v in values], span=1, step=1)
    finally:
        bw.close()
    return path


# ---------------------------------------------------------------------------
# Session / config files
# ---------------------------------------------------------------------------


def _track_items(
    track_paths: list[Any] | dict[str, Any],
) -> list[tuple[str, str]]:
    """Normalise ``track_paths`` (list or ``{name: path}``) → ``[(name, path)]``."""
    if isinstance(track_paths, dict):
        return [(name, str(p)) for name, p in track_paths.items()]
    return [(Path(p).name, str(p)) for p in track_paths]


_TRACK_FORMAT = {
    ".bb": "bigBed",
    ".bigbed": "bigBed",
    ".bw": "bigWig",
    ".bigwig": "bigWig",
    ".bed": "bed",
    ".bed12": "bed",
    ".bam": "bam",
}


def write_igv_session(
    track_paths: list[Any] | dict[str, Any],
    genome_fa: str | Path,
    path: str | Path,
) -> Path:
    """Write a minimal IGV session XML referencing the genome + tracks; return Path."""
    import xml.etree.ElementTree as ET

    session = ET.Element("Session", {"genome": str(genome_fa), "version": "8"})
    resources = ET.SubElement(session, "Resources")
    for _, p in _track_items(track_paths):
        ET.SubElement(resources, "Resource", {"path": p})
    for name, p in _track_items(track_paths):
        ET.SubElement(session, "Track", {"id": p, "name": name})

    path = Path(path)
    with atomic_write(path) as fh:
        ET.ElementTree(session).write(fh, encoding="unicode", xml_declaration=True)
    return path


def write_jbrowse_config(
    track_paths: list[Any] | dict[str, Any],
    path: str | Path,
    assembly_name: str = "HelixForge",
    genome_fa: str | Path | None = None,
) -> Path:
    """Write a JBrowse2 ``config.json`` (assembly + tracks); return Path.

    Each track gets a type inferred from its suffix (bigWig→QuantitativeTrack,
    everything else→FeatureTrack) and a ``uri`` adapter pointing at the emitted
    file. Standards-valid enough for JBrowse2 to load the emitted tracks.
    """
    tracks = []
    for name, p in _track_items(track_paths):
        fmt = _TRACK_FORMAT.get(Path(p).suffix.lower(), "bed")
        if fmt == "bigWig":
            track_type, adapter_type = "QuantitativeTrack", "BigWigAdapter"
            adapter: dict[str, Any] = {
                "type": adapter_type,
                "bigWigLocation": {"uri": p},
            }
        elif fmt == "bigBed":
            track_type, adapter_type = "FeatureTrack", "BigBedAdapter"
            adapter = {"type": adapter_type, "bigBedLocation": {"uri": p}}
        else:
            track_type, adapter_type = "FeatureTrack", "BedAdapter"
            adapter = {"type": adapter_type, "bedLocation": {"uri": p}}
        tracks.append(
            {
                "type": track_type,
                "trackId": name,
                "name": name,
                "assemblyNames": [assembly_name],
                "adapter": adapter,
            }
        )

    config = {
        "assemblies": [
            {
                "name": assembly_name,
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": f"{assembly_name}-ref",
                    "adapter": {
                        "type": "IndexedFastaAdapter",
                        "fastaLocation": {"uri": str(genome_fa) if genome_fa else ""},
                    },
                },
            }
        ],
        "tracks": tracks,
    }
    path = Path(path)
    with atomic_write(path) as fh:
        fh.write(json.dumps(config, indent=2))
    return path
