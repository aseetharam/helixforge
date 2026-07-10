"""Structured-ncRNA hook: tRNAscan-SE / Infernal(Rfam)."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

from helixforge.constants import NCRNA_ID_BASE
from helixforge.prep._subprocess import output_is_fresh, run_tool
from helixforge.qc.flags import NOVEL_LOCUS
from helixforge.utils.logging import get_logger

if TYPE_CHECKING:
    from helixforge.reconcile.models import ReconciledGene

_log = get_logger(__name__)

# tRNAscan-SE isotype/note → biotype. Everything tRNAscan reports is a tRNA; the
# pseudogene note maps to a pseudogene biotype (a tRNA-derived pseudogene).
_TRNASCAN_PSEUDO_NOTE = "pseudo"


def run_trnascan_se(
    genome_fasta: str | Path,
    out_path: str | Path,
    *,
    trnascan_bin: str = "tRNAscan-SE",
    mode: str = "-E",
    threads: int = 1,
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
) -> Path:
    """``tRNAscan-SE <mode> --thread <n> -o <out> <genome>`` → ``out_path``.

    ``mode`` selects the search model (``-E`` eukaryotic, ``-B`` bacterial,
    ``-A`` archaeal, ``-O`` organellar). Skips when ``out_path`` is already fresh
    unless ``force``. Returns the output path.
    """
    out_path = Path(out_path)
    if not force and output_is_fresh(out_path, [genome_fasta]):
        _log.info("skip tRNAscan-SE: %s up to date", out_path.name)
        return out_path
    argv: list[str | Path | int] = [
        trnascan_bin,
        mode,
        "--thread",
        threads,
        "-o",
        out_path,
        genome_fasta,
    ]
    if extra_args:
        argv += list(extra_args)
    run_tool(argv, log=_log)
    return out_path


def run_infernal_cmscan(
    genome_fasta: str | Path,
    rfam_cm: str | Path,
    out_tblout: str | Path,
    *,
    cmscan_bin: str = "cmscan",
    clanin: str | Path | None = None,
    threads: int = 1,
    extra_args: Iterable[str | Path] | None = None,
    force: bool = False,
) -> Path:
    """``cmscan --cut_ga --rfam [--clanin <f>] --cpu <n> --tblout <out> <cm> <genome>``.

    Scans the genome against an Rfam covariance-model database (``rfam_cm``) and
    writes the tabular hit table to ``out_tblout`` (returned). ``--cut_ga`` /
    ``--rfam`` are the Rfam-recommended thresholds; ``clanin`` enables clan-based
    overlap resolution. Skips when fresh unless ``force``.
    """
    out_tblout = Path(out_tblout)
    if not force and output_is_fresh(out_tblout, [genome_fasta, rfam_cm]):
        _log.info("skip cmscan: %s up to date", out_tblout.name)
        return out_tblout
    argv: list[str | Path | int] = [cmscan_bin, "--cut_ga", "--rfam"]
    if clanin is not None:
        argv += ["--clanin", clanin]
    argv += [
        "--cpu",
        threads,
        "--tblout",
        out_tblout,
        "--fmt",
        "2",
        rfam_cm,
        genome_fasta,
    ]
    if extra_args:
        argv += list(extra_args)
    run_tool(argv, log=_log)
    return out_tblout


# ---------------------------------------------------------------------------
# Parsers → (seqid, start, end, strand, biotype) tuples (0-based half-open)
# ---------------------------------------------------------------------------


def parse_trnascan(path: str | Path) -> list[tuple[str, int, int, str, str]]:
    """Parse a tRNAscan-SE tabular ``-o`` output into ncRNA loci.

    tRNAscan columns: ``Name  tRNA#  Begin  End  Type  Anticodon  ...  Note``.
    Coordinates are 1-based inclusive and **strand-encoded** (Begin > End on the
    minus strand); they are normalised to internal 0-based half-open ``(low, high)``
    with an explicit strand. A ``pseudo`` note maps to ``pseudogene``; everything
    else is ``tRNA``.
    """
    out: list[tuple[str, int, int, str, str]] = []
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t") if "\t" in line else line.split()
        if len(fields) < 4:
            continue
        seqid = fields[0].strip()
        # Skip the two header rows + the dashed separator row tRNAscan emits.
        if not seqid or seqid.lower() in ("name", "sequence") or seqid.startswith("-"):
            continue
        try:
            begin = int(fields[2])
            end = int(fields[3])
        except ValueError:
            continue
        strand = "+" if begin <= end else "-"
        low, high = (begin - 1, end) if begin <= end else (end - 1, begin)
        note = fields[-1].lower() if fields else ""
        biotype = "pseudogene" if _TRNASCAN_PSEUDO_NOTE in note else "tRNA"
        out.append((seqid, low, high, strand, biotype))
    return out


def build_ncrna_genes(
    records: Iterable[tuple[str, int, int, str, str]],
    *,
    id_base: int = NCRNA_ID_BASE,
) -> list["ReconciledGene"]:
    """Build single-exon ``ReconciledGene`` ncRNA loci from parsed records.

    Each record is ``(seqid, start, end, strand, biotype)`` (internal 0-based
    half-open). Genes are origin ``novel`` (no Helixer anchor), Tier 4 (no
    CDS/expression evidence, the biotype carries the meaning), flagged
    ``NOVEL_LOCUS``, and ids ``HFG_<id_base+i>``. Sorted by ``(seqid, start)``.
    """
    from helixforge.reconcile.models import (
        Exon,
        LocusClassification,
        ReconciledGene,
        TranscriptCandidate,
    )

    ordered = sorted(records, key=lambda r: (r[0], r[1], r[2]))
    genes: list[ReconciledGene] = []
    for i, (seqid, start, end, strand, biotype) in enumerate(ordered):
        gene_id = f"HFG_{id_base + i:05d}"
        tid = f"{gene_id}.1"
        tx = TranscriptCandidate(
            transcript_id=tid,
            locus_id=gene_id,
            source="ncrna",
            seqid=seqid,
            start=start,
            end=end,
            strand=strand,
            exons=[Exon(start, end)],
            is_primary=True,
            biotype=biotype,
        )
        genes.append(
            ReconciledGene(
                gene_id=gene_id,
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                tier=4,
                transcripts=[tx],
                primary_transcript_id=tid,
                classification=LocusClassification(
                    gene_id, "SILENT", evidence_source="none"
                ),
                origin="novel",
                flags=[NOVEL_LOCUS],
                biotype=biotype,
            )
        )
    return genes


def scan_structured_ncrna(
    genome_fasta: str | Path,
    work_dir: str | Path,
    *,
    enabled: bool = False,
    tool: str = "trnascan",
    trnascan_bin: str = "tRNAscan-SE",
    trnascan_mode: str = "-E",
    cmscan_bin: str = "cmscan",
    rfam_cm: str | Path | None = None,
    clanin: str | Path | None = None,
    threads: int = 1,
    id_base: int = NCRNA_ID_BASE,
    force: bool = False,
) -> list["ReconciledGene"]:
    """Run the structured-ncRNA hook; return the ncRNA loci (``[]`` when disabled).

    **Off by default** (``enabled=False`` → returns ``[]`` and runs no subprocess),
    so the default pipeline path is unaffected. When enabled, runs ``tool`` ∈
    ``{trnascan, infernal}`` and parses its output into ncRNA genes. Infernal
    requires ``rfam_cm``. The Infernal tabular parser is intentionally minimal
    (it relies on the external scan output, which is mocked in tests).
    """
    if not enabled:
        return []
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    if tool == "trnascan":
        out = run_trnascan_se(
            genome_fasta,
            work_dir / "trnascan.out",
            trnascan_bin=trnascan_bin,
            mode=trnascan_mode,
            threads=threads,
            force=force,
        )
        records = parse_trnascan(out)
    elif tool == "infernal":
        if rfam_cm is None:
            raise ValueError("infernal scan requires rfam_cm (the Rfam CM database)")
        out = run_infernal_cmscan(
            genome_fasta,
            rfam_cm,
            work_dir / "cmscan.tblout",
            cmscan_bin=cmscan_bin,
            clanin=clanin,
            threads=threads,
            force=force,
        )
        records = parse_infernal_tblout(out)
    else:
        raise ValueError(f"tool must be 'trnascan' or 'infernal', got {tool!r}")
    genes = build_ncrna_genes(records, id_base=id_base)
    _log.info("structured-ncRNA hook (%s): %d loci", tool, len(genes))
    return genes


# Rfam/Infernal biotype mapping by hit type-name prefix (cmscan --fmt 2 'tblout').
_INFERNAL_BIOTYPE_PREFIX = (
    ("tRNA", "tRNA"),
    ("5S_rRNA", "rRNA"),
    ("5_8S_rRNA", "rRNA"),
    ("SSU_rRNA", "rRNA"),
    ("LSU_rRNA", "rRNA"),
    ("rRNA", "rRNA"),
    ("SNOR", "snoRNA"),
    ("snoR", "snoRNA"),
    ("U1", "snRNA"),
    ("U2", "snRNA"),
    ("U4", "snRNA"),
    ("U5", "snRNA"),
    ("U6", "snRNA"),
    ("MIR", "miRNA"),
    ("mir-", "miRNA"),
)


def _infernal_biotype(target_name: str) -> str:
    name = target_name or ""
    for prefix, biotype in _INFERNAL_BIOTYPE_PREFIX:
        if name.upper().startswith(prefix.upper()):
            return biotype
    return "ncRNA_undetermined"


def parse_infernal_tblout(path: str | Path) -> list[tuple[str, int, int, str, str]]:
    """Parse an Infernal ``cmscan --fmt 2 --tblout`` table into ncRNA loci.

    Reads the model name (col 2), the genome seqid (col 4), the strand (col 12),
    and the sequence-from/to (cols 10/11; strand-encoded). Comment lines (``#``)
    are skipped. The model name maps to a biotype via :func:`_infernal_biotype`.
    """
    out: list[tuple[str, int, int, str, str]] = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 12:
            continue
        model_name, seqid = f[1], f[3]
        try:
            seq_from, seq_to = int(f[9]), int(f[10])
        except ValueError:
            continue
        strand = f[11] if f[11] in ("+", "-") else ("+" if seq_from <= seq_to else "-")
        low, high = (
            (seq_from - 1, seq_to) if seq_from <= seq_to else (seq_to - 1, seq_from)
        )
        out.append((seqid, low, high, strand, _infernal_biotype(model_name)))
    return out
