"""GFF3 parsing and writing."""

from __future__ import annotations

import os
from typing import Any, IO

import gffutils

from helixforge.reconcile.models import CDSSegment, Exon, HelixerLocus
from helixforge.utils.atomic import atomic_write
from helixforge.utils.regions import gff3_to_internal, internal_to_gff3

_TRANSCRIPT_TYPES = ("mRNA", "transcript")


def _phase_from_frame(frame: str) -> int:
    """gffutils frame string → phase int ∈ {0,1,2} (default 0 if absent)."""
    return int(frame) if frame in ("0", "1", "2") else 0


def build_gffutils_db(
    source_path: str | os.PathLike[str],
    dbfn: str | None = None,
    keep_db: bool = False,
) -> Any:
    """Build (or reuse) a gffutils DB for ``source_path``.

    The single home for the ``create_db(..., dbfn=":memory:", force=True)`` call
    that ``GFF3Parser``, ``MiniprotParser`` and ``parse_loci_gff3`` each made on
    every construction — a full re-parse + index build with no caching.

    * ``dbfn=None`` (default): the historical in-memory DB, rebuilt every time —
      byte-for-byte unchanged from before, so existing tests/timings hold.
    * ``dbfn`` given + the file already exists + ``keep_db=True``: open and
      **reuse** the on-disk DB via ``FeatureDB`` — no re-parse of the GFF3.
    * ``dbfn`` given + the file exists + ``keep_db=False``: rebuild it
      (``force=True``), treating the existing file as stale.
    * ``dbfn`` given + the file is absent: build it on disk (``force=False``).

    The pipeline may pass a persistent ``.gffutils.db`` path with ``keep_db=True``
    in a later phase (Phase 22/24) so the same GFF3 is parsed once and memory-
    mapped on reuse.
    """
    if dbfn is None:
        return gffutils.create_db(
            str(source_path),
            dbfn=":memory:",
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    if keep_db and os.path.exists(dbfn):
        return gffutils.FeatureDB(str(dbfn), keep_order=True)
    return gffutils.create_db(
        str(source_path),
        dbfn=str(dbfn),
        force=not keep_db,  # absent file → harmless; stale file → rebuild
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


class GFF3Parser:
    """Parse a GFF3 file into HelixForge models (gffutils-backed).

    ``dbfn``/``keep_db``: pass an on-disk path to build the gffutils
    DB once and reuse it instead of rebuilding ``:memory:`` each construction; the
    default (``None``) is the historical in-memory behavior. See
    :func:`build_gffutils_db`.
    """

    def __init__(
        self,
        gff_path: str | os.PathLike[str],
        dbfn: str | None = None,
        keep_db: bool = False,
    ) -> None:
        if not os.path.exists(gff_path):
            raise FileNotFoundError(f"GFF3 not found: {gff_path}")
        self.gff_path = str(gff_path)
        self._db: Any = build_gffutils_db(self.gff_path, dbfn=dbfn, keep_db=keep_db)

    def _first_transcript(self, gene: Any) -> list[Any]:
        for ttype in _TRANSCRIPT_TYPES:
            kids = list(self._db.children(gene, featuretype=ttype, order_by="start"))
            if kids:
                return kids
        return []

    def parse_helixer_genes(self) -> list[HelixerLocus]:
        """Parse a Helixer gene→mRNA→exon,CDS hierarchy into ``HelixerLocus``.

        One mRNA per gene (Helixer's design); if several are present the first
        is used. CDS is optional. Returns a list sorted by (seqid, start).
        """
        loci: list[HelixerLocus] = []
        for gene in self._db.features_of_type("gene", order_by=("seqid", "start")):
            transcripts = self._first_transcript(gene)
            if not transcripts:
                continue
            mrna = transcripts[0]

            exons = [
                Exon(*gff3_to_internal(e.start, e.end))
                for e in self._db.children(mrna, featuretype="exon", order_by="start")
            ]
            cds_feats = list(
                self._db.children(mrna, featuretype="CDS", order_by="start")
            )
            cds: list[CDSSegment] | None = (
                [
                    CDSSegment(
                        *gff3_to_internal(c.start, c.end), _phase_from_frame(c.frame)
                    )
                    for c in cds_feats
                ]
                if cds_feats
                else None
            )

            g_start, g_end = gff3_to_internal(gene.start, gene.end)
            loci.append(
                HelixerLocus(
                    gene_id=gene.id,
                    seqid=gene.seqid,
                    start=g_start,
                    end=g_end,
                    strand=gene.strand,
                    exons=exons,
                    cds=cds,
                )
            )
        loci.sort(key=lambda g: (g.seqid, g.start))
        return loci

    def parse_genes_generic(self) -> list[dict[str, Any]]:
        """Parse an arbitrary GFF3 into gene dicts (0-based half-open).

        Each gene dict: ``{gene_id, seqid, start, end, strand, transcripts}``
        where each transcript is ``{transcript_id, exons: list[Exon],
        cds: list[CDSSegment] | None}``. Multiple transcripts per gene are
        preserved. Used later for reference annotations in benchmarking.
        """
        genes: list[dict[str, Any]] = []
        for gene in self._db.features_of_type("gene", order_by=("seqid", "start")):
            g_start, g_end = gff3_to_internal(gene.start, gene.end)
            transcripts: list[dict[str, Any]] = []
            for ttype in _TRANSCRIPT_TYPES:
                for tx in self._db.children(gene, featuretype=ttype, order_by="start"):
                    exons = [
                        Exon(*gff3_to_internal(e.start, e.end))
                        for e in self._db.children(
                            tx, featuretype="exon", order_by="start"
                        )
                    ]
                    cds_feats = list(
                        self._db.children(tx, featuretype="CDS", order_by="start")
                    )
                    cds: list[CDSSegment] | None = (
                        [
                            CDSSegment(
                                *gff3_to_internal(c.start, c.end),
                                _phase_from_frame(c.frame),
                            )
                            for c in cds_feats
                        ]
                        if cds_feats
                        else None
                    )
                    transcripts.append(
                        {"transcript_id": tx.id, "exons": exons, "cds": cds}
                    )
            genes.append(
                {
                    "gene_id": gene.id,
                    "seqid": gene.seqid,
                    "start": g_start,
                    "end": g_end,
                    "strand": gene.strand,
                    "transcripts": transcripts,
                }
            )
        return genes

    def get_features_in_region(
        self,
        seqid: str,
        start: int,
        end: int,
        featuretype: str = "gene",
    ) -> list[Any]:
        """Return gffutils Features of ``featuretype`` overlapping ``[start, end)``.

        The query is internal (0-based half-open); it is converted to GFF3
        coordinates for the gffutils region query. The caller converts the
        returned Features back to internal coordinates.
        """
        g_start, g_end = internal_to_gff3(start, end)
        return list(
            self._db.region(
                seqid=seqid, start=g_start, end=g_end, featuretype=featuretype
            )
        )


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

# Characters that must be escaped within GFF3 column-9 attribute values.
_ESCAPE: list[tuple[str, str]] = [
    ("%", "%25"),  # must be first
    (";", "%3B"),
    ("=", "%3D"),
    ("&", "%26"),
    (",", "%2C"),
    ("\t", "%09"),
    ("\n", "%0A"),
    ("\r", "%0D"),
]


def _encode(value: Any) -> str:
    s = str(value)
    for ch, rep in _ESCAPE:
        s = s.replace(ch, rep)
    return s


def _fmt_value(value: Any) -> str:
    """Scalars are encoded; lists/tuples become comma-joined encoded tokens."""
    if isinstance(value, (list, tuple)):
        return ",".join(_encode(v) for v in value)
    return _encode(value)


def _attrs(pairs: list[tuple[str, Any]]) -> str:
    """Build a column-9 string from (key, value) pairs, dropping empties."""
    parts = []
    for key, value in pairs:
        if value is None or value == "" or value == []:
            continue
        parts.append(f"{key}={_fmt_value(value)}")
    return ";".join(parts)


class GFF3Writer:
    """Hand-rolled GFF3 writer for ``ReconciledGene`` objects.

    Spec-completeness knobs — all
    optional; **omitting every one keeps the output byte-for-byte identical** to
    the historical ``##gff-version 3``-only header, so existing callers/tests are
    unchanged:

    * ``sequence_regions`` — ``{seqid: length}``; emits one ``##sequence-region
      <seqid> 1 <len>`` directive per contig (1-based inclusive per the spec).
    * ``provenance`` — an object exposing ``to_gff3_lines()`` (see
      :mod:`helixforge.provenance`); its ``#!`` lines are written as comments.
    * ``embed_fasta`` — a genome FASTA path; its sequence is appended after a
      ``##FASTA`` directive so the GFF3 is fully self-describing.
    """

    def __init__(
        self,
        output_path: str | os.PathLike[str],
        *,
        sequence_regions: dict[str, int] | None = None,
        provenance: Any | None = None,
        embed_fasta: str | os.PathLike[str] | None = None,
    ) -> None:
        self.output_path = str(output_path)
        self.sequence_regions = sequence_regions
        self.provenance = provenance
        self.embed_fasta = embed_fasta

    def write_header(self, fh: IO[str]) -> None:
        # ``##gff-version 3`` MUST be the first line (GFF3 spec); provenance
        # comments and ``##sequence-region`` directives follow.
        fh.write("##gff-version 3\n")
        if self.provenance is not None:
            for line in self.provenance.to_gff3_lines():
                fh.write(line + "\n")
        if self.sequence_regions:
            for seqid in sorted(self.sequence_regions):
                length = self.sequence_regions[seqid]
                # GFF3 ##sequence-region is 1-based inclusive: a length-N contig
                # spans 1..N.
                fh.write(f"##sequence-region {seqid} 1 {length}\n")

    def write_footer(self, fh: IO[str]) -> None:
        # Optionally embed the genome so the GFF3 is self-contained (##FASTA must
        # be the last directive; everything after it is sequence).
        if self.embed_fasta is not None:
            fh.write("##FASTA\n")
            with open(self.embed_fasta) as src:
                for line in src:
                    fh.write(line)
        return None

    def write_genes(
        self,
        genes: Any,
        source: str = "HelixForge",
        functional: dict[str, Any] | None = None,
    ) -> None:
        """Write the gene→mRNA→exon/CDS hierarchy.

        ``functional`` optionally maps ``transcript_id`` → a
        ``FunctionalRecord`` (``go_terms`` / ``dbxrefs``); when present, each mRNA
        gains GFF3-reserved ``Ontology_term=`` (GO) and ``Dbxref=`` (Pfam/InterPro/
        UniProt) attributes. The default ``None`` keeps the output byte-identical
        (count-neutral: never adds, drops, or reorders genes).
        """
        # Atomic: never leave a partial GFF3 on a crash/kill.
        with atomic_write(self.output_path) as fh:
            self.write_header(fh)
            for gene in genes:
                self._write_gene(fh, gene, source, functional)
                fh.write("###\n")
            self.write_footer(fh)

    # --- internals ---
    def _line(
        self,
        seqid: str,
        source: str,
        ftype: str,
        start: int,
        end: int,
        strand: str,
        phase: str,
        attrs_str: str,
    ) -> str:
        g_start, g_end = internal_to_gff3(start, end)
        return (
            f"{seqid}\t{source}\t{ftype}\t{g_start}\t{g_end}\t.\t"
            f"{strand}\t{phase}\t{attrs_str}\n"
        )

    def _write_gene(
        self,
        fh: IO[str],
        gene: Any,
        source: str,
        functional: dict[str, Any] | None = None,
    ) -> None:
        strand = gene.strand
        gene_attrs = _attrs(
            [
                ("ID", gene.gene_id),
                ("Name", gene.gene_id),
                ("gene_biotype", getattr(gene, "biotype", None)),
                ("tier", gene.tier),
                ("origin", gene.origin),
                ("status", gene.classification.status),
                ("merged_from", list(gene.merged_from) or None),
            ]
        )
        fh.write(
            self._line(
                gene.seqid,
                source,
                "gene",
                gene.start,
                gene.end,
                strand,
                ".",
                gene_attrs,
            )
        )

        flag_names: list[str] | None = [f.name for f in gene.flags] or None
        as_summary: list[str] | None = [
            f"{e.kind}@{internal_to_gff3(e.start, e.end)[0]}-{internal_to_gff3(e.start, e.end)[1]}"
            for e in gene.as_events
        ] or None

        for tx in gene.transcripts:
            # GFF3-reserved functional attributes from the optional
            # InterProScan/eggNOG hook. ``Ontology_term`` carries GO,
            # ``Dbxref`` carries Pfam/InterPro/UniProt. Absent → omitted (None).
            func = functional.get(tx.transcript_id) if functional else None
            go_terms = list(func.go_terms) if func else None
            dbxrefs = list(func.dbxrefs) if func else None
            mrna_attrs = _attrs(
                [
                    ("ID", tx.transcript_id),
                    ("Parent", gene.gene_id),
                    (
                        "transcript_biotype",
                        getattr(tx, "biotype", None) or getattr(gene, "biotype", None),
                    ),
                    ("tier", gene.tier),
                    ("origin", gene.origin),
                    ("structure_source", tx.source),
                    ("protein_id", tx.protein_id),
                    ("tpm", tx.tpm),
                    ("junction_support", tx.junction_support_fraction),
                    ("combined_score", tx.combined_score),
                    ("primary", "true" if tx.is_primary else None),
                    ("Ontology_term", go_terms or None),
                    ("Dbxref", dbxrefs or None),
                    ("flags", flag_names),
                    ("as_events", as_summary),
                ]
            )
            fh.write(
                self._line(
                    tx.seqid, source, "mRNA", tx.start, tx.end, strand, ".", mrna_attrs
                )
            )

            # exon/CDS feature IDs: 1-based genomic order (low→high) regardless of strand.
            for i, ex in enumerate(tx.exons, start=1):
                exon_attrs = _attrs(
                    [
                        ("ID", f"{tx.transcript_id}.exon{i}"),
                        ("Parent", tx.transcript_id),
                    ]
                )
                fh.write(
                    self._line(
                        tx.seqid,
                        source,
                        "exon",
                        ex.start,
                        ex.end,
                        strand,
                        ".",
                        exon_attrs,
                    )
                )
            if tx.cds:
                for i, c in enumerate(tx.cds, start=1):
                    cds_attrs = _attrs(
                        [
                            ("ID", f"{tx.transcript_id}.CDS{i}"),
                            ("Parent", tx.transcript_id),
                        ]
                    )
                    fh.write(
                        self._line(
                            tx.seqid,
                            source,
                            "CDS",
                            c.start,
                            c.end,
                            strand,
                            str(c.phase),
                            cds_attrs,
                        )
                    )
