"""GFF3/GTF file handling.

This module provides readers and writers for GFF3 and GTF format files,
which are the standard formats for genome annotation.

Features:
    - Parse GFF3/GTF files into GeneModel objects
    - Write GeneModel objects to GFF3/GTF format
    - Validate GFF3 format compliance
    - Handle parent-child relationships
    - Streaming iteration for large files
    - Region-based queries

Example:
    >>> from helixforge.io.gff import GFF3Parser, GFF3Writer
    >>> parser = GFF3Parser("annotations.gff3")
    >>> for gene in parser.iter_genes():
    ...     print(gene.gene_id, len(gene.transcripts))
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator, Literal

import attrs

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

# GFF3 column indices
COL_SEQID = 0
COL_SOURCE = 1
COL_TYPE = 2
COL_START = 3
COL_END = 4
COL_SCORE = 5
COL_STRAND = 6
COL_PHASE = 7
COL_ATTRIBUTES = 8

# Standard feature types
FEATURE_GENE = "gene"
FEATURE_MRNA = "mRNA"
FEATURE_TRANSCRIPT = "transcript"
FEATURE_EXON = "exon"
FEATURE_CDS = "CDS"
FEATURE_UTR5 = "five_prime_UTR"
FEATURE_UTR3 = "three_prime_UTR"

# Helixer-specific feature types
FEATURE_TYPES_GENE = {"gene"}
FEATURE_TYPES_TRANSCRIPT = {"mRNA", "transcript", "ncRNA", "lnc_RNA"}
FEATURE_TYPES_EXON = {"exon"}
FEATURE_TYPES_CDS = {"CDS"}

Strand = Literal["+", "-"]


# =============================================================================
# Data Models
# =============================================================================


@attrs.define(slots=True)
class TranscriptModel:
    """Represents a transcript/mRNA with its features.

    Attributes:
        transcript_id: Unique transcript identifier.
        parent_gene: Parent gene ID.
        seqid: Scaffold/chromosome name.
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
        strand: Strand (+ or -).
        source: Annotation source.
        exons: List of (start, end) tuples for exons.
        cds: List of (start, end, phase) tuples for CDS.
        attributes: Additional attributes from GFF3.
    """

    transcript_id: str
    parent_gene: str
    seqid: str
    start: int
    end: int
    strand: str
    source: str = "helixer"
    exons: list[tuple[int, int]] = attrs.Factory(list)
    cds: list[tuple[int, int, int]] = attrs.Factory(list)  # (start, end, phase)
    attributes: dict[str, str] = attrs.Factory(dict)

    @property
    def n_exons(self) -> int:
        """Number of exons."""
        return len(self.exons)

    @property
    def n_cds(self) -> int:
        """Number of CDS segments."""
        return len(self.cds)

    @property
    def cds_length(self) -> int:
        """Total CDS length in nucleotides."""
        return sum(end - start for start, end, _ in self.cds)

    @property
    def transcript_length(self) -> int:
        """Total exon length (spliced transcript length)."""
        return sum(end - start for start, end in self.exons)

    @property
    def introns(self) -> list[tuple[int, int]]:
        """Get list of intron coordinates."""
        if len(self.exons) < 2:
            return []

        sorted_exons = sorted(self.exons)
        introns = []
        for i in range(len(sorted_exons) - 1):
            intron_start = sorted_exons[i][1]
            intron_end = sorted_exons[i + 1][0]
            if intron_start < intron_end:
                introns.append((intron_start, intron_end))
        return introns


@attrs.define(slots=True)
class GeneModel:
    """Represents a gene with its child features.

    Attributes:
        gene_id: Unique gene identifier.
        seqid: Scaffold/chromosome name.
        start: Start position (0-based, inclusive).
        end: End position (0-based, exclusive).
        strand: Strand (+ or -).
        source: Annotation source.
        transcripts: List of transcript models.
        attributes: Additional attributes from GFF3.
        confidence: Confidence score (added by HelixForge).
        qc_flags: QC flags (added by HelixForge).
    """

    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    source: str = "helixer"
    transcripts: list[TranscriptModel] = attrs.Factory(list)
    attributes: dict[str, str] = attrs.Factory(dict)
    confidence: float | None = None
    qc_flags: list[str] = attrs.Factory(list)

    @property
    def n_transcripts(self) -> int:
        """Number of transcripts."""
        return len(self.transcripts)

    @property
    def gene_length(self) -> int:
        """Gene length (end - start)."""
        return self.end - self.start

    @property
    def primary_transcript(self) -> TranscriptModel | None:
        """Get the primary/longest transcript."""
        if not self.transcripts:
            return None
        return max(self.transcripts, key=lambda t: t.cds_length or t.transcript_length)


# =============================================================================
# Attribute Parsing
# =============================================================================


def parse_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF3 attribute string into dictionary.

    Args:
        attr_string: Semicolon-separated key=value pairs.

    Returns:
        Dictionary of attribute key-value pairs.
    """
    attributes = {}
    if not attr_string or attr_string == ".":
        return attributes

    for item in attr_string.split(";"):
        item = item.strip()
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            # URL decode
            value = value.replace("%3B", ";").replace("%3D", "=").replace("%26", "&")
            value = value.replace("%2C", ",")
            attributes[key] = value

    return attributes


def format_attributes(attributes: dict[str, str]) -> str:
    """Format attribute dictionary as GFF3 string.

    Args:
        attributes: Dictionary of attributes.

    Returns:
        Semicolon-separated key=value string.
    """
    if not attributes:
        return "."

    parts = []
    for key, value in attributes.items():
        # URL encode special characters
        value = str(value).replace(";", "%3B").replace("=", "%3D")
        value = value.replace("&", "%26").replace(",", "%2C")
        parts.append(f"{key}={value}")

    return ";".join(parts)


# =============================================================================
# GFF3 Parser
# =============================================================================


class GFF3Parser:
    """Parse Helixer GFF3 output into structured gene models.

    Handles:
    - Parent-child relationships
    - Multiple transcripts per gene
    - Attribute parsing
    - Streaming iteration for large files

    Attributes:
        path: Path to the GFF3 file.
        gene_count: Total number of genes (after parsing).

    Example:
        >>> parser = GFF3Parser("annotations.gff3")
        >>> for gene in parser.iter_genes():
        ...     print(gene.gene_id, gene.n_transcripts)
    """

    def __init__(self, gff_path: Path | str) -> None:
        """Initialize the parser.

        Args:
            gff_path: Path to GFF3 file.

        Raises:
            FileNotFoundError: If file doesn't exist.
        """
        self.path = Path(gff_path)
        if not self.path.exists():
            raise FileNotFoundError(f"GFF3 file not found: {self.path}")

        self._genes: dict[str, GeneModel] | None = None
        self._index: dict[str, list[str]] | None = None  # seqid -> gene_ids

    def _parse_line(self, line: str) -> dict[str, Any] | None:
        """Parse a single GFF3 line.

        Args:
            line: Raw GFF3 line.

        Returns:
            Parsed feature dictionary or None for comments/empty.
        """
        line = line.strip()
        if not line or line.startswith("#"):
            return None

        parts = line.split("\t")
        if len(parts) < 9:
            logger.warning(f"Malformed GFF3 line (expected 9 columns): {line[:50]}...")
            return None

        try:
            # Parse coordinates (GFF3 is 1-based, convert to 0-based)
            start = int(parts[COL_START]) - 1  # Convert to 0-based
            end = int(parts[COL_END])  # Keep as exclusive end

            feature = {
                "seqid": parts[COL_SEQID],
                "source": parts[COL_SOURCE],
                "type": parts[COL_TYPE],
                "start": start,
                "end": end,
                "score": None if parts[COL_SCORE] == "." else float(parts[COL_SCORE]),
                "strand": parts[COL_STRAND] if parts[COL_STRAND] in ("+", "-") else "+",
                "phase": None if parts[COL_PHASE] == "." else int(parts[COL_PHASE]),
                "attributes": parse_attributes(parts[COL_ATTRIBUTES]),
            }
            return feature

        except (ValueError, IndexError) as e:
            logger.warning(f"Error parsing GFF3 line: {e}")
            return None

    def _build_genes(self) -> dict[str, GeneModel]:
        """Build gene models from GFF3 file."""
        genes: dict[str, GeneModel] = {}
        transcripts: dict[str, TranscriptModel] = {}

        # First pass: collect all features
        gene_features: dict[str, dict] = {}
        transcript_features: dict[str, dict] = {}
        exon_features: list[dict] = []
        cds_features: list[dict] = []

        with open(self.path) as f:
            for line in f:
                feature = self._parse_line(line)
                if feature is None:
                    continue

                ftype = feature["type"]
                attrs = feature["attributes"]

                if ftype in FEATURE_TYPES_GENE:
                    gene_id = attrs.get("ID", attrs.get("gene_id", f"gene_{len(gene_features)}"))
                    gene_features[gene_id] = feature

                elif ftype in FEATURE_TYPES_TRANSCRIPT:
                    tx_id = attrs.get("ID", attrs.get("transcript_id", f"tx_{len(transcript_features)}"))
                    transcript_features[tx_id] = feature

                elif ftype in FEATURE_TYPES_EXON:
                    exon_features.append(feature)

                elif ftype in FEATURE_TYPES_CDS:
                    cds_features.append(feature)

        # Build gene models
        for gene_id, gf in gene_features.items():
            gene = GeneModel(
                gene_id=gene_id,
                seqid=gf["seqid"],
                start=gf["start"],
                end=gf["end"],
                strand=gf["strand"],
                source=gf["source"],
                attributes=gf["attributes"],
            )
            genes[gene_id] = gene

        # Build transcript models and associate with genes
        for tx_id, tf in transcript_features.items():
            parent = tf["attributes"].get("Parent", "")
            # Handle multiple parents
            parent_ids = parent.split(",") if parent else []

            transcript = TranscriptModel(
                transcript_id=tx_id,
                parent_gene=parent_ids[0] if parent_ids else "",
                seqid=tf["seqid"],
                start=tf["start"],
                end=tf["end"],
                strand=tf["strand"],
                source=tf["source"],
                attributes=tf["attributes"],
            )
            transcripts[tx_id] = transcript

            # Add to parent gene
            for parent_id in parent_ids:
                if parent_id in genes:
                    genes[parent_id].transcripts.append(transcript)

        # Add exons to transcripts
        for ef in exon_features:
            parent = ef["attributes"].get("Parent", "")
            parent_ids = parent.split(",") if parent else []

            for parent_id in parent_ids:
                if parent_id in transcripts:
                    transcripts[parent_id].exons.append((ef["start"], ef["end"]))

        # Add CDS to transcripts
        for cf in cds_features:
            parent = cf["attributes"].get("Parent", "")
            parent_ids = parent.split(",") if parent else []
            phase = cf["phase"] if cf["phase"] is not None else 0

            for parent_id in parent_ids:
                if parent_id in transcripts:
                    transcripts[parent_id].cds.append((cf["start"], cf["end"], phase))

        # Sort exons and CDS within each transcript
        for tx in transcripts.values():
            tx.exons.sort()
            tx.cds.sort()

        # Build index
        self._index = defaultdict(list)
        for gene_id, gene in genes.items():
            self._index[gene.seqid].append(gene_id)

        logger.info(f"Parsed {len(genes)} genes, {len(transcripts)} transcripts")
        return genes

    def _ensure_parsed(self) -> None:
        """Ensure the GFF3 file has been parsed."""
        if self._genes is None:
            self._genes = self._build_genes()

    def iter_genes(self) -> Iterator[GeneModel]:
        """Iterate over genes with their child features.

        Yields:
            GeneModel objects.
        """
        self._ensure_parsed()
        assert self._genes is not None
        yield from self._genes.values()

    def get_gene(self, gene_id: str) -> GeneModel | None:
        """Retrieve specific gene by ID.

        Args:
            gene_id: Gene identifier.

        Returns:
            GeneModel or None if not found.
        """
        self._ensure_parsed()
        assert self._genes is not None
        return self._genes.get(gene_id)

    def get_genes_in_region(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> list[GeneModel]:
        """Get all genes overlapping a region.

        Args:
            seqid: Scaffold name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            List of overlapping GeneModel objects.
        """
        self._ensure_parsed()
        assert self._genes is not None
        assert self._index is not None

        result = []
        for gene_id in self._index.get(seqid, []):
            gene = self._genes[gene_id]
            # Check overlap
            if gene.start < end and gene.end > start:
                result.append(gene)

        return result

    @property
    def gene_count(self) -> int:
        """Total number of genes."""
        self._ensure_parsed()
        assert self._genes is not None
        return len(self._genes)

    @property
    def scaffold_gene_counts(self) -> dict[str, int]:
        """Genes per scaffold."""
        self._ensure_parsed()
        assert self._index is not None
        return {seqid: len(genes) for seqid, genes in self._index.items()}

    @property
    def gene_ids(self) -> list[str]:
        """List of all gene IDs."""
        self._ensure_parsed()
        assert self._genes is not None
        return list(self._genes.keys())


# =============================================================================
# GFF3 Writer
# =============================================================================


class GFF3Writer:
    """Write refined gene models to GFF3 format.

    Features:
    - Preserve original attributes
    - Add HelixForge-specific attributes (confidence, flags)
    - Provenance tracking

    Example:
        >>> writer = GFF3Writer("output.gff3")
        >>> writer.write_header(genome_path="genome.fa")
        >>> for gene in genes:
        ...     writer.write_gene(gene)
        >>> writer.close()
    """

    def __init__(
        self,
        output_path: Path | str,
        source: str = "HelixForge",
    ) -> None:
        """Initialize the writer.

        Args:
            output_path: Output file path.
            source: Source field value for GFF3.
        """
        self.path = Path(output_path)
        self.source = source
        self._file = open(self.path, "w")
        self._header_written = False

    def __enter__(self) -> GFF3Writer:
        """Context manager entry."""
        return self

    def __exit__(self, *args: Any) -> None:
        """Context manager exit."""
        self.close()

    def close(self) -> None:
        """Close the output file."""
        if self._file:
            self._file.close()
            self._file = None

    def write_header(
        self,
        genome_path: Path | str | None = None,
        helixer_gff: Path | str | None = None,
        helixforge_version: str = "0.1.0",
    ) -> None:
        """Write GFF3 header with provenance.

        Args:
            genome_path: Path to reference genome.
            helixer_gff: Path to original Helixer GFF.
            helixforge_version: HelixForge version.
        """
        self._file.write("##gff-version 3\n")
        self._file.write(f"#!processor HelixForge v{helixforge_version}\n")

        if genome_path:
            self._file.write(f"#!genome-build-file {genome_path}\n")
        if helixer_gff:
            self._file.write(f"#!original-file {helixer_gff}\n")

        self._header_written = True

    def _format_line(
        self,
        seqid: str,
        feature_type: str,
        start: int,
        end: int,
        strand: str,
        attributes: dict[str, str],
        score: float | None = None,
        phase: int | None = None,
    ) -> str:
        """Format a GFF3 line.

        Args:
            seqid: Scaffold name.
            feature_type: Feature type.
            start: Start (0-based).
            end: End (0-based, exclusive).
            strand: Strand.
            attributes: Attribute dictionary.
            score: Feature score.
            phase: CDS phase.

        Returns:
            Formatted GFF3 line.
        """
        # Convert to 1-based coordinates for GFF3
        gff_start = start + 1
        gff_end = end

        score_str = "." if score is None else f"{score:.4f}"
        phase_str = "." if phase is None else str(phase)
        attr_str = format_attributes(attributes)

        return f"{seqid}\t{self.source}\t{feature_type}\t{gff_start}\t{gff_end}\t{score_str}\t{strand}\t{phase_str}\t{attr_str}\n"

    def write_gene(self, gene: GeneModel) -> None:
        """Write a gene and its child features.

        Args:
            gene: GeneModel to write.
        """
        if not self._header_written:
            self.write_header()

        # Gene attributes
        gene_attrs = {"ID": gene.gene_id}
        gene_attrs.update(gene.attributes)

        # Add HelixForge attributes
        if gene.confidence is not None:
            gene_attrs["confidence"] = f"{gene.confidence:.4f}"
        if gene.qc_flags:
            gene_attrs["qc_flags"] = ",".join(gene.qc_flags)

        # Write gene line
        self._file.write(
            self._format_line(
                gene.seqid,
                FEATURE_GENE,
                gene.start,
                gene.end,
                gene.strand,
                gene_attrs,
            )
        )

        # Write transcripts and their features
        for transcript in gene.transcripts:
            self._write_transcript(transcript)

    def _write_transcript(self, transcript: TranscriptModel) -> None:
        """Write a transcript and its features."""
        # Transcript attributes
        tx_attrs = {
            "ID": transcript.transcript_id,
            "Parent": transcript.parent_gene,
        }
        tx_attrs.update(transcript.attributes)

        # Write mRNA line
        self._file.write(
            self._format_line(
                transcript.seqid,
                FEATURE_MRNA,
                transcript.start,
                transcript.end,
                transcript.strand,
                tx_attrs,
            )
        )

        # Write exons
        for i, (start, end) in enumerate(transcript.exons, 1):
            exon_attrs = {
                "ID": f"{transcript.transcript_id}.exon{i}",
                "Parent": transcript.transcript_id,
            }
            self._file.write(
                self._format_line(
                    transcript.seqid,
                    FEATURE_EXON,
                    start,
                    end,
                    transcript.strand,
                    exon_attrs,
                )
            )

        # Write CDS
        for i, (start, end, phase) in enumerate(transcript.cds, 1):
            cds_attrs = {
                "ID": f"{transcript.transcript_id}.CDS{i}",
                "Parent": transcript.transcript_id,
            }
            self._file.write(
                self._format_line(
                    transcript.seqid,
                    FEATURE_CDS,
                    start,
                    end,
                    transcript.strand,
                    cds_attrs,
                    phase=phase,
                )
            )

    def write_genes(self, genes: list[GeneModel]) -> None:
        """Write multiple genes.

        Args:
            genes: List of GeneModel objects.
        """
        for gene in genes:
            self.write_gene(gene)


# =============================================================================
# Convenience Functions
# =============================================================================


def read_gff(path: Path | str) -> list[GeneModel]:
    """Read gene models from a GFF3 file.

    Args:
        path: Path to the GFF3 file.

    Returns:
        List of GeneModel objects.
    """
    parser = GFF3Parser(path)
    return list(parser.iter_genes())


def iter_gff(path: Path | str) -> Iterator[GeneModel]:
    """Iterate over gene models from a GFF3 file.

    Args:
        path: Path to the GFF3 file.

    Yields:
        GeneModel objects.
    """
    parser = GFF3Parser(path)
    yield from parser.iter_genes()


def write_gff(
    genes: list[GeneModel],
    path: Path | str,
    source: str = "HelixForge",
) -> None:
    """Write gene models to a GFF3 file.

    Args:
        genes: List of GeneModel objects.
        path: Output file path.
        source: Source field value.
    """
    with GFF3Writer(path, source=source) as writer:
        writer.write_header()
        writer.write_genes(genes)


def format_gff_line(
    seqid: str,
    source: str,
    feature_type: str,
    start: int,
    end: int,
    score: float | None = None,
    strand: str = ".",
    phase: int | None = None,
    attributes: dict[str, str] | None = None,
) -> str:
    """Format a single GFF3 line.

    Args:
        seqid: Sequence identifier.
        source: Source of the annotation.
        feature_type: Type of feature.
        start: Start position (0-based).
        end: End position (0-based, exclusive).
        score: Feature score.
        strand: Strand.
        phase: CDS phase.
        attributes: Feature attributes.

    Returns:
        Formatted GFF3 line.
    """
    # Convert to 1-based for GFF3
    gff_start = start + 1
    gff_end = end

    score_str = "." if score is None else f"{score:.4f}"
    phase_str = "." if phase is None else str(phase)
    attr_str = format_attributes(attributes or {})

    return f"{seqid}\t{source}\t{feature_type}\t{gff_start}\t{gff_end}\t{score_str}\t{strand}\t{phase_str}\t{attr_str}"
