"""Homology search using Diamond or MMseqs2.

This module provides wrappers for running protein homology searches
against reference databases like Swiss-Prot.

Features:
    - Run Diamond or MMseqs2 blastp searches
    - Parse tabular search results into structured objects
    - Extract and translate protein sequences from gene models
    - Format databases for searching
    - Support for parallel execution

Example:
    >>> from helixforge.homology.search import HomologySearch, SearchTool
    >>> searcher = HomologySearch(
    ...     tool=SearchTool.DIAMOND,
    ...     database="uniprot.dmnd",
    ...     threads=4,
    ... )
    >>> results = searcher.search_and_parse("proteins.fa")
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Iterator

import attrs

if TYPE_CHECKING:
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GFF3Parser, GeneModel

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

# Standard genetic code
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Diamond output format fields
DIAMOND_OUTFMT = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore",
    "qlen", "slen",
]


# =============================================================================
# Enums
# =============================================================================


class SearchTool(Enum):
    """Supported homology search tools."""

    DIAMOND = "diamond"
    MMSEQS2 = "mmseqs2"


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class HomologyHit:
    """Single homology search hit.

    Attributes:
        query_id: Gene/transcript ID.
        subject_id: Database protein ID.
        identity: Percent identity (0-100).
        alignment_length: Length of alignment in residues.
        mismatches: Number of mismatches.
        gap_opens: Number of gap openings.
        query_start: Start position in query (1-based).
        query_end: End position in query (1-based).
        subject_start: Start position in subject (1-based).
        subject_end: End position in subject (1-based).
        evalue: E-value of hit.
        bitscore: Bit score of hit.
        query_length: Total query length (if available).
        subject_length: Total subject length (if available).
    """

    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bitscore: float
    query_length: int | None = None
    subject_length: int | None = None

    @property
    def query_coverage(self) -> float | None:
        """Fraction of query aligned (0-1)."""
        if self.query_length and self.query_length > 0:
            aligned = self.query_end - self.query_start + 1
            return aligned / self.query_length
        return None

    @property
    def subject_coverage(self) -> float | None:
        """Fraction of subject aligned (0-1)."""
        if self.subject_length and self.subject_length > 0:
            aligned = self.subject_end - self.subject_start + 1
            return aligned / self.subject_length
        return None

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "query_id": self.query_id,
            "subject_id": self.subject_id,
            "identity": self.identity,
            "alignment_length": self.alignment_length,
            "mismatches": self.mismatches,
            "gap_opens": self.gap_opens,
            "query_start": self.query_start,
            "query_end": self.query_end,
            "subject_start": self.subject_start,
            "subject_end": self.subject_end,
            "evalue": self.evalue,
            "bitscore": self.bitscore,
            "query_length": self.query_length,
            "subject_length": self.subject_length,
            "query_coverage": self.query_coverage,
            "subject_coverage": self.subject_coverage,
        }

    @classmethod
    def from_diamond_line(
        cls,
        line: str,
        query_lengths: dict[str, int] | None = None,
    ) -> HomologyHit:
        """Parse Diamond tabular output line.

        Expects format 6 with columns: qseqid sseqid pident length
        mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

        Args:
            line: Tab-separated line from Diamond output.
            query_lengths: Optional {query_id: length} for fallback.

        Returns:
            HomologyHit instance.
        """
        fields = line.strip().split("\t")
        if len(fields) < 12:
            raise ValueError(f"Expected at least 12 fields, got {len(fields)}")

        query_id = fields[0]
        subject_id = fields[1]
        identity = float(fields[2])
        alignment_length = int(fields[3])
        mismatches = int(fields[4])
        gap_opens = int(fields[5])
        query_start = int(fields[6])
        query_end = int(fields[7])
        subject_start = int(fields[8])
        subject_end = int(fields[9])
        evalue = float(fields[10])
        bitscore = float(fields[11])

        # Optional length fields
        query_length = None
        subject_length = None

        if len(fields) >= 14:
            query_length = int(fields[12]) if fields[12] else None
            subject_length = int(fields[13]) if fields[13] else None

        # Fallback to provided lengths
        if query_length is None and query_lengths:
            query_length = query_lengths.get(query_id)

        return cls(
            query_id=query_id,
            subject_id=subject_id,
            identity=identity,
            alignment_length=alignment_length,
            mismatches=mismatches,
            gap_opens=gap_opens,
            query_start=query_start,
            query_end=query_end,
            subject_start=subject_start,
            subject_end=subject_end,
            evalue=evalue,
            bitscore=bitscore,
            query_length=query_length,
            subject_length=subject_length,
        )

    @classmethod
    def from_mmseqs_line(
        cls,
        line: str,
        query_lengths: dict[str, int] | None = None,
    ) -> HomologyHit:
        """Parse MMseqs2 tabular output line.

        Args:
            line: Tab-separated line from MMseqs2 output.
            query_lengths: Optional {query_id: length} for coverage.

        Returns:
            HomologyHit instance.
        """
        # MMseqs2 default format is similar to BLAST format 6
        return cls.from_diamond_line(line, query_lengths)


@attrs.define(slots=True)
class GeneHomology:
    """Aggregated homology results for a gene.

    Attributes:
        gene_id: Gene identifier.
        hits: List of all homology hits.
        best_hit_id: ID of best hit (by e-value).
        best_hit_evalue: E-value of best hit.
        best_hit_identity: Percent identity of best hit.
        best_hit_coverage: Query coverage of best hit.
        status: Validation status string.
        flags: List of issue flags.
    """

    gene_id: str
    hits: list[HomologyHit] = attrs.Factory(list)
    best_hit_id: str | None = None
    best_hit_evalue: float | None = None
    best_hit_identity: float | None = None
    best_hit_coverage: float | None = None
    status: str = "unknown"
    flags: list[str] = attrs.Factory(list)

    @classmethod
    def from_hits(cls, gene_id: str, hits: list[HomologyHit]) -> GeneHomology:
        """Create from list of hits, computing best hit stats.

        Args:
            gene_id: Gene identifier.
            hits: List of HomologyHit objects.

        Returns:
            GeneHomology with best hit statistics.
        """
        if not hits:
            return cls(gene_id=gene_id, hits=[], status="no_hit")

        # Sort by e-value (lower is better)
        sorted_hits = sorted(hits, key=lambda h: h.evalue)
        best = sorted_hits[0]

        return cls(
            gene_id=gene_id,
            hits=sorted_hits,
            best_hit_id=best.subject_id,
            best_hit_evalue=best.evalue,
            best_hit_identity=best.identity,
            best_hit_coverage=best.query_coverage,
        )

    @property
    def n_hits(self) -> int:
        """Number of hits."""
        return len(self.hits)

    @property
    def has_homology(self) -> bool:
        """Whether any hits were found."""
        return len(self.hits) > 0

    @property
    def unique_subjects(self) -> set[str]:
        """Set of unique subject protein IDs."""
        return {h.subject_id for h in self.hits}

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "gene_id": self.gene_id,
            "n_hits": self.n_hits,
            "has_homology": self.has_homology,
            "best_hit_id": self.best_hit_id,
            "best_hit_evalue": self.best_hit_evalue,
            "best_hit_identity": self.best_hit_identity,
            "best_hit_coverage": self.best_hit_coverage,
            "unique_subjects": list(self.unique_subjects),
            "status": self.status,
            "flags": self.flags,
        }


# =============================================================================
# Search Class
# =============================================================================


class HomologySearch:
    """Run homology searches against protein databases.

    Supports Diamond (recommended) and MMseqs2.

    Example:
        >>> searcher = HomologySearch(
        ...     tool=SearchTool.DIAMOND,
        ...     database="swissprot.dmnd",
        ...     threads=8,
        ... )
        >>> results = searcher.search_and_parse("proteins.fa")
        >>> for gene_id, hits in results.items():
        ...     print(f"{gene_id}: {len(hits)} hits")
    """

    def __init__(
        self,
        tool: SearchTool | str = SearchTool.DIAMOND,
        database: Path | str | None = None,
        threads: int = 4,
        evalue: float = 1e-5,
        max_target_seqs: int = 10,
        sensitive: bool = False,
        tmp_dir: Path | str | None = None,
    ) -> None:
        """Initialize search.

        Args:
            tool: Search tool to use (diamond or mmseqs2).
            database: Path to formatted database.
            threads: Number of threads.
            evalue: E-value threshold.
            max_target_seqs: Maximum hits per query.
            sensitive: Use sensitive search mode.
            tmp_dir: Temporary directory for intermediate files.
        """
        self.tool = SearchTool(tool) if isinstance(tool, str) else tool
        self.database = Path(database) if database else None
        self.threads = threads
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs
        self.sensitive = sensitive
        self.tmp_dir = Path(tmp_dir) if tmp_dir else Path(tempfile.gettempdir())

        self._check_tool_available()

    def _check_tool_available(self) -> None:
        """Verify search tool is installed.

        Raises:
            RuntimeError: If tool is not found in PATH.
        """
        if self.tool == SearchTool.DIAMOND:
            executable = "diamond"
        else:
            executable = "mmseqs"

        if shutil.which(executable) is None:
            raise RuntimeError(
                f"{self.tool.value} not found in PATH. "
                f"Please install {self.tool.value} and ensure it's in your PATH."
            )

    def format_database(
        self,
        fasta_path: Path | str,
        output_path: Path | str | None = None,
    ) -> Path:
        """Format protein database for searching.

        Args:
            fasta_path: Input protein FASTA.
            output_path: Output database path. Auto-generated if None.

        Returns:
            Path to formatted database.

        Raises:
            FileNotFoundError: If FASTA doesn't exist.
            RuntimeError: If formatting fails.
        """
        fasta_path = Path(fasta_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        if output_path is None:
            output_path = fasta_path.with_suffix("")
        output_path = Path(output_path)

        # Handle case where user passes a directory instead of a filename
        if output_path.is_dir() or str(output_path) in (".", "./", ".."):
            # Use input filename as base, place in the specified directory
            base_name = fasta_path.stem
            output_path = output_path / base_name
        elif output_path.name == "" or output_path.name == ".":
            raise ValueError(
                f"Invalid output path: '{output_path}'. "
                "Please specify a filename, not just a directory. "
                f"Example: --output {fasta_path.stem}"
            )

        logger.info(f"Formatting database from {fasta_path}")

        if self.tool == SearchTool.DIAMOND:
            cmd = [
                "diamond", "makedb",
                "--in", str(fasta_path),
                "--db", str(output_path),
                "--threads", str(self.threads),
            ]
            # Diamond adds .dmnd extension automatically
            expected_output = output_path.with_suffix(".dmnd")
        else:
            # MMseqs2
            cmd = [
                "mmseqs", "createdb",
                str(fasta_path),
                str(output_path),
            ]
            expected_output = output_path

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
            logger.debug(f"Database formatting output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Database formatting failed: {e.stderr}"
            ) from e

        if self.tool == SearchTool.DIAMOND:
            self.database = expected_output
        else:
            self.database = output_path

        logger.info(f"Database created: {self.database}")
        return self.database

    def search(
        self,
        query_fasta: Path | str,
        output_path: Path | str | None = None,
    ) -> Path:
        """Run homology search.

        Args:
            query_fasta: Query protein sequences.
            output_path: Output path for results. Auto-generated if None.

        Returns:
            Path to tabular output file.

        Raises:
            FileNotFoundError: If query or database doesn't exist.
            ValueError: If database not set.
            RuntimeError: If search fails.
        """
        query_fasta = Path(query_fasta)
        if not query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")

        if self.database is None:
            raise ValueError("Database not set. Call format_database() first.")

        if not self.database.exists():
            raise FileNotFoundError(f"Database not found: {self.database}")

        if output_path is None:
            output_path = self.tmp_dir / f"search_results_{query_fasta.stem}.tsv"
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Running {self.tool.value} search")
        logger.info(f"  Query: {query_fasta}")
        logger.info(f"  Database: {self.database}")

        if self.tool == SearchTool.DIAMOND:
            self._run_diamond(query_fasta, output_path)
        else:
            self._run_mmseqs(query_fasta, output_path)

        logger.info(f"Search results written to: {output_path}")
        return output_path

    def _run_diamond(self, query: Path, output: Path) -> None:
        """Execute Diamond blastp.

        Args:
            query: Query FASTA path.
            output: Output file path.
        """
        # Build command with outfmt fields as separate arguments
        # Newer Diamond versions require: --outfmt 6 field1 field2 ...
        # (not: --outfmt "6 field1 field2 ...")
        cmd = [
            "diamond", "blastp",
            "--db", str(self.database),
            "--query", str(query),
            "--out", str(output),
            "--outfmt", "6",
        ] + DIAMOND_OUTFMT + [
            "--threads", str(self.threads),
            "--evalue", str(self.evalue),
            "--max-target-seqs", str(self.max_target_seqs),
        ]

        if self.sensitive:
            cmd.append("--sensitive")

        logger.debug(f"Diamond command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
            logger.debug(f"Diamond output: {result.stderr}")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Diamond search failed: {e.stderr}") from e

    def _run_mmseqs(self, query: Path, output: Path) -> None:
        """Execute MMseqs2 search.

        Args:
            query: Query FASTA path.
            output: Output file path.
        """
        # MMseqs2 requires a temporary directory for intermediate files
        with tempfile.TemporaryDirectory(dir=self.tmp_dir) as tmp:
            tmp_path = Path(tmp)
            query_db = tmp_path / "query_db"
            result_db = tmp_path / "result_db"

            # Create query database
            subprocess.run(
                ["mmseqs", "createdb", str(query), str(query_db)],
                check=True,
                capture_output=True,
            )

            # Run search
            search_cmd = [
                "mmseqs", "search",
                str(query_db),
                str(self.database),
                str(result_db),
                str(tmp_path),
                "--threads", str(self.threads),
                "-e", str(self.evalue),
                "--max-seqs", str(self.max_target_seqs),
            ]

            if self.sensitive:
                search_cmd.extend(["-s", "7.5"])

            subprocess.run(search_cmd, check=True, capture_output=True)

            # Convert to tabular format
            subprocess.run(
                [
                    "mmseqs", "convertalis",
                    str(query_db),
                    str(self.database),
                    str(result_db),
                    str(output),
                    "--format-output",
                    "query,target,pident,alnlen,mismatch,gapopen,"
                    "qstart,qend,tstart,tend,evalue,bits,qlen,tlen",
                ],
                check=True,
                capture_output=True,
            )

    def parse_results(
        self,
        results_path: Path | str,
        query_lengths: dict[str, int] | None = None,
    ) -> dict[str, list[HomologyHit]]:
        """Parse tabular search results.

        Args:
            results_path: Path to search output.
            query_lengths: Optional {query_id: length} for coverage.

        Returns:
            {query_id: [HomologyHit, ...]}
        """
        results_path = Path(results_path)
        if not results_path.exists():
            raise FileNotFoundError(f"Results file not found: {results_path}")

        hits_by_query: dict[str, list[HomologyHit]] = {}

        with open(results_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                try:
                    if self.tool == SearchTool.DIAMOND:
                        hit = HomologyHit.from_diamond_line(line, query_lengths)
                    else:
                        hit = HomologyHit.from_mmseqs_line(line, query_lengths)

                    if hit.query_id not in hits_by_query:
                        hits_by_query[hit.query_id] = []
                    hits_by_query[hit.query_id].append(hit)

                except (ValueError, IndexError) as e:
                    logger.warning(f"Failed to parse line: {line[:50]}... ({e})")
                    continue

        logger.info(f"Parsed hits for {len(hits_by_query)} queries")
        return hits_by_query

    def search_and_parse(
        self,
        query_fasta: Path | str,
        query_lengths: dict[str, int] | None = None,
    ) -> dict[str, list[HomologyHit]]:
        """Convenience method: search and parse in one call.

        Args:
            query_fasta: Query protein FASTA.
            query_lengths: Optional {query_id: length} for coverage.

        Returns:
            {query_id: [HomologyHit, ...]}
        """
        output_path = self.search(query_fasta)
        return self.parse_results(output_path, query_lengths)


# =============================================================================
# Protein Extraction
# =============================================================================


def translate_sequence(
    cds_sequence: str,
    genetic_code: int = 1,
    stop_symbol: str = "*",
) -> str:
    """Translate CDS to protein sequence.

    Args:
        cds_sequence: Nucleotide sequence (CDS).
        genetic_code: NCBI genetic code number (1 = standard).
        stop_symbol: Character for stop codons.

    Returns:
        Protein sequence string.

    Note:
        Currently only supports standard genetic code (1).
        Partial codons at end are ignored.
        Internal stops are kept as specified symbol.
    """
    if genetic_code != 1:
        logger.warning(f"Genetic code {genetic_code} not fully supported, using standard")

    seq = cds_sequence.upper().replace("U", "T")
    protein = []

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if "N" in codon:
            protein.append("X")  # Unknown amino acid
        else:
            aa = CODON_TABLE.get(codon, "X")
            protein.append(aa if aa != "*" else stop_symbol)

    return "".join(protein)


def _get_cds_sequence(
    gene: GeneModel,
    genome: GenomeAccessor,
) -> str | None:
    """Extract CDS sequence for a gene.

    Args:
        gene: Gene model with CDS features.
        genome: Genome accessor for sequence extraction.

    Returns:
        CDS sequence string, or None if no CDS.
    """
    if not gene.transcripts:
        return None

    # Use primary (first) transcript
    transcript = gene.transcripts[0]

    if not transcript.cds:
        return None

    # Sort CDS by position
    cds_regions = sorted(transcript.cds, key=lambda x: x[0])

    # Concatenate CDS sequences
    cds_seqs = []
    for start, end, _phase in cds_regions:
        seq = genome.get_sequence(gene.seqid, start, end, gene.strand)
        cds_seqs.append(seq)

    # Handle strand
    if gene.strand == "-":
        # CDS already extracted in correct orientation by get_sequence
        # but we need to reverse the order of exons for minus strand
        cds_seqs = cds_seqs[::-1]

    return "".join(cds_seqs)


def extract_proteins_from_gff(
    gff_parser: GFF3Parser,
    genome: GenomeAccessor,
    output_fasta: Path | str,
    longest_isoform: bool = True,
) -> dict[str, int]:
    """Extract protein sequences from gene models.

    Args:
        gff_parser: Parsed GFF3 with gene models.
        genome: Genome accessor for sequence extraction.
        output_fasta: Output FASTA path.
        longest_isoform: If True, output only longest isoform per gene.

    Returns:
        {protein_id: length} mapping.
    """
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    protein_lengths: dict[str, int] = {}
    n_genes = 0
    n_extracted = 0
    n_no_cds = 0

    with open(output_fasta, "w") as f:
        for gene in gff_parser.iter_genes():
            n_genes += 1

            cds_seq = _get_cds_sequence(gene, genome)

            if cds_seq is None or len(cds_seq) < 3:
                n_no_cds += 1
                continue

            protein_seq = translate_sequence(cds_seq)

            # Remove terminal stop if present
            if protein_seq.endswith("*"):
                protein_seq = protein_seq[:-1]

            if not protein_seq:
                n_no_cds += 1
                continue

            # Use gene ID as protein ID
            protein_id = gene.gene_id

            f.write(f">{protein_id}\n")
            # Write in 60-character lines
            for i in range(0, len(protein_seq), 60):
                f.write(protein_seq[i:i + 60] + "\n")

            protein_lengths[protein_id] = len(protein_seq)
            n_extracted += 1

    logger.info(
        f"Extracted {n_extracted} proteins from {n_genes} genes "
        f"({n_no_cds} had no CDS)"
    )

    return protein_lengths


def iter_fasta_sequences(fasta_path: Path | str) -> Iterator[tuple[str, str]]:
    """Iterate over sequences in a FASTA file.

    Args:
        fasta_path: Path to FASTA file.

    Yields:
        (sequence_id, sequence) tuples.
    """
    fasta_path = Path(fasta_path)

    with open(fasta_path) as f:
        current_id = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            yield current_id, "".join(current_seq)


def get_sequence_lengths(fasta_path: Path | str) -> dict[str, int]:
    """Get lengths of all sequences in a FASTA file.

    Args:
        fasta_path: Path to FASTA file.

    Returns:
        {sequence_id: length}
    """
    lengths = {}
    for seq_id, seq in iter_fasta_sequences(fasta_path):
        lengths[seq_id] = len(seq)
    return lengths


# =============================================================================
# Legacy Function Wrappers
# =============================================================================


def run_diamond(
    query: Path | str,
    database: Path | str,
    output: Path | str | None = None,
    evalue: float = 1e-5,
    max_target_seqs: int = 5,
    threads: int = 1,
    sensitive: bool = False,
    block_size: float = 2.0,
) -> list[HomologyHit]:
    """Run DIAMOND BLASTP search (legacy wrapper).

    Args:
        query: Path to query FASTA file.
        database: Path to DIAMOND database (.dmnd).
        output: Path for output file. If None, uses temp file.
        evalue: E-value threshold.
        max_target_seqs: Maximum hits per query.
        threads: Number of threads to use.
        sensitive: Use sensitive mode.
        block_size: Memory block size in GB (unused).

    Returns:
        List of HomologyHit objects.
    """
    searcher = HomologySearch(
        tool=SearchTool.DIAMOND,
        database=database,
        threads=threads,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        sensitive=sensitive,
    )

    results = searcher.search_and_parse(query)

    # Flatten to list
    all_hits = []
    for hits in results.values():
        all_hits.extend(hits)

    return all_hits


def run_blast(
    query: Path | str,
    database: Path | str,
    output: Path | str | None = None,
    program: str = "blastp",
    evalue: float = 1e-5,
    max_target_seqs: int = 5,
    threads: int = 1,
) -> list[HomologyHit]:
    """Run NCBI BLAST search (stub - use Diamond instead).

    Args:
        query: Path to query FASTA file.
        database: Path to BLAST database.
        output: Path for output file.
        program: BLAST program.
        evalue: E-value threshold.
        max_target_seqs: Maximum hits per query.
        threads: Number of threads.

    Returns:
        List of HomologyHit objects.

    Note:
        BLAST is not fully implemented. Use Diamond instead.
    """
    raise NotImplementedError(
        "BLAST support not implemented. Use Diamond instead with run_diamond() "
        "or HomologySearch(tool=SearchTool.DIAMOND)"
    )


def parse_diamond_output(output_path: Path | str) -> list[HomologyHit]:
    """Parse DIAMOND output file (legacy wrapper).

    Args:
        output_path: Path to DIAMOND output file.

    Returns:
        List of HomologyHit objects.
    """
    searcher = HomologySearch(tool=SearchTool.DIAMOND)
    results = searcher.parse_results(output_path)

    all_hits = []
    for hits in results.values():
        all_hits.extend(hits)

    return all_hits


def build_diamond_database(
    fasta: Path | str,
    output: Path | str,
    threads: int = 1,
) -> Path:
    """Build a DIAMOND database from a FASTA file.

    Args:
        fasta: Path to input FASTA file.
        output: Path for output database.
        threads: Number of threads to use.

    Returns:
        Path to the created database.
    """
    searcher = HomologySearch(
        tool=SearchTool.DIAMOND,
        threads=threads,
    )
    return searcher.format_database(fasta, output)
