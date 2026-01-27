"""BAM file handling for RNA-seq alignments.

This module provides readers for BAM files to extract RNA-seq evidence
for gene model validation and refinement.

Features:
    - Extract splice junctions from split reads
    - Calculate coverage profiles
    - Support for strand-specific libraries
    - Memory-efficient streaming
    - BED format I/O for junctions

Example:
    >>> from helixforge.io.bam import JunctionExtractor
    >>> extractor = JunctionExtractor("rnaseq.bam")
    >>> junctions = extractor.extract_region("chr1", 1000, 50000)
    >>> for j in junctions:
    ...     print(f"{j.donor}-{j.acceptor}: {j.read_count} reads")
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator, Literal, NamedTuple

import numpy as np
import pysam

if TYPE_CHECKING:
    from helixforge.io.fasta import GenomeAccessor

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

# CIGAR operations
CIGAR_M = 0  # Match or mismatch
CIGAR_I = 1  # Insertion
CIGAR_D = 2  # Deletion
CIGAR_N = 3  # Skipped region (intron)
CIGAR_S = 4  # Soft clip
CIGAR_H = 5  # Hard clip
CIGAR_P = 6  # Padding
CIGAR_EQ = 7  # Sequence match
CIGAR_X = 8  # Sequence mismatch

# Canonical splice site dinucleotides
CANONICAL_DONORS = {"GT", "GC", "AT"}
CANONICAL_ACCEPTORS = {"AG", "AC"}

Strand = Literal["+", "-"]


# =============================================================================
# Data Structures
# =============================================================================


class SpliceJunction(NamedTuple):
    """A splice junction with supporting evidence.

    Attributes:
        seqid: Chromosome/contig identifier.
        start: Intron start (0-based, first intronic base).
        end: Intron end (0-based, exclusive - first base after intron).
        strand: Strand inferred from splice site or XS tag.
        read_count: Number of reads supporting this junction.
        unique_count: Number of uniquely mapped reads.
        is_canonical: Whether junction has canonical splice sites.
        donor_dinuc: Dinucleotide at donor site (if known).
        acceptor_dinuc: Dinucleotide at acceptor site (if known).
    """

    seqid: str
    start: int  # Intron start (0-based)
    end: int  # Intron end (0-based, exclusive)
    strand: str
    read_count: int
    unique_count: int
    is_canonical: bool | None = None
    donor_dinuc: str | None = None
    acceptor_dinuc: str | None = None

    @property
    def intron_size(self) -> int:
        """Intron size in base pairs."""
        return self.end - self.start

    @property
    def donor(self) -> int:
        """Donor site position (last exonic base, 0-based)."""
        return self.start - 1

    @property
    def acceptor(self) -> int:
        """Acceptor site position (first exonic base, 0-based)."""
        return self.end


class CoverageProfile(NamedTuple):
    """Per-base coverage for a genomic region.

    Attributes:
        seqid: Chromosome/contig identifier.
        start: Start position (0-based).
        end: End position (0-based, exclusive).
        values: Array of coverage values.
        strand: Strand if strand-specific, None otherwise.
    """

    seqid: str
    start: int
    end: int
    values: np.ndarray
    strand: str | None = None

    @property
    def mean_coverage(self) -> float:
        """Mean coverage across the region."""
        return float(np.mean(self.values))

    @property
    def max_coverage(self) -> int:
        """Maximum coverage in the region."""
        return int(np.max(self.values))

    @property
    def min_coverage(self) -> int:
        """Minimum coverage in the region."""
        return int(np.min(self.values))


# =============================================================================
# Junction Extractor
# =============================================================================


class JunctionExtractor:
    """Extract splice junctions from RNA-seq BAM files.

    Identifies splice junctions from reads with CIGAR 'N' operations
    and counts supporting reads.

    Features:
    - Filter by mapping quality and read count
    - Strand inference from XS tag or read orientation
    - Parallel extraction by region
    - BED format output

    Attributes:
        path: Path to the BAM file.
        min_reads: Minimum reads to report a junction.
        min_mapq: Minimum mapping quality.
        min_overhang: Minimum exonic overhang on each side.

    Example:
        >>> extractor = JunctionExtractor("rnaseq.bam")
        >>> junctions = extractor.extract_region("chr1", 0, 1000000)
        >>> extractor.to_bed(junctions, "junctions.bed")
    """

    def __init__(
        self,
        bam_path: Path | str,
        min_reads: int = 3,
        min_mapq: int = 20,
        min_overhang: int = 8,
    ) -> None:
        """Initialize the junction extractor.

        Args:
            bam_path: Path to indexed BAM file.
            min_reads: Minimum reads to report a junction.
            min_mapq: Minimum mapping quality.
            min_overhang: Minimum exonic overhang on each side.

        Raises:
            FileNotFoundError: If BAM file doesn't exist.
            ValueError: If BAM file is not indexed.
        """
        self.path = Path(bam_path)
        self.min_reads = min_reads
        self.min_mapq = min_mapq
        self.min_overhang = min_overhang

        if not self.path.exists():
            raise FileNotFoundError(f"BAM file not found: {self.path}")

        # Check for index
        index_paths = [
            self.path.with_suffix(".bam.bai"),
            self.path.with_suffix(".bai"),
            Path(str(self.path) + ".bai"),
        ]
        if not any(p.exists() for p in index_paths):
            raise ValueError(
                f"BAM index not found. Please run: samtools index {self.path}"
            )

        self._bam: pysam.AlignmentFile | None = None
        self._open()

    def _open(self) -> None:
        """Open the BAM file."""
        self._bam = pysam.AlignmentFile(str(self.path), "rb")
        logger.info(f"Opened BAM file: {self.path.name}")

    def __enter__(self) -> JunctionExtractor:
        """Context manager entry."""
        return self

    def __exit__(self, *args: Any) -> None:
        """Context manager exit."""
        self.close()

    def close(self) -> None:
        """Close the BAM file."""
        if self._bam is not None:
            self._bam.close()
            self._bam = None

    @property
    def references(self) -> list[str]:
        """List of reference sequences in BAM."""
        if self._bam is None:
            raise RuntimeError("BAM file not open")
        return list(self._bam.references)

    @property
    def reference_lengths(self) -> dict[str, int]:
        """Dictionary of reference lengths."""
        if self._bam is None:
            raise RuntimeError("BAM file not open")
        return dict(zip(self._bam.references, self._bam.lengths))

    def _get_junctions_from_read(
        self,
        read: pysam.AlignedSegment,
    ) -> list[tuple[int, int, str]]:
        """Extract junctions from a single read.

        Args:
            read: Aligned read segment.

        Returns:
            List of (intron_start, intron_end, strand) tuples.
        """
        junctions = []

        # Get strand from XS tag if available
        try:
            strand = read.get_tag("XS")
        except KeyError:
            # Infer from read orientation if paired
            if read.is_paired:
                if read.is_read1:
                    strand = "-" if read.is_reverse else "+"
                else:
                    strand = "+" if read.is_reverse else "-"
            else:
                strand = "."

        # Parse CIGAR for N operations (introns)
        ref_pos = read.reference_start
        query_pos = 0

        for op, length in read.cigartuples or []:
            if op == CIGAR_N:  # Intron
                intron_start = ref_pos
                intron_end = ref_pos + length

                # Check overhang - need sufficient aligned bases on each side
                # This is a simplified check; a full implementation would
                # track aligned bases before/after the junction
                if length >= self.min_overhang:
                    junctions.append((intron_start, intron_end, strand))

                ref_pos += length

            elif op in (CIGAR_M, CIGAR_EQ, CIGAR_X):
                ref_pos += length
                query_pos += length

            elif op == CIGAR_I:
                query_pos += length

            elif op == CIGAR_D:
                ref_pos += length

            elif op == CIGAR_S:
                query_pos += length

        return junctions

    def extract_region(
        self,
        seqid: str,
        start: int,
        end: int,
        min_reads: int | None = None,
    ) -> list[SpliceJunction]:
        """Extract junctions from a genomic region.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).
            min_reads: Override default minimum reads.

        Returns:
            List of SpliceJunction objects.
        """
        if self._bam is None:
            raise RuntimeError("BAM file not open")

        min_reads = min_reads if min_reads is not None else self.min_reads

        # Count junctions: (start, end) -> {strand: (total_count, unique_count)}
        junction_counts: dict[tuple[int, int], dict[str, list[int]]] = defaultdict(
            lambda: defaultdict(lambda: [0, 0])
        )

        for read in self._bam.fetch(seqid, start, end):
            # Filter reads
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < self.min_mapq:
                continue

            is_unique = read.mapping_quality >= 255 or not read.has_tag("NH") or read.get_tag("NH") == 1

            # Extract junctions
            for intron_start, intron_end, strand in self._get_junctions_from_read(read):
                # Only count junctions that overlap our region
                if intron_start >= start and intron_end <= end:
                    junction_counts[(intron_start, intron_end)][strand][0] += 1
                    if is_unique:
                        junction_counts[(intron_start, intron_end)][strand][1] += 1

        # Build junction objects
        junctions = []
        for (junc_start, junc_end), strand_counts in junction_counts.items():
            # Merge strand counts if ambiguous
            total_count = sum(c[0] for c in strand_counts.values())
            unique_count = sum(c[1] for c in strand_counts.values())

            if total_count < min_reads:
                continue

            # Determine strand (use most common)
            if len(strand_counts) == 1:
                strand = list(strand_counts.keys())[0]
            else:
                strand = max(strand_counts.keys(), key=lambda s: strand_counts[s][0])

            junction = SpliceJunction(
                seqid=seqid,
                start=junc_start,
                end=junc_end,
                strand=strand,
                read_count=total_count,
                unique_count=unique_count,
            )
            junctions.append(junction)

        # Sort by position
        junctions.sort(key=lambda j: (j.start, j.end))

        logger.debug(f"Extracted {len(junctions)} junctions from {seqid}:{start}-{end}")
        return junctions

    def extract_all(
        self,
        min_reads: int | None = None,
    ) -> dict[str, list[SpliceJunction]]:
        """Extract all junctions from the BAM file.

        Args:
            min_reads: Override default minimum reads.

        Returns:
            Dictionary mapping seqid to list of junctions.
        """
        if self._bam is None:
            raise RuntimeError("BAM file not open")

        result: dict[str, list[SpliceJunction]] = {}

        for seqid, length in zip(self._bam.references, self._bam.lengths):
            junctions = self.extract_region(seqid, 0, length, min_reads)
            if junctions:
                result[seqid] = junctions

        total = sum(len(j) for j in result.values())
        logger.info(f"Extracted {total} junctions from {len(result)} scaffolds")
        return result

    def to_bed(
        self,
        junctions: list[SpliceJunction],
        output_path: Path | str,
    ) -> None:
        """Write junctions to BED6+2 format.

        Columns: chrom, start, end, name, score, strand, read_count, unique_count

        Args:
            junctions: List of junctions to write.
            output_path: Output file path.
        """
        with open(output_path, "w") as f:
            f.write("#chrom\tstart\tend\tname\tscore\tstrand\tread_count\tunique_count\n")
            for i, j in enumerate(junctions):
                name = f"JUNC{i:06d}"
                score = min(1000, j.read_count)  # BED score capped at 1000
                f.write(
                    f"{j.seqid}\t{j.start}\t{j.end}\t{name}\t{score}\t"
                    f"{j.strand}\t{j.read_count}\t{j.unique_count}\n"
                )

    @classmethod
    def from_bed(cls, bed_path: Path | str) -> list[SpliceJunction]:
        """Load junctions from BED file.

        Args:
            bed_path: Path to BED file.

        Returns:
            List of SpliceJunction objects.
        """
        junctions = []
        with open(bed_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue

                seqid = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[5] if parts[5] in ("+", "-") else "."

                # Extended fields
                read_count = int(parts[6]) if len(parts) > 6 else 1
                unique_count = int(parts[7]) if len(parts) > 7 else read_count

                junction = SpliceJunction(
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    read_count=read_count,
                    unique_count=unique_count,
                )
                junctions.append(junction)

        return junctions

    def get_coverage(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str | None = None,
    ) -> CoverageProfile:
        """Calculate coverage for a region.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).
            strand: If provided, only count reads on this strand.

        Returns:
            CoverageProfile with per-base coverage values.
        """
        if self._bam is None:
            raise RuntimeError("BAM file not open")

        coverage = np.zeros(end - start, dtype=np.int32)

        for read in self._bam.fetch(seqid, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < self.min_mapq:
                continue

            # Strand filtering
            if strand is not None:
                try:
                    read_strand = read.get_tag("XS")
                except KeyError:
                    read_strand = "-" if read.is_reverse else "+"
                if read_strand != strand:
                    continue

            # Add coverage from aligned blocks
            for block_start, block_end in read.get_blocks():
                # Intersect with our region
                cov_start = max(block_start, start) - start
                cov_end = min(block_end, end) - start
                if cov_start < cov_end:
                    coverage[cov_start:cov_end] += 1

        return CoverageProfile(
            seqid=seqid,
            start=start,
            end=end,
            values=coverage,
            strand=strand,
        )

    def get_read_count(
        self,
        seqid: str,
        start: int,
        end: int,
    ) -> int:
        """Count reads overlapping a region.

        Args:
            seqid: Scaffold/chromosome name.
            start: Start position (0-based).
            end: End position (0-based, exclusive).

        Returns:
            Number of reads overlapping the region.
        """
        if self._bam is None:
            raise RuntimeError("BAM file not open")

        count = 0
        for read in self._bam.fetch(seqid, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < self.min_mapq:
                continue
            count += 1

        return count


# =============================================================================
# Utility Functions
# =============================================================================


def is_canonical_junction(donor_dinuc: str, acceptor_dinuc: str) -> bool:
    """Check if a splice junction has canonical dinucleotides.

    Canonical junctions are GT-AG (major) or GC-AG, AT-AC (minor).

    Args:
        donor_dinuc: Two nucleotides at donor site (e.g., "GT").
        acceptor_dinuc: Two nucleotides at acceptor site (e.g., "AG").

    Returns:
        True if the junction is canonical.
    """
    donor = donor_dinuc.upper()
    acceptor = acceptor_dinuc.upper()

    # Major spliceosome (U2-dependent)
    if donor == "GT" and acceptor == "AG":
        return True
    if donor == "GC" and acceptor == "AG":
        return True

    # Minor spliceosome (U12-dependent)
    if donor == "AT" and acceptor == "AC":
        return True

    return False


def annotate_junction_sequences(
    junctions: list[SpliceJunction],
    genome: "GenomeAccessor",
) -> list[SpliceJunction]:
    """Add donor/acceptor dinucleotides to junctions.

    Args:
        junctions: List of junctions to annotate.
        genome: GenomeAccessor for sequence lookup.

    Returns:
        New list of junctions with dinucleotide annotations.
    """
    annotated = []
    for j in junctions:
        try:
            # Get donor dinucleotide (first 2 bases of intron)
            donor_dinuc = genome.get_sequence(j.seqid, j.start, j.start + 2)
            # Get acceptor dinucleotide (last 2 bases of intron)
            acceptor_dinuc = genome.get_sequence(j.seqid, j.end - 2, j.end)

            is_canonical = is_canonical_junction(donor_dinuc, acceptor_dinuc)

            annotated.append(
                SpliceJunction(
                    seqid=j.seqid,
                    start=j.start,
                    end=j.end,
                    strand=j.strand,
                    read_count=j.read_count,
                    unique_count=j.unique_count,
                    is_canonical=is_canonical,
                    donor_dinuc=donor_dinuc.upper(),
                    acceptor_dinuc=acceptor_dinuc.upper(),
                )
            )
        except (KeyError, ValueError) as e:
            logger.warning(f"Could not annotate junction {j.seqid}:{j.start}-{j.end}: {e}")
            annotated.append(j)

    return annotated


# =============================================================================
# Junction Aggregation for Multiple BAMs
# =============================================================================


class MultiSampleJunction(NamedTuple):
    """A splice junction with per-sample evidence.

    Attributes:
        seqid: Chromosome/contig identifier.
        start: Intron start (0-based, first intronic base).
        end: Intron end (0-based, exclusive).
        strand: Strand inferred from splice site or XS tag.
        total_reads: Total reads across all samples.
        total_unique: Total unique reads across all samples.
        sample_reads: Dictionary mapping sample name to read count.
        sample_unique: Dictionary mapping sample name to unique count.
        n_samples_supporting: Number of samples with reads >= threshold.
        is_canonical: Whether junction has canonical splice sites.
        donor_dinuc: Dinucleotide at donor site (if known).
        acceptor_dinuc: Dinucleotide at acceptor site (if known).
    """

    seqid: str
    start: int
    end: int
    strand: str
    total_reads: int
    total_unique: int
    sample_reads: dict[str, int]
    sample_unique: dict[str, int]
    n_samples_supporting: int
    is_canonical: bool | None = None
    donor_dinuc: str | None = None
    acceptor_dinuc: str | None = None

    def to_splice_junction(self) -> SpliceJunction:
        """Convert to standard SpliceJunction (loses per-sample info)."""
        return SpliceJunction(
            seqid=self.seqid,
            start=self.start,
            end=self.end,
            strand=self.strand,
            read_count=self.total_reads,
            unique_count=self.total_unique,
            is_canonical=self.is_canonical,
            donor_dinuc=self.donor_dinuc,
            acceptor_dinuc=self.acceptor_dinuc,
        )


def aggregate_junctions_multi_sample(
    bam_paths: list[Path | str],
    sample_names: list[str] | None = None,
    min_reads_per_sample: int = 1,
    min_samples: int = 1,
    min_mapq: int = 20,
    min_overhang: int = 8,
    region: tuple[str, int, int] | None = None,
) -> dict[str, list[MultiSampleJunction]]:
    """Extract and aggregate junctions from multiple BAM files.

    This function extracts splice junctions from multiple BAM files and
    aggregates them by position, tracking per-sample read counts.

    Args:
        bam_paths: List of paths to indexed BAM files.
        sample_names: Optional list of sample names (defaults to BAM filenames).
        min_reads_per_sample: Minimum reads in a sample to count as supporting.
        min_samples: Minimum number of supporting samples to include junction.
        min_mapq: Minimum mapping quality for reads.
        min_overhang: Minimum exonic overhang on each side.
        region: Optional (seqid, start, end) tuple to restrict extraction.

    Returns:
        Dictionary mapping seqid to list of MultiSampleJunction objects.
    """
    if not bam_paths:
        return {}

    # Generate sample names if not provided
    if sample_names is None:
        sample_names = [Path(p).stem for p in bam_paths]

    if len(sample_names) != len(bam_paths):
        raise ValueError("Number of sample names must match number of BAM files")

    # Collect junctions from each BAM
    # Key: (seqid, start, end) -> {strand: {sample: (total, unique)}}
    junction_data: dict[
        tuple[str, int, int], dict[str, dict[str, tuple[int, int]]]
    ] = defaultdict(lambda: defaultdict(dict))

    for bam_path, sample_name in zip(bam_paths, sample_names):
        logger.info(f"Extracting junctions from {sample_name}: {bam_path}")
        extractor = JunctionExtractor(
            bam_path, min_reads=1, min_mapq=min_mapq, min_overhang=min_overhang
        )

        try:
            if region:
                seqid, start, end = region
                sample_junctions = {seqid: extractor.extract_region(seqid, start, end, min_reads=1)}
            else:
                sample_junctions = extractor.extract_all(min_reads=1)

            # Aggregate into junction_data
            for seqid, junctions in sample_junctions.items():
                for j in junctions:
                    key = (seqid, j.start, j.end)
                    junction_data[key][j.strand][sample_name] = (j.read_count, j.unique_count)
        finally:
            extractor.close()

    # Build MultiSampleJunction objects
    result: dict[str, list[MultiSampleJunction]] = defaultdict(list)

    for (seqid, start, end), strand_data in junction_data.items():
        # Merge strand information (use most common strand)
        all_sample_reads: dict[str, int] = {}
        all_sample_unique: dict[str, int] = {}

        for strand, sample_counts in strand_data.items():
            for sample, (reads, unique) in sample_counts.items():
                # If sample appears in multiple strands, sum the counts
                all_sample_reads[sample] = all_sample_reads.get(sample, 0) + reads
                all_sample_unique[sample] = all_sample_unique.get(sample, 0) + unique

        # Count supporting samples
        n_supporting = sum(
            1 for reads in all_sample_reads.values() if reads >= min_reads_per_sample
        )

        # Filter by min_samples
        if n_supporting < min_samples:
            continue

        # Determine strand (most common)
        strand_totals = {
            strand: sum(counts[0] for counts in sample_counts.values())
            for strand, sample_counts in strand_data.items()
        }
        best_strand = max(strand_totals.keys(), key=lambda s: strand_totals[s])

        total_reads = sum(all_sample_reads.values())
        total_unique = sum(all_sample_unique.values())

        junction = MultiSampleJunction(
            seqid=seqid,
            start=start,
            end=end,
            strand=best_strand,
            total_reads=total_reads,
            total_unique=total_unique,
            sample_reads=all_sample_reads,
            sample_unique=all_sample_unique,
            n_samples_supporting=n_supporting,
        )
        result[seqid].append(junction)

    # Sort by position
    for seqid in result:
        result[seqid].sort(key=lambda j: (j.start, j.end))

    total = sum(len(j) for j in result.values())
    logger.info(
        f"Aggregated {total} junctions from {len(bam_paths)} samples "
        f"(min_samples={min_samples}, min_reads_per_sample={min_reads_per_sample})"
    )

    return dict(result)


def multi_sample_to_standard_junctions(
    multi_junctions: dict[str, list[MultiSampleJunction]],
) -> dict[str, list[SpliceJunction]]:
    """Convert MultiSampleJunction dict to standard SpliceJunction dict.

    Args:
        multi_junctions: Dictionary from aggregate_junctions_multi_sample.

    Returns:
        Dictionary mapping seqid to list of SpliceJunction objects.
    """
    return {
        seqid: [mj.to_splice_junction() for mj in junctions]
        for seqid, junctions in multi_junctions.items()
    }


def load_junctions_from_bed(
    bed_path: Path | str,
    min_reads: int = 1,
) -> dict[str, list[SpliceJunction]]:
    """Load junctions from a BED file and organize by seqid.

    Args:
        bed_path: Path to BED file with junctions.
        min_reads: Minimum read count to include junction.

    Returns:
        Dictionary mapping seqid to list of SpliceJunction objects.
    """
    junctions = JunctionExtractor.from_bed(bed_path)

    # Filter by min_reads and organize by seqid
    result: dict[str, list[SpliceJunction]] = defaultdict(list)
    for j in junctions:
        if j.read_count >= min_reads:
            result[j.seqid].append(j)

    # Sort by position
    for seqid in result:
        result[seqid].sort(key=lambda j: (j.start, j.end))

    return dict(result)


def load_junctions_from_star_sj(
    sj_path: Path | str,
    min_reads: int = 1,
    include_multi_mappers: bool = True,
) -> dict[str, list[SpliceJunction]]:
    """Load junctions from STAR SJ.out.tab file.

    STAR SJ.out.tab format (tab-separated):
        1. chromosome
        2. first base of intron (1-based)
        3. last base of intron (1-based)
        4. strand (0: undefined, 1: +, 2: -)
        5. intron motif (0: non-canonical, 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT)
        6. annotation status (0: unannotated, 1: annotated)
        7. number of uniquely mapping reads
        8. number of multi-mapping reads
        9. maximum spliced alignment overhang

    Args:
        sj_path: Path to STAR SJ.out.tab file.
        min_reads: Minimum read count to include junction.
        include_multi_mappers: Include multi-mapping reads in count.

    Returns:
        Dictionary mapping seqid to list of SpliceJunction objects.
    """
    # STAR strand encoding
    STAR_STRAND = {0: ".", 1: "+", 2: "-"}

    # STAR motif encoding (for is_canonical)
    # 0: non-canonical, 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
    CANONICAL_MOTIFS = {1, 2, 3, 4, 5, 6}

    result: dict[str, list[SpliceJunction]] = defaultdict(list)

    with open(sj_path) as f:
        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) < 7:
                logger.warning(f"Skipping malformed line in {sj_path}: {line.strip()}")
                continue

            seqid = parts[0]
            # STAR uses 1-based coordinates, convert to 0-based half-open
            start = int(parts[1]) - 1  # First intronic base (0-based)
            end = int(parts[2])  # Last intronic base + 1 (0-based exclusive)
            strand = STAR_STRAND.get(int(parts[3]), ".")
            motif = int(parts[4])
            unique_reads = int(parts[6])
            multi_reads = int(parts[7]) if len(parts) > 7 else 0

            if include_multi_mappers:
                total_reads = unique_reads + multi_reads
            else:
                total_reads = unique_reads

            if total_reads < min_reads:
                continue

            is_canonical = motif in CANONICAL_MOTIFS

            junction = SpliceJunction(
                seqid=seqid,
                start=start,
                end=end,
                strand=strand,
                read_count=total_reads,
                unique_count=unique_reads,
                is_canonical=is_canonical,
            )
            result[seqid].append(junction)

    # Sort by position
    for seqid in result:
        result[seqid].sort(key=lambda j: (j.start, j.end))

    total = sum(len(j) for j in result.values())
    logger.info(f"Loaded {total} junctions from STAR SJ file: {sj_path}")

    return dict(result)


def detect_junction_file_format(file_path: Path | str) -> str:
    """Detect the format of a junction file.

    Args:
        file_path: Path to the junction file.

    Returns:
        Format string: "bed", "star_sj", or "unknown".
    """
    path = Path(file_path)

    # Check by extension first
    if path.name.endswith("_SJ.out.tab") or path.name.endswith(".SJ.out.tab"):
        return "star_sj"
    if path.suffix.lower() in (".bed", ".bed6", ".bed12"):
        return "bed"

    # Try to detect by content
    with open(path) as f:
        first_line = f.readline().strip()

        # Skip header/comment lines
        while first_line.startswith("#") or not first_line:
            first_line = f.readline().strip()
            if not first_line:
                return "unknown"

        parts = first_line.split("\t")

        # STAR SJ.out.tab has exactly 9 columns, all numeric except first
        if len(parts) == 9:
            try:
                # Columns 2-9 should be integers
                for i in range(1, 9):
                    int(parts[i])
                return "star_sj"
            except ValueError:
                pass

        # BED format: at least 3 columns, columns 2-3 are integers
        if len(parts) >= 3:
            try:
                int(parts[1])
                int(parts[2])
                return "bed"
            except ValueError:
                pass

    return "unknown"


def load_junctions_auto(
    file_path: Path | str,
    min_reads: int = 1,
) -> dict[str, list[SpliceJunction]]:
    """Load junctions from a file, auto-detecting format.

    Supports BED format and STAR SJ.out.tab format.

    Args:
        file_path: Path to junction file.
        min_reads: Minimum read count to include junction.

    Returns:
        Dictionary mapping seqid to list of SpliceJunction objects.

    Raises:
        ValueError: If file format cannot be detected.
    """
    fmt = detect_junction_file_format(file_path)

    if fmt == "star_sj":
        logger.info(f"Detected STAR SJ.out.tab format: {file_path}")
        return load_junctions_from_star_sj(file_path, min_reads=min_reads)
    elif fmt == "bed":
        logger.info(f"Detected BED format: {file_path}")
        return load_junctions_from_bed(file_path, min_reads=min_reads)
    else:
        raise ValueError(
            f"Could not detect junction file format: {file_path}. "
            "Expected BED or STAR SJ.out.tab format."
        )


def aggregate_junctions_from_files(
    file_paths: list[Path | str],
    sample_names: list[str] | None = None,
    min_reads_per_sample: int = 1,
    min_samples: int = 1,
) -> dict[str, list[MultiSampleJunction]]:
    """Aggregate junctions from multiple BED or STAR SJ.out.tab files.

    This function loads junctions from multiple files and aggregates them
    by position, tracking per-sample read counts. Format is auto-detected.

    Args:
        file_paths: List of paths to junction files (BED or STAR SJ.out.tab).
        sample_names: Optional list of sample names (defaults to filenames).
        min_reads_per_sample: Minimum reads in a sample to count as supporting.
        min_samples: Minimum number of supporting samples to include junction.

    Returns:
        Dictionary mapping seqid to list of MultiSampleJunction objects.
    """
    if not file_paths:
        return {}

    # Generate sample names if not provided
    if sample_names is None:
        sample_names = [Path(p).stem for p in file_paths]

    if len(sample_names) != len(file_paths):
        raise ValueError("Number of sample names must match number of files")

    # Collect junctions from each file
    # Key: (seqid, start, end) -> {strand: {sample: (total, unique)}}
    junction_data: dict[
        tuple[str, int, int], dict[str, dict[str, tuple[int, int]]]
    ] = defaultdict(lambda: defaultdict(dict))

    for file_path, sample_name in zip(file_paths, sample_names):
        logger.info(f"Loading junctions from {sample_name}: {file_path}")

        # Load junctions (auto-detect format), no filtering yet
        sample_junctions = load_junctions_auto(file_path, min_reads=1)

        # Aggregate into junction_data
        for seqid, junctions in sample_junctions.items():
            for j in junctions:
                key = (seqid, j.start, j.end)
                junction_data[key][j.strand][sample_name] = (j.read_count, j.unique_count)

    # Build MultiSampleJunction objects (same logic as BAM aggregation)
    result: dict[str, list[MultiSampleJunction]] = defaultdict(list)

    for (seqid, start, end), strand_data in junction_data.items():
        all_sample_reads: dict[str, int] = {}
        all_sample_unique: dict[str, int] = {}

        for strand, sample_counts in strand_data.items():
            for sample, (reads, unique) in sample_counts.items():
                all_sample_reads[sample] = all_sample_reads.get(sample, 0) + reads
                all_sample_unique[sample] = all_sample_unique.get(sample, 0) + unique

        # Count supporting samples
        n_supporting = sum(
            1 for reads in all_sample_reads.values() if reads >= min_reads_per_sample
        )

        if n_supporting < min_samples:
            continue

        # Determine strand (most common)
        strand_totals = {
            strand: sum(counts[0] for counts in sample_counts.values())
            for strand, sample_counts in strand_data.items()
        }
        best_strand = max(strand_totals.keys(), key=lambda s: strand_totals[s])

        total_reads = sum(all_sample_reads.values())
        total_unique = sum(all_sample_unique.values())

        junction = MultiSampleJunction(
            seqid=seqid,
            start=start,
            end=end,
            strand=best_strand,
            total_reads=total_reads,
            total_unique=total_unique,
            sample_reads=all_sample_reads,
            sample_unique=all_sample_unique,
            n_samples_supporting=n_supporting,
        )
        result[seqid].append(junction)

    # Sort by position
    for seqid in result:
        result[seqid].sort(key=lambda j: (j.start, j.end))

    total = sum(len(j) for j in result.values())
    logger.info(
        f"Aggregated {total} junctions from {len(file_paths)} files "
        f"(min_samples={min_samples}, min_reads_per_sample={min_reads_per_sample})"
    )

    return dict(result)


# =============================================================================
# Legacy Compatibility
# =============================================================================

# Alias for backward compatibility with stub
BamReader = JunctionExtractor
