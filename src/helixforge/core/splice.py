"""Splice site refinement using RNA-seq evidence and PWM scoring.

This module provides tools for refining splice sites in Helixer gene predictions
using empirical RNA-seq junction evidence and position weight matrix scoring.

Key components:
- SpliceSiteType: Classification of splice site dinucleotides
- SpliceSite: Individual splice site representation
- IntronModel: Intron with donor/acceptor sites
- PositionWeightMatrix: PWM for scoring splice site motifs
- SpliceRefiner: Main refinement engine
- SpliceCorrection: Records corrections made to splice sites
- GeneSpliceReport: Per-gene refinement report
- SpliceReportWriter: Output report generation

Example:
    >>> from helixforge.core.splice import SpliceRefiner, PositionWeightMatrix
    >>> from helixforge.io.fasta import GenomeAccessor
    >>> from helixforge.io.bam import JunctionExtractor
    >>>
    >>> genome = GenomeAccessor("genome.fa")
    >>> extractor = JunctionExtractor("rnaseq.bam")
    >>> junctions = extractor.extract_all()
    >>>
    >>> refiner = SpliceRefiner(genome, junctions)
    >>> refined_gene, report = refiner.refine_gene(gene)
"""

from __future__ import annotations

import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Iterator, Literal

import attrs
import numpy as np

if TYPE_CHECKING:
    from helixforge.io.bam import SpliceJunction
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GeneModel, TranscriptModel

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

# Canonical splice site dinucleotides
DONOR_CANONICAL = {"GT", "GC"}  # 5' splice site (U2 spliceosome)
ACCEPTOR_CANONICAL = {"AG"}  # 3' splice site

# Minor spliceosome (U12-dependent)
DONOR_U12 = {"AT"}
ACCEPTOR_U12 = {"AC"}

# Nucleotide to index mapping for PWM
BASE_TO_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
INDEX_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

# Default background frequencies (approximate plant GC content ~35%)
DEFAULT_BACKGROUND = np.array([0.33, 0.17, 0.17, 0.33])

# Default thresholds
DEFAULT_MAX_SHIFT = 15
DEFAULT_MIN_JUNCTION_READS = 3
DEFAULT_PWM_THRESHOLD = 0.0


# =============================================================================
# Enums
# =============================================================================


class SpliceSiteType(Enum):
    """Splice site dinucleotide classification."""

    GT_AG = "GT-AG"  # Canonical (~98.5% of plant introns)
    GC_AG = "GC-AG"  # Minor canonical (~1%)
    AT_AC = "AT-AC"  # U12-type minor spliceosome (~0.5%)
    OTHER = "other"  # Non-canonical (likely error)


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class SpliceSite:
    """Represents a single splice site (donor or acceptor).

    Attributes:
        seqid: Scaffold/chromosome name.
        position: 0-based, points to first base of dinucleotide.
        site_type: "donor" or "acceptor".
        strand: Strand (+ or -).
        dinucleotide: GT, GC, AT for donor; AG, AC for acceptor.
    """

    seqid: str
    position: int
    site_type: Literal["donor", "acceptor"]
    strand: str
    dinucleotide: str

    @property
    def splice_site_type(self) -> SpliceSiteType:
        """Classify based on dinucleotide."""
        dinuc = self.dinucleotide.upper()

        if self.site_type == "donor":
            if dinuc == "GT":
                return SpliceSiteType.GT_AG
            elif dinuc == "GC":
                return SpliceSiteType.GC_AG
            elif dinuc == "AT":
                return SpliceSiteType.AT_AC
            else:
                return SpliceSiteType.OTHER
        else:  # acceptor
            if dinuc == "AG":
                # Could be GT-AG or GC-AG; return GT-AG as default
                return SpliceSiteType.GT_AG
            elif dinuc == "AC":
                return SpliceSiteType.AT_AC
            else:
                return SpliceSiteType.OTHER

    @property
    def is_canonical(self) -> bool:
        """Check if this is a canonical splice site."""
        return self.splice_site_type in (
            SpliceSiteType.GT_AG,
            SpliceSiteType.GC_AG,
        )


@attrs.define(slots=True)
class IntronModel:
    """Represents an intron with its splice sites.

    Attributes:
        seqid: Scaffold/chromosome name.
        start: 0-based, first base of intron (donor).
        end: 0-based, last base + 1 (acceptor).
        strand: Strand (+ or -).
        parent_transcript: Parent transcript ID.
        donor_dinucleotide: Dinucleotide at donor site.
        acceptor_dinucleotide: Dinucleotide at acceptor site.
        rnaseq_support: Read count from junction extraction.
        pwm_score_donor: PWM score for donor site.
        pwm_score_acceptor: PWM score for acceptor site.
    """

    seqid: str
    start: int
    end: int
    strand: str
    parent_transcript: str

    donor_dinucleotide: str | None = None
    acceptor_dinucleotide: str | None = None

    # Evidence
    rnaseq_support: int = 0
    pwm_score_donor: float | None = None
    pwm_score_acceptor: float | None = None

    @property
    def length(self) -> int:
        """Intron length in base pairs."""
        return self.end - self.start

    @property
    def splice_type(self) -> SpliceSiteType:
        """Classify intron by donor-acceptor pair."""
        donor = (self.donor_dinucleotide or "").upper()
        acceptor = (self.acceptor_dinucleotide or "").upper()

        if donor == "GT" and acceptor == "AG":
            return SpliceSiteType.GT_AG
        elif donor == "GC" and acceptor == "AG":
            return SpliceSiteType.GC_AG
        elif donor == "AT" and acceptor == "AC":
            return SpliceSiteType.AT_AC
        else:
            return SpliceSiteType.OTHER

    @property
    def is_canonical(self) -> bool:
        """Check if intron has canonical splice sites."""
        return self.splice_type in (SpliceSiteType.GT_AG, SpliceSiteType.GC_AG)


@attrs.define(slots=True)
class SpliceCorrection:
    """Records a correction made to a splice site.

    Attributes:
        gene_id: Gene identifier.
        transcript_id: Transcript identifier.
        intron_index: Index of intron within transcript.
        site: "donor" or "acceptor".
        original_position: Original position (0-based).
        corrected_position: Corrected position (0-based).
        shift: Positive = downstream, negative = upstream.
        reason: Reason for correction.
        confidence: Confidence in correction (0-1).
        rnaseq_reads: Supporting read count (0 if PWM-only).
    """

    gene_id: str
    transcript_id: str
    intron_index: int
    site: Literal["donor", "acceptor"]
    original_position: int
    corrected_position: int
    shift: int
    reason: str
    confidence: float
    rnaseq_reads: int = 0

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "gene_id": self.gene_id,
            "transcript_id": self.transcript_id,
            "intron_index": self.intron_index,
            "site": self.site,
            "original_position": self.original_position,
            "corrected_position": self.corrected_position,
            "shift": self.shift,
            "reason": self.reason,
            "confidence": self.confidence,
            "rnaseq_reads": self.rnaseq_reads,
        }


@attrs.define(slots=True)
class GeneSpliceReport:
    """Splice refinement report for a gene.

    Attributes:
        gene_id: Gene identifier.
        n_introns: Total number of introns.
        n_supported: Introns with RNA-seq support.
        n_corrected: Number of corrected introns.
        n_canonical: Number of canonical introns.
        n_noncanonical: Number of non-canonical introns.
        corrections: List of corrections made.
        unsupported_introns: Indices of introns without evidence.
        flags: List of flags/warnings.
    """

    gene_id: str
    n_introns: int
    n_supported: int
    n_corrected: int
    n_canonical: int
    n_noncanonical: int
    corrections: list[SpliceCorrection] = attrs.Factory(list)
    unsupported_introns: list[int] = attrs.Factory(list)
    flags: list[str] = attrs.Factory(list)

    @property
    def support_ratio(self) -> float:
        """Fraction of introns with RNA-seq support."""
        return self.n_supported / self.n_introns if self.n_introns > 0 else 1.0

    @property
    def fully_supported(self) -> bool:
        """Check if all introns have RNA-seq support."""
        return self.n_supported == self.n_introns

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return {
            "gene_id": self.gene_id,
            "n_introns": self.n_introns,
            "n_supported": self.n_supported,
            "n_corrected": self.n_corrected,
            "n_canonical": self.n_canonical,
            "n_noncanonical": self.n_noncanonical,
            "support_ratio": self.support_ratio,
            "fully_supported": self.fully_supported,
            "corrections_summary": len(self.corrections),
            "unsupported_introns": ",".join(map(str, self.unsupported_introns)),
            "flags": ",".join(self.flags),
        }


# =============================================================================
# Position Weight Matrix
# =============================================================================


@attrs.define(slots=True)
class PositionWeightMatrix:
    """PWM for scoring splice site motifs.

    Donor site: -3 to +6 around GT (9 positions)
    Acceptor site: -14 to +1 around AG (16 positions)

    Attributes:
        matrix: Shape (n_positions, 4) for A, C, G, T.
        site_type: "donor" or "acceptor".
        offset: Position of splice site within window.
        background: Background nucleotide frequencies.
        positions: Position labels relative to splice site.
    """

    matrix: np.ndarray
    site_type: Literal["donor", "acceptor"]
    offset: int
    background: np.ndarray = attrs.Factory(lambda: DEFAULT_BACKGROUND.copy())
    positions: list[int] = attrs.Factory(list)

    @property
    def window_size(self) -> int:
        """Size of the motif window."""
        return len(self.matrix)

    def score(self, sequence: str) -> float:
        """Score a sequence against the PWM.

        Returns log-likelihood ratio vs background.

        Args:
            sequence: DNA sequence of length == window_size.

        Returns:
            Log-likelihood ratio score.

        Raises:
            ValueError: If sequence length doesn't match PWM.
        """
        sequence = sequence.upper()

        if len(sequence) != self.window_size:
            raise ValueError(
                f"Sequence length {len(sequence)} doesn't match "
                f"PWM window size {self.window_size}"
            )

        score = 0.0
        for i, base in enumerate(sequence):
            if base not in BASE_TO_INDEX:
                # Handle N or ambiguous bases - use background
                continue

            base_idx = BASE_TO_INDEX[base]
            prob = self.matrix[i, base_idx]
            bg = self.background[base_idx]

            # Log-likelihood ratio
            if prob > 0 and bg > 0:
                score += np.log2(prob / bg)

        return float(score)

    def score_all_positions(
        self,
        sequence: str,
        step: int = 1,
    ) -> np.ndarray:
        """Score all valid positions in a sequence.

        Args:
            sequence: DNA sequence longer than window_size.
            step: Step size for sliding window.

        Returns:
            Array of scores for each position.
        """
        sequence = sequence.upper()
        n_positions = (len(sequence) - self.window_size) // step + 1

        if n_positions <= 0:
            return np.array([])

        scores = np.zeros(n_positions)
        for i in range(n_positions):
            start = i * step
            end = start + self.window_size
            scores[i] = self.score(sequence[start:end])

        return scores

    @classmethod
    def from_sequences(
        cls,
        sequences: list[str],
        site_type: Literal["donor", "acceptor"],
        pseudocount: float = 0.1,
    ) -> PositionWeightMatrix:
        """Build PWM from aligned sequences.

        Args:
            sequences: List of aligned sequences (same length).
            site_type: "donor" or "acceptor".
            pseudocount: Pseudocount for smoothing (avoid log(0)).

        Returns:
            PositionWeightMatrix instance.

        Raises:
            ValueError: If sequences are empty or unequal length.
        """
        if not sequences:
            raise ValueError("Cannot build PWM from empty sequences")

        # Validate all sequences same length
        seq_len = len(sequences[0])
        if not all(len(s) == seq_len for s in sequences):
            raise ValueError("All sequences must have same length")

        # Count nucleotides at each position
        counts = np.zeros((seq_len, 4)) + pseudocount

        for seq in sequences:
            seq = seq.upper()
            for i, base in enumerate(seq):
                if base in BASE_TO_INDEX:
                    counts[i, BASE_TO_INDEX[base]] += 1

        # Normalize to frequencies
        matrix = counts / counts.sum(axis=1, keepdims=True)

        # Determine offset based on site type
        if site_type == "donor":
            # Standard: -3 to +6 around GT, GT at positions 3,4
            offset = 3
            positions = list(range(-3, 7))
        else:
            # Standard: -14 to +1 around AG, AG at positions 13,14
            offset = 13
            positions = list(range(-14, 2))

        # Adjust positions if sequence length differs from standard
        if len(positions) != seq_len:
            positions = list(range(seq_len))

        return cls(
            matrix=matrix,
            site_type=site_type,
            offset=offset,
            positions=positions,
        )

    @classmethod
    def load_plant_defaults(cls) -> tuple[PositionWeightMatrix, PositionWeightMatrix]:
        """Load pre-computed plant splice site PWMs.

        Derived from Arabidopsis + rice + maize canonical junctions.

        Returns:
            Tuple of (donor_pwm, acceptor_pwm).
        """
        # Try to load from package data
        try:
            # Python 3.9+ style
            data_path = resources.files("helixforge.data").joinpath("pwm_plant.json")
            with resources.as_file(data_path) as path:
                with open(path) as f:
                    data = json.load(f)
        except (FileNotFoundError, TypeError, AttributeError):
            # Fallback to embedded defaults
            logger.warning("Could not load PWM data file, using embedded defaults")
            data = _get_embedded_pwm_data()

        donor_pwm = cls.from_dict(data["donor"])
        acceptor_pwm = cls.from_dict(data["acceptor"])

        return donor_pwm, acceptor_pwm

    def to_dict(self) -> dict[str, Any]:
        """Serialize for storage."""
        return {
            "matrix": self.matrix.tolist(),
            "site_type": self.site_type,
            "offset": self.offset,
            "background": self.background.tolist(),
            "positions": self.positions,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PositionWeightMatrix:
        """Deserialize from dictionary."""
        return cls(
            matrix=np.array(data["matrix"], dtype=np.float64),
            site_type=data["site_type"],
            offset=data["offset"],
            background=np.array(data.get("background", DEFAULT_BACKGROUND)),
            positions=data.get("positions", []),
        )


def _get_embedded_pwm_data() -> dict[str, Any]:
    """Return embedded default PWM data for plant splice sites.

    Based on consensus from plant literature.
    """
    # Donor site: -3 to +6 (9 positions)
    # Consensus: M-A-G | G-T-R-A-G-T
    # M = A/C, R = A/G
    donor_matrix = [
        [0.35, 0.35, 0.15, 0.15],  # -3: M (A/C)
        [0.60, 0.10, 0.20, 0.10],  # -2: A strong
        [0.10, 0.05, 0.80, 0.05],  # -1: G
        [0.00, 0.00, 1.00, 0.00],  # +1: G (invariant)
        [0.00, 0.00, 0.00, 1.00],  # +2: T (invariant)
        [0.55, 0.05, 0.25, 0.15],  # +3: A/R
        [0.65, 0.05, 0.15, 0.15],  # +4: A
        [0.15, 0.05, 0.60, 0.20],  # +5: G
        [0.20, 0.10, 0.20, 0.50],  # +6: T
    ]

    # Acceptor site: -14 to +1 (16 positions)
    # Consensus: Y(12)-N-C-A-G | G
    # Y = pyrimidine tract (C/T rich)
    acceptor_matrix = [
        [0.15, 0.35, 0.15, 0.35],  # -14: Y
        [0.15, 0.35, 0.15, 0.35],  # -13: Y
        [0.15, 0.35, 0.15, 0.35],  # -12: Y
        [0.15, 0.35, 0.15, 0.35],  # -11: Y
        [0.15, 0.35, 0.15, 0.35],  # -10: Y
        [0.15, 0.35, 0.15, 0.35],  # -9: Y
        [0.15, 0.35, 0.15, 0.35],  # -8: Y
        [0.12, 0.38, 0.12, 0.38],  # -7: Y (stronger)
        [0.12, 0.38, 0.12, 0.38],  # -6: Y
        [0.12, 0.38, 0.12, 0.38],  # -5: Y
        [0.15, 0.40, 0.10, 0.35],  # -4: Y/C
        [0.25, 0.25, 0.25, 0.25],  # -3: N
        [0.15, 0.60, 0.10, 0.15],  # -2: C
        [0.70, 0.10, 0.10, 0.10],  # -1: A
        [0.00, 0.00, 1.00, 0.00],  # 0: G (invariant)
        [0.15, 0.15, 0.50, 0.20],  # +1: G
    ]

    return {
        "donor": {
            "matrix": donor_matrix,
            "site_type": "donor",
            "offset": 3,
            "background": [0.27, 0.23, 0.23, 0.27],
            "positions": [-3, -2, -1, 1, 2, 3, 4, 5, 6],
        },
        "acceptor": {
            "matrix": acceptor_matrix,
            "site_type": "acceptor",
            "offset": 14,
            "background": [0.27, 0.23, 0.23, 0.27],
            "positions": list(range(-14, 2)),
        },
    }


# =============================================================================
# Splice Refiner
# =============================================================================


class SpliceRefiner:
    """Refine splice sites using RNA-seq evidence and PWM scoring.

    Strategy:
    1. Match Helixer introns to RNA-seq junctions (exact or nearby)
    2. For unmatched introns, search for nearby junctions
    3. Score candidates with PWM
    4. Apply corrections within allowed window

    Attributes:
        genome: GenomeAccessor for sequence lookups.
        junctions: Pre-extracted RNA-seq junctions by scaffold.
        donor_pwm: PWM for donor site scoring.
        acceptor_pwm: PWM for acceptor site scoring.
        max_shift: Maximum correction distance.
        min_junction_reads: Minimum reads to trust a junction.
        pwm_threshold: Minimum PWM score for non-RNA-seq corrections.

    Example:
        >>> refiner = SpliceRefiner(genome, junctions)
        >>> refined_gene, report = refiner.refine_gene(gene)
    """

    def __init__(
        self,
        genome: GenomeAccessor,
        junctions: dict[str, list[SpliceJunction]],
        donor_pwm: PositionWeightMatrix | None = None,
        acceptor_pwm: PositionWeightMatrix | None = None,
        max_shift: int = DEFAULT_MAX_SHIFT,
        min_junction_reads: int = DEFAULT_MIN_JUNCTION_READS,
        pwm_threshold: float = DEFAULT_PWM_THRESHOLD,
    ) -> None:
        """Initialize refiner.

        Args:
            genome: GenomeAccessor for sequence lookups.
            junctions: Pre-extracted RNA-seq junctions by scaffold.
            donor_pwm: PWM for donor sites (uses plant default if None).
            acceptor_pwm: PWM for acceptor sites.
            max_shift: Maximum bases to shift a splice site.
            min_junction_reads: Minimum reads to trust a junction.
            pwm_threshold: Minimum PWM score for non-RNA-seq corrections.
        """
        self.genome = genome
        self.max_shift = max_shift
        self.min_junction_reads = min_junction_reads
        self.pwm_threshold = pwm_threshold

        # Load default PWMs if not provided
        if donor_pwm is None or acceptor_pwm is None:
            default_donor, default_acceptor = PositionWeightMatrix.load_plant_defaults()
            self.donor_pwm = donor_pwm or default_donor
            self.acceptor_pwm = acceptor_pwm or default_acceptor
        else:
            self.donor_pwm = donor_pwm
            self.acceptor_pwm = acceptor_pwm

        # Store raw junctions
        self._junctions = junctions

        # Build interval index for fast junction lookup
        self._junction_index = self._build_junction_index(junctions)

    def _build_junction_index(
        self,
        junctions: dict[str, list[SpliceJunction]],
    ) -> dict[str, dict[tuple[int, int], SpliceJunction]]:
        """Build index for fast junction lookup by coordinates.

        Returns dict mapping seqid -> {(start, end): junction}
        """
        index: dict[str, dict[tuple[int, int], SpliceJunction]] = {}

        for seqid, junction_list in junctions.items():
            index[seqid] = {}
            for j in junction_list:
                index[seqid][(j.start, j.end)] = j

        return index

    def refine_gene(
        self,
        gene: GeneModel,
    ) -> tuple[GeneModel, GeneSpliceReport]:
        """Refine all splice sites in a gene.

        Args:
            gene: GeneModel to refine.

        Returns:
            Tuple of (refined_gene, splice_report).
        """
        from copy import deepcopy

        # Work on a copy
        refined_gene = deepcopy(gene)

        all_corrections: list[SpliceCorrection] = []
        n_introns_total = 0
        n_supported = 0
        n_corrected = 0
        n_canonical = 0
        n_noncanonical = 0
        unsupported_indices: list[int] = []
        flags: list[str] = []

        # Process each transcript
        for transcript in refined_gene.transcripts:
            introns = self._get_introns(transcript)
            n_introns_total += len(introns)

            for idx, intron in enumerate(introns):
                # Refine this intron
                refined_intron, corrections = self._refine_intron(
                    intron,
                    gene.gene_id,
                    transcript.transcript_id,
                    idx,
                )

                # Track statistics
                if refined_intron.rnaseq_support >= self.min_junction_reads:
                    n_supported += 1
                else:
                    unsupported_indices.append(idx)

                if corrections:
                    n_corrected += 1
                    all_corrections.extend(corrections)

                if refined_intron.is_canonical:
                    n_canonical += 1
                else:
                    n_noncanonical += 1

                # Apply corrections to transcript if any
                if corrections:
                    self._apply_corrections_to_transcript(
                        transcript, idx, refined_intron
                    )

        # Update gene coordinates to match refined transcripts
        self._update_gene_bounds(refined_gene)

        # Generate flags
        if n_noncanonical > 0:
            flags.append(f"noncanonical_introns:{n_noncanonical}")
        if unsupported_indices:
            flags.append(f"unsupported_introns:{len(unsupported_indices)}")
        if n_corrected > 0:
            flags.append(f"corrected_introns:{n_corrected}")

        report = GeneSpliceReport(
            gene_id=gene.gene_id,
            n_introns=n_introns_total,
            n_supported=n_supported,
            n_corrected=n_corrected,
            n_canonical=n_canonical,
            n_noncanonical=n_noncanonical,
            corrections=all_corrections,
            unsupported_introns=unsupported_indices,
            flags=flags,
        )

        return refined_gene, report

    def refine_genes_parallel(
        self,
        genes: Iterable[GeneModel],
        n_workers: int = 1,
    ) -> Iterator[tuple[GeneModel, GeneSpliceReport]]:
        """Refine multiple genes in parallel.

        Args:
            genes: Iterable of GeneModel objects.
            n_workers: Number of parallel workers.

        Yields:
            Tuples of (refined_gene, splice_report).
        """
        genes_list = list(genes)

        if n_workers <= 1:
            for gene in genes_list:
                yield self.refine_gene(gene)
        else:
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(self.refine_gene, gene): gene
                    for gene in genes_list
                }

                for future in as_completed(futures):
                    try:
                        yield future.result()
                    except Exception as e:
                        gene = futures[future]
                        logger.error(f"Error refining gene {gene.gene_id}: {e}")
                        # Return unchanged gene with error report
                        report = GeneSpliceReport(
                            gene_id=gene.gene_id,
                            n_introns=0,
                            n_supported=0,
                            n_corrected=0,
                            n_canonical=0,
                            n_noncanonical=0,
                            flags=["refinement_error"],
                        )
                        yield gene, report

    def _get_introns(self, transcript: TranscriptModel) -> list[IntronModel]:
        """Extract introns from transcript exons."""
        introns = []

        sorted_exons = sorted(transcript.exons)
        for i in range(len(sorted_exons) - 1):
            intron_start = sorted_exons[i][1]  # End of current exon
            intron_end = sorted_exons[i + 1][0]  # Start of next exon

            if intron_end > intron_start:
                # Get dinucleotides
                donor_dinuc, acceptor_dinuc = self._get_dinucleotides(
                    transcript.seqid,
                    intron_start,
                    intron_end,
                    transcript.strand,
                )

                intron = IntronModel(
                    seqid=transcript.seqid,
                    start=intron_start,
                    end=intron_end,
                    strand=transcript.strand,
                    parent_transcript=transcript.transcript_id,
                    donor_dinucleotide=donor_dinuc,
                    acceptor_dinucleotide=acceptor_dinuc,
                )
                introns.append(intron)

        return introns

    def _refine_intron(
        self,
        intron: IntronModel,
        gene_id: str,
        transcript_id: str,
        intron_idx: int,
    ) -> tuple[IntronModel, list[SpliceCorrection]]:
        """Refine a single intron.

        Steps:
        1. Check for exact junction match
        2. If no match, search within max_shift window
        3. If multiple candidates, score with PWM
        4. Apply best correction if above threshold

        Args:
            intron: IntronModel to refine.
            gene_id: Parent gene ID.
            transcript_id: Parent transcript ID.
            intron_idx: Index of intron in transcript.

        Returns:
            Tuple of (refined_intron, list of corrections).
        """
        corrections: list[SpliceCorrection] = []

        # Try exact match first
        exact_junction = self._find_matching_junction(
            intron.seqid, intron.start, intron.end, intron.strand, tolerance=0
        )

        if exact_junction is not None:
            # Perfect match - just record support
            return IntronModel(
                seqid=intron.seqid,
                start=intron.start,
                end=intron.end,
                strand=intron.strand,
                parent_transcript=intron.parent_transcript,
                donor_dinucleotide=intron.donor_dinucleotide,
                acceptor_dinucleotide=intron.acceptor_dinucleotide,
                rnaseq_support=exact_junction.read_count,
                pwm_score_donor=self._score_donor_site(
                    intron.seqid, intron.start, intron.strand
                ),
                pwm_score_acceptor=self._score_acceptor_site(
                    intron.seqid, intron.end, intron.strand
                ),
            ), corrections

        # Search for nearby junctions
        nearby = self._find_nearby_junctions(
            intron.seqid,
            intron.start,
            intron.end,
            intron.strand,
            self.max_shift,
        )

        if nearby:
            # Score candidates and pick best
            best_junction = self._select_best_junction(
                intron, nearby
            )

            if best_junction is not None:
                # Create corrections
                new_start = best_junction.start
                new_end = best_junction.end

                if new_start != intron.start:
                    donor_shift = new_start - intron.start
                    corrections.append(
                        SpliceCorrection(
                            gene_id=gene_id,
                            transcript_id=transcript_id,
                            intron_index=intron_idx,
                            site="donor",
                            original_position=intron.start,
                            corrected_position=new_start,
                            shift=donor_shift,
                            reason="rnaseq_junction",
                            confidence=self._calculate_correction_confidence(
                                best_junction.read_count,
                                abs(donor_shift),
                            ),
                            rnaseq_reads=best_junction.read_count,
                        )
                    )

                if new_end != intron.end:
                    acceptor_shift = new_end - intron.end
                    corrections.append(
                        SpliceCorrection(
                            gene_id=gene_id,
                            transcript_id=transcript_id,
                            intron_index=intron_idx,
                            site="acceptor",
                            original_position=intron.end,
                            corrected_position=new_end,
                            shift=acceptor_shift,
                            reason="rnaseq_junction",
                            confidence=self._calculate_correction_confidence(
                                best_junction.read_count,
                                abs(acceptor_shift),
                            ),
                            rnaseq_reads=best_junction.read_count,
                        )
                    )

                # Get new dinucleotides
                donor_dinuc, acceptor_dinuc = self._get_dinucleotides(
                    intron.seqid, new_start, new_end, intron.strand
                )

                return IntronModel(
                    seqid=intron.seqid,
                    start=new_start,
                    end=new_end,
                    strand=intron.strand,
                    parent_transcript=intron.parent_transcript,
                    donor_dinucleotide=donor_dinuc,
                    acceptor_dinucleotide=acceptor_dinuc,
                    rnaseq_support=best_junction.read_count,
                    pwm_score_donor=self._score_donor_site(
                        intron.seqid, new_start, intron.strand
                    ),
                    pwm_score_acceptor=self._score_acceptor_site(
                        intron.seqid, new_end, intron.strand
                    ),
                ), corrections

        # No RNA-seq support - try to rescue non-canonical sites
        if not intron.is_canonical:
            rescued_intron, rescue_correction = self._rescue_noncanonical(
                intron, gene_id, transcript_id, intron_idx
            )
            if rescue_correction is not None:
                corrections.append(rescue_correction)
                return rescued_intron, corrections

        # Return original with PWM scores
        return IntronModel(
            seqid=intron.seqid,
            start=intron.start,
            end=intron.end,
            strand=intron.strand,
            parent_transcript=intron.parent_transcript,
            donor_dinucleotide=intron.donor_dinucleotide,
            acceptor_dinucleotide=intron.acceptor_dinucleotide,
            rnaseq_support=0,
            pwm_score_donor=self._score_donor_site(
                intron.seqid, intron.start, intron.strand
            ),
            pwm_score_acceptor=self._score_acceptor_site(
                intron.seqid, intron.end, intron.strand
            ),
        ), corrections

    def _find_matching_junction(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str,
        tolerance: int = 0,
    ) -> SpliceJunction | None:
        """Find RNA-seq junction matching intron coordinates.

        Args:
            seqid: Scaffold name.
            start: Intron start.
            end: Intron end.
            strand: Strand.
            tolerance: Allowed position difference.

        Returns:
            Matching SpliceJunction or None.
        """
        if seqid not in self._junction_index:
            return None

        junc_dict = self._junction_index[seqid]

        if tolerance == 0:
            return junc_dict.get((start, end))

        # Search within tolerance
        for junc_coords, junction in junc_dict.items():
            if (
                abs(junc_coords[0] - start) <= tolerance
                and abs(junc_coords[1] - end) <= tolerance
            ):
                # Check strand compatibility
                if junction.strand in (strand, "."):
                    return junction

        return None

    def _find_nearby_junctions(
        self,
        seqid: str,
        start: int,
        end: int,
        strand: str,
        window: int,
    ) -> list[SpliceJunction]:
        """Find all junctions within window of intron.

        Args:
            seqid: Scaffold name.
            start: Intron start.
            end: Intron end.
            strand: Strand.
            window: Search window size.

        Returns:
            List of nearby SpliceJunction objects.
        """
        if seqid not in self._junction_index:
            return []

        nearby = []
        junc_dict = self._junction_index[seqid]

        for (junc_start, junc_end), junction in junc_dict.items():
            # Check if within window
            start_diff = abs(junc_start - start)
            end_diff = abs(junc_end - end)

            if start_diff <= window and end_diff <= window:
                # Check strand compatibility
                if junction.strand in (strand, "."):
                    if junction.read_count >= self.min_junction_reads:
                        nearby.append(junction)

        return nearby

    def _select_best_junction(
        self,
        intron: IntronModel,
        candidates: list[SpliceJunction],
    ) -> SpliceJunction | None:
        """Select best junction from candidates.

        Scoring:
        1. Higher read count
        2. Lower distance to original
        3. Better PWM score

        Args:
            intron: Original intron.
            candidates: Candidate junctions.

        Returns:
            Best junction or None.
        """
        if not candidates:
            return None

        def score_junction(j: SpliceJunction) -> tuple[float, float, float]:
            # Distance penalty
            dist = abs(j.start - intron.start) + abs(j.end - intron.end)
            dist_score = -dist  # Lower distance is better

            # Read count (log scale)
            read_score = np.log1p(j.read_count)

            # PWM score
            pwm_donor = self._score_donor_site(j.seqid, j.start, intron.strand)
            pwm_acceptor = self._score_acceptor_site(j.seqid, j.end, intron.strand)
            pwm_score = (pwm_donor or 0) + (pwm_acceptor or 0)

            return (read_score, dist_score, pwm_score)

        # Sort by composite score
        scored = [(score_junction(j), j) for j in candidates]
        scored.sort(key=lambda x: x[0], reverse=True)

        return scored[0][1] if scored else None

    def _score_donor_site(
        self,
        seqid: str,
        position: int,
        strand: str,
    ) -> float | None:
        """Score a potential donor site with PWM.

        Args:
            seqid: Scaffold name.
            position: Position of first intronic base (0-based).
            strand: Strand.

        Returns:
            PWM score or None if sequence unavailable.
        """
        try:
            # Donor window: -3 to +6 around the GT
            # Position is first intronic base (G of GT)
            window_start = position - 3
            window_end = position + 6

            if window_start < 0:
                return None

            seq = self.genome.get_sequence(seqid, window_start, window_end, strand)

            if len(seq) < self.donor_pwm.window_size:
                return None

            return self.donor_pwm.score(seq[: self.donor_pwm.window_size])

        except (KeyError, ValueError):
            return None

    def _score_acceptor_site(
        self,
        seqid: str,
        position: int,
        strand: str,
    ) -> float | None:
        """Score a potential acceptor site with PWM.

        Args:
            seqid: Scaffold name.
            position: Position past last intronic base (0-based, exclusive).
            strand: Strand.

        Returns:
            PWM score or None if sequence unavailable.
        """
        try:
            # Acceptor window: -14 to +1 around AG
            # Position is first exonic base (after AG)
            # AG is at position-2, position-1
            window_start = position - 16
            window_end = position

            if window_start < 0:
                return None

            seq = self.genome.get_sequence(seqid, window_start, window_end, strand)

            if len(seq) < self.acceptor_pwm.window_size:
                return None

            return self.acceptor_pwm.score(seq[: self.acceptor_pwm.window_size])

        except (KeyError, ValueError):
            return None

    def _get_dinucleotides(
        self,
        seqid: str,
        intron_start: int,
        intron_end: int,
        strand: str,
    ) -> tuple[str, str]:
        """Get donor and acceptor dinucleotides for an intron.

        Args:
            seqid: Scaffold name.
            intron_start: First base of intron (0-based).
            intron_end: Position past last base of intron (0-based).
            strand: Strand.

        Returns:
            Tuple of (donor_dinuc, acceptor_dinuc).
        """
        try:
            # Donor: first 2 bases of intron
            donor = self.genome.get_sequence(
                seqid, intron_start, intron_start + 2, strand
            )

            # Acceptor: last 2 bases of intron
            acceptor = self.genome.get_sequence(
                seqid, intron_end - 2, intron_end, strand
            )

            return donor.upper(), acceptor.upper()

        except (KeyError, ValueError):
            return "NN", "NN"

    def _rescue_noncanonical(
        self,
        intron: IntronModel,
        gene_id: str,
        transcript_id: str,
        intron_idx: int,
    ) -> tuple[IntronModel, SpliceCorrection | None]:
        """Attempt to rescue non-canonical splice sites.

        Search nearby for canonical GT-AG sites with good PWM scores.
        Only apply if substantial improvement.

        Args:
            intron: Non-canonical intron.
            gene_id: Gene ID.
            transcript_id: Transcript ID.
            intron_idx: Intron index.

        Returns:
            Tuple of (intron, correction or None).
        """
        # Search for nearby canonical sites
        best_donor_pos = None
        best_donor_score = float("-inf")
        best_acceptor_pos = None
        best_acceptor_score = float("-inf")

        # Search window (smaller for rescue)
        rescue_window = min(5, self.max_shift // 2)

        # Search for donor sites
        for offset in range(-rescue_window, rescue_window + 1):
            pos = intron.start + offset
            try:
                dinuc = self.genome.get_sequence(
                    intron.seqid, pos, pos + 2, intron.strand
                ).upper()

                if dinuc in DONOR_CANONICAL:
                    score = self._score_donor_site(intron.seqid, pos, intron.strand)
                    if score is not None and score > best_donor_score:
                        best_donor_score = score
                        best_donor_pos = pos
            except (KeyError, ValueError):
                continue

        # Search for acceptor sites
        for offset in range(-rescue_window, rescue_window + 1):
            pos = intron.end + offset
            try:
                dinuc = self.genome.get_sequence(
                    intron.seqid, pos - 2, pos, intron.strand
                ).upper()

                if dinuc in ACCEPTOR_CANONICAL:
                    score = self._score_acceptor_site(intron.seqid, pos, intron.strand)
                    if score is not None and score > best_acceptor_score:
                        best_acceptor_score = score
                        best_acceptor_pos = pos
            except (KeyError, ValueError):
                continue

        # Only rescue if we found good canonical sites
        if best_donor_pos is not None and best_acceptor_pos is not None:
            if best_donor_score > self.pwm_threshold and best_acceptor_score > self.pwm_threshold:
                donor_shift = best_donor_pos - intron.start
                acceptor_shift = best_acceptor_pos - intron.end

                # Get new dinucleotides
                donor_dinuc, acceptor_dinuc = self._get_dinucleotides(
                    intron.seqid, best_donor_pos, best_acceptor_pos, intron.strand
                )

                new_intron = IntronModel(
                    seqid=intron.seqid,
                    start=best_donor_pos,
                    end=best_acceptor_pos,
                    strand=intron.strand,
                    parent_transcript=intron.parent_transcript,
                    donor_dinucleotide=donor_dinuc,
                    acceptor_dinucleotide=acceptor_dinuc,
                    rnaseq_support=0,
                    pwm_score_donor=best_donor_score,
                    pwm_score_acceptor=best_acceptor_score,
                )

                # Create single correction for both if shifted
                if donor_shift != 0 or acceptor_shift != 0:
                    correction = SpliceCorrection(
                        gene_id=gene_id,
                        transcript_id=transcript_id,
                        intron_index=intron_idx,
                        site="donor" if abs(donor_shift) >= abs(acceptor_shift) else "acceptor",
                        original_position=intron.start if donor_shift != 0 else intron.end,
                        corrected_position=best_donor_pos if donor_shift != 0 else best_acceptor_pos,
                        shift=donor_shift if donor_shift != 0 else acceptor_shift,
                        reason="canonical_rescue",
                        confidence=min(1.0, (best_donor_score + best_acceptor_score) / 20),
                        rnaseq_reads=0,
                    )
                    return new_intron, correction

        return intron, None

    def _calculate_correction_confidence(
        self,
        read_count: int,
        shift_distance: int,
    ) -> float:
        """Calculate confidence score for a correction.

        Higher read count and lower shift distance = higher confidence.

        Args:
            read_count: Supporting read count.
            shift_distance: Correction distance.

        Returns:
            Confidence score (0-1).
        """
        # Read count component (log scale, saturates around 100 reads)
        read_conf = min(1.0, np.log1p(read_count) / np.log1p(100))

        # Distance penalty (linear decay)
        dist_penalty = max(0, 1.0 - shift_distance / self.max_shift)

        # Combined score
        return read_conf * 0.7 + dist_penalty * 0.3

    def _apply_corrections_to_transcript(
        self,
        transcript: TranscriptModel,
        intron_idx: int,
        refined_intron: IntronModel,
    ) -> None:
        """Apply intron corrections to transcript exons and CDS.

        Modifies transcript in place.

        Args:
            transcript: Transcript to modify.
            intron_idx: Index of corrected intron.
            refined_intron: Refined intron with new coordinates.
        """
        sorted_exons = sorted(transcript.exons)

        if intron_idx >= len(sorted_exons) - 1:
            return

        # Update exon boundaries
        # Exon before intron: end changes to intron.start
        # Exon after intron: start changes to intron.end
        exon_before_idx = intron_idx
        exon_after_idx = intron_idx + 1

        # Update exons list
        new_exons = []
        for i, (start, end) in enumerate(sorted_exons):
            if i == exon_before_idx:
                new_end = refined_intron.start
                # Validate: don't create invalid exon (start >= end)
                if new_end <= start:
                    # Skip this correction - would create invalid exon
                    new_exons.append((start, end))
                else:
                    new_exons.append((start, new_end))
            elif i == exon_after_idx:
                new_start = refined_intron.end
                # Validate: don't create invalid exon (start >= end)
                if new_start >= end:
                    # Skip this correction - would create invalid exon
                    new_exons.append((start, end))
                else:
                    new_exons.append((new_start, end))
            else:
                new_exons.append((start, end))

        transcript.exons = new_exons

        # Update CDS coordinates based on coordinate overlap with affected exons
        # (not by index - CDS indices don't necessarily match exon indices due to UTRs)
        if transcript.cds:
            sorted_cds = sorted(transcript.cds)
            new_cds = []

            # Get the original exon boundaries (before correction)
            orig_exon_before = sorted_exons[exon_before_idx]
            orig_exon_after = sorted_exons[exon_after_idx]

            for cds_start, cds_end, phase in sorted_cds:
                new_cds_start = cds_start
                new_cds_end = cds_end

                # Check if this CDS overlaps with the exon BEFORE the intron
                # (its end might need to be adjusted to the new intron start)
                if (cds_start >= orig_exon_before[0] and cds_start < orig_exon_before[1] and
                        cds_end > orig_exon_before[0] and cds_end <= orig_exon_before[1] + 10):
                    # CDS is within or near the exon before the intron
                    # Adjust end if it extends past the new intron start
                    if cds_end > refined_intron.start:
                        new_cds_end = refined_intron.start

                # Check if this CDS overlaps with the exon AFTER the intron
                # (its start might need to be adjusted to the new intron end)
                if (cds_start >= orig_exon_after[0] - 10 and cds_start < orig_exon_after[1] and
                        cds_end > orig_exon_after[0] and cds_end <= orig_exon_after[1]):
                    # CDS is within or near the exon after the intron
                    # Adjust start if it's before the new intron end
                    if cds_start < refined_intron.end:
                        new_cds_start = refined_intron.end

                # Validate: don't create invalid CDS (start >= end)
                if new_cds_start >= new_cds_end:
                    # Keep original coordinates if correction would be invalid
                    new_cds.append((cds_start, cds_end, phase))
                else:
                    new_cds.append((new_cds_start, new_cds_end, phase))

            transcript.cds = new_cds

    def _update_gene_bounds(self, gene: GeneModel) -> None:
        """Update gene bounds to match refined transcripts.

        Args:
            gene: Gene to update.
        """
        if not gene.transcripts:
            return

        all_exons = []
        for tx in gene.transcripts:
            all_exons.extend(tx.exons)

        if all_exons:
            gene.start = min(e[0] for e in all_exons)
            gene.end = max(e[1] for e in all_exons)


# =============================================================================
# Report Writers
# =============================================================================


class SpliceReportWriter:
    """Write splice refinement reports."""

    @staticmethod
    def write_tsv(
        reports: Iterable[GeneSpliceReport],
        output_path: Path | str,
    ) -> None:
        """Write gene-level splice report.

        Columns:
        gene_id, n_introns, n_supported, n_corrected, n_canonical,
        n_noncanonical, support_ratio, corrections_summary, flags

        Args:
            reports: Iterable of GeneSpliceReport objects.
            output_path: Output file path.
        """
        output_path = Path(output_path)

        columns = [
            "gene_id",
            "n_introns",
            "n_supported",
            "n_corrected",
            "n_canonical",
            "n_noncanonical",
            "support_ratio",
            "fully_supported",
            "corrections_summary",
            "unsupported_introns",
            "flags",
        ]

        with open(output_path, "w") as f:
            f.write("\t".join(columns) + "\n")

            for report in reports:
                data = report.to_dict()
                row = [str(data.get(col, "")) for col in columns]
                f.write("\t".join(row) + "\n")

        logger.info(f"Wrote splice report to {output_path}")

    @staticmethod
    def write_corrections_detail(
        reports: Iterable[GeneSpliceReport],
        output_path: Path | str,
    ) -> None:
        """Write detailed per-correction report.

        Columns:
        gene_id, transcript_id, intron_index, site, original_pos,
        corrected_pos, shift, reason, confidence, rnaseq_reads

        Args:
            reports: Iterable of GeneSpliceReport objects.
            output_path: Output file path.
        """
        output_path = Path(output_path)

        columns = [
            "gene_id",
            "transcript_id",
            "intron_index",
            "site",
            "original_position",
            "corrected_position",
            "shift",
            "reason",
            "confidence",
            "rnaseq_reads",
        ]

        with open(output_path, "w") as f:
            f.write("\t".join(columns) + "\n")

            for report in reports:
                for correction in report.corrections:
                    data = correction.to_dict()
                    row = [str(data.get(col, "")) for col in columns]
                    f.write("\t".join(row) + "\n")

        logger.info(f"Wrote corrections detail to {output_path}")

    @staticmethod
    def write_unsupported_introns(
        reports: Iterable[GeneSpliceReport],
        genes: dict[str, GeneModel],
        output_path: Path | str,
    ) -> None:
        """Write BED of introns lacking RNA-seq support.

        Useful for identifying potential false positive introns.

        Args:
            reports: Iterable of GeneSpliceReport objects.
            genes: Dict mapping gene_id to GeneModel.
            output_path: Output file path.
        """
        output_path = Path(output_path)

        with open(output_path, "w") as f:
            f.write("#chrom\tstart\tend\tname\tscore\tstrand\n")

            for report in reports:
                if not report.unsupported_introns:
                    continue

                gene = genes.get(report.gene_id)
                if gene is None or not gene.transcripts:
                    continue

                tx = gene.transcripts[0]
                sorted_exons = sorted(tx.exons)

                for idx in report.unsupported_introns:
                    if idx < len(sorted_exons) - 1:
                        intron_start = sorted_exons[idx][1]
                        intron_end = sorted_exons[idx + 1][0]

                        name = f"{report.gene_id}_intron{idx}"
                        f.write(
                            f"{gene.seqid}\t{intron_start}\t{intron_end}\t"
                            f"{name}\t0\t{gene.strand}\n"
                        )

        logger.info(f"Wrote unsupported introns BED to {output_path}")

    @staticmethod
    def summary_statistics(
        reports: list[GeneSpliceReport],
    ) -> dict[str, Any]:
        """Compute summary statistics.

        Args:
            reports: List of GeneSpliceReport objects.

        Returns:
            Dict with summary statistics.
        """
        total_genes = len(reports)
        total_introns = sum(r.n_introns for r in reports)
        supported_introns = sum(r.n_supported for r in reports)
        corrected_introns = sum(r.n_corrected for r in reports)
        canonical_introns = sum(r.n_canonical for r in reports)
        noncanonical_introns = sum(r.n_noncanonical for r in reports)

        fully_supported_genes = sum(1 for r in reports if r.fully_supported)

        # Splice type distribution
        splice_types = {
            "GT-AG": canonical_introns,
            "non-canonical": noncanonical_introns,
        }

        return {
            "total_genes": total_genes,
            "total_introns": total_introns,
            "supported_introns": supported_introns,
            "support_rate": supported_introns / total_introns if total_introns > 0 else 0,
            "corrections_made": corrected_introns,
            "correction_rate": corrected_introns / total_introns if total_introns > 0 else 0,
            "canonical_introns": canonical_introns,
            "canonical_rate": canonical_introns / total_introns if total_introns > 0 else 0,
            "fully_supported_genes": fully_supported_genes,
            "fully_supported_rate": fully_supported_genes / total_genes if total_genes > 0 else 0,
            "splice_type_distribution": splice_types,
        }


# =============================================================================
# Utility Functions
# =============================================================================


def get_splice_dinucleotides(
    sequence: str,
    intron_start: int,
    intron_end: int,
) -> tuple[str, str]:
    """Extract donor and acceptor dinucleotides from sequence.

    Args:
        sequence: The genomic sequence (forward strand).
        intron_start: Start of intron (0-based, first intronic base).
        intron_end: End of intron (0-based, last intronic base + 1).

    Returns:
        Tuple of (donor_dinuc, acceptor_dinuc).
    """
    donor = sequence[intron_start : intron_start + 2]
    acceptor = sequence[intron_end - 2 : intron_end]
    return donor.upper(), acceptor.upper()


def is_canonical_splice_site(donor: str, acceptor: str) -> bool:
    """Check if a splice site pair is canonical.

    Args:
        donor: Donor dinucleotide.
        acceptor: Acceptor dinucleotide.

    Returns:
        True if the pair is canonical (GT-AG, GC-AG, or AT-AC).
    """
    donor = donor.upper()
    acceptor = acceptor.upper()

    # U2-dependent (major spliceosome)
    if donor in DONOR_CANONICAL and acceptor in ACCEPTOR_CANONICAL:
        return True

    # U12-dependent (minor spliceosome)
    if donor in DONOR_U12 and acceptor in ACCEPTOR_U12:
        return True

    return False


def classify_splice_type(donor: str, acceptor: str) -> SpliceSiteType:
    """Classify splice type based on dinucleotides.

    Args:
        donor: Donor dinucleotide.
        acceptor: Acceptor dinucleotide.

    Returns:
        SpliceSiteType classification.
    """
    donor = donor.upper()
    acceptor = acceptor.upper()

    if donor == "GT" and acceptor == "AG":
        return SpliceSiteType.GT_AG
    elif donor == "GC" and acceptor == "AG":
        return SpliceSiteType.GC_AG
    elif donor == "AT" and acceptor == "AC":
        return SpliceSiteType.AT_AC
    else:
        return SpliceSiteType.OTHER


# =============================================================================
# Legacy Compatibility
# =============================================================================

# Aliases for backward compatibility
SpliceSiteAnalyzer = SpliceRefiner

# Re-export for convenience
__all__ = [
    "SpliceSiteType",
    "SpliceSite",
    "IntronModel",
    "SpliceCorrection",
    "GeneSpliceReport",
    "PositionWeightMatrix",
    "SpliceRefiner",
    "SpliceReportWriter",
    "get_splice_dinucleotides",
    "is_canonical_splice_site",
    "classify_splice_type",
]
