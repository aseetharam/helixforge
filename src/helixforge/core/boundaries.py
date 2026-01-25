"""Gene boundary refinement.

This module provides tools for refining the boundaries of gene models,
including:

- Start codon optimization
- Stop codon verification
- UTR extension using RNA-seq data
- Phase consistency verification

Example:
    >>> from helixforge.core.boundaries import BoundaryAdjuster
    >>> from helixforge.io.fasta import GenomeAccessor
    >>> genome = GenomeAccessor("genome.fa")
    >>> adjuster = BoundaryAdjuster(genome)
    >>> adjusted_gene, info = adjuster.adjust_start_codon(gene)
"""

from __future__ import annotations

import logging
from copy import deepcopy
from typing import TYPE_CHECKING, Any, Literal, NamedTuple

if TYPE_CHECKING:
    from helixforge.io.bam import CoverageProfile
    from helixforge.io.fasta import GenomeAccessor
    from helixforge.io.gff import GeneModel, TranscriptModel

logger = logging.getLogger(__name__)

# =============================================================================
# Constants
# =============================================================================

# Start codons by frequency (standard genetic code)
START_CODONS = {
    "ATG": 1.0,  # Canonical start
    "CTG": 0.3,  # Alternative (Leu)
    "GTG": 0.2,  # Alternative (Val)
    "TTG": 0.1,  # Alternative (Leu)
}

# Stop codons (standard genetic code)
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Context preferences for start codon (Kozak-like)
# Position -3 to -1 before ATG, and +4 after ATG
KOZAK_PREFERENCE = {
    -3: {"A": 0.5, "G": 0.3, "C": 0.1, "T": 0.1},
    -2: {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},  # Neutral
    -1: {"C": 0.4, "A": 0.3, "G": 0.2, "T": 0.1},
    +4: {"G": 0.5, "A": 0.3, "C": 0.1, "T": 0.1},
}

# Codon tables (standard genetic code = 1)
CODON_TABLES = {
    1: {  # Standard
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
    },
}


# =============================================================================
# Data Structures
# =============================================================================


class StartCodonCandidate(NamedTuple):
    """A potential start codon for a transcript.

    Attributes:
        position: Genomic position of the start codon (0-based).
        codon: The start codon sequence.
        score: Combined score for this candidate.
        context_score: Kozak context score.
        coverage_support: Coverage evidence score.
        in_frame: Whether it's in frame with the existing CDS.
    """

    position: int
    codon: str
    score: float
    context_score: float
    coverage_support: float
    in_frame: bool


class BoundaryRefinement(NamedTuple):
    """Results of boundary refinement for a transcript.

    Attributes:
        transcript_id: ID of the refined transcript.
        original_start: Original CDS start.
        refined_start: Refined CDS start.
        original_end: Original CDS end.
        refined_end: Refined CDS end.
        utr5_extended: Whether 5' UTR was extended.
        utr3_extended: Whether 3' UTR was extended.
        changes: List of change descriptions.
    """

    transcript_id: str
    original_start: int
    refined_start: int
    original_end: int
    refined_end: int
    utr5_extended: bool
    utr3_extended: bool
    changes: list[str]


# =============================================================================
# BoundaryAdjuster Class
# =============================================================================


class BoundaryAdjuster:
    """Fine-tune exon boundaries after splice correction.

    Handles:
    - Start codon positioning
    - Stop codon positioning
    - UTR boundary refinement
    - Phase consistency

    Attributes:
        genome: GenomeAccessor for sequence lookups.
        codon_table: Genetic code table to use.

    Example:
        >>> adjuster = BoundaryAdjuster(genome)
        >>> adjusted, info = adjuster.adjust_start_codon(gene)
    """

    def __init__(
        self,
        genome: GenomeAccessor,
        codon_table: int = 1,
    ) -> None:
        """Initialize the boundary adjuster.

        Args:
            genome: GenomeAccessor for sequence lookups.
            codon_table: Genetic code table (default 1 = standard).
        """
        self.genome = genome
        self.codon_table = codon_table
        self._codons = CODON_TABLES.get(codon_table, CODON_TABLES[1])

    def adjust_start_codon(
        self,
        gene: GeneModel,
        search_window: int = 30,
    ) -> tuple[GeneModel, dict[str, Any]]:
        """Find optimal start codon position.

        Search upstream/downstream for ATG in frame.
        Prefer position with better Kozak-like context.

        Args:
            gene: GeneModel to adjust.
            search_window: Number of bases to search.

        Returns:
            Tuple of (adjusted_gene, adjustment_info).
        """
        adjusted_gene = deepcopy(gene)
        adjustments: list[dict[str, Any]] = []

        for transcript in adjusted_gene.transcripts:
            if not transcript.cds:
                continue

            # Get current CDS start
            sorted_cds = sorted(transcript.cds)
            if transcript.strand == "+":
                cds_start_pos = sorted_cds[0][0]
            else:
                cds_start_pos = sorted_cds[-1][1]

            # Find start codons in window
            candidates = self._find_start_codons(
                transcript.seqid,
                cds_start_pos - search_window if transcript.strand == "+" else cds_start_pos,
                cds_start_pos + search_window if transcript.strand == "+" else cds_start_pos + search_window,
                transcript.strand,
            )

            if not candidates:
                continue

            # Score candidates
            scored = []
            for pos, context in candidates:
                context_score = self._score_kozak_context(context)
                # Distance penalty (prefer positions closer to original)
                dist = abs(pos - cds_start_pos)
                dist_penalty = 1.0 - (dist / (search_window * 2))
                total_score = context_score * 0.6 + dist_penalty * 0.4
                scored.append((total_score, pos, context, context_score))

            # Pick best
            scored.sort(reverse=True)
            best_score, best_pos, best_context, best_context_score = scored[0]

            # Apply if different from current
            if best_pos != cds_start_pos:
                shift = best_pos - cds_start_pos
                self._apply_cds_start_shift(transcript, shift)

                adjustments.append({
                    "transcript_id": transcript.transcript_id,
                    "original_position": cds_start_pos,
                    "new_position": best_pos,
                    "shift": shift,
                    "codon": best_context[3:6],  # ATG part
                    "context_score": best_context_score,
                    "total_score": best_score,
                })

        # Update gene bounds
        self._update_gene_bounds(adjusted_gene)

        info = {
            "n_adjusted": len(adjustments),
            "adjustments": adjustments,
        }

        return adjusted_gene, info

    def adjust_stop_codon(
        self,
        gene: GeneModel,
        search_window: int = 30,
    ) -> tuple[GeneModel, dict[str, Any]]:
        """Find optimal stop codon position.

        Ensure CDS ends with valid stop (TAA, TAG, TGA).
        Extend or truncate as needed.

        Args:
            gene: GeneModel to adjust.
            search_window: Number of bases to search.

        Returns:
            Tuple of (adjusted_gene, adjustment_info).
        """
        adjusted_gene = deepcopy(gene)
        adjustments: list[dict[str, Any]] = []

        for transcript in adjusted_gene.transcripts:
            if not transcript.cds:
                continue

            # Get current CDS end
            sorted_cds = sorted(transcript.cds)
            if transcript.strand == "+":
                cds_end_pos = sorted_cds[-1][1]
            else:
                cds_end_pos = sorted_cds[0][0]

            # Check current stop codon
            current_has_stop, current_stop = self._check_stop_codon(
                transcript.seqid,
                cds_end_pos,
                transcript.strand,
            )

            if current_has_stop:
                continue  # Already has valid stop

            # Search for stop codons
            frame = self._get_cds_frame(transcript)
            stop_positions = self._find_stop_codons(
                transcript.seqid,
                cds_end_pos - search_window if transcript.strand == "-" else cds_end_pos,
                cds_end_pos + search_window if transcript.strand == "+" else cds_end_pos,
                transcript.strand,
                frame,
            )

            if not stop_positions:
                adjustments.append({
                    "transcript_id": transcript.transcript_id,
                    "original_position": cds_end_pos,
                    "new_position": cds_end_pos,
                    "shift": 0,
                    "codon": None,
                    "warning": "no_stop_found",
                })
                continue

            # Pick closest in-frame stop
            closest = min(stop_positions, key=lambda p: abs(p - cds_end_pos))
            shift = closest - cds_end_pos

            # Apply adjustment
            self._apply_cds_end_shift(transcript, shift)

            # Get the stop codon
            if transcript.strand == "+":
                stop_codon = self.genome.get_sequence(
                    transcript.seqid, closest, closest + 3
                )
            else:
                stop_codon = self.genome.get_sequence(
                    transcript.seqid, closest - 3, closest, "-"
                )

            adjustments.append({
                "transcript_id": transcript.transcript_id,
                "original_position": cds_end_pos,
                "new_position": closest,
                "shift": shift,
                "codon": stop_codon,
            })

        # Update gene bounds
        self._update_gene_bounds(adjusted_gene)

        info = {
            "n_adjusted": len([a for a in adjustments if a.get("shift", 0) != 0]),
            "adjustments": adjustments,
        }

        return adjusted_gene, info

    def verify_phase_consistency(
        self,
        gene: GeneModel,
    ) -> list[str]:
        """Check that CDS phases are consistent across exons.

        For each transcript:
        - Verify phase annotations match calculated values
        - Check that total CDS length is divisible by 3
        - Verify start and stop codons are present

        Args:
            gene: GeneModel to check.

        Returns:
            List of phase errors found.
        """
        errors: list[str] = []

        for transcript in gene.transcripts:
            if not transcript.cds:
                continue

            # Sort CDS by position
            sorted_cds = sorted(transcript.cds)
            if transcript.strand == "-":
                sorted_cds = sorted_cds[::-1]

            # Calculate expected phases
            cumulative_length = 0
            for i, (cds_start, cds_end, annotated_phase) in enumerate(sorted_cds):
                expected_phase = cumulative_length % 3
                if annotated_phase != expected_phase:
                    errors.append(
                        f"{transcript.transcript_id}: CDS {i} has phase {annotated_phase}, "
                        f"expected {expected_phase}"
                    )
                cds_length = cds_end - cds_start
                cumulative_length += cds_length

            # Check total length
            total_cds_length = sum(e - s for s, e, _ in sorted_cds)
            if total_cds_length % 3 != 0:
                errors.append(
                    f"{transcript.transcript_id}: CDS length {total_cds_length} not divisible by 3"
                )

            # Check start codon
            if transcript.strand == "+":
                first_cds = sorted_cds[0]
                start_seq = self.genome.get_sequence(
                    transcript.seqid, first_cds[0], first_cds[0] + 3
                )
            else:
                first_cds = sorted_cds[0]
                start_seq = self.genome.get_sequence(
                    transcript.seqid, first_cds[1] - 3, first_cds[1], "-"
                )

            if start_seq.upper() not in START_CODONS:
                errors.append(
                    f"{transcript.transcript_id}: No valid start codon "
                    f"(found {start_seq.upper()})"
                )

            # Check stop codon
            if transcript.strand == "+":
                last_cds = sorted_cds[-1]
                stop_seq = self.genome.get_sequence(
                    transcript.seqid, last_cds[1], last_cds[1] + 3
                )
            else:
                last_cds = sorted_cds[-1]
                stop_seq = self.genome.get_sequence(
                    transcript.seqid, last_cds[0] - 3, last_cds[0], "-"
                )

            if stop_seq.upper() not in STOP_CODONS:
                errors.append(
                    f"{transcript.transcript_id}: No valid stop codon "
                    f"(found {stop_seq.upper()})"
                )

        return errors

    def _find_start_codons(
        self,
        seqid: str,
        region_start: int,
        region_end: int,
        strand: str,
    ) -> list[tuple[int, str]]:
        """Find all ATG positions in region with context.

        Args:
            seqid: Scaffold name.
            region_start: Start of search region.
            region_end: End of search region.
            strand: Strand.

        Returns:
            List of (position, context) tuples.
            Context is the 7bp sequence around ATG (-3 to +3).
        """
        candidates = []

        # Ensure valid region
        region_start = max(0, region_start)
        try:
            region_end = min(region_end, self.genome.get_length(seqid))
        except KeyError:
            return candidates

        # Need 3 bp before and 4 bp after for context
        search_start = max(3, region_start)
        search_end = region_end - 4

        if search_end <= search_start:
            return candidates

        try:
            # Get sequence with context
            seq = self.genome.get_sequence(
                seqid, search_start - 3, search_end + 4, strand
            )
        except (KeyError, ValueError):
            return candidates

        # Scan for ATG
        seq_upper = seq.upper()
        for i in range(3, len(seq_upper) - 4):
            codon = seq_upper[i : i + 3]
            if codon in START_CODONS:
                # Position in genomic coordinates
                if strand == "+":
                    pos = search_start - 3 + i
                else:
                    # For minus strand, position is at the end
                    pos = search_end + 4 - i - 3

                # Context: 3 bp before, ATG, 1 bp after
                context = seq_upper[i - 3 : i + 4]
                candidates.append((pos, context))

        return candidates

    def _find_stop_codons(
        self,
        seqid: str,
        region_start: int,
        region_end: int,
        strand: str,
        frame: int,
    ) -> list[int]:
        """Find all in-frame stop codons in region.

        Args:
            seqid: Scaffold name.
            region_start: Start of search region.
            region_end: End of search region.
            strand: Strand.
            frame: Reading frame (0, 1, or 2).

        Returns:
            List of stop codon positions.
        """
        positions = []

        # Ensure valid region
        region_start = max(0, region_start)
        try:
            region_end = min(region_end, self.genome.get_length(seqid))
        except KeyError:
            return positions

        if region_end <= region_start + 3:
            return positions

        try:
            seq = self.genome.get_sequence(seqid, region_start, region_end, strand)
        except (KeyError, ValueError):
            return positions

        seq_upper = seq.upper()

        # Scan in frame
        start_offset = (3 - frame) % 3
        for i in range(start_offset, len(seq_upper) - 2, 3):
            codon = seq_upper[i : i + 3]
            if codon in STOP_CODONS:
                if strand == "+":
                    pos = region_start + i + 3  # Position after stop codon
                else:
                    pos = region_end - i  # Position before stop codon
                positions.append(pos)

        return positions

    def _score_kozak_context(self, context: str) -> float:
        """Score the Kozak-like context around a start codon.

        Args:
            context: 7bp sequence with ATG in positions 3-5.

        Returns:
            Context score (0.0 to 1.0).
        """
        if len(context) < 7:
            return 0.0

        score = 0.0
        n_positions = 0

        # Score positions -3, -1, +4 (relative to A of ATG)
        for offset, pos in [(-3, 0), (-1, 2), (+4, 6)]:
            if offset in KOZAK_PREFERENCE and pos < len(context):
                base = context[pos].upper()
                if base in KOZAK_PREFERENCE[offset]:
                    score += KOZAK_PREFERENCE[offset][base]
                    n_positions += 1

        # Bonus for canonical ATG
        codon = context[3:6].upper()
        if codon == "ATG":
            score += 0.5
            n_positions += 1

        return score / n_positions if n_positions > 0 else 0.0

    def _check_stop_codon(
        self,
        seqid: str,
        position: int,
        strand: str,
    ) -> tuple[bool, str | None]:
        """Check if position has a stop codon.

        Args:
            seqid: Scaffold name.
            position: CDS end position.
            strand: Strand.

        Returns:
            Tuple of (has_stop, stop_codon or None).
        """
        try:
            if strand == "+":
                codon = self.genome.get_sequence(seqid, position, position + 3)
            else:
                codon = self.genome.get_sequence(seqid, position - 3, position, "-")

            codon = codon.upper()
            return codon in STOP_CODONS, codon if codon in STOP_CODONS else None
        except (KeyError, ValueError):
            return False, None

    def _get_cds_frame(self, transcript: TranscriptModel) -> int:
        """Get the reading frame of the CDS.

        Args:
            transcript: Transcript to check.

        Returns:
            Frame (0, 1, or 2).
        """
        if not transcript.cds:
            return 0

        sorted_cds = sorted(transcript.cds)
        if transcript.strand == "-":
            sorted_cds = sorted_cds[::-1]

        # First CDS has phase annotation
        return sorted_cds[0][2] if len(sorted_cds[0]) > 2 else 0

    def _apply_cds_start_shift(
        self,
        transcript: TranscriptModel,
        shift: int,
    ) -> None:
        """Apply shift to CDS start position.

        Args:
            transcript: Transcript to modify.
            shift: Shift amount (positive = downstream).
        """
        if not transcript.cds:
            return

        sorted_cds = sorted(transcript.cds)

        if transcript.strand == "+":
            # Modify first CDS start
            old_start, old_end, phase = sorted_cds[0]
            new_cds = [(old_start + shift, old_end, 0)]  # Reset phase to 0
            new_cds.extend(sorted_cds[1:])
        else:
            # Modify last CDS end
            old_start, old_end, phase = sorted_cds[-1]
            new_cds = sorted_cds[:-1]
            new_cds.append((old_start, old_end - shift, 0))

        transcript.cds = new_cds

        # Also adjust first exon if needed
        sorted_exons = sorted(transcript.exons)
        if transcript.strand == "+":
            if sorted_exons[0][0] > sorted_cds[0][0] + shift:
                # Extend first exon
                new_exons = [(sorted_cds[0][0] + shift, sorted_exons[0][1])]
                new_exons.extend(sorted_exons[1:])
                transcript.exons = new_exons
        else:
            if sorted_exons[-1][1] < sorted_cds[-1][1] - shift:
                # Extend last exon
                new_exons = sorted_exons[:-1]
                new_exons.append((sorted_exons[-1][0], sorted_cds[-1][1] - shift))
                transcript.exons = new_exons

    def _apply_cds_end_shift(
        self,
        transcript: TranscriptModel,
        shift: int,
    ) -> None:
        """Apply shift to CDS end position.

        Args:
            transcript: Transcript to modify.
            shift: Shift amount (positive = downstream).
        """
        if not transcript.cds:
            return

        sorted_cds = sorted(transcript.cds)

        if transcript.strand == "+":
            # Modify last CDS end
            old_start, old_end, phase = sorted_cds[-1]
            new_cds = sorted_cds[:-1]
            new_cds.append((old_start, old_end + shift, phase))
        else:
            # Modify first CDS start
            old_start, old_end, phase = sorted_cds[0]
            new_cds = [(old_start - shift, old_end, phase)]
            new_cds.extend(sorted_cds[1:])

        transcript.cds = new_cds

        # Also adjust exon if needed
        sorted_exons = sorted(transcript.exons)
        if transcript.strand == "+":
            if sorted_exons[-1][1] < sorted_cds[-1][1] + shift:
                # Extend last exon
                new_exons = sorted_exons[:-1]
                new_exons.append((sorted_exons[-1][0], sorted_cds[-1][1] + shift))
                transcript.exons = new_exons
        else:
            if sorted_exons[0][0] > sorted_cds[0][0] - shift:
                # Extend first exon
                new_exons = [(sorted_cds[0][0] - shift, sorted_exons[0][1])]
                new_exons.extend(sorted_exons[1:])
                transcript.exons = new_exons

    def _update_gene_bounds(self, gene: GeneModel) -> None:
        """Update gene bounds to match transcripts.

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
# Legacy BoundaryRefiner (for backward compatibility)
# =============================================================================


class BoundaryRefiner:
    """Refines gene model boundaries using sequence and evidence.

    Optimizes start codons, verifies stop codons, and extends UTRs
    based on sequence context and RNA-seq evidence.

    Note: This class wraps BoundaryAdjuster for backward compatibility.
    New code should use BoundaryAdjuster directly.
    """

    def __init__(
        self,
        genome: GenomeAccessor,
        max_scan_distance: int = 300,
        min_utr_coverage: float = 2.0,
    ) -> None:
        """Initialize the boundary refiner.

        Args:
            genome: GenomeAccessor for sequence access.
            max_scan_distance: Maximum bp to scan for alternative starts.
            min_utr_coverage: Minimum coverage to extend UTR.
        """
        self._adjuster = BoundaryAdjuster(genome)
        self.genome = genome
        self.max_scan_distance = max_scan_distance
        self.min_utr_coverage = min_utr_coverage

    def refine_gene(
        self,
        gene: GeneModel,
        coverage: CoverageProfile | None = None,
    ) -> GeneModel:
        """Refine boundaries for all transcripts in a gene.

        Args:
            gene: The gene model to refine.
            coverage: Optional coverage profile for the region.

        Returns:
            Gene model with refined boundaries.
        """
        # Adjust start codons
        refined, _ = self._adjuster.adjust_start_codon(
            gene, search_window=self.max_scan_distance
        )

        # Adjust stop codons
        refined, _ = self._adjuster.adjust_stop_codon(
            refined, search_window=self.max_scan_distance
        )

        return refined

    def refine_transcript(
        self,
        transcript: TranscriptModel,
        coverage: CoverageProfile | None = None,
    ) -> tuple[TranscriptModel, BoundaryRefinement]:
        """Refine boundaries for a single transcript.

        Args:
            transcript: The transcript to refine.
            coverage: Optional coverage profile.

        Returns:
            Tuple of (refined transcript, refinement details).
        """
        from helixforge.io.gff import GeneModel

        # Create temporary gene wrapper
        temp_gene = GeneModel(
            gene_id="temp",
            seqid=transcript.seqid,
            start=transcript.start,
            end=transcript.end,
            strand=transcript.strand,
            transcripts=[transcript],
            source="helixer",
        )

        refined = self.refine_gene(temp_gene, coverage)
        refined_tx = refined.transcripts[0]

        # Get CDS bounds
        old_cds = sorted(transcript.cds) if transcript.cds else []
        new_cds = sorted(refined_tx.cds) if refined_tx.cds else []

        orig_start = old_cds[0][0] if old_cds else transcript.start
        orig_end = old_cds[-1][1] if old_cds else transcript.end
        new_start = new_cds[0][0] if new_cds else refined_tx.start
        new_end = new_cds[-1][1] if new_cds else refined_tx.end

        changes = []
        if new_start != orig_start:
            changes.append(f"start_shifted:{new_start - orig_start}")
        if new_end != orig_end:
            changes.append(f"end_shifted:{new_end - orig_end}")

        refinement = BoundaryRefinement(
            transcript_id=transcript.transcript_id,
            original_start=orig_start,
            refined_start=new_start,
            original_end=orig_end,
            refined_end=new_end,
            utr5_extended=False,  # TODO: implement UTR extension
            utr3_extended=False,
            changes=changes,
        )

        return refined_tx, refinement

    def find_start_codons(
        self,
        seqid: str,
        region_start: int,
        region_end: int,
        strand: str,
        current_start: int,
    ) -> list[StartCodonCandidate]:
        """Find potential start codons in a region."""
        raw_candidates = self._adjuster._find_start_codons(
            seqid, region_start, region_end, strand
        )

        candidates = []
        for pos, context in raw_candidates:
            context_score = self._adjuster._score_kozak_context(context)
            in_frame = (pos - current_start) % 3 == 0

            candidates.append(
                StartCodonCandidate(
                    position=pos,
                    codon=context[3:6],
                    score=context_score,
                    context_score=context_score,
                    coverage_support=0.0,
                    in_frame=in_frame,
                )
            )

        return sorted(candidates, key=lambda x: x.score, reverse=True)

    def verify_stop_codon(
        self,
        seqid: str,
        position: int,
        strand: str,
    ) -> bool:
        """Verify that a position contains a stop codon."""
        has_stop, _ = self._adjuster._check_stop_codon(seqid, position, strand)
        return has_stop

    def extend_utr(
        self,
        transcript: TranscriptModel,
        coverage: CoverageProfile,
        direction: Literal["5", "3"],
    ) -> int:
        """Extend UTR based on coverage evidence.

        TODO: Implement UTR extension logic.
        """
        if direction == "5":
            return transcript.start
        else:
            return transcript.end


# =============================================================================
# Utility Functions
# =============================================================================


def score_kozak_context(sequence: str, start_pos: int) -> float:
    """Score the Kozak-like context around a start codon.

    Args:
        sequence: Genomic sequence.
        start_pos: Position of the A in ATG (0-based).

    Returns:
        Context score (0.0 to 1.0).
    """
    if start_pos < 3 or start_pos + 4 > len(sequence):
        return 0.0

    context = sequence[start_pos - 3 : start_pos + 4]

    # Create temporary adjuster logic inline
    score = 0.0
    n_positions = 0

    for offset, pos in [(-3, 0), (-1, 2), (+4, 6)]:
        if offset in KOZAK_PREFERENCE and pos < len(context):
            base = context[pos].upper()
            if base in KOZAK_PREFERENCE[offset]:
                score += KOZAK_PREFERENCE[offset][base]
                n_positions += 1

    codon = context[3:6].upper()
    if codon == "ATG":
        score += 0.5
        n_positions += 1

    return score / n_positions if n_positions > 0 else 0.0


def find_in_frame_stops(
    sequence: str,
    start_pos: int,
    strand: str = "+",
) -> list[int]:
    """Find all in-frame stop codons downstream of a start.

    Args:
        sequence: Genomic sequence (forward strand).
        start_pos: Position of the start codon (0-based).
        strand: Strand (+ or -).

    Returns:
        List of stop codon positions.
    """
    stops = []
    seq_upper = sequence.upper()

    if strand == "+":
        for i in range(start_pos + 3, len(seq_upper) - 2, 3):
            codon = seq_upper[i : i + 3]
            if codon in STOP_CODONS:
                stops.append(i)
    else:
        for i in range(start_pos - 3, 2, -3):
            codon = seq_upper[i - 2 : i + 1]
            # Reverse complement
            rc_codon = "".join(
                {"A": "T", "T": "A", "C": "G", "G": "C"}.get(b, "N")
                for b in reversed(codon)
            )
            if rc_codon in STOP_CODONS:
                stops.append(i)

    return stops


# =============================================================================
# Exports
# =============================================================================

__all__ = [
    "BoundaryAdjuster",
    "BoundaryRefiner",
    "BoundaryRefinement",
    "StartCodonCandidate",
    "START_CODONS",
    "STOP_CODONS",
    "score_kozak_context",
    "find_in_frame_stops",
]
