"""Homology-based validation of gene predictions.

This module validates gene models by comparing predicted proteins
against reference protein databases and detecting issues like:
- Incomplete predictions (partial coverage)
- Orphan genes (no homology)
- Chimeric genes (fusions of multiple genes)
- Fragmented genes (single gene split into multiple predictions)
- TE overlap (genes overlapping transposable elements)

Example:
    >>> from helixforge.homology.validate import HomologyValidator, ValidationThresholds
    >>> validator = HomologyValidator(thresholds=ValidationThresholds.strict())
    >>> results = validator.validate_genes(gene_homology_map)
    >>> for gene_id, result in results.items():
    ...     print(f"{gene_id}: {result.status.value}")
"""

from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple

import attrs

if TYPE_CHECKING:
    from helixforge.homology.search import GeneHomology, HomologyHit
    from helixforge.io.gff import GeneModel

logger = logging.getLogger(__name__)


# =============================================================================
# Enums
# =============================================================================


class HomologyStatus(Enum):
    """Validation status based on homology evidence.

    Status values indicate the confidence and completeness
    of a gene prediction based on homology support.
    """

    COMPLETE = "complete"
    """Good coverage of a known protein (>=80% query coverage, >=50% subject coverage)."""

    PARTIAL = "partial"
    """Incomplete coverage (<80% query or <50% subject coverage)."""

    NO_HIT = "no_hit"
    """No significant homology found."""

    TE_OVERLAP = "te_overlap"
    """Gene overlaps transposable element annotation."""

    CHIMERIC = "chimeric"
    """Gene appears to be fusion of multiple unrelated genes."""

    FRAGMENTED = "fragmented"
    """Gene appears to be part of a fragmented prediction (fragment of larger gene)."""


# =============================================================================
# Configuration
# =============================================================================


@attrs.define(slots=True, frozen=True)
class ValidationThresholds:
    """Thresholds for homology-based validation.

    Attributes:
        min_identity: Minimum percent identity (0-100).
        min_query_coverage: Minimum query coverage (0-1).
        min_subject_coverage: Minimum subject coverage (0-1).
        max_evalue: Maximum e-value for significant hit.
        chimera_min_gap: Minimum gap between non-overlapping hits (aa).
        chimera_max_overlap: Maximum overlap for chimera detection (aa).
        fragment_max_overlap: Maximum overlap between fragment hits (aa).
        fragment_min_coverage: Minimum combined coverage for fragment group.
    """

    min_identity: float = 30.0
    min_query_coverage: float = 0.5
    min_subject_coverage: float = 0.5
    max_evalue: float = 1e-5
    chimera_min_gap: int = 20
    chimera_max_overlap: int = 10
    fragment_max_overlap: int = 50
    fragment_min_coverage: float = 0.8
    complete_query_coverage: float = 0.8
    complete_subject_coverage: float = 0.5

    @classmethod
    def default(cls) -> ValidationThresholds:
        """Create default thresholds for standard validation."""
        return cls()

    @classmethod
    def strict(cls) -> ValidationThresholds:
        """Create strict thresholds for high-confidence validation.

        Use for well-annotated reference databases like Swiss-Prot.
        """
        return cls(
            min_identity=50.0,
            min_query_coverage=0.7,
            min_subject_coverage=0.7,
            max_evalue=1e-10,
            complete_query_coverage=0.9,
            complete_subject_coverage=0.7,
        )

    @classmethod
    def relaxed(cls) -> ValidationThresholds:
        """Create relaxed thresholds for divergent species.

        Use for species without close relatives in databases.
        """
        return cls(
            min_identity=20.0,
            min_query_coverage=0.3,
            min_subject_coverage=0.3,
            max_evalue=1e-3,
            complete_query_coverage=0.6,
            complete_subject_coverage=0.3,
        )


# =============================================================================
# Data Structures
# =============================================================================


@attrs.define(slots=True)
class ValidationResult:
    """Result of homology validation for a gene.

    Attributes:
        gene_id: ID of the validated gene.
        status: Validation status (complete, partial, etc.).
        best_hit: Best homology hit (if any).
        n_hits: Total number of significant hits.
        query_coverage: Coverage of query by best hit.
        subject_coverage: Coverage of subject by best hit.
        identity: Percent identity of best hit.
        evalue: E-value of best hit.
        flags: Issue flags detected.
        chimeric_evidence: Evidence if chimera detected.
        fragment_group: Group ID if gene is part of fragment group.
        te_overlap_fraction: Fraction of gene overlapping TEs.
        notes: Additional notes about validation.
    """

    gene_id: str
    status: HomologyStatus
    best_hit: HomologyHit | None = None
    n_hits: int = 0
    query_coverage: float | None = None
    subject_coverage: float | None = None
    identity: float | None = None
    evalue: float | None = None
    flags: list[str] = attrs.Factory(list)
    chimeric_evidence: ChimericEvidence | None = None
    fragment_group: str | None = None
    te_overlap_fraction: float = 0.0
    notes: list[str] = attrs.Factory(list)

    @property
    def is_valid(self) -> bool:
        """Whether the gene passes validation."""
        return self.status in (HomologyStatus.COMPLETE, HomologyStatus.PARTIAL)

    @property
    def has_homology(self) -> bool:
        """Whether any homology was found."""
        return self.status != HomologyStatus.NO_HIT

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "gene_id": self.gene_id,
            "status": self.status.value,
            "best_hit_id": self.best_hit.subject_id if self.best_hit else None,
            "n_hits": self.n_hits,
            "query_coverage": self.query_coverage,
            "subject_coverage": self.subject_coverage,
            "identity": self.identity,
            "evalue": self.evalue,
            "flags": self.flags,
            "is_chimeric": self.chimeric_evidence is not None,
            "fragment_group": self.fragment_group,
            "te_overlap_fraction": self.te_overlap_fraction,
            "notes": self.notes,
        }


class ChimericEvidence(NamedTuple):
    """Evidence for a chimeric (fusion) gene.

    Chimeric genes show distinct regions with homology to
    different proteins, suggesting improper gene merging.

    Attributes:
        gene_id: ID of the potential chimeric gene.
        hit_a: First homology hit (N-terminal region).
        hit_b: Second homology hit (C-terminal region).
        gap: Gap between hits in query coordinates (aa).
        overlap: Overlap between hits (negative gap).
        confidence: Confidence in chimeric call (0-1).
    """

    gene_id: str
    hit_a: HomologyHit
    hit_b: HomologyHit
    gap: int
    overlap: int
    confidence: float

    @property
    def subject_a(self) -> str:
        """Subject ID of first hit."""
        return self.hit_a.subject_id

    @property
    def subject_b(self) -> str:
        """Subject ID of second hit."""
        return self.hit_b.subject_id

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "gene_id": self.gene_id,
            "subject_a": self.subject_a,
            "subject_b": self.subject_b,
            "gap": self.gap,
            "overlap": self.overlap,
            "confidence": self.confidence,
        }


@attrs.define(slots=True)
class FragmentGroup:
    """Group of genes that appear to be fragments of one gene.

    Fragment groups contain multiple predicted genes that all
    hit the same reference protein in non-overlapping regions,
    suggesting a single gene was incorrectly split.

    Attributes:
        group_id: Unique group identifier.
        subject_id: Reference protein ID they all hit.
        gene_ids: List of gene IDs in this group.
        subject_length: Length of the reference protein.
        total_coverage: Combined coverage of reference protein.
        hits: Mapping of gene_id to their best hit for this subject.
    """

    group_id: str
    subject_id: str
    gene_ids: list[str] = attrs.Factory(list)
    subject_length: int = 0
    total_coverage: float = 0.0
    hits: dict[str, HomologyHit] = attrs.Factory(dict)

    @property
    def n_fragments(self) -> int:
        """Number of fragments in this group."""
        return len(self.gene_ids)

    def add_gene(self, gene_id: str, hit: HomologyHit) -> None:
        """Add a gene to this fragment group."""
        if gene_id not in self.gene_ids:
            self.gene_ids.append(gene_id)
        self.hits[gene_id] = hit

    def calculate_coverage(self) -> float:
        """Calculate combined coverage of reference protein."""
        if not self.hits or self.subject_length == 0:
            return 0.0

        # Collect all covered regions
        covered_positions: set[int] = set()
        for hit in self.hits.values():
            for pos in range(hit.subject_start, hit.subject_end + 1):
                covered_positions.add(pos)

        self.total_coverage = len(covered_positions) / self.subject_length
        return self.total_coverage

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "group_id": self.group_id,
            "subject_id": self.subject_id,
            "gene_ids": self.gene_ids,
            "n_fragments": self.n_fragments,
            "subject_length": self.subject_length,
            "total_coverage": self.total_coverage,
        }


# =============================================================================
# TE Overlap Detection
# =============================================================================


@attrs.define(slots=True)
class TEInterval:
    """Transposable element interval.

    Attributes:
        seqid: Scaffold/chromosome name.
        start: Start position (0-based).
        end: End position (exclusive).
        te_family: TE family name (if known).
        te_class: TE class (if known).
    """

    seqid: str
    start: int
    end: int
    te_family: str = ""
    te_class: str = ""

    @property
    def length(self) -> int:
        """Length of the interval."""
        return self.end - self.start


def load_te_annotations(
    bed_path: Path | str,
) -> dict[str, list[TEInterval]]:
    """Load TE annotations from BED file.

    Args:
        bed_path: Path to BED file with TE annotations.
            Expected format: seqid, start, end, name (optional), score (optional), strand (optional)

    Returns:
        {seqid: [TEInterval, ...]} mapping.
    """
    bed_path = Path(bed_path)
    if not bed_path.exists():
        raise FileNotFoundError(f"TE annotation file not found: {bed_path}")

    te_by_seqid: dict[str, list[TEInterval]] = {}

    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                continue

            seqid = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # Optional fields
            te_family = fields[3] if len(fields) > 3 else ""
            te_class = ""

            # Parse TE class from name if format is "Class/Family"
            if "/" in te_family:
                te_class, te_family = te_family.split("/", 1)

            te = TEInterval(
                seqid=seqid,
                start=start,
                end=end,
                te_family=te_family,
                te_class=te_class,
            )

            if seqid not in te_by_seqid:
                te_by_seqid[seqid] = []
            te_by_seqid[seqid].append(te)

    # Sort intervals by start position
    for seqid in te_by_seqid:
        te_by_seqid[seqid].sort(key=lambda x: x.start)

    n_total = sum(len(v) for v in te_by_seqid.values())
    logger.info(f"Loaded {n_total} TE intervals from {bed_path}")

    return te_by_seqid


def calculate_te_overlap(
    gene_seqid: str,
    gene_start: int,
    gene_end: int,
    te_intervals: dict[str, list[TEInterval]],
) -> tuple[float, list[TEInterval]]:
    """Calculate fraction of gene overlapping TEs.

    Args:
        gene_seqid: Scaffold/chromosome of gene.
        gene_start: Gene start position (0-based).
        gene_end: Gene end position (exclusive).
        te_intervals: TE intervals by seqid.

    Returns:
        Tuple of (overlap_fraction, overlapping_tes).
    """
    if gene_seqid not in te_intervals:
        return 0.0, []

    gene_length = gene_end - gene_start
    if gene_length <= 0:
        return 0.0, []

    overlapping_tes = []
    overlapped_bases = 0

    for te in te_intervals[gene_seqid]:
        # Check for overlap
        if te.end <= gene_start or te.start >= gene_end:
            continue

        # Calculate overlap
        overlap_start = max(gene_start, te.start)
        overlap_end = min(gene_end, te.end)
        overlap_length = overlap_end - overlap_start

        if overlap_length > 0:
            overlapping_tes.append(te)
            overlapped_bases += overlap_length

    # Calculate fraction (could be > 1 if TEs overlap each other)
    overlap_fraction = min(1.0, overlapped_bases / gene_length)

    return overlap_fraction, overlapping_tes


# =============================================================================
# Validator Class
# =============================================================================


class HomologyValidator:
    """Validates gene predictions using protein homology.

    Analyzes homology search results to classify gene predictions
    and detect issues like chimeras, fragments, and TE overlaps.

    Example:
        >>> from helixforge.homology.validate import HomologyValidator
        >>> from helixforge.homology.search import HomologySearch, GeneHomology
        >>>
        >>> validator = HomologyValidator()
        >>> results = validator.validate_genes(gene_homology_map)
        >>> for gene_id, result in results.items():
        ...     print(f"{gene_id}: {result.status.value}")
    """

    def __init__(
        self,
        thresholds: ValidationThresholds | None = None,
        te_annotations: dict[str, list[TEInterval]] | None = None,
        te_overlap_threshold: float = 0.5,
    ) -> None:
        """Initialize validator.

        Args:
            thresholds: Validation thresholds. Uses default if None.
            te_annotations: TE intervals by seqid for overlap detection.
            te_overlap_threshold: Fraction of TE overlap to flag gene.
        """
        self.thresholds = thresholds or ValidationThresholds.default()
        self.te_intervals = te_annotations or {}
        self.te_overlap_threshold = te_overlap_threshold

    def validate_genes(
        self,
        gene_homology: dict[str, GeneHomology],
        gene_info: dict[str, tuple[str, int, int]] | None = None,
    ) -> dict[str, ValidationResult]:
        """Validate multiple genes.

        Args:
            gene_homology: {gene_id: GeneHomology} mapping.
            gene_info: Optional {gene_id: (seqid, start, end)} for TE overlap.

        Returns:
            {gene_id: ValidationResult} mapping.
        """
        results: dict[str, ValidationResult] = {}

        # First pass: validate individual genes
        for gene_id, homology in gene_homology.items():
            result = self._validate_single(gene_id, homology)

            # Check TE overlap if gene info provided
            if gene_info and gene_id in gene_info:
                seqid, start, end = gene_info[gene_id]
                overlap_frac, overlapping_tes = calculate_te_overlap(
                    seqid, start, end, self.te_intervals
                )
                result.te_overlap_fraction = overlap_frac

                if overlap_frac >= self.te_overlap_threshold:
                    result.status = HomologyStatus.TE_OVERLAP
                    result.flags.append("te_overlap")
                    te_names = [t.te_family or "unknown" for t in overlapping_tes[:3]]
                    result.notes.append(f"Overlaps TEs: {', '.join(te_names)}")

            results[gene_id] = result

        # Second pass: detect fragment groups
        fragment_groups = self._detect_fragment_groups(gene_homology)
        for group in fragment_groups:
            for gene_id in group.gene_ids:
                if gene_id in results:
                    results[gene_id].status = HomologyStatus.FRAGMENTED
                    results[gene_id].fragment_group = group.group_id
                    results[gene_id].flags.append("fragmented")
                    results[gene_id].notes.append(
                        f"Part of fragment group {group.group_id} "
                        f"({group.n_fragments} fragments)"
                    )

        logger.info(f"Validated {len(results)} genes")
        return results

    def _validate_single(
        self,
        gene_id: str,
        homology: GeneHomology,
    ) -> ValidationResult:
        """Validate a single gene.

        Args:
            gene_id: Gene identifier.
            homology: GeneHomology object with hits.

        Returns:
            ValidationResult for this gene.
        """
        # Filter hits by thresholds
        valid_hits = self._filter_hits(homology.hits)

        if not valid_hits:
            return ValidationResult(
                gene_id=gene_id,
                status=HomologyStatus.NO_HIT,
                n_hits=0,
                flags=["no_homology"],
                notes=["No significant homology found"],
            )

        # Get best hit (by e-value)
        best_hit = min(valid_hits, key=lambda h: h.evalue)

        # Determine status based on coverage
        query_cov = best_hit.query_coverage or 0.0
        subject_cov = best_hit.subject_coverage or 0.0

        if (
            query_cov >= self.thresholds.complete_query_coverage
            and subject_cov >= self.thresholds.complete_subject_coverage
        ):
            status = HomologyStatus.COMPLETE
        else:
            status = HomologyStatus.PARTIAL

        result = ValidationResult(
            gene_id=gene_id,
            status=status,
            best_hit=best_hit,
            n_hits=len(valid_hits),
            query_coverage=query_cov,
            subject_coverage=subject_cov,
            identity=best_hit.identity,
            evalue=best_hit.evalue,
        )

        # Check for chimera
        chimeric_evidence = self._detect_chimera(gene_id, valid_hits)
        if chimeric_evidence:
            result.status = HomologyStatus.CHIMERIC
            result.chimeric_evidence = chimeric_evidence
            result.flags.append("chimeric")
            result.notes.append(
                f"Chimeric: hits {chimeric_evidence.subject_a} and "
                f"{chimeric_evidence.subject_b}"
            )

        # Add coverage flags
        if query_cov < self.thresholds.min_query_coverage:
            result.flags.append("low_query_coverage")
            result.notes.append(f"Query coverage {query_cov:.1%} below threshold")

        if subject_cov < self.thresholds.min_subject_coverage:
            result.flags.append("low_subject_coverage")
            result.notes.append(f"Subject coverage {subject_cov:.1%} below threshold")

        return result

    def _filter_hits(self, hits: list[HomologyHit]) -> list[HomologyHit]:
        """Filter hits by thresholds.

        Args:
            hits: List of HomologyHit objects.

        Returns:
            Filtered list of significant hits.
        """
        valid = []
        for hit in hits:
            if hit.evalue > self.thresholds.max_evalue:
                continue
            if hit.identity < self.thresholds.min_identity:
                continue
            valid.append(hit)
        return valid

    def _detect_chimera(
        self,
        gene_id: str,
        hits: list[HomologyHit],
    ) -> ChimericEvidence | None:
        """Detect if gene is chimeric (fusion).

        A gene is considered chimeric if it has non-overlapping hits
        to different proteins in distinct regions of the query.

        Args:
            gene_id: Gene identifier.
            hits: Filtered homology hits.

        Returns:
            ChimericEvidence if detected, None otherwise.
        """
        if len(hits) < 2:
            return None

        # Group hits by unique subject
        hits_by_subject: dict[str, list[HomologyHit]] = {}
        for hit in hits:
            if hit.subject_id not in hits_by_subject:
                hits_by_subject[hit.subject_id] = []
            hits_by_subject[hit.subject_id].append(hit)

        # Need at least 2 different subjects
        if len(hits_by_subject) < 2:
            return None

        # Get best hit per subject
        best_per_subject = []
        for subject_id, subject_hits in hits_by_subject.items():
            best = min(subject_hits, key=lambda h: h.evalue)
            best_per_subject.append(best)

        # Sort by query start position
        best_per_subject.sort(key=lambda h: h.query_start)

        # Look for non-overlapping hits to different subjects
        for i in range(len(best_per_subject) - 1):
            hit_a = best_per_subject[i]
            hit_b = best_per_subject[i + 1]

            # Skip if same subject
            if hit_a.subject_id == hit_b.subject_id:
                continue

            # Calculate gap/overlap in query coordinates
            gap = hit_b.query_start - hit_a.query_end
            overlap = -gap if gap < 0 else 0

            # Check chimera criteria
            if gap >= self.thresholds.chimera_min_gap or (
                0 <= overlap <= self.thresholds.chimera_max_overlap
            ):
                # Calculate confidence based on hit quality
                avg_identity = (hit_a.identity + hit_b.identity) / 2
                confidence = min(1.0, avg_identity / 100.0)

                return ChimericEvidence(
                    gene_id=gene_id,
                    hit_a=hit_a,
                    hit_b=hit_b,
                    gap=max(0, gap),
                    overlap=overlap,
                    confidence=confidence,
                )

        return None

    def _detect_fragment_groups(
        self,
        gene_homology: dict[str, GeneHomology],
    ) -> list[FragmentGroup]:
        """Detect genes that appear to be fragments of the same gene.

        Fragments are detected when multiple genes have non-overlapping
        hits to the same reference protein, suggesting the reference
        gene was incorrectly split during prediction.

        Args:
            gene_homology: {gene_id: GeneHomology} mapping.

        Returns:
            List of FragmentGroup objects.
        """
        # Collect all hits grouped by subject
        hits_by_subject: dict[str, list[tuple[str, HomologyHit]]] = {}

        for gene_id, homology in gene_homology.items():
            valid_hits = self._filter_hits(homology.hits)
            for hit in valid_hits:
                if hit.subject_id not in hits_by_subject:
                    hits_by_subject[hit.subject_id] = []
                hits_by_subject[hit.subject_id].append((gene_id, hit))

        # Find subjects with hits from multiple genes
        fragment_groups = []
        group_counter = 0

        for subject_id, gene_hits in hits_by_subject.items():
            # Need at least 2 different genes
            unique_genes = set(gene_id for gene_id, _ in gene_hits)
            if len(unique_genes) < 2:
                continue

            # Get best hit per gene
            best_by_gene: dict[str, HomologyHit] = {}
            for gene_id, hit in gene_hits:
                if gene_id not in best_by_gene:
                    best_by_gene[gene_id] = hit
                elif hit.evalue < best_by_gene[gene_id].evalue:
                    best_by_gene[gene_id] = hit

            # Check if hits are non-overlapping in subject coordinates
            hits_sorted = sorted(
                best_by_gene.items(),
                key=lambda x: x[1].subject_start,
            )

            # Find genes with minimal overlap
            fragment_genes = []
            for i, (gene_id, hit) in enumerate(hits_sorted):
                is_fragment = True

                # Check overlap with next hit
                if i < len(hits_sorted) - 1:
                    next_gene, next_hit = hits_sorted[i + 1]
                    overlap = hit.subject_end - next_hit.subject_start
                    if overlap > self.thresholds.fragment_max_overlap:
                        # Too much overlap, probably same region
                        is_fragment = False

                if is_fragment:
                    fragment_genes.append((gene_id, hit))

            # Only create group if we have non-overlapping fragments
            if len(fragment_genes) >= 2:
                group_counter += 1
                group = FragmentGroup(
                    group_id=f"frag_{group_counter:04d}",
                    subject_id=subject_id,
                )

                # Get subject length from first hit
                first_hit = fragment_genes[0][1]
                if first_hit.subject_length:
                    group.subject_length = first_hit.subject_length

                for gene_id, hit in fragment_genes:
                    group.add_gene(gene_id, hit)

                # Calculate combined coverage
                group.calculate_coverage()

                # Only report if combined coverage is good
                if group.total_coverage >= self.thresholds.fragment_min_coverage:
                    fragment_groups.append(group)
                    logger.debug(
                        f"Fragment group {group.group_id}: {group.gene_ids} "
                        f"-> {subject_id} (coverage: {group.total_coverage:.1%})"
                    )

        logger.info(f"Detected {len(fragment_groups)} fragment groups")
        return fragment_groups

    def validate_from_search(
        self,
        hits_by_gene: dict[str, list[HomologyHit]],
        gene_info: dict[str, tuple[str, int, int]] | None = None,
    ) -> dict[str, ValidationResult]:
        """Validate genes directly from search results.

        Convenience method that wraps hits in GeneHomology objects.

        Args:
            hits_by_gene: {gene_id: [HomologyHit, ...]} from search.
            gene_info: Optional {gene_id: (seqid, start, end)} for TE overlap.

        Returns:
            {gene_id: ValidationResult} mapping.
        """
        from helixforge.homology.search import GeneHomology

        gene_homology = {
            gene_id: GeneHomology.from_hits(gene_id, hits)
            for gene_id, hits in hits_by_gene.items()
        }

        return self.validate_genes(gene_homology, gene_info)


# =============================================================================
# Summary Statistics
# =============================================================================


def summarize_validation(
    results: dict[str, ValidationResult],
) -> dict[str, int | float]:
    """Summarize validation results.

    Args:
        results: {gene_id: ValidationResult} mapping.

    Returns:
        Summary statistics dictionary.
    """
    total = len(results)
    if total == 0:
        return {
            "total": 0,
            "complete": 0,
            "partial": 0,
            "no_hit": 0,
            "chimeric": 0,
            "fragmented": 0,
            "te_overlap": 0,
            "pct_with_homology": 0.0,
            "pct_valid": 0.0,
        }

    counts = {status: 0 for status in HomologyStatus}
    for result in results.values():
        counts[result.status] += 1

    with_homology = total - counts[HomologyStatus.NO_HIT]
    valid = counts[HomologyStatus.COMPLETE] + counts[HomologyStatus.PARTIAL]

    return {
        "total": total,
        "complete": counts[HomologyStatus.COMPLETE],
        "partial": counts[HomologyStatus.PARTIAL],
        "no_hit": counts[HomologyStatus.NO_HIT],
        "chimeric": counts[HomologyStatus.CHIMERIC],
        "fragmented": counts[HomologyStatus.FRAGMENTED],
        "te_overlap": counts[HomologyStatus.TE_OVERLAP],
        "pct_with_homology": 100.0 * with_homology / total,
        "pct_valid": 100.0 * valid / total,
    }


def get_genes_by_status(
    results: dict[str, ValidationResult],
    status: HomologyStatus,
) -> list[str]:
    """Get gene IDs with a specific status.

    Args:
        results: {gene_id: ValidationResult} mapping.
        status: Status to filter by.

    Returns:
        List of gene IDs with that status.
    """
    return [
        gene_id
        for gene_id, result in results.items()
        if result.status == status
    ]


def get_flagged_genes(
    results: dict[str, ValidationResult],
    flag: str,
) -> list[str]:
    """Get gene IDs with a specific flag.

    Args:
        results: {gene_id: ValidationResult} mapping.
        flag: Flag to filter by.

    Returns:
        List of gene IDs with that flag.
    """
    return [
        gene_id
        for gene_id, result in results.items()
        if flag in result.flags
    ]
