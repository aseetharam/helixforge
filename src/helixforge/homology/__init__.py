"""Homology-based validation for HelixForge.

This module provides tools for validating gene predictions using
protein homology searches and domain annotations:

- Diamond/MMseqs2 searches against protein databases
- Validation scoring based on homology evidence
- Chimera and fragment detection
- TE overlap analysis
- Database download and management utilities

Example:
    >>> from helixforge.homology import HomologySearch, HomologyValidator
    >>> from helixforge.homology import SearchTool, HomologyStatus
    >>>
    >>> # Run homology search
    >>> searcher = HomologySearch(tool=SearchTool.DIAMOND, database="swissprot.dmnd")
    >>> hits = searcher.search_and_parse("proteins.fa")
    >>>
    >>> # Validate gene predictions
    >>> validator = HomologyValidator()
    >>> results = validator.validate_from_search(hits)
    >>> for gene_id, result in results.items():
    ...     print(f"{gene_id}: {result.status.value}")
"""

from helixforge.homology.search import (
    GeneHomology,
    HomologyHit,
    HomologySearch,
    SearchTool,
    extract_proteins_from_gff,
    get_sequence_lengths,
    translate_sequence,
)
from helixforge.homology.validate import (
    ChimericEvidence,
    FragmentGroup,
    HomologyStatus,
    HomologyValidator,
    TEInterval,
    ValidationResult,
    ValidationThresholds,
    get_flagged_genes,
    get_genes_by_status,
    load_te_annotations,
    summarize_validation,
)
from helixforge.homology.databases import (
    DatabaseInfo,
    DatabaseManager,
    DatabaseType,
    download_swissprot,
    download_swissprot_plants,
    download_uniref,
    list_databases,
    setup_database,
)

__all__ = [
    # Search
    "HomologySearch",
    "HomologyHit",
    "GeneHomology",
    "SearchTool",
    "extract_proteins_from_gff",
    "get_sequence_lengths",
    "translate_sequence",
    # Validation
    "HomologyValidator",
    "HomologyStatus",
    "ValidationResult",
    "ValidationThresholds",
    "ChimericEvidence",
    "FragmentGroup",
    "TEInterval",
    "load_te_annotations",
    "summarize_validation",
    "get_genes_by_status",
    "get_flagged_genes",
    # Database management
    "DatabaseManager",
    "DatabaseInfo",
    "DatabaseType",
    "download_swissprot",
    "download_swissprot_plants",
    "download_uniref",
    "setup_database",
    "list_databases",
]
