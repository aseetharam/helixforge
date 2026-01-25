"""Homology-based validation for HelixForge.

This module provides tools for validating gene predictions using
protein homology searches and domain annotations:

- BLAST/DIAMOND searches against protein databases
- Validation scoring based on homology evidence
- InterProScan integration for domain annotation

Example:
    >>> from helixforge.homology import HomologyValidator
    >>> validator = HomologyValidator(protein_db="uniprot.fa")
    >>> results = validator.validate(gene_models)
"""

# TODO: Import and expose main classes once implemented
# from helixforge.homology.search import run_diamond, run_blast
# from helixforge.homology.validate import HomologyValidator
# from helixforge.homology.domains import annotate_domains

__all__: list[str] = []
