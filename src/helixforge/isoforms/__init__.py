"""Alternative splicing and isoform reconstruction.

This module provides tools for analyzing alternative splicing and
reconstructing transcript isoforms:

- Collect splice junction evidence from RNA-seq
- Reconstruct alternative transcripts
- Select representative transcripts

Example:
    >>> from helixforge.isoforms import IsoformReconstructor
    >>> reconstructor = IsoformReconstructor(genome)
    >>> isoforms = reconstructor.reconstruct(gene, junctions)
"""

# TODO: Import and expose main classes once implemented
# from helixforge.isoforms.evidence import collect_junction_evidence
# from helixforge.isoforms.reconstruct import IsoformReconstructor
# from helixforge.isoforms.select import select_representative

__all__: list[str] = []
