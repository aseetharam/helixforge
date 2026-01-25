"""HelixForge: A modular toolkit for refining Helixer gene predictions.

HelixForge transforms raw Helixer predictions into publication-quality
genome annotations using RNA-seq evidence, homology validation,
confidence scoring, and isoform reconstruction.

Example:
    >>> import helixforge
    >>> helixforge.__version__
    '0.1.0-alpha'

Modules:
    io: Input/output handlers for HDF5, GFF, FASTA, and BAM files
    core: Core refinement logic (confidence, splice, boundaries)
    homology: Homology-based validation
    isoforms: Alternative splicing and isoform reconstruction
    qc: Quality control flags, reports, and filters
    viz: Visualization tools
    parallel: Parallelization utilities
    utils: General utilities
"""

__version__ = "0.1.0-alpha"
__author__ = "Arun Seetharam"

# TODO: Import and expose main classes once implemented
# from helixforge.core.models import GeneModel, Transcript, Exon, GenomicRegion
# from helixforge.core.confidence import ConfidenceScore
# from helixforge.io.gff import read_gff, write_gff
# from helixforge.io.fasta import FastaReader

__all__ = [
    "__version__",
    "__author__",
]
