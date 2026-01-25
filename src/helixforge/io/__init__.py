"""Input/output handlers for HelixForge.

This module provides readers and writers for various file formats
used in genome annotation:

- HDF5: Helixer prediction files
- GFF3/GTF: Gene annotation files
- FASTA: Genome sequence files
- BAM: RNA-seq alignment files

Example:
    >>> from helixforge.io import read_gff, read_fasta
    >>> genes = read_gff("annotations.gff3")
    >>> genome = read_fasta("genome.fa")
"""

# TODO: Import and expose main I/O functions once implemented
# from helixforge.io.gff import read_gff, write_gff
# from helixforge.io.fasta import FastaReader
# from helixforge.io.hdf5 import HelixerReader
# from helixforge.io.bam import BamReader

__all__: list[str] = [
    # "read_gff",
    # "write_gff",
    # "FastaReader",
    # "HelixerReader",
    # "BamReader",
]
