"""Sequence manipulation utilities.

This module provides utilities for working with nucleotide and protein
sequences:

- Reverse complement
- Translation
- Codon tables
- Sequence statistics

Example:
    >>> from helixforge.utils.sequences import reverse_complement, translate
    >>> rc = reverse_complement("ATGCATGC")
    >>> protein = translate("ATGAAATAG")

TODO:
    - Implement translation with multiple genetic codes
    - Add sequence validation
    - Add ORF finding
    - Support for ambiguous bases
"""

from typing import Literal

# =============================================================================
# Constants
# =============================================================================

# Standard complement mapping
COMPLEMENT = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "n": "n",
}

# Standard genetic code (NCBI Table 1)
CODON_TABLE_STANDARD = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

# Start codons
START_CODONS = {"ATG", "CTG", "GTG", "TTG"}

# Stop codons
STOP_CODONS = {"TAA", "TAG", "TGA"}


# =============================================================================
# Complement and Reverse Complement
# =============================================================================


def complement(sequence: str) -> str:
    """Get the complement of a DNA sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        Complement sequence.
    """
    return "".join(COMPLEMENT.get(base, "N") for base in sequence)


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        Reverse complement sequence.
    """
    return complement(sequence)[::-1]


# =============================================================================
# Translation
# =============================================================================


def translate(
    sequence: str,
    table: int = 1,
    to_stop: bool = False,
) -> str:
    """Translate a DNA sequence to protein.

    Args:
        sequence: DNA coding sequence.
        table: NCBI genetic code table number.
        to_stop: If True, stop at first stop codon.

    Returns:
        Amino acid sequence.

    Raises:
        ValueError: If sequence length is not a multiple of 3.
    """
    if len(sequence) % 3 != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) is not a multiple of 3")

    # Get codon table
    if table != 1:
        raise NotImplementedError(f"Genetic code table {table} not yet implemented")

    codon_table = CODON_TABLE_STANDARD
    protein = []

    for i in range(0, len(sequence), 3):
        codon = sequence[i : i + 3].upper()
        aa = codon_table.get(codon, "X")

        if aa == "*" and to_stop:
            break

        protein.append(aa)

    return "".join(protein)


def find_orfs(
    sequence: str,
    min_length: int = 100,
    strand: Literal["+", "-", "both"] = "both",
) -> list[tuple[int, int, str, str]]:
    """Find open reading frames in a sequence.

    Args:
        sequence: DNA sequence.
        min_length: Minimum ORF length in nucleotides.
        strand: Which strand(s) to search.

    Returns:
        List of (start, end, strand, protein) tuples.
    """
    # TODO: Implement ORF finding
    raise NotImplementedError("find_orfs not yet implemented")


# =============================================================================
# Sequence Statistics
# =============================================================================


def gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence.

    Args:
        sequence: DNA sequence.

    Returns:
        GC content as a fraction (0.0 to 1.0).
    """
    if not sequence:
        return 0.0

    sequence = sequence.upper()
    gc = sum(1 for base in sequence if base in "GC")
    return gc / len(sequence)


def count_nucleotides(sequence: str) -> dict[str, int]:
    """Count nucleotides in a sequence.

    Args:
        sequence: DNA sequence.

    Returns:
        Dictionary mapping nucleotide to count.
    """
    counts = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
    for base in sequence.upper():
        if base in counts:
            counts[base] += 1
        else:
            counts["N"] += 1
    return counts


def is_valid_dna(sequence: str) -> bool:
    """Check if a sequence contains only valid DNA bases.

    Args:
        sequence: Sequence to check.

    Returns:
        True if sequence is valid DNA.
    """
    valid_bases = set("ATGCNatgcn")
    return all(base in valid_bases for base in sequence)


def is_valid_protein(sequence: str) -> bool:
    """Check if a sequence contains only valid amino acids.

    Args:
        sequence: Sequence to check.

    Returns:
        True if sequence is valid protein.
    """
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY*X")
    return all(aa in valid_aa for aa in sequence.upper())
