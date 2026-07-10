"""Sequence utilities: complement, translation, codon checks."""

from __future__ import annotations

from helixforge.constants import DEFAULT_TRANSL_TABLE

# Strict nucleotide alphabet (unambiguous bases).
_ACGT = frozenset("ACGT")

# Full IUPAC nucleotide complement map (uppercase). Keys are the complete IUPAC
# degenerate alphabet, so this doubles as the set of characters the tolerant
# sequence path accepts. ``N`` -> ``N``; degenerate codes complement to their
# base-complement set (e.g. ``R`` = A/G -> ``Y`` = T/C).
COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "R": "Y",
    "Y": "R",  # R=A/G  Y=C/T
    "S": "S",
    "W": "W",  # S=G/C  W=A/T  (self-complementary)
    "K": "M",
    "M": "K",  # K=G/T  M=A/C
    "B": "V",
    "V": "B",  # B=C/G/T  V=A/C/G
    "D": "H",
    "H": "D",  # D=A/G/T  H=A/C/T
}

# Set of characters the IUPAC-tolerant path accepts (the COMPLEMENT keys).
IUPAC_ALPHABET = frozenset(COMPLEMENT)

# Standard genetic code (NCBI translation table 1). 64 entries; '*' = stop.
CODON_TABLE = {
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

STANDARD_START_CODONS = ("ATG",)
STANDARD_STOP_CODONS = ("TAA", "TAG", "TGA")

# --- NCBI alternative translation tables ---
# Each table is expressed as the set of amino-acid reassignments relative to
# table 1 (the standard code); ``_build_table`` applies them. Only tables a
# plant pipeline plausibly needs for organellar contigs are provided, extend
# this map to add more. Start-codon variation is intentionally NOT encoded
# Start codons stay ATG by default; the start check is conservative.
# ``DEFAULT_TRANSL_TABLE`` is sourced from helixforge.constants (the single knob).
_TABLE_OVERRIDES: dict[int, dict[str, str]] = {
    2: {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},  # vertebrate mito
    4: {"TGA": "W"},  # mold/protozoan/coelenterate mito; Mycoplasma/Spiroplasma
    11: {},  # bacterial/archaeal/plastid: AA identical to table 1
}


def _build_table(table_id: int) -> dict[str, str]:
    tbl = dict(CODON_TABLE)
    tbl.update(_TABLE_OVERRIDES[table_id])
    return tbl


# Public registry: NCBI transl_table id -> codon dict. Table 1 is the canonical
# ``CODON_TABLE`` object (identity-preserved for backward compatibility).
CODON_TABLES: dict[int, dict[str, str]] = {1: CODON_TABLE}
CODON_TABLES.update({tid: _build_table(tid) for tid in _TABLE_OVERRIDES})


def get_codon_table(transl_table: int = DEFAULT_TRANSL_TABLE) -> dict[str, str]:
    """Return the codon→amino-acid dict for an NCBI ``transl_table`` id.

    Defaults to table 1 (standard code). Raises ``ValueError`` for an unknown id
    so a typo fails loudly rather than silently mis-translating.
    """
    try:
        return CODON_TABLES[transl_table]
    except KeyError:
        raise ValueError(
            f"unsupported NCBI transl_table id {transl_table!r}; "
            f"known tables: {sorted(CODON_TABLES)}"
        ) from None


def stop_codons_for_table(transl_table: int = DEFAULT_TRANSL_TABLE) -> tuple[str, ...]:
    """Return the stop codons (codons mapping to ``'*'``) for ``transl_table``."""
    tbl = get_codon_table(transl_table)
    return tuple(c for c, aa in tbl.items() if aa == "*")


# ---------------------------------------------------------------------------
# Cleaning / validation
# ---------------------------------------------------------------------------


def _clean(seq: str) -> str:
    """Uppercase and validate against the IUPAC alphabet (the default path).

    Tolerates ``N``/IUPAC degenerate codes (common in draft plant assemblies);
    raises ``ValueError`` only on a character outside the IUPAC nucleotide
    alphabet entirely (e.g. ``'Z'``, ``'1'``).
    """
    s = seq.upper()
    for ch in s:
        if ch not in IUPAC_ALPHABET:
            raise ValueError(f"non-IUPAC character {ch!r} in sequence")
    return s


def require_acgt(seq: str) -> str:
    """Strict path: uppercase and reject any non-ACGT character.

    Retained for callers that genuinely require unambiguous bases. The default
    sequence helpers use the IUPAC-tolerant :func:`_clean` instead.
    """
    s = seq.upper()
    for ch in s:
        if ch not in _ACGT:
            raise ValueError(f"non-ACGT character {ch!r} in sequence")
    return s


def has_ambiguous_base(seq: str) -> bool:
    """True if ``seq`` contains any base outside strict ACGT (after upper-casing)."""
    return any(ch not in _ACGT for ch in seq.upper())


def reverse_complement(seq: str) -> str:
    """Reverse complement of an IUPAC nucleotide sequence.

    Complements the full IUPAC alphabet (``N`` -> ``N``, ``R`` -> ``Y``, …);
    raises ``ValueError`` only on a non-IUPAC character.
    """
    s = _clean(seq)
    return "".join(COMPLEMENT[ch] for ch in reversed(s))


def translate(
    seq: str, phase: int = 0, transl_table: int = DEFAULT_TRANSL_TABLE
) -> str:
    """Translate ``seq`` to a protein string starting at offset ``phase``.

    Trailing 1–2 nt that do not form a full codon are ignored. Stop codons
    translate to ``'*'``. Any codon containing an ``N``/IUPAC ambiguous base
    (or otherwise absent from the chosen table) translates to ``'X'``, never an
    exception. ``transl_table`` selects the NCBI genetic code.
    """
    if phase not in (0, 1, 2):
        raise ValueError(f"phase must be 0, 1 or 2, got {phase}")
    table = get_codon_table(transl_table)
    s = _clean(seq)
    protein = []
    for i in range(phase, len(s) - 2, 3):
        protein.append(table.get(s[i : i + 3], "X"))
    return "".join(protein)


def is_start_codon(
    codon: str, start_codons: tuple[str, ...] = STANDARD_START_CODONS
) -> bool:
    """True if ``codon`` is a start codon (default: ATG).

    An ambiguous (``N``/IUPAC) codon is *not* a start codon (returns ``False``,
    never raises); callers that need to distinguish "ambiguous" from "absent"
    should test :func:`has_ambiguous_base` first.
    """
    return _clean(codon) in start_codons


def is_stop_codon(codon: str, transl_table: int = DEFAULT_TRANSL_TABLE) -> bool:
    """True if ``codon`` is a stop codon under ``transl_table``.

    An ambiguous codon translates to ``'X'`` and is therefore *not* a stop
    (returns ``False``, never raises).
    """
    return get_codon_table(transl_table).get(_clean(codon)) == "*"


def check_internal_stops(
    seq: str, phase: int = 0, transl_table: int = DEFAULT_TRANSL_TABLE
) -> list[int]:
    """Return nucleotide positions of stop codons that are *not* the final codon.

    A premature (internal) stop indicates a broken ORF.
    Positions are 0-based indices into ``seq`` of the stop codon's first base.
    An ``N``/IUPAC codon translates to ``'X'`` and is treated as "unknown," not a
    premature stop, never raises on ambiguity.
    """
    table = get_codon_table(transl_table)
    s = _clean(seq)
    codon_starts = list(range(phase, len(s) - 2, 3))
    last = codon_starts[-1] if codon_starts else None
    positions = []
    for i in codon_starts:
        if table.get(s[i : i + 3], "X") == "*" and i != last:
            positions.append(i)
    return positions


def extract_start_codon(genome_seq: str, cds_start: int, strand: str) -> str:
    """Extract the start codon (3 nt, coding direction) from a contig sequence.

    ``genome_seq`` is the plus-strand contig sequence. For ``+`` strand,
    ``cds_start`` is the low genomic coordinate of the CDS 5' end and the codon
    is ``genome_seq[cds_start:cds_start+3]``. For ``-`` strand, ``cds_start`` is
    the *high* genomic coordinate (the 5' end in coding direction) and the codon
    is ``reverse_complement(genome_seq[cds_start-3:cds_start])``. Tolerates
    ``N``/IUPAC bases.
    """
    _validate_strand_arg(strand)
    if strand == "+":
        return _clean(genome_seq[cds_start : cds_start + 3])
    return reverse_complement(genome_seq[cds_start - 3 : cds_start])


def extract_stop_codon_after_cds(genome_seq: str, cds_end: int, strand: str) -> str:
    """Extract the stop codon immediately following the CDS 3' end.

    ``genome_seq`` is the plus-strand contig sequence. For ``+`` strand,
    ``cds_end`` is the half-open high coordinate of the CDS and the stop codon is
    ``genome_seq[cds_end:cds_end+3]``. For ``-`` strand, ``cds_end`` is the low
    genomic coordinate of the CDS and the stop (downstream = lower coords) is
    ``reverse_complement(genome_seq[cds_end-3:cds_end])``. Tolerates ``N``/IUPAC.
    """
    _validate_strand_arg(strand)
    if strand == "+":
        return _clean(genome_seq[cds_end : cds_end + 3])
    return reverse_complement(genome_seq[cds_end - 3 : cds_end])


def _validate_strand_arg(strand: str) -> None:
    if strand not in ("+", "-"):
        raise ValueError(f"strand must be '+' or '-', got {strand!r}")
