"""Phase 25 D2 — property-based codon / sequence invariants (M10).

Hypothesis round-trip + symmetry properties for ``utils/sequences.py`` — the
highest-risk, strand-sensitive code in the project (CLAUDE.md §5/§6; the v1
catastrophe was a minus-strand codon-handling error). Every property that can
involve strand is tested on **both** strands with concrete seeds (``@example``)
plus generated cases.

Coding-direction convention (CLAUDE.md §5): plus strand 5' = lowest-coordinate
CDS start; minus strand 5' = highest-coordinate CDS end (reverse-complemented).
"""

from __future__ import annotations

from hypothesis import assume, example, given
from hypothesis import strategies as st

from helixforge.utils.sequences import (
    check_internal_stops,
    extract_start_codon,
    extract_stop_codon_after_cds,
    reverse_complement,
    translate,
)

_dna = st.text(alphabet="ACGT", min_size=0, max_size=300)
_dna_nonempty = st.text(alphabet="ACGT", min_size=1, max_size=300)
# Lengths divisible by 3 keep "no trailing partial codon" tidy where it matters.
_dna_codon_aligned = st.integers(min_value=0, max_value=100).map(
    lambda n: n * 3
).flatmap(lambda L: st.text(alphabet="ACGT", min_size=L, max_size=L))


@given(_dna)
@example("")
@example("ATGAAACCCGGGTTTTAA")          # the conftest plus-strand ORF
@example("A")                            # odd length
@example("AC")                           # length-2 edge
def test_revcomp_is_an_involution(seq):
    assert reverse_complement(reverse_complement(seq)) == seq


@given(_dna)
@example("ATGAAACCCGGGTTTTAA")
def test_revcomp_preserves_length(seq):
    assert len(reverse_complement(seq)) == len(seq)


@given(_dna, st.integers(min_value=0, max_value=2))
@example("ATGAAACCCGGGTTTTAA", 0)
@example("ATGAAACCCGGGTTTTAA", 1)        # phase shifts the start offset
def test_translate_length_parity(seq, phase):
    """``translate`` emits exactly the number of *complete* codons at ``phase``;
    trailing 1-2 nt are dropped (CLAUDE.md §5.3 / docstring)."""
    protein = translate(seq, phase)
    expected = max(0, (len(seq) - phase)) // 3
    assert len(protein) == expected
    # Every residue is a single-char amino acid or a stop marker.
    assert all(c == "*" or c.isalpha() for c in protein)


@given(_dna_codon_aligned)
@example("ATGAAACCCGGG")                 # no stop at all
def test_translate_codon_aligned_length(seq):
    assert len(translate(seq, 0)) == len(seq) // 3


# --- start/stop extraction symmetry, BOTH strands ----------------------------

# A genome long enough to slice a codon on either side of any interior position.
_genome = st.text(alphabet="ACGT", min_size=6, max_size=200)


@given(_genome, st.data())
@example("ATGAAACCCGGGTTTTAA" + "A" * 42, None)
def test_start_codon_plus_is_slice_minus_is_revcomp(genome, data):
    """``extract_start_codon`` reads the plus-strand slice on ``+`` and the
    reverse-complement of the symmetric slice on ``-`` (CLAUDE.md §5 coding
    direction). Verified at a concrete and a generated interior position."""
    n = len(genome)
    # plus: codon at [pos, pos+3); needs pos+3 <= n
    pos_plus = data.draw(st.integers(min_value=0, max_value=n - 3)) if data else 0
    assert extract_start_codon(genome, pos_plus, "+") == genome[pos_plus:pos_plus + 3].upper()
    # minus: cds_start is the HIGH coord; codon = revcomp(genome[pos-3:pos]); needs pos >= 3
    pos_minus = data.draw(st.integers(min_value=3, max_value=n)) if data else 3
    assert extract_start_codon(genome, pos_minus, "-") == reverse_complement(
        genome[pos_minus - 3:pos_minus]
    )


@given(_genome, st.data())
@example("ATGAAACCCGGGTTTTAA" + "A" * 42, None)
def test_stop_codon_after_cds_both_strands(genome, data):
    """``extract_stop_codon_after_cds``: on ``+`` the stop is the slice just past
    the half-open CDS high coord; on ``-`` it is revcomp of the slice just below
    the CDS low coord (downstream = lower coords)."""
    n = len(genome)
    cds_end_plus = data.draw(st.integers(min_value=0, max_value=n - 3)) if data else 0
    assert extract_stop_codon_after_cds(genome, cds_end_plus, "+") == (
        genome[cds_end_plus:cds_end_plus + 3].upper()
    )
    cds_end_minus = data.draw(st.integers(min_value=3, max_value=n)) if data else 3
    assert extract_stop_codon_after_cds(genome, cds_end_minus, "-") == reverse_complement(
        genome[cds_end_minus - 3:cds_end_minus]
    )


@given(_genome, st.integers(min_value=3))
def test_minus_start_equals_plus_start_on_revcomp_genome(genome, pos):
    """Coding-direction consistency: the minus-strand start codon read at high
    coord ``pos`` equals the plus-strand start codon read on the reverse-
    complemented genome at the mirror position. This is the exact symmetry whose
    violation destroyed 96% of v1 genes (CLAUDE.md §6)."""
    n = len(genome)
    assume(0 < pos <= n)
    rc = reverse_complement(genome)
    minus_codon = extract_start_codon(genome, pos, "-")
    # Mirror position on the revcomp'd genome: the low coord of the same codon.
    plus_codon = extract_start_codon(rc, n - pos, "+")
    assert minus_codon == plus_codon


@given(_dna)
@example("ATGAAACCCGGGTTTTAA")          # only a terminal stop → no internal stops
@example("ATGTAAAAATTT")                 # TAA at index 3 IS internal
def test_check_internal_stops_consistent_with_translation(seq):
    """A position is flagged internal-stop iff its codon translates to ``*`` and
    it is not the final codon. Symmetric with ``translate``."""
    positions = check_internal_stops(seq, 0)
    codon_starts = list(range(0, max(0, len(seq) - 2), 3))
    last = codon_starts[-1] if codon_starts else None
    expected = [
        i for i in codon_starts
        if translate(seq[i:i + 3], 0) == "*" and i != last
    ]
    assert positions == expected
