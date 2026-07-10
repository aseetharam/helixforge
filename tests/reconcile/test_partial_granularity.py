"""Phase 18 D2 — 5'/3' partial-ORF granularity.

A single ``cds_partial`` bool used to make a 5'-partial-but-3'-complete CDS skip
the **stop** check too. With ``cds_partial_5prime`` / ``cds_partial_3prime``
split out, the start check skips only on 5'-partial and the stop check skips only
on 3'-partial — the gate gets *more* precise, never looser. ``cds_partial``
remains a derived property (= 5p OR 3p).

Concrete literal coordinates only (CLAUDE.md §12); both strands mandatory.
"""

import attrs
import pytest

from helixforge.reconcile.models import CDSSegment, Exon, TranscriptCandidate
from helixforge.reconcile.validate import check_start_codon, check_stop_codon
from helixforge.utils.sequences import reverse_complement

CONTIG_LEN = 3000


class MockGenome:
    def __init__(self, sequences):
        self.sequences = dict(sequences)

    def get_sequence(self, seqid, start, end, strand="+"):
        seq = self.sequences[seqid][start:end]
        return reverse_complement(seq) if strand == "-" else seq


def genome_with(planted):
    arr = ["A"] * CONTIG_LEN
    for pos, seq in planted.items():
        for i, ch in enumerate(seq):
            arr[pos + i] = ch
    return MockGenome({"chr1": "".join(arr)})


def make_tx(exons, strand, cds, p5=False, p3=False):
    return TranscriptCandidate(
        transcript_id="HFG_00001.1",
        locus_id="HFG_00001",
        source="mikado",
        seqid="chr1",
        start=exons[0][0],
        end=exons[-1][1],
        strand=strand,
        exons=[Exon(s, e) for s, e in exons],
        cds=[CDSSegment(*c) for c in cds],
        cds_partial_5prime=p5,
        cds_partial_3prime=p3,
        is_primary=True,
    )


# plus: exon (1000,1100), CDS [1000,1099); start [1000,1003), stop [1096,1099)
def PLUS(p5=False, p3=False):
    return make_tx([(1000, 1100)], "+", [(1000, 1099, 0)], p5=p5, p3=p3)


# minus: exon (2000,2100), CDS [2001,2100); start [2097,2100), stop [2001,2004)
def MINUS(p5=False, p3=False):
    return make_tx([(2000, 2100)], "-", [(2001, 2100, 0)], p5=p5, p3=p3)


# --- the property identity (= 5p OR 3p) ---

def test_cds_partial_property_equals_5p_or_3p():
    assert PLUS(p5=False, p3=False).cds_partial is False
    assert PLUS(p5=True, p3=False).cds_partial is True
    assert PLUS(p5=False, p3=True).cds_partial is True
    assert PLUS(p5=True, p3=True).cds_partial is True


def test_cds_partial_compat_alias_sets_both():
    t = TranscriptCandidate(
        transcript_id="x.1", locus_id="x", source="mikado", seqid="chr1",
        start=1000, end=1100, strand="+", exons=[Exon(1000, 1100)],
        cds=[CDSSegment(1000, 1099, 0)], cds_partial=True,
    )
    assert t.cds_partial_5prime is True
    assert t.cds_partial_3prime is True
    assert t.cds_partial is True


def test_evolve_of_one_flag_not_clobbered_by_compat():
    # The compat alias must not resurrect on attrs.evolve of a granular flag.
    t = PLUS(p5=True, p3=False)
    t2 = attrs.evolve(t, cds_partial_3prime=True)
    assert t2.cds_partial_5prime is True
    assert t2.cds_partial_3prime is True


# --- 5'-partial-only: start skipped, STOP STILL CHECKED (the D2 fix) ---

def test_5prime_partial_plus_skips_start_checks_stop():
    g = genome_with({1000: "CCC", 1096: "CCC"})  # bad start AND bad stop
    tx = PLUS(p5=True, p3=False)
    assert check_start_codon(tx, g) is None              # skipped (5'-partial)
    assert check_stop_codon(tx, g).name == "NO_STOP"     # STILL verified


def test_5prime_partial_minus_skips_start_checks_stop():
    g = genome_with({2097: "AAA", 2001: "GGG"})  # RC(AAA)=TTT: bad start; RC(GGG)=CCC: bad stop
    tx = MINUS(p5=True, p3=False)
    assert check_start_codon(tx, g) is None
    assert check_stop_codon(tx, g).name == "NO_STOP"


def test_5prime_partial_plus_good_stop_passes():
    g = genome_with({1000: "CCC", 1096: "TAA"})  # good stop present
    tx = PLUS(p5=True, p3=False)
    assert check_stop_codon(tx, g) is None


# --- 3'-partial-only: stop skipped, START still checked ---

def test_3prime_partial_plus_skips_stop_checks_start():
    g = genome_with({1000: "CCC", 1096: "CCC"})
    tx = PLUS(p5=False, p3=True)
    assert check_stop_codon(tx, g) is None               # skipped (3'-partial)
    assert check_start_codon(tx, g).name == "NO_START"   # STILL verified


def test_3prime_partial_minus_skips_stop_checks_start():
    g = genome_with({2097: "AAA", 2001: "GGG"})
    tx = MINUS(p5=False, p3=True)
    assert check_stop_codon(tx, g) is None
    assert check_start_codon(tx, g).name == "NO_START"


def test_3prime_partial_plus_good_start_passes():
    g = genome_with({1000: "ATG", 1096: "CCC"})  # good start present
    tx = PLUS(p5=False, p3=True)
    assert check_start_codon(tx, g) is None


# --- both-partial: both skipped + mod-3 exempt ---

def test_both_partial_plus_skips_both():
    g = genome_with({1000: "CCC", 1096: "CCC"})
    tx = PLUS(p5=True, p3=True)
    assert check_start_codon(tx, g) is None
    assert check_stop_codon(tx, g) is None


def test_both_partial_minus_skips_both():
    g = genome_with({2097: "AAA", 2001: "GGG"})
    tx = MINUS(p5=True, p3=True)
    assert check_start_codon(tx, g) is None
    assert check_stop_codon(tx, g) is None


def test_partial_either_end_exempts_mod3():
    # 50 bp CDS (not mod-3) is accepted when either end is partial, rejected when
    # complete — exemption keys on the OR, both strands.
    make_tx([(1000, 1100)], "+", [(1000, 1050, 0)], p5=True, p3=False)
    make_tx([(2000, 2100)], "-", [(2050, 2100, 0)], p5=False, p3=True)
    with pytest.raises(ValueError):
        make_tx([(1000, 1100)], "+", [(1000, 1050, 0)], p5=False, p3=False)
