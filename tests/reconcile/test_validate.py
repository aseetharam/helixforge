"""Phase 7 D3 — structural codon gate (read-only).

Concrete literal coordinates only (CLAUDE.md §12); both strands mandatory.
Codons are planted at exact genomic positions in a synthetic ``MockGenome``.
"""

import pytest

from helixforge.reconcile.models import (
    CDSSegment,
    Exon,
    LocusClassification,
    ReconciledGene,
    TranscriptCandidate,
)
from helixforge.reconcile.validate import (
    check_cds_mod3,
    check_cds_within_exons,
    check_internal_stops,
    check_long_introns,
    check_short_cds,
    check_short_exons,
    check_start_codon,
    check_stop_codon,
    extract_cds_sequence,
    validate_all,
    validate_gene,
)
from helixforge.utils.sequences import reverse_complement

CONTIG_LEN = 3000


class MockGenome:
    """Tiny in-memory genome; reverse-complements on the minus strand."""

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


def make_tx(exons, strand, cds=None, partial=False, tid="HFG_00001.1"):
    return TranscriptCandidate(
        transcript_id=tid,
        locus_id="HFG_00001",
        source="mikado",
        seqid="chr1",
        start=exons[0][0],
        end=exons[-1][1],
        strand=strand,
        exons=[Exon(s, e) for s, e in exons],
        cds=[CDSSegment(*c) for c in cds] if cds else None,
        cds_partial=partial,
        is_primary=True,
    )


def make_gene(transcript, origin="mikado_1to1", tier=2, transcripts=None, primary=None):
    txs = transcripts or [transcript]
    return ReconciledGene(
        gene_id="HFG_00001",
        seqid="chr1",
        start=min(t.start for t in txs),
        end=max(t.end for t in txs),
        strand=transcript.strand,
        tier=tier,
        transcripts=txs,
        primary_transcript_id=primary or transcript.transcript_id,
        classification=LocusClassification("HFG_00001", "EXPRESSED"),
        origin=origin,
    )


def names(flags):
    return {f.name for f in flags}


# A plus-strand single-CDS transcript: exon (1000,1100), CDS [1000,1099).
PLUS_TX = lambda partial=False, cds=((1000, 1099, 0),): make_tx(  # noqa: E731
    [(1000, 1100)], "+", cds=cds, partial=partial
)
# A minus-strand single-CDS transcript: exon (2000,2100), CDS [2001,2100).
MINUS_TX = lambda partial=False, cds=((2001, 2100, 0),): make_tx(  # noqa: E731
    [(2000, 2100)], "-", cds=cds, partial=partial
)


# ---------------------------------------------------------------------------
# start codon — both strands
# ---------------------------------------------------------------------------

def test_start_codon_plus_good():
    g = genome_with({1000: "ATG", 1096: "TAA"})
    assert check_start_codon(PLUS_TX(), g) is None


def test_start_codon_plus_bad():
    g = genome_with({1000: "CCC", 1096: "TAA"})
    assert check_start_codon(PLUS_TX(), g).name == "NO_START"


def test_start_codon_minus_good():
    g = genome_with({2097: "CAT", 2001: "TTA"})  # RC(CAT) = ATG
    assert check_start_codon(MINUS_TX(), g) is None


def test_start_codon_minus_bad():
    g = genome_with({2097: "AAA", 2001: "TTA"})  # RC(AAA) = TTT
    assert check_start_codon(MINUS_TX(), g).name == "NO_START"


def test_start_codon_skipped_without_genome():
    assert check_start_codon(PLUS_TX(), None) is None


def test_start_codon_skipped_for_partial():
    g = genome_with({1000: "CCC"})  # bad start, but partial
    assert check_start_codon(PLUS_TX(partial=True), g) is None


# ---------------------------------------------------------------------------
# stop codon — both strands (stop-inclusive: last 3 CDS bases are the stop)
# ---------------------------------------------------------------------------

def test_stop_codon_plus_good():
    # CDS [1000,1099): last 3 bases = [1096,1099)
    g = genome_with({1000: "ATG", 1096: "TAA"})
    assert check_stop_codon(PLUS_TX(), g) is None


def test_stop_codon_plus_bad():
    g = genome_with({1000: "ATG", 1096: "CCC"})
    assert check_stop_codon(PLUS_TX(), g).name == "NO_STOP"


def test_stop_codon_minus_good():
    # CDS [2001,2100): 3' coding end = low genomic [2001,2004); RC(TTA)=TAA
    g = genome_with({2097: "CAT", 2001: "TTA"})
    assert check_stop_codon(MINUS_TX(), g) is None


def test_stop_codon_minus_bad():
    g = genome_with({2097: "CAT", 2001: "GGG"})  # RC(GGG) = CCC
    assert check_stop_codon(MINUS_TX(), g).name == "NO_STOP"


def test_stop_codon_skipped_for_partial():
    g = genome_with({1000: "ATG", 1096: "CCC"})
    assert check_stop_codon(PLUS_TX(partial=True), g) is None


# ---------------------------------------------------------------------------
# internal stop — both strands
# ---------------------------------------------------------------------------

def test_internal_stop_plus_clean():
    # Terminal stop at [1096,1099) is excluded from internal-stop check
    g = genome_with({1000: "ATG", 1096: "TAA"})
    assert check_internal_stops(PLUS_TX(), g) is None


def test_internal_stop_plus_present():
    g = genome_with({1000: "ATG", 1003: "TAA", 1096: "TAA"})
    assert check_internal_stops(PLUS_TX(), g).name == "INTERNAL_STOP"


def test_internal_stop_minus_clean():
    # Terminal stop at [2001,2004) is excluded from internal-stop check
    g = genome_with({2097: "CAT", 2001: "TTA"})
    assert check_internal_stops(MINUS_TX(), g) is None


def test_internal_stop_minus_present():
    # coding codon #1 = RC(genome[2094:2097]); RC(TTA) = TAA -> internal stop
    g = genome_with({2097: "CAT", 2094: "TTA", 2001: "TTA"})
    assert check_internal_stops(MINUS_TX(), g).name == "INTERNAL_STOP"


def test_internal_stop_skipped_without_genome():
    assert check_internal_stops(PLUS_TX(), None) is None


# ---------------------------------------------------------------------------
# extract_cds_sequence — both strands + None
# ---------------------------------------------------------------------------

def test_extract_cds_sequence_plus():
    g = genome_with({1000: "ATGAAATTT"})
    tx = make_tx([(1000, 1010)], "+", cds=((1000, 1009, 0),))
    assert extract_cds_sequence(tx, g) == "ATGAAATTT"


def test_extract_cds_sequence_minus():
    g = genome_with({2000: "ATGAAATTT"})
    tx = make_tx([(2000, 2010)], "-", cds=((2000, 2009, 0),))
    assert extract_cds_sequence(tx, g) == reverse_complement("ATGAAATTT")
    assert extract_cds_sequence(tx, g) == "AAATTTCAT"


def test_extract_cds_sequence_multi_exon_plus():
    g = genome_with({1000: "ATG", 1100: "GGG"})
    tx = make_tx([(1000, 1030), (1100, 1130)], "+", cds=((1000, 1030, 0), (1100, 1112, 0)))
    seq = extract_cds_sequence(tx, g)
    assert seq.startswith("ATG")
    assert len(seq) == 30 + 12


def test_extract_cds_sequence_none_when_no_cds():
    tx = make_tx([(1000, 1100)], "+", cds=None)
    assert extract_cds_sequence(tx, None) is None


# ---------------------------------------------------------------------------
# short CDS / short exon / long intron (coordinate-only)
# ---------------------------------------------------------------------------

def test_short_cds_flagged():
    assert check_short_cds(PLUS_TX(), threshold=300).name == "SHORT_CDS"


def test_short_cds_not_flagged():
    assert check_short_cds(PLUS_TX(), threshold=10) is None


def test_short_cds_skipped_without_cds():
    assert check_short_cds(make_tx([(1000, 1100)], "+"), threshold=300) is None


def test_short_exon_flagged():
    tx = make_tx([(1000, 1005)], "+")
    assert check_short_exons(tx, threshold=10).name == "SHORT_EXON"


def test_short_exon_not_flagged():
    assert check_short_exons(PLUS_TX(), threshold=10) is None


def test_long_intron_flagged():
    tx = make_tx([(1000, 1100), (200000, 200100)], "+")
    assert check_long_introns(tx, threshold=100_000).name == "LONG_INTRON"


def test_long_intron_not_flagged():
    tx = make_tx([(1000, 1100), (1300, 1400)], "+")
    assert check_long_introns(tx, threshold=100_000) is None


def test_long_intron_minus_flagged():
    tx = make_tx([(1000, 1100), (200000, 200100)], "-")
    assert check_long_introns(tx, threshold=100_000).name == "LONG_INTRON"


# ---------------------------------------------------------------------------
# defense-in-depth helpers (return bool, emit no flag)
# ---------------------------------------------------------------------------

def test_check_cds_mod3_true_for_valid_complete():
    assert check_cds_mod3(PLUS_TX()) is True


def test_check_cds_mod3_true_for_partial():
    assert check_cds_mod3(PLUS_TX(partial=True, cds=((1000, 1098, 0),))) is True


def test_check_cds_mod3_true_when_no_cds():
    assert check_cds_mod3(make_tx([(1000, 1100)], "+")) is True


def test_check_cds_within_exons_true():
    tx = make_tx([(1000, 1200), (1300, 1500)], "+", cds=((1050, 1200, 0), (1300, 1450, 0)))
    assert check_cds_within_exons(tx) is True


def test_check_cds_within_exons_true_no_cds():
    assert check_cds_within_exons(make_tx([(1000, 1100)], "+")) is True


# ---------------------------------------------------------------------------
# validate_gene integration
# ---------------------------------------------------------------------------

def test_validate_gene_plus_clean_no_structural_flags():
    g = genome_with({1000: "ATG", 1096: "TAA"})
    gene = make_gene(PLUS_TX())
    flags = names(validate_gene(gene, g, short_cds_threshold=0))
    assert "NO_START" not in flags
    assert "NO_STOP" not in flags
    assert "INTERNAL_STOP" not in flags


def test_validate_gene_minus_clean_no_structural_flags():
    g = genome_with({2097: "CAT", 2001: "TTA"})
    gene = make_gene(MINUS_TX())
    flags = names(validate_gene(gene, g, short_cds_threshold=0))
    assert {"NO_START", "NO_STOP", "INTERNAL_STOP"}.isdisjoint(flags)


def test_validate_gene_plus_all_broken():
    g = genome_with({1000: "CCC", 1003: "TAA", 1096: "CCC"})
    gene = make_gene(PLUS_TX())
    flags = names(validate_gene(gene, g, short_cds_threshold=0))
    assert {"NO_START", "NO_STOP", "INTERNAL_STOP"} <= flags


def test_validate_gene_genome_none_skips_sequence_checks():
    gene = make_gene(PLUS_TX())  # would be bad-start if checked, but genome=None
    flags = names(validate_gene(gene, None, short_cds_threshold=300))
    assert "NO_START" not in flags
    assert "SHORT_CDS" in flags  # coordinate check still runs


def test_validate_gene_no_cds_skips_cds_checks():
    tx = make_tx([(1000, 1005)], "+")  # no CDS, tiny exon
    gene = make_gene(tx)
    flags = names(validate_gene(gene, None))
    assert "SHORT_EXON" in flags
    assert "NO_START" not in flags
    assert "SHORT_CDS" not in flags


def test_validate_gene_partial_orf_flag_and_skipped_termini():
    g = genome_with({1000: "CCC"})  # bad start, but partial
    tx = PLUS_TX(partial=True, cds=((1000, 1098, 0),))
    gene = make_gene(tx)
    flags = names(validate_gene(gene, g, short_cds_threshold=0))
    assert "PARTIAL_ORF" in flags
    assert "NO_START" not in flags
    assert "NO_STOP" not in flags


def test_validate_gene_partial_still_flags_internal_stop():
    g = genome_with({1000: "ATG", 1003: "TAA"})
    tx = PLUS_TX(partial=True)
    gene = make_gene(tx)
    flags = names(validate_gene(gene, g, short_cds_threshold=0))
    assert "INTERNAL_STOP" in flags
    assert "PARTIAL_ORF" in flags


def test_validate_gene_validates_primary_transcript():
    primary = make_tx([(1000, 1100)], "+", cds=((1000, 1099, 0),), tid="HFG_00001.1")
    alt = make_tx([(1500, 1600)], "+", cds=((1500, 1599, 0),), tid="HFG_00001.2")
    g = genome_with({1000: "CCC", 1500: "ATG"})  # primary bad, alt good
    gene = make_gene(primary, transcripts=[primary, alt], primary="HFG_00001.1")
    assert "NO_START" in names(validate_gene(gene, g, short_cds_threshold=0))


def test_validate_gene_ignores_non_primary_transcript():
    primary = make_tx([(1500, 1600)], "+", cds=((1500, 1599, 0),), tid="HFG_00001.2")
    alt = make_tx([(1000, 1100)], "+", cds=((1000, 1099, 0),), tid="HFG_00001.1")
    g = genome_with({1000: "CCC", 1500: "ATG", 1596: "TAA"})  # alt bad, primary good
    gene = make_gene(primary, transcripts=[alt, primary], primary="HFG_00001.2")
    assert "NO_START" not in names(validate_gene(gene, g, short_cds_threshold=0))


def test_validate_gene_flags_deduplicated():
    g = genome_with({1000: "ATG", 1096: "TAA"})
    flags = validate_gene(make_gene(PLUS_TX()), g, short_cds_threshold=300)
    assert len(flags) == len({f.name for f in flags})


# ---------------------------------------------------------------------------
# validate_all
# ---------------------------------------------------------------------------

def test_validate_all_keyed_by_gene_id():
    g = genome_with({1000: "ATG", 1096: "TAA"})
    gene1 = make_gene(PLUS_TX())
    gene2 = ReconciledGene(
        gene_id="HFG_00002",
        seqid="chr1",
        start=1000,
        end=1005,
        strand="+",
        tier=4,
        transcripts=[make_tx([(1000, 1005)], "+", tid="HFG_00002.1")],
        primary_transcript_id="HFG_00002.1",
        classification=LocusClassification("HFG_00002", "SILENT"),
        origin="helixer_backstop",
    )
    result = validate_all([gene1, gene2], g)
    assert set(result) == {"HFG_00001", "HFG_00002"}
    assert "SHORT_EXON" in {f.name for f in result["HFG_00002"]}


def test_validate_all_empty():
    assert validate_all([]) == {}


# ===========================================================================
# Phase 27 D2 — ambiguous (N/IUPAC) start/stop codons → AMBIGUOUS_CODON
# (assessment §1.1). Both strands. An N in a terminal codon is indeterminate,
# never a false NO_START/NO_STOP; an N internal codon is not a premature stop.
# ===========================================================================

def test_start_codon_n_is_ambiguous_not_no_start_plus():
    # plus start codon = genome[1000:1003] = "ANG" (masked middle base).
    g = genome_with({1000: "ANG", 1096: "TAA"})
    flag = check_start_codon(PLUS_TX(), g)
    assert flag is not None and flag.name == "AMBIGUOUS_CODON"


def test_start_codon_n_is_ambiguous_not_no_start_minus():
    # minus start codon = RC(genome[2097:2100]); plant "CNT" -> RC = "ANG".
    g = genome_with({2097: "CNT", 2001: "TTA"})
    flag = check_start_codon(MINUS_TX(), g)
    assert flag is not None and flag.name == "AMBIGUOUS_CODON"


def test_stop_codon_n_is_ambiguous_not_no_stop_plus():
    # Stop-inclusive: last 3 CDS bases = [1096,1099); plant "TNA" -> ambiguous.
    g = genome_with({1000: "ATG", 1096: "TNA"})
    flag = check_stop_codon(PLUS_TX(), g)
    assert flag is not None and flag.name == "AMBIGUOUS_CODON"


def test_stop_codon_n_is_ambiguous_not_no_stop_minus():
    # Stop-inclusive: 3' coding end [2001,2004); RC("TNA") = "TNA" -> ambiguous.
    g = genome_with({2097: "CAT", 2001: "TNA"})
    flag = check_stop_codon(MINUS_TX(), g)
    assert flag is not None and flag.name == "AMBIGUOUS_CODON"


def test_internal_stop_with_n_codon_not_premature_plus():
    # CDS [1000,1099). codon #2 (1003) = "NNN" (unknown, not a stop); terminal
    # codon at 1096 should not be flagged. No INTERNAL_STOP despite the N codon.
    g = genome_with({1000: "ATG", 1003: "NNN", 1096: "TAA"})
    assert check_internal_stops(PLUS_TX(), g) is None


def test_internal_stop_real_stop_survives_n_elsewhere_plus():
    # A genuine internal stop at codon #2 is still detected even when a *later*
    # codon contains N (previously the whole check was skipped on any N).
    g = genome_with({1000: "ATG", 1003: "TAA", 1006: "NNN", 1096: "TAA"})
    flag = check_internal_stops(PLUS_TX(), g)
    assert flag is not None and flag.name == "INTERNAL_STOP"


def test_clean_gene_unchanged_no_ambiguous_flag():
    # A clean ACGT gene gets no AMBIGUOUS_CODON and no NO_START/NO_STOP.
    g = genome_with({1000: "ATG", 1096: "TAA"})
    gene = make_gene(PLUS_TX())
    flags = names(validate_gene(gene, g))
    assert "AMBIGUOUS_CODON" not in flags
    assert "NO_START" not in flags
    assert "NO_STOP" not in flags


def test_validate_gene_ambiguous_start_emits_flag_not_no_start():
    # End-to-end through validate_gene: N start codon -> AMBIGUOUS_CODON only.
    g = genome_with({1000: "ANG", 1096: "TAA"})
    gene = make_gene(PLUS_TX())
    flags = names(validate_gene(gene, g))
    assert "AMBIGUOUS_CODON" in flags
    assert "NO_START" not in flags


def test_per_scaffold_table_applies_in_validate_gene():
    # Phase 27 D3: a {seqid: table_id} map changes the genetic code for the
    # scaffold. With table 4 (TGA->W) a "TGA" internal codon is no longer a
    # premature stop. CDS [1000,1099): codon #2 at 1003 = "TGA".
    g = genome_with({1000: "ATG", 1003: "TGA", 1096: "TAA"})
    gene = make_gene(PLUS_TX())
    # default table 1 -> internal stop
    assert "INTERNAL_STOP" in names(validate_gene(gene, g))
    # per-scaffold table 4 for chr1 -> no internal stop
    flags = names(validate_gene(gene, g, transl_table_map={"chr1": 4}))
    assert "INTERNAL_STOP" not in flags
