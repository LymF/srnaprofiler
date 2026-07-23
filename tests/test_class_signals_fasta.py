from samio import Aln
from class_signals import strand_sequences


def test_strand_sequences_orient_trim_uniq():
    alns = [
        Aln("chr", 0, "+", "TACGTA", 1),        # + -> TACGTA, trim4 -> UACG
        Aln("chr", 0, "+", "TACGTT", 1),        # + -> TACGTT, trim4 -> UACG (dup)
        Aln("chr", 0, "-", "AACGTA", 1),         # revcomp = TACGTT, trim4 -> UACG (dup)
    ]
    seqs = strand_sequences(alns, "+", size_lo=6, size_hi=6, trim=4)
    assert seqs == ["UACG"]                      # collapsed to one distinct sequence


def test_strand_sequences_minus_uses_revcomp():
    alns = [Aln("chr", 0, "-", "AAATTT", 1)]     # revcomp = AAATTT -> trim3 -> AAA
    seqs = strand_sequences(alns, "-", size_lo=6, size_hi=6, trim=3)
    assert seqs == ["AAA"]
