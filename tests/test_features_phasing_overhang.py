"""FIX 3 (TDD): 2-nt overhang correction in phasing.

Dicer's phased products have a 2-nt 3' overhang, so the 5' end of the - strand
is shifted (+2) relative to the + strand's register. Without the correction, a
double-stranded phased locus splits the abundance across two registers (ratio ~0.5).
With the correction, sense and antisense fall into the SAME dominant register
(ratio == 1.0).
"""
from samio import read_alignments
from features import build_record

CONFIG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def test_overhang_pools_both_strands_in_one_register():
    alns = list(read_alignments("tests/data/phasing_overhang.sam"))
    # region_start=1000, detected pl = 21 (all reads are 21 nt)
    rec = build_record(("chr1", 1000, 1200), alns, lib_size=1_000_000, config=CONFIG)
    assert rec.dicer_call == "21"
    # 3 phase cycles, covered by + and - together after the 2-nt correction
    assert rec.phase_number == 3
    # without the correction, +/- would fall into different registers -> ratio == 0.5;
    # with the correction both pool into the same register -> ratio == 1.0
    assert rec.phase_ratio == 1.0
