"""phase_score must not depend on how deep the library was sequenced.

It is `ratio * number * ln(mass)` on raw 1/NH mass, so it correlates 0.572 with
log(reads) -- median 0.89 in loci under 50 reads against 10.88 above 5000. That
score feeds a sigmoid whose midpoint (`phase_x0: 5.0`) was calibrated on a
19.9M-read library, so the same locus in a deeper run scores higher for no
biological reason. sRNAminer normalises first (`PHASMinerCLI --normLibSize
20000000`), and 20M is kept here because it is where our calibration already
sits: the plant library is 19.93M, so the shipped thresholds keep their meaning.
"""
import math

from signals import NORM_LIB_SIZE, phase_metrics, phase_test


def _phased(n_cycles, abundance=4.0, phase=21):
    return [(c * phase, abundance) for c in range(n_cycles)]


def test_the_reference_library_size_is_the_one_we_calibrated_on():
    assert NORM_LIB_SIZE == 20_000_000


def test_the_same_locus_scores_the_same_in_a_deeper_library():
    """Twice the depth, twice the mass, same biology -- same score."""
    shallow = phase_metrics(_phased(8, abundance=4.0), 0, 21, lib_size=10_000_000)
    deep = phase_metrics(_phased(8, abundance=8.0), 0, 21, lib_size=20_000_000)
    assert abs(shallow[0] - deep[0]) < 1e-6


def test_raw_mass_still_used_when_the_library_size_is_unknown():
    """Without lib_size there is nothing to normalise by; do not invent one."""
    positions = _phased(8, abundance=4.0)
    assert phase_metrics(positions, 0, 21)[0] == phase_metrics(
        positions, 0, 21, lib_size=None)[0]


def test_a_library_at_the_reference_size_is_left_alone():
    """The factor is 1.0 there, so the calibration library is untouched."""
    positions = _phased(8, abundance=4.0)
    assert (phase_metrics(positions, 0, 21, lib_size=NORM_LIB_SIZE)[0]
            == phase_metrics(positions, 0, 21)[0])


def test_the_normalised_mass_is_what_enters_the_logarithm():
    positions = _phased(6, abundance=5.0)
    lib = 5_000_000
    ratio, number = 1.0, 6
    mass = 6 * 5.0 * (NORM_LIB_SIZE / lib)
    expected = round(ratio * number * math.log(mass), 3)
    assert phase_metrics(positions, 0, 21, lib_size=lib)[0] == expected


def test_number_and_ratio_do_not_move_with_the_library():
    """Only the abundance term is normalised; the register geometry is not."""
    a = phase_metrics(_phased(7, abundance=4.0), 0, 21, lib_size=5_000_000)
    b = phase_metrics(_phased(7, abundance=4.0), 0, 21)
    assert a[1] == b[1] and a[2] == b[2]


def test_phase_test_passes_the_library_size_through():
    """The wiring, not the formula: a deeper library must reach phase_metrics."""
    shallow = phase_test(_phased(8, abundance=4.0), 0, 21, lib_size=10_000_000)
    deep = phase_test(_phased(8, abundance=8.0), 0, 21, lib_size=20_000_000)
    assert abs(shallow.score - deep.score) < 1e-6


def test_build_record_normalises_by_the_library_it_was_given():
    """The wiring through features, not the formula in signals.

    The genomic 5' null was threaded through every layer except the final call
    and stayed off while both validations came back byte-identical. What the
    record ends up carrying is the only thing that proves the size arrived.
    """
    from features import build_record
    from samio import Aln

    cfg = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}
    tails = ["ACGT", "TGCA", "GTAC", "CATG", "AGCT", "TCGA", "GATC", "CTAG"]
    reads = [Aln("c", i * 21, "+", (tails[i % len(tails)] * 6)[:21], 1)
             for i in range(8)]

    shallow = build_record(("c", 0, 400), reads, 5_000_000, cfg)
    deep = build_record(("c", 0, 400), reads, NORM_LIB_SIZE, cfg)

    assert shallow.phase_number == deep.phase_number
    assert shallow.phase_ratio == deep.phase_ratio
    assert shallow.phase_score > deep.phase_score
