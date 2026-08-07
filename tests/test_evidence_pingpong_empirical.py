"""Ping-pong tested against the locus's own background, not a uniform null.

The 1/20 null assumed every 5'-5' offset is equally likely. Offsets are not
uniform -- coverage structure alone puts pairs at some offsets more than
others -- which is why the reference implementation estimates the background
from the other offsets instead of assuming it.

Here the rate comes from the same locus's offsets 1-20 excluding 10, so a locus
whose pairs are spread thin and one whose pairs pile up everywhere are held to
their own standard.
"""
from evidence import pingpong_score_empirical


def _flat(value, target=0):
    """Counts for offsets 1..20, all `value`, with `target` extra at offset 10."""
    c = {o: float(value) for o in range(1, 21)}
    c[10] = float(value + target)
    return c


def test_a_peak_far_above_its_own_background_scores_high():
    assert pingpong_score_empirical(_flat(5, target=45)) > 5


def test_no_peak_stays_below_the_evidence_threshold():
    # a flat spectrum gives p ~ 0.5, i.e. a score near zero -- not exactly zero,
    # which is what a p-value should do when the observation is unremarkable
    assert pingpong_score_empirical(_flat(5)) < 1.0


def test_the_same_peak_counts_less_on_a_noisy_locus():
    """50 pairs at offset 10 is strong against a background of 5, weak against 45."""
    quiet = pingpong_score_empirical(_flat(5, target=45))
    noisy = pingpong_score_empirical(_flat(45, target=5))
    assert quiet > noisy


def test_the_uniform_null_could_not_tell_them_apart():
    """Both loci have the same 10-nt fraction; only the background differs."""
    from evidence import pingpong_score
    quiet_frac = 50 / (5 * 19 + 50)
    noisy_frac = 50 / (45 * 19 + 50)
    assert quiet_frac > noisy_frac          # the old score sees only this
    # and it has no way to know the second locus's background is nine times higher


def test_empty_and_degenerate_inputs():
    assert pingpong_score_empirical({}) == 0.0
    assert pingpong_score_empirical({10: 0.0}) == 0.0
    assert pingpong_score_empirical({o: 0.0 for o in range(1, 21)}) == 0.0
