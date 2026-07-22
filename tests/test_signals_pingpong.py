from signals import pingpong_zscore


def test_pingpong_peak_at_10_gives_high_z():
    counts = {o: (o % 4) + 1 for o in range(1, 21)}  # background with real variance (1..4)
    counts[10] = 100
    z, pairs = pingpong_zscore(counts)
    assert pairs == 100
    assert z > 3.0


def test_pingpong_flat_gives_zero_z():
    counts = {o: 5 for o in range(1, 21)}
    z, pairs = pingpong_zscore(counts)
    assert z == 0.0
    assert pairs == 5


def test_pingpong_insufficient_background():
    z, pairs = pingpong_zscore({10: 50})
    assert z == 0.0
    assert pairs == 50
