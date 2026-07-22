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

def test_pingpong_sparse_dense_locus_returns_zero():
    # ultra-dense narrow locus: huge peak, background almost all zero (1 non-zero)
    z, pairs = pingpong_zscore({10: 1_000_000, 5: 3})
    assert z == 0.0            # < min_bg non-zero background offsets -> not meaningful
    assert pairs == 1_000_000

def test_pingpong_zscore_is_capped():
    # 3 non-zero background offsets but a huge peak -> capped, not thousands
    z, _ = pingpong_zscore({10: 1_000_000, 1: 1, 2: 1, 3: 1})
    assert z == 100.0
