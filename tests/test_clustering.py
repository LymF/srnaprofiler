from clustering import call_islands

def test_call_islands_merges_within_pad():
    depth = {10: 5, 11: 5, 12: 5,   15: 5, 16: 5}  # gap 12->15 = 2 vazios
    # pad=2 funde os dois blocos; pad=1 separa
    assert call_islands(depth, mincov=1, pad=2) == [(10, 17)]
    assert call_islands(depth, mincov=1, pad=1) == [(10, 13), (15, 17)]

def test_call_islands_respects_mincov():
    depth = {10: 1, 11: 3, 12: 1}
    assert call_islands(depth, mincov=2, pad=0) == [(11, 12)]

def test_call_islands_empty():
    assert call_islands({}, mincov=1, pad=0) == []
