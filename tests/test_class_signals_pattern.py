import statistics
from samio import Aln
from class_signals import pattern_row


def test_pattern_row_layout_and_zscore():
    # one 21-nt + read (size 21 -> index 6) and one 24-nt - read (size 24 -> +21+9=30)
    alns = [Aln("chr", 0, "+", "A" * 21, 1), Aln("chr", 0, "-", "A" * 24, 2)]
    row = pattern_row(alns, minlen=15, maxlen=35)
    assert len(row) == 42
    # raw counts before z-score: idx 6 = 1.0 (+), idx 30 = 0.5 (-)
    # z-score is monotonic, so the two non-zero cells are the two largest values
    top2 = sorted(range(42), key=lambda i: row[i], reverse=True)[:2]
    assert set(top2) == {6, 30}


def test_pattern_row_all_zero_is_flat():
    row = pattern_row([], minlen=15, maxlen=35)
    assert len(row) == 42
    assert all(v == 0.0 for v in row)
