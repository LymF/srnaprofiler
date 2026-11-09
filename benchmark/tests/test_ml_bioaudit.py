import gzip

import numpy as np
import pandas as pd

from benchmark.ml import bioaudit as ba
from benchmark.ml import evaluate as ev
from benchmark.scripts import build_te


def _table():
    return pd.DataFrame({
        "organism": "Drosophila melanogaster",
        "chrom": "2L",
        "start": [100, 100, 100, 100, 100, 100, 100, 100],
        "end":   [130, 130, 130, 130, 130, 130, 130, 130],
        "label": [None, None, None, "piRNA", None, None, None, None],
        "carved": [False, False, False, False, True, False, False, False],
        "repro_frac": [0.5, 0.49, 1.0, 1.0, 1.0, 1.0, np.nan, 1.0],
    })


def test_supported_repro_and_te_cases():
    t = _table()
    # row 0: piRNA repro 0.5 but no TE -> False; use a TE touching by exactly 1 bp
    te = ba.TEIndex([("2L", 129, 500)])            # overlaps [100,130) by 1 bp
    pred = ["piRNA", "piRNA", "siRNA", "piRNA", "miRNA", "ambiguous", "siRNA", None]
    out = ba.supported(t, pred, te)
    assert out.tolist() == [True,    # piRNA, repro 0.5, 1 bp TE
                            False,   # repro 0.49
                            True,    # siRNA: repro only
                            False,   # annotated row excluded
                            False,   # carved excluded
                            False,   # ambiguous
                            False,   # NaN repro
                            False]   # no call
    far = ba.TEIndex([("2L", 130, 500)])           # adjacent, 0 bp
    assert not ba.supported(t, pred, far)[0]       # 0 bp (adjacent) is not an overlap
    assert ba.supported(t, pred, None)[0]           # no TE configured: repro only
    assert ba.supported(t, pred, {"Other sp": far})[0]   # organism without TE: repro only
    assert not ba.supported(t, pred, {"Drosophila melanogaster": far})[0]
    rate, n = ba.supported_rate(t, pred, te)
    assert n == 4 and rate == 0.5                  # unannotated calls: rows 0, 1, 2, 6
    assert ba.supported_rate(t, [None] * 8, te) == (None, 0)


def test_fold_metrics_reports_supported_rates():
    t = _table()
    te = ba.TEIndex([("2L", 129, 500)])
    pred = np.array(["piRNA", "piRNA", "siRNA", "piRNA", "miRNA", "ambiguous", "siRNA", None],
                    dtype=object)
    rule = np.array(["ambiguous", "piRNA", "piRNA", "piRNA", "miRNA", "ambiguous", "siRNA",
                     "ambiguous"], dtype=object)
    m = ev.fold_metrics(t, pred, np.full(8, 0.9), rule, None, te)
    assert m["unann_supported_rate"] == 0.5 and m["n_unann_calls"] == 4
    assert m["rules_n_unann_calls"] == 3           # rows 1, 2, 6
    assert m["rules_unann_supported_rate"] == 1 / 3   # only row 2 (row 1 repro 0.49, row 6 NaN)


def test_te_build_and_load(tmp_path):
    rows = [
        ["585", "1", "0", "0", "0", "chr2L", "100", "200", "-1", "+", "Gypsy", "LTR", "Gypsy", "0", "0", "0", "1"],
        ["585", "1", "0", "0", "0", "chrX", "5", "50", "-1", "+", "Jockey", "LINE", "Jockey", "0", "0", "0", "2"],
        ["585", "1", "0", "0", "0", "chr2L", "10", "20", "-1", "+", "(CA)n", "Simple_repeat", "x", "0", "0", "0", "3"],
        ["585", "1", "0", "0", "0", "chr3R", "1", "9", "-1", "+", "Foo", "LTR?", "x", "0", "0", "0", "4"],
        ["585", "1", "0", "0", "0", "chr3R", "30", "90", "-1", "+", "Bar", "Unknown", "x", "0", "0", "0", "5"],
    ]
    src = tmp_path / "rmsk.txt.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("\n".join("\t".join(r) for r in rows) + "\n")
    out = tmp_path / "rmsk_te.bed"
    build_te.main([str(src), "-o", str(out)])
    assert out.read_text().splitlines() == ["2L\t100\t200", "3R\t30\t90", "X\t5\t50"]
    te = ba.load_te(out)
    assert te.overlaps("2L", 199, 300) and not te.overlaps("2L", 200, 300)
    import pytest
    with pytest.raises(FileNotFoundError, match="Drosophila melanogaster.*d/rmsk_te.bed"):
        ba.load_te_indexes({"te": {"Drosophila melanogaster": "d/rmsk_te.bed"}}, tmp_path)
    # no entry, or entry for an organism not in the table: no error, repro only
    assert ba.load_te_indexes({"te": {}}, tmp_path) == {}
    assert ba.load_te_indexes({"te": {"Drosophila melanogaster": "d/x.bed"}}, tmp_path,
                              organisms={"Bombyx mori"}) == {}
    (tmp_path / "d").mkdir()
    out.rename(tmp_path / "d" / "rmsk_te.bed")
    assert set(ba.load_te_indexes({"te": {"Drosophila melanogaster": "d/rmsk_te.bed"}},
                                  tmp_path)) == {"Drosophila melanogaster"}
