"""Per-strand size histogram and read-start shape features."""
from collections import namedtuple

import features
from reporter import write_tsv, _COLS
from records import ClusterRecord

CAln = namedtuple("CAln", "chrom pos seq nh strand")
CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _build(alns):
    return features.build_record(("chr1", 0, 1000), alns, 1e6, CFG)


def _sz(r):
    return {f"sz{s}_{L}": getattr(r, f"sz{s}_{L}")
            for s in "pm" for L in range(18, 36)}


def test_size_histogram_by_strand():
    alns = [CAln("chr1", 100 + 30 * i, "A" * 21, 1, "+") for i in range(3)]
    alns.append(CAln("chr1", 300, "A" * 24, 2, "-"))
    r = _build(alns)
    sz = _sz(r)
    assert sz["szp_21"] == round(3 / 3.5, 4)
    assert sz["szm_24"] == round(0.5 / 3.5, 4)
    for k, v in sz.items():
        if k not in ("szp_21", "szm_24"):
            assert v == 0
    assert abs(sum(sz.values()) - r.frac_in_window) < 1e-3


def test_shape_features_single_position():
    alns = [CAln("chr1", 100, "A" * 22, 1, "+") for _ in range(5)]
    r = _build(alns)
    assert r.start_entropy == 0.0
    assert r.top_start_frac == 1.0
    assert all(getattr(r, f"szm_{L}") == 0 for L in range(18, 36))


def test_shape_features_uniform():
    alns = [CAln("chr1", 100 + 30 * i, "A" * 22, 1, "+") for i in range(4)]
    r = _build(alns)
    assert r.start_entropy == 1.0
    assert r.top_start_frac == 0.25
    assert r.starts_per_kb is not None


def test_tsv_has_new_columns(tmp_path):
    p = tmp_path / "c.tsv"
    write_tsv([ClusterRecord(chrom="chr1", start=0, end=10, name="x")], str(p))
    header = p.read_text().splitlines()[0].split("\t")
    i = header.index("read_span")
    expect = ([f"szp_{L}" for L in range(18, 36)]
              + [f"szm_{L}" for L in range(18, 36)]
              + ["start_entropy", "top_start_frac", "starts_per_kb"])
    assert header[i + 1:i + 1 + len(expect)] == expect
    row = p.read_text().splitlines()[1].split("\t")
    assert row[header.index("start_entropy")] == "NA"
