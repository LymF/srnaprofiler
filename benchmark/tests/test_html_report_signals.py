"""Per-cluster figure data for the HTML report (single pass over the BAM)."""
import json
import math
import random
import weakref
from collections import namedtuple

import pysam
import pytest

import bamutil
import class_signals
import report_signals as rs
from samio import Aln

Cluster = namedtuple("Cluster", "name chrom start end")


def _open(tmp_path, reads, chroms=(("c1", 5000),)):
    bam = tmp_path / "t.bam"
    bamutil.make_seq_bam(bam, reads, chroms)
    return pysam.AlignmentFile(str(bam), "rb")


def _sig(tmp_path, reads, start=0, end=5000, chroms=(("c1", 5000),), **kw):
    with _open(tmp_path, reads, chroms) as fh:
        return rs.cluster_signals(fh, "c1", start, end, xm_as_nh=False, **kw)


def _alns(reads):
    return [Aln(c, p, s, seq, nh or 1) for c, p, seq, s, nh, *_ in reads]


def _numbers(x):
    if isinstance(x, dict):
        for v in x.values():
            yield from _numbers(v)
    elif isinstance(x, list):
        for v in x:
            yield from _numbers(v)
    elif isinstance(x, (int, float)) and not isinstance(x, bool):
        yield x


def _seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _random_reads(seed=1, n=400):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        L = rng.randint(16, 30)
        out.append(("c1", rng.randint(0, 600), _seq(rng, L), rng.choice("+-"),
                    rng.choice([1, 1, 2, 4])))
    return out


def test_matches_class_signals_reference(tmp_path):
    reads = _random_reads()
    # reference input: the reads a fetch of [50, 500) returns
    alns = [a for a in _alns(reads) if a.pos < 500 and a.pos + len(a.seq) > 50]
    d = _sig(tmp_path, reads, 50, 500)
    assert 0 < d["n_alns"] == len(alns) < len(reads)
    comp = class_signals.base_composition(alns, 18, 35)
    for L in range(18, 36):
        for s in "+-":
            for b in "ACGT":
                got = d["sizes"]["U" if b == "T" else b][s][L - 18]
                assert got == pytest.approx(round(comp[L][s][b], 3), abs=1e-3)
    for s in "+-":
        ref = class_signals.phasing_series(alns, s)
        assert len(d["phasing"][s]) == 100
        assert d["phasing"][s] == pytest.approx(
            [round(ref.get(k, 0.0), 3) for k in range(1, 101)], abs=1e-3)
    pp = class_signals.pingpong_offsets(alns, hi=30)
    assert d["pingpong"] == pytest.approx([round(pp.get(k, 0.0), 3) for k in range(1, 31)],
                                          abs=1e-3)
    cov = d["cov"]
    plus, minus = class_signals.coverage_track(alns, 50, 500)
    widths = [min(cov["bin"], 450 - i * cov["bin"]) for i in range(len(cov["+"]))]
    assert sum(v * w for v, w in zip(cov["+"], widths)) == pytest.approx(sum(plus), abs=0.5)
    assert sum(v * w for v, w in zip(cov["-"], widths)) == pytest.approx(sum(minus), abs=0.5)
    p5, m5 = class_signals.five_prime_ends(alns)
    assert sum(d["ends"]["+"]) == pytest.approx(
        sum(v for p, v in p5.items() if 50 <= p < 500), abs=0.01)
    assert sum(d["ends"]["-"]) == pytest.approx(
        sum(v for p, v in m5.items() if 50 <= p < 500), abs=0.01)


def test_nh_weighting(tmp_path):
    reads = [("c1", 100, "T" + "A" * 20, "+", 1),
             ("c1", 100, "T" + "A" * 20, "+", 2),
             ("c1", 200, "A" * 22, "-", 4),
             ("c1", 300, "G" * 21, "+", None)]  # no NH tag: unique
    d = _sig(tmp_path, reads)
    assert d["sizes"]["U"]["+"][21 - 18] == 1.5
    assert d["sizes"]["U"]["-"][22 - 18] == 0.25  # revcomp of A* starts with U
    assert d["sizes"]["G"]["+"][21 - 18] == 1.0
    assert d["n_alns"] == 4
    rows = {(r[0], r[2]): r[3] for r in d["seqs"]}
    assert rows[("U" + "A" * 20, "+")] == 1.5
    assert rows[("U" * 22, "-")] == 0.25


def test_xm_as_nh(tmp_path):
    bam = tmp_path / "t.bam"
    bamutil.make_seq_bam(bam, [("c1", 100, "A" * 21, "+", None)])
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        # rewrite with an XM tag (bowtie1 multiplicity)
        out = tmp_path / "x.bam"
        with pysam.AlignmentFile(str(out), "wb", template=fh) as w:
            for rec in fh:
                rec.set_tag("XM", 4)
                w.write(rec)
    pysam.index(str(out))
    with pysam.AlignmentFile(str(out), "rb") as fh:
        a = rs.cluster_signals(fh, "c1", 0, 500, xm_as_nh=True)
        b = rs.cluster_signals(fh, "c1", 0, 500, xm_as_nh=False)
    assert a["sizes"]["A"]["+"][3] == 0.25
    assert b["sizes"]["A"]["+"][3] == 1.0


def test_filter_matches_samio(tmp_path):
    reads = [("c1", 100, "A" * 21, "+", 1),
             ("c1", 110, "A" * 21, "+", 1, 256),    # secondary
             ("c1", 120, "A" * 21, "+", 1, 2048)]   # supplementary
    d = _sig(tmp_path, reads)
    assert d["n_alns"] == 1


def test_minus_strand_is_reverse_complemented_in_logo_and_seqs(tmp_path):
    # stored (reference orientation) T + A*21 -> RNA U*21 + A
    reads = [("c1", 100, "T" + "A" * 21, "-", 1)]
    d = _sig(tmp_path, reads)
    assert d["dominant"] == 22
    lm = d["logo"]["-"]
    assert len(lm) == 26 and lm[0]["U"] == 1.0 and lm[1]["U"] == 1.0
    assert lm[21]["A"] == 1.0
    assert lm[22] == {"A": 0, "C": 0, "G": 0, "U": 0}
    assert d["seqs"] == [["U" * 21 + "A", 22, "-", 1.0, "U"]]
    assert d["sizes"]["U"]["-"][22 - 18] == 1.0
    # 5' end of a - read is its last reference base
    assert d["ends"]["-"][121 // d["ends"]["bin"]] == 1.0


def test_sense_only_cluster_has_zero_antisense(tmp_path):
    reads = [("c1", 100 + i, "T" + "C" * 20, "+", 1) for i in range(20)]
    d = _sig(tmp_path, reads, 0, 1000)
    for b in "ACGU":
        assert not any(d["sizes"][b]["-"])
    assert not any(d["phasing"]["-"]) and not any(d["cov"]["-"]) and not any(d["ends"]["-"])
    assert all(sum(p.values()) == 0 for p in d["logo"]["-"])
    assert d["logo"]["+"][0]["U"] == 1.0
    assert all(math.isfinite(x) for x in _numbers(d))
    json.dumps(d, allow_nan=False)


def test_bins_bounded(tmp_path):
    d = _sig(tmp_path, [("c1", 10, "A" * 20, "+", 1)], 0, 37)
    assert d["cov"]["bin"] == 1 and len(d["cov"]["+"]) == 37 <= 37
    assert len(d["ends"]["+"]) == 37 and d["ends"]["bin"] == 1
    assert d["cov"]["+"][10] == 1.0 and d["cov"]["+"][29] == 1.0 and d["cov"]["+"][30] == 0.0
    big = [("c1", 500_000, "A" * 20, "+", 1), ("c1", 999_990, "A" * 20, "-", 1)]
    d = _sig(tmp_path, big, 0, 1_000_000, chroms=(("c1", 1_000_000),))
    assert len(d["cov"]["+"]) <= 400 and len(d["cov"]["-"]) <= 400
    assert len(d["ends"]["+"]) == len(d["cov"]["+"])
    assert d["cov"]["bin"] == 2500
    # mean depth: 20 bases of depth 1 in a 2500-bp bin
    assert d["cov"]["+"][200] == round(20 / 2500, 3)
    assert d["cov"]["-"][-1] == round(10 / 2500, 3)  # clipped at the window end


def test_partial_last_bin_is_a_mean_over_its_width(tmp_path):
    reads = [("c1", 1000, "A" * 3, "+", 1)]
    d = _sig(tmp_path, reads, 0, 1003, max_bins=100)
    assert d["cov"]["bin"] == 11 and len(d["cov"]["+"]) == 92
    # last bin spans 1001-1002 (2 bases), both covered
    assert d["cov"]["+"][-1] == 1.0
    assert d["cov"]["+"][-2] == round(1 / 11, 3)


def test_seqs_top_and_fewer_than_ten(tmp_path):
    rng = random.Random(3)
    reads = []
    seqs = [_seq(rng, 21) for _ in range(15)]
    for i, s in enumerate(seqs):
        reads += [("c1", 100 + i, s, "+", 1)] * (i + 1)
    d = _sig(tmp_path, reads)
    assert len(d["seqs"]) == 10
    assert [r[3] for r in d["seqs"]] == [float(n) for n in range(15, 5, -1)]
    assert d["seqs"][0][0] == seqs[14].replace("T", "U")
    assert d["seqs"][0][4] == seqs[14].replace("T", "U")[0]
    few = [("c1", 100, "G" * 21, "+", 1), ("c1", 200, "C" * 21, "+", 1)]
    d = _sig(tmp_path, few)
    # ties broken by sequence
    assert [r[0] for r in d["seqs"]] == ["C" * 21, "G" * 21]
    assert len(d["seqs"]) == 2


def test_dominant_inside_window_and_given(tmp_path):
    reads = [("c1", 0, "G" * 15, "+", 1)] * 10 + [("c1", 0, "T" + "A" * 20, "+", 1)] * 2
    reads += [("c1", 0, "C" + "A" * 25, "+", 1)]
    d = _sig(tmp_path, reads)
    assert d["dominant"] == 21 and d["logo"]["+"][0]["U"] == 1.0
    d = _sig(tmp_path, reads, dominant=26)
    assert d["dominant"] == 26 and d["logo"]["+"][0]["C"] == 1.0


def test_zero_reads(tmp_path):
    d = _sig(tmp_path, [("c1", 4000, "A" * 20, "+", 1)], 0, 500)
    assert d["n_alns"] == 0 and d["dominant"] is None and d["seqs"] == []
    assert len(d["pingpong"]) == 30 and not any(d["pingpong"])
    for s in "+-":
        assert len(d["phasing"][s]) == 100 and not any(d["phasing"][s])
        assert len(d["cov"][s]) > 0 and not any(d["cov"][s]) and not any(d["ends"][s])
        assert len(d["logo"][s]) == 26 and all(sum(p.values()) == 0 for p in d["logo"][s])
        for b in "ACGU":
            assert len(d["sizes"][b][s]) == 18 and not any(d["sizes"][b][s])
    json.dumps(d, allow_nan=False)


class _Rec:
    """Minimal stand-in for a pysam record (weak-referenceable)."""
    flag = 0
    is_reverse = False
    reference_start = 100
    query_sequence = "T" + "A" * 20

    def has_tag(self, t):
        return t == "NH"

    def get_tag(self, t):
        return 2


class _Fh:
    def __init__(self, n):
        self.n, self.alive_max = n, 0

    def fetch(self, chrom, start, end):
        # Structural check: every yielded record is weakly tracked; if the
        # consumer kept a list, earlier records would stay alive.
        live = weakref.WeakSet()
        for _ in range(self.n):
            r = _Rec()
            live.add(r)
            self.alive_max = max(self.alive_max, len(live))
            yield r


def test_memory_is_bounded_by_distinct_keys():
    fh = _Fh(200_000)
    n, keys, seqs = rs._accumulate(fh, "c1", 0, 500, xm_as_nh=False)
    assert n == 200_000
    assert len(keys) == 1 and len(seqs) == 1  # one distinct key, not 200k entries
    assert fh.alive_max <= 2  # no record outlives the next fetch step
    d = rs._build(n, keys, seqs, 0, 500)
    assert d["sizes"]["U"]["+"][3] == 100_000.0


def test_missing_chromosome(tmp_path, capsys):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False), ("c1", 110, False)])
    out = rs.signals_for(str(bam), [Cluster("a", "chrZ", 0, 500),
                                    Cluster("b", "c1", 0, 500)], xm_as_nh=False)
    assert out["a"]["n_alns"] == 0 and len(out["a"]["logo"]["+"]) == 26
    assert not out["a"].get("failed")
    assert out["b"]["n_alns"] == 2
    assert "WARN" not in capsys.readouterr().err


def test_fetch_error_is_zero_filled_and_warned(tmp_path, capsys, monkeypatch):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False)])

    def boom(*a, **k):
        raise ValueError("bad record")
    monkeypatch.setattr(rs, "_accumulate", boom)
    out = rs.signals_for(str(bam), [Cluster("a", "c1", 0, 500),
                                    Cluster("b", "c2", 0, 500)], xm_as_nh=False)
    assert out["a"]["n_alns"] == 0 and out["a"]["failed"]
    assert len(out["a"]["logo"]["+"]) == 26
    err = capsys.readouterr().err
    assert "WARN" in err and "2 cluster(s)" in err and "ValueError: bad record" in err


def test_none_coordinate_is_zero_filled_and_warned(tmp_path, capsys):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False)])
    out = rs.signals_for(str(bam), [Cluster("a", "c1", None, 500),
                                    Cluster("b", "c1", 0, 500)], xm_as_nh=False)
    assert out["a"]["n_alns"] == 0 and out["a"].get("failed")
    assert out["b"]["n_alns"] == 1 and not out["b"].get("failed")
    assert "1 cluster(s)" in capsys.readouterr().err


def test_threads_argument_accepted(tmp_path):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False)])
    cl = [Cluster("b", "c1", 0, 500)]
    assert rs.signals_for(str(bam), cl, xm_as_nh=False, threads=4) == \
        rs.signals_for(str(bam), cl, xm_as_nh=False)


def _binned_reference(lane, B):
    return [round(sum(lane[i:i + B]) / len(lane[i:i + B]), 3) for i in range(0, len(lane), B)]


@pytest.mark.parametrize("start,end,max_bins", [
    (37, 1040, 100),   # 1003 bp: bin 11, partial last bin, clipped window
    (90, 213, 400),    # shorter than max_bins: bin 1
    (5, 1205, 7),      # few large bins
])
def test_coverage_bins_equal_per_base_reference(tmp_path, start, end, max_bins):
    rng = random.Random(start + end)
    reads = [("c1", rng.randint(0, 1250), _seq(rng, rng.randint(15, 40)), rng.choice("+-"),
              rng.choice([1, 3])) for _ in range(500)]
    reads += [("c1", start - 5, "A" * 30, "+", 1), ("c1", end - 3, "C" * 25, "-", 1)]
    alns = [a for a in _alns(reads) if a.pos < end and a.pos + len(a.seq) > start]
    d = _sig(tmp_path, reads, start, end, max_bins=max_bins)
    B = d["cov"]["bin"]
    assert B == -(-(end - start) // max_bins) and len(d["cov"]["+"]) <= max_bins
    plus, minus = class_signals.coverage_track(alns, start, end)
    assert d["cov"]["+"] == _binned_reference(plus, B)
    assert d["cov"]["-"] == _binned_reference(minus, B)
    p5, m5 = class_signals.five_prime_ends(alns)
    for s, ends in (("+", p5), ("-", m5)):
        ref = [0.0] * len(d["cov"][s])
        for p, w in ends.items():
            if start <= p < end:
                ref[(p - start) // B] += w
        assert d["ends"][s] == pytest.approx([round(v, 3) for v in ref], abs=1e-3)


def test_logo_and_seqs_equal_reference(tmp_path):
    rng = random.Random(7)
    pool = [_seq(rng, rng.choice([20, 21, 21, 22, 24])) for _ in range(40)]
    reads = [("c1", rng.randint(0, 800), rng.choice(pool), rng.choice("+-"),
              rng.choice([1, 2, 3, 5])) for _ in range(600)]
    alns = _alns(reads)
    d = _sig(tmp_path, reads, 0, 1000, top_seqs=15)
    # straightforward reference from class_signals._oriented
    w = {}
    for a in alns:
        k = (class_signals._oriented(a).replace("T", "U"), a.strand)
        w[k] = w.get(k, 0.0) + 1.0 / a.nh
    ref_rows = sorted(w.items(), key=lambda kv: (-round(kv[1], 9), kv[0][0], kv[0][1]))[:15]
    assert [r[:3] + [r[4]] for r in d["seqs"]] == \
        [[o, len(o), s, o[0]] for (o, s), _ in ref_rows]
    assert [r[3] for r in d["seqs"]] == pytest.approx([round(v, 3) for _, v in ref_rows],
                                                      abs=1e-3)
    dom = d["dominant"]
    by_len = {}
    for a in alns:
        by_len[len(a.seq)] = by_len.get(len(a.seq), 0.0) + 1.0 / a.nh
    assert dom == max(by_len, key=lambda k: (by_len[k], -k))
    for s in "+-":
        cnt = [dict.fromkeys("ACGU", 0.0) for _ in range(26)]
        for a in alns:
            if a.strand == s and len(a.seq) == dom:
                for i, b in enumerate(class_signals._oriented(a).replace("T", "U")[:26]):
                    cnt[i][b] += 1.0 / a.nh
        for i, p in enumerate(cnt):
            t = sum(p.values())
            exp = {b: (v / t if t else 0.0) for b, v in p.items()}
            assert d["logo"][s][i] == pytest.approx(exp, abs=1e-3)
