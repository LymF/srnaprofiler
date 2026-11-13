"""miRNA hairpin data of the report: matching, segments, reads, fold, payload."""
import os
import random
import stat
import sys
from collections import namedtuple

import pysam
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))
import bamutil
import report_mirna as rm
import report_signals as rs
from features import revcomp
from report_data import Cluster  # noqa: F401  (documents the real type)

Cl = namedtuple("Cl", "name chrom start end cls")

COLS = ["name", "chrom", "start", "end", "strand", "engine", "mature_seq", "mature_len",
        "star_seq", "star_len", "precursor_seq", "precursor_len", "mature_start", "mature_end",
        "star_start", "star_end", "result", "flags", "mismatches", "bulges", "mature_reads",
        "star_reads", "total_reads", "precision", "mirdeep_score", "mirdeep_prob", "rfam_alert",
        "randfold", "mirbase_hit", "mirbase_same_seed"]

# hairpin 200..260 (0-based), flank 25: precursor window 175..285, 110 nt
START, END, LO, HI = 200, 260, 175, 285
M_AB, S_AB = (28, 50), (75, 97)  # 0-based half-open in the precursor


def precursor(seed=3):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGU") for _ in range(HI - LO))


def row_dict(strand="+", prec=None, **kw):
    prec = prec or precursor()
    d = dict.fromkeys(COLS, "NA")
    d.update(name="mir1", chrom="c1", start=START, end=END, strand=strand, engine="mirdeep2",
             mature_seq=prec[M_AB[0]:M_AB[1]].replace("U", "T"), mature_len=22,
             star_seq=prec[S_AB[0]:S_AB[1]].replace("U", "T"), star_len=22,
             precursor_seq=prec.replace("U", "T"), precursor_len=len(prec),
             mature_start=M_AB[0] + 1, mature_end=M_AB[1], star_start=S_AB[0] + 1, star_end=S_AB[1],
             result="Pass", flags="NA", mismatches=2, bulges=1, mature_reads=120, star_reads=15,
             total_reads=140, precision=96.4, mirdeep_score="12.5", mirdeep_prob="82%",
             rfam_alert="-", randfold="0.01", mirbase_hit="ath-miR166a", mirbase_same_seed="NA")
    d.update(kw)
    return d


def write_tsv(path, rows):
    path.write_text("\n".join("\t".join(str(r[c]) for c in COLS) for r in [dict(zip(COLS, COLS))] + rows) + "\n")


def load_one(tmp_path, **kw):
    p = tmp_path / "x.mirna.tsv"
    write_tsv(p, [row_dict(**kw)])
    return rm.load(str(p))[0]


# ---- load and match ----

def test_load_parses_and_maps_na(tmp_path):
    r = load_one(tmp_path, star_seq="NA", star_start="NA", star_end="NA", precision="NA")
    assert r.name == "mir1" and (r.start, r.end) == (START, END) and r.strand == "+"
    assert r.star_seq is None and r.star_start is None and r.precision is None
    assert "T" not in r.precursor_seq and r.mature_reads == 120.0 and r.mismatches == 2
    assert r.mirdeep_prob == "82%" and r.mirbase_same_seed is None


def test_load_without_locus_columns_is_an_error(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("a\tb\n1\t2\n")
    with pytest.raises(ValueError, match="lacks the column"):
        rm.load(str(p))


def mkrow(name, chrom, start, end, result="Pass"):
    return rm.Row(name=name, chrom=chrom, start=start, end=end, strand="+", result=result)


def test_match_overlap_rule_chrom_and_class():
    rows = [mkrow("a", "c1", 100, 200), mkrow("b", "c2", 100, 200)]
    cl = [Cl("hit", "c1", 110, 190, "miRNA"),        # inside
          Cl("half", "c1", 150, 250, "miRNA"),        # exactly 50% of its span covered
          Cl("less", "c1", 160, 260, "miRNA"),        # 40%
          Cl("chrom", "c3", 110, 190, "miRNA"),
          Cl("pi", "c1", 110, 190, "piRNA"),          # not a miRNA cluster
          Cl("none", None, None, None, "miRNA")]
    m = rm.match(cl, rows)
    assert {k: v.name for k, v in m.items()} == {"hit": "a", "half": "a"}


def test_match_prefers_largest_overlap_then_pass():
    rows = [mkrow("small", "c1", 100, 160), mkrow("big", "c1", 100, 190),
            mkrow("failed", "c1", 100, 190, "Fail")]
    c = Cl("c", "c1", 100, 200, "miRNA")
    assert rm.match([c], rows)["c"].name == "big"        # larger overlap; Pass before Fail on a tie
    assert rm.match([c], rows[::-1])["c"].name == "big"


# ---- geometry ----

def test_window_flanks_and_clipping():
    r = rm.Row("m", "c1", START, END, "+", precursor_seq="A" * 110, precursor_len=110)
    assert rm.window(r) == (LO, HI)
    # hairpin near position 0: left flank clipped to start
    r = rm.Row("m", "c1", 10, 70, "+", precursor_seq="A" * 95, precursor_len=95)
    assert rm.window(r) == (0, 95)
    # right flank clipped (chromosome end)
    r = rm.Row("m", "c1", 200, 260, "+", precursor_seq="A" * 100, precursor_len=100)
    assert rm.window(r) == (175, 275)
    # impossible lengths or no sequence: no window
    assert rm.window(rm.Row("m", "c1", 200, 260, "+", precursor_seq="A" * 200, precursor_len=200)) is None
    assert rm.window(rm.Row("m", "c1", 200, 260, "+", precursor_seq="A" * 40, precursor_len=40)) is None
    assert rm.window(rm.Row("m", "c1", 200, 260, "+")) is None
    assert rm.window(rm.Row("m", "c1", 200, 260, "+", precursor_seq="A" * 100, precursor_len=99)) is None


def test_to_precursor_plus_and_minus():
    assert rm.to_precursor(LO + 28, LO + 50, "+", LO, HI) == (28, 50)
    assert rm.to_precursor(HI - 50, HI - 28, "-", LO, HI) == (28, 50)
    assert rm.to_precursor(LO, LO + 1, "-", LO, HI) == (HI - LO - 1, HI - LO)   # lowest coordinate is the last base


@pytest.mark.parametrize("strand", "+-")
def test_segments_match_the_sequences(tmp_path, strand):
    r = load_one(tmp_path, strand=strand)
    seg = rm.segments(r)
    assert seg["m"] == M_AB and seg["s"] == S_AB and seg["loop"] == (M_AB[1], S_AB[0])
    p = r.precursor_seq
    assert p[slice(*seg["m"])] == r.mature_seq and p[slice(*seg["s"])] == r.star_seq


def test_segments_star_first_without_star_and_zero_based_start(tmp_path):
    # 3p mature: the star comes first, the loop is between
    prec = precursor()
    r = load_one(tmp_path, mature_seq=prec[75:97], mature_start=76, mature_end=97,
                 star_seq=prec[28:50], star_start=29, star_end=50)
    seg = rm.segments(r)
    assert seg["m"] == (75, 97) and seg["s"] == (28, 50) and seg["loop"] == (50, 75)
    # no star: only the mature
    r = load_one(tmp_path, star_seq="NA", star_start="NA", star_end="NA")
    assert rm.segments(r) == {"m": M_AB, "s": None, "loop": None}
    # a 0-based start (miRScore's failure rows) and a wrong start are both corrected by the sequence
    assert rm.segments(load_one(tmp_path, mature_start=M_AB[0], mature_end=M_AB[1] + 0))["m"] == M_AB
    assert rm.segments(load_one(tmp_path, mature_start=3, mature_end=25))["m"] == M_AB


# ---- reads inside the precursor, on a tiny BAM ----

def plus_reads(prec):
    dna = prec.replace("U", "T")
    mat = dna[M_AB[0]:M_AB[1]]
    mm = mat[:5] + ("A" if mat[5] != "A" else "C") + mat[6:]      # one mismatch at read position 5
    star = dna[S_AB[0]:S_AB[1]]
    g = lambda a: LO + a
    return [("c1", g(28), mat, "+", 1), ("c1", g(28), mat, "+", 1), ("c1", g(28), mat, "+", 2),
            ("c1", g(28), mm, "+", 1),
            ("c1", g(75), star, "+", 1),
            ("c1", g(52), dna[52:72], "+", 1),                       # in the loop: kind ""
            ("c1", g(28), revcomp(mat), "-", 1),                     # antisense, not shown
            ("c1", LO - 5, dna[:10] + "AAAAA", "+", 1)]              # overhangs the window start


def hp_of(tmp_path, reads, strand="+", prec=None):
    row = rm.Row("m", "c1", START, END, strand, mature_seq=None, precursor_seq=prec or precursor(),
                 precursor_len=HI - LO)
    spec = (LO, HI, strand, row.precursor_seq, M_AB, S_AB)
    bam = tmp_path / "h.bam"
    bamutil.make_seq_bam(bam, reads)
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        hp = rs.HairpinReads(spec)
        rs._accumulate(fh, "c1", START, END, xm_as_nh=False, hp=hp)
    return hp.build()


def test_reads_collapsed_weighted_and_mismatches_marked(tmp_path):
    prec = precursor()
    out = hp_of(tmp_path, plus_reads(prec), prec=prec)
    rows = {(a, b, q): (w, k) for a, b, q, w, k in out["rd"]}
    mat = prec[M_AB[0]:M_AB[1]].lower()
    assert rows[M_AB[0], M_AB[1], mat] == (2.5, "m")             # 2 reads + one NH=2 read (0.5), collapsed
    mm = [r for r in out["rd"] if r[2] != mat and r[0] == M_AB[0]]
    assert len(mm) == 1 and mm[0][3] == 1.0 and mm[0][4] == "m"
    assert sum(1 for c in mm[0][2] if c != c.lower()) == 1 and mm[0][2][5].isupper()
    assert mm[0][2].lower() == mat[:5] + mm[0][2][5].lower() + mat[6:]
    star = [r for r in out["rd"] if r[0] == S_AB[0]]
    assert star[0][4] == "s" and star[0][2] == prec[S_AB[0]:S_AB[1]].lower()
    other = [r for r in out["rd"] if r[0] == 52]
    assert other[0][4] == ""
    assert out["anti"] == 1.0 and out["ovh"] == 1.0              # opposite strand; past the window start
    assert out["more"] == 0 and out["mr"] == 0
    assert len(out["cov"]) == HI - LO
    assert out["cov"][M_AB[0]] == pytest.approx(3.5)             # 2.5 + the mismatched read (1.0)
    assert out["rd"][0][3] >= out["rd"][-1][3]                   # most abundant first


def test_coverage_is_the_sum_of_collapsed_reads(tmp_path):
    prec = precursor()
    out = hp_of(tmp_path, plus_reads(prec), prec=prec)
    total = sum(r[3] * (r[1] - r[0]) for r in out["rd"])
    assert sum(out["cov"]) == pytest.approx(total, abs=0.01)


def test_minus_strand_reads_are_reverse_complemented_and_flipped(tmp_path):
    prec = precursor()
    dna = prec.replace("U", "T")
    # on '-' the precursor's position 0 is genomic HI-1: read [a,b) lies at genomic [HI-b, HI-a)
    reads = [("c1", HI - M_AB[1], revcomp(dna[M_AB[0]:M_AB[1]]), "-", 1),
             ("c1", HI - M_AB[1], revcomp(dna[M_AB[0]:M_AB[1]]), "-", 1),
             ("c1", HI - S_AB[1], revcomp(dna[S_AB[0]:S_AB[1]]), "-", 1),
             ("c1", HI - M_AB[1], dna[M_AB[0]:M_AB[1]], "+", 1)]       # antisense for a '-' precursor
    out = hp_of(tmp_path, reads, strand="-", prec=prec)
    assert out["rd"][0] == [M_AB[0], M_AB[1], prec[M_AB[0]:M_AB[1]].lower(), 2.0, "m"]
    assert out["rd"][1][0:2] == list(S_AB) and out["rd"][1][4] == "s"
    assert out["anti"] == 1.0
    # mapping a precursor interval back to the genome is the inverse
    assert rm.to_precursor(HI - M_AB[1], HI - M_AB[0], "-", LO, HI) == M_AB


def test_row_cap_and_more(tmp_path):
    prec = precursor()
    dna = prec.replace("U", "T")
    # 30 distinct 20-mers, i+1 reads each
    reads = [("c1", LO + 10 + i, dna[10 + i:30 + i], "+", 1) for i in range(30) for _ in range(i + 1)]
    out = hp_of(tmp_path, reads, prec=prec)
    assert len(out["rd"]) == rs.READ_ROWS == 25
    assert out["more"] == 5 and out["mr"] == sum(range(1, 6))    # the five least abundant: 1..5 reads
    assert out["rd"][0][3] == 30.0 and out["rd"][-1][3] == 6.0


def test_cluster_figures_unchanged_by_the_precursor_pass(tmp_path):
    prec = precursor()
    bam = tmp_path / "h.bam"
    bamutil.make_seq_bam(bam, plus_reads(prec))
    C = namedtuple("C", "name chrom start end")
    cl = [C("a", "c1", START, END)]
    spec = (LO, HI, "+", prec, M_AB, S_AB)
    plain = rs.signals_for(str(bam), cl, xm_as_nh=False)
    with_hp = rs.signals_for(str(bam), cl, xm_as_nh=False, hairpins={"a": spec})
    hp = with_hp["a"].pop("hp")
    assert with_hp == plain and hp["rd"]
    # the precursor window is read even where the cluster window holds no read at all
    far = rs.signals_for(str(bam), [C("a", "c1", 3000, 3010)], xm_as_nh=False, hairpins={"a": spec})
    assert far["a"]["n_alns"] == 0 and far["a"]["hp"]["rd"]
    assert rs.signals_for(str(bam), cl, xm_as_nh=False, threads=2, hairpins={"a": spec})["a"]["hp"] == hp


def test_hairpin_failure_is_zero_filled(tmp_path):
    bam = tmp_path / "h.bam"
    bamutil.make_seq_bam(bam, plus_reads(precursor()))
    C = namedtuple("C", "name chrom start end")
    out = rs.signals_for(str(bam), [C("a", "nochrom", START, END)], xm_as_nh=False,
                         hairpins={"a": (LO, HI, "+", precursor(), M_AB, S_AB)})
    assert out["a"]["n_alns"] == 0 and "hp" not in out["a"]


# ---- RNAfold ----

def test_fold_absent_path(monkeypatch):
    monkeypatch.setattr(rm.shutil, "which", lambda name: None)
    monkeypatch.setattr(rm, "_import_rna", lambda: None)
    assert rm.fold("GGGGAAAACCCC") == (None, None)
    assert rm.fold_tool_missing()


def _fake_rnafold(tmp_path, body):
    exe = tmp_path / "RNAfold"
    exe.write_text("#!/bin/sh\n" + body)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return str(exe)


def test_fold_parses_rnafold_output_and_writes_nothing(tmp_path, monkeypatch):
    seq = "GGGGAAAACCCC"
    exe = _fake_rnafold(tmp_path, f'read s\necho "$s"\necho "((((....)))) ( -4.30)"\n')
    monkeypatch.setattr(rm.shutil, "which", lambda name: exe)
    cwd = tmp_path / "w"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    assert rm.fold(seq) == ("((((....))))", -4.3)
    assert not rm.fold_tool_missing()
    assert list(cwd.iterdir()) == []


def test_fold_failure_is_none(tmp_path, monkeypatch):
    monkeypatch.setattr(rm.shutil, "which", lambda name: _fake_rnafold(tmp_path, "exit 1\n"))
    assert rm.fold("GGGGAAAACCCC") == (None, None)
    monkeypatch.setattr(rm.shutil, "which", lambda name: _fake_rnafold(tmp_path, 'echo x\necho "(((" \n'))
    assert rm.fold("GGGGAAAACCCC") == (None, None)             # wrong length


def test_fold_uses_the_python_module_when_there_is_no_executable(monkeypatch):
    class FakeRNA:
        @staticmethod
        def fold(s):
            return "(" * 4 + "." * (len(s) - 8) + ")" * 4, -3.0
    monkeypatch.setattr(rm.shutil, "which", lambda name: None)
    monkeypatch.setattr(rm, "_import_rna", lambda: FakeRNA)
    db, mfe = rm.fold("GGGGAAAACCCC")
    assert db == "((((....))))" and mfe == -3.0 and not rm.fold_tool_missing()


@pytest.mark.skipif(rm.fold_tool_missing(), reason="RNAfold (ViennaRNA) not installed")
def test_real_rnafold_folds_a_hairpin():
    db, mfe = rm.fold("GGGGGAAAAAAACCCCC")
    assert db.count("(") == db.count(")") and len(db) == 17 and mfe <= 0


# ---- payload ----

def test_criteria_statuses_and_wording(tmp_path):
    r = load_one(tmp_path)
    c = {k: (st, obs) for k, st, obs, rule in rm.criteria(r)}
    assert c == {"len": (1, "60 nt"), "mm": (1, "2"), "prec": (1, "96.4%")}
    r = load_one(tmp_path, mismatches=6, precision=74.9)
    r.end = END + 300                                           # a 360 nt hairpin locus
    c = {k: st for k, st, obs, rule in rm.criteria(r)}
    assert c == {"len": 0, "mm": 0, "prec": 0}
    r = load_one(tmp_path, mismatches="NA", precision="NA")
    assert {k: st for k, st, _, _ in rm.criteria(r)}["mm"] == 2
    assert {k: st for k, st, _, _ in rm.criteria(r)}["prec"] == 2
    boundary = load_one(tmp_path, mismatches=5, precision=75)
    assert [st for _, st, _, _ in rm.criteria(boundary)] == [1, 1, 1]


def test_numbers_cover_the_requested_fields(tmp_path):
    keys = [k for k, _ in rm.numbers(load_one(tmp_path))]
    for k in ("mature_reads", "star_reads", "precision", "mismatches", "bulges", "precursor_len",
              "mirdeep_score", "mirdeep_prob", "randfold", "rfam_alert", "mirbase_hit",
              "mirbase_same_seed"):
        assert k in keys


def test_figure_data_shapes(tmp_path):
    r = load_one(tmp_path)
    hp = {"rd": [[28, 50, "a" * 22, 2.0, "m"]], "cov": [0.0] * 110, "more": 0, "mr": 0.0,
          "anti": 0.0, "ovh": 0.0}
    d = rm.figure_data(r, ("((..))", -2.0, False), hp)
    assert d["L"] == 110 and d["m"] == [28, 50] and d["s"] == [75, 97] and d["lp"] == [50, 75]
    assert d["db"] == "((..))" and d["mfe"] == -2.0 and "ns" not in d
    assert d["lo"] == LO and d["rd"] == hp["rd"] and d["st"] == "+"
    # no structure: the reason is carried, and no db
    d = rm.figure_data(r, (None, None, True), hp)
    assert "db" not in d and d["ns"] == "m"
    assert rm.figure_data(r, (None, None, False), None)["ns"] == "f"
    # without a precursor the numbers still come
    d = rm.figure_data(load_one(tmp_path, precursor_seq="NA", precursor_len="NA"), None, None)
    assert "seq" not in d and d["t"] and d["pc"]
    # a window that cannot be placed keeps the sequence, drops the reads
    d = rm.figure_data(load_one(tmp_path, end=END + 100), (None, None, False), hp)
    assert d["nowin"] == 1 and "rd" not in d


def test_figure_data_failed_read_and_unknown_strand(tmp_path):
    r = load_one(tmp_path)
    assert rm.figure_data(r, (None, None, False), None, failed=True)["err"] == 1
    assert "err" not in rm.figure_data(r, (None, None, False), {"rd": [], "cov": []})
    u = load_one(tmp_path, strand=".")
    assert rm.hairpin_spec(u) is None
    d = rm.figure_data(u, (None, None, False), None, failed=True)
    assert d["nostrand"] == 1 and "err" not in d and "rd" not in d and d["st"] == "+"
    assert rm.hairpin_spec(r)[2] == "+"


def test_f_returns_none_for_non_finite():
    for bad in (float("inf"), float("-inf"), float("nan"), "inf", "nan", "x", None):
        assert rm._f(bad) is None
    assert rm._f("1.5") == 1.5 and rm._f(0) == 0.0


# ---- 2D layout ----

PS = """%!PS
/coor [
[92.50000000 92.50000000]
[92.50000000 77.50000000]
[92.50000000 62.50000000]
[80.40541840 22.58762932]
] def
/arcs [
[]
] def
"""


def test_parse_coor():
    assert rm.parse_coor(PS) == [(92.5, 92.5), (92.5, 77.5), (92.5, 62.5), (80.40541840, 22.58762932)]
    assert rm.parse_coor("no coordinates") == [] and rm.parse_coor(None) == []


def test_quantise_orients_five_prime_top_left_and_scales():
    # a straight 5-base line drawn upward-right in tool coordinates (y up), 20 units apart
    xy = rm.quantise_xy([(0, 0), (20, 0), (40, 0), (60, 0), (80, 0)])
    assert len(xy) == 10 and all(isinstance(v, int) for v in xy)
    xs, ys = xy[0::2], xy[1::2]
    assert min(xs) == 0 and min(ys) == 0 and (xs[0], ys[0]) == (0, 0)
    assert max(xs) - min(xs) == 4 * rm.UNIT or max(ys) - min(ys) == 4 * rm.UNIT
    # a hairpin whose 5' end is at the bottom right is rotated so it ends up near the top-left
    pts = [(100, 0), (100, 10), (100, 20), (95, 30), (90, 30), (85, 20), (85, 10), (85, 0)]
    q = rm.quantise_xy(pts)
    x0, y0 = q[0], q[1]
    assert x0 + y0 == min(a + b for a, b in zip(q[0::2], q[1::2])) or x0 + y0 <= 20
    assert rm.quantise_xy([(0, 0)]) == [] and rm.quantise_xy([(1, 1), (1, 1)]) == []


class _P:
    def __init__(self, x, y):
        self.X, self.Y = x, y


def test_layout_uses_the_module_and_drops_its_extra_point(monkeypatch):
    class Mod:
        @staticmethod
        def naview_xy_coordinates(db):
            return [_P(i * 15.0, 0.0) for i in range(len(db) + 1)]
    monkeypatch.setattr(rm, "_import_rna", lambda: Mod)
    xy = rm.layout("((..))")
    assert len(xy) == 12
    # naview failing falls back to the simple layout
    class Mod2:
        @staticmethod
        def naview_xy_coordinates(db):
            raise RuntimeError("x")

        @staticmethod
        def simple_xy_coordinates(db):
            return [_P(i * 15.0, 0.0) for i in range(len(db) + 1)]
    monkeypatch.setattr(rm, "_import_rna", lambda: Mod2)
    assert len(rm.layout("((..))")) == 12


def test_layout_falls_back_to_rnaplot_and_parses_its_postscript(tmp_path, monkeypatch):
    exe = tmp_path / "RNAplot"
    body = "\n".join(f"[{i * 15}.0 0.0]" for i in range(6))
    exe.write_text(f"#!/bin/sh\ncat > /dev/null\nprintf '/coor [\\n{body}\\n] def\\n' > rna.ps\n")
    exe.chmod(0o755)
    monkeypatch.setattr(rm, "_import_rna", lambda: None)
    monkeypatch.setattr(rm.shutil, "which", lambda name: str(exe) if name == "RNAplot" else None)
    assert len(rm.layout("((..))")) == 12
    exe.write_text("#!/bin/sh\nexit 1\n")
    assert rm.layout("((..))") == []


def test_layout_without_a_tool_is_empty(monkeypatch):
    monkeypatch.setattr(rm, "_import_rna", lambda: None)
    monkeypatch.setattr(rm.shutil, "which", lambda name: None)
    assert rm.layout("((..))") == [] and rm.layout(None) == []


@pytest.mark.skipif(rm._import_rna() is None and not rm.shutil.which("RNAplot"),
                    reason="ViennaRNA layout not installed")
def test_real_layout():
    db = "(((((.......)))))"
    xy = rm.layout(db)
    assert len(xy) == 2 * len(db)


def test_xy_in_the_payload_is_small_and_checked(tmp_path):
    import json
    import random
    r = load_one(tmp_path)
    rng = random.Random(1)
    xy = []
    for i in range(110):                                     # about the extent of a 110-nt hairpin
        xy += [rng.randint(0, 400), rng.randint(0, 700)]
    db = "(" * 40 + "." * 30 + ")" * 40
    d = rm.figure_data(r, (db, -30.0, False, xy), None)
    assert d["xy"] == xy and len(json.dumps(d["xy"], separators=(",", ":"))) <= 4096
    assert "xy" not in rm.figure_data(r, (db, -30.0, False, xy[:-2]), None)   # wrong length is dropped
    assert "xy" not in rm.figure_data(r, (db, -30.0, False), None)
    assert "xy" not in rm.figure_data(r, (None, None, True, xy), None)


class _Rec:
    def __init__(self, start, length, reverse):
        self.reference_start, self.reference_end = start, start + length
        self.is_reverse = reverse
        self.cigartuples = [(0, length)]
        self.query_sequence = "A" * length


def test_reads_outside_the_precursor_window_are_not_counted():
    hp = rs.HairpinReads((1000, 1100, "+", "A" * 100, None, None))
    hp.add(_Rec(1500, 22, False), 1)      # sense, elsewhere in the cluster
    hp.add(_Rec(1600, 22, True), 1)       # antisense, elsewhere in the cluster
    assert hp.anti == 0 and hp.ovh == 0 and not hp.rows
    hp.add(_Rec(1050, 22, True), 1)       # antisense inside the window still counts
    hp.add(_Rec(990, 22, False), 1)       # sense overhanging the window start still counts
    assert hp.anti == 1 and hp.ovh == 1
