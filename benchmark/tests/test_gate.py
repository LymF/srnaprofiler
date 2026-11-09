"""Biology gate: spec v2 s2.1. A class is allowed only where the locus shows its biology."""
import copy

import pytest

import decider
import reporter
from decider import RuleScorer, load_config
from records import ClusterRecord

ROOT = decider.__file__.rsplit("/", 1)[0]


def _scorer(name):
    return RuleScorer(load_config(f"{ROOT}/{name}"))


@pytest.fixture(scope="module")
def plant():
    return _scorer("rules.plant.yaml")


@pytest.fixture(scope="module")
def arth():
    return _scorer("rules.yaml")


def _rec(call, **kw):
    counts = {int(call): 900, 21: 50, 22: 50} if call != "N" else {20: 10, 26: 10}
    base = dict(chrom="c", start=0, end=5000, name="L", reads=500.0,
                dicer_call=call, distinct_seqs=200, size_counts=counts,
                read_span=5000, density_pvalue=0.001, phase_pvalue=0.0001,
                phase_number=10, e_1u=0.9, e_5a=0.9, e_fractop_bid=0.9,
                e_pingpong=0.9, e_phase=0.9, frac_top=0.5,
                s_1u=5.0, s_5a=5.0, s_fractop_bid=5.0, s_pingpong=5.0, s_phase=5.0)
    base.update(kw)
    return ClusterRecord(**base)


def test_clean_records_pass(plant, arth):
    assert arth.passes_gate("piRNA", _rec("27"))
    assert plant.passes_gate("phasiRNA", _rec("21"))
    assert plant.passes_gate("phasiRNA24", _rec("24"))
    assert plant.passes_gate("hcsiRNA", _rec("24", phase_pvalue=1.0, phase_number=0))


@pytest.mark.parametrize("scorer,cls,rec", [
    ("arth", "piRNA", _rec("N")),
    ("arth", "piRNA", _rec("21")),                       # outside window
    ("arth", "piRNA", _rec("27", read_span=500)),        # short span
    ("arth", "piRNA", _rec("27", density_pvalue=0.5)),   # density p
    ("plant", "phasiRNA24", _rec("24", phase_pvalue=0.2)),  # phase p
    ("arth", "piRNA", _rec("27", s_1u=1.0, s_pingpong=1.9)),  # neither OR member p <= 0.01
    ("arth", "siRNA", _rec("21", s_fractop_bid=1.5)),
    ("plant", "hcsiRNA", _rec("24", phase_pvalue=1.0, s_5a=1.0)),
])
def test_gate_rejects(request, scorer, cls, rec):
    assert not request.getfixturevalue(scorer).passes_gate(cls, rec)


def test_gate_rejects_disabled_or_unknown_class(arth):
    assert not arth.passes_gate("hcsiRNA", _rec("24"))
    assert not arth.passes_gate("other", _rec("27"))


def test_gate_does_not_mutate_record(arth):
    r = _rec("27")
    before = copy.deepcopy(r)
    arth.passes_gate("piRNA", r)
    assert r == before


def test_gated_choice_falls_through(plant):
    r = _rec("24", phase_pvalue=1.0, phase_number=0)
    P = {"phasiRNA24": 0.6, "hcsiRNA": 0.3, "siRNA": 0.1}
    assert plant.gated_choice(P, r) == "hcsiRNA"


def test_gated_choice_other_on_top_is_ambiguous(plant):
    r = _rec("24")
    assert plant.gated_choice({"other": 0.5, "phasiRNA24": 0.4, "hcsiRNA": 0.1}, r) \
        == RuleScorer.AMBIGUOUS


def test_gated_choice_skips_mirna_and_ambiguous_when_none_pass(arth):
    r = _rec("N")
    assert arth.gated_choice({"miRNA": 0.6, "piRNA": 0.3, "siRNA": 0.1}, r) \
        == RuleScorer.AMBIGUOUS
    assert arth.gated_choice({"miRNA": 0.6, "piRNA": 0.4}, _rec("27")) == "piRNA"


def test_gate_threshold_in_both_yaml():
    for f in ("rules.yaml", "rules.plant.yaml"):
        assert load_config(f"{ROOT}/{f}")["decision"]["ml_gate_max_pvalue"] == 0.01


def test_rules_decide_unchanged(arth):
    r = arth.decide(_rec("27"))
    assert r.cls == "piRNA"
    # gate helpers leave the rules path alone
    arth.gated_choice({"piRNA": 1.0}, _rec("27"))
    assert arth.decide(_rec("27")).cls == r.cls


def test_reporter_writes_evidence_columns(tmp_path):
    cols = reporter._COLS
    i = cols.index("s_phase")
    assert cols[i + 1:i + 6] == ["e_1u", "e_5a", "e_fractop_bid", "e_pingpong", "e_phase"]
    p = tmp_path / "c.tsv"
    reporter.write_tsv([_rec("27", e_1u=0.25)], str(p))
    head, row = [l.split("\t") for l in p.read_text().splitlines()]
    assert row[head.index("e_1u")] == "0.25"


def test_nan_defining_evidence_fails_the_gate(arth, plant):
    nan = float("nan")
    assert not arth.passes_gate("piRNA", _rec("27", s_1u=nan, s_pingpong=nan))
    assert not arth.passes_gate("siRNA", _rec("21", s_fractop_bid=nan))
    assert not plant.passes_gate("phasiRNA", _rec("21", s_phase=nan))
    # one significant member of an OR slot is enough, even if the other is NaN
    assert arth.passes_gate("piRNA", _rec("27", s_1u=nan, s_pingpong=3.0))


def test_siRNA_passes_with_perfect_duplex_strength(arth):
    # R27: e_fractop_bid cannot exceed 0.44, so a strength floor of 0.5 blocked every siRNA
    assert arth.passes_gate("siRNA", _rec("21", e_fractop_bid=0.44, s_fractop_bid=10.0))


def test_gate_flags_cover_enabled_classes(arth, plant):
    f = arth.gate_flags(_rec("27"))
    assert set(f) == set(RuleScorer.enabled_classes(arth.cfg)) and f["piRNA"] is True
    assert "hcsiRNA" not in f and "hcsiRNA" in plant.gate_flags(_rec("24"))


def test_reporter_writes_gate_columns(tmp_path):
    r = _rec("27")
    r.gates = {"piRNA": True, "siRNA": False}
    p = tmp_path / "c.tsv"
    reporter.write_tsv([r], str(p), ["piRNA", "siRNA"])
    head, row = [l.split("\t") for l in p.read_text().splitlines()]
    assert head[-2:] == ["gate_piRNA", "gate_siRNA"] and row[-2:] == ["1", "0"]
    assert head[:-2] == reporter._COLS
    reporter.write_tsv([r], str(p))                       # no classes: columns unchanged
    assert p.read_text().splitlines()[0].split("\t") == reporter._COLS


# --- round trip: records -> gate columns in clusters.tsv -> benchmark decide == production

def _real_records(n=60, seed=3):
    import random
    from collections import namedtuple
    import features
    aln = namedtuple("A", "chrom pos seq nh strand")
    rnd = random.Random(seed)
    cfg = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}
    out = []
    for i in range(n):
        lens = rnd.choice([(21,), (24,), (24, 21), (27, 26), (22, 23, 24, 25)])
        first = rnd.choice("TTTA")       # 1U- or 5'A-biased loci, so the gate passes sometimes
        alns = []
        for k in range(rnd.randint(30, 120)):
            L = rnd.choice(lens)
            seq = first + "".join(rnd.choice("ACGT") for _ in range(L - 1))
            alns.append(aln("c", 10000 * i + rnd.randint(0, 3000), seq,
                            rnd.choice([1, 1, 2, 3]), rnd.choice("+-")))
        r = features.build_record(("c", 10000 * i, 10000 * i + 3500), alns, 1e6, cfg)
        # size_counts and dicer_call are the real ones; the evidence that would need
        # real ping-pong/phasing reads is drawn, so the gate passes and fails by class
        for e in ("1u", "5a", "fractop_bid", "pingpong", "phase"):
            setattr(r, f"e_{e}", rnd.choice([0.2, 0.7, 0.9]))
            setattr(r, f"s_{e}", rnd.choice([0.5, 1.9, 2.0, 6.0]))
        r.density_pvalue = rnd.choice([0.001, 0.5])
        r.phase_pvalue = rnd.choice([0.001, 0.4])
        r.phase_number = rnd.randint(0, 10)
        out.append(r)
    return out


class _Fixed:
    def __init__(self, P, classes):
        import numpy as np
        self.P, self.classes = np.asarray(P, float), classes

    def predict_proba(self, table):
        import pandas as pd
        return pd.DataFrame(self.P, columns=self.classes, index=table.index)


@pytest.mark.parametrize("name", ["plant", "arth"])
def test_benchmark_decide_matches_production_on_real_rows(request, tmp_path, name):
    import numpy as np
    import pandas as pd
    from benchmark.metrics_class import read_loci
    from benchmark.ml import evaluate
    scorer = request.getfixturevalue(name)
    recs = _real_records()
    for r in recs:
        scorer.decide(r)
        r.gates = scorer.gate_flags(r)
    classes = list(RuleScorer.enabled_classes(scorer.cfg))
    assert any(any(r.gates.values()) for r in recs) and not all(all(r.gates.values()) for r in recs)
    path = tmp_path / "c.clusters.tsv"
    reporter.write_tsv(recs, str(path), classes)
    loci = read_loci(path)
    table = pd.DataFrame([{**l, "carved": False,
                           **{c: float(l[c]) for c in l if c.startswith("gate_")}}
                          for l in loci])
    cols = ["other", *[c for c in classes if c != "miRNA"]]
    P = np.random.default_rng(0).dirichlet(np.ones(len(cols)), len(recs))
    got = evaluate.decide(_Fixed(P, cols), table, [], 0.0).tolist()
    want = [scorer.gated_choice(dict(zip(cols, p)), r) for p, r in zip(P, recs)]
    assert got == want
    assert "other" not in got and "miRNA" not in got
    for r, g in zip(recs, table.to_dict("records")):   # the columns are the gate itself
        assert all(bool(g[f"gate_{c}"]) == scorer.passes_gate(c, r) for c in classes)


def test_unphased_24nt_phasirna_falls_through_to_hcsirna_in_benchmark(plant):
    import pandas as pd
    from benchmark.ml import evaluate
    r = _rec("24", phase_pvalue=1.0, phase_number=0)
    flags = plant.gate_flags(r)
    table = pd.DataFrame([{"carved": False, **{f"gate_{k}": float(v) for k, v in flags.items()}}])
    m = _Fixed([[0.6, 0.3, 0.1]], ["phasiRNA24", "hcsiRNA", "siRNA"])
    assert evaluate.decide(m, table, [], 0.0).tolist() == ["hcsiRNA"]
