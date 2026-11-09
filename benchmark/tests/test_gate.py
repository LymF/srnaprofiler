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
                e_pingpong=0.9, e_phase=0.9, frac_top=0.5)
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
    ("arth", "piRNA", _rec("27", e_1u=0.2, e_pingpong=0.3)),  # low defining
    ("arth", "siRNA", _rec("21", e_fractop_bid=0.4)),
    ("arth", "siRNA", _rec("21", size_counts={21: 10, 27: 990})),  # size strength low
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
        assert load_config(f"{ROOT}/{f}")["decision"]["ml_gate_min_score"] == 0.5


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
