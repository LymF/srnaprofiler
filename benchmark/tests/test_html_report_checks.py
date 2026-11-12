"""Class checklist ("Why this class"): rows, thresholds and statuses mirror the decider."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import report_checks as rc
import report_data
import report_text as rt
from test_html_report_html import mk

PASS, FAIL, NA, UNTESTED = rc.PASS, rc.FAIL, rc.NA, rc.UNTESTED


def statuses(ck):
    return [r[0] for r in ck[1:]]


def test_pirna_rows_follow_the_rules():
    plan = rc.skeleton("animal")["piRNA"]
    assert [r[0] for r in plan] == ["check.size.title", "check.1u_pingpong.title",
                                    "check.length.title", "check.density.title"]
    assert plan[0][1] == "Dominant size of 24 to 31 nt"
    assert "1000" in plan[2][1] and "0.01" in plan[3][1]


def test_pirna_all_pass():
    c = mk(cls="piRNA", dicer_call=27, pct_1u=0.8, start=0, end=5000,
           e={"1u": 0.6, "pingpong": 0.0, "5a": 0.0, "fractop_bid": 0.0, "phase": 0.0})
    ck = rc.checklist(c, "animal")
    assert ck[0] == "piRNA"
    assert statuses(ck) == [PASS, PASS, PASS, PASS]
    assert ck[1][1] == "27 nt"
    assert "80%" in ck[2][1] and "0.60" in ck[2][1]


def test_pirna_or_slot_either_member_passes_neither_fails():
    base = dict(cls="piRNA", dicer_call=27, start=0, end=5000)
    only_pp = mk(e={"1u": 0.0, "pingpong": 0.4, "5a": 0, "fractop_bid": 0, "phase": 0}, **base)
    neither = mk(e={"1u": 0.0, "pingpong": 0.0, "5a": 0, "fractop_bid": 0, "phase": 0}, **base)
    assert rc.checklist(only_pp, "animal")[2][0] == PASS
    assert rc.checklist(neither, "animal")[2][0] == FAIL


def test_old_tsv_evidence_is_not_available():
    c = mk(cls="piRNA", dicer_call=27, start=0, end=5000,
           e={k: None for k in ("1u", "5a", "fractop_bid", "pingpong", "phase")})
    ck = rc.checklist(c, "animal")
    assert ck[2][0] == NA
    assert "80%" in ck[2][1]           # the raw column still shows


def test_sirna_strand_and_size():
    c = mk(cls="siRNA", dicer_call=21, frac_top=0.62,
           e={"1u": 0, "5a": 0, "fractop_bid": 0.45, "pingpong": 0, "phase": 0})
    ck = rc.checklist(c, "animal")
    plan = rc.skeleton("animal")["siRNA"]
    assert [r[0] for r in plan] == ["check.size.title", "check.strand.title"]
    assert statuses(ck) == [PASS, PASS]
    assert "38%" in ck[2][1]
    assert "10%" in plan[1][1]


def test_mirna_and_other_have_no_checklist():
    assert rc.checklist(mk(cls="miRNA"), "animal") is None
    assert rc.checklist(mk(cls="other"), "animal") is None
    assert "miRNA" not in rc.skeleton("animal")


def test_ambiguous_size_outside_window_fails_size_and_skips_criteria():
    c = mk(cls="ambiguous", dicer_call=19, start=0, end=200,
           evidence="nearest=piRNA; size=0.3",
           e={"1u": 0.0, "pingpong": 0.0, "5a": 0, "fractop_bid": 0, "phase": 0})
    ck = rc.checklist(c, "animal")
    assert ck[0] == "piRNA"
    assert statuses(ck) == [FAIL, FAIL, UNTESTED, UNTESTED]


def test_ambiguous_criterion_failure_is_read_from_the_trail():
    c = mk(cls="ambiguous", dicer_call=27, start=0, end=500,
           evidence="nearest=piRNA; below_cluster_criterion(length); size=1.0",
           e={"1u": 0.5, "pingpong": 0.0, "5a": 0, "fractop_bid": 0, "phase": 0})
    ck = rc.checklist(c, "animal")
    assert statuses(ck) == [PASS, PASS, FAIL, PASS]
    both = mk(cls="ambiguous", dicer_call=27, start=0, end=500,
              evidence="nearest=piRNA; below_cluster_criterion(length,density)")
    assert statuses(rc.checklist(both, "animal"))[2:] == [FAIL, FAIL]


def test_ambiguous_without_size_call():
    c = mk(cls="ambiguous", dicer_call=None, evidence="nearest=siRNA; x=1")
    ck = rc.checklist(c, "animal")
    assert ck[0] == "siRNA" and ck[1][0] == FAIL and "no size call" in ck[1][1]


def test_ambiguous_without_trail_or_ml_has_none():
    assert rc.checklist(mk(cls="ambiguous", evidence=""), "animal") is None
    assert rc.checklist(mk(cls="ambiguous", evidence="ml:rf; other"), "animal") is None


def test_plant_phasirna_rows():
    plan = rc.skeleton("plant")["phasiRNA"]
    assert [r[0] for r in plan][:2] == ["check.size.title", "check.phasing.title"]
    c = mk(cls="phasiRNA", dicer_call=21, phase_pvalue=0.3, phase_number=7,
           e={"1u": 0, "5a": 0, "fractop_bid": 0, "pingpong": 0, "phase": 0.3}, start=0, end=300)
    ck = rc.checklist(c, "plant")
    assert statuses(ck)[:2] == [PASS, PASS]       # phasiRNA has no p limit
    assert "7 cycles" in ck[2][1] and "0.3" in ck[2][1]


def test_plant_phasirna24_p_limit_gates_phasing():
    plan = rc.skeleton("plant")["phasiRNA24"]
    assert "0.01" in plan[1][1]
    ok = mk(cls="phasiRNA24", dicer_call=24, phase_pvalue=0.001, phase_number=6,
            e={"1u": 0, "5a": 0, "fractop_bid": 0, "pingpong": 0, "phase": 0.2}, start=0, end=500)
    assert statuses(rc.checklist(ok, "plant"))[:2] == [PASS, PASS]
    bad = mk(cls="ambiguous", dicer_call=24, phase_pvalue=0.2, phase_number=6,
             evidence="nearest=phasiRNA24; x=1",
             e={"1u": 0, "5a": 0, "fractop_bid": 0, "pingpong": 0, "phase": 0.2}, start=0, end=500)
    assert rc.checklist(bad, "plant")[2][0] == FAIL
    # p missing: the decider treats it as no phasing signal
    nop = mk(cls="ambiguous", dicer_call=24, phase_pvalue=None, evidence="nearest=phasiRNA24",
             e={"1u": 0, "5a": 0, "fractop_bid": 0, "pingpong": 0, "phase": 0.2}, start=0, end=500)
    assert rc.checklist(nop, "plant")[2][0] == FAIL


def test_plant_hcsirna_5a_and_size():
    c = mk(cls="hcsiRNA", dicer_call=24, pct_5a=0.7,
           e={"1u": 0, "5a": 0.5, "fractop_bid": 0, "pingpong": 0, "phase": 0})
    c.pct_5a = 0.7
    ck = rc.checklist(c, "plant")
    assert [r[0] for r in rc.skeleton("plant")["hcsiRNA"]] == ["check.size.title", "check.5a.title"]
    assert statuses(ck) == [PASS, PASS] and "70%" in ck[2][1]


def test_every_row_has_name_and_threshold_strings():
    for kingdom in ("animal", "plant"):
        for cls, plan in rc.skeleton(kingdom).items():
            assert plan, cls
            for name_key, thr in plan:
                assert name_key in rt.STRINGS, name_key
                assert thr and "{" not in thr


def test_checklist_follows_class_specs_slots():
    for kingdom in ("animal", "plant"):
        for cls, spec in report_data.class_specs(kingdom).items():
            extra = (spec["min_length"] is not None) + (spec["density_p"] is not None)
            assert len(rc.skeleton(kingdom)[cls]) == len(spec["slots"]) + extra


def test_payload_budget_per_page():
    c = mk(cls="piRNA", dicer_call=27, start=0, end=5000, pingpong_z=12.3,
           e={"1u": 0.6, "pingpong": 0.4, "5a": 0, "fractop_bid": 0, "phase": 0})
    assert len(json.dumps(rc.checklist(c, "animal"), separators=(",", ":"))) <= 300
    p = mk(cls="phasiRNA24", dicer_call=24, phase_pvalue=0.0012, phase_number=8, start=0, end=500,
           e={"1u": 0, "5a": 0, "fractop_bid": 0, "pingpong": 0, "phase": 0.3})
    assert len(json.dumps(rc.checklist(p, "plant"), separators=(",", ":"))) <= 300


def test_payload_carries_checklist_and_skeleton(tmp_path):
    from test_html_report_html import build, data_of
    cl = [mk("a", cls="piRNA", dicer_call=27, start=0, end=5000), mk("b", cls="miRNA"),
          mk("c", cls="ambiguous", dicer_call=19, evidence="nearest=piRNA; x=1")]
    sel, html, _ = build(tmp_path, cl)
    d = data_of(html)
    assert set(d["meta"]["cl"]) == {"piRNA", "siRNA"}
    by_name = {r["n"]: str(i) for i, r in enumerate(d["rows"])}
    assert d["details"][by_name["a"]]["ck"][0] == "piRNA"
    assert "ck" not in d["details"][by_name["b"]]
    assert d["details"][by_name["c"]]["ck"][0] == "piRNA"
    assert len(d["details"][by_name["a"]]["ck"]) - 1 == len(d["meta"]["cl"]["piRNA"])


def test_plant_payload_skeleton(tmp_path):
    from test_html_report_html import build, data_of
    cl = [mk("a", cls="hcsiRNA", dicer_call=24, start=0, end=300)]
    _, html, _ = build(tmp_path, cl)
    assert "hcsiRNA" in data_of(html)["meta"]["cl"]
