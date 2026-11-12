"""The miRNA hairpin figure under JavaScriptCore with the stub DOM."""
import json
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import report_html
import report_mirna as rm
from test_html_report_html import SUMMARY, data_of, full_sig, mk
from test_html_report_js import run_js
from test_html_report_mirna import (END, HI, LO, M_AB, S_AB, START, load_one, precursor)

pytestmark = pytest.mark.skipif(shutil.which("osascript") is None,
                                reason="needs osascript (JavaScriptCore)")

PROBE = """
var res = {};
H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
var det = H.id("sec-detail");
res.figs = H.byClass(det, "rpt-srna-fig", "section").map(function (f) { return f.getAttribute("data-fig"); });
var fig = H.byClass(det, "rpt-srna-fig", "section").filter(function (f) { return f.getAttribute("data-fig") === "hairpin"; })[0];
res.has = !!fig;
if (fig) {
  var svg = H.all(fig, "svg")[0];
  res.width = svg.getAttribute("width");
  res.height = svg.getAttribute("height");
  var cls = function (c) { return H.byClass(svg, c, "rect").length; };
  res.bands = cls("mk-band"); res.tints = cls("mk-tint"); res.reads = cls("mk-read");
  res.mis = cls("mk-mis"); res.bars = cls("mk-bar");
  res.texts = H.all(svg, "text").map(H.text);
  res.buttons = H.all(fig, "button").map(H.text);
  res.legend = H.all(H.byClass(fig, "rpt-srna-legend", "ul")[0], "li").map(H.text);
  res.notes = H.byClass(fig, "rpt-srna-note", "p").map(H.text);
  res.scroll = H.byClass(fig, "rpt-srna-fig__plot--scroll", "div").length;
  H.click(H.buttons(fig, TBL)[0]);
  var tbls = H.all(H.byClass(fig, "rpt-srna-fig__table", "div")[0], "table");
  res.tables = tbls.map(function (t) { return H.all(t, "tr").length; });
  res.tblText = H.text(tbls[0]).slice(0, 80);
  H.click(H.buttons(fig, TBL_BACK)[0]);
  H.click(H.buttons(fig, DL)[0]);
  res.download = H.stats().downloads.length ? H.stats().downloads[0].name : null;
  res.blob = H.stats().lastBlob || "";
}
res.detailText = H.text(det);
var secs = H.byClass(det, "rpt-srna-ck", "section");
res.ck = secs.map(function (s) { return H.text(H.all(s, "h3")[0]); });
res.rowsOf = secs.map(function (s) { return H.all(s, "tr").length; });
res.script = H.doc("script").length;
return res;
"""


def build(tmp_path, mi_args, strand="+", cls="miRNA", fold=("((((..))))", -3.2, False), hp=None):
    row = load_one(tmp_path, strand=strand, **(mi_args or {}))
    prec = row.precursor_seq or precursor()
    rd = [[M_AB[0], M_AB[1], prec[M_AB[0]:M_AB[1]].lower(), 12.0, "m"],
          [M_AB[0], M_AB[1], prec[M_AB[0]:M_AB[1]].lower()[:5] + prec[M_AB[0] + 5].upper() +
           prec[M_AB[0]:M_AB[1]].lower()[6:], 3.0, "m"],
          [S_AB[0], S_AB[1], prec[S_AB[0]:S_AB[1]].lower(), 4.5, "s"],
          [52, 72, prec[52:72].lower(), 1.0, ""]]
    cov = [0.0] * (HI - LO)
    for a, b, _, w, _k in rd:
        for i in range(a, b):
            cov[i] += w
    hp = hp if hp is not None else {"rd": rd, "cov": cov, "more": 3, "mr": 2.0, "anti": 1.0, "ovh": 0.0}
    c = mk("mir1", reads=900, cls=cls, start=START, end=END)
    from report_data import select
    sel = select([c], 300.0)
    sig = {c.name: dict(full_sig(), hp=hp)}
    out = tmp_path / "r.html"
    report_html.build_report(sel, SUMMARY, sig, "t", str(out), mirna={"mir1": (row, fold)})
    d = data_of(out.read_text())
    return d, json.dumps(d)


def run(tmp_path, text, d, probe=PROBE, width=720):
    S = d["strings"]
    pre = "var TBL = %s; var TBL_BACK = %s; var DL = %s;\n" % (
        json.dumps(S["ui.fig.table"]), json.dumps(S["ui.fig.chart"]), json.dumps(S["ui.fig.svg"]))
    return run_js(tmp_path, text, pre + probe, width=width)


def test_hairpin_figure_renders(tmp_path):
    d, text = build(tmp_path, None)
    assert "mi" in d["details"]["0"]
    r = run(tmp_path, text, d)
    S = d["strings"]
    assert r["has"] and r["figs"][0] == "hairpin" and r["figs"][1:] == list("ABCDEFGH")
    assert r["bands"] == 3 and r["tints"] >= 3                    # mature, loop, star
    assert r["reads"] == 4 and r["bars"] == 4 and r["mis"] == 1   # one mismatch cell
    mi = d["details"]["0"]["mi"]
    assert int(r["width"]) == 48 + mi["L"] * 8 + 196 and r["scroll"] == 1
    assert mi["seq"] in "".join(r["texts"])                       # the runs, in order, spell the precursor
    assert mi["db"] in r["texts"]
    for k in ("mi.seg.mature", "mi.seg.loop", "mi.seg.star"):     # segments named in text, not colour alone
        assert S[k] in r["texts"]
    assert S["mi.more"].format(n=3, reads=2) in r["texts"]
    assert r["buttons"] == [S["ui.fig.table"], S["ui.fig.svg"]]
    assert S["mi.seg.mature"] in r["legend"] and S["mi.lg.mismatch"] in r["legend"]
    assert r["tables"] == [mi["L"] + 1, 4 + 1]                    # one row per position; reads table
    assert r["download"].endswith("_hairpin.svg") and r["blob"].startswith("<svg")
    assert "hit" not in r["blob"]
    assert S["mi.reads.note"] in r["notes"]
    assert S["mi.struct"].format(mfe=-3.2) in r["notes"]
    assert S["mi.anti"].format(anti=1, ovh=0) in r["notes"]
    assert S["mi.num.title"] in r["ck"] and S["mi.crit.title"] in r["ck"]
    assert S["ck.mirna"] in r["detailText"]                       # the existing line is unchanged
    for need in ("96.4%", "ath-miR166a", "82%", "Pass", "120", "15"):
        assert need in r["detailText"], need
    assert r["script"] == 0


def test_criteria_rows_have_icon_and_text(tmp_path):
    d, text = build(tmp_path, {"mismatches": 7})
    r = run(tmp_path, text, d)
    S = d["strings"]
    assert S["ck.st.0"] in r["detailText"] and S["ck.st.1"] in r["detailText"]   # Not met and Met, as words
    assert S["mi.crit.note"] in r["detailText"]
    assert S["mi.rule.mm"].format(n=5) in r["detailText"]


def test_no_structure_says_why(tmp_path):
    d, text = build(tmp_path, None, fold=(None, None, True))
    r = run(tmp_path, text, d)
    S = d["strings"]
    assert S["mi.nostruct.m"] == "Structure not available: RNAfold not installed."
    assert S["mi.nostruct.m"] in r["notes"]
    assert not any(set(t) <= set("().") and len(t) > 20 for t in r["texts"])


def test_minus_strand_and_no_reads_and_no_precursor(tmp_path):
    d, text = build(tmp_path, None, strand="-")
    r = run(tmp_path, text, d)
    assert r["has"] and "on the − strand" in " ".join(r["notes"])
    # no reads fall in the precursor: no profile, no read rows, a note
    d, text = build(tmp_path, None, hp={"rd": [], "cov": [0.0] * (HI - LO), "more": 0, "mr": 0.0,
                                        "anti": 0.0, "ovh": 0.0})
    r = run(tmp_path, text, d)
    assert r["has"] and r["reads"] == 0 and d["strings"]["mi.noreads"] in r["notes"]
    assert r["tables"] == [HI - LO + 1]
    # no precursor sequence (miRScore did not run): a message, the numbers still show
    d, text = build(tmp_path, {"precursor_seq": "NA", "precursor_len": "NA"})
    r = run(tmp_path, text, d)
    assert not r["has"] and d["strings"]["mi.noseq"] in r["detailText"]
    assert d["strings"]["mi.num.title"] in r["ck"]


def test_unplaceable_window_keeps_the_sequence(tmp_path):
    d, text = build(tmp_path, {"end": END + 100})
    r = run(tmp_path, text, d)
    assert r["has"] and r["reads"] == 0 and d["strings"]["mi.nowin"] in r["notes"]


def test_pages_without_a_hairpin_are_unchanged(tmp_path):
    c = mk("p", reads=900)
    from report_data import select
    out = tmp_path / "p.html"
    report_html.build_report(select([c], 300.0), SUMMARY, {c.name: full_sig()}, "t", str(out))
    d2 = data_of(out.read_text())
    assert "mi" not in d2["details"]["0"]
    r = run(tmp_path, json.dumps(d2), d2)
    assert not r["has"] and r["figs"] == list("ABCDEFGH")
