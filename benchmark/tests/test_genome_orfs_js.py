"""Task 7: the EXPERIMENTAL predicted-ORF track in the genomic-context figure (stub DOM, JavaScriptCore)."""
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import test_genome_figure_js as t5  # noqa: E402
from test_genome_figure_js import data_with, gm_payload, run  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("osascript") is None,
                                reason="needs osascript (JavaScriptCore)")

O = [[1500, 1800, "+", 99], [3000, 3600, "-", 199], [4000, 4300, ".", 99]]
TE_COLORS = {"#c46a2a", "#6b4d99", "#a8413a", "#2c5f8d", "#137d6f", "#8493a3"}


def orf_data(tmp_path, o=O, min_aa=100, **over):
    gm = gm_payload(**over)
    if o is not None:
        gm["o"] = o
    d = data_with(tmp_path, gm)
    if min_aa:
        d["meta"]["orf_min_aa"] = min_aa
    return d


GUTTER = """
function gutter() {
  return H.all(gmEl, "text").filter(function (t) { return /rpt-gm-(tl|ts)/.test(H.cls(t)); })
    .map(H.text).join(" ").replace(/\\s+/g, " ");
}
"""


def test_no_track_without_o_or_flag(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = GUTTER + """
    return {g: gutter(), bars: H.all(gmEl, "g").filter(function (g) { return /rpt-gm-orf/.test(H.cls(g)); }).length,
            html: H.serialize(gmEl), empties: H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-empty"; }).length};
    """
    r = run(tmp_path, d, body)
    S = d["strings"]
    assert "Predicted ORFs" not in r["g"] and S["ui.gm.orf.note"] not in r["g"]
    assert r["bars"] == 0 and "orf" not in r["html"].lower()
    assert r["empties"] == 2      # te layer given with rows: only ncrna and phas (no layer) ... see below


def test_track_with_o_after_phas_with_note_and_legend(tmp_path):
    d = orf_data(tmp_path)
    S = d["strings"]
    body = GUTTER + """
    var tl = H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-tl"; }).map(H.text);
    return {g: gutter(), tl: tl, leg: H.text(H.byClass(gmEl, "rpt-srna-legend", "ul")[0])};
    """
    r = run(tmp_path, d, body)
    assert S["ui.gm.track.orf"].split(" (")[0] in r["g"]
    assert "(experimental)" in r["g"]
    assert " ".join(S["ui.gm.orf.note"].split()) in r["g"]
    assert r["tl"].index("PHAS loci") < len(r["tl"]) - 1      # something comes after PHAS
    assert S["ui.gm.orf.legend"] in r["leg"]


def test_note_text_is_the_agreed_wording(tmp_path):
    S = orf_data(tmp_path)["strings"]
    assert S["ui.gm.orf.note"] == "experimental: predicted from sequence, not annotation"
    assert S["ui.gm.orf.empty"] == "No ORFs of at least {n} aa outside annotated CDS"
    assert S["ui.gm.orf.tip"] == "Predicted, not annotated"


def test_bars_strand_chevrons_and_own_colour(tmp_path):
    body = """
    var bars = H.all(gmEl, "g").filter(function (g) { return /rpt-gm-orf/.test(H.cls(g)); });
    return bars.map(function (b) {
      var rect = H.all(b, "rect")[0], ch = H.all(b, "polyline");
      var pts = ch.length ? ch[0].getAttribute("points").split(" ").map(function (p) { return parseFloat(p.split(",")[0]); }) : [];
      return {n: ch.length, dir: pts.length ? (pts[1] > pts[0] ? 1 : -1) : 0, fill: rect.getAttribute("fill"),
              stroke: rect.getAttribute("stroke"), track: b.getAttribute("data-track")};
    });
    """
    r = run(tmp_path, orf_data(tmp_path), body)
    assert len(r) == 3
    assert r[0]["n"] > 0 and r[0]["dir"] == 1          # +
    assert r[1]["n"] > 0 and r[1]["dir"] == -1         # -
    assert r[2]["n"] == 0                              # unstranded: no chevron
    for b in r:
        assert b["track"] == "o"
        assert b["fill"] not in TE_COLORS and b["stroke"] not in TE_COLORS


def test_empty_state_uses_the_min_length(tmp_path):
    d = orf_data(tmp_path, o=None, min_aa=150)
    S = d["strings"]
    body = """
    return H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-empty"; }).map(H.text);
    """
    notes = run(tmp_path, d, body)
    assert S["ui.gm.orf.empty"].replace("{n}", "150") in notes
    body2 = GUTTER + "return gutter();"
    assert " ".join(S["ui.gm.orf.note"].split()) in run(tmp_path, d, body2)


def test_table_rows_have_the_track_name(tmp_path):
    d = orf_data(tmp_path)
    S = d["strings"]
    body = """
    H.click(H.buttons(gmEl, %s)[0]);
    var tb = H.byClass(gmEl, "rpt-srna-fig__table", "div")[0];
    return H.all(tb, "tr").slice(1).map(function (tr) { return H.all(tr, "td").map(H.text); })
      .filter(function (r) { return r[0] === %s; });
    """ % (repr(S["ui.fig.table"]).replace("'", '"'), repr(S["ui.gm.th.orf"]).replace("'", '"'))
    rows = run(tmp_path, d, body)
    assert S["ui.gm.th.orf"] == "Predicted ORFs"
    assert len(rows) == 3
    assert rows[0][3] == "+" and rows[0][4:6] == ["1,501", "1,800"] and rows[1][3] == "−"
    assert S["ui.gm.orf.len"].replace("{n}", "99") in rows[0]


def test_tooltip_has_length_strand_coordinates_and_the_caveat(tmp_path):
    d = orf_data(tmp_path)
    S = d["strings"]
    body = """
    return H.all(gmEl, "rect").filter(function (r) { return r.getAttribute("data-track") === "o"; })
      .map(function (r) { return r.getAttribute("aria-label"); });
    """
    tips = run(tmp_path, d, body)
    assert len(tips) == 3
    t = tips[0]
    assert S["ui.gm.orf.len"].replace("{n}", "99") in t and "1,501" in t and "1,800" in t
    assert S["ui.gm.orf.tip"] in t and "frame" not in t.lower()


def test_more_note_from_mo_o(tmp_path):
    d = orf_data(tmp_path, mo={"o": 5})
    body = """
    return H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-more"; }).map(H.text);
    """
    assert run(tmp_path, d, body) == [d["strings"]["ui.gm.more"].replace("{n}", "5")]


def test_orfs_follow_the_view(tmp_path):
    d = orf_data(tmp_path)
    body = """
    function n() { return H.all(gmEl, "g").filter(function (g) { return /rpt-gm-orf/.test(H.cls(g)); }).length; }
    var before = n();
    H.click(H.buttons(gmEl, %s)[0]); H.frames();     // the locus view (3000..3100): only the minus ORF at 3000..3600 shows
    return {before: before, after: n()};
    """ % ('"' + d["strings"]["ui.gm.locus"] + '"')
    r = run(tmp_path, d, body)
    assert r["before"] == 3 and r["after"] == 1


def test_o_text_is_never_html(tmp_path):
    d = orf_data(tmp_path)
    body = "return H.doc('b').length;"
    assert run(tmp_path, d, body) == 0
