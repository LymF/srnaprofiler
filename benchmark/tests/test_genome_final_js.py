"""Final-review fixes in the genomic-context figure (stub DOM, JavaScriptCore): M4, M5, M7, M8, M9, I-3."""
import json
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import test_genome_figure_js as t5  # noqa: E402
from test_genome_figure_js import data_with, gm_payload  # noqa: E402
from test_html_report_js import APP  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("osascript") is None, reason="needs osascript (JavaScriptCore)")

HOOK = "window.__gm = {gmBelowBin: gmBelowBin, gmZoom: gmZoom};\n  buildTabs();\n  renderAbout();"
APPX = APP.replace("  buildTabs();\n  renderAbout();", HOOK, 1)
APP_FILE = os.environ.get("GM_APP")
if APP_FILE:
    APPX = open(APP_FILE, encoding="utf-8").read()


def run(tmp_path, d, body="return res;", width=900):
    return t5.run_js(tmp_path, json.dumps(d), t5.OPEN + body, app=APPX, width=width)


def test_below_bin_helper(tmp_path):
    d = data_with(tmp_path, gm_payload())
    r = run(tmp_path, d, "return [window.__gm.gmBelowBin([0, 39], 20), window.__gm.gmBelowBin([0, 40], 20), "
                         "window.__gm.gmBelowBin([0, 10], 1), window.__gm.gmBelowBin([0, 99], 50), "
                         "window.__gm.gmBelowBin([0, 100], 50)];")
    assert r == [True, False, False, True, False]


def test_zoomed_below_bin_shows_the_note_and_zoom_out_removes_it(tmp_path):
    d = data_with(tmp_path, gm_payload(L=[3000, 3030]))     # 50 bp bins; the locus view is 40 bp
    S = d["strings"]
    body = """
    var rng = H.byClass(gmEl, "rpt-srna-fig__n", "span")[0];
    var lbl = %s;
    var loc = H.all(gmEl, "button").filter(function (b) { return b.getAttribute("aria-label") === lbl; })[0];
    res.before = H.text(rng);
    H.click(loc);
    res.after = H.text(rng);
    return res;
    """ % json.dumps(S["ui.gm.locus"])
    r = run(tmp_path, d, body)
    note = S["ui.gm.below"].replace("{n}", "50")
    assert note not in r["before"] and r["after"].endswith(note)


def test_end_unit_is_per_bin(tmp_path):
    assert t5.run_js is not None
    from report_text import STRINGS
    assert "base" not in STRINGS["ui.gm.ends.unit"] and "bin" in STRINGS["ui.gm.ends.unit"]


def test_layer_without_contig_has_its_own_message(tmp_path):
    d = data_with(tmp_path, gm_payload(g=[], t=[], ly=["gene"], nm=["te"]))
    S = d["strings"]
    body = """
    return H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-empty"; }).map(H.text);
    """
    notes = run(tmp_path, d, body)
    assert notes == [S["ui.gm.empty"], S["ui.gm.nocontig"], S["ui.gm.nolayer"].replace("{kind}", "ncrna"),
                     S["ui.gm.nolayer"].replace("{kind}", "phas")]


def test_unnamed_feature_tooltip_and_aria_fall_back_to_the_type(tmp_path):
    gm = gm_payload(g=[["", "+", 1500, 2700, 0, [["E", 1500, 1800]]]],
                    t=[["", "LTR/Gypsy", "+", 3000, 3200, 1]], r=[["", "tRNA", "+", 3300, 3400, 0]])
    d = data_with(tmp_path, gm)
    body = """
    return H.byClass(gmEl, "rpt-gm-hit", "rect").map(function (e) { return e.getAttribute("aria-label"); });
    """
    r = run(tmp_path, d, body)
    assert r and all(a and not a.startswith(";") and not a.startswith(" ") for a in r)
    assert any(a.startswith("Gene; ") for a in r) and any(a.startswith("LTR/Gypsy; ") for a in r)
    assert any(a.startswith("tRNA; ") for a in r)


def test_download_name_has_no_stray_letter(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    H.click(H.buttons(gmEl, %s)[0]);
    return H.stats().downloads[0].name;
    """ % json.dumps(d["strings"]["ui.fig.svg"])
    name = run(tmp_path, d, body)
    assert name.endswith("_context.svg") and not name.endswith("_I.svg")


def test_precursor_marker_wording(tmp_path):
    d = data_with(tmp_path, gm_payload(), mi={"n": "x", "lo": 2900, "L": 120, "nowin": 0, "t": [], "pc": []})
    S = d["strings"]
    body = """
    res.txt = H.all(gmEl, "text").map(H.text).join("|");
    res.leg = H.text(H.byClass(gmEl, "rpt-srna-legend", "ul")[0]);
    res.aria = H.byClass(gmEl, "rpt-gm-hit", "rect").map(function (e) { return e.getAttribute("aria-label"); }).join("|");
    return res;
    """
    r = run(tmp_path, d, body)
    for k in ("txt", "leg", "aria"):
        assert S["ui.gm.prec"] in r[k], k
        assert "miRNA hairpin" not in r[k], k
