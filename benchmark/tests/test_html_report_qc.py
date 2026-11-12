"""Library quality strip: flag logic, fastp / qc_summary parsing, payload."""
import json
import os
import sys

import pytest
import yaml

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import report_qc as q
import report_text as rt

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")


def counts(**by_len):
    """27 slots from {"l22": 10, ...}; "lt" and "gt" are the overflow bins."""
    c = [0.0] * 27
    for k, v in by_len.items():
        c[0 if k == "lt" else 26 if k == "gt" else int(k[1:]) - 15] = v
    return c


def flag(res, key):
    return next(f for f in res["fl"] if f[0] == key)


def test_share_20_25_threshold_and_status():
    ok = q.build({"counts": counts(l22=25, l30=75), "unmapped": 0})
    assert flag(ok, "share")[1:] == [0.25, 1]               # exactly 25% passes
    low = q.build({"counts": counts(l22=24, l30=76), "unmapped": 0})
    assert flag(low, "share")[1:] == [0.24, 0]


def test_share_counts_20_to_25_inclusive():
    r = q.build({"counts": counts(l19=10, l20=10, l25=10, l26=10), "unmapped": 0})
    assert flag(r, "share")[1] == 0.5


def test_dominant_length_and_range_status():
    r = q.build({"counts": counts(l22=50, l24=30, l30=20), "unmapped": 0})
    assert flag(r, "dominant")[1:] == [22, 1]
    r = q.build({"counts": counts(l17=60, l22=40), "unmapped": 0})
    assert flag(r, "dominant")[1:] == [17, 0]               # below the 18 nt the classifier keeps
    r = q.build({"counts": counts(l36=60, l22=40), "unmapped": 0})
    assert flag(r, "dominant")[1:] == [36, 0]
    assert flag(q.build({"counts": counts(l35=1), "unmapped": 0}), "dominant")[2] == 1
    assert flag(q.build({"counts": counts(l18=1), "unmapped": 0}), "dominant")[2] == 1


def test_dominant_in_overflow_bin_is_named_and_flagged():
    r = q.build({"counts": counts(lt=70, l22=30), "unmapped": 0})
    assert flag(r, "dominant")[1:] == ["<16", 0]
    r = q.build({"counts": counts(gt=70, l22=30), "unmapped": 0})
    assert flag(r, "dominant")[1:] == [">40", 0]


def test_dominant_ties_take_the_shorter_length():
    assert flag(q.build({"counts": counts(l22=5, l24=5), "unmapped": 0}), "dominant")[1] == 22


def test_short_share_counts_below_18_including_overflow():
    r = q.build({"counts": counts(lt=10, l16=10, l17=5, l18=75), "unmapped": 0})
    assert flag(r, "short")[1] == 0.25 and flag(r, "short")[2] == 1      # 25% is not above 25%
    r = q.build({"counts": counts(lt=10, l16=10, l17=6, l18=74), "unmapped": 0})
    assert flag(r, "short")[2] == 0


def test_empty_library_gives_nothing():
    assert q.build({"counts": [0.0] * 27, "unmapped": 5}) is None
    assert q.build(None) is None


def test_payload_is_compact_and_json_safe():
    r = q.build({"counts": counts(l22=1234.5678, l24=10), "unmapped": 7})
    assert r["u"] == 7 and len(r["c"]) == 27 and r["c"][7] == 1234.6
    assert len(json.dumps(r, allow_nan=False)) < 600
    assert "tr" not in r


def test_extractor_range_matches_the_shipped_rules():
    for f in ("rules.yaml", "rules.plant.yaml"):
        ext = yaml.safe_load(open(os.path.join(ROOT, f)))["extractor"]
        assert (ext["dicermin"], ext["dicermax"]) == q.KEPT_RANGE


def test_strings_exist_for_every_flag_and_status():
    for key in ("share", "dominant", "short"):
        for st in (0, 1):
            assert f"lq.{key}.{st}" in rt.STRINGS
        assert f"lq.{key}.label" in rt.STRINGS


# ---- trimming numbers ----

FASTP = {"summary": {"before_filtering": {"total_reads": 1000}, "after_filtering": {"total_reads": 900}},
         "adapter_cutting": {"adapter_trimmed_reads": 850}}


def test_fastp_json(tmp_path):
    p = tmp_path / "x.fastp.json"
    p.write_text(json.dumps(FASTP))
    assert q.load_trim(fastp=str(p)) == {"b": 1000, "a": 900, "ad": 850}


def test_qc_summary_only_has_before_and_after(tmp_path):
    p = tmp_path / "x.qc_summary.txt"
    p.write_text("QC Report\nTrimming Statistics:\n------\n  Reads before: 1,000\n"
                 "  Reads after:  900 (90.0%)\n")
    assert q.load_trim(qc_summary=str(p)) == {"b": 1000, "a": 900}


def test_fastp_wins_and_qc_summary_fills_gaps(tmp_path):
    f = tmp_path / "x.fastp.json"
    f.write_text(json.dumps({"adapter_cutting": {"adapter_trimmed_reads": 5}}))
    t = tmp_path / "x.qc_summary.txt"
    t.write_text("  Reads before: 20\n  Reads after:  10 (50.0%)\n")
    assert q.load_trim(fastp=str(f), qc_summary=str(t)) == {"b": 20, "a": 10, "ad": 5}


def test_absent_or_broken_files_omit_the_part(tmp_path):
    assert q.load_trim() is None
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    assert q.load_trim(fastp=str(bad)) is None
    junk = tmp_path / "j.txt"
    junk.write_text("nothing here\n")
    assert q.load_trim(qc_summary=str(junk)) is None


def test_trim_goes_into_the_payload():
    r = q.build({"counts": counts(l22=10), "unmapped": 0}, {"b": 1000, "a": 900, "ad": 850})
    assert r["tr"] == {"b": 1000, "a": 900, "ad": 850}


def test_trim_paths_beside_clusters(tmp_path):
    (tmp_path / "run.fastp.json").write_text("{}")
    c = str(tmp_path / "run.clusters.tsv")
    assert q.trim_paths(c, None, None) == (str(tmp_path / "run.fastp.json"), None)
    assert q.trim_paths(c, "/e/f.json", "/e/q.txt") == ("/e/f.json", "/e/q.txt")
    (tmp_path / "run.qc_summary.txt").write_text("")
    assert q.trim_paths(c, None, None)[1] == str(tmp_path / "run.qc_summary.txt")
    assert q.trim_paths(str(tmp_path / "other.tsv"), None, None) == (None, None)


def test_flag_strings_quote_the_constants():
    S = rt.STRINGS
    assert f"{SHARE_PCT}%" in S["lq.share.1"] and f"{SHARE_PCT}%" in S["lq.share.0"]
    assert f"{q.SHARE_LO} to {q.SHARE_HI} nt" in S["lq.share.label"]
    assert f"{q.KEPT_RANGE[0]} to {q.KEPT_RANGE[1]} nt" in S["lq.dominant.1"]
    assert f"{q.SHORT_BELOW} nt" in S["lq.short.label"]
    assert f"{round(q.SHORT_MAX * 100)}%" in S["lq.short.1"] + S["lq.short.0"]


SHARE_PCT = round(q.SHARE_MIN * 100)
