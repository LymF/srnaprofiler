import hashlib
import json
import os
import re
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))
import bamutil
import report_cli

COLS = ["name", "chrom", "start", "end", "reads", "reads_rpm", "cls", "confidence",
        "dicer_call", "size_purity", "frac_top", "pct_1u", "pct_10a", "pingpong_z",
        "phase_pvalue", "multimap_frac", "unique_frac", "review", "major_rna", "evidence"]


def write_clusters(path):
    rows = [["cl1", "c1", 100, 200, 500, 50, "piRNA", 0.9, 26, 0.8, 0.5, 0.8, 0.5, 3, 0.5, 0.1, 0.9, "False", "UGA", ""],
            ["cl2", "c2", 300, 400, 50, 5, "siRNA", 0.6, 21, 0.7, 0.5, 0.3, 0.2, 0, 0.5, 0.1, 0.9, "False", "UGA", ""],
            ["cl3", "chr1", 10, 90, 5, 0.5, "other", 0.2, 22, 0.5, 0.5, 0.3, 0.2, 0, 0.5, 0.1, 0.9, "False", "UGA", ""]]
    path.write_text("\n".join("\t".join(map(str, r)) for r in [COLS] + rows) + "\n")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def args_for(tmp_path, **kw):
    clusters = tmp_path / "run.clusters.tsv"
    write_clusters(clusters)
    bam = tmp_path / "s.bam"
    bamutil.make_bam(bam, [("c1", 100, False), ("c1", 110, False), ("c2", 300, False)])
    a = dict(clusters=str(clusters), bam=str(bam), outdir=str(tmp_path / "out" / "x"),
             prefix="report", min_reads=0.0, max_detail=3000, fragments=None, title=None,
             threads=1)
    a.update(kw)
    return SimpleNamespace(**a)


def test_report_end_to_end(tmp_path):
    a = args_for(tmp_path)
    report_cli.run_report(a)
    out = tmp_path / "out" / "x" / "report.report.html"
    assert out.exists()
    m = re.search(r'<script type="application/json" id="report-data">(.*?)</script>',
                  out.read_text(), re.S)
    assert len(json.loads(m.group(1))["rows"]) == 3


def test_min_reads_hides(tmp_path):
    a = args_for(tmp_path, min_reads=300.0)
    report_cli.run_report(a)
    m = re.search(r'<script type="application/json" id="report-data">(.*?)</script>',
                  (tmp_path / "out" / "x" / "report.report.html").read_text(), re.S)
    assert len(json.loads(m.group(1))["rows"]) == 1


def test_parser_accepts_flags():
    import srnaprofiler
    ns = srnaprofiler.build_parser().parse_args(
        ["report", "--clusters", "a.clusters.tsv", "--bam", "b.bam", "-o", "o", "-p", "P",
         "--min-reads", "10", "--max-detail", "5", "--fragments", "f.tsv", "--title", "T"])
    assert ns.func is report_cli.run_report
    assert (ns.min_reads, ns.max_detail, ns.prefix, ns.title) == (10.0, 5, "P", "T")
    d = srnaprofiler.build_parser().parse_args(
        ["report", "--clusters", "a", "--bam", "b", "-o", "o"])
    assert (d.prefix, d.min_reads, d.max_detail, d.fragments) == ("report", 300.0, 3000, None)


def test_inputs_untouched(tmp_path):
    a = args_for(tmp_path)
    before = (sha(a.clusters), sha(a.bam))
    report_cli.run_report(a)
    assert (sha(a.clusters), sha(a.bam)) == before


def test_fragments_default_beside_clusters(tmp_path):
    a = args_for(tmp_path)
    (tmp_path / "run.fragments.tsv").write_text("reads\n5\n7\n")
    assert report_cli._fragments_path(a) == str(tmp_path / "run.fragments.tsv")
    assert report_cli._fragments_path(args_for(tmp_path, fragments="x.tsv")) == "x.tsv"


def test_missing_inputs_clear_error(tmp_path, capsys):
    a = args_for(tmp_path, bam=str(tmp_path / "nope.bam"))
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert e.value.code != 0
    assert "nope.bam" in capsys.readouterr().err + str(e.value)
    a = args_for(tmp_path, clusters=str(tmp_path / "nope.clusters.tsv"))
    with pytest.raises(SystemExit):
        report_cli.run_report(a)


def test_missing_bam_index_is_a_clear_error(tmp_path):
    a = args_for(tmp_path)
    os.remove(a.bam + ".bai")
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert "index" in str(e.value)


@pytest.mark.parametrize("kw", [{"max_detail": -1}, {"min_reads": float("nan")},
                                {"min_reads": float("inf")}, {"min_reads": -1.0}])
def test_bad_numbers_are_clear_errors(tmp_path, kw):
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(args_for(tmp_path, **kw))
    assert "report:" in str(e.value)


@pytest.mark.parametrize("bad", ["0", "-1", "abc", "1.5"])
def test_threads_must_be_positive_int(bad, capsys):
    import srnaprofiler
    with pytest.raises(SystemExit) as e:
        srnaprofiler.build_parser().parse_args(
            ["report", "--clusters", "a", "--bam", "b", "-o", "o", "--threads", bad])
    assert e.value.code != 0
    assert "--threads" in capsys.readouterr().err


def test_threads_default_and_run_report_guard(tmp_path):
    import srnaprofiler
    d = srnaprofiler.build_parser().parse_args(
        ["report", "--clusters", "a", "--bam", "b", "-o", "o"])
    assert d.threads == 1
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(args_for(tmp_path, threads=0))
    assert "--threads" in str(e.value)


def test_report_with_threads(tmp_path):
    report_cli.run_report(args_for(tmp_path, threads=2))
    assert (tmp_path / "out" / "x" / "report.report.html").exists()


def test_version_single_source():
    import subprocess
    import version
    root = os.path.join(os.path.dirname(__file__), "..", "..")
    r = subprocess.run([sys.executable, os.path.join(root, "srnaprofiler.py"), "--version"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and version.__version__ in (r.stdout + r.stderr)
    assert "Traceback" not in r.stderr
