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
    # more than 50 clusters, so the worker pool (not the single-process path) runs
    a = args_for(tmp_path, threads=2)
    rows = ["\t".join(map(str, ["cl%d" % i, "c1", i * 10, i * 10 + 200, 500, 50, "piRNA", 0.9,
                                26, 0.8, 0.5, 0.8, 0.5, 3, 0.5, 0.1, 0.9, "False", "UGA", ""]))
            for i in range(120)]
    open(a.clusters, "w").write("\n".join(["\t".join(COLS)] + rows) + "\n")
    report_cli.run_report(a)
    out = tmp_path / "out" / "x" / "report.report.html"
    m = re.search(r'<script type="application/json" id="report-data">(.*?)</script>',
                  out.read_text(), re.S)
    assert len(json.loads(m.group(1))["rows"]) == 120


def test_warns_on_very_large_max_detail(tmp_path, capsys):
    report_cli.run_report(args_for(tmp_path, max_detail=5001))
    err = capsys.readouterr().err
    assert "WARN" in err and "--max-detail" in err and "7 MB" in err
    report_cli.run_report(args_for(tmp_path, max_detail=5000))
    assert "7 MB" not in capsys.readouterr().err


def test_max_detail_default_is_3000():
    import srnaprofiler
    d = srnaprofiler.build_parser().parse_args(
        ["report", "--clusters", "a", "--bam", "b", "-o", "o"])
    assert d.max_detail == 3000


def test_readme_mentions_threads_and_no_glossary():
    text = open(os.path.join(os.path.dirname(__file__), "..", "..", "README.md")).read()
    assert "--threads" in text
    assert "glossary" not in text.lower()


def test_banner_uses_shared_version(capsys):
    import srnaprofiler
    import version
    srnaprofiler.print_header()
    out = capsys.readouterr().out
    assert f"v{version.__version__}" in out


def test_version_single_source():
    import subprocess
    import version
    root = os.path.join(os.path.dirname(__file__), "..", "..")
    r = subprocess.run([sys.executable, os.path.join(root, "srnaprofiler.py"), "--version"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and version.__version__ in (r.stdout + r.stderr)
    assert "Traceback" not in r.stderr


# ---- library quality strip flags ----

def _data(tmp_path):
    out = tmp_path / "out" / "x" / "report.report.html"
    m = re.search(r'<script type="application/json" id="report-data">(.*?)</script>',
                  out.read_text(), re.S)
    return json.loads(m.group(1))


def test_library_qc_parser_flags():
    import srnaprofiler
    p = srnaprofiler.build_parser()
    d = p.parse_args(["report", "--clusters", "a", "--bam", "b", "-o", "o"])
    assert (d.no_library_qc, d.fastp, d.qc_summary) == (False, None, None)
    n = p.parse_args(["report", "--clusters", "a", "--bam", "b", "-o", "o", "--no-library-qc",
                      "--fastp", "f.json", "--qc-summary", "q.txt"])
    assert (n.no_library_qc, n.fastp, n.qc_summary) == (True, "f.json", "q.txt")


def test_report_has_library_qc_by_default(tmp_path):
    report_cli.run_report(args_for(tmp_path))
    qc = _data(tmp_path)["library"]["qc"]
    assert qc["c"][5] == 3 and sum(qc["c"]) == 3            # the three 20-nt test reads
    assert qc["fl"][1][:2] == ["dominant", 20] and "tr" not in qc


def test_no_library_qc_skips_the_pass(tmp_path, monkeypatch):
    import report_signals

    def boom(*a, **k):
        raise AssertionError("BAM pass must not run")
    monkeypatch.setattr(report_signals, "library_lengths", boom)
    report_cli.run_report(args_for(tmp_path, no_library_qc=True))
    assert "qc" not in _data(tmp_path)["library"]


def test_trim_files_beside_clusters_and_explicit(tmp_path):
    (tmp_path / "run.fastp.json").write_text(json.dumps(
        {"summary": {"before_filtering": {"total_reads": 100}, "after_filtering": {"total_reads": 90}},
         "adapter_cutting": {"adapter_trimmed_reads": 80}}))
    report_cli.run_report(args_for(tmp_path))
    assert _data(tmp_path)["library"]["qc"]["tr"] == {"b": 100, "a": 90, "ad": 80}
    other = tmp_path / "other.txt"
    other.write_text("  Reads before: 50\n  Reads after:  40 (80.0%)\n")
    (tmp_path / "run.fastp.json").unlink()
    report_cli.run_report(args_for(tmp_path, qc_summary=str(other)))
    assert _data(tmp_path)["library"]["qc"]["tr"] == {"b": 50, "a": 40}


def test_missing_explicit_trim_file_is_an_error(tmp_path):
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(args_for(tmp_path, fastp=str(tmp_path / "nope.json")))
    assert "nope.json" in str(e.value)


def test_readme_documents_library_qc():
    text = open(os.path.join(os.path.dirname(__file__), "..", "..", "README.md")).read()
    for flag in ("--no-library-qc", "--fastp", "--qc-summary"):
        assert flag in text
