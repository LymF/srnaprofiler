import json
from pathlib import Path

from benchmark.manifest import Library
from benchmark.run import pipeline_sha256
from benchmark.report import LibraryResult, build_report, collect_results


def lib(acc, role, organism="Bombyx mori", group="arthropod"):
    return Library(acc, organism, role, group, "ovary", 1, "http://x", ["t"])


def metrics(precision=0.9, fdr=0.1):
    return {
        "classification": {
            "per_class": {"piRNA": {"precision": precision, "fdr": fdr, "recall": 0.5, "tp": 1, "fp": 1},
                          "miRNA": {"precision": None, "fdr": None, "recall": None, "tp": 0, "fp": 0}},
            "ambiguous_coverage": 0.75,
            "unknown": {"count": 7},
            "ece": 0.0123,
        },
        "delimitation": {
            "global": {"tiny_counts": {"50": 3, "100": 5},
                       "low_read_counts": {"20": 8, "50": 9, "100": 10}},
            "by_source": {
                "SrcA": {"summary": {"n_clusters": 4, "mean_recovered_fraction": 0.6,
                                     "mean_fragments": 1.5}, "merged_distinct": 2},
                "SrcB": {"summary": {"n_clusters": 9, "mean_recovered_fraction": 0.7,
                                     "mean_fragments": 2.5}, "merged_distinct": 6},
            },
        },
        "by_source": {"SrcA": {"piRNA": {"precision": 0.8, "fdr": 0.2, "recall": 0.4, "tp": 1, "fp": 1}}},
        "n_loci": 10,
    }


def make_rules(tmp_path, text="a: 1\n"):
    rules = tmp_path / "rules.yaml"
    rules.write_text(text)
    (tmp_path / "pipeline.py").write_text("x = 1\n")
    return rules, pipeline_sha256(tmp_path)


def res(acc, role, sha, organism="Bombyx mori", **kw):
    meta = {"pipeline_sha256": sha, "rules_path": "/elsewhere/rules.yaml", "role": role}
    return LibraryResult(lib(acc, role, organism), metrics(**kw), meta)


def test_report_keeps_tuning_and_test_separate(tmp_path):
    _, sha = make_rules(tmp_path)
    text = build_report([
        res("SRRTUNE", "tuning", sha, "Drosophila melanogaster", precision=0.111),
        res("SRRTEST", "test", sha, precision=0.222),
    ], tmp_path)
    assert text.index("## Ajuste") < text.index("## Teste")
    ajuste, teste = text.split("## Teste")
    assert "SRRTUNE" in ajuste and "SRRTEST" not in ajuste
    assert "SRRTEST" in teste and "SRRTUNE" not in teste
    assert "0.111" in ajuste and "0.111" not in teste
    assert "0.222" in teste and "0.222" not in ajuste


def test_report_flags_rules_changed_after_test_run(tmp_path):
    rules, sha = make_rules(tmp_path)
    results = [res("SRRTEST", "test", sha, precision=0.777),
               res("SRRTUNE", "tuning", sha, "Drosophila melanogaster")]
    ok = build_report(results, tmp_path)
    assert "INVALID" not in ok and "stale" not in ok
    assert "mean over 1 valid" in ok

    rules.write_text("a: 2\n")
    text = build_report(results, tmp_path)
    ajuste, teste = text.split("## Teste")
    assert "INVALID" in teste
    assert "stale" in ajuste and "INVALID" not in ajuste
    assert "no valid test libraries" in teste
    aggregate = teste.split("### Aggregate")[1]
    assert "0.777" not in aggregate


def test_report_flags_pipeline_code_change_and_missing_hash(tmp_path):
    _, sha = make_rules(tmp_path)
    results = [res("SRRTEST", "test", sha), res("SRRTUNE", "tuning", sha, "Drosophila melanogaster")]
    assert "INVALID" not in build_report(results, tmp_path)
    (tmp_path / "pipeline.py").write_text("x = 2\n")  # code, not rules, changed
    text = build_report(results, tmp_path)
    ajuste, teste = text.split("## Teste")
    assert "INVALID" in teste and "stale" in ajuste
    # a result without pipeline_sha256 is never valid
    old = res("SRROLD", "test", sha)
    del old.meta["pipeline_sha256"]
    (tmp_path / "pipeline.py").write_text("x = 1\n")
    assert "INVALID" in build_report([old], tmp_path)


def test_report_shows_short_pipeline_hash_per_row(tmp_path):
    _, sha = make_rules(tmp_path)
    text = build_report([res("SRRTEST", "test", sha)], tmp_path)
    assert sha[:8] in text and sha[:9] not in text


def test_report_shows_mirna_loci_used_per_row(tmp_path):
    _, sha = make_rules(tmp_path)
    used, unused = res("SRRUSED", "test", sha), res("SRRNONE", "test", sha)
    used.meta["mirna_loci_used"] = True
    unused.meta["mirna_loci_used"] = False
    text = build_report([used, unused], tmp_path)
    row = {}
    for l in text.splitlines():  # the classification table comes first in the section
        if l.startswith("| SRR"):
            row.setdefault(l.split("|")[1].strip(), l)
    header = next(l for l in text.splitlines() if "mirna_loci_used" in l)
    col = [c.strip() for c in header.split("|")].index("mirna_loci_used")
    assert row["SRRUSED"].split("|")[col].strip() == "yes"
    assert row["SRRNONE"].split("|")[col].strip() == "no"


def test_report_none_rendered_as_dash(tmp_path):
    _, sha = make_rules(tmp_path)
    text = build_report([res("SRRTEST", "test", sha)], tmp_path)
    assert "None" not in text
    assert "—" in text


def test_report_has_delimitation_and_source_tables(tmp_path):
    _, sha = make_rules(tmp_path)
    text = build_report([res("SRRTEST", "test", sha)], tmp_path)
    assert "mean_recovered_fraction" in text
    assert "tiny<=50" in text and "reads<100" in text
    assert "SrcA" in text
    # one delimitation table per source, not pooled
    assert "Delimitation: SrcA" in text and "Delimitation: SrcB" in text
    sec_a = text.split("Delimitation: SrcA")[1].split("Delimitation: SrcB")[0]
    assert "| 4 | 0.600 | 1.500 | 2 |" in sec_a and "| 9 |" not in sec_a


def test_collect_results_reads_metrics_and_meta(tmp_path):
    d = tmp_path / "results" / "SRRTEST"
    d.mkdir(parents=True)
    (d / "metrics.json").write_text(json.dumps(metrics() | {"accession": "SRRTEST"}))
    (d / "benchmark_run.json").write_text(json.dumps({"role": "test"}))
    (tmp_path / "results" / "empty").mkdir()
    out = collect_results(tmp_path / "results", [lib("SRRTEST", "test"), lib("OTHER", "test")])
    assert [r.lib.accession for r in out] == ["SRRTEST"]
    assert out[0].meta["role"] == "test"
