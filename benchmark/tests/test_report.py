import hashlib
import json
from pathlib import Path

from benchmark.manifest import Library
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
            "summary": {"n_clusters": 4, "mean_recovered_fraction": 0.6, "mean_fragments": 1.5},
            "merged_distinct": 2,
            "tiny_counts": {"50": 3, "100": 5},
            "low_read_counts": {"20": 8, "50": 9, "100": 10},
        },
        "by_source": {"SrcA": {"piRNA": {"precision": 0.8, "fdr": 0.2, "recall": 0.4, "tp": 1, "fp": 1}}},
        "n_loci": 10,
    }


def make_rules(tmp_path, text="a: 1\n"):
    rules = tmp_path / "rules.yaml"
    rules.write_text(text)
    return rules, hashlib.sha256(rules.read_bytes()).hexdigest()


def res(acc, role, sha, organism="Bombyx mori", **kw):
    meta = {"rules_sha256": sha, "rules_path": "/elsewhere/rules.yaml", "role": role}
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


def test_collect_results_reads_metrics_and_meta(tmp_path):
    d = tmp_path / "results" / "SRRTEST"
    d.mkdir(parents=True)
    (d / "metrics.json").write_text(json.dumps(metrics() | {"accession": "SRRTEST"}))
    (d / "benchmark_run.json").write_text(json.dumps({"role": "test"}))
    (tmp_path / "results" / "empty").mkdir()
    out = collect_results(tmp_path / "results", [lib("SRRTEST", "test"), lib("OTHER", "test")])
    assert [r.lib.accession for r in out] == ["SRRTEST"]
    assert out[0].meta["role"] == "test"
