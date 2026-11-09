"""--classifier ml: optional mode; rules stay the default and need no scikit-learn."""
import json
import sys
from collections import namedtuple
from types import SimpleNamespace

import pytest

import classify
import features
import mlfeatures
from records import ClusterRecord

CAln = namedtuple("CAln", "chrom pos seq nh strand")
CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _records(n=24):
    out = []
    for i in range(n):
        length = (21, 24, 26)[i % 3]
        alns = [CAln("chr1", 1000 * i + 30 * k + (i % 5), "TA" * (length // 2) + "A" * (length % 2),
                     1, "+" if (k + i) % 4 else "-") for k in range(8 + i)]
        r = features.build_record(("chr1", 1000 * i, 1000 * i + 600), alns, 1e6, CFG)
        r.cls = ("siRNA", "piRNA", "ambiguous")[i % 3]
        out.append(r)
    out[0].cls = "miRNA"
    out[0].confidence = 1.0
    return out


def _scorer():
    from pathlib import Path
    from decider import RuleScorer, load_config
    return RuleScorer(load_config(Path(classify.__file__).resolve().parent / "rules.yaml"))


def _patch_pipeline(monkeypatch, records):
    from pathlib import Path
    from decider import load_config
    cfg = load_config(Path(classify.__file__).resolve().parent / "rules.yaml")
    monkeypatch.setattr(classify, "load_config", lambda p: cfg)
    monkeypatch.setattr(classify, "class_size_windows", lambda cfg: [])
    monkeypatch.setattr(classify, "ensure_sorted_indexed_bam", lambda *a: "x.bam")
    monkeypatch.setattr(classify, "prefers_xm_as_nh", lambda b: False)
    monkeypatch.setattr(classify, "mapped_chromosomes", lambda b: ["chr1"])
    monkeypatch.setattr(classify, "reference_alignment_counts", lambda b: [4])
    monkeypatch.setattr(classify, "available_memory_mb", lambda: 10 ** 6)
    monkeypatch.setattr(classify, "extract_by_chrom", lambda *a, **k: (records, []))


def _args(tmp_path, **kw):
    base = dict(outdir=str(tmp_path), config=None, prefix="t", sam="x.sam", norm=1e6,
                genome="g.fa", mincov=1, mirna_loci=None, organelle_contigs=None,
                threads=1, pad=50)
    base.update(kw)
    return SimpleNamespace(**base)


def _tiny_model(tmp_path, records, classes=("siRNA", "piRNA")):
    """Fit a logreg on the records' own feature table; return the .joblib path."""
    import joblib
    import pandas as pd
    from mlmodels import REGISTRY
    rows = [mlfeatures.derive({n: getattr(r, n, None) for n in
                               [*mlfeatures.FEATURES, "dicer_call"]}) for r in records]
    table = pd.DataFrame(rows)
    table["rule_cls"] = "ambiguous"
    table["carved"] = False
    y = pd.Series([classes[i % len(classes)] for i in range(len(table))])
    model = REGISTRY["logreg"](0, list(mlfeatures.FEATURES)).fit(table, y)
    path = tmp_path / "tiny.joblib"
    joblib.dump(model, path)
    schema = {"features": list(mlfeatures.FEATURES), "classes": list(model.classes_),
              "threshold": 0.5, "algorithm": "logreg", "variant": "main", "seed": 0,
              "train_commit": "unknown", "pipeline_hash": "x"}
    (tmp_path / "tiny.schema.json").write_text(json.dumps(schema))
    return path, schema


def test_rules_mode_without_sklearn(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "sklearn", None)
    monkeypatch.setitem(sys.modules, "joblib", None)
    monkeypatch.delitem(sys.modules, "mlclassify", raising=False)
    recs = _records(6)
    _patch_pipeline(monkeypatch, recs)
    real = classify.RuleScorer

    class Stub(real):
        def __init__(self, cfg):
            pass

        def decide(self, r):
            r.cls = "siRNA"

        def gate_flags(self, r):
            return {}

    monkeypatch.setattr(classify, "RuleScorer", Stub)
    classify.run_classify(_args(tmp_path))
    lines = (tmp_path / "t.clusters.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    assert header.index("p_real") == header.index("score") + 1
    assert {l.split("\t")[header.index("p_real")] for l in lines[1:]} == {"NA"}
    assert "mlclassify" not in sys.modules


def test_ml_mode_schema_mismatch(tmp_path):
    import mlclassify
    recs = _records()
    path, schema = _tiny_model(tmp_path, recs)
    schema["features"] = [*schema["features"], "bogus_feature"]
    (tmp_path / "tiny.schema.json").write_text(json.dumps(schema))
    dec = mlclassify.load(path)
    with pytest.raises(ValueError, match="bogus_feature"):
        dec.decide_all(recs, _scorer())


def test_ml_mode_sets_fields(tmp_path):
    import mlclassify
    recs = _records()
    path, schema = _tiny_model(tmp_path, recs)
    mlclassify.load(path).decide_all(recs, _scorer())
    assert recs[0].cls == "miRNA" and recs[0].p_real is None
    for r in recs[1:]:
        assert r.cls in [*schema["classes"], "ambiguous"]
        assert 0.0 <= r.p_real <= 1.0
        assert r.evidence == "ml:logreg"
        assert r.review == (r.cls == "ambiguous")


def test_ml_mode_end_to_end_writes_p_real(monkeypatch, tmp_path):
    recs = _records()
    path, _ = _tiny_model(tmp_path, recs)
    _patch_pipeline(monkeypatch, recs)
    classify.run_classify(_args(tmp_path, classifier="ml", ml_model=str(path)))
    lines = (tmp_path / "t.clusters.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    vals = [l.split("\t")[header.index("p_real")] for l in lines[1:]]
    assert vals[0] == "NA" and all(0 <= float(v) <= 1 for v in vals[1:])


def test_default_is_rules():
    import srnaprofiler
    args = srnaprofiler.build_parser().parse_args(["classify", "-p", "a", "-o", "o", "--sam", "x"])
    assert args.classifier == "rules" and args.ml_model is None


def test_resolve_model_path_default_from_config():
    from pathlib import Path
    p = classify.resolve_ml_model(None, {"ml": {"model": "models/arthropod.joblib"}})
    assert p == Path(classify.__file__).resolve().parent / "models" / "arthropod.joblib"
    assert classify.resolve_ml_model("/a/b.joblib", {}) == Path("/a/b.joblib")
    with pytest.raises(ValueError):
        classify.resolve_ml_model(None, {})


def test_train_final_writes_model_and_schema(tmp_path):
    import joblib
    import pandas as pd
    from benchmark.ml import train_final
    recs = _records(60)
    rows = [mlfeatures.derive({n: getattr(r, n, None) for n in [*mlfeatures.FEATURES, "dicer_call"]})
            for r in recs]
    table = pd.DataFrame(rows)
    table["rule_cls"] = "ambiguous"
    table["carved"] = False
    table["acc"] = "A"
    table["organism"] = "Bombyx mori"
    table["chrom"] = ["c%d" % (i % 6) for i in range(len(table))]
    table["label"] = [("siRNA", "piRNA", None)[i % 3] for i in range(len(table))]
    table["gate_siRNA"] = table["gate_piRNA"] = 1.0
    out = tmp_path / "m.joblib"
    sch = train_final.train_and_save(table, "arthropod", "logreg", "main", 0, out)
    saved = json.loads((tmp_path / "m.schema.json").read_text())
    assert saved == sch
    assert saved["features"] == list(mlfeatures.FEATURES)
    assert set(saved) >= {"features", "classes", "threshold", "algorithm", "variant", "seed",
                          "train_commit", "pipeline_hash"}
    assert joblib.load(out).classes_ == saved["classes"]


def test_loaded_model_needs_no_benchmark_module(tmp_path):
    import subprocess
    from pathlib import Path
    recs = _records()
    path, _ = _tiny_model(tmp_path, recs)
    repo = Path(classify.__file__).resolve().parent
    code = (
        "import sys, mlclassify\n"
        "from types import SimpleNamespace as N\n"
        "import mlfeatures\n"
        "d = mlclassify.load(sys.argv[1])\n"
        "names = [n for n in mlfeatures.FEATURES if n != 'dicer_len'] + ['dicer_call']\n"
        "r = N(cls='', confidence=0, review=False, p_real=None, evidence='', score=0, class_scores={},\n"
        "      **{n: 1.0 for n in names if n != 'dicer_call'}, dicer_call='21')\n"
        "sc = N(gated_choice=lambda p, rec: max(p, key=p.get))\n"
        "d.decide_all([r], sc)\n"
        "bad = [m for m in sys.modules if m.split('.')[0] == 'benchmark']\n"
        "assert not bad, bad\n"
        "assert r.p_real is not None\n")
    res = subprocess.run([sys.executable, "-c", code, str(path)], cwd=repo,
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr


def test_loaded_model_refuses_refit(tmp_path):
    import joblib
    path, _ = _tiny_model(tmp_path, _records())
    model = joblib.load(path)
    with pytest.raises(RuntimeError, match="loaded from an artifact"):
        model.fit(None, None)


def test_incompatible_artifact_is_runtime_error(tmp_path, monkeypatch):
    import mlclassify
    path, _ = _tiny_model(tmp_path, _records())
    import joblib
    (tmp_path / "gonemod.py").write_text("class Cls:\n    pass\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    import gonemod
    joblib.dump(gonemod.Cls(), path)
    del sys.modules["gonemod"]
    (tmp_path / "gonemod.py").unlink()
    with pytest.raises(RuntimeError, match="incompatible mlmodels"):
        mlclassify.load(path)


def test_records_missing_attribute_named(tmp_path):
    import mlclassify
    recs = _records()
    path, _ = _tiny_model(tmp_path, recs)
    bare = [SimpleNamespace(cls="", reads=1.0)]
    with pytest.raises(ValueError, match="pingpong_z"):
        mlclassify.load(path).decide_all(bare, _scorer())


@pytest.mark.parametrize("cfg", ["rules.yaml", "rules.plant.yaml"])
def test_shipped_configs_with_ml_key_still_score(cfg):
    from pathlib import Path
    from decider import RuleScorer, load_config
    config = load_config(Path(classify.__file__).resolve().parent / cfg)
    assert config["ml"]["model"].startswith("models/")
    r = _records(3)[1]
    RuleScorer(config).decide(r)
    assert r.cls


def test_production_decide_all_matches_benchmark_decide(tmp_path):
    """Same trained model, same rows: MLDecider on records == evaluate.decide on gate columns."""
    import mlclassify
    import pandas as pd
    from benchmark.ml import evaluate
    recs = _records()
    path, schema = _tiny_model(tmp_path, recs, classes=("siRNA", "piRNA", "other"))
    scorer = _scorer()
    for r in recs[1:]:
        scorer.decide(r)               # sets the evidence the gate reads
    flags = [scorer.gate_flags(r) for r in recs]
    dec = mlclassify.load(path)
    rows = [mlfeatures.derive({n: getattr(r, n, None) for n in
                               [*mlfeatures.FEATURES, "dicer_call"]}) for r in recs[1:]]
    table = pd.DataFrame(rows)
    table["carved"] = False
    for c in scorer.enabled_classes(scorer.cfg):
        table[f"gate_{c}"] = [float(f[c]) for f in flags[1:]]
    want = evaluate.decide(dec.model, table, [], dec.threshold).tolist()
    dec.decide_all(recs, scorer)
    assert [r.cls for r in recs[1:]] == want
    assert not {"other", "miRNA"} & set(want)
