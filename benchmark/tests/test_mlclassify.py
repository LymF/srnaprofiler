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


def _patch_pipeline(monkeypatch, records):
    monkeypatch.setattr(classify, "load_config", lambda p: {})
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
    from benchmark.ml.models import REGISTRY
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
    monkeypatch.setattr(classify, "RuleScorer",
                        lambda cfg: SimpleNamespace(decide=lambda r: setattr(r, "cls", "siRNA")))
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
        dec.decide_all(recs)


def test_ml_mode_sets_fields(tmp_path):
    import mlclassify
    recs = _records()
    path, schema = _tiny_model(tmp_path, recs)
    mlclassify.load(path).decide_all(recs)
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
    out = tmp_path / "m.joblib"
    sch = train_final.train_and_save(table, "arthropod", "logreg", "main", 0, out)
    saved = json.loads((tmp_path / "m.schema.json").read_text())
    assert saved == sch
    assert saved["features"] == list(mlfeatures.FEATURES)
    assert set(saved) >= {"features", "classes", "threshold", "algorithm", "variant", "seed",
                          "train_commit", "pipeline_hash"}
    assert joblib.load(out).classes_ == saved["classes"]
