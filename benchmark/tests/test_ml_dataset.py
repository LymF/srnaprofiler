import math

import pytest

import mlfeatures
from benchmark.ml import dataset
from benchmark.ml.dataset import FEATURES, build_table
from benchmark.truth import Region
from reporter import _COLS

MANIFEST_HEADER = ["accession", "organism", "role", "taxon_group", "tissue", "replicate",
                   "genome_url", "truth_sources"]


def _write_clusters(path, loci):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(_COLS) + "\n")
        for chrom, start, end, cls in loci:
            row = {c: "NA" for c in _COLS}
            row.update(chrom=chrom, start=start, end=end, name=f"{chrom}:{start}", reads=10,
                       dicer_call="21", cls=cls, confidence="0.9", score="3", pct_1u="0.5")
            fh.write("\t".join(str(row[c]) for c in _COLS) + "\n")


@pytest.fixture
def env(tmp_path, monkeypatch):
    manifest = tmp_path / "manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("\t".join(MANIFEST_HEADER) + "\n")
        fh.write("\t".join(["T1", "Dmel", "tuning", "arthropod", "ovary", "1", "u", "x.gff"]) + "\n")
        fh.write("\t".join(["X1", "Bmor", "test", "arthropod", "ovary", "1", "u", "x.gff"]) + "\n")
    results = tmp_path / "results"
    _write_clusters(results / "T1" / "work" / "T1.clusters.tsv",
                    [("chr1", 1000, 1600, "piRNA"), ("chr1", 50000, 50200, "miRNA")])
    _write_clusters(results / "X1" / "work" / "X1.clusters.tsv", [("chr1", 0, 100, "siRNA")])
    monkeypatch.setattr(dataset, "load_truth",
                        lambda lib: [Region("chr1", 500, 1500, "piRNA", "src", "r1")])
    return results, manifest


def test_build_table_labels_by_match(env):
    results, manifest = env
    df = build_table(["T1"], results, manifest)
    assert len(df) == 2
    assert df.loc[0, "label"] == "piRNA"
    assert df.loc[1, "label"] is None
    assert df.loc[0, "profile"] == "arthropod"
    assert df.loc[0, "dicer_len"] == 21.0
    assert df.loc[0, "score"] == 3.0 and math.isnan(df.loc[0, "repro_frac"])


def test_build_table_refuses_test_library(env):
    results, manifest = env
    with pytest.raises(ValueError):
        build_table(["X1"], results, manifest)
    assert len(build_table(["X1"], results, manifest, final_test=True)) == 1


def test_build_table_features_exclude_rule_outputs():
    assert set(FEATURES) & {"cls", "confidence", "score", "evidence"} == set()
    assert not any(f.startswith("s_") for f in FEATURES)


def test_carved_flag(env):
    results, manifest = env
    df = build_table(["T1"], results, manifest)
    assert list(df["carved"]) == [False, True]


def test_derive_dicer_len_and_na():
    row = {c: "NA" for c in _COLS}
    row.update(dicer_call="21", reads="5", gc=0.4)
    d = mlfeatures.derive(row)
    assert d["dicer_len"] == 21.0 and d["reads"] == 5.0 and d["gc"] == 0.4
    assert math.isnan(d["pct_1u"]) and math.isnan(d["szp_21"])
    assert all(math.isnan(d[k]) for k in ("ctx_te", "ctx_trna", "ctx_rrna", "ctx_gene"))
    row["dicer_call"] = "N"
    assert math.isnan(mlfeatures.derive(row)["dicer_len"])
    assert list(d) == mlfeatures.FEATURES and len(FEATURES) == 22 + 1 + 36 + 3 + 4


def test_build_table_index_unique_across_libraries(env):
    results, manifest = env
    with open(manifest, "a") as fh:
        fh.write("\t".join(["T2", "Dmel", "tuning", "arthropod", "ovary", "2", "u", "x.gff"]) + "\n")
    _write_clusters(results / "T2" / "work" / "T2.clusters.tsv",
                    [("chr2", 10, 200, "siRNA"), ("chr2", 900, 1200, "siRNA")])
    df = build_table(["T1", "T2"], results, manifest)
    assert len(df) == 4 and df.index.is_unique
    assert list(df.index) == list(range(4))
