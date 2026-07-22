"""FIX 1 (+FIX 4) integration: ~27 nt piRNA locus correctly classified.

The extractor's MEASUREMENT window must be a superset of the class windows
(piRNA uses 24-31). With dicermax=24 fixed, every 27 nt piRNA turned into
dicer_call="N" and was forced to review=True. Here a 27 nt locus, high 1U,
double-stranded with ping-pong (5'-5' overlap = 10) must come out as piRNA
with no review.
"""
import types
from clustering import cluster_regions
from features import extract
from decider import RuleScorer, load_config


def test_pirna_locus_classified_without_review():
    config = load_config("rules.yaml")
    ext_cfg = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}
    ext_cfg.update(config.get("extractor", {}))

    regions = cluster_regions("tests/data/pirna_locus.sam", mincov=1, pad=5)
    recs = extract("tests/data/pirna_locus.sam", regions, lib_size=1_000_000, config=ext_cfg)
    rec = max(recs, key=lambda r: r.reads)

    scorer = RuleScorer(config)
    scorer.decide(rec)

    assert rec.dicer_call == "27"
    assert rec.cls == "piRNA"
    assert rec.review is False


def test_pirna_via_run_classify(tmp_path):
    from classify import run_classify
    args = types.SimpleNamespace(
        sam="tests/data/pirna_locus.sam", genome=None, prefix="pi",
        outdir=str(tmp_path), mincov=1, pad=5, mirna_loci=None,
        config="rules.yaml", norm=1_000_000)
    run_classify(args)
    tsv = (tmp_path / "pi.clusters.tsv").read_text().splitlines()
    header, rows = tsv[0].split("\t"), [l.split("\t") for l in tsv[1:]]
    ci = {c: i for i, c in enumerate(header)}
    # the most abundant locus is the piRNA one
    row = max(rows, key=lambda r: float(r[ci["reads"]]))
    assert row[ci["dicer_call"]] == "27"
    assert row[ci["cls"]] == "piRNA"
    assert row[ci["review"]] == "False"
