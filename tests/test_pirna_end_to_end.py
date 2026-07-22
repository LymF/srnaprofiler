"""FIX 1 (+FIX 4) integração: locus de piRNA ~27 nt classificado corretamente.

A janela de MEDIÇÃO do extractor tem que ser superset das janelas de classe
(piRNA usa 24-31). Com dicermax=24 fixo, todo piRNA de 27 nt virava
dicer_call="N" e era forçado a review=True. Aqui um locus de 27 nt, 1U alto,
bifita com ping-pong (overlap 5'-5' = 10) deve sair como piRNA sem review.
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
    # o locus mais abundante é o de piRNA
    row = max(rows, key=lambda r: int(r[ci["reads"]]))
    assert row[ci["dicer_call"]] == "27"
    assert row[ci["cls"]] == "piRNA"
    assert row[ci["review"]] == "False"
