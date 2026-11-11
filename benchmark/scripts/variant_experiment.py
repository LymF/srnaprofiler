"""Feature-set experiment for one ML algorithm: main vs bio vs ctx vs bio_ctx.

  main     the production features without ctx_*
  bio      read size, strand, nucleotide, ping-pong and phasing only (VARIANTS["bio"])
  ctx      main plus ctx_gene / ctx_trna / ctx_rrna filled from the genome annotation
           (fraction of the locus covered by CDS+UTR, tRNA, rRNA); ctx_te stays NaN so the
           transposon overlap of the "supported" metric stays independent of the model
  bio_ctx  bio plus the three filled ctx columns

CV: leave-one-unit-out on the tuning libraries, threshold 0.0 (as in the final models).
Arthropod also gets an offline Bombyx test (fit on all tuning rows, predict the Bombyx table,
locus-level metrics like the CV, same code path as fold_metrics).

Run from the repo root:
  python benchmark/scripts/variant_experiment.py <profile> <algorithm> <results_dir> <annotation_dir> <out.json> [n_jobs] [bombyx_rules_results_dir]
"""
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import joblib
import numpy as np

sys.path.insert(0, ".")
from benchmark.ml.cv import VARIANTS, _te_indexes, folds_for, tuning_accessions
from benchmark.ml.dataset import build_table
from benchmark.ml.evaluate import TableFit, decide, fold_metrics
from mlmodels import REGISTRY

GFF = {"Drosophila melanogaster": "Drosophila_melanogaster.BDGP6.54.63.gff3.gz",
       "Arabidopsis thaliana": "Arabidopsis_thaliana.TAIR10.63.gff3.gz",
       "Oryza sativa": "Oryza_sativa.IRGSP-1.0.63.gff3.gz",
       "Zea mays": "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gff3.gz",
       "Bombyx mori": "GCF_030269925.1_ASM3026992v2_genomic.gff.gz"}
CTX_TYPES = {"ctx_gene": {"CDS", "five_prime_UTR", "three_prime_UTR"},
             "ctx_trna": {"tRNA"}, "ctx_rrna": {"rRNA"}}
CTX = list(CTX_TYPES)


def read_intervals(path, alias=None):
    out = {c: defaultdict(list) for c in CTX_TYPES}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            for c, types in CTX_TYPES.items():
                if f[2] in types:
                    chrom = alias.get(f[0], f[0]) if alias else f[0]
                    out[c][chrom].append((int(f[3]) - 1, int(f[4])))
    merged = {}
    for c, d in out.items():
        merged[c] = {}
        for chrom, iv in d.items():
            iv.sort()
            m = [list(iv[0])]
            for a, b in iv[1:]:
                if a <= m[-1][1]:
                    m[-1][1] = max(m[-1][1], b)
                else:
                    m.append([a, b])
            arr = np.array(m)
            merged[c][chrom] = (arr[:, 0], arr[:, 1], np.concatenate([[0], np.cumsum(arr[:, 1] - arr[:, 0])]))
    return merged


def covered_fraction(iv, s, e):
    if iv is None or e <= s:
        return 0.0
    st, en, cum = iv
    lo = int(np.searchsorted(en, s, side="right"))   # first interval ending after s
    hi = int(np.searchsorted(st, e, side="left"))    # intervals starting before e
    if hi <= lo:
        return 0.0
    tot = cum[hi] - cum[lo]
    tot -= max(0, s - st[lo])
    tot -= max(0, en[hi - 1] - e)
    return max(0.0, tot) / (e - s)


def add_ctx(table, annotation_dir, aliases=None):
    table = table.copy()
    for c in CTX:
        table[c] = np.nan
    for org in sorted(set(table["organism"])):
        rows = table["organism"] == org
        merged = read_intervals(Path(annotation_dir) / GFF[org], (aliases or {}).get(org))
        sub = table[rows]
        for c in CTX:
            iv = merged[c]
            table.loc[rows, c] = [covered_fraction(iv.get(str(ch)), int(s), int(e))
                                  for ch, s, e in zip(sub.chrom, sub.start, sub.end)]
        print(org, {c: round(float(table.loc[rows, c].mean()), 4) for c in CTX}, flush=True)
    return table


def variants_for():
    # the table carries the filled ctx_* columns, so "main" must leave them out explicitly
    base = {"main": [f for f in VARIANTS["main"] if not f.startswith("ctx_")],
            "bio": VARIANTS["bio"]}
    base["ctx"] = [f for f in VARIANTS["main"] if f != "ctx_te"]   # ctx_te stays out (independence)
    base["bio_ctx"] = VARIANTS["bio"] + CTX
    return base


def fit_predict(table_train, table_val, feats, algorithm):
    model = TableFit(REGISTRY[algorithm](0, feats)).fit(table_train)
    pred = decide(model, table_val, feats, 0.0).to_numpy(dtype=object)
    return pred, model.p_real(table_val)


def cv_fold(table, algorithm, name, feats, unit, tr, va, te):
    pred, preal = fit_predict(table.loc[tr], table.loc[va], feats, algorithm)
    val = table.loc[va]
    m = fold_metrics(val, pred, preal, val["rule_cls"].to_numpy(dtype=object), None, te)
    return {"variant": name, "unit": unit, "metrics": m}


def main(profile, algorithm, results_dir, annotation_dir, out, n_jobs=6, bombyx_results=None):
    manifest = Path("benchmark/manifest.tsv")
    res = Path(results_dir)
    table = build_table(tuning_accessions(manifest, profile), res, manifest)
    table = table[table["profile"] == profile].reset_index(drop=True)
    table = add_ctx(table, annotation_dir)
    te = _te_indexes(tuple(sorted(set(table["organism"]))))
    folds = [f for f in folds_for(table, profile) if not str(f[0]).startswith("xspecies")]
    vs = variants_for()
    jobs = [joblib.delayed(cv_fold)(table, algorithm, n, f, u, tr, va, te)
            for n, f in vs.items() for u, tr, va in folds]
    cv = joblib.Parallel(n_jobs=int(n_jobs))(jobs)
    result = {"cv": cv, "bombyx": {}}
    json.dump(result, open(out, "w"))     # keep the CV even if the Bombyx step fails
    if profile == "arthropod" and bombyx_results:
        alias = {}
        for line in open(Path(annotation_dir) / "bmor_refseq_to_genome.tsv"):
            a = line.rstrip("\n").split("\t")
            if len(a) == 2:
                alias[a[0]] = a[1]
        bm = build_table(["DRR377498"], Path(bombyx_results), manifest, final_test=True).reset_index(drop=True)
        bm = add_ctx(bm, annotation_dir, {"Bombyx mori": alias})
        rule = bm["rule_cls"].to_numpy(dtype=object)
        result["bombyx"]["rules"] = fold_metrics(bm, rule, np.zeros(len(bm)), rule, None, None)
        for n, f in vs.items():
            pred, preal = fit_predict(table, bm, f, algorithm)
            result["bombyx"][n] = fold_metrics(bm, pred, preal, rule, None, None)
    json.dump(result, open(out, "w"))
    print("wrote", out)


if __name__ == "__main__":
    main(*sys.argv[1:8])
