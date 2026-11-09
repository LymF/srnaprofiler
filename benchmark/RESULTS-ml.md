# ML classifier: one-shot test on held-out maize

One run of the plant random forest on the four held-out maize libraries, compared with the rules classifier on the same loci. Selection of the model: `benchmark/ML-CV.md` (tuning libraries only). Format follows `benchmark/RESULTS-run6.md`. Numbers are reported as obtained; nothing in the model, the code or the metrics was changed after seeing them.

## 1. What was tested

- Model: plant random forest (Breiman, 2001), variant `main`, seed 0, chosen by tuning-only cross-validation (`ML-CV.md` section 4: macro-F1 0.753 ± 0.123 against 0.574 ± 0.089 for the rules, p = 5.56e-07, eligible). Classes: hcsiRNA, phasiRNA, phasiRNA24; threshold 0.0 (the model never abstains).
- Run: `--classifier ml`, once, on the four maize libraries SRR15530410, SRR15530411, SRR15530418 and SRR15530419. The user approved the run on 2026-11-08.
- Bombyx was **not** run: no arthropod model passed selection (`ML-CV.md` section 5, winner: none).
- The miRNA class is not an ML output: miRNA loci come only from miRDeep-P2 (ruling R23, user decision 2026-11-07), so the miRNA rows are identical by construction.

## 2. Provenance

| item | value |
|---|---|
| model file | `models/plant.joblib` on the server, sha256 `f6780f1c3716772d8ff46ac937b433f233b8f6c70f5b248261532fbace0708f9` |
| model schema | algorithm rf, variant main, seed 0; 65 features; classes hcsiRNA, phasiRNA, phasiRNA24; threshold 0.0 |
| training code | commit 44c190d (`feat(ml): classifier never outputs miRNA (miRNA only from miRDeep2)`) |
| training pipeline hash | `99724418152630e93b82f6170abb7c7f28194b1f4bbc324e2175297f82d3b436` |
| upstream outputs | `results_ml` reuses the run #7 qc, mapping and miRDeep-P2 outputs; only the classify/ML step was run (ruling R25) |

On R25: the copied `benchmark_run.json` has its `pipeline_sha256` rewritten to the HEAD hash, and the original is kept as `reused_from_pipeline_sha256`, so the file carries both hashes. The justification is that `git diff 99e849e..HEAD` touches only classify and ML code (`classify.py`, `mlclassify.py`, `mlmodels.py`, `records.p_real`, the reporter column, the yaml `ml` key, the CLI flags). If that were wrong, the ML numbers would rest on stale upstream outputs. The 94 MB model file is not in git (ruling R24); only `schema.json` and the sha256 are.

## 3. Results per library

Values are rules -> ML. P precision, R recall, "called" the number of loci called in the class, "unannotated" the called loci that match no truth region. "-" means not defined (no truth) or no call; "none" means the classifier made no call.

### SRR15530410

| class | P | R | called | unannotated |
|---|---|---|---|---|
| hcsiRNA | 0.000 -> 0.000 | - -> - | 3173 -> 3901 | 2893 -> 3708 |
| miRNA | 1.000 -> 1.000 | 0.186 -> 0.186 | 42 -> 42 | 11 -> 11 |
| phasiRNA | 0.975 -> 0.974 | 0.734 -> 0.787 | 1034 -> 2709 | 481 -> 2064 |
| phasiRNA24 | 1.000 -> 0.990 | 0.125 -> 0.317 | 174 -> 2127 | 96 -> 1935 |
| siRNA | 0.000 -> - | - -> - | 429 -> none | 369 -> none |

### SRR15530411

| class | P | R | called | unannotated |
|---|---|---|---|---|
| hcsiRNA | 0.000 -> 0.000 | - -> - | 3888 -> 6456 | 3703 -> 6386 |
| miRNA | 1.000 -> 1.000 | 0.419 -> 0.419 | 83 -> 83 | 12 -> 12 |
| phasiRNA | 0.956 -> 0.940 | 0.809 -> 0.811 | 1443 -> 3087 | 918 -> 2537 |
| phasiRNA24 | 1.000 -> 0.995 | 0.096 -> 0.369 | 280 -> 3175 | 232 -> 2989 |
| siRNA | 0.000 -> - | - -> - | 1006 -> none | 995 -> none |

### SRR15530418

| class | P | R | called | unannotated |
|---|---|---|---|---|
| hcsiRNA | 0.000 -> 0.000 | - -> - | 5850 -> 6714 | 5575 -> 6501 |
| miRNA | 1.000 -> 1.000 | 0.198 -> 0.198 | 46 -> 46 | 12 -> 12 |
| phasiRNA | 0.976 -> 0.965 | 0.739 -> 0.768 | 1160 -> 2428 | 622 -> 1822 |
| phasiRNA24 | 1.000 -> 1.000 | 0.079 -> 0.254 | 218 -> 2055 | 177 -> 1920 |
| siRNA | 0.000 -> - | - -> - | 404 -> none | 354 -> none |

### SRR15530419

| class | P | R | called | unannotated |
|---|---|---|---|---|
| hcsiRNA | 0.000 -> 0.000 | - -> - | 3550 -> 3773 | 3259 -> 3570 |
| miRNA | 1.000 -> 1.000 | 0.134 -> 0.134 | 36 -> 36 | 13 -> 13 |
| phasiRNA | 0.977 -> 0.983 | 0.697 -> 0.727 | 905 -> 2382 | 388 -> 1798 |
| phasiRNA24 | 1.000 -> 0.988 | 0.090 -> 0.300 | 148 -> 1637 | 98 -> 1471 |
| siRNA | 0.000 -> - | - -> - | 309 -> none | 258 -> none |

## 4. Summary

- **phasiRNA**: precision about unchanged (rules 0.956-0.977, ML 0.940-0.983; the largest drop is 0.016, in SRR15530411) and recall slightly up (+0.002 to +0.053 across the four libraries).
- **phasiRNA24**: recall about 2.5 to 3.8 times higher (0.079-0.125 -> 0.254-0.369), at precision 0.988-1.000 (rules 1.000 in all four).
- **miRNA**: identical in every library, as expected (miRDeep-P2, ruling R23). Recall is 0.134-0.419 and says nothing about the ML.
- **hcsiRNA**: maize has no hcsiRNA truth, so precision 0.000 in both columns is not informative; every call counts as a false positive against another annotated class or as unannotated. The same artefact is documented for the rice folds in `ML-CV.md` section 8. Recall is undefined. Only the call counts are informative (3173-5850 -> 3773-6714).
- **siRNA**: the ML never calls siRNA, because the plant tuning data has no siRNA labels to learn from. The rules called 309-1006 siRNA loci per library, of which 258-995 were unannotated, with precision 0.000 in all four.

## 5. Caveat that governs the reading of the table

The ML calls many more loci than the rules, almost all of them on unannotated loci:

| class | called, ML / rules (four libraries) | share of the extra calls that are unannotated |
|---|---|---|
| phasiRNA | 2.1x to 2.6x (2382-3087 against 905-1443) | 94% to 98% |
| phasiRNA24 | 9.4x to 12.2x (1637-3175 against 148-280) | 92% to 95% |

Precision and recall are measured only on annotated loci. The extra calls fall outside the annotation, so **this table does not validate them**. They could be real phasiRNA loci missing from the annotation or false positives, and the table cannot tell which. A biological audit (phasing evidence on the extra loci) is pending.

Maize is now spent as an independent test. Any further change to the ML (features, model, threshold, class set) needs new held-out species before it can be evaluated again; reusing maize would make the test no longer independent.

## 6. Context: why annotation-based metrics are unreliable here

The audit on the fly tuning data (out-of-fold, rf main, seed 0) shows the same pattern in a second species. 24369 ML-piRNA loci with no annotation carry the piRNA signature (24-31 nt fraction 0.88, 1U 0.74, ping-pong z 7.3, TE context), so the annotation is incomplete and an annotation-based precision understates the model. In the other direction, the piRNA precision loss comes mostly from 109 miRNA-like loci (22 nt, 1U 0.94, no ping-pong or TE) that are forced into piRNA because the ML cannot output miRNA (R23); this points to the need for an abstain or "other" class. The threshold is 0 in every fold, so the model never abstains. Finally, `p_real` is unreliable: it falls below 0.5 on 17246 piRNA-like unannotated loci, which violates the assumption of positive-unlabeled learning that labelled loci are a random sample of the real loci (Elkan & Noto, 2008; Bekker & Davis, 2020). Annotated loci are a biased sample, so a score trained to separate them from the rest inherits that bias. The maize result should be read with the same limits: annotation-based precision may understate the ML on unannotated loci, and may also hide false positives there.

## References

- Bekker, J., & Davis, J. (2020). Learning from positive and unlabeled data: a survey. *Machine Learning*, 109(4), 719–760. https://doi.org/10.1007/s10994-020-05877-5
- Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
- Elkan, C., & Noto, K. (2008). Learning classifiers from only positive and unlabeled data. In *Proceedings of the 14th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 213–220). https://doi.org/10.1145/1401890.1401920
