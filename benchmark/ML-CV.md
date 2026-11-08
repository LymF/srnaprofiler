# ML classifier: cross-validation on tuning data

Cross-validation of 12 machine-learning algorithms against the rules classifier, on the tuning libraries only, for the two rules profiles (plant, arthropod). Design and success criterion: `docs/superpowers/specs/2026-11-07-ml-classifier-design.md` (sections 6 to 8). The test libraries (Bombyx, maize) were not used for anything below. Format follows `benchmark/RESULTS-run6.md`.

## 1. Provenance

| item | value |
|---|---|
| code | commit 44c190d (branch feat/benchmark) |
| tuning libraries | fly SRR30544370 / SRR30544371 / SRR30544372; Arabidopsis SRR24236124 / SRR24236125 / SRR24236126; rice DRR076318 / DRR076319 / DRR076320 |
| test libraries | Bombyx and maize: **not used** (held out; no test evaluation has been run) |
| pipeline outputs | run #7; metrics equal run #6 within 0.001 for every tuning library (the only change in run #7 is a miRNA locus count on a Bombyx test library, DRR377498 `n_called` 100 to 99, no rate moved by more than 0.001) |
| cross-validation | chromosome-grouped leave-one-unit-out; a unit is a chromosome holding at least 5% of a species' loci, the rest pooled into one unit; the unit is held out in all replicates, so a locus never sits in train and validation through another replicate |
| seeds | 5 |
| algorithms | 12: k-means, Gaussian mixture, HDBSCAN, hierarchical (unsupervised); label propagation, label spreading, self-training (semi-supervised); logistic regression, random forest, XGBoost, LightGBM, SVM (supervised); each in variants `main`, `stack`, `repro` |
| plant | 18 chromosome folds plus 2 cross-species folds (Arabidopsis to rice and rice to Arabidopsis), reported separately and not used for selection; 3620 fold lines, 0 errors |
| arthropod | 7 units (fly only; piRNA has one tuning species); 1267 fold lines, 0 errors |

Implementations: scikit-learn (Pedregosa et al., 2011), XGBoost (Chen & Guestrin, 2016), LightGBM (Ke et al., 2017); random forest is Breiman (2001). Grouping folds by chromosome follows the argument that random splits overstate performance on structured data (Roberts et al., 2017). The four O(n^2) methods (hierarchical, label propagation, label spreading, SVM) fit on a seeded subsample capped at 20000 rows, labelled rows always kept, so they are slightly disadvantaged. The "is real" score uses positive-unlabeled learning (annotated loci positive, the rest unlabeled; Elkan & Noto, 2008; Bekker & Davis, 2020), with each algorithm's own estimator family as base.

## 2. Labelled loci

A locus takes the class of the truth region it matches under the benchmark rule (reciprocal overlap at least 0.5); unmatched loci are unlabeled. Approximate labelled loci per library (from the run ledger, not recounted here):

| species | class | labelled loci per library |
|---|---|---|
| fly | piRNA | about 2160 |
| fly | miRNA | about 115 |
| fly | siRNA | about 33 |
| Arabidopsis | hcsiRNA | about 19200 |
| Arabidopsis | miRNA | about 150 |
| Arabidopsis | phasiRNA | about 82 |
| Arabidopsis | phasiRNA24 | about 38 |
| rice | phasiRNA | about 2340 |
| rice | phasiRNA24 | about 440 |
| rice | miRNA | about 78 |

The three replicates of a species are not independent. Caveats on the labels:

- Arabidopsis phasiRNA24 (PHAS24) labels are sRNAanno algorithmic predictions, in a species without DCL5; they are noisy truth.
- hcsiRNA is labelled only in Arabidopsis; the rice hc-siRNA file was excluded as truth.
- The user decided to keep both as they are (2026-11-07). Annotated truth is used as label as is, as in the benchmark.

## 3. Selection rule (spec section 7.3), as applied

Fixed before any run. On tuning data, under the cross-validation above:

- no class loses more than 0.01 precision or recall versus the rules;
- the share of ambiguous loci falls;
- on both "is real" proxies (reproducibility across replicates, negative-control rate over rRNA and protein-coding exons), the model does better than the rules.

Among algorithms that pass, the winner has the best mean macro-F1 over classes on annotated loci across folds and seeds; if two are within the fold-to-fold variation, the simpler wins. The comparison uses the corrected resampled t-test (Nadeau & Bengio, 2003) on per-fold, seed-averaged macro-F1. Eligibility is judged on variant `main`. Per-class metrics and the 0.01 constraint are computed only where the class has truth in the fold. Cross-species folds are excluded from eligibility, ranking and the t-test.

Per-class values are locus-level (predicted class versus label), not the region-level recall of `RESULTS-run6.md`, so they differ in definition from that report. miRNA loci come from miRDeep2 and are identical across algorithms (section 7), so the miRNA columns say nothing about the ML.

## 4. Plant results (variant main)

macro_f1, ambiguous_share and control_real_rate are mean ± sd over folds x seeds (n = folds x seeds; rules n = folds). Per-class values are means. "p vs rules" is the corrected resampled t-test on macro_f1.

| algorithm | n | macro_f1 | ambiguous_share | control_real_rate | p vs rules | eligible |
|---|---|---|---|---|---|---|
| rules | 18 | 0.574 ± 0.089 | 0.387 ± 0.248 | 0.326 ± 0.272 | - | baseline |
| gmm | 90 | 0.465 ± 0.065 | 0.000 ± 0.000 | 0.203 ± 0.259 | 3.38e-05 | no |
| hdbscan | 90 | 0.513 ± 0.050 | 0.013 ± 0.079 | 0.265 ± 0.239 | 0.0364 | no |
| hier | 90 | 0.515 ± 0.049 | 0.000 ± 0.000 | 0.203 ± 0.259 | 0.0181 | no |
| kmeans | 90 | 0.517 ± 0.047 | 0.000 ± 0.000 | 0.203 ± 0.259 | 0.047 | no |
| labelprop | 90 | 0.723 ± 0.124 | 0.000 ± 0.000 | 0.206 ± 0.250 | 8.54e-06 | yes |
| labelspread | 90 | 0.723 ± 0.124 | 0.000 ± 0.000 | 0.209 ± 0.252 | 7.95e-06 | yes |
| lgbm | 90 | 0.760 ± 0.130 | 0.000 ± 0.000 | 0.180 ± 0.248 | 6.03e-07 | no |
| logreg | 90 | 0.746 ± 0.139 | 0.000 ± 0.000 | 0.203 ± 0.259 | 1.14e-05 | no |
| rf | 90 | 0.753 ± 0.123 | 0.000 ± 0.000 | 0.184 ± 0.248 | 5.56e-07 | yes |
| selftrain | 90 | 0.752 ± 0.130 | 0.000 ± 0.000 | 0.184 ± 0.250 | 1.2e-06 | yes |
| svm | 90 | 0.746 ± 0.129 | 0.000 ± 0.000 | 0.222 ± 0.275 | 4.19e-06 | yes |
| xgb | 90 | 0.765 ± 0.128 | 0.000 ± 0.000 | 0.172 ± 0.238 | 4.1e-07 | no |

Per-class precision / recall (P / R):

| algorithm | hcsiRNA P / R | miRNA P / R | phasiRNA P / R | phasiRNA24 P / R |
|---|---|---|---|---|
| rules | 0.293 / 0.975 | 0.982 / 0.446 | 0.752 / 0.866 | 0.706 / 0.124 |
| gmm | 0.292 / 0.988 | 0.982 / 0.446 | 0.718 / 0.636 | n/a / 0.000 |
| hdbscan | 0.292 / 0.999 | 0.982 / 0.446 | 0.860 / 0.752 | 1.000 / 0.003 |
| hier | 0.292 / 0.998 | 0.982 / 0.446 | 0.789 / 0.851 | n/a / 0.000 |
| kmeans | 0.292 / 1.000 | 0.982 / 0.446 | 0.794 / 0.853 | n/a / 0.000 |
| labelprop | 0.293 / 0.999 | 0.982 / 0.446 | 0.817 / 0.929 | 0.827 / 0.510 |
| labelspread | 0.293 / 0.999 | 0.982 / 0.446 | 0.814 / 0.935 | 0.818 / 0.510 |
| lgbm | 0.293 / 1.000 | 0.982 / 0.446 | 0.845 / 0.978 | 0.905 / 0.657 |
| logreg | 0.293 / 0.998 | 0.982 / 0.446 | 0.819 / 0.933 | 0.693 / 0.630 |
| rf | 0.293 / 0.999 | 0.982 / 0.446 | 0.843 / 0.981 | 0.975 / 0.608 |
| selftrain | 0.293 / 0.999 | 0.982 / 0.446 | 0.837 / 0.968 | 0.758 / 0.619 |
| svm | 0.292 / 0.999 | 0.982 / 0.446 | 0.836 / 0.930 | 0.799 / 0.598 |
| xgb | 0.311 / 0.999 | 0.982 / 0.446 | 0.851 / 0.978 | 0.782 / 0.676 |

Variants `stack` and `repro` do not change the ranking: for every algorithm the macro-F1 of the three variants spans at most 0.008 (plant) (`select-plant.md`). Eligibility is evaluated on `main` only.

### 4.1 Winner and generalization

- **Winner: rf** (macro-F1 0.753 vs rules 0.574, p = 5.56e-07; ambiguous share 0.000 vs 0.387; control_real_rate 0.184 vs 0.326). Eligible algorithms: labelprop, labelspread, rf, selftrain, svm. Runner-up selftrain (0.752), p = 0.847 against rf, i.e. indistinguishable; the simpler model (rf) is chosen.
- xgb (0.765) and lgbm (0.760) have a higher macro-F1 but are not eligible: phasiRNA24 precision 0.782 (xgb) and 0.905 (lgbm), below the fold-aligned rules reference printed in `select-plant.md` (0.800 and 0.923) minus 0.01. logreg fails the same criterion (0.693 against 0.706 - 0.01).
- Clustering models (gmm, hdbscan, hier, kmeans) lose phasiRNA recall and phasiRNA24 recall.
- **Complete_rules** (rule class kept where not ambiguous, ML only where the rules abstain; not used for eligibility): rf 0.579, rules 0.574; every algorithm 0.573 to 0.582.

### 4.2 Cross-species folds (main variant; not used for selection)

| algorithm | n | macro_f1 | control_real_rate |
|---|---|---|---|
| rules | 2 | 0.563 ± 0.084 | 0.418 ± 0.395 |
| gmm | 10 | 0.187 ± 0.014 | 0.102 ± 0.061 |
| hdbscan | 10 | 0.197 ± 0.028 | 0.228 ± 0.093 |
| hier | 10 | 0.237 ± 0.122 | 0.102 ± 0.061 |
| kmeans | 10 | 0.182 ± 0.012 | 0.102 ± 0.061 |
| labelprop | 10 | 0.402 ± 0.106 | 0.146 ± 0.099 |
| labelspread | 10 | 0.395 ± 0.114 | 0.149 ± 0.103 |
| lgbm | 10 | 0.368 ± 0.121 | 0.153 ± 0.143 |
| logreg | 10 | 0.413 ± 0.113 | 0.102 ± 0.061 |
| rf | 10 | 0.404 ± 0.120 | 0.199 ± 0.192 |
| selftrain | 10 | 0.360 ± 0.147 | 0.165 ± 0.154 |
| svm | 10 | 0.218 ± 0.007 | 0.035 ± 0.019 |
| xgb | 10 | 0.382 ± 0.116 | 0.149 ± 0.133 |

The rf macro-F1 on the two cross-species folds is 0.404 against 0.563 for the rules; the best ML in this table (logreg) reaches 0.413. In these folds the ML models score 0 for hcsiRNA. This is the main generalization caveat: the chosen model does better than the rules within species but clearly worse across species. n = 10 per ML row (2 folds x 5 seeds) and 2 for the rules, so these estimates are very coarse.

## 5. Arthropod results (variant main)

Fly is the only tuning arthropod, so there are 7 units and no cross-species check; piRNA generalization can only be measured on the final test (fly to Bombyx).

| algorithm | n | macro_f1 | ambiguous_share | control_real_rate | p vs rules | eligible |
|---|---|---|---|---|---|---|
| rules | 7 | 0.614 ± 0.097 | 0.475 ± 0.153 | 0.358 ± 0.169 | - | baseline |
| gmm | 35 | 0.695 ± 0.167 | 0.000 ± 0.000 | 0.404 ± 0.176 | 0.52 | no |
| hdbscan | 35 | 0.739 ± 0.158 | 0.000 ± 0.000 | 0.767 ± 0.202 | 0.327 | no |
| hier | 35 | 0.695 ± 0.167 | 0.000 ± 0.000 | 0.404 ± 0.176 | 0.52 | no |
| kmeans | 35 | 0.695 ± 0.167 | 0.000 ± 0.000 | 0.404 ± 0.176 | 0.52 | no |
| labelprop | 35 | 0.887 ± 0.057 | 0.000 ± 0.000 | 0.590 ± 0.252 | 0.00855 | no |
| labelspread | 35 | 0.886 ± 0.054 | 0.000 ± 0.000 | 0.585 ± 0.249 | 0.00883 | no |
| lgbm | 35 | 0.883 ± 0.051 | 0.000 ± 0.000 | 0.202 ± 0.099 | 0.0112 | no |
| logreg | 35 | 0.868 ± 0.065 | 0.000 ± 0.000 | 0.404 ± 0.176 | 0.0192 | no |
| rf | 35 | 0.881 ± 0.058 | 0.000 ± 0.000 | 0.276 ± 0.136 | 0.0123 | no |
| selftrain | 35 | 0.845 ± 0.081 | 0.000 ± 0.000 | 0.217 ± 0.104 | 0.0273 | no |
| svm | 35 | 0.882 ± 0.058 | 0.000 ± 0.000 | 0.214 ± 0.092 | 0.00989 | no |
| xgb | 35 | 0.879 ± 0.057 | 0.000 ± 0.000 | 0.170 ± 0.082 | 0.0136 | no |

Per-class precision / recall (P / R):

| algorithm | miRNA P / R | piRNA P / R | siRNA P / R |
|---|---|---|---|
| rules | 0.982 / 0.677 | 0.998 / 0.558 | 0.131 / 0.807 |
| gmm | 0.982 / 0.677 | 0.973 / 1.000 | n/a / 0.000 |
| hdbscan | 0.982 / 0.677 | 0.975 / 0.999 | 0.819 / 0.145 |
| hier | 0.982 / 0.677 | 0.973 / 1.000 | n/a / 0.000 |
| kmeans | 0.982 / 0.677 | 0.973 / 1.000 | n/a / 0.000 |
| labelprop | 0.982 / 0.677 | 0.982 / 0.997 | 0.715 / 0.758 |
| labelspread | 0.982 / 0.677 | 0.982 / 0.996 | 0.715 / 0.753 |
| lgbm | 0.982 / 0.677 | 0.983 / 0.998 | 0.820 / 0.776 |
| logreg | 0.982 / 0.677 | 0.984 / 0.995 | 0.560 / 0.805 |
| rf | 0.982 / 0.677 | 0.981 / 0.999 | 0.927 / 0.674 |
| selftrain | 0.982 / 0.677 | 0.981 / 0.996 | 0.589 / 0.615 |
| svm | 0.982 / 0.677 | 0.981 / 0.998 | 0.741 / 0.698 |
| xgb | 0.982 / 0.677 | 0.982 / 0.998 | 0.704 / 0.715 |

Variants `stack` and `repro` do not change the ranking: macro-F1 of the three variants spans at most 0.010 per algorithm (`select-arthropod.md`).

### 5.1 No eligible algorithm

**Winner: none.** Every algorithm fails the piRNA-precision constraint: 0.973 to 0.984 against 0.998 - 0.01 = 0.988 (closest: logreg 0.984). Most also fail siRNA recall (rules 0.807, threshold 0.797): the clustering models 0.000 to 0.145, labelprop 0.758, labelspread 0.753, lgbm 0.776, rf 0.674, selftrain 0.615, svm 0.698, xgb 0.715. Several also fail the control criterion (control_real_rate not below the rules' 0.358): gmm, hier, kmeans 0.404, hdbscan 0.767, labelprop 0.590, labelspread 0.585, logreg 0.404.

What the ML would have gained, for the record: siRNA precision from 0.131 to between 0.560 (logreg) and 0.927 (rf) among the supervised and semi-supervised models; piRNA recall from 0.558 to between 0.995 and 0.999; macro-F1 from 0.614 to between 0.845 and 0.887; ambiguous share from 0.475 to 0. The price is 0.014 to 0.025 of piRNA precision (0.973 to 0.984 against 0.998).

Per the pre-registered rule the arthropod profile stays **rules-only**. The rule is not relaxed after seeing the numbers. No arthropod model is shipped.

## 6. Plant versus arthropod, in one line each

- Plant: the ML wins on macro-F1, abstains never, calls fewer negative controls real, and no class loses more than 0.01 (rf), but it generalizes worse than the rules across species.
- Arthropod: the ML wins on macro-F1 and siRNA/piRNA coverage but fails the pre-registered precision constraint on piRNA.

## 7. The miRNA decision

The first plant run (code ddf0ab4, ML allowed to output miRNA; `select-plant-v1.md`) showed the ML adding miRNA calls on loci that miRDeep2 missed: miRNA recall rose from 0.446 (rules) to 0.754 (rf) and 0.779 (xgb), while precision fell from 0.982 to 0.915 (rf) and 0.856 (xgb). Every supervised and semi-supervised model lost miRNA precision (0.802 to 0.927 against 0.982) and none was eligible; the plant profile had no winner.

The user then decided (2026-11-07) that **miRNA comes only from miRDeep2** and the ML never outputs it. Loci labelled miRNA that miRDeep2 did not carve stay positives for the "is real" score but are not a class target. The code was changed (commit 44c190d) and the cross-validation re-run; all results above are from that re-run. This is a design change made after seeing tuning-data results, and it is stated here openly. The arthropod run on the first code was stopped part-way and archived; the arthropod numbers above are from the re-run.

Future work: an ML layer around miRDeep2 (the v1 numbers show room: recall 0.446 to 0.754 at precision 0.982 to 0.915).

## 8. Known metric artefact: hcsiRNA precision

hcsiRNA precision in the plant table is about 0.29 ± 0.46 for every algorithm and for the rules (0.293 ± 0.468). That is not a measure of hcsiRNA quality. Folds without hcsiRNA truth (the rice folds, since the rice hc-siRNA file was excluded) score precision 0 for any hcsiRNA call, so the mean is a mixture of folds near 1 and folds at 0, hence the large standard deviation. Checked per fold for the rules: in the five Arabidopsis chromosome folds hcsiRNA precision is 0.995–0.998 (e.g. chromosome 1: 12423 TP, 20 FP); in all twelve rice folds n_truth is 0, so every hcsiRNA call is a false positive against another annotated class (43–138 per fold) and precision is 0; the pooled Arabidopsis "other" fold has no hcsiRNA calls (undefined, excluded from the mean). Recall (0.975 for the rules, 0.988 to 1.000 for the ML) is the informative hcsiRNA number, and it hardly differs between algorithms.

## 9. Limitations and next step

- Tuning data only: one fly, one Arabidopsis and one rice experiment; replicates are not independent.
- Labels are annotation-derived and noisy (section 2); recall is locus-level (section 3).
- Rare classes (fly siRNA, Arabidopsis phasiRNA/phasiRNA24) have few loci; per-class standard deviations are large (phasiRNA24 recall for rf is 0.608 ± 0.396).
- Cross-species performance of the plant winner is below the rules (section 4.2).
- The miRNA design change (section 7) was made after seeing tuning results; the selection thresholds themselves were not changed.
- The fixed number of clusters (classes + 2) was not tuned in the folds; the four O(n^2) methods use a 20000-row subsample.

Next step: one-shot evaluation of the winner on the test libraries (Bombyx and maize), only after explicit user approval. **This has not been done.**

## References

- Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324
- Bekker, J., & Davis, J. (2020). Learning from positive and unlabeled data: a survey. *Machine Learning*, 109(4), 719–760. https://doi.org/10.1007/s10994-020-05877-5
- Chen, T., & Guestrin, C. (2016). XGBoost: A scalable tree boosting system. In *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 785–794). https://doi.org/10.1145/2939672.2939785
- Elkan, C., & Noto, K. (2008). Learning classifiers from only positive and unlabeled data. In *Proceedings of the 14th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 213–220). https://doi.org/10.1145/1401890.1401920
- Ke, G., Meng, Q., Finley, T., Wang, T., Chen, W., Ma, W., Ye, Q., & Liu, T.-Y. (2017). LightGBM: A highly efficient gradient boosting decision tree. *Advances in Neural Information Processing Systems*, 30, 3146–3154.
- Nadeau, C., & Bengio, Y. (2003). Inference for the generalization error. *Machine Learning*, 52(3), 239–281. https://doi.org/10.1023/A:1024068626366
- Pedregosa, F., Varoquaux, G., Gramfort, A., et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830.
- Roberts, D. R., Bahn, V., Ciuti, S., et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40(8), 913–929. https://doi.org/10.1111/ecog.02881
