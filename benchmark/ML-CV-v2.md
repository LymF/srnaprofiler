# ML classifier v2 (biology-gated): cross-validation on tuning data and maize report

Cross-validation of the 12 algorithms of `ML-CV.md` against the rules classifier, rerun with the biology gate (spec `docs/superpowers/specs/2026-11-09-ml-v2-biology-gated-design.md`), selection criterion v3 and the R28 gate (rules' winner checks plus phasing for phase-defined classes). Tuning libraries only. The maize section at the end is a **report-only** run on a test set already spent on the v1 model (`RESULTS-ml.md`); it was not used to choose anything. Bombyx was run once at the end (section 6).

## 1. Provenance

| item | value |
|---|---|
| code | commit 0c4246c (R28), deployed at `~/srnaprofiler-ml` |
| pipeline hash (HEAD) | `e0c10a7c164e73e5da7d99969d089e8747c8d92b08e473723195781888c4fa6a` |
| tuning libraries | fly SRR30544370/71/72; Arabidopsis SRR24236124/25/26; rice DRR076318/19/20 |
| upstream outputs | run #7 qc, mapping and miRDeep outputs; classify rerun in rules mode so `clusters.tsv` carries `e_*` |
| rules check | rules baseline equals v1 exactly: plant macro-F1 0.574 ± 0.089 (18 folds), arthropod 0.614 ± 0.097 (7 folds) |
| cross-validation | chromosome-grouped leave-one-unit-out, 5 seeds, variants main/stack/repro; plant 3620 fold lines, arthropod 1267, 0 errors |
| selection | criterion v3 (class precision floor 0.02, class F1 floor 0.05, alpha 0.05; t-test reported, not gated) |
| raw outputs (server) | `~/mlcv3/plant/folds.jsonl`, `~/mlcv3/arthropod/folds.jsonl`, `select.md` beside each |
| maize model | plant lgbm, main, seed 0, trained on all tuning rows; `~/mlcv3/final/plant_lgbm.joblib`, sha256 `f822c56bab43388eef2b47d613b016dd7d17143902593b2369daf55fe2104490`; classes hcsiRNA, phasiRNA, phasiRNA24, other; threshold 0.0 |

For the maize run the upstream outputs were reused by the R25 procedure: `git diff 99e849e..HEAD` touches only classify and ML code (`classify.py`, `decider.py`, `mlclassify.py`, `mlgate.py`, `mlmodels.py`, `records.py`, `reporter.py`, the two rules yaml files) and adds the `--classifier` and `--ml-model` flags to `srnaprofiler.py`; `work/benchmark_run.json` keeps the old hash in `reused_from_pipeline_sha256`.

The miRNA class is not an ML output (R23): miRNA rows are identical to the rules by construction.

## 2. Selection result

Implementations: scikit-learn, XGBoost and LightGBM; random forest after Breiman (2001). Chromosome-grouped folds follow Roberts et al. (2017). The t-test is the corrected resampled test of Nadeau and Bengio (2003).

| profile | winner | runner-up | macro-F1 rules | macro-F1 winner | p vs rules |
|---|---|---|---|---|---|
| plant | lgbm | xgb (p = 0.944) | 0.574 | 0.587 | 0.401 |
| arthropod | labelprop | labelspread (p = 0.306) | 0.614 | 0.697 | 0.141 |

Eligible (main variant): plant lgbm, xgb, selftrain; arthropod labelprop, labelspread, lgbm, logreg, rf, selftrain, svm, xgb. The four unsupervised methods are never eligible (macro-F1 below the rules; plant phasiRNA24 and fly siRNA F1 collapse to 0). In plants, rf, logreg, svm and the label-propagation methods fail only the independent-biology condition (supported share of unannotated calls below the rules' 0.921).

## 3. Plants (18 chromosome units)

Pooled over units (seed-averaged), rules to lgbm:

| class | precision | recall | F1 |
|---|---|---|---|
| hcsiRNA | 0.979 -> 0.984 | 0.975 -> 0.980 | 0.977 -> 0.982 |
| phasiRNA | 0.935 -> 0.998 | 0.948 -> 0.882 | 0.941 -> 0.936 |
| phasiRNA24 | 0.784 -> 1.000 | 0.169 -> 0.187 | 0.279 -> 0.315 |
| miRNA | 0.982 -> 0.982 | 0.479 -> 0.479 | 0.644 -> 0.644 |

Paired by unit, ML minus rules, with a bootstrap 95% interval over units (4000 resamples):

| metric | lgbm | xgb | rf |
|---|---|---|---|
| macro-F1 | +0.013 [-0.004, +0.037]; better in 9/17 | +0.014 [-0.004, +0.037]; 11/17 | +0.013 [-0.003, +0.034]; 9/17 |
| controls called real (lower is better) | -0.127 [-0.155, -0.108]; worse in 0/18 | -0.127; 0/18 | -0.117; 0/18 |
| supported share of unannotated calls | +0.034 [+0.017, +0.052]; 16/18 | +0.036 [+0.020, +0.053]; 17/18 | -0.002 [-0.017, +0.011]; 9/18 |
| ambiguous share | +0.110 [+0.079, +0.149]; 18/18 | +0.115; 18/18 | +0.095; 18/18 |

Reading: macro-F1 does not move. What the ML buys is precision (fewer controls called real, better supported unannotated calls) and it pays with recall of phasiRNA and with more loci left ambiguous.

## 4. Arthropods (fly, 7 chromosome units)

| class | rules P / R / F1 | labelprop | lgbm |
|---|---|---|---|
| miRNA | 0.979 / 0.682 / 0.804 | same | same |
| piRNA | 0.998 / 0.597 / 0.747 | 0.999 / 0.568 / 0.724 | 0.999 / 0.579 / 0.733 |
| siRNA | 0.147 / 0.800 / 0.248 | 0.799 / 0.596 / 0.683 | 0.955 / 0.420 / 0.583 |

The macro-F1 gain is entirely siRNA, with about 100 truth loci in the whole tuning set. Paired by unit, controls called real: labelprop -0.006 [-0.039, +0.028] (no gain), lgbm -0.096 [-0.144, -0.042]; ambiguous share +0.10 for every algorithm (7/7 units).

Biology audit (out-of-fold predictions, fly, RepeatMasker dm6; scripts in `~/mlcv/audit_fly*.py`, output `~/mlaudit3/`): lgbm calls 8998 unannotated piRNA loci, of which 82% overlap a transposon, 92% have ping-pong z > 2 and the reproducibility fraction is 1.0 (median); the annotated piRNA loci it keeps have 98% TE overlap and 86% ping-pong z > 2. Annotated siRNA or miRNA loci called piRNA: 3 for lgbm, 4 for labelprop. Unannotated siRNA calls: 19 (lgbm), 51 (labelprop).

## 5. Maize (report only, spent test set)

lgbm model of section 1 against the rules on SRR15530410/11/18/19. Rules -> ML.

| library | phasiRNA P | phasiRNA R | phasiRNA called (unannotated) | phasiRNA24 R | phasiRNA24 called (unannotated) | hcsiRNA called |
|---|---|---|---|---|---|---|
| SRR15530410 | 0.975 -> 0.997 | 0.734 -> 0.588 | 1034 (481) -> 426 (46) | 0.125 -> 0.135 | 174 (96) -> 211 (132) | 3173 -> 1564 |
| SRR15530411 | 0.956 -> 0.998 | 0.809 -> 0.753 | 1443 (918) -> 579 (130) | 0.096 -> 0.096 | 280 (232) -> 264 (216) | 3888 -> 2839 |
| SRR15530418 | 0.976 -> 0.998 | 0.739 -> 0.635 | 1160 (622) -> 460 (54) | 0.079 -> 0.058 | 218 (177) -> 155 (123) | 5850 -> 4752 |
| SRR15530419 | 0.977 -> 1.000 | 0.697 -> 0.542 | 905 (388) -> 377 (26) | 0.090 -> 0.098 | 148 (98) -> 158 (104) | 3550 -> 2156 |

miRNA rows are identical. For comparison, the v1 model (`RESULTS-ml.md`) called 2000 to 3000 unannotated phasiRNA24 loci per library; the v2 model calls about 100 to 200, and phasiRNA24 recall no longer jumps (0.25-0.37 in v1, 0.06-0.14 now). Precision of phasiRNA rises to 0.997-1.000 at the price of 6 to 15 points of recall. hcsiRNA has no truth in maize (precision undefined), so its drop in calls is not scored.

## 6. Bombyx mori one-shot (held-out test, user-approved 2026-11-11)

One library (DRR377498), run once per model with `--classifier ml`, upstream outputs reused by the R25 procedure (same diff as section 1). The pre-fixed selection rule names **labelprop** (official result); **lgbm** was run as well, as exploratory, because it is the model that reduces controls called real in the tuning data. Nothing was changed after seeing the numbers.

| item | labelprop (official) | lgbm (exploratory) |
|---|---|---|
| model sha256 | `e8e3112f0f7bf2f55f0c85dc3005ac8e392b3995a4237d8c6139c11e9a3ac496` | `a7bcde28b2d229c91afc811663e66133c7ce8ae81bc89b635dc3a910b891e230` |
| piRNA precision (rules 0.843) | 0.850 | 0.861 |
| piRNA recall (rules 0.949) | 0.691 | 0.449 |
| piRNA loci called (rules 11339; unannotated 10957) | 7537 (7290) | 3926 (3732) |
| ambiguous coverage (rules 0.535) | 0.356 | 0.188 |
| miRNA P / R / called | 1.000 / 0.196 / 99, identical to rules | identical |
| siRNA | no calls (rules: 49 calls, no truth) | no calls |

Reading: on a species the model never saw, the ML removes 33% (labelprop) or 65% (lgbm) of the piRNA calls and gains 0.007 to 0.018 of precision, at a cost of 26 to 50 points of recall. In this library the precision gain does not pay for the recall loss. The loci removed are mostly unannotated (a drop of 3667 and 7225 unannotated calls), but unannotated is not false: no TE file or reproducibility group exists for Bombyx here, so whether the removed calls were real piRNA loci was not audited. The fly-trained piRNA models transfer poorly in recall.

## 7. What the ML adds beyond the gate, and what it uses

Two controls on the tuning folds (scripts `benchmark/scripts/ablate_gate.py` and `benchmark/scripts/feature_importance.py`; raw output on the server in `~/mlcv3/ablation/` and `~/mlcv3/importance/`).

### 7.1 Ablation: rules + gate, no ML

The rules' class is kept when its `gate_<class>` column passes and turned ambiguous otherwise. Differences from the rules, paired by unit (bootstrap 95% interval):

| profile | metric | rules + gate | lgbm |
|---|---|---|---|
| plant | macro-F1 | +0.007 [-0.008, +0.027] | +0.013 [-0.004, +0.037] |
| plant | controls called real | -0.031 [-0.038, -0.025] | -0.127 [-0.155, -0.108] |
| plant | supported share of unannotated calls | -0.015 [-0.020, -0.009] | +0.034 [+0.017, +0.052] |
| plant | ambiguous share | +0.043 | +0.110 |
| arthropod | every metric | +0.000 (the gate changes 0.00% of the loci) | see section 4 |

In plants the gate alone changes 2.4% of the loci. It explains the phasiRNA effect (precision 0.935 -> 0.996, recall 0.948 -> 0.882, the same as lgbm: 0.998 / 0.882). It explains about a quarter of the drop in controls called real and none of the gain in supported calls. phasiRNA24 (precision 0.784 -> 1.000 with the ML) is not touched by the gate. In arthropods the gate does nothing, so the whole effect is the ML.

### 7.2 Feature groups (lgbm, main, seed 0, t = 0; change against the full model)

Permutation of each group on the validation rows (3 repeats), and refit without the group for context and abundance.

| profile | group | macro-F1 | controls called real | supported share |
|---|---|---|---|---|
| plant | size (read-length profile) | -0.204 | +0.275 | +0.002 |
| plant | phasing | -0.017 | +0.051 | -0.001 |
| plant | abundance (permute / refit) | -0.000 / -0.001 | +0.002 / +0.002 | -0.035 / -0.040 |
| plant | ping-pong, nucleotide, shape, mapping, strand | about 0 | about 0 | about 0 |
| arthropod | size | -0.133 | +0.051 | -0.004 |
| arthropod | shape (span, start entropy, density) | -0.110 | +0.091 | -0.001 |
| arthropod | mapping (multimap and unique fraction) | -0.108 | +0.059 | -0.005 |
| arthropod | nucleotide (1U, 10A, 5A, GC) | -0.029 | +0.033 | -0.005 |
| arthropod | ping-pong | -0.023 | +0.058 | -0.005 |
| arthropod | phasing, strand, abundance | -0.014 to +0.001 | small | about 0 |

What this says:

- Plants: the model rests on read size and phasing, the two signals that define the plant classes. That is biologically coherent.
- Plants: the gain in the supported share of unannotated calls disappears when abundance is removed (+0.034 -> about -0.006). It comes from the model preferring abundant loci, which recur across replicates more easily; it is not independent biology. Macro-F1 and controls called real do not depend on abundance, so those gains are not affected.
- Arthropods: size matters, but shape and mapping (multimap fraction) weigh as much, while ping-pong and 1U/10A, the textbook piRNA signature, weigh little. Part of the fly piRNA decision therefore rests on mapping behaviour, which is a plausible proxy (piRNAs come from repeats) and a possible species-specific shortcut. This fits the poor recall on Bombyx (section 6) but does not prove it. Permutation importance understates correlated features (ping-pong with mapping and size).
- The four genomic-context features (`ctx_te`, `ctx_trna`, `ctx_rrna`, `ctx_gene`) are 100% missing in the benchmark tables (checked on both profiles), so the model never used context. Consequence: the TE overlap of the fly supported-rate metric is independent of the model, as intended. It also means the context signals that could remove tRNA/rRNA-derived false positives have never been tried.

### 7.3 Feature-set experiments (lgbm, t = 0; script `benchmark/scripts/variant_experiment.py`)

Four feature sets, leave-one-unit-out on the tuning data, and for arthropods an offline Bombyx test (fit on all fly tuning rows, locus-level metrics, genome annotation as context): `main` (production features, no context), `bio` (read size, strand, 1U/10A/5A, ping-pong, phasing only), `ctx` (`main` plus the fraction of the locus covered by CDS/UTR, tRNA, rRNA; `ctx_te` stays empty), `bio_ctx`.

| arthropod CV | macro-F1 | controls called real | supported share |
|---|---|---|---|
| main | 0.674 | 0.261 | 0.822 |
| bio | 0.657 (-0.017 [-0.039, +0.001]) | 0.287 | 0.820 |
| ctx | 0.673 (-0.001) | 0.096 (see below) | 0.835 (+0.012 [-0.002, +0.036]) |
| bio_ctx | 0.670 (-0.004) | 0.102 | 0.832 |

| Bombyx offline, locus level | piRNA precision | piRNA recall | ambiguous share | macro-F1 |
|---|---|---|---|---|
| rules | 0.843 | 0.695 | 0.465 | 0.702 |
| main | 0.861 | 0.361 | 0.812 | 0.575 |
| bio | 0.873 | 0.268 | 0.815 | 0.526 |
| ctx | 0.836 | 0.594 | 0.501 | 0.668 |
| bio_ctx | 0.838 | 0.624 | 0.485 | 0.678 |

Plants: `bio` leaves macro-F1 unchanged (+0.001) but lowers the supported share (-0.042 [-0.062, -0.024]); `ctx` lowers controls called real (-0.023) and the supported share (-0.031). Neither helps.

What this says:

- **The controls metric is circular for `ctx`.** The negative controls are rRNA features and protein-coding exons of the same annotation (`benchmark/scripts/build_controls.py`), which is what `ctx_gene` and `ctx_rrna` measure. The drop of controls called real with context (arthropod -0.165) is therefore not evidence. Only the metrics against the independent truth count.
- **The mapping/shape shortcut hypothesis is rejected.** Removing mapping, shape and abundance (`bio`) does not repair Bombyx; it makes it worse (recall 0.268). The model is not failing because of those features.
- **The collapse on Bombyx comes from the `other` class.** Without context the model abstains on 81% of the Bombyx loci (rules 47%): the fly-learned boundary between real loci and exon/rRNA controls does not transfer. With genomic context the abstention falls to 49-50% and piRNA recall returns to 0.59-0.62.
- **Even then the ML does not beat the rules on Bombyx:** precision 0.836-0.838 against 0.843, recall 0.594-0.624 against 0.695, macro-F1 0.668-0.678 against 0.702.
- **Reading for production:** the context features are what makes the arthropod model usable on a new species, but they require genome annotation as an input (the tool does not take it today) and, at best, they bring the ML back to the rules' level on Bombyx. One test species and one library.

## 8. Conclusions

1. With the biology gate the ML is a precision tool, not a recall tool. Macro-F1 gains are not significant (plant p = 0.40; fly p = 0.14); the v1 gain (macro-F1 0.753, p = 5.6e-7) did not survive the gate.
2. What holds: fewer controls called real in plants (0/18 units worse), better supported unannotated calls (lgbm, xgb), phasiRNA and phasiRNA24 precision near 1.0 on tuning and on maize, and far fewer unsupported phasiRNA24 calls on maize.
3. What it costs: phasiRNA recall (-0.07 pooled on tuning, -0.06 to -0.16 on maize) and 0.10-0.11 more ambiguous loci.
4. Bombyx (new species, section 6) does not confirm the gain: precision +0.007 (labelprop) or +0.018 (lgbm), piRNA recall down from 0.949 to 0.691 or 0.449. The arthropod ML does not generalise well to an unseen species on this library.
5. Controls (section 7): in plants about three quarters of the drop in controls called real is the ML, not the gate; the gate explains the phasiRNA precision/recall trade. The ML relies on read size and phasing (coherent), but its gain in supported unannotated calls is an abundance effect. In arthropods the gate is inert and the model leans on mapping and shape more than on ping-pong.
6. Feature-set experiments (section 7.3): restricting to textbook features does not help; genomic context repairs the Bombyx collapse (abstention 81% -> 50%) but only reaches the rules' level there, and its effect on controls is circular by construction. A cross-species `other` class is the weak point of the arthropod model.
7. The rules stay the default. `--classifier ml` is an opt-in layer for users who want fewer false positives, best used on plants (phasiRNA) where the gain held on tuning and on maize; for arthropod piRNA it is not recommended.

## 9. Caveats

- The arthropod winner under the pre-fixed rule (labelprop) does not reduce controls called real (-0.006). lgbm does (-0.096) and is also the plant winner. Choosing lgbm for the arthropod profile would be a decision taken after seeing the results; it is not applied here.
- The context features were computed offline from the annotation GFF for this experiment; the production tool does not compute them (`mlfeatures.derive` returns NaN for `ctx_*`).
- Bombyx is a single library and one species; precision there is measured against an incomplete truth.
- Fly siRNA rests on about 100 truth loci in one species; one tuning species for piRNA.
- Arabidopsis phasiRNA24 labels are sRNAanno predictions in a species without DCL5 (noisy truth). Plants have no TE file, so supported means reproducibility only.
- Maize was used for the v1 one-shot; the numbers in section 5 are descriptive and were not used to choose anything.
- The classify v3 logs on the server print a missing-tools message from the dependency check although the steps exit 0; the rules metrics match v1 to three decimals, so the outputs are consistent, but the logs were not audited line by line.
- The suite passes on the server (475 passed, 1 skipped) with the new `bio` variant; it was not run locally.

## References

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324

Chen, T., & Guestrin, C. (2016). XGBoost: A scalable tree boosting system. In *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 785–794). https://doi.org/10.1145/2939672.2939785

Ke, G., Meng, Q., Finley, T., Wang, T., Chen, W., Ma, W., Ye, Q., & Liu, T.-Y. (2017). LightGBM: A highly efficient gradient boosting decision tree. In *Advances in Neural Information Processing Systems 30* (pp. 3146–3154). (No DOI.)

Nadeau, C., & Bengio, Y. (2003). Inference for the generalization error. *Machine Learning*, 52(3), 239–281. https://doi.org/10.1023/A:1024068626366

Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., & Duchesnay, E. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825–2830. (No DOI.)

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schröder, B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F., & Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40(8), 913–929. https://doi.org/10.1111/ecog.02881
