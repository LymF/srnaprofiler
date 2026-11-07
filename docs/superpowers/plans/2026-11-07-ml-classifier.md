# ML Classifier Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an ML classifier that assigns every locus a class (with abstention) and a probability of being real, compare all algorithm families against the rules on tuning data, and ship the winner as `classify --classifier ml` (default stays `rules`).

**Architecture:**
- **Pipeline (repo root, covered by the pipeline hash):**
  - New per-locus features computed in `features.build_record` and written to `clusters.tsv`.
  - A runtime module `mlclassify.py` that loads a frozen model and decides records.
- **Research code (`benchmark/ml/`):**
  - Builds a feature table from finished benchmark runs.
  - Runs the chromosome-grouped CV over an algorithm registry, scores with the existing benchmark metrics, and picks and freezes the winner.

**Tech Stack:** Python 3, numpy, pandas, scikit-learn ≥ 1.3 (includes `HDBSCAN`), xgboost, lightgbm, joblib, scipy. pytest is run in the scratchpad venv. Real runs happen on `litrp3`.

**Spec:** `docs/superpowers/specs/2026-11-07-ml-classifier-design.md`

## Global Constraints

- Test libraries (Bombyx DRR377498; maize SRR15530410/11/18/19) are never used to choose an algorithm, feature, hyperparameter, threshold or class mapping. The selected model runs on test **once**, and that result is reported as is.
- An annotation used as a yardstick for a task is never a feature of the model evaluated on that task.
- `classify --classifier` defaults to `rules`. ML dependencies are optional: rules mode must import and run without scikit-learn, xgboost or lightgbm.
- Main models exclude rule outputs (`cls`, `confidence`, `score`, `s_*`, `e_*`, `review`, `evidence`, `class_scores`) and replicate features. Variants that include them are named `stack` and `repro` and are reported separately.
- Classes per profile:
  - arthropod (`rules.yaml`): miRNA, siRNA, piRNA;
  - plant (`rules.plant.yaml`): miRNA, siRNA, hcsiRNA, phasiRNA, phasiRNA24.
  - A class is learned only if at least one tuning library of that profile has truth for it.
- Labels: a locus takes the class of the truth region chosen by `benchmark.metrics_class.match_loci` with `min_overlap_frac` 0.5; a locus with no match is unlabeled. Unknown is never "false".
- Loci with rule `cls == "miRNA"` (carved from `--mirna-loci`) are excluded from training and keep class miRNA.
- Everything that is fitted is fitted inside the training fold: scaler, imputer, hyperparameters, PU class prior, abstention threshold, and cluster→class mapping.
- At least 5 seeds per algorithm.
- Fixes and features go into the pipeline `.py`/YAML, not into benchmark post-processing. Do not push or merge without asking. `tests/` and `extras/` stay gitignored. Tests live in `benchmark/tests/`.

## Review Focus

1. **A test accession reaches the ML table during development.** Expected: `build_table` raises unless `final_test=True` is passed explicitly (Task 3 test `test_build_table_refuses_test_library`).
2. **A class with no truth in a library.** Expected: its precision and recall are `None` (not 0), and macro-F1 averages only over classes with truth in that library (Task 6 test `test_macro_f1_skips_classes_without_truth`).
3. **Rules mode on a machine without scikit-learn, or ml mode with a model whose feature schema differs.** Expected:
   - rules mode runs unchanged;
   - ml mode exits with a message naming the missing package or the mismatched columns (Task 9 tests `test_rules_mode_without_sklearn`, `test_ml_mode_schema_mismatch`).
4. **A locus with all reads on one strand, or all reads starting at one position.** Expected: histogram fractions are 0 for the empty strand, `start_entropy` is 0.0, and nothing is NaN or raises (Task 1 test `test_shape_features_single_position`).
5. **NaN features (Bombyx replicate features; reserved context columns) reaching models that do not accept NaN.** Expected: every registered model is a pipeline with an imputer before such estimators, so fit and predict succeed on a table with NaN columns (Task 5 test `test_every_model_fits_with_nan_columns`).

---

### Task 1: New per-locus features in the pipeline

**Files:**
- Modify: `records.py` (`ClusterRecord`), `features.py` (`build_record`), `reporter.py` (`_COLS`)
- Test: `benchmark/tests/test_ml_features.py`

**Interfaces:**
- Produces: in `clusters.tsv`, after `read_span`, 36 columns `szp_18 … szp_35` and `szm_18 … szm_35`. Each is the fraction of the locus's 1/NH mass with that length on the + and − strands; together they sum to `frac_in_window`. They are followed by `start_entropy`, `top_start_frac` and `starts_per_kb`.
  - `ClusterRecord` fields have the same names; the `sz*` default is 0.0, and the three shape fields default to `None` (written `NA`).
- Definitions (5′ start = `a.pos` on +, `a.pos + len(a.seq) - 1` on −; weights 1/NH):
  - `start_entropy` = Shannon entropy (base 2) of the weighted 5′-start distribution, divided by log2(number of distinct starts); 0.0 when there is one distinct start.
  - `top_start_frac` = weight of the heaviest start / total weight.
  - `starts_per_kb` = distinct starts / (read_span / 1000).

- [ ] **Step 1: Write the failing tests** (build alignments with the helper in `benchmark/tests/bamutil.py` or the same light alignment objects that `test_locus_floor.py` uses):
  - `test_size_histogram_by_strand`:
    - Input: three 21-nt reads on +, each NH 1, and one 24-nt read on − with NH 2.
    - Assert `szp_21 == round(3/3.5, 4)`, `szm_24 == round(0.5/3.5, 4)`, and every other `sz*` is 0.
    - Assert the sum of `sz*` equals `frac_in_window` within 1e-3.
  - `test_shape_features_single_position`:
    - Input: five reads, all on + and all starting at position 100.
    - Assert `start_entropy == 0.0`, `top_start_frac == 1.0`, and every `szm_*` == 0.
  - `test_shape_features_uniform`:
    - Input: four reads, NH 1, starting at four distinct positions.
    - Assert `start_entropy == 1.0` and `top_start_frac == 0.25`.
  - `test_tsv_has_new_columns`:
    - Call `write_tsv` on one record.
    - Assert the header contains `szp_18`, `szm_35` and `start_entropy`, in the order defined above, right after `read_span`.
- [ ] **Step 2: Run** `pytest benchmark/tests/test_ml_features.py -v`. Expected: FAIL (missing attributes).
- [ ] **Step 3: Implement** in `build_record`, inside the existing per-read loop: accumulate per-strand length weights and 5′-start weights, then round to 4 decimals like the existing fields.
- [ ] **Step 4: Run** the new test file and the full suite `pytest benchmark/tests -q`. Expected: all pass. Update any existing test that pins the TSV header.
- [ ] **Step 5: Commit** `feat(features): per-strand size histogram and read-start shape per locus`.

### Task 2: Re-run classify on all 14 libraries with the new features (litrp3)

**Files:**
- Create: `benchmark/scripts/rerun_classify.sh` (driver)

**Interfaces:**
- Consumes: Task 1 commit.
- Produces: `~/srnaprofiler/data/benchmark/results/<ACC>/work/<ACC>.clusters.tsv` with the new columns for all 14 libraries, metrics rebuilt, and report `~/bench_report7.md`.

- [ ] **Step 1:** Write the driver. For each accession in `benchmark/manifest.tsv`:
  - delete `<ACC>.clusters.tsv`, so the step is not skipped as fresh (`_fresh` looks only at mtimes);
  - then run `python -m benchmark run --library <ACC> --results-dir <results> --threads 4`;
  - run up to 14 in parallel, appending `<ACC> exit=<code>` to a status file and writing `ALL_DONE` at the end.
  - qc, mapping and mirna are reused because their outputs are fresh.
- [ ] **Step 2:** Ship the commit to `litrp3` as before (`git bundle` → `scp` → `git clone -b feat/benchmark ~/srnaprofiler-run7`). Run the driver with `nohup`. Expected: 14 × `exit=0`.
- [ ] **Step 3:** Run `python -m benchmark report > ~/bench_report7.md`. Check that every classification metric equals run #6 within 0.001: the new features do not touch the decision, so any change is a bug to stop on.
- [ ] **Step 4:** Commit the driver: `chore(benchmark): classify-only rerun driver`. Add a ledger line with the run #7 hash.

### Task 3: Feature table builder

**Files:**
- Create: `benchmark/ml/__init__.py`, `benchmark/ml/dataset.py`, `benchmark/ml/config.yaml`
- Test: `benchmark/tests/test_ml_dataset.py`

**Interfaces:**
- Consumes:
  - `benchmark.manifest.load_manifest`, `load_config`;
  - `benchmark.orchestrate._load_truth` (promote it to the public name `load_truth`, keeping the old name as an alias);
  - `benchmark.metrics_class.read_loci`, `match_loci`.
- Produces:
  - `FEATURES: list[str]`: the main-model feature columns. These are the numeric existing columns `reads, reads_rpm, distinct_seqs, size_purity, frac_in_window, frac_top, pct_1u, pct_10a, pct_5a, gc, pingpong_z, pingpong_pairs10, pingpong_frac, overlap_dup_frac, phase_score, phase_pvalue, phase_number, phase_ratio, multimap_frac, unique_frac, density_pvalue, read_span`; plus `dicer_len`, which is `dicer_call` as a float, NaN for `N`; plus the 39 new Task 1 columns; plus the reserved context columns `ctx_te, ctx_trna, ctx_rrna, ctx_gene`, always NaN for now.
  - `RULE_COLS: list[str]`: the columns of the `stack` variant, which are `score, confidence` and the `s_*` columns.
  - `REPRO_COLS = ["repro_frac", "repro_rpm_ratio"]`.
  - `build_table(accessions: list[str], results_dir: Path, manifest: Path, final_test: bool = False) -> pandas.DataFrame`. One row per locus, with these columns: `acc, organism, tissue, replicate, profile` (`arthropod`/`plant`, from `taxon_group`), `chrom, start, end, rule_cls` (the rules' `cls`), `label` (truth class or `None`), `is_control` (bool, from Task 4; False until then), then the `FEATURES`, `RULE_COLS` and `REPRO_COLS` (NaN until Task 4).
  - It raises `ValueError` if any accession has `role == "test"` and `final_test` is False.
  - Rows with `rule_cls == "miRNA"` are kept but flagged with `carved = True`.

- [ ] **Step 1: Write the failing tests** against a tiny synthetic results dir (two loci, one truth region) built in `tmp_path`:
  - `test_build_table_labels_by_match`: a locus overlapping a 1 kb `piRNA` region by 60% gets `label == "piRNA"`; a far locus gets `label is None`.
  - `test_build_table_refuses_test_library`: `pytest.raises(ValueError)` for a test accession with the default `final_test`.
  - `test_build_table_features_exclude_rule_outputs`: `set(FEATURES) & {"cls", "confidence", "score", "evidence"} == set()`, and nothing starting with `s_` is in `FEATURES`.
  - `test_carved_flag`: a locus with `cls == "miRNA"` has `carved` True.
- [ ] **Step 2: Run** the tests. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** `build_table`. The truth comes from `load_truth(lib)`, the same regions the benchmark scores against, seeds included for tuning. `benchmark/ml/config.yaml` holds the `min_overlap_frac` value 0.5 and the fold rule of Task 5.
- [ ] **Step 4: Run** `pytest benchmark/tests -q`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(benchmark/ml): locus feature table from benchmark runs`.

### Task 4: Negative controls and replicate features

**Files:**
- Create: `benchmark/scripts/build_controls.py`, `benchmark/ml/controls.py`, `benchmark/ml/replicates.py`
- Modify: `benchmark/ml/dataset.py` (fill `is_control`, `repro_*`), `benchmark/PROVENANCE.md` (annotation sources)
- Test: `benchmark/tests/test_ml_controls.py`

**Interfaces:**
- Produces:
  - `build_controls.py <gff3> <truth_bed_or_gff3 ...> --aliases <tsv> -o <controls.bed>`. It writes the union of `rRNA` gene features and `exon` features of protein-coding transcripts, minus every interval that overlaps any listed truth region of any class. The output is BED, with column 4 = `rRNA` or `exon`.
  - Sources (pinned):
    - Ensembl release 116 GFF3 for `drosophila_melanogaster` (BDGP6.54), `arabidopsis_thaliana` (TAIR10), `oryza_sativa` (IRGSP-1.0) and `zea_mays` (Zm-B73-REFERENCE-NAM-5.0), on the same FTP tree as the genome URLs in the manifest;
    - for Bombyx, the NCBI RefSeq annotation of GCF_030269925.1. Its contigs are renamed to the GCA names of the benchmark genome through the assembly report's `GenBank-Accn`↔`RefSeq-Accn` columns, written as an aliases TSV.
  - Output paths are `benchmark/truth/<dir>/controls.bed` (`dmel`, `atha`, `osat`, `zmay`, `bmor`). They are generated on `litrp3` and not committed (size). Their sha256 is recorded in `PROVENANCE.md`.
  - `controls.is_control(table, controls_bed) -> pandas.Series[bool]`: a locus is a control if it overlaps a control interval by ≥ 0.5 of the locus or of the interval (the same rule as `match_loci`).
  - `replicates.repro_features(table) -> DataFrame[repro_frac, repro_rpm_ratio]`. The replicate group is (`organism`, `tissue`). `repro_frac` is the fraction of the other replicates in the group that have a locus overlapping this one under the same 0.5 rule. `repro_rpm_ratio` is the median RPM of the matched loci divided by this locus's RPM (NaN if there is no match). Both are NaN when the group has one library.

- [ ] **Step 1: Write the failing tests:**
  - `test_controls_drop_truth_overlap`: an exon overlapping a miRNA truth region is absent from the output; an exon that overlaps no truth is present.
  - `test_is_control_rule`: a locus 50% inside an exon is a control, and one 40% inside is not.
  - `test_repro_three_replicates`: a locus present in 1 of 2 other replicates has `repro_frac == 0.5`.
  - `test_repro_single_replicate_nan`: Bombyx-like single library gives NaN for both columns.
- [ ] **Step 2: Run** the tests. Expected: FAIL.
- [ ] **Step 3: Implement**, reusing `benchmark.metrics_class.IntervalIndex` for overlaps.
- [ ] **Step 4: Run** `pytest benchmark/tests -q`. Expected: all pass.
- [ ] **Step 5:** On `litrp3`, download the five annotations, build the five `controls.bed` files, record the sha256 values in `PROVENANCE.md` together with the URLs and dates, and commit: `feat(benchmark/ml): negative-control regions and replicate reproducibility`.

### Task 5: Algorithm registry and CV folds

**Files:**
- Create: `benchmark/ml/models.py`, `benchmark/ml/folds.py`
- Test: `benchmark/tests/test_ml_models.py`

**Interfaces:**
- Produces:
  - `folds.chromosome_units(table: DataFrame, min_share: float = 0.05) -> pandas.Series[str]`. It gives each row a fold unit: a chromosome holding at least `min_share` of its species' loci is its own unit, and all other contigs are pooled into `"other"`. The value 0.05 lives in `benchmark/ml/config.yaml`.
  - `folds.split(table) -> Iterator[(train_idx, val_idx, unit)]`: leave-one-unit-out within each profile, over tuning rows only.
  - `models.REGISTRY: dict[str, Callable[[int], Model]]` (the argument is the seed), with these entries:
    - unsupervised: `kmeans`, `gmm`, `hdbscan` (`sklearn.cluster.HDBSCAN`), `hier` (`AgglomerativeClustering`);
    - semi-supervised: `labelprop`, `labelspread`, `selftrain` (`SelfTrainingClassifier` over a gradient-boosting base);
    - supervised: `logreg`, `rf`, `xgb`, `lgbm`, `svm`.
    - Plus `rules`, a passthrough that returns `rule_cls`. It is the baseline and has no fitting.
  - `Model` protocol:
    - `fit(X: DataFrame, y: Series[str | None]) -> Model`; `y` holds `None` for unlabeled rows, and supervised models drop those rows.
    - `predict_proba(X) -> DataFrame` (columns = learned classes; rows sum to 1).
    - `p_real(X) -> ndarray`.
  - Every model is an sklearn `Pipeline` with `SimpleImputer(strategy="median")` and `StandardScaler` before estimators that need them.
  - Real score:
    - supervised and semi-supervised models use a PU classifier, positive = any label, with the Elkan–Noto correction (class prior c estimated on a held-out 10% of the training-fold positives);
    - for unsupervised models it is 1 − the share of the HDBSCAN noise label among the k nearest training loci (k = 15), and for the other clustering models it is the PU classifier.
  - Cluster→class mapping (unsupervised): the majority label among the labeled training rows of each cluster; a cluster with no labeled rows maps to no class (its loci are predicted ambiguous).
  - Clustering models use the number of learned classes + 2 as k / n_components. **[builder may tune this inside the fold]**

- [ ] **Step 1: Write the failing tests** on a synthetic table: three Gaussian blobs, one per class, 10% labeled, plus two all-NaN context columns and one NaN `repro_frac` column.
  - `test_every_model_fits_with_nan_columns`: for each name in `REGISTRY`, `fit` then `predict_proba` returns rows that sum to 1 ± 1e-6, and `p_real` lies in [0, 1].
  - `test_supervised_beats_chance_on_blobs`: `rf` accuracy > 0.9 on the labeled rows held out.
  - `test_chromosome_units_pool_small_contigs`: three contigs holding 60%, 38% and 2% of loci give units `{c1, c2, "other"}`.
  - `test_folds_never_mix_units`: no unit appears in both train and val of a split.
  - `test_mapping_fitted_on_train_only`: labels present only in the val rows never appear as a mapped class.
- [ ] **Step 2: Run** the tests. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `pytest benchmark/tests -q`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(benchmark/ml): algorithm registry and chromosome-grouped folds`.

### Task 6: Abstention and evaluation

**Files:**
- Create: `benchmark/ml/evaluate.py`
- Test: `benchmark/tests/test_ml_evaluate.py`

**Interfaces:**
- Consumes: `benchmark.metrics_class.classification_metrics` and `_ece`; Task 5 `Model`.
- Produces:
  - `choose_threshold(model, X_train, y_train, rule_cls_train, seed) -> float`.
    - It runs an inner 3-fold split by chromosome unit inside the training rows.
    - It returns the smallest t in {0.00, 0.05, …, 0.95} such that, on annotated loci, the macro precision of predictions with max probability ≥ t is ≥ the rules' macro precision on the same loci.
    - If no t qualifies, it returns 0.95.
  - `decide(model, table: DataFrame, features: list[str], t: float) -> Series[str]`: the argmax class over `table[features]`, or `"ambiguous"` when max proba < t; rows with `table.carved` keep `"miRNA"`.
  - `fold_metrics(table_val, pred_cls, p_real, rule_cls) -> dict` with these keys:
    - `per_class` (from `classification_metrics`, on the val rows turned into locus dicts with `cls = pred`);
    - `macro_f1` (mean F1 over classes whose truth count is > 0 in these rows);
    - `ambiguous_share`;
    - `real_repro_gap` (mean `repro_frac` of loci with p_real ≥ 0.5 minus that of loci < 0.5; None if `repro_frac` is all NaN);
    - `control_real_rate` (share of `is_control` loci with p_real ≥ 0.5) and `rules_control_rate` (share of `is_control` loci with a non-ambiguous rule class);
    - `ece`.
  - `complete_rules(rule_cls, pred_cls) -> Series[str]`: the "complete the rules" variant. It keeps the rule class where it is not ambiguous and takes the ML class elsewhere.

- [ ] **Step 1: Write the failing tests:**
  - `test_macro_f1_skips_classes_without_truth`: a val set with piRNA truth only gives `per_class["siRNA"]["precision"] is None`, and `macro_f1` equals the piRNA F1.
  - `test_decide_abstains_below_threshold`: proba rows (0.5, 0.5) with t = 0.6 give `"ambiguous"`.
  - `test_carved_stays_mirna`.
  - `test_complete_rules_keeps_rule_calls`: rule `["piRNA", "ambiguous"]` and ML `["siRNA", "piRNA"]` give `["piRNA", "piRNA"]`.
  - `test_control_rates`: two control loci, one with p_real 0.9 and rule class piRNA, the other with p_real 0.1 and rule ambiguous, give `control_real_rate == 0.5` and `rules_control_rate == 0.5`.
- [ ] **Step 2: Run** the tests. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `pytest benchmark/tests -q`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(benchmark/ml): abstention threshold and fold metrics`.

### Task 7: CV runner, statistics and selection

**Files:**
- Create: `benchmark/ml/cv.py`, `benchmark/ml/stats.py`, `benchmark/ml/select.py`
- Modify: `benchmark/__main__.py` (subcommands `ml-cv`, `ml-select`)
- Test: `benchmark/tests/test_ml_cv.py`

**Interfaces:**
- Produces:
  - `python -m benchmark ml-cv --profile {arthropod,plant} --algorithms <names|all> --variants main,stack,repro --seeds 5 --results-dir <dir> --out <dir>`.
    - It builds the table for the tuning accessions of the profile (`final_test` False).
    - For each (variant, algorithm, seed, fold) it fits, chooses t, and writes one JSON line to `<out>/folds.jsonl`: `{profile, variant, algorithm, seed, unit, t, metrics, complete_rules_metrics}`.
    - Plant profiles also run `cross_species` folds (train Arabidopsis → validate rice, and the reverse) with `unit = "xspecies:<org>"`.
    - Runs are parallel across folds (joblib, `n_jobs` flag, default 8).
  - `stats.corrected_resampled_ttest(a: list[float], b: list[float], n_train: int, n_test: int) -> (t, p)`, as in Nadeau & Bengio (2003): variance factor (1/k + n_test/n_train).
  - `python -m benchmark ml-select --folds <out>/folds.jsonl -o <out>/selection.md` applies spec §7.3 in this order:
    1. Eligible algorithms are those of variant `main` where, in mean over folds and seeds:
       - no class has precision or recall more than 0.01 below `rules`;
       - `ambiguous_share` < that of rules;
       - `real_repro_gap` > 0;
       - `control_real_rate` < `rules_control_rate`.
    2. The winner has the highest mean `macro_f1`.
    3. If the winner and the runner-up do not differ (corrected t-test p ≥ 0.05), take the simpler one by the fixed order `logreg < kmeans < gmm < hier < hdbscan < labelprop < labelspread < svm < rf < selftrain < lgbm < xgb`.
    4. Write a markdown table of every algorithm × variant: mean ± sd of each metric, the p-value against rules, eligibility, and the winner. The `complete_rules` variant goes in its own column group.
- [ ] **Step 1: Write the failing tests:**
  - `test_corrected_ttest_known_value`: k = 10, a − b differences all 0.01 with sd 0.005, n_train/n_test = 9 gives the hand-computed t (write the expected value in the test, computed from the formula).
  - `test_select_rejects_class_loss`: a fake jsonl where `xgb` beats rules on macro-F1 but loses 0.02 precision on one class means `xgb` is not eligible.
  - `test_select_tie_prefers_simpler`: `rf` and `xgb` have equal means and p = 0.9, so `rf` wins.
  - `test_cv_end_to_end_tiny`: the synthetic table from Task 5 with two algorithms, one seed and two units gives a `folds.jsonl` with 2 × 1 × 2 lines for the variant `main`.
- [ ] **Step 2: Run** the tests. Expected: FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `pytest benchmark/tests -q`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(benchmark/ml): CV runner, corrected t-test and selection rule`.

### Task 8: Run the comparison on litrp3 (tuning only)

**Files:**
- Create: `benchmark/ML-CV.md` (results, written from `selection.md`)

- [ ] **Step 1:** In the `srnaprofiler` env on `litrp3`, install `scikit-learn>=1.3 xgboost lightgbm joblib` with pip, and record the versions in `benchmark/ML-CV.md`.
- [ ] **Step 2:** Run `ml-cv` for `arthropod` and `plant`, with all algorithms, all three variants and 5 seeds, under `nohup` with a status file. Expected: `folds.jsonl` for each profile, with no failed fold. A failed fold is an error to fix, not to drop.
- [ ] **Step 3:** Run `ml-select` for each profile.
- [ ] **Step 4:** Write `benchmark/ML-CV.md`:
  - provenance (commit, pipeline hash, run #7 report);
  - the selection table;
  - the winner per profile, or "none eligible" if no algorithm passes. In that case stop and report to the user; do not relax the criterion;
  - the `complete_rules` comparison;
  - caveats: few piRNA folds, the "real" proxies are proxies;
  - `## Referências`: Elkan & Noto 2008, Bekker & Davis 2020, Nadeau & Bengio 2003, and the scikit-learn, XGBoost and LightGBM papers, all with DOIs.
- [ ] **Step 5: Commit** `docs(benchmark): ML cross-validation results (tuning only)`.

### Task 9: Freeze the winner and integrate `--classifier ml`

**Files:**
- Create: `benchmark/ml/train_final.py`, `mlclassify.py` (repo root), `models/arthropod.joblib`, `models/plant.joblib`, `models/*.schema.json`
- Modify: `classify.py`, `srnaprofiler.py` (classify args), `reporter.py` (`p_real` column), `rules.yaml`, `rules.plant.yaml` (`ml.model` key), `README.md`
- Test: `benchmark/tests/test_mlclassify.py`

**Interfaces:**
- Produces:
  - `python -m benchmark ml-train --profile <p> --algorithm <winner> --variant main --seed 0 -o models/<p>.joblib`. It trains on all tuning rows of the profile (`final_test` False) and chooses t by the Task 6 rule on all tuning rows. It writes the joblib and `<p>.schema.json` = `{"features": FEATURES, "classes": [...], "threshold": t, "algorithm": ..., "train_commit": ..., "pipeline_hash": ...}`.
  - `mlclassify.load(path) -> MLDecider` imports sklearn lazily. If sklearn is missing it raises `RuntimeError("--classifier ml needs scikit-learn (pip install scikit-learn xgboost lightgbm)")`.
  - `MLDecider.decide_all(records) -> None`:
    - builds the feature matrix from the record attributes in schema order and raises `ValueError` listing missing or extra columns on mismatch;
    - sets `cls`, `confidence` (max proba), `review` (True when ambiguous), `p_real`, and `evidence = "ml:<algorithm>"`;
    - skips records with `cls == "miRNA"`.
  - `classify --classifier {rules,ml}`, default `rules`, and `--ml-model PATH`, default from YAML `ml: {model: models/arthropod.joblib}` (plant: `models/plant.joblib`), resolved relative to the repo root.
  - `ClusterRecord.p_real: Optional[float] = None`; `_COLS` gains `p_real` after `score`, written `NA` in rules mode.
- [ ] **Step 1: Write the failing tests:**
  - `test_rules_mode_without_sklearn`: monkeypatch `sys.modules["sklearn"] = None`; `run_classify` in rules mode on the small BAM fixture writes `clusters.tsv`, with `p_real` all `NA`.
  - `test_ml_mode_schema_mismatch`: a schema listing a feature the record lacks gives `ValueError` naming it.
  - `test_ml_mode_sets_fields`: a tiny model trained on the fixture's own table, after `decide_all`, leaves every non-miRNA record with `cls` in the schema classes or `"ambiguous"` and with `0 <= p_real <= 1`.
  - `test_default_is_rules`: parse `classify` args without `--classifier`, which gives `"rules"`.
- [ ] **Step 2: Run** the tests. Expected: FAIL.
- [ ] **Step 3: Implement.** Train the two shipped models with `ml-train` on `litrp3`, and copy the artifacts and schema into `models/`. Document `--classifier`, `--ml-model`, `p_real` and the optional dependencies in `README.md`.
- [ ] **Step 4: Run** `pytest benchmark/tests -q`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(classify): optional ML classifier (--classifier ml), rules stay default`.

### Task 10: One-shot test evaluation (ask the user before launching)

**Files:**
- Modify: `benchmark/run.py` (pass `--classifier`/`--ml-model` through), `benchmark/__main__.py` (`run --classifier`)
- Create: `benchmark/RESULTS-ml.md`

- [ ] **Step 1:** Thread `--classifier` through `benchmark run`, with default `rules`. Add a test asserting that the classify argv contains `--classifier ml` when it is requested. Commit.
- [ ] **Step 2: Stop and ask the user** before running. This is the single allowed look at test data.
- [ ] **Step 3:** On approval, run classify in ml mode on the 14 libraries, with tuning included for the report: results go to a separate results dir, `~/srnaprofiler/data/benchmark/results_ml`, reusing the BAMs. Run `benchmark report`. Run `ml-cv`'s `fold_metrics` once on the test tables with `final_test=True`, for the "real" proxies.
- [ ] **Step 4:** Write `benchmark/RESULTS-ml.md`:
  - rules (run #7) against ML, for tuning (in-sample, labeled as such) and test;
  - the `complete_rules` variant;
  - the "real" proxies on test, which for Bombyx means the controls only;
  - caveats;
  - `## Referências` with DOIs.
  - No change to any model or threshold after seeing these numbers.
- [ ] **Step 5: Commit** `docs(benchmark): ML vs rules on tuning and test`.
