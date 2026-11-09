# ML v2 (biology-gated) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** the ML classifier only assigns a class whose defining biology the locus shows, can abstain through an explicit `other` target, has a supervised "is real" score, and is selected with the v2 criterion (annotation + independent biology).

**Architecture:** the gate is a method of the rules decider (`decider.RuleScorer`), driven by the existing rules YAML, so the pipeline (`mlclassify.py`) and the benchmark CV (`benchmark/ml/evaluate.py`) call the same code. The evidence strengths the gate needs (`e_*`) become `clusters.tsv` columns, so the benchmark can rebuild a record from a table row.

**Tech Stack:** Python 3.13, numpy/pandas/scikit-learn (+xgboost/lightgbm on litrp3), pytest.

**Spec:** `docs/superpowers/specs/2026-11-09-ml-v2-biology-gated-design.md` (extends `docs/superpowers/specs/2026-11-07-ml-classifier-design.md`).

## Global Constraints

- Test data is never used for tuning. Maize is spent; Bombyx stays unseen until the user approves its one-shot run.
- Rules stay the default classifier (`--classifier rules`); rules-mode class calls and metrics must not change (new columns only).
- Fixes go in the pipeline `.py`/YAML, not in benchmark post-processing.
- miRNA comes only from miRDeep2 (R23); the ML never outputs miRNA.
- No push or merge without asking. `tests/` and `extras/` stay gitignored. Commits end with the session's Co-Authored-By / Claude-Session lines.
- Local test command: `/private/tmp/claude-501/-Users-lucasmelo/041e0fc2-0125-4b54-b1b9-cb47570b8137/scratchpad/venv/bin/python -m pytest -q -W error::FutureWarning` (baseline at 44c190d: 402 passed, 5 skipped; xgb/lgbm skip locally).
- Subagents do not dispatch subagents.

## Review Focus

1. A locus whose dicer_call is "N" or whose size fits no class window must end ambiguous in ML mode, never the first class in config order.
2. Rules-mode output must be byte-identical in its existing columns after Task 1 (only `e_*` columns added).
3. A table row rebuilt into a record must gate exactly as the original ClusterRecord did (round-trip test on real rows).
4. `other` must never appear as an output class (TSV, metrics); it maps to ambiguous.
5. Control loci in training must not leak into the control metric: the metric only counts held-out (validation) loci.

---

### Task 1: evidence columns + gate in the decider

**Files:** Modify `reporter.py` (`_COLS`), `records.py` (only if a TSV reader parses columns back), `decider.py`, `rules.yaml`, `rules.plant.yaml`. Test: `benchmark/tests/test_gate.py` (new) or the existing decider test module if one exists.

**Interfaces — Produces:**
- `clusters.tsv` gains `e_1u, e_5a, e_fractop_bid, e_pingpong, e_phase` right after `s_phase`.
- YAML key `decision.ml_gate_min_score: 0.5` in both rules files (comment: spec v2 §2.1, midpoint of the evidence sigmoids; fixed, not tuned).
- `RuleScorer.passes_gate(name: str, record) -> bool`: False if the class is not enabled; if `record.dicer_call == "N"`; if `_outside_window(spec, record)`; if `_cluster_failures(spec, record)` is non-empty; if `not _phase_supported(spec, record)`; if the class has `defining_evidence` and `min(defining_scores(spec, record)) < ml_gate_min_score`. True otherwise.
- `RuleScorer.gated_choice(proba: Mapping[str, float], record) -> str`: walk classes by descending probability; the first that is `"other"` returns `RuleScorer.AMBIGUOUS`; the first real class that passes the gate is returned; `"miRNA"` is skipped; nothing passes → `AMBIGUOUS`.

- [ ] Tests: gate rejects each failure kind (N size, outside window, short span, density p, phase p, low defining score) and accepts a clean piRNA / phasiRNA / hcsiRNA record built with real-looking values; `gated_choice` falls through to the next class, returns ambiguous for `other` on top and when nothing passes; rules-mode `decide` unchanged on the existing decider tests; reporter writes the five new columns.
- [ ] Implement, run the full suite, commit.

### Task 2: gated ML decisions and the `other` target

**Files:** Modify `mlfeatures.py` (or a new root module `mlgate.py` if cleaner), `mlmodels.py`, `mlclassify.py`, `benchmark/ml/dataset.py`, `benchmark/ml/evaluate.py`, `benchmark/ml/cv.py`, `benchmark/ml/train_final.py`. Tests: `benchmark/tests/test_ml_*.py`.

**Interfaces:**
- Consumes Task 1's `passes_gate` / `gated_choice` and the `e_*` columns.
- `GATE_COLS` (root module): the raw columns the gate reads (at least `dicer_call`, `read_span`, `density_pvalue`, `phase_pvalue`, `phase_number`, `frac_top`, `start`, `end`, `e_*`, and whatever `_size` reads). `build_table` keeps them in the table (not as model features).
- `record_from_row(row) -> ClusterRecord` (root module) so the benchmark and the pipeline gate the same object type.
- `mlmodels`: rows labelled with a class in `ML_EXCLUDED_CLASSES` and non-carved control rows (`is_control`) get the class target `"other"`. `classes_` may contain `"other"`.
- `evaluate.decide(model, table, features, t, scorer)`: carved → miRNA; max prob < t → ambiguous; else `scorer.gated_choice(row_proba, record_from_row(row))`. `choose_threshold` uses the same decide. `cv._run_fold` builds the scorer from the profile's rules YAML (`rules.yaml` for arthropod, `rules.plant.yaml` for plant) once per fold.
- `MLDecider.decide_all(records)`: uses `gated_choice` with the records directly; never writes `other`.

- [ ] Tests: round-trip — rows of a real `clusters.tsv` fixture rebuilt with `record_from_row` give the same `passes_gate` as the original records; a model trained with miRNA-labelled and control rows has `"other"` in `classes_` and `decide` never outputs `"other"` or `"miRNA"` on non-carved rows; an unphased 24-nt locus predicted phasiRNA24 falls through to hcsiRNA when it passes hcsiRNA's gate; production `decide_all` on records matches benchmark `decide` on the same rows.
- [ ] Implement, run the full suite, commit.

### Task 3: supervised "is real"

**Files:** `mlmodels.py`, tests.

**Interfaces:**
- `p_real`: binary model per algorithm family (same family choice as R20 / `_pu()` builders), positives = rows with any class label (miRNA and others, non-carved), negatives = control rows; fitted inside the fold. If the training slice has no control rows, fall back to the existing PU score and set `self.real_mode_ = "pu"` (else `"supervised"`).
- `p_real_pu(X)` keeps the old PU score for reporting.

- [ ] Tests: with controls present the supervised mode is used and scores controls lower than labelled rows on a separable synthetic set; without controls the PU fallback is used; carved rows never enter either fit.
- [ ] Implement, run the full suite, commit.

### Task 4: independent biology metric and selection v2

**Files:** Create `benchmark/ml/bioaudit.py`; modify `benchmark/ml/evaluate.py` (`fold_metrics`), `benchmark/ml/cv.py`, `benchmark/ml/select.py`, `benchmark/ml/config.yaml`, `benchmark/PROVENANCE.md`, `.gitignore`. Tests.

**Interfaces:**
- `config.yaml`: `te: {Drosophila melanogaster: dmel/rmsk_te.bed}` (relative to the truth dir, gitignored like controls); `selection: {class_floor: 0.02, alpha: 0.05}`.
- A build step (CLI or documented command, like `build_controls.py`) turns UCSC dm6 `rmsk.txt.gz` (https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz) into `rmsk_te.bed`: repClass in {LTR, LINE, DNA, RC, Other, Unknown}, `chr` prefix stripped. PROVENANCE gets URL, date, sha256, counts.
- `bioaudit.supported(table_val, pred_cls, te_index) -> np.ndarray[bool]`: for unannotated (label null), non-carved, non-ambiguous calls: `repro_frac >= 0.5` and, for class `piRNA`, overlap ≥ 1 bp with a TE interval. Rows that are not unannotated calls are False and excluded from the denominator.
- `fold_metrics` adds `unann_supported_rate` (ML) and `rules_unann_supported_rate` (rules on the same val rows), plus `n_unann_calls` for both. `repro_frac` must be present on every table used by CV (computed by `repro_features`, never added to the main/stack feature lists).
- `select.eligibility` v2: (1) macro-F1 mean > rules and corrected t-test p < alpha; (2) no class loses > class_floor precision or recall (fold-aligned, as now); (3) `unann_supported_rate` ≥ rules', fold-aligned means; (4) control rate < rules'; (5) ambiguous share < rules'. `real_repro_gap` stays reported, no longer gates. Output states the criterion version (v2) and the config values.

- [ ] Tests: `supported` on a hand-built table (repro and TE cases, annotated rows excluded); eligibility fails each condition separately and passes when all hold; a fly row overlapping a TE interval by 1 bp counts.
- [ ] Implement, run the full suite, commit.

### Task 5 (controller, litrp3): rerun and report

- [ ] Deploy HEAD to `~/srnaprofiler-ml`; build `dmel/rmsk_te.bed`.
- [ ] Re-run classify (rules mode) for the 9 tuning libraries so `clusters.tsv` carries `e_*`; reuse qc/mapping/miRDeep outputs by the R25 procedure only after `git diff` confirms Task 1–4 touch no qc/mapping/mirna code; check rules metrics equal run #7 within 0.001.
- [ ] `ml-cv` both profiles, all algorithms, main/stack/repro, 5 seeds; `ml-select` v2.
- [ ] Re-run the fly and maize biology audits with v2 predictions (maize report-only).
- [ ] Write `benchmark/ML-CV-v2.md` (provenance, tables, winner or none, audits, caveats, `## References` with DOIs); commit.
- [ ] Stop and ask the user before any Bombyx (or new-species) one-shot test.
