# ML classifier for srnaprofiler: design

Status: approved in conversation 2026-11-07 (user: "está certo, pode escrever"). Items marked **[builder's call / Claude's call]** were not decided by the user.

## 1. Intent

The paper makes two claims: srnaprofiler classifies small-RNA loci better than other tools, and machine learning improves that classification. Today a rule-based decider (`decider.py`, `rules.yaml` / `rules.plant.yaml`) assigns each locus a class; the run #6 benchmark (`benchmark/RESULTS-run6.md`) shows it is precise on annotated loci, but it gives up on many loci (ambiguous: fly 51%, rice 56%, Bombyx 46%, maize 45%) and it has no notion of whether a locus is a real small-RNA cluster or noise/degradation.

The ML stage must learn, from per-locus features, two things:

1. **Class** of each locus (the classes of the active rules profile: miRNA, siRNA, piRNA for arthropods; miRNA, siRNA, hcsiRNA, phasiRNA, phasiRNA24 for plants).
2. **Is it real**: a probability that the locus is a genuine small-RNA-producing cluster rather than noise or degradation.

**The ML decides every locus** (user decision 2026-11-07), not only those the rules leave ambiguous: it is a full classifier that may change a class the rules assign. It must be able to **abstain** (report ambiguous) when not confident, so that a lower ambiguous share is earned, not forced; the abstention threshold is chosen on tuning data. The "complete the rules" variant (rules keep their calls, ML decides only the rule-ambiguous loci) is reported alongside as a comparison. **[variant reporting: Claude's call]**

The final aim is to improve srnaprofiler itself. The ML enters as an **optional mode** (`--classifier ml`); the rules stay the default until the user explicitly decides otherwise.

## 2. Why features, not only annotations

Annotated loci are a small share of calls (fly piRNA: 28% of calls fall on any annotation; Bombyx piRNA: 3%). No public set of "false" loci exists. So the work is built on the features of **all** loci; annotations are used as labels where they exist (semi-supervised and supervised families) and as the yardstick for evaluation. Unknown is never treated as false.

## 3. Scope

- **Organisms now:** the five benchmark species. Tuning (development): fly SRR30544370/71/72, Arabidopsis SRR24236124/25/26, rice DRR076318/19/20. Test: Bombyx DRR377498, maize SRR15530410/11/18/19. Wild-type, public data only.
- **Later:** more organisms. The design must let a new species be added by adding manifest rows and truth files, with no code change.
- **Non-goal now:** genomic-context features (genes, transposons, rRNA as features). They come from a separate annotation module, built after this one. This stage must leave room for them (§5.3).

## 4. Hard constraints (anti-goals)

1. **Test data (Bombyx, maize) is never used** to choose an algorithm, feature, hyperparameter, threshold or class mapping. The selected model is run on test **once**, and that number is reported as is.
2. **No label leakage:** an annotation used as a yardstick for a task is never a feature of the model evaluated on that task (e.g. Rfam miRNA families vs miRBase; rRNA annotation used as a negative control vs any rRNA feature).
3. **The default classifier stays the rules** until the user says otherwise.
4. Fixes and features go into the pipeline (`.py`/YAML), not into benchmark post-processing.
5. `tests/` and `extras/` stay gitignored; no push or merge without asking.

## 5. Features

### 5.1 Already computed (per locus, `clusters.tsv`)

reads, reads_rpm, distinct_seqs, major_rna_reads, multimap_frac, unique_frac, dicer_call, dicer_call_2, size_purity, frac_in_window, read_span, gc, frac_top, strand_call, pct_1u, pct_10a, pct_5a, pingpong_z, pingpong_pairs10, pingpong_frac, overlap_dup_frac, phase_score, phase_pvalue, phase_number, phase_ratio, density_pvalue.

The rule-derived outputs (s_*, e_*, score, confidence, class_scores, cls, review) are **not** features of the main models, so that "features alone" is a clean comparison against the rules. A separate, clearly labeled variant may stack ML on top of the rule scores. **[Claude's call]**

### 5.2 New, from the existing BAMs (first round)

1. **Full size profile:** fraction of 1/NH read mass per length 18–35 nt, separately for each strand (36 values), plus the same for 5′ nucleotide by length where cheap. **[exact encoding: builder's call]**
2. **Coverage shape inside the locus:** e.g. entropy of read 5′ start positions, fraction of mass in the top position, number of distinct start positions per kb. **[exact set: builder's call]**
3. **Replicate reproducibility:** for each locus, in how many other replicates of the same condition an overlapping locus is called, and the correlation of its expression. Replicate groups: fly 3, Arabidopsis 3, rice 3, maize 2 + 2 (1.5 mm and 2.0 mm anthers), Bombyx 1 (single replicate: the feature is missing). Because it is also the yardstick of the "is real" task (§7.2), it is used **either** as a feature **or** as an evaluation proxy in a given model, never both; the main models exclude it, and a variant with it is reported separately. **[Claude's call]**

### 5.3 Reserved for the annotation module (later)

Overlap with transposons/repeats, tRNA, rRNA/snoRNA (non-miRNA Rfam families), and genes when the user supplies a GFF. Feature columns are reserved now; models must run when they are missing (models that accept missing values, or separate with/without-context variants).

## 6. Labels and data split

- **Unit:** one locus of a run #6-style output (`clusters.tsv`, after the locus floor), per library. miRNA loci carved from miRDeep2/miRDeep-P2 are excluded from training and keep the miRNA class as today. **[Claude's call]**
- **Class labels:** a locus takes the class of the truth region it matches under the benchmark rule (`benchmark/metrics_class.match_loci`, reciprocal overlap ≥ `min_overlap_frac` = 0.5); loci with no match are unlabeled.
- **One model per rules profile** (arthropod, plant), predicting that profile's classes. **[Claude's call]**
- **Cross-validation on tuning data:** grouped by **chromosome** within a species, holding out whole chromosomes across all replicates, so the same genomic locus never sits in train and validation through another replicate. Plants also get a cross-species check (train Arabidopsis → validate rice and vice versa). piRNA has only one tuning species (fly), so its cross-species generalization is measured only on the final test (fly → Bombyx). **[Claude's call]**
- **Seeds:** at least 5 per algorithm.
- **Folds:** a chromosome-level unit (e.g. chromosome arm in fly: 2L, 2R, 3L, 3R, X, with small chromosomes and scaffolds pooled into one fold), leave-one-unit-out. For piRNA this gives few folds and low test power; that is accepted and reported, not worked around. **[Claude's call]**
- **Everything fitted is fitted inside the fold**: scalers, hyperparameters, PU class prior, and the cluster→class mapping of unsupervised models. **[Claude's call]**
- **Classes without labels:** a class is learned only if at least one tuning species has truth for it. Annotated truth is used as label as is, even where the source is noisy (sRNAanno), as in the benchmark. Per-class metrics and the 0.01 constraint are computed only where that class has truth in the library; macro-F1 averages over those classes. **[Claude's call]**

## 7. Evaluation

The baseline in every comparison is the rules classifier at the run #6 commit (a5c7615, pipeline hash b89dedd6), same libraries, same truth, same benchmark code.

### 7.1 Class

Per class, per library: precision and recall on annotated loci, computed with the existing benchmark functions; plus the share of ambiguous loci.

### 7.2 Is real

No ground truth exists, so two proxies, neither used for training the evaluated model:

1. **Reproducibility:** loci the model calls real should be called in more of the other replicates than loci it rejects.
2. **Negative controls:** loci over rRNA and over protein-coding exons are treated as degradation, after removing any control region that overlaps sRNA truth of any class (no mRNA-seq exists, so "highly expressed" is not used; the control is all protein-coding exons). A better model calls fewer of them real than the rules call any non-ambiguous class (for the rules, "real" = a non-ambiguous call). The coordinates come from public annotation (Ensembl GFF3 for fly, Arabidopsis, rice, maize; NCBI for Bombyx), downloaded for the benchmark only, like the other truth files. Control loci are excluded from training. **[Claude's call]**

Reproducibility (proxy 1) uses the benchmark overlap rule (≥ 0.5 of either locus) to match loci between replicates. Bombyx has one replicate, so its test report uses the negative-control proxy only. The "is real" probability is judged on ranking (both proxies) and on calibration (ECE against annotated positives); its operating threshold is chosen on tuning data. **[Claude's call]**

### 7.3 Success criterion (fixed before any run)

On tuning data, under the CV of §6:

- no class loses more than 0.01 precision or recall versus the rules;
- the share of ambiguous loci falls;
- on both "is real" proxies, the model does better than the rules.

Among algorithms that pass, the winner has the best mean macro-F1 over classes on annotated loci across folds and seeds. If two are within the fold-to-fold variation, the simpler wins. The comparison uses a paired test over folds that accounts for overlapping training sets (corrected resampled t-test; Nadeau & Bengio, 2003). **[test choice: Claude's call]**

Then the winner, trained on all tuning data, runs once on Bombyx and maize.

## 8. Algorithms

All families, one protocol:

- **Unsupervised (features only; annotation only to name groups and to evaluate):** k-means, Gaussian mixture, HDBSCAN, hierarchical clustering. Groups are mapped to classes by majority annotated label on tuning data; HDBSCAN noise feeds the "is real" score.
- **Semi-supervised (few labels + unlabeled loci):** label propagation/spreading, self-training, PU learning for "is real" (annotated = positive, rest = unlabeled; Elkan & Noto, 2008; Bekker & Davis, 2020).
- **Supervised (annotated loci only), as reference:** logistic regression, random forest, gradient boosting (XGBoost/LightGBM), SVM.

The list may grow; the protocol must make adding an algorithm a small change. **[implementation layout: builder's call]**

## 9. Deliverables

1. `--classifier {rules,ml}` in `classify` (default `rules`). The ml mode loads a model shipped with the package, versioned with its feature schema; it refuses to run on a mismatched schema.
2. Output gains a column with the probability of being real; the class column keeps its meaning. **[column name: builder's call]**
3. A comparison table of every algorithm against the rules, in the style of `benchmark/RESULTS-run6.md`, with references.
4. ML dependencies (scikit-learn, xgboost/lightgbm, hdbscan) are optional; the rules mode must work without them. On `litrp3` (56 cores, 1 TB RAM) they are not installed yet.
5. A retrain script in the repository (feature table → model artifact), so the model is reproducible and can be retrained when species or context features are added. One artifact per rules profile; a context variant is added with the annotation module. The ml mode must not make `classify` noticeably slower than the feature extraction it already does. **[Claude's call]**

## 10. Order of work

1. This ML stage, with BAM-derived features and the reserved context columns empty.
2. The annotation module (tRNAscan-SE, Infernal + Rfam without miRNA families, Red or RepeatModeler2 + RepeatMasker; user GFF when given).
3. ML again with context features, reporting the gain with versus without context.

## Referências

Bekker, J., & Davis, J. (2020). Learning from positive and unlabeled data: a survey. *Machine Learning*, 109(4), 719–760. https://doi.org/10.1007/s10994-020-05877-5

Elkan, C., & Noto, K. (2008). Learning classifiers from only positive and unlabeled data. In *Proceedings of the 14th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* (pp. 213–220). https://doi.org/10.1145/1401890.1401920

Nadeau, C., & Bengio, Y. (2003). Inference for the generalization error. *Machine Learning*, 52(3), 239–281. https://doi.org/10.1023/A:1024068626366
