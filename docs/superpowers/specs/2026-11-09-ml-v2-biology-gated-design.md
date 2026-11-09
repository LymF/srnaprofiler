# ML classifier v2: biology-gated decisions — design

Status: draft 2026-11-09, following the user's direction "temos que pegar o que classificou e ver a biologia dele pra ver se faz sentido a classificacao e ir melhorando o ML" and "segue com tudo". Extends `2026-11-07-ml-classifier-design.md`; everything there holds unless changed here. Items marked **[Claude's call]** were not decided by the user.

## 1. Why a v2

The v1 audits showed that the annotation-based metrics hide two opposite problems:

1. **The annotation is incomplete, so v1 is under-credited.** Fly (tuning, out-of-fold random forest): 24,369 unannotated loci called piRNA have the piRNA signature of the annotated ones (24–31 nt share 0.88, 1U 0.74, ping-pong z 7.3, transposon overlap). The piRNA precision loss that failed selection comes mostly from 109 miRNA-like loci (22 nt, 1U 0.94, no ping-pong) that v1 had to put somewhere once it could no longer output miRNA (R23).
2. **v1 makes biologically impossible calls on unannotated loci, so it is over-credited.** Maize (spent test, audited for diagnosis only): of ~10,600 ML phasiRNA calls, 5,113 sit on loci the rules leave ambiguous with no phasing at all (median phase p-value 1.0, 0 phase cycles, 21-nt share 0.22); of ~9,000 ML phasiRNA24 calls, 3,300 are unphased 24-nt hcsiRNA-like loci and 4,197 have no 24-nt majority. Annotation precision never sees these, because they are unannotated.
3. **v1 never abstains:** the threshold rule picked t = 0 in every fold.
4. **The PU "is real" score is not usable:** it scores 17,246 piRNA-like fly loci below 0.5. PU learning assumes the annotated loci are a random sample of the real ones (Bekker & Davis, 2020); here the annotation favours large, well-known clusters.

## 2. Changes

### 2.1 Biology gate (the core change)

The ML may only give a locus a class whose **defining biology the locus shows**. After the model's argmax, the class goes through the same checks the rules already apply to their own winner, from the same rules YAML, with no new per-class thresholds:

- the dominant size is inside the class `size` window (`_outside_window`);
- `cluster_min_length` / `cluster_max_density_pvalue` when the class has them (`_cluster_failures`);
- `phase_max_pvalue` when the class has it (`_phase_supported`);
- the class's `defining_evidence` holds: min over the defining evidence strengths (the rules' `score`) ≥ **0.5**. **[Claude's call: 0.5 is the midpoint of the evidence sigmoids, i.e. "the evidence is more present than absent"; fixed now, not tuned.]**

A locus that fails the gate for the predicted class gets the **next class by probability that passes**; if none passes, it is **ambiguous**. **[Claude's call: fall through rather than straight to ambiguous, so a 24-nt unphased locus predicted phasiRNA24 can still become hcsiRNA.]**

The gate lives in the pipeline (`decider.py` exposes it, `mlclassify.py` and the CV use it), not in benchmark post-processing.

### 2.2 An explicit "other" target, so the ML can abstain

Non-carved loci labelled miRNA (miRDeep2 missed them) and negative-control loci (rRNA / protein-coding exon, §7.2 of v1) become one training target, `other`. An `other` prediction is reported as ambiguous. Control loci therefore enter training; the control proxy is then evaluated only on held-out chromosomes (the CV already holds out whole chromosomes). **[Claude's call]**

The confidence threshold is kept (§6 of v1) but it is no longer the only way to abstain.

### 2.3 "Is real" without the PU assumption

`p_real` becomes a supervised score: positives = annotated loci of any class; negatives = control loci; fitted inside each fold like everything else. It is judged, as before, by replicate reproducibility (never a feature of the main model) and by control loci on held-out chromosomes. The PU version is kept only as a reported comparison. **[Claude's call]**

### 2.4 Biological audit as a reported axis

Every CV run also reports, per class, the share of **unannotated** calls that pass evidence **not used by the gate or the features**:

- arthropod piRNA: transposon overlap (RepeatMasker; UCSC dm6 for fly) and replicate reproducibility;
- all classes: replicate reproducibility; negative-control rate.

Signature summaries (size, 1U/5'A, ping-pong, phasing, strand) are reported too, but labelled as **consistency, not validation**, because the gate and the features use them.

### 2.5 Not in v2

- A `mixed` siRNA+piRNA class (fly: 1,178 loci with both signatures on transposons): reported in the audit, not modelled yet.
- ML recovery of miRNAs missed by miRDeep2: future work (user, 2026-11-07).
- tRNA/snoRNA/rRNA classes: wait for the annotation module (Ensembl maize GFF3 has no ncRNA genes).

## 3. Selection criterion v2 — fixed before the v2 run

User decision 2026-11-09 ("Anotação + biologia"). On tuning data, under the v1 CV, an algorithm is eligible when all hold:

1. **Annotation, overall:** mean macro-F1 is higher than the rules', corrected resampled t-test p < 0.05 (Nadeau & Bengio, 2003).
2. **Annotation, per class:** no class loses more than **0.02** precision or recall versus the rules (v1 used 0.01; the change is declared in the paper as a v2 decision taken after the v1 audit showed the annotation is incomplete).
3. **Independent biology of unannotated calls:** the share of the model's unannotated non-ambiguous calls that are *supported* is at least the rules' share, computed the same way on the rules' unannotated calls. *Supported* is fixed now: the locus recurs in at least half of the other replicates of its group (`repro_frac` ≥ 0.5), and, for arthropod piRNA, it also overlaps a RepeatMasker transposon by ≥ 1 bp. **[Claude's call: the exact definition of "supported".]**
4. **Controls:** the model calls fewer held-out control loci real than the rules do (unchanged).
5. **Ambiguous share falls** (unchanged).

Among eligible algorithms the winner has the best mean macro-F1; ties within the fold-to-fold variation go to the simpler one (unchanged). No relaxing after the run: if none is eligible, that is the result.

## 4. Data

- Tuning: unchanged (fly, Arabidopsis, rice). CV, folds, seeds, algorithms: unchanged.
- **Maize is spent** (used once by v1, then audited). It cannot be the independent test of v2. Bombyx was never run with an ML model and stays a valid held-out test for the arthropod profile. A plant held-out test needs a new species (from the annotation survey, memory `srna-annotation-sources`) chosen before v2 results are seen.

## Referências

Bekker, J., & Davis, J. (2020). Learning from positive and unlabeled data: a survey. *Machine Learning*, 109(4), 719–760. https://doi.org/10.1007/s10994-020-05877-5

Brennecke, J., Aravin, A. A., Stark, A., Dus, M., Kellis, M., Sachidanandam, R., & Hannon, G. J. (2007). Discrete small RNA-generating loci as master regulators of transposon activity in *Drosophila*. *Cell*, 128(6), 1089–1103. https://doi.org/10.1016/j.cell.2007.01.043

Rosenkranz, D., & Zischler, H. (2012). proTRAC — a software for probabilistic piRNA cluster detection, visualization and analysis. *BMC Bioinformatics*, 13, 5. https://doi.org/10.1186/1471-2105-13-5

Nadeau, C., & Bengio, Y. (2003). Inference for the generalization error. *Machine Learning*, 52(3), 239–281. https://doi.org/10.1023/A:1024068626366
