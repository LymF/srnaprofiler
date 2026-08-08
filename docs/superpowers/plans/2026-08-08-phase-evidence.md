# Phase evidence by effect size — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `e_phase` measure how many of the window's phase cycles are in register, so the `score` of a `phasiRNA` locus reads the quantity that actually separates real phasing from background.

**Architecture:** `signals.phase_test` already computes the number of cycles in the dominant register (`PhaseResult.number`) and is not touched. Only the line in `features.build_record` that turns the `PhaseResult` into `r.e_phase` changes, plus a named constant for the window size so the denominator is not a magic number in two files. The p-value path (`phase_pvalue`, `s_phase`) is untouched by design — one acceptance criterion is that `phase_pvalue` comes out byte for byte identical.

**Tech Stack:** Python 3, pytest, pysam, pandas/numpy (validation scripts only). Environment: `~/miniconda3/envs/srnaprofiler/bin/python`.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-08-08-phase-evidence-design.md`. Every number quoted below comes from it.
- New formula, verbatim: `e_phase = effect_strength(k = min(phase_number, CYCLES), n = CYCLES, p0 = 1/phase_length)` with `CYCLES = 11`.
- **Do not modify** `phase_test`, `phase_metrics`, `hypergeometric_phase_pvalue`, `_phase_window`, or `tests/test_signals_phase_null.py`.
- **Do not modify** the reported columns `phase_pvalue`, `s_phase`, `phase_score`, `phase_ratio`, `phase_number`, nor `w_phase` or any class decision.
- All code, comments, docstrings and log messages in English. Documentation in `docs/` stays in Portuguese, matching the files already there.
- Commit messages: imperative mood, no AI attribution, no `Co-Authored-By`.
- Run tests with `~/miniconda3/envs/srnaprofiler/bin/python -m pytest`. The whole suite is 334 tests and must stay green.
- Outputs of validation runs go under `/mnt/d/srnaprofiler-teste/`.

---

## File Structure

| file | responsibility | change |
|---|---|---|
| `signals.py` | phasing statistics | add the module constant `PHASE_WINDOW_CYCLES = 11` and use it as `phase_test`'s default, so the denominator has one source of truth |
| `features.py` | per-locus record building | the `r.e_phase` assignment (currently `features.py:204`) |
| `tests/test_features_phase_evidence.py` | new | what `e_phase` must mean |
| `tests/test_features_effect_strength.py` | existing | one assertion tightened |
| `docs/FORMULAS.md` | formula reference | the `e_phase` row and its caveat |
| `extras/validate_phase_evidence.py` | new | AUC of `score` against the truth sets, so acceptance criterion 1 is a command, not prose |

---

### Task 1: One source of truth for the window size

`phase_test` takes `cycles=11` as a default argument. `e_phase` needs the same 11 as its denominator, and a second literal in another file is how the two silently drift apart.

**Files:**
- Modify: `signals.py:83` (the `phase_test` signature)
- Test: `tests/test_signals_phase_test.py`

**Interfaces:**
- Produces: `signals.PHASE_WINDOW_CYCLES` (int, 11) — Task 2 imports it.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_signals_phase_test.py`:

```python
def test_the_window_size_is_a_named_constant():
    """e_phase divides by it, so it cannot be a literal in two files."""
    import inspect

    from signals import PHASE_WINDOW_CYCLES, phase_test

    assert PHASE_WINDOW_CYCLES == 11
    default = inspect.signature(phase_test).parameters["cycles"].default
    assert default == PHASE_WINDOW_CYCLES
```

- [ ] **Step 2: Run test to verify it fails**

Run: `~/miniconda3/envs/srnaprofiler/bin/python -m pytest tests/test_signals_phase_test.py::test_the_window_size_is_a_named_constant -v`
Expected: FAIL with `ImportError: cannot import name 'PHASE_WINDOW_CYCLES' from 'signals'`

- [ ] **Step 3: Write minimal implementation**

In `signals.py`, above `def phase_metrics` (around line 50), add:

```python
# Windows of this many phase cycles. e_phase divides by it: "of the cycles the
# test can see, how many are in register".
PHASE_WINDOW_CYCLES = 11
```

and change the signature at `signals.py:83` from

```python
def phase_test(positions, region_start, phase_length, cycles=11):
```

to

```python
def phase_test(positions, region_start, phase_length, cycles=PHASE_WINDOW_CYCLES):
```

Change nothing else in the function.

- [ ] **Step 4: Run the whole file to verify nothing shifted**

Run: `~/miniconda3/envs/srnaprofiler/bin/python -m pytest tests/test_signals_phase_test.py tests/test_signals_phase_null.py tests/test_signals_phasing.py -v`
Expected: PASS, all of them.

- [ ] **Step 5: Commit**

```bash
git add signals.py tests/test_signals_phase_test.py
git commit -m "refactor(signals): name the phase window size, e_phase divides by it"
```

---

### Task 2: e_phase counts cycles in register

**Files:**
- Modify: `features.py:204` (the `r.e_phase` assignment) and the `signals` import at `features.py:7`
- Test: `tests/test_features_phase_evidence.py` (create)

**Interfaces:**
- Consumes: `signals.PHASE_WINDOW_CYCLES` from Task 1; `signals.phase_test` returning `PhaseResult(score, number, ratio, pvalue, draws, observed, null_frac)`; `evidence.effect_strength(k, n, p0, conf=0.95) -> float` in `[0, 1]`.
- Produces: `ClusterRecord.e_phase` on the new definition. No other field changes.

- [ ] **Step 1: Write the failing test**

Create `tests/test_features_phase_evidence.py`:

```python
"""e_phase measures how far the phasing runs, not how surprising it is.

Measured on the 8 experimental TAS (12 loci) and 68 loci on annotated PHAS21,
against controls matched by distinct_seqs: the count of cycles in register is
the best discriminator available (AUC 0.826 / 0.769), the current p-value is at
chance (0.560 / 0.481) and the abundance ratio is worse than a coin (0.415).
Dividing the count by the locus's own occupied cycles destroys it (0.647),
because the positives simply span more cycles than their controls.
"""
from evidence import effect_strength
from features import build_record
from samio import Aln
from signals import PHASE_WINDOW_CYCLES

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}

_TAILS = ["ACGT", "TGCA", "GTAC", "CATG", "AGCT", "TCGA", "GATC", "CTAG",
          "AACC", "GGTT", "CCAA", "TTGG", "ACAC", "GTGT", "CACA", "TGTG"]


def _seqs(n, length=21):
    return [(_TAILS[i % len(_TAILS)] * ((length // 4) + 1))[:length]
            for i in range(n)]


def _in_register(n_cycles, phase=21, start=0):
    """One 21-nt read per phase cycle, all in the same register."""
    return [Aln("c", start + i * phase, "+", seq, 1)
            for i, seq in enumerate(_seqs(n_cycles, phase))]


def test_a_full_window_in_register_is_the_maximum():
    r = build_record(("c", 0, 400), _in_register(PHASE_WINDOW_CYCLES),
                     1_000_000, CFG)
    assert r.dicer_call == "21"
    assert r.e_phase == effect_strength(PHASE_WINDOW_CYCLES,
                                        PHASE_WINDOW_CYCLES, 1 / 21)
    assert r.e_phase > 0.6


def test_more_cycles_in_register_is_more_evidence():
    """The whole point: the statistic has to grow with the run's length."""
    short = build_record(("c", 0, 400), _in_register(3), 1_000_000, CFG)
    long = build_record(("c", 0, 400), _in_register(9), 1_000_000, CFG)
    assert long.e_phase > short.e_phase


def test_the_denominator_is_the_window_not_the_locus():
    """A 3-cycle locus must not score like a 3-cycle run inside a long one.

    Normalising by the locus's own cycles was measured and rejected: it drops
    the AUC from 0.826 to 0.647, because a TAS spans 10 cycles and its matched
    control spans 5.
    """
    r = build_record(("c", 0, 400), _in_register(3), 1_000_000, CFG)
    assert r.e_phase == effect_strength(3, PHASE_WINDOW_CYCLES, 1 / 21)
    assert r.e_phase < 0.3


def test_abundance_alone_does_not_make_evidence():
    """One hot position with 40 reads is stacking, not processivity."""
    stacked = [Aln("c", 100, "+", seq, 1) for seq in _seqs(40)]
    r = build_record(("c", 0, 400), stacked, 1_000_000, CFG)
    assert r.e_phase < 0.3


def test_an_empty_locus_has_no_phase_evidence():
    r = build_record(("c", 0, 400), [], 1_000_000, CFG)
    assert r.e_phase == 0.0


def test_a_locus_with_no_dicer_call_has_no_phase_evidence():
    """Sizes all over the window: there is no phase length to speak of."""
    mixed = [Aln("c", i * 30, "+", seq[:18 + (i % 12)], 1)
             for i, seq in enumerate(_seqs(20, 35))]
    r = build_record(("c", 0, 900), mixed, 1_000_000, CFG)
    assert r.dicer_call == "N"
    assert r.e_phase == 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `~/miniconda3/envs/srnaprofiler/bin/python -m pytest tests/test_features_phase_evidence.py -v`
Expected: FAIL. `test_a_full_window_in_register_is_the_maximum` and `test_the_denominator_is_the_window_not_the_locus` fail on the equality, because `e_phase` is still `effect_strength(res.observed, res.draws, res.null_frac)`.

- [ ] **Step 3: Write minimal implementation**

In `features.py:7`, extend the import:

```python
from signals import (offset_spectrum, pingpong_zscore, phase_test,
                     PHASE_WINDOW_CYCLES)
```

Replace `features.py:204`:

```python
        r.e_phase = effect_strength(res.observed, res.draws, res.null_frac)
```

with:

```python
        # Cycles in register, out of the cycles the window can hold. The run's
        # length is the processive cut's signature; the p-value of the same
        # count is at chance against size-matched controls, and the abundance
        # in the register is worse than that.
        r.e_phase = effect_strength(min(res.number, PHASE_WINDOW_CYCLES),
                                    PHASE_WINDOW_CYCLES, 1 / pl)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `~/miniconda3/envs/srnaprofiler/bin/python -m pytest tests/test_features_phase_evidence.py -v`
Expected: PASS, 6 tests.

- [ ] **Step 5: Run the whole suite**

Run: `~/miniconda3/envs/srnaprofiler/bin/python -m pytest tests/ -q`
Expected: PASS. If `tests/test_features_effect_strength.py::test_a_clean_phase_register_gets_a_high_phase_strength` fails, do **not** weaken it — it builds 11 reads in perfect register and asserts `e_phase > 0.3`, which the new formula satisfies (`effect_strength(11, 11, 1/21)` ≈ 0.70). A failure there means the implementation is wrong.

- [ ] **Step 6: Commit**

```bash
git add features.py tests/test_features_phase_evidence.py
git commit -m "fix(features): phase evidence is the length of the run, not the p-value

Measured against the 8 experimental TAS (12 loci) and 68 loci on annotated
PHAS21, with controls matched by distinct_seqs: the count of cycles in
register discriminates at AUC 0.826/0.769, the p-value of the same data at
0.560/0.481 -- chance -- and e_phase as a proportion of the occupied positions
at 0.618/0.639.

The denominator is the window's 11 cycles, not the locus's own: normalising by
the locus drops the AUC to 0.647, because a TAS spans 10 cycles where its
matched control spans 5."
```

---

### Task 3: Tighten the existing effect-strength test

`tests/test_features_effect_strength.py:83` asserts only `e_phase > 0.3` for a perfect 11-cycle register. Under the new definition that case is the maximum the statistic can reach, and the test should say so — otherwise it passes for a formula that gets it half right.

**Files:**
- Modify: `tests/test_features_effect_strength.py:82-89`

**Interfaces:**
- Consumes: `signals.PHASE_WINDOW_CYCLES`, `evidence.effect_strength` (Task 1, Task 2).

- [ ] **Step 1: Rewrite the assertion**

Replace the body of `test_a_clean_phase_register_gets_a_high_phase_strength` with:

```python
def test_a_clean_phase_register_gets_a_high_phase_strength():
    # 11 reads of 21 nt, one per cycle, in perfect register -> dicer_call "21".
    from signals import PHASE_WINDOW_CYCLES

    reads = [Aln("c", i * 21, "+", seq, 1)
             for i, seq in enumerate(_distinct_seqs(11, 21))]
    r = build_record(("c", 0, 400), reads, 1_000_000, CFG)
    assert r.dicer_call == "21"
    # A full window in register is as much phase evidence as there is.
    assert r.e_phase == effect_strength(PHASE_WINDOW_CYCLES,
                                        PHASE_WINDOW_CYCLES, 1 / 21)
```

Add `from evidence import effect_strength` to the file's imports if it is not already there.

- [ ] **Step 2: Run the file**

Run: `~/miniconda3/envs/srnaprofiler/bin/python -m pytest tests/test_features_effect_strength.py -v`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_features_effect_strength.py
git commit -m "test(features): a full window in register is the ceiling of e_phase"
```

---

### Task 4: The acceptance criterion as a command

Criterion 1 of the spec is an AUC against truth sets. Prose cannot be re-run; this script can.

**Files:**
- Create: `extras/validate_phase_evidence.py`
- Test: none — this is a validation script over real data, in the same family as `extras/validate_plant_tas.py` and `extras/validate_structural.py`, none of which are unit-tested.

**Interfaces:**
- Consumes: a `.clusters.tsv` from `classify`, `extras/phas_experimental_at_TAIR10.bed`, `extras/phas21_ref_at_TAIR10.gff3`.
- Produces: nothing importable; prints AUCs.

- [ ] **Step 1: Write the script**

Create `extras/validate_phase_evidence.py`:

```python
#!/usr/bin/env python3
"""Does the score separate real phasing from background of the same size?

Positives: the 8 experimentally characterised TAS loci, and the loci on PHAS21
annotated by sRNAminer. Controls: 21-nt loci on neither, matched to each
positive by distinct_seqs -- unmatched controls would be won by locus size
alone, since the dominant register's mass fraction is 0.954 in loci with <=20
distinct sequences against 1/21 expected.

There is no threshold here that counts as passing. The number is for comparing
formulations: the p-value scored 0.560 (TAS) and 0.481 (PHAS), i.e. chance.

Usage:  python3 extras/validate_phase_evidence.py run.clusters.tsv
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TAS_BED = os.path.join(HERE, "phas_experimental_at_TAIR10.bed")
PHAS_GFF = os.path.join(HERE, "phas21_ref_at_TAIR10.gff3")

# The PHAS annotation is in TAIR names, the runs are on the RefSeq assembly.
TAIR_TO_REFSEQ = {"Chr1": "NC_003070.9", "Chr2": "NC_003071.7",
                  "Chr3": "NC_003074.8", "Chr4": "NC_003075.7",
                  "Chr5": "NC_003076.8"}


def load_intervals():
    tas = pd.read_csv(TAS_BED, sep="\t", comment="#",
                      names=["chrom", "start", "end", "name"])
    phas = []
    with open(PHAS_GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                phas.append((TAIR_TO_REFSEQ.get(f[0], f[0]),
                             int(f[3]) - 1, int(f[4])))
    return tas, pd.DataFrame(phas, columns=["chrom", "start", "end"])


def overlaps(d, iv):
    hit = np.zeros(len(d), dtype=bool)
    for chrom, grp in iv.groupby("chrom"):
        same = (d.chrom == chrom).values
        if not same.any():
            continue
        for _, r in grp.iterrows():
            hit |= same & (d.end.values > r.start) & (d.start.values < r.end)
    return hit


def auc(pos, neg):
    """Probability a positive outranks a control, ties counted as half."""
    a = np.asarray(pos, dtype=float)
    b = np.asarray(neg, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    wins = (a[:, None] > b[None, :]).sum()
    ties = (a[:, None] == b[None, :]).sum()
    return (wins + 0.5 * ties) / (len(a) * len(b))


def main(tsv_path):
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    d = d[d.dicer_call.astype(str) == "21"].reset_index(drop=True)
    tas, phas = load_intervals()
    d["in_tas"] = overlaps(d, tas)
    d["in_phas"] = overlaps(d, phas) | d.in_tas

    pool = d[~d.in_phas]
    negs = []
    for i in d.index[d.in_phas]:
        n = d.at[i, "distinct_seqs"]
        cand = pool[(pool.distinct_seqs >= n * 0.5)
                    & (pool.distinct_seqs <= n * 1.5)]
        cand = cand[~cand.index.isin(negs)]
        if len(cand):
            negs.append(int(cand.sample(1, random_state=len(negs)).index[0]))

    ctrl = d.loc[negs]
    print(f"{len(d)} loci 21-nt; {d.in_tas.sum()} on experimental TAS, "
          f"{d.in_phas.sum()} on annotated PHAS; {len(ctrl)} matched controls")
    print(f"\n{'column':14} {'AUC TAS':>8} {'AUC PHAS':>9}")
    for col in ("score", "s_phase", "phase_number", "phase_ratio"):
        if col not in d.columns:
            continue
        print(f"{col:14} {auc(d[d.in_tas][col], ctrl[col]):8.3f} "
              f"{auc(d[d.in_phas][col], ctrl[col]):9.3f}")
    print("\nmedian score  TAS / PHAS / control: "
          f"{d[d.in_tas].score.median():.3f} / "
          f"{d[d.in_phas].score.median():.3f} / {ctrl.score.median():.3f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])
```

- [ ] **Step 2: Run it against the pre-change run, to record the baseline**

Run:
```bash
~/miniconda3/envs/srnaprofiler/bin/python extras/validate_phase_evidence.py \
    /mnt/d/srnaprofiler-teste/plant/classify/at_ncrna_org.clusters.tsv
```
Expected: it prints AUCs. `score` should come out near **0.618 (TAS)** and **0.639 (PHAS)** — this file was produced before the change. If it does not, stop: either the script disagrees with the experiment it reproduces, or the file is not the one the spec measured.

- [ ] **Step 3: Commit**

```bash
git add extras/validate_phase_evidence.py
git commit -m "feat(extras): AUC of the score against phasing truth sets

The spec's first acceptance criterion was a number in prose. This runs it."
```

---

### Task 5: Rerun both organisms and check every acceptance criterion

**Files:**
- No source changes. Produces `/mnt/d/srnaprofiler-teste/plant/classify/at_ephase.clusters.tsv` and `/mnt/d/srnaprofiler-teste/dros/classify/dros_ephase.clusters.tsv`.

**Interfaces:**
- Consumes: the sorted BAMs from the ncRNA-filtered runs, which already exist.

- [ ] **Step 1: Rerun Arabidopsis, on the same BAM and the same denominator as the baseline**

```bash
cd /home/lucas/srnaprofiler
export PATH=~/miniconda3/envs/srnaprofiler/bin:$PATH
~/miniconda3/envs/srnaprofiler/bin/python srnaprofiler.py classify \
    --sam /mnt/d/srnaprofiler-teste/plant/classify/at_ncrna.sorted.bam \
    -g /mnt/d/srnaprofiler-teste/plant/GCF_000001735.4_TAIR10.1_genomic.fna \
    -o /mnt/d/srnaprofiler-teste/plant/classify -p at_ephase \
    --config rules.plant.yaml \
    --mirna-loci /mnt/d/srnaprofiler-teste/refs/ath_mirbase.bed \
    --organelle-contigs NC_000932.1,NC_037304.1 \
    --norm 19934660 --threads 1
```
Takes about 10 minutes. The comparison baseline is `at_ncrna_org.clusters.tsv`, produced with exactly these arguments.

- [ ] **Step 2: Criterion 4 — the p-value must not have moved**

```bash
~/miniconda3/envs/srnaprofiler/bin/python - <<'EOF'
import pandas as pd
b = "/mnt/d/srnaprofiler-teste/plant/classify/"
old = pd.read_csv(b + "at_ncrna_org.clusters.tsv", sep="\t", low_memory=False)
new = pd.read_csv(b + "at_ephase.clusters.tsv", sep="\t", low_memory=False)
assert len(old) == len(new), (len(old), len(new))
for col in ("phase_pvalue", "s_phase", "phase_score", "phase_ratio",
            "phase_number", "cls"):
    same = (old[col].fillna(-1).values == new[col].fillna(-1).values).mean()
    print(f"{col:14} identical in {same*100:.2f}% of loci")
    assert same == 1.0, f"{col} moved -- only e_phase was supposed to change"
print("OK: only the score changed")
EOF
```
Expected: every column identical in 100% of loci, including `cls` (criterion 3, zero class changes).

- [ ] **Step 3: Criterion 1 — discrimination**

```bash
~/miniconda3/envs/srnaprofiler/bin/python extras/validate_phase_evidence.py \
    /mnt/d/srnaprofiler-teste/plant/classify/at_ephase.clusters.tsv
```
Expected: `score` AUC **≥ 0.80 on TAS** (was 0.618) and **≥ 0.75 on PHAS** (was 0.639).

If it lands below, do not adjust the threshold: report the number. The spec's own denominator experiments are in `extras/experiment_phase_statistic.py` and the next step would be to rerun those, not to move the goalposts.

- [ ] **Step 4: Criterion 2 — the plant floor**

```bash
P=~/miniconda3/envs/srnaprofiler/bin/python
T=/mnt/d/srnaprofiler-teste/plant/classify/at_ephase.clusters.tsv
$P extras/validate_plant_tas.py $T \
    /mnt/d/srnaprofiler-teste/plant/GCF_000001735.4_TAIR10.1_genomic.gff | grep -E "/6"
$P extras/validate_mirna.py $T ath.gff3 | sed -n '3,7p'
```
Expected: `TAS loci with at least one called locus: 6/6`, and 242/242 precursors labelled miRNA.

- [ ] **Step 5: Rerun Drosophila and check its floor**

```bash
cd /home/lucas/srnaprofiler
~/miniconda3/envs/srnaprofiler/bin/python srnaprofiler.py classify \
    --sam /mnt/d/srnaprofiler-teste/dros/classify/dros_ncrna.sorted.bam \
    -g /mnt/d/srnaprofiler-teste/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna \
    -o /mnt/d/srnaprofiler-teste/dros/classify -p dros_ephase \
    --config rules.yaml \
    --mirna-loci /mnt/d/srnaprofiler-teste/refs/dme_mirbase.bed \
    --norm 9658002 --threads 2

P=~/miniconda3/envs/srnaprofiler/bin/python
T=/mnt/d/srnaprofiler-teste/dros/classify/dros_ephase.clusters.tsv
$P extras/validate_piclusters.py $T | grep -E "RECALL"
$P extras/validate_mirna.py $T dme.gff3 | sed -n '3,7p'
```
Expected: piC 85/90 and miRBase 134/134, unchanged from `dros_ncrna_matched.clusters.tsv`.

Drosophila's `defining_evidence` for piRNA is `[[1u, pingpong]]`, so no locus there should change score at all — a difference here means `e_phase` leaked into a class that does not use it.

- [ ] **Step 6: Criterion 5 — calibration untouched**

```bash
~/miniconda3/envs/srnaprofiler/bin/python -m pytest tests/test_signals_phase_null.py -v
```
Expected: PASS, with no edit to that file. It tests `phase_test`, which this plan does not touch.

- [ ] **Step 7: Commit the measured result into the docs**

Update `docs/FORMULAS.md`:

- In the table at line 143, replace the `e_phase` row

```
| `e_phase` | `observados / sorteios` | `sorteios` | `sucessos / população` |
```

with

```
| `e_phase` | `ciclos em registro / 11` | `11` (a janela) | `1 / phase_length` |
```

- Replace the caveat at line 215

```
- **`e_phase` é otimista**: a janela é escolhida por maximizar evidência, e o Bonferroni corrige
  significância, não tamanho de efeito.
```

with

```
- **`e_phase` mede a extensão da fase, não a surpresa.** É a fração dos 11 ciclos da janela que
  está no registro dominante — a corrida de cortes processivos é a assinatura do mecanismo.
  Medido contra controles casados por `distinct_seqs`: AUC 0,826 (TAS experimental) e 0,769
  (PHAS anotado), contra 0,560/0,481 do p-valor dos mesmos dados e 0,415 da razão de
  abundância. Ressalvas: a janela é escolhida por maximizar evidência, e o `p₀ = 1/phase_length`
  é anticonservador porque o registro reportado é o máximo entre 21.
```

Then update `docs/ROADMAP.md`: in the section **A0b**, replace the "Precisa de spec própria" sentence with the outcome — the spec is `docs/superpowers/specs/2026-08-08-phase-evidence-design.md`, the premise that abundance was missing was measured and rejected, and the AUCs before and after.

```bash
git add docs/FORMULAS.md docs/ROADMAP.md
git commit -m "docs: e_phase measures how far the phasing runs, with the numbers"
```

---

## Self-Review

**Spec coverage:**

| spec requirement | task |
|---|---|
| `e_phase = effect_strength(min(number, CYCLES), CYCLES, 1/phase_length)` | Task 2 |
| `CYCLES = 11`, one source of truth | Task 1 |
| `phase_test` and the p-value path untouched | constraint; verified in Task 5 step 2 |
| `docs/FORMULAS.md` updated | Task 5 step 7 |
| Criterion 1, AUC ≥ 0.80 / 0.75 | Task 4 (tool) + Task 5 step 3 |
| Criterion 2, floors | Task 5 steps 4 and 5 |
| Criterion 3, zero class change | Task 5 step 2 (`cls` column) |
| Criterion 4, `phase_pvalue` identical | Task 5 step 2 |
| Criterion 5, permutation null passing unedited | Task 5 step 6 |
| Consecutive cycles — out of scope | not planned, by design |

**Placeholders:** none. Every step has the command or the code to paste.

**Type consistency:** `PHASE_WINDOW_CYCLES` is defined in Task 1 and used under that name in Tasks 2 and 3. `effect_strength(k, n, p0)` matches `evidence.py:24`. `res.number` is field 2 of `PhaseResult` as declared at `signals.py:76`. `pl` is the local `int(r.dicer_call)` already in scope at `features.py:197`.
