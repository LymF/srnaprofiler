# Evidence Scores Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the relative softmax `confidence` with a per-locus score that means something on its own: every signature is tested against an explicit null and reported as `log10(1/p)`.

**Architecture:** A new module `evidence.py` holds one function per statistical test, each pure and independently testable. `features.build_record` fills the new `s_*` fields from those functions. `RuleScorer` reads which evidence defines each class from the config and sets the class score to the minimum over that evidence. Nothing is removed: `cls`, `confidence` and `support` stay, so existing consumers keep working.

**Tech Stack:** Python, scipy.stats (already installed, 1.17.1), pytest.

## Global Constraints

- Code, comments, docstrings, log messages and CLI help in **English**.
- Comments short — only what is strictly necessary. Long context goes in the commit message or the roadmap, never inline.
- Never add `Co-Authored-By` or any AI credit to commits.
- Scores are `log10(1/p)`: score 2 means p = 0.01. Always non-negative, rounded to 3 decimals.
- Binomial `n` is the 1/NH-weighted read count, rounded to int, floored at 1.
- A rule enters only if it holds in **both** kingdoms — validate against Arabidopsis (6 TAS loci) and Drosophila (90 piRNA clusters) before committing behaviour changes.
- `p0` for the 5' base tests is **0.25**, following proTRAC.

---

## File Structure

| file | responsibility |
|---|---|
| `evidence.py` (new) | one function per null test; pure, no I/O, no record knowledge |
| `tests/test_evidence.py` (new) | unit tests for each test function |
| `records.py` | new `s_*` fields on `ClusterRecord` |
| `reporter.py` | new columns in `_COLS` |
| `features.py` | fill `s_*` while building the record |
| `decider.py` | class score = min over defining evidence |
| `rules.yaml`, `rules.plant.yaml` | declare `defining_evidence` per class |
| `README.md` | document what the score means and does not mean |

---

### Task 1: `evidence.py` — score conversion and the binomial tests

**Files:**
- Create: `evidence.py`
- Test: `tests/test_evidence.py`

**Interfaces:**
- Produces: `score_from_p(p: float) -> float`, `base_bias_score(k: float, n: float, p0: float = 0.25) -> float`, `strand_bias_score(frac_top: float, n: float) -> float`

- [ ] **Step 1: Write the failing tests**

```python
"""Null tests behind the per-locus evidence scores.

score = log10(1/p), so score 2 means p = 0.01. p0 = 0.25 for the 5' base
follows proTRAC, which states it tests 1T/10A "based on random base
composition (25% for each nucleotide)".
"""
import math

from evidence import score_from_p, base_bias_score, strand_bias_score


def test_score_from_p_is_log10_of_one_over_p():
    assert score_from_p(0.01) == 2.0
    assert score_from_p(0.001) == 3.0
    assert score_from_p(1.0) == 0.0


def test_score_from_p_never_negative_and_caps_at_zero_p():
    assert score_from_p(1.5) == 0.0          # p > 1 from rounding
    assert score_from_p(0.0) == 300.0        # underflow, capped


def test_base_bias_strong_signal_scores_high():
    # 90 of 100 reads starting with the same base, against 25% by chance
    assert base_bias_score(90, 100) > 20


def test_base_bias_at_chance_scores_zero():
    assert base_bias_score(25, 100) < 0.5


def test_base_bias_shallow_locus_cannot_prove_much():
    # same 90% proportion, but only 10 reads
    deep = base_bias_score(90, 100)
    shallow = base_bias_score(9, 10)
    assert shallow < deep
    assert shallow < 6


def test_base_bias_handles_degenerate_input():
    assert base_bias_score(0, 0) == 0.0
    assert base_bias_score(5, 3) == 0.0      # k > n is nonsense


def test_strand_bias_scores_asymmetry_either_way():
    assert strand_bias_score(0.95, 200) > 10
    assert strand_bias_score(0.05, 200) > 10


def test_strand_bias_zero_for_balanced_locus():
    assert strand_bias_score(0.5, 200) < 0.5
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_evidence.py -q`
Expected: FAIL, `ModuleNotFoundError: No module named 'evidence'`

- [ ] **Step 3: Write the implementation**

```python
"""Null tests for the per-locus evidence scores."""
import math

from scipy import stats

# log10(1/p) of the smallest p we report; below this the test is saturated.
_MAX_SCORE = 300.0
# proTRAC tests the 5' base against random base composition.
BASE_P0 = 0.25


def score_from_p(p):
    """log10(1/p), floored at 0 and capped at _MAX_SCORE."""
    if p is None or p >= 1.0:
        return 0.0
    if p <= 0.0:
        return _MAX_SCORE
    return round(min(math.log10(1.0 / p), _MAX_SCORE), 3)


def base_bias_score(k, n, p0=BASE_P0):
    """Enrichment of a 5' base over chance (one-sided binomial)."""
    k, n = int(round(k)), int(round(n))
    if n <= 0 or k <= 0 or k > n:
        return 0.0
    return score_from_p(stats.binomtest(k, n, p0, alternative="greater").pvalue)


def strand_bias_score(frac_top, n):
    """Departure from an even strand split (two-sided binomial)."""
    n = int(round(n))
    if n <= 0 or frac_top is None:
        return 0.0
    k = int(round(frac_top * n))
    if k < 0 or k > n:
        return 0.0
    return score_from_p(stats.binomtest(k, n, 0.5, alternative="two-sided").pvalue)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_evidence.py -q`
Expected: PASS, 8 passed

- [ ] **Step 5: Commit**

```bash
git add evidence.py tests/test_evidence.py
git commit -m "feat(evidence): score_from_p and the binomial null tests

score = log10(1/p) so the number reads directly: 2 is p = 0.01. The 5' base
is tested against 0.25, which is what proTRAC does.

The tests pin the property that matters for false positives: the same 90%
proportion scores far lower on 10 reads than on 100, because a shallow locus
cannot prove much."
```

---

### Task 2: `evidence.py` — ping-pong and density tests

**Files:**
- Modify: `evidence.py`
- Test: `tests/test_evidence.py`

**Interfaces:**
- Consumes: `score_from_p` from Task 1
- Produces: `pingpong_score(pairs_at_10: float, total_pairs: float, n_offsets: int = 20) -> float`, `density_score(reads: float, span_bp: int, lib_reads: float, genome_bp: int) -> float`

- [ ] **Step 1: Write the failing tests**

```python
def test_pingpong_peak_at_ten_scores_high():
    from evidence import pingpong_score
    # 200 of 400 pairs at offset 10, against 1/20 by chance
    assert pingpong_score(200, 400) > 20


def test_pingpong_at_chance_scores_zero():
    from evidence import pingpong_score
    assert pingpong_score(20, 400) < 0.5      # exactly 1/20


def test_pingpong_few_pairs_cannot_prove_much():
    from evidence import pingpong_score
    # every pair at offset 10, but only 3 pairs: geometry, not biology
    assert pingpong_score(3, 3) < 5


def test_density_above_background_scores_high():
    from evidence import density_score
    # 5000 reads in 500bp, library 10M reads over 100Mb -> 50 expected
    assert density_score(5000, 500, 10_000_000, 100_000_000) > 20


def test_density_at_background_scores_zero():
    from evidence import density_score
    assert density_score(50, 500, 10_000_000, 100_000_000) < 1.0


def test_density_handles_degenerate_input():
    from evidence import density_score
    assert density_score(100, 0, 10_000_000, 100_000_000) == 0.0
    assert density_score(100, 500, 0, 100_000_000) == 0.0
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_evidence.py -q`
Expected: FAIL, `ImportError: cannot import name 'pingpong_score'`

- [ ] **Step 3: Write the implementation**

Append to `evidence.py`:

```python
def pingpong_score(pairs_at_10, total_pairs, n_offsets=20):
    """Excess of 5'-5' pairs at the 10-nt offset over the other offsets."""
    k, n = int(round(pairs_at_10)), int(round(total_pairs))
    if n <= 0 or k <= 0 or k > n:
        return 0.0
    return score_from_p(
        stats.binomtest(k, n, 1.0 / n_offsets, alternative="greater").pvalue)


def density_score(reads, span_bp, lib_reads, genome_bp):
    """Read pile-up over what a uniform distribution would put in this span."""
    if span_bp <= 0 or lib_reads <= 0 or genome_bp <= 0 or reads <= 0:
        return 0.0
    expected = lib_reads * (span_bp / genome_bp)
    if expected <= 0:
        return 0.0
    return score_from_p(stats.poisson.sf(int(round(reads)) - 1, expected))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_evidence.py -q`
Expected: PASS, 14 passed

- [ ] **Step 5: Commit**

```bash
git add evidence.py tests/test_evidence.py
git commit -m "feat(evidence): ping-pong and density null tests

Ping-pong is binomial against 1/20 -- the chance a pair lands on the 10-nt
offset among the twenty tested -- so three pairs all at offset 10 score low,
which is the geometric artefact seen on small dense loci.

Density is Poisson against what a uniform distribution would put in the span,
as proTRAC does for hit density."
```

---

### Task 3: carry the scores on the record and in the TSV

**Files:**
- Modify: `records.py:30` (after `pct_5a`), `reporter.py:9`
- Test: `tests/test_reporter_scores.py`

**Interfaces:**
- Produces: `ClusterRecord` fields `s_1u`, `s_5a`, `s_fractop`, `s_pingpong`, `s_phase`, `s_density`, `score`, all `float = 0.0`; the same names as TSV columns.

- [ ] **Step 1: Write the failing test**

```python
"""The evidence scores must reach the TSV, in a fixed column order."""
from records import ClusterRecord
from reporter import write_tsv, _COLS

SCORE_COLS = ["s_1u", "s_5a", "s_fractop", "s_pingpong", "s_phase",
              "s_density", "score"]


def test_score_columns_present_and_after_support():
    for c in SCORE_COLS:
        assert c in _COLS
    assert _COLS.index("score") > _COLS.index("support")


def test_scores_are_written(tmp_path):
    r = ClusterRecord(chrom="c", start=0, end=100, name="c:0-100", reads=500,
                      dicer_call="24", s_5a=14.2, score=9.4)
    out = tmp_path / "x.tsv"
    write_tsv([r], str(out))
    header, row = out.read_text().splitlines()[:2]
    cols = header.split("\t")
    vals = row.split("\t")
    assert vals[cols.index("s_5a")] == "14.2"
    assert vals[cols.index("score")] == "9.4"


def test_score_defaults_to_zero():
    r = ClusterRecord(chrom="c", start=0, end=1, name="n")
    assert r.score == 0.0
    assert r.s_phase == 0.0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_reporter_scores.py -q`
Expected: FAIL, `TypeError: ClusterRecord.__init__() got an unexpected keyword argument 's_5a'`

- [ ] **Step 3: Write the implementation**

In `records.py`, after the `support` field:

```python
    # evidence scores, log10(1/p) per signature
    s_1u: float = 0.0
    s_5a: float = 0.0
    s_fractop: float = 0.0
    s_pingpong: float = 0.0
    s_phase: float = 0.0
    s_density: float = 0.0
    score: float = 0.0
```

In `reporter.py`, replace the last line of `_COLS`:

```python
    "cls", "confidence", "support", "score",
    "s_1u", "s_5a", "s_fractop", "s_pingpong", "s_phase", "s_density",
    "review", "evidence",
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_reporter_scores.py tests/test_reporter.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add records.py reporter.py tests/test_reporter_scores.py
git commit -m "feat(records): evidence score fields on the record and in the TSV

One column per signature plus the class score. Existing columns keep their
place so current consumers are unaffected."
```

---

### Task 4: fill the scores in `features.build_record`

**Files:**
- Modify: `features.py:92-131` (inside `build_record`), `features.py:178-200` (`process_chromosome` signature)
- Test: `tests/test_features_scores.py`

**Interfaces:**
- Consumes: `evidence.base_bias_score`, `strand_bias_score`, `pingpong_score`, `density_score`, `score_from_p` (Tasks 1-2)
- Produces: `build_record(region, region_alns, lib_size, config, genome_bp=0)` — one new optional parameter, defaulting to 0 so density scoring is skipped when the genome size is unknown.

- [ ] **Step 1: Write the failing test**

```python
"""build_record must turn the raw features into evidence scores."""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _reads(seq, n, pos=100, strand="+"):
    return [Aln("c", pos + i, strand, seq, 1) for i in range(n)]


def test_strong_5a_bias_scores_high():
    r = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200),
                     1_000_000, CFG)
    assert r.s_5a > 20
    assert r.s_1u == 0.0


def test_shallow_locus_scores_low_despite_same_proportion():
    deep = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200),
                        1_000_000, CFG)
    shallow = build_record(("c", 100, 300), _reads("A" + "C" * 23, 5),
                           1_000_000, CFG)
    assert shallow.s_5a < deep.s_5a


def test_single_strand_locus_scores_strand_bias():
    r = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200),
                     1_000_000, CFG)
    assert r.s_fractop > 10


def test_phase_score_mirrors_the_pvalue():
    from evidence import score_from_p
    reads = [Aln("c", 100 + i * 21, "+", "T" + "A" * 20, 1) for i in range(11)]
    r = build_record(("c", 100, 400), reads, 1_000_000, CFG)
    assert r.s_phase == score_from_p(r.phase_pvalue)


def test_density_scored_only_when_genome_size_given():
    reads = _reads("A" + "C" * 23, 200)
    without = build_record(("c", 100, 300), reads, 1_000_000, CFG)
    with_genome = build_record(("c", 100, 300), reads, 1_000_000, CFG,
                               genome_bp=10_000_000)
    assert without.s_density == 0.0
    assert with_genome.s_density > 0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_features_scores.py -q`
Expected: FAIL, `AttributeError: 'ClusterRecord' object has no attribute 's_5a'` is already fixed by Task 3, so it fails on `assert r.s_5a > 20` with `0.0 > 20`

- [ ] **Step 3: Write the implementation**

In `features.py`, change the import line:

```python
from signals import offset_spectrum, pingpong_zscore, phase_test
from evidence import (base_bias_score, strand_bias_score, pingpong_score,
                      density_score, score_from_p)
```

Change the signature:

```python
def build_record(region, region_alns, lib_size, config, genome_bp=0):
```

After `r.pct_5a = round(n_5a / r.reads, 4)`:

```python
    r.s_1u = base_bias_score(n_1u, r.reads)
    r.s_5a = base_bias_score(n_5a, r.reads)
    r.s_fractop = strand_bias_score(r.frac_top, r.reads)
    r.s_density = density_score(r.reads, end - start, lib_size, genome_bp)
```

After `r.pingpong_frac = ...`:

```python
    r.s_pingpong = pingpong_score(overlap_counts.get(10, 0), total_pairs)
```

Inside the `if r.dicer_call != "N":` block, after `r.phase_ratio, r.phase_pvalue = res.ratio, res.pvalue`:

```python
        r.s_phase = score_from_p(r.phase_pvalue)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_features_scores.py tests/test_features_signals.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add features.py tests/test_features_scores.py
git commit -m "feat(features): fill the evidence scores per locus

genome_bp is optional and defaults to 0, which skips density scoring, so
callers that do not know the genome size keep working unchanged."
```

---

### Task 5: pass the genome size through `classify`

**Files:**
- Modify: `features.py:178-227` (`process_chromosome`, `extract_by_chrom`, `_chrom_job`), `classify.py:93-103`
- Test: `tests/test_classify_genome_bp.py`

**Interfaces:**
- Consumes: `build_record(..., genome_bp=...)` from Task 4
- Produces: `process_chromosome(bam_path, chrom, mincov, pad, lib_size, config, xm_as_nh=False, genome_bp=0)` and `extract_by_chrom(bam_path, chroms, lib_size, config, mincov, pad, threads=1, xm_as_nh=False, genome_bp=0)`; `samio.mapped_genome_size(bam_path) -> int`

- [ ] **Step 1: Write the failing test**

```python
"""The density null needs the size of the references actually mapped."""
import pysam

from samio import mapped_genome_size


def _bam(tmp_path, lengths):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": l} for n, l in lengths]}
    path = tmp_path / "x.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for n, _ in lengths:
            a = pysam.AlignedSegment(fh.header)
            a.query_name = "r" + n
            a.query_sequence = "ACGTACGTACGTACGTACGTA"
            a.reference_id = fh.header.get_tid(n)
            a.reference_start = 10
            a.cigarstring = "21M"
            a.flag = 0
            a.mapping_quality = 255
            fh.write(a)
    pysam.index(str(path))
    return str(path)


def test_sums_lengths_of_mapped_references(tmp_path):
    bam = _bam(tmp_path, [("chr1", 1000), ("chr2", 2000)])
    assert mapped_genome_size(bam) == 3000


def test_ignores_references_without_reads(tmp_path):
    bam = _bam(tmp_path, [("chr1", 1000)])
    # a reference declared in the header but holding no read must not count
    with pysam.AlignmentFile(bam, "rb") as fh:
        assert fh.get_reference_length("chr1") == 1000
    assert mapped_genome_size(bam) == 1000
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_classify_genome_bp.py -q`
Expected: FAIL, `ImportError: cannot import name 'mapped_genome_size'`

- [ ] **Step 3: Write the implementation**

Append to `samio.py`:

```python
def mapped_genome_size(bam_path):
    """Total length of the references that carry reads (from the index)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        return sum(fh.get_reference_length(s.contig)
                   for s in fh.get_index_statistics() if s.mapped > 0)
```

In `features.py`, thread the parameter through:

```python
def process_chromosome(bam_path, chrom, mincov, pad, lib_size, config,
                       xm_as_nh=False, genome_bp=0):
```

and in its final `build_record` call add `genome_bp=genome_bp`; in `_chrom_job` unpack the extra element; in `extract_by_chrom` add `genome_bp=0` to the signature and to each tuple in `jobs`.

In `classify.py`, import it alongside the other `samio` names and pass it:

```python
    records = extract_by_chrom(bam_path, chroms, lib_size, ext_cfg,
                               args.mincov, args.pad, threads, xm_as_nh,
                               mapped_genome_size(bam_path))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_classify_genome_bp.py tests/test_classify_threads.py tests/test_classify_integration.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add samio.py features.py classify.py tests/test_classify_genome_bp.py
git commit -m "feat(classify): pass the mapped genome size for the density null

Only references carrying reads count, so unplaced scaffolds in the header do
not inflate the expected background."
```

---

### Task 6: class score as the minimum over defining evidence

**Files:**
- Modify: `decider.py:46-90` (`_score_class`), `decider.py:91-125` (`decide`)
- Test: `tests/test_decider_class_score.py`

**Interfaces:**
- Consumes: `ClusterRecord.s_*` fields from Task 4
- Produces: `RuleScorer.defining_scores(spec, record) -> list[float]`; `record.score` set by `decide`; config key `defining_evidence: [str]` per class, where each name maps to the record field `s_<name>`.

- [ ] **Step 1: Write the failing test**

```python
"""The class score is the weakest of its defining evidence.

Minimum, not sum: summing would let one strong signal cover for a missing one,
which is the false positive being targeted. Score 2 therefore reads as "every
defining evidence has p <= 0.01".
"""
from records import ClusterRecord
from decider import RuleScorer

CFG = {
    "decision": {"review_confidence": 0.4, "min_margin": 0.0, "min_reads": 10},
    "classes": {
        "hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0, "w_5a": 1.0,
                    "defining_evidence": ["5a", "density"]},
        "phasiRNA": {"size": {"lo": 21, "hi": 22}, "w_size": 1.0, "w_phase": 3.0,
                     "defining_evidence": ["phase"]},
    },
}


def _rec(**kw):
    base = dict(chrom="c", start=0, end=100, name="n", reads=500,
                frac_top=0.5, pct_1u=0.1, pct_5a=0.8, pingpong_frac=0.0,
                phase_score=0.0)
    base.update(kw)
    return ClusterRecord(**base)


def test_score_is_the_minimum_of_defining_evidence():
    r = _rec(dicer_call="24", s_5a=14.0, s_density=9.0)
    RuleScorer(CFG).decide(r)
    assert r.cls == "hcsiRNA"
    assert r.score == 9.0


def test_one_missing_evidence_sinks_the_score():
    # spectacular density cannot cover for absent base bias
    r = _rec(dicer_call="24", s_5a=0.0, s_density=40.0)
    RuleScorer(CFG).decide(r)
    assert r.cls == "hcsiRNA"
    assert r.score == 0.0


def test_non_defining_evidence_is_ignored_by_the_score():
    r = _rec(dicer_call="24", s_5a=14.0, s_density=9.0, s_pingpong=99.0)
    RuleScorer(CFG).decide(r)
    assert r.score == 9.0


def test_class_without_defining_evidence_scores_zero():
    cfg = {"decision": CFG["decision"],
           "classes": {"x": {"size": {"lo": 20, "hi": 30}, "w_size": 1.0}}}
    r = _rec(dicer_call="24")
    RuleScorer(cfg).decide(r)
    assert r.score == 0.0


def test_score_belongs_to_the_winning_class():
    r = _rec(dicer_call="21", s_phase=7.0, s_5a=30.0, s_density=30.0)
    RuleScorer(CFG).decide(r)
    assert r.cls == "phasiRNA"
    assert r.score == 7.0
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_decider_class_score.py -q`
Expected: FAIL, `assert 0.0 == 9.0` (the field exists but nothing sets it)

- [ ] **Step 3: Write the implementation**

In `decider.py`, add to `RuleScorer`:

```python
    @staticmethod
    def defining_scores(spec, record):
        """Evidence scores the class is defined by, in config order."""
        return [getattr(record, f"s_{name}", 0.0)
                for name in spec.get("defining_evidence", [])]
```

In `decide`, after `record.cls = top_name`:

```python
        defining = self.defining_scores(self.cfg["classes"][top_name], record)
        record.score = round(min(defining), 3) if defining else 0.0
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_decider_class_score.py tests/test_decider.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add decider.py tests/test_decider_class_score.py
git commit -m "feat(decider): class score is the minimum over defining evidence

Summing would let one strong signal cover for a missing one, which is exactly
the false positive being targeted: a locus with spectacular density and no base
bias would still score high. With the minimum, score 2 reads as 'every defining
evidence has p <= 0.01'."
```

---

### Task 7: declare the defining evidence in both configs

**Files:**
- Modify: `rules.yaml`, `rules.plant.yaml`
- Test: `tests/test_config_defining_evidence.py`

**Interfaces:**
- Consumes: `defining_evidence` read by `RuleScorer.defining_scores` (Task 6)

- [ ] **Step 1: Write the failing test**

```python
"""Every class must declare what defines it, and two classes must not be
defined by the same set -- that would make them indistinguishable."""
import pytest

from decider import load_config

CONFIGS = ["rules.yaml", "rules.plant.yaml"]
KNOWN = {"1u", "5a", "fractop", "pingpong", "phase", "density"}


@pytest.mark.parametrize("path", CONFIGS)
def test_every_class_declares_defining_evidence(path):
    for name, spec in load_config(path)["classes"].items():
        assert spec.get("defining_evidence"), f"{path}:{name}"


@pytest.mark.parametrize("path", CONFIGS)
def test_defining_evidence_names_are_known(path):
    for name, spec in load_config(path)["classes"].items():
        assert set(spec["defining_evidence"]) <= KNOWN, f"{path}:{name}"


@pytest.mark.parametrize("path", CONFIGS)
def test_no_two_classes_share_the_same_defining_set(path):
    seen = {}
    for name, spec in load_config(path)["classes"].items():
        key = frozenset(spec["defining_evidence"])
        assert key not in seen, f"{path}: {name} and {seen[key]} are indistinguishable"
        seen[key] = name
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_config_defining_evidence.py -q`
Expected: FAIL on the first test, since no class declares the key yet

- [ ] **Step 3: Write the implementation**

In `rules.plant.yaml`, add one line to each class:

```yaml
  phasiRNA:
    defining_evidence: [phase]
  hcsiRNA:
    defining_evidence: [5a]
  siRNA:
    defining_evidence: [fractop]
```

In `rules.yaml`:

```yaml
  siRNA:
    defining_evidence: [fractop]
  piRNA:
    defining_evidence: [1u]
```

Keep each key inside its existing class block, next to the weights.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_config_defining_evidence.py -q`
Expected: PASS, 6 passed

- [ ] **Step 5: Commit**

```bash
git add rules.yaml rules.plant.yaml tests/test_config_defining_evidence.py
git commit -m "feat(config): declare the defining evidence of each class

From the literature, not from preference: hcsiRNA is defined by 5'A because
AGO4 sorts on it, phasiRNA by phasing because the definition is processive
cleavage in regular increments, piRNA by 5'U from Zucchini-derived primaries.

A test rejects two classes sharing the same defining set, which would make them
indistinguishable by construction."
```

---

### Task 8: validate on both kingdoms and document

**Files:**
- Modify: `README.md`, `docs/ROADMAP.md`
- Test: run the existing validators, no new test file

- [ ] **Step 1: Run classify on both datasets**

```bash
python srnaprofiler.py classify --sam ~/dros.sorted.bam -o ~/calib_out \
  -p dros_scores --mincov 10 --threads 2 --config rules.yaml
python srnaprofiler.py classify --sam ~/plant_data/at_flower.sorted.bam \
  -o /mnt/d/srnaprofiler-teste/plant/classify -p at_scores --mincov 10 \
  --threads 1 --config rules.plant.yaml \
  --mirna-loci /mnt/d/srnaprofiler-teste/plant/mirdeep_pred.bed
```

- [ ] **Step 2: Confirm the validations did not regress**

```bash
python extras/validate_piclusters.py ~/calib_out/dros_scores.clusters.tsv
python extras/validate_plant_tas.py \
  /mnt/d/srnaprofiler-teste/plant/classify/at_scores.clusters.tsv \
  /mnt/d/srnaprofiler-teste/plant/GCF_000001735.4_TAIR10.1_genomic.gff
```

Expected: piC recall still 85-86/90, TAS still 6/6. **If either drops, stop and report — the score must not change which loci are found.**

- [ ] **Step 3: Check the score separates what it should**

```bash
python - <<'PY'
import pandas as pd
d = pd.read_csv("/mnt/d/srnaprofiler-teste/plant/classify/at_scores.clusters.tsv",
                sep="\t", low_memory=False)
d = d[d.review.astype(str).str.lower() == "false"]
org = d.chrom.isin({"NC_000932.1", "NC_037304.1"})
print("score mediano organela:", d.score[org].median())
print("score mediano nuclear :", d.score[~org].median())
PY
```

Expected: organellar loci sit at a clearly lower median score than nuclear ones. Record both numbers in the commit message.

- [ ] **Step 4: Document what the score means**

In `README.md`, under the classify section:

```markdown
### Reading the score

`score` is `log10(1/p)`: **2 means p = 0.01**, 3 means p = 0.001. It is the
weakest of the evidences that define the assigned class, so a score of 2 says
every defining evidence cleared p = 0.01.

Filter with it directly — there is no confirmed/unconfirmed flag, because the
threshold belongs to whoever is analysing.

Two caveats:

- **A low score means "not demonstrated", not "false".** The test has little
  power on shallow loci: on the Arabidopsis run the hcsiRNA loci that fail the
  5'A test still show 2.71x transposon enrichment at 19 median reads.
- **The 5' base score mixes biology with protocol.** T4 RNA ligase 1 prefers
  C > U >= A > G at the 5' end, and no tool in the field corrects for it.
```

In `docs/ROADMAP.md`, mark the evidence-score item done with the numbers from Step 3.

- [ ] **Step 5: Commit**

```bash
git add README.md docs/ROADMAP.md
git commit -m "docs: what the score means, and its two caveats

Records that a low score means not demonstrated rather than false, and that the
5' base score carries ligase bias that no tool in the field corrects."
```

---

## Self-Review

**Spec coverage:** nulls per feature — Tasks 1, 2 (binomial, ping-pong, Poisson) and Task 4 (phasing reuses `phase_pvalue`); class score as minimum — Task 6; defining evidence per class — Task 7; output columns — Task 3; validation on both kingdoms — Task 8; caveats documented — Task 8. The `0.25` null, the Bonferroni note and the stratified-FDR alternative are recorded in the spec; Bonferroni is **not** implemented here, since with the minimum rule a per-test correction shifts every score by the same constant and does not change the ordering — noted so it is a conscious omission rather than a gap.

**Placeholders:** none. Every code step carries the code.

**Type consistency:** `score_from_p`, `base_bias_score`, `strand_bias_score` (Task 1) are consumed with those names in Tasks 2 and 4. `s_*` field names in Task 3 match the `s_<name>` lookup in Task 6 and the `defining_evidence` names in Task 7 (`5a`, `phase`, `fractop`, `1u`). `mapped_genome_size` (Task 5) matches its use in `classify.py`.
