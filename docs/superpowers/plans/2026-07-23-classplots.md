# classplots Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `classplots` subcommand that re-aggregates alignments grouped by the class each locus was assigned in `.clusters.tsv`, producing per-class publication signal plots (R + WebLogo) and a matplotlib QC summary.

**Architecture:** A standalone subcommand reads the classifier's `.clusters.tsv` (coords + `cls`) and the indexed BAM. Pure accumulators in `class_signals.py` build per-locus signal matrices from `samio.Aln` (1/NH-weighted, `-F 256`, `-`-strand revcomp'd), pooled per class and summed; the orchestrator in `classplots.py` writes the intermediate tables each plotter expects and shells out to R / WebLogo, plus a pure-Python QC summary in `class_summary.py`.

**Tech Stack:** Python 3 (pysam, pandas, matplotlib), R (ggplot2 / base), external `weblogo` CLI (WebLogo 3). Tests: pytest.

## Global Constraints

- Code, comments, docstrings, logs, and help text in **English**.
- Author is **Lucas**; never credit AI, including in commit messages (no `Co-Authored-By`).
- Reuse the engine's read model via `samio`: primary only (`-F 256`), weight `w = 1/NH`
  (NH via `XM` for bowtie1, `samio._nh`), reverse-complement `-`-strand reads before
  reading first base / building logos. The **only** exception is the WebLogo, which counts
  **distinct** sequences (classic convention).
- Do not modify `classify.py`, `records.py`, `decider.py`, or `rules.yaml`.
- Design doc: `docs/superpowers/specs/2026-07-23-classplots-design.md`.

---

### Task 1: Extract the ping-pong offset spectrum into `signals.py`

Shared 5'-5' offset accumulation, reused by the classifier and by the per-class plots.

**Files:**
- Modify: `signals.py` (add `offset_spectrum`)
- Modify: `features.py:113-124` (call `offset_spectrum` instead of the inline loop)
- Test: `tests/test_signals_offset_spectrum.py` (create)

**Interfaces:**
- Produces: `signals.offset_spectrum(a5: dict[int,float], b5: dict[int,float], lo: int = 1, hi: int = 20) -> collections.Counter` — `Counter[d] = Σ_p a5[p] * b5[p + d - 1]` for `d` in `lo..hi` (overlap length `d`; ping-pong peak at `d = 10`). Same-strand phasing passes `a5 == b5`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_signals_offset_spectrum.py
from collections import Counter
from signals import offset_spectrum


def test_offset_spectrum_peak_at_10():
    plus5 = {100: 2.0}
    minus5 = {109: 3.0}          # 100 + 10 - 1 = 109 -> overlap length 10
    counts = offset_spectrum(plus5, minus5, lo=1, hi=20)
    assert counts[10] == 6.0
    assert sum(counts.values()) == 6.0


def test_offset_spectrum_same_strand_distance():
    a = {10: 1.0, 31: 1.0, 52: 1.0}   # phased ladder, period 21
    counts = offset_spectrum(a, a, lo=1, hi=60)
    assert counts[21] == 2.0          # (10,31) and (31,52)
    assert counts[42] == 1.0          # (10,52)
    assert counts.get(0, 0) == 0      # self-pairs excluded (lo=1)


def test_offset_spectrum_respects_hi():
    counts = offset_spectrum({0: 1.0}, {29: 1.0}, lo=1, hi=20)
    assert counts.get(30, 0) == 0     # overlap length 30 is beyond hi
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_signals_offset_spectrum.py -v`
Expected: FAIL with `ImportError: cannot import name 'offset_spectrum'`

- [ ] **Step 3: Implement `offset_spectrum` in `signals.py`**

Add at the top of `signals.py` (after the existing imports):

```python
from collections import Counter


def offset_spectrum(a5, b5, lo=1, hi=20):
    """5'-5' offset spectrum: Counter[d] = sum over p of a5[p] * b5[p + d - 1].

    `d` is the overlap length in `lo..hi` (ping-pong peak at d=10). Cross-strand
    ping-pong passes plus5/minus5; same-strand phasing passes the same dict twice,
    which yields the 5'-5' distance distribution (peaks at multiples of the phase
    length) and is phase-invariant, so it pools correctly across loci.
    """
    counts = Counter()
    for p, ac in a5.items():
        for d in range(lo, hi + 1):
            bc = b5.get(p + d - 1, 0)
            if bc:
                counts[d] += ac * bc
    return counts
```

- [ ] **Step 4: Rewrite the inline loop in `features.py`**

In `features.py`, replace lines 113-118 (the `overlap_counts = Counter()` loop) with:

```python
    from signals import offset_spectrum
    overlap_counts = offset_spectrum(plus5, minus5, lo=1, hi=20)
```

Add `offset_spectrum` to the existing `from signals import ...` line at the top instead of the local import if preferred; keep behavior identical. Note the previous inline code used `range(0, 20)` with `d + 1` (overlap 1..20) and `minus5.get(a5 + d)` — `offset_spectrum(..., lo=1, hi=20)` reproduces exactly `minus5[p + d - 1]` for d in 1..20.

- [ ] **Step 5: Run the new and existing signal/feature tests**

Run: `pytest tests/test_signals_offset_spectrum.py tests/test_features_pingpong_frac.py tests/test_features_signals.py tests/test_features_core.py -v`
Expected: PASS (existing ping-pong values unchanged)

- [ ] **Step 6: Commit**

```bash
git add signals.py features.py tests/test_signals_offset_spectrum.py
git commit -m "refactor(signals): extract offset_spectrum shared by classifier and plots"
```

---

### Task 2: `class_signals.base_composition` (weighted first-base by size/strand)

**Files:**
- Create: `class_signals.py`
- Test: `tests/test_class_signals_basecomp.py` (create)

**Interfaces:**
- Consumes: `samio.Aln`, `features.revcomp`.
- Produces: `class_signals.base_composition(alns, minlen=15, maxlen=35) -> dict` shaped
  `{size: {'+': {'A':w,'C':w,'G':w,'T':w}, '-': {...}}}` (first base of the RNA-oriented
  sequence, weighted by `1/NH`). Matches `plots.write_base_composition_tsv`'s input.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_class_signals_basecomp.py
from samio import Aln
from class_signals import base_composition


def test_base_composition_weights_and_orients():
    alns = [
        Aln("chr", 0, "+", "TACGT", 1),       # len 5, + -> first base T, w=1
        Aln("chr", 0, "-", "AACGT", 2),        # revcomp = ACGTT -> first base A, w=0.5
    ]
    comp = base_composition(alns, minlen=5, maxlen=5)
    assert comp[5]['+']['T'] == 1.0
    assert comp[5]['-']['A'] == 0.5
    assert comp[5]['+']['A'] == 0.0


def test_base_composition_size_filter():
    alns = [Aln("chr", 0, "+", "AAA", 1)]      # len 3, outside 5..35
    comp = base_composition(alns, minlen=5, maxlen=35)
    assert all(comp[s]['+']['A'] == 0.0 for s in comp)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_class_signals_basecomp.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'class_signals'`

- [ ] **Step 3: Create `class_signals.py` with `base_composition`**

```python
"""Per-class signal accumulators over samio.Aln (1/NH-weighted, RNA 5'->3')."""
from collections import Counter, defaultdict
from features import revcomp

_BASES = ("A", "C", "G", "T")


def _oriented(a):
    """Sequence in the RNA 5'->3' orientation (revcomp on the - strand)."""
    return a.seq if a.strand == "+" else revcomp(a.seq)


def base_composition(alns, minlen=15, maxlen=35):
    """First-base composition by size and strand, weighted by 1/NH."""
    comp = {s: {"+": {b: 0.0 for b in _BASES}, "-": {b: 0.0 for b in _BASES}}
            for s in range(minlen, maxlen + 1)}
    for a in alns:
        seq = _oriented(a)
        L = len(seq)
        if L < minlen or L > maxlen or not seq:
            continue
        first = seq[0]
        if first in comp[L][a.strand]:
            comp[L][a.strand][first] += 1.0 / a.nh
    return comp
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_class_signals_basecomp.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add class_signals.py tests/test_class_signals_basecomp.py
git commit -m "feat(class_signals): weighted first-base composition by size/strand"
```

---

### Task 3: `class_signals.pingpong_offsets` and `phasing_series`

**Files:**
- Modify: `class_signals.py`
- Test: `tests/test_class_signals_biogenesis.py` (create)

**Interfaces:**
- Produces:
  - `class_signals.five_prime_ends(alns) -> tuple[Counter, Counter]` — `(plus5, minus5)`,
    weighted 5' end position of `+` reads (`pos`) and `-` reads (`pos + len - 1`).
  - `class_signals.pingpong_offsets(alns, hi=30) -> Counter` — cross-strand overlap spectrum
    `1..hi` (via `signals.offset_spectrum`), peak at 10.
  - `class_signals.phasing_series(alns, strand, maxdist=100) -> Counter` — same-strand 5'-5'
    **distance** spectrum for one strand (phase-invariant; peaks at multiples of the phase
    length). NOTE: `offset_spectrum` keys are 1-indexed overlap length = distance + 1, so
    `phasing_series` remaps the key by −1 to report true distance (a period-21 ladder peaks
    at 21, aligning with the R reference lines).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_class_signals_biogenesis.py
from samio import Aln
from class_signals import pingpong_offsets, phasing_series


def test_pingpong_peak_at_10():
    # + read 5' at 100 (len 20 -> spans 100..119); - read whose 5' end is at 109
    plus = Aln("chr", 100, "+", "A" * 20, 1)
    minus = Aln("chr", 90, "-", "A" * 20, 1)   # 5' end = 90 + 20 - 1 = 109
    counts = pingpong_offsets([plus, minus], hi=30)
    assert counts[10] == 1.0


def test_phasing_series_peaks_at_period():
    # three + reads with 5' ends 21 nt apart -> phasing period 21
    reads = [Aln("chr", p, "+", "A" * 21, 1) for p in (0, 21, 42)]
    counts = phasing_series(reads, "+", maxdist=100)
    assert counts[21] == 2.0
    assert counts[42] == 1.0


def test_phasing_series_strand_filter():
    reads = [Aln("chr", 0, "-", "A" * 21, 1), Aln("chr", 21, "-", "A" * 21, 1)]
    plus = phasing_series(reads, "+", maxdist=100)
    assert sum(plus.values()) == 0.0     # no + reads
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_class_signals_biogenesis.py -v`
Expected: FAIL with `ImportError: cannot import name 'pingpong_offsets'`

- [ ] **Step 3: Implement in `class_signals.py`**

Append:

```python
from signals import offset_spectrum


def five_prime_ends(alns):
    """Weighted 5' end positions: (plus5, minus5) Counters."""
    plus5, minus5 = Counter(), Counter()
    for a in alns:
        w = 1.0 / a.nh
        if a.strand == "+":
            plus5[a.pos] += w
        else:
            minus5[a.pos + len(a.seq) - 1] += w
    return plus5, minus5


def pingpong_offsets(alns, hi=30):
    """Cross-strand 5'-5' overlap spectrum (ping-pong peak at 10)."""
    plus5, minus5 = five_prime_ends(alns)
    return offset_spectrum(plus5, minus5, lo=1, hi=hi)


def phasing_series(alns, strand, maxdist=100):
    """Same-strand 5'-5' distance spectrum for `strand` (phase-invariant).

    offset_spectrum keys are 1-indexed overlap length (= distance + 1); we remap
    by -1 so the result is keyed by true 5'-5' distance, peaking at multiples of
    the phase length (aligning with the R reference lines).
    """
    plus5, minus5 = five_prime_ends(alns)
    ends = plus5 if strand == "+" else minus5
    raw = offset_spectrum(ends, ends, lo=2, hi=maxdist + 1)
    return Counter({k - 1: v for k, v in raw.items()})
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_class_signals_biogenesis.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add class_signals.py tests/test_class_signals_biogenesis.py
git commit -m "feat(class_signals): ping-pong and phase-invariant phasing spectra"
```

---

### Task 4: `class_signals.pattern_matrix` (per-locus z-scored size×strand rows)

Ported from `clusters.py::_generate_zscore_matrix`, one row per locus, weighted.
The 42-value layout requires exactly 21 sizes; this plot uses the fixed 15–35 window.

**Files:**
- Modify: `class_signals.py`
- Test: `tests/test_class_signals_pattern.py` (create)

**Interfaces:**
- Produces: `class_signals.pattern_row(alns, minlen=15, maxlen=35) -> list[float]` — 42
  z-scored values: indices 0–20 = `+` strand counts for sizes 15..35, 21–41 = `-` strand,
  z-scored across the 42 values (population like the original). Rows are written by the
  orchestrator with a leading locus-name column, `#`-separated, header `Chromosome#0#..#41`,
  matching `plotPattern_publication.R`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_class_signals_pattern.py
import statistics
from samio import Aln
from class_signals import pattern_row


def test_pattern_row_layout_and_zscore():
    # one 21-nt + read (size 21 -> index 6) and one 24-nt - read (size 24 -> +21+9=30)
    alns = [Aln("chr", 0, "+", "A" * 21, 1), Aln("chr", 0, "-", "A" * 24, 2)]
    row = pattern_row(alns, minlen=15, maxlen=35)
    assert len(row) == 42
    # raw counts before z-score: idx 6 = 1.0 (+), idx 30 = 0.5 (-)
    # z-score is monotonic, so the two non-zero cells are the two largest values
    top2 = sorted(range(42), key=lambda i: row[i], reverse=True)[:2]
    assert set(top2) == {6, 30}


def test_pattern_row_all_zero_is_flat():
    row = pattern_row([], minlen=15, maxlen=35)
    assert len(row) == 42
    assert all(v == 0.0 for v in row)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_class_signals_pattern.py -v`
Expected: FAIL with `ImportError: cannot import name 'pattern_row'`

- [ ] **Step 3: Implement `pattern_row` in `class_signals.py`**

```python
import statistics


def pattern_row(alns, minlen=15, maxlen=35):
    """42 z-scored size x strand counts (0-20 = +, 21-41 = -), weighted by 1/NH."""
    span = maxlen - minlen + 1
    values = [0.0] * (span * 2)
    for a in alns:
        L = len(_oriented(a))
        if L < minlen or L > maxlen:
            continue
        idx = L - minlen if a.strand == "+" else (L - minlen) + span
        values[idx] += 1.0 / a.nh
    nonzero = [v for v in values if v != 0.0]
    if len(nonzero) < 2:
        return [0.0] * len(values)
    mean = statistics.fmean(values)
    sd = statistics.pstdev(values)
    if sd == 0:
        return [0.0] * len(values)
    return [(v - mean) / sd for v in values]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_class_signals_pattern.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add class_signals.py tests/test_class_signals_pattern.py
git commit -m "feat(class_signals): per-locus z-scored size x strand pattern row"
```

---

### Task 5: `class_signals.coverage_profile` (per-strand normalized metaplot)

**Files:**
- Modify: `class_signals.py`
- Test: `tests/test_class_signals_coverage.py` (create)

**Interfaces:**
- Produces: `class_signals.coverage_profile(alns, region_start, region_end, bins=100) -> dict`
  shaped `{'+': [float]*bins, '-': [float]*bins}` — weighted coverage over the locus span,
  each locus rescaled to `[0,1)` into `bins` bins (so loci of different lengths are
  comparable). The orchestrator averages these across loci of a class.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_class_signals_coverage.py
from samio import Aln
from class_signals import coverage_profile


def test_coverage_profile_bins_and_strand():
    # a + read covering the first half of a length-100 region
    alns = [Aln("chr", 0, "+", "A" * 50, 1)]
    prof = coverage_profile(alns, region_start=0, region_end=100, bins=10)
    assert len(prof["+"]) == 10
    assert prof["+"][0] > 0.0        # first bin covered
    assert prof["+"][9] == 0.0       # last bin not covered
    assert all(v == 0.0 for v in prof["-"])


def test_coverage_profile_weighting():
    alns = [Aln("chr", 0, "+", "A" * 100, 2)]   # w = 0.5 over whole region
    prof = coverage_profile(alns, region_start=0, region_end=100, bins=10)
    assert all(abs(v - 0.5) < 1e-9 for v in prof["+"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_class_signals_coverage.py -v`
Expected: FAIL with `ImportError: cannot import name 'coverage_profile'`

- [ ] **Step 3: Implement `coverage_profile` in `class_signals.py`**

```python
def coverage_profile(alns, region_start, region_end, bins=100):
    """Per-strand weighted coverage over the locus span, rescaled to `bins` bins."""
    length = region_end - region_start
    prof = {"+": [0.0] * bins, "-": [0.0] * bins}
    if length <= 0:
        return prof
    for a in alns:
        w = 1.0 / a.nh
        lane = prof[a.strand]
        read_start = max(a.pos, region_start)
        read_end = min(a.pos + len(a.seq), region_end)
        for p in range(read_start, read_end):
            b = int((p - region_start) / length * bins)
            if b >= bins:
                b = bins - 1
            lane[b] += w
    return prof
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_class_signals_coverage.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add class_signals.py tests/test_class_signals_coverage.py
git commit -m "feat(class_signals): per-strand normalized coverage profile"
```

---

### Task 6: `class_signals.strand_fasta` (distinct sequences for WebLogo)

**Files:**
- Modify: `class_signals.py`
- Test: `tests/test_class_signals_fasta.py` (create)

**Interfaces:**
- Produces: `class_signals.strand_sequences(alns, strand, size_lo, size_hi, trim) -> list[str]`
  — sorted list of **distinct** RNA-oriented sequences on `strand`, size-filtered to
  `[size_lo, size_hi]`, trimmed to the first `trim` nt from the RNA 5' end, `T→U`.
  `-` strand reads are reverse-complemented first (SAM stores genome-forward SEQ), so
  position 1 is the RNA 5' end for both strands.
- Produces: `class_signals.write_fasta(sequences, path) -> None` — writes `>seq_N` records.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_class_signals_fasta.py
from samio import Aln
from class_signals import strand_sequences


def test_strand_sequences_orient_trim_uniq():
    alns = [
        Aln("chr", 0, "+", "TACGTA", 1),        # + -> TACGTA, trim4 -> UACG
        Aln("chr", 0, "+", "TACGTT", 1),        # + -> TACGTT, trim4 -> UACG (dup)
        Aln("chr", 0, "-", "AACGTA", 1),         # revcomp = TACGTT, trim4 -> UACG (dup)
    ]
    seqs = strand_sequences(alns, "+", size_lo=6, size_hi=6, trim=4)
    assert seqs == ["UACG"]                      # collapsed to one distinct sequence


def test_strand_sequences_minus_uses_revcomp():
    alns = [Aln("chr", 0, "-", "AAATTT", 1)]     # revcomp = AAATTT -> trim3 -> AAA
    seqs = strand_sequences(alns, "-", size_lo=6, size_hi=6, trim=3)
    assert seqs == ["AAA"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_class_signals_fasta.py -v`
Expected: FAIL with `ImportError: cannot import name 'strand_sequences'`

- [ ] **Step 3: Implement in `class_signals.py`**

```python
def strand_sequences(alns, strand, size_lo, size_hi, trim):
    """Distinct RNA-oriented, size-filtered, trimmed, T->U sequences on `strand`."""
    seen = set()
    for a in alns:
        if a.strand != strand:
            continue
        seq = _oriented(a)
        L = len(seq)
        if L < size_lo or L > size_hi:
            continue
        trimmed = seq[:trim].replace("T", "U").replace("t", "u")
        if len(trimmed) == trim:
            seen.add(trimmed)
    return sorted(seen)


def write_fasta(sequences, path):
    with open(path, "w") as fh:
        for i, seq in enumerate(sequences, 1):
            fh.write(f">seq_{i}\n{seq}\n")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_class_signals_fasta.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add class_signals.py tests/test_class_signals_fasta.py
git commit -m "feat(class_signals): distinct-sequence FASTA prep for WebLogo"
```

---

### Task 7: `classplots.load_clusters` + `Locus`

**Files:**
- Create: `classplots.py`
- Test: `tests/test_classplots_load.py` (create)
- Test data: `tests/data/mini.clusters.tsv` (create)

**Interfaces:**
- Produces: `classplots.Locus` dataclass with `chrom:str, start:int, end:int, name:str,
  cls:str, reads:float, dicer_call:str`.
- Produces: `classplots.load_clusters(tsv_path) -> list[Locus]` — parses the reporter TSV
  (header from `reporter._COLS`), skipping rows with a blank `cls`. Raises `ValueError` if a
  required column is missing.

- [ ] **Step 1: Create the fixture TSV**

`tests/data/mini.clusters.tsv` (tab-separated; only the columns used here need realism, but include the full header so parsing is realistic):

```
chrom	start	end	name	reads	reads_rpm	distinct_seqs	major_rna	dicer_call	frac_in_window	frac_top	strand_call	pct_1u	pct_10a	gc	pingpong_z	pingpong_pairs10	pingpong_frac	overlap19_frac	phase_score	phase_pvalue	phase_number	phase_ratio	multimap_frac	cls	confidence	review	evidence
chr1	100	400	chr1:100-400	500	1000	40	ACGT	24	0.9	0.95	+	0.8	0.2	0.45	3.1	50	0.4	0.05	2.0	0.01	4	0.6	0.9	piRNA	0.8	False	sz;1u
chr1	800	1000	chr1:800-1000	300	600	25	TTGC	21	0.95	0.5	.	0.3	0.1	0.5	0.2	2	0.0	0.02	5.0	0.001	6	0.7	0.2	siRNA	0.75	False	sz;fractop
chr2	50	200	chr2:50-200	20	40	5	AAAA	N	0.4	0.9	+	0.5	0.1	0.6	0.0	0	0.0	0.0	0.0	1.0	0	0.0	0.5		0.0	True	lowcov
```

- [ ] **Step 2: Write the failing test**

```python
# tests/test_classplots_load.py
import pytest
from classplots import load_clusters, Locus

TSV = "tests/data/mini.clusters.tsv"


def test_load_clusters_parses_and_skips_blank_cls():
    loci = load_clusters(TSV)
    assert len(loci) == 2                       # blank-cls row skipped
    assert isinstance(loci[0], Locus)
    piRNA = [l for l in loci if l.cls == "piRNA"][0]
    assert piRNA.chrom == "chr1" and piRNA.start == 100 and piRNA.end == 400
    assert piRNA.reads == 500.0 and piRNA.dicer_call == "24"


def test_load_clusters_missing_column_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("chrom\tstart\tend\n1\t0\t9\n")
    with pytest.raises(ValueError):
        load_clusters(str(bad))
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest tests/test_classplots_load.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'classplots'`

- [ ] **Step 4: Implement `classplots.py` with `Locus` and `load_clusters`**

```python
"""Per-class plotting: read the classifier TSV + BAM, plot signals grouped by class."""
from dataclasses import dataclass

_REQUIRED = ("chrom", "start", "end", "name", "reads", "dicer_call", "cls")


@dataclass
class Locus:
    chrom: str
    start: int
    end: int
    name: str
    cls: str
    reads: float
    dicer_call: str


def load_clusters(tsv_path):
    """Parse the reporter TSV into Locus rows, skipping rows with a blank `cls`."""
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in _REQUIRED if c not in idx]
        if missing:
            raise ValueError(f"clusters TSV missing columns: {missing}")
        loci = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cls = f[idx["cls"]]
            if not cls:
                continue
            loci.append(Locus(
                chrom=f[idx["chrom"]], start=int(f[idx["start"]]),
                end=int(f[idx["end"]]), name=f[idx["name"]],
                cls=cls, reads=float(f[idx["reads"]]),
                dicer_call=f[idx["dicer_call"]]))
    return loci
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_classplots_load.py -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add classplots.py tests/test_classplots_load.py tests/data/mini.clusters.tsv
git commit -m "feat(classplots): Locus dataclass and clusters TSV loader"
```

---

### Task 8: `classplots.group_by_class`, `top_n_loci`, `dominant_dicer`

**Files:**
- Modify: `classplots.py`
- Test: `tests/test_classplots_group.py` (create)

**Interfaces:**
- Produces:
  - `classplots.group_by_class(loci, classes=None) -> dict[str, list[Locus]]` — groups by
    `cls`; if `classes` is given, restricts to that set.
  - `classplots.top_n_loci(loci, n) -> list[Locus]` — the `n` loci with the largest `reads`
    (descending); `n <= 0` returns `[]`.
  - `classplots.dominant_dicer(loci, default=21) -> int` — the most common integer
    `dicer_call` among the loci (ignoring `"N"`); `default` if none.
  - `classplots.parse_region(spec, name=None) -> Locus` — parses `"chrom:start-end"` into a
    `Locus` (0-based half-open `start`/`end`, `cls=""`, `reads=0.0`, `dicer_call="N"`);
    raises `ValueError` on a malformed spec.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_classplots_group.py
from classplots import Locus, group_by_class, top_n_loci, dominant_dicer


def _l(name, cls, reads, dicer="24"):
    return Locus("chr", 0, 9, name, cls, reads, dicer)


def test_group_by_class_and_restrict():
    loci = [_l("a", "piRNA", 10), _l("b", "siRNA", 5), _l("c", "piRNA", 20)]
    grouped = group_by_class(loci)
    assert set(grouped) == {"piRNA", "siRNA"}
    assert len(grouped["piRNA"]) == 2
    only = group_by_class(loci, classes={"piRNA"})
    assert set(only) == {"piRNA"}


def test_top_n_loci_orders_by_reads():
    loci = [_l("a", "piRNA", 10), _l("c", "piRNA", 20), _l("b", "piRNA", 5)]
    top = top_n_loci(loci, 2)
    assert [l.name for l in top] == ["c", "a"]
    assert top_n_loci(loci, 0) == []


def test_dominant_dicer_ignores_N():
    loci = [_l("a", "piRNA", 1, "24"), _l("b", "piRNA", 1, "24"), _l("c", "piRNA", 1, "N")]
    assert dominant_dicer(loci) == 24
    assert dominant_dicer([_l("a", "piRNA", 1, "N")], default=21) == 21


def test_parse_region():
    from classplots import parse_region
    import pytest
    loc = parse_region("chr1:100-400")
    assert (loc.chrom, loc.start, loc.end) == ("chr1", 100, 400)
    assert loc.name == "chr1:100-400"
    with pytest.raises(ValueError):
        parse_region("chr1-100-400")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_classplots_group.py -v`
Expected: FAIL with `ImportError: cannot import name 'group_by_class'`

- [ ] **Step 3: Implement in `classplots.py`**

```python
from collections import Counter, defaultdict


def group_by_class(loci, classes=None):
    grouped = defaultdict(list)
    for l in loci:
        if classes is None or l.cls in classes:
            grouped[l.cls].append(l)
    return dict(grouped)


def top_n_loci(loci, n):
    if n <= 0:
        return []
    return sorted(loci, key=lambda l: l.reads, reverse=True)[:n]


def dominant_dicer(loci, default=21):
    calls = Counter(int(l.dicer_call) for l in loci if l.dicer_call not in ("N", ""))
    return calls.most_common(1)[0][0] if calls else default


def parse_region(spec, name=None):
    """Parse 'chrom:start-end' (0-based half-open) into a Locus."""
    try:
        chrom, coords = spec.rsplit(":", 1)
        start_s, end_s = coords.split("-")
        start, end = int(start_s), int(end_s)
    except ValueError:
        raise ValueError(f"malformed region (expected chrom:start-end): {spec!r}")
    return Locus(chrom=chrom, start=start, end=end, name=name or spec,
                 cls="", reads=0.0, dicer_call="N")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_classplots_group.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add classplots.py tests/test_classplots_group.py
git commit -m "feat(classplots): grouping, top-N, and dominant-DicerCall helpers"
```

---

### Task 9: New R scripts (ping-pong, phasing, coverage) and retire the old one

**Files:**
- Create: `plotPingPong.R`, `plotPhasing.R`, `plotCoverageProfile.R`
- Delete: `plotOffsetPhasingSeries.R`
- Test: `tests/test_r_scripts_smoke.py` (create)

**Interfaces:**
- `plotPingPong.R <tsv> <out_prefix>` — `tsv` is `overlap<TAB>pairs` (header); bar chart over
  overlap 1–30, the 10-nt bar highlighted.
- `plotPhasing.R <plus_tsv> <minus_tsv> <period> <out_prefix>` — each tsv `distance<TAB>count`
  (header); two stacked line panels (+ / −) with dashed reference lines at multiples of `period`.
- `plotCoverageProfile.R <tsv> <out_prefix>` — `tsv` is `bin<TAB>plus<TAB>minus` (header);
  two lines over normalized position 0–1.

- [ ] **Step 1: Write `plotPingPong.R`**

```r
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
d <- read.table(args[1], header = TRUE, sep = "\t")
d$fill <- ifelse(d$overlap == 10, "ping-pong (10)", "other")
pdf(paste0(args[2], ".pingpong.pdf"), width = 6, height = 4)
ggplot(d, aes(x = overlap, y = pairs, fill = fill)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("ping-pong (10)" = "#e41a1c", "other" = "grey60")) +
  theme_classic(base_size = 16) +
  xlab("5'-5' overlap (nt)") + ylab("Weighted pairs") +
  theme(legend.title = element_blank())
dev.off()
```

- [ ] **Step 2: Write `plotPhasing.R`**

```r
args <- commandArgs(trailingOnly = TRUE)
plus_file <- args[1]; minus_file <- args[2]
period <- as.numeric(args[3]); out <- args[4]
panel <- function(file, title, period) {
  if (file.exists(file) && file.info(file)$size > 0) {
    d <- read.table(file, header = TRUE, sep = "\t")
    d <- d[order(d$distance), ]
    plot(d$distance, d$count, type = "l", lwd = 1, col = "black",
         main = title, xlab = "5'-5' distance (nt)", ylab = "Weighted pairs",
         panel.first = grid(lty = "dotted", col = "lightgray"))
    for (m in seq(period, max(d$distance, na.rm = TRUE), by = period))
      abline(v = m, lty = 2, col = "red")
  } else {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
    text(1, 1, "No data")
  }
}
pdf(paste0(out, ".phasing.pdf"), width = 5, height = 6)
par(mfrow = c(2, 1), mar = c(4.5, 4.5, 3, 2))
panel(plus_file, "Sense (+)", period)
panel(minus_file, "Antisense (-)", period)
dev.off()
```

- [ ] **Step 3: Write `plotCoverageProfile.R`**

```r
args <- commandArgs(trailingOnly = TRUE)
d <- read.table(args[1], header = TRUE, sep = "\t")
x <- d$bin / max(d$bin)
pdf(paste0(args[2], ".coverage.pdf"), width = 6, height = 4)
plot(x, d$plus, type = "l", col = "#377eb8", lwd = 2,
     ylim = range(c(d$plus, d$minus), 0), xlab = "Normalized position",
     ylab = "Mean weighted coverage", main = "Coverage profile")
lines(x, d$minus, col = "#e41a1c", lwd = 2)
legend("topright", legend = c("+", "-"), col = c("#377eb8", "#e41a1c"), lwd = 2, bty = "n")
dev.off()
```

- [ ] **Step 4: Delete the retired script**

```bash
git rm plotOffsetPhasingSeries.R
```

- [ ] **Step 5: Write a smoke test that skips without Rscript**

```python
# tests/test_r_scripts_smoke.py
import shutil, subprocess
from pathlib import Path
import pytest

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript absent")
ROOT = Path(__file__).resolve().parent.parent


def test_plot_pingpong_runs(tmp_path):
    tsv = tmp_path / "pp.tsv"
    lines = ["overlap\tpairs"] + [f"{o}\t{5 if o != 10 else 90}" for o in range(1, 31)]
    tsv.write_text("\n".join(lines) + "\n")
    out = tmp_path / "x"
    r = subprocess.run(["Rscript", str(ROOT / "plotPingPong.R"), str(tsv), str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "x.pingpong.pdf").exists()
```

- [ ] **Step 6: Run the smoke test**

Run: `pytest tests/test_r_scripts_smoke.py -v`
Expected: PASS (or SKIP if `Rscript` is not installed)

- [ ] **Step 7: Commit**

```bash
git add plotPingPong.R plotPhasing.R plotCoverageProfile.R tests/test_r_scripts_smoke.py
git commit -m "feat(plots): correct ping-pong/phasing/coverage R scripts; retire offset script"
```

---

### Task 10: WebLogo invocation wrapper

**Files:**
- Modify: `classplots.py`
- Test: `tests/test_classplots_weblogo.py` (create)

**Interfaces:**
- Produces: `classplots.run_weblogo(fasta_path, out_pdf, title) -> bool` — runs the external
  `weblogo` CLI (WebLogo 3) with `-F pdf -A rna -U bits -s large --resolution 300
  --errorbars NO -c classic -P <title>`; returns `False` (with a warning) if the `weblogo`
  binary is absent from PATH or the FASTA is empty.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_classplots_weblogo.py
import shutil
from classplots import run_weblogo


def test_run_weblogo_absent_binary_returns_false(tmp_path, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    fasta = tmp_path / "s.fa"
    fasta.write_text(">seq_1\nUACG\n")
    assert run_weblogo(str(fasta), str(tmp_path / "o.pdf"), "sense") is False


def test_run_weblogo_empty_fasta_returns_false(tmp_path):
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")
    assert run_weblogo(str(fasta), str(tmp_path / "o.pdf"), "sense") is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_classplots_weblogo.py -v`
Expected: FAIL with `ImportError: cannot import name 'run_weblogo'`

- [ ] **Step 3: Implement `run_weblogo` in `classplots.py`**

```python
import shutil
import subprocess
import os
from utils import log


def run_weblogo(fasta_path, out_pdf, title):
    """Run the external WebLogo 3 CLI; skip (return False) if unavailable/empty."""
    if shutil.which("weblogo") is None:
        log("weblogo CLI not found in PATH; skipping logo", "WARN")
        return False
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        log(f"empty FASTA, skipping logo: {fasta_path}", "WARN")
        return False
    cmd = ["weblogo", "-F", "pdf", "-A", "rna", "-U", "bits", "-s", "large",
           "--resolution", "300", "--errorbars", "NO", "-c", "classic",
           "-P", title, "-f", fasta_path, "-o", out_pdf]
    r = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if r.returncode != 0:
        log(f"weblogo failed: {r.stderr}", "ERROR")
        return False
    return True
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_classplots_weblogo.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add classplots.py tests/test_classplots_weblogo.py
git commit -m "feat(classplots): WebLogo CLI wrapper with graceful skip"
```

---

### Task 11: `class_summary` QC panels (matplotlib)

**Files:**
- Create: `class_summary.py`
- Test: `tests/test_class_summary.py` (create)

**Interfaces:**
- Consumes: the reporter `.clusters.tsv` via pandas.
- Produces: `class_summary.write_summary(tsv_path, out_dir, prefix) -> list[str]` — writes
  PDF panels (`{prefix}.summary.counts.pdf`, `.dicer.pdf`, `.signals.pdf`, `.review.pdf`)
  and returns the list of written paths. Uses the non-interactive `Agg` backend. Skips a
  panel if its column is entirely missing but never raises on empty groups.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_class_summary.py
from pathlib import Path
from class_summary import write_summary

TSV = "tests/data/mini.clusters.tsv"


def test_write_summary_creates_panels(tmp_path):
    paths = write_summary(TSV, str(tmp_path), "mini")
    assert paths, "expected at least one panel"
    for p in paths:
        assert Path(p).exists() and Path(p).stat().st_size > 0
    names = {Path(p).name for p in paths}
    assert "mini.summary.counts.pdf" in names
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_class_summary.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'class_summary'`

- [ ] **Step 3: Implement `class_summary.py`**

```python
"""Matplotlib QC summary panels from the classifier .clusters.tsv (no R needed)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def write_summary(tsv_path, out_dir, prefix):
    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_csv(tsv_path, sep="\t")
    df = df[df["cls"].notna() & (df["cls"].astype(str) != "")]
    written = []

    def save(fig, tag):
        path = os.path.join(out_dir, f"{prefix}.summary.{tag}.pdf")
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    # counts per class
    fig, ax = plt.subplots(figsize=(5, 4))
    df["cls"].value_counts().plot.bar(ax=ax, color="grey")
    ax.set_ylabel("loci"); ax.set_title("Loci per class")
    save(fig, "counts")

    # dicer_call distribution per class
    if "dicer_call" in df:
        fig, ax = plt.subplots(figsize=(6, 4))
        (df.groupby("cls")["dicer_call"].value_counts()
           .unstack(fill_value=0).T.plot.bar(ax=ax))
        ax.set_xlabel("DicerCall"); ax.set_ylabel("loci")
        save(fig, "dicer")

    # scalar signals per class (box)
    signals = [c for c in ("pingpong_frac", "phase_score", "pct_1u", "gc", "multimap_frac")
               if c in df]
    if signals:
        fig, axes = plt.subplots(1, len(signals), figsize=(3 * len(signals), 4))
        if len(signals) == 1:
            axes = [axes]
        for ax, col in zip(axes, signals):
            groups = [g[col].dropna().values for _, g in df.groupby("cls")]
            ax.boxplot(groups, labels=list(df.groupby("cls").groups))
            ax.set_title(col)
        save(fig, "signals")

    # confidence + review rate per class
    if {"confidence", "review"} <= set(df.columns):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4))
        for cls, g in df.groupby("cls"):
            a1.hist(g["confidence"].dropna(), bins=10, alpha=0.5, label=cls)
        a1.set_title("Confidence"); a1.legend()
        rate = df.groupby("cls")["review"].apply(
            lambda s: (s.astype(str) == "True").mean())
        rate.plot.bar(ax=a2, color="salmon"); a2.set_title("Review rate")
        save(fig, "review")

    return written
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_class_summary.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add class_summary.py tests/test_class_summary.py
git commit -m "feat(class_summary): matplotlib QC panels from clusters TSV"
```

---

### Task 12: `classplots.run_classplots` orchestrator

Ties everything together: read loci, group, per-class pooled signals + top-N, write tables,
invoke R / WebLogo, and the QC summary. Gated by `--no-signals` / `--no-summary`.

**Files:**
- Modify: `classplots.py`
- Test: `tests/test_classplots_run.py` (create)

**Interfaces:**
- Consumes: `samio.read_region`, `samio.prefers_xm_as_nh`, all `class_signals.*`,
  `plots.generate_base_distribution_plot`, `plots.generate_pattern_plot`,
  `plots.run_r_script`, `plots.write_base_composition_tsv`, `class_summary.write_summary`,
  `run_weblogo`.
- Produces: `classplots.run_classplots(args)` where `args` has `clusters, bam, outdir,
  prefix, top_n, classes, region, min_length, max_length, logo_trim, no_signals, no_summary`.
  `args.region` is an optional list of `chrom:start-end` strings (accessed via
  `getattr(args, "region", None)`), each plotted as its own single-locus unit. Writes into
  `outdir/tables/` and `outdir/plots/`; returns `None`.

- [ ] **Step 1: Write the failing test (summary-only path, no R/weblogo needed)**

```python
# tests/test_classplots_run.py
import types
from pathlib import Path
import pysam
from classplots import run_classplots


def _make_bam(path):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 2000}],
              "PG": [{"ID": "bowtie", "PN": "bowtie"}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for pos in range(100, 380, 21):
            a = pysam.AlignedSegment()
            a.query_name = f"r{pos}"; a.query_sequence = "T" + "A" * 23
            a.flag = 0; a.reference_id = 0; a.reference_start = pos
            a.mapping_quality = 255
            a.cigar = [(0, 24)]
            a.set_tag("XM", 1)
            bam.write(a)
    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


def _args(**kw):
    base = dict(clusters="tests/data/mini.clusters.tsv", bam=None, outdir=None,
               prefix="mini", top_n=1, classes=None, min_length=15, max_length=35,
               logo_trim=None, no_signals=True, no_summary=False)
    base.update(kw)
    return types.SimpleNamespace(**base)


def test_run_classplots_summary_only(tmp_path):
    bam = tmp_path / "reads.bam"; _make_bam(bam)
    out = tmp_path / "out"
    run_classplots(_args(bam=str(bam), outdir=str(out)))
    pdfs = list((out / "plots").glob("mini.summary.*.pdf"))
    assert pdfs, "summary panels should be written even with --no-signals"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_classplots_run.py -v`
Expected: FAIL with `ImportError: cannot import name 'run_classplots'`

- [ ] **Step 3: Implement `run_classplots` in `classplots.py`**

```python
import plots
import class_signals as cs
import class_summary
from samio import read_region, prefers_xm_as_nh


def _fetch(bam, locus, xm):
    return list(read_region(bam, locus.chrom, locus.start, locus.end, xm_as_nh=xm))


def _write_pingpong_tsv(counts, path, hi=30):
    with open(path, "w") as fh:
        fh.write("overlap\tpairs\n")
        for o in range(1, hi + 1):
            fh.write(f"{o}\t{counts.get(o, 0.0)}\n")


def _write_distance_tsv(counts, path, maxdist=100):
    with open(path, "w") as fh:
        fh.write("distance\tcount\n")
        for d in range(1, maxdist + 1):
            fh.write(f"{d}\t{counts.get(d, 0.0)}\n")


def _write_coverage_tsv(profile, path):
    with open(path, "w") as fh:
        fh.write("bin\tplus\tminus\n")
        for i, (p, m) in enumerate(zip(profile["+"], profile["-"])):
            fh.write(f"{i}\t{p}\t{m}\n")


def _write_pattern_matrix(rows, path):
    with open(path, "w") as fh:
        fh.write("Chromosome#" + "#".join(str(i) for i in range(42)) + "\n")
        for name, row in rows:
            fh.write(name + "#" + "#".join(f"{z:.6f}" for z in row) + "\n")


def _signal_plots_for(name, alns_list, period, size_lo, size_hi, trim,
                      tables, plotsdir, prefix, rpath):
    """Write tables + render all signal plots for one pooled unit (class or locus)."""
    from collections import Counter
    all_alns = [a for alns in alns_list for a in alns]

    # size distribution (reuse the existing base-composition writer + R)
    comp = cs.base_composition(all_alns)
    sd_tsv = tables / f"{prefix}.{name}.sizedist.tsv"
    plots.write_base_composition_tsv(comp, sd_tsv)
    plots.generate_base_distribution_plot(
        sd_tsv, plotsdir / f"{prefix}.{name}.sizedist", rpath, "per_base")

    # pattern matrix (one row per unit member)
    rows = [(f"{name}_{i}", cs.pattern_row(alns)) for i, alns in enumerate(alns_list)]
    pm = tables / f"{prefix}.{name}.pattern.matrix"
    _write_pattern_matrix(rows, pm)
    plots.generate_pattern_plot(pm, plotsdir / f"{prefix}.{name}.pattern", rpath)

    # ping-pong
    pp = Counter()
    for alns in alns_list:
        pp.update(cs.pingpong_offsets(alns, hi=30))
    pp_tsv = tables / f"{prefix}.{name}.pingpong.tsv"
    _write_pingpong_tsv(pp, pp_tsv)
    plots.run_r_script("plotPingPong.R", [str(pp_tsv), str(plotsdir / f"{prefix}.{name}")],
                       rpath)

    # phasing (per strand)
    for strand, tag in (("+", "plus"), ("-", "minus")):
        ph = Counter()
        for alns in alns_list:
            ph.update(cs.phasing_series(alns, strand, maxdist=100))
        _write_distance_tsv(ph, tables / f"{prefix}.{name}.phasing.{tag}.tsv")
    plots.run_r_script("plotPhasing.R", [
        str(tables / f"{prefix}.{name}.phasing.plus.tsv"),
        str(tables / f"{prefix}.{name}.phasing.minus.tsv"),
        str(period), str(plotsdir / f"{prefix}.{name}")], rpath)

    # coverage profile (average normalized coverage across members)
    acc = {"+": [0.0] * 100, "-": [0.0] * 100}
    n = 0
    for alns, region in alns_list_with_regions(alns_list):
        pass  # replaced below
```

Note: the coverage step needs each member's region; refactor `_signal_plots_for` to accept
`members` as `(alns, region_start, region_end)` tuples. Use this corrected version instead of
the placeholder coverage block above:

```python
def _signal_plots_for(name, members, period, size_lo, size_hi, trim,
                      tables, plotsdir, prefix, rpath):
    """members: list of (alns, region_start, region_end)."""
    from collections import Counter
    alns_list = [m[0] for m in members]
    all_alns = [a for alns in alns_list for a in alns]

    comp = cs.base_composition(all_alns)
    sd_tsv = tables / f"{prefix}.{name}.sizedist.tsv"
    plots.write_base_composition_tsv(comp, sd_tsv)
    plots.generate_base_distribution_plot(
        sd_tsv, plotsdir / f"{prefix}.{name}.sizedist", rpath, "per_base")

    rows = [(f"{name}_{i}", cs.pattern_row(alns)) for i, alns in enumerate(alns_list)]
    _write_pattern_matrix(rows, tables / f"{prefix}.{name}.pattern.matrix")
    plots.generate_pattern_plot(
        tables / f"{prefix}.{name}.pattern.matrix",
        plotsdir / f"{prefix}.{name}.pattern", rpath)

    pp = Counter()
    for alns in alns_list:
        pp.update(cs.pingpong_offsets(alns, hi=30))
    _write_pingpong_tsv(pp, tables / f"{prefix}.{name}.pingpong.tsv")
    plots.run_r_script("plotPingPong.R",
        [str(tables / f"{prefix}.{name}.pingpong.tsv"), str(plotsdir / f"{prefix}.{name}")],
        rpath)

    for strand, tag in (("+", "plus"), ("-", "minus")):
        ph = Counter()
        for alns in alns_list:
            ph.update(cs.phasing_series(alns, strand, maxdist=100))
        _write_distance_tsv(ph, tables / f"{prefix}.{name}.phasing.{tag}.tsv")
    plots.run_r_script("plotPhasing.R", [
        str(tables / f"{prefix}.{name}.phasing.plus.tsv"),
        str(tables / f"{prefix}.{name}.phasing.minus.tsv"),
        str(period), str(plotsdir / f"{prefix}.{name}")], rpath)

    acc = {"+": [0.0] * 100, "-": [0.0] * 100}
    for alns, rs, re in members:
        prof = cs.coverage_profile(alns, rs, re, bins=100)
        for s in ("+", "-"):
            for i in range(100):
                acc[s][i] += prof[s][i]
    if members:
        for s in ("+", "-"):
            acc[s] = [v / len(members) for v in acc[s]]
    _write_coverage_tsv(acc, tables / f"{prefix}.{name}.coverage.tsv")
    plots.run_r_script("plotCoverageProfile.R",
        [str(tables / f"{prefix}.{name}.coverage.tsv"), str(plotsdir / f"{prefix}.{name}")],
        rpath)

    # logos (sense + antisense) via WebLogo
    for strand, title in (("+", "sense"), ("-", "antisense")):
        seqs = cs.strand_sequences(all_alns, strand, size_lo, size_hi, trim)
        fa = tables / f"{prefix}.{name}.logo.{title}.fasta"
        cs.write_fasta(seqs, fa)
        run_weblogo(str(fa), str(plotsdir / f"{prefix}.{name}.logo.{title}.pdf"), title)


def run_classplots(args):
    from pathlib import Path
    out = Path(args.outdir)
    tables = out / "tables"; plotsdir = out / "plots"
    tables.mkdir(parents=True, exist_ok=True)
    plotsdir.mkdir(parents=True, exist_ok=True)

    if not args.no_summary:
        class_summary.write_summary(args.clusters, str(plotsdir), args.prefix)

    if args.no_signals:
        return

    loci = load_clusters(args.clusters)
    classes = set(args.classes.split(",")) if args.classes else None
    grouped = group_by_class(loci, classes)
    xm = prefers_xm_as_nh(args.bam)
    rpath = plots.get_r_scripts_path()

    for cls, cls_loci in grouped.items():
        period = dominant_dicer(cls_loci)
        trim = args.logo_trim or period
        members = [(_fetch(args.bam, l, xm), l.start, l.end) for l in cls_loci]
        _signal_plots_for(cls, members, period, args.min_length, args.max_length,
                          trim, tables, plotsdir, args.prefix, rpath)
        for l in top_n_loci(cls_loci, args.top_n):
            m = [(_fetch(args.bam, l, xm), l.start, l.end)]
            _signal_plots_for(f"{cls}.locus.{l.name}", m, period, args.min_length,
                              args.max_length, trim, tables, plotsdir, args.prefix, rpath)

    # explicit user-requested coordinates (regardless of class)
    default_period = dominant_dicer(loci)
    for spec in (getattr(args, "region", None) or []):
        reg = parse_region(spec)
        trim = args.logo_trim or default_period
        m = [(_fetch(args.bam, reg, xm), reg.start, reg.end)]
        safe = reg.name.replace(":", "_").replace("-", "_")
        _signal_plots_for(f"region.{safe}", m, default_period, args.min_length,
                          args.max_length, trim, tables, plotsdir, args.prefix, rpath)
```

Remove the first (placeholder) `_signal_plots_for` definition entirely; keep only the
corrected `members`-based version.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_classplots_run.py -v`
Expected: PASS

- [ ] **Step 5: Run the whole suite**

Run: `pytest -q`
Expected: PASS (existing 37 + new tests)

- [ ] **Step 6: Commit**

```bash
git add classplots.py tests/test_classplots_run.py
git commit -m "feat(classplots): run_classplots orchestrator (per-class + top-N + summary)"
```

---

### Task 13: Wire the `classplots` subcommand into the CLI

**Files:**
- Modify: `srnaprofiler.py:595-608` region (add a new subparser after `classify`)
- Test: `tests/test_cli_classplots.py` (create)

**Interfaces:**
- Consumes: `classplots.run_classplots`.
- Produces: `srnaprofiler classplots ...` subcommand whose parsed args match the
  `run_classplots(args)` attribute names in Task 12.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cli_classplots.py
from srnaprofiler import build_parser


def test_classplots_subcommand_parses():
    parser = build_parser()
    args = parser.parse_args([
        "classplots", "--clusters", "c.tsv", "--bam", "x.bam", "-o", "out",
        "--prefix", "run", "--top-n", "3", "--region", "chr1:100-400",
        "--region", "chr2:5-9", "--no-summary"])
    assert args.clusters == "c.tsv"
    assert args.bam == "x.bam"
    assert args.outdir == "out"
    assert args.top_n == 3
    assert args.region == ["chr1:100-400", "chr2:5-9"]
    assert args.no_summary is True
    assert args.func.__name__ == "run_classplots"
```

Note: the parser factory is `build_parser` (verified at `srnaprofiler.py:53`).

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cli_classplots.py -v`
Expected: FAIL (unknown subcommand `classplots`)

- [ ] **Step 3: Add the subparser in `srnaprofiler.py`**

Immediately after `cl.set_defaults(func=run_classify)` (line 608), add:

```python
    from classplots import run_classplots
    cp = subparsers.add_parser("classplots",
        help="Per-class plots from a classify .clusters.tsv + BAM")
    cp.add_argument("--clusters", required=True, help="classify .clusters.tsv")
    cp.add_argument("--bam", required=True, help="Sorted+indexed BAM used by classify")
    cp.add_argument("-o", "--outdir", required=True)
    cp.add_argument("-p", "--prefix", default="classplots")
    cp.add_argument("--top-n", dest="top_n", type=int, default=5,
                    help="Per-locus plots for the N most abundant loci per class (0 disables)")
    cp.add_argument("--classes", help="Comma-separated subset of classes to plot")
    cp.add_argument("--region", action="append", metavar="CHROM:START-END",
                    help="Plot one explicit coordinate span (repeatable), regardless of class")
    cp.add_argument("--min-length", dest="min_length", type=int, default=15)
    cp.add_argument("--max-length", dest="max_length", type=int, default=35)
    cp.add_argument("--logo-trim", dest="logo_trim", type=int, default=None,
                    help="Trim reads to N nt for the logo (default: class dominant DicerCall)")
    cp.add_argument("--no-signals", dest="no_signals", action="store_true",
                    help="Skip R/WebLogo signal plots (QC summary only)")
    cp.add_argument("--no-summary", dest="no_summary", action="store_true",
                    help="Skip the matplotlib QC summary")
    cp.set_defaults(func=run_classplots)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cli_classplots.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add srnaprofiler.py tests/test_cli_classplots.py
git commit -m "feat(cli): add classplots subcommand"
```

---

### Task 14: Docs — README + ROADMAP

**Files:**
- Modify: `README.md`, `docs/ROADMAP.md`

- [ ] **Step 1: Add a `classplots` usage block to `README.md`**

Add under the subcommands section:

```markdown
### classplots — per-class plots

Re-aggregates a classify run's alignments grouped by the assigned class and produces
publication signal plots (size distribution, pattern, ping-pong, phasing, coverage,
sequence logos) plus a matplotlib QC summary.

    srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \
        -o classplots_out --prefix run --top-n 5

Plot a specific coordinate span (regardless of class) with one or more `--region`:

    srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \
        -o classplots_out --region chr1:100-400 --region X:9000-9500

Signal plots need `Rscript` (skipped with a warning if absent) and the `weblogo` CLI for
logos (`pip install weblogo`); the QC summary is pure Python.
```

- [ ] **Step 2: Mark ROADMAP item #8 done**

In `docs/ROADMAP.md`, change the Tier 3 item 8 line to note it is implemented on the
`classplots` branch (edit the "Plots por classe" bullet to past tense with the subcommand name).

- [ ] **Step 3: Commit**

```bash
git add README.md docs/ROADMAP.md
git commit -m "docs: document classplots subcommand; mark roadmap item done"
```

---

## Self-Review Notes (addressed)

- **Phasing overhang −2:** the spec's phasing plot named a −2 overhang; during planning the
  pooled metric was changed to the phase-invariant same-strand distance spectrum
  (`offset_spectrum(ends, ends)`), which pools across loci with different absolute phase
  without smearing. The −2 correction applied only to cross-strand register pooling and is
  therefore not used here; the synthetic test asserts the peak at the phase length instead.
  This deviation was raised with the user before writing the plan.
- **Weighted floats into R:** the size-distribution and pattern R scripts receive float
  counts; base R / ggplot handle floats. Verified conceptually; the smoke test (Task 9) and
  integration run (Task 12) exercise real rendering.
- **Coverage vs frac_top overlap:** kept intentionally (positional shape is not a scalar).
- **`create_parser` name:** Task 13 flags verifying the parser factory name before importing.
