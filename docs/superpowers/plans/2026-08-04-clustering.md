# Clustering Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Define loci the way the field does — coverage in RPM, bilateral extension, merging conditioned on the size and strand profile — so that a locus holds enough reads for its signatures to be measurable.

**Architecture:** `clustering.py` keeps its role (turn per-position depth into regions) and gains the new merge rules as separate, pure functions. `classify.py` converts the RPM threshold once per run using the weighted library size it already computes. The per-locus size composition is measured in `features.build_record`, where the reads are already in hand.

**Tech Stack:** Python, numpy (already a dependency), pytest.

## Global Constraints

- Code, comments, docstrings, log messages and CLI help in **English**.
- Comments short — only what is strictly necessary. Context belongs in the commit message.
- Never add `Co-Authored-By` or any AI credit to commits.
- A change enters only if it holds in **both** kingdoms: Drosophila (90 piRNA clusters, `extras/validate_piclusters.py`) and Arabidopsis (6 TAS loci, `extras/validate_plant_tas.py`). Recall must not drop.
- `--mincov` becomes RPM (float), `--pad` stays an int in bp.
- Existing callers must keep working: every new parameter gets a default that reproduces current behaviour unless explicitly set.

---

## File Structure

| file | responsibility |
|---|---|
| `clustering.py` | islands, bilateral extension, profile-conditioned merge |
| `tests/test_clustering.py` | existing island tests, extended |
| `tests/test_clustering_merge.py` (new) | extension and merge rules |
| `classify.py` | RPM → depth conversion, wiring |
| `features.py` | size composition per locus; pass profiles to the merge |
| `records.py`, `reporter.py` | `size_purity`, `dicer_call_2` |
| `srnaprofiler.py` | `--mincov` as float RPM, help text |
| `rules.yaml`, `rules.plant.yaml` | `clustering.pad` per organism |

---

### Task 1: bilateral extension

**Files:**
- Modify: `clustering.py`
- Test: `tests/test_clustering_merge.py`

**Interfaces:**
- Produces: `extend_and_merge(islands, pad, chrom_len=None) -> list[tuple[int,int]]`

- [ ] **Step 1: Write the failing test**

```python
"""Bilateral extension, as ShortStack does it (bedtools slop + merge).

Each island grows by `pad` on both sides and overlapping results are merged, so
the effective tolerated gap is 2*pad -- not pad, which is what the old
gap-based rule gave.
"""
from clustering import extend_and_merge


def test_islands_closer_than_two_pads_merge():
    # gap of 150 between them, pad 100 -> each grows 100, so they overlap
    assert extend_and_merge([(1000, 1100), (1250, 1350)], pad=100) == [(900, 1450)]


def test_islands_further_than_two_pads_stay_separate():
    # gap of 250, pad 100 -> extended ranges do not touch
    assert extend_and_merge([(1000, 1100), (1350, 1450)], pad=100) == [(900, 1200), (1250, 1550)]


def test_extension_is_clamped_at_zero_and_chrom_end():
    assert extend_and_merge([(10, 50)], pad=100, chrom_len=200) == [(0, 150)]
    assert extend_and_merge([(150, 190)], pad=100, chrom_len=200) == [(50, 200)]


def test_pad_zero_returns_islands_untouched():
    assert extend_and_merge([(10, 20), (40, 50)], pad=0) == [(10, 20), (40, 50)]


def test_empty_input():
    assert extend_and_merge([], pad=100) == []


def test_touching_islands_merge_into_one():
    assert extend_and_merge([(100, 200), (200, 300)], pad=0) == [(100, 300)]
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_clustering_merge.py -q`
Expected: FAIL, `ImportError: cannot import name 'extend_and_merge'`

- [ ] **Step 3: Write the implementation**

Append to `clustering.py`:

```python
def extend_and_merge(islands, pad, chrom_len=None):
    """Grow each island by `pad` on both sides, then merge what overlaps.

    Effective tolerated gap is 2*pad. Mirrors ShortStack's slop+merge.
    """
    if not islands:
        return []
    grown = []
    for s, e in sorted(islands):
        s = max(0, s - pad)
        e = e + pad if chrom_len is None else min(chrom_len, e + pad)
        grown.append((s, e))
    merged = [grown[0]]
    for s, e in grown[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/test_clustering_merge.py -q`
Expected: PASS, 6 passed

- [ ] **Step 5: Commit**

```bash
git add clustering.py tests/test_clustering_merge.py
git commit -m "feat(clustering): bilateral extension before merging

Each island grows by pad on both sides and overlaps are merged, so the
effective gap tolerated is 2*pad. This is what ShortStack does with bedtools
slop + merge, and it treats a peak as having a zone of influence rather than
only measuring the hole between two peaks."
```

---

### Task 2: profile-conditioned merge

**Files:**
- Modify: `clustering.py`
- Test: `tests/test_clustering_merge.py`

**Interfaces:**
- Consumes: `extend_and_merge` from Task 1
- Produces: `merge_compatible(islands, profiles, pad, chrom_len=None, max_strand_diff=0.5) -> list[tuple[int,int]]` where `profiles` is a list parallel to `islands`, each entry a dict `{"size": int|None, "frac_top": float|None}`

- [ ] **Step 1: Write the failing test**

```python
"""Merging also requires a compatible profile (YASMA's criterion).

Distance alone would fuse a 24-nt neighbour with a 21-nt one, which are DCL3 and
DCL4 products -- different pathways that must stay apart.
"""
from clustering import merge_compatible


def _p(size, frac_top):
    return {"size": size, "frac_top": frac_top}


def test_same_profile_and_close_merges():
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(24, 0.5), _p(24, 0.55)]
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1450)]


def test_different_dominant_size_does_not_merge():
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(24, 0.5), _p(21, 0.5)]
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1200), (1150, 1450)]


def test_opposite_strand_bias_does_not_merge():
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(24, 0.05), _p(24, 0.95)]     # difference 0.90 > 0.5
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1200), (1150, 1450)]


def test_unknown_profile_never_blocks_merging():
    # a locus with no size call must not prevent a merge it knows nothing about
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(None, None), _p(24, 0.5)]
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1450)]


def test_chain_merges_only_through_compatible_neighbours():
    # A(24) B(24) C(21): A+B merge, C stays out even though it is near B
    islands = [(1000, 1100), (1250, 1350), (1500, 1600)]
    profiles = [_p(24, 0.5), _p(24, 0.5), _p(21, 0.5)]
    out = merge_compatible(islands, profiles, pad=100)
    assert out[0] == (900, 1450)
    assert out[-1] == (1400, 1700)


def test_empty_input():
    assert merge_compatible([], [], pad=100) == []
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_clustering_merge.py -q`
Expected: FAIL, `ImportError: cannot import name 'merge_compatible'`

- [ ] **Step 3: Write the implementation**

Append to `clustering.py`:

```python
def _compatible(a, b, max_strand_diff):
    """Whether two neighbouring profiles may be merged."""
    if a["size"] is not None and b["size"] is not None and a["size"] != b["size"]:
        return False
    ta, tb = a["frac_top"], b["frac_top"]
    if ta is not None and tb is not None and abs(ta - tb) > max_strand_diff:
        return False
    return True


def merge_compatible(islands, profiles, pad, chrom_len=None, max_strand_diff=0.5):
    """Extend and merge, but only across neighbours with a compatible profile.

    Unknown values (None) never block a merge: absence of a size call is not
    evidence of incompatibility.
    """
    if not islands:
        return []
    order = sorted(range(len(islands)), key=lambda i: islands[i])
    out = []
    cur_s, cur_e = None, None
    cur_prof = None
    for i in order:
        s, e = islands[i]
        s = max(0, s - pad)
        e = e + pad if chrom_len is None else min(chrom_len, e + pad)
        prof = profiles[i]
        if cur_s is None:
            cur_s, cur_e, cur_prof = s, e, prof
            continue
        if s <= cur_e and _compatible(cur_prof, prof, max_strand_diff):
            cur_e = max(cur_e, e)
            if cur_prof["size"] is None:
                cur_prof = prof
        else:
            out.append((cur_s, cur_e))
            cur_s, cur_e, cur_prof = s, e, prof
    out.append((cur_s, cur_e))
    return out
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `python -m pytest tests/test_clustering_merge.py -q`
Expected: PASS, 12 passed

- [ ] **Step 5: Commit**

```bash
git add clustering.py tests/test_clustering_merge.py
git commit -m "feat(clustering): condition merging on the size and strand profile

Distance alone would fuse a 24-nt neighbour with a 21-nt one -- DCL3 and DCL4
products, different pathways. Following YASMA, a merge now also requires the
same dominant size and a strand proportion within 0.5.

An unknown profile never blocks a merge: no size call is not evidence of
incompatibility."
```

---

### Task 3: profiles from the reads, and wiring the new merge

**Files:**
- Modify: `features.py` (`process_chromosome`)
- Test: `tests/test_clustering_profiles.py`

**Interfaces:**
- Consumes: `clustering.call_islands`, `clustering.merge_compatible` (Task 2)
- Produces: `features.island_profiles(islands, alns) -> list[dict]` with the same `{"size", "frac_top"}` shape `merge_compatible` expects

- [ ] **Step 1: Write the failing test**

```python
"""The merge needs a profile per island, computed from the reads in it."""
from features import island_profiles
from samio import Aln


def _reads(seq, positions, strand="+"):
    return [Aln("c", p, strand, seq, 1) for p in positions]


def test_profile_reports_dominant_size_and_strand():
    alns = _reads("A" * 24, [100, 110, 120])
    prof = island_profiles([(100, 150)], alns)
    assert prof[0]["size"] == 24
    assert prof[0]["frac_top"] == 1.0


def test_profile_of_dual_strand_island():
    alns = _reads("A" * 24, [100, 110]) + _reads("A" * 24, [120, 130], strand="-")
    prof = island_profiles([(100, 160)], alns)
    assert prof[0]["frac_top"] == 0.5


def test_island_without_reads_has_unknown_profile():
    prof = island_profiles([(5000, 5100)], _reads("A" * 24, [100]))
    assert prof[0]["size"] is None
    assert prof[0]["frac_top"] is None


def test_dominant_size_wins_over_minority():
    alns = _reads("A" * 24, [100, 110, 120]) + _reads("A" * 21, [130])
    prof = island_profiles([(100, 160)], alns)
    assert prof[0]["size"] == 24
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_clustering_profiles.py -q`
Expected: FAIL, `ImportError: cannot import name 'island_profiles'`

- [ ] **Step 3: Write the implementation**

Add to `features.py`:

```python
def island_profiles(islands, alns):
    """Dominant read length and strand fraction per island, for the merge rule."""
    sizes = [Counter() for _ in islands]
    plus = [0.0] * len(islands)
    total = [0.0] * len(islands)
    starts = [s for s, _ in islands]
    for a in alns:
        j = bisect.bisect_right(starts, a.pos) - 1
        if j < 0 or a.pos >= islands[j][1]:
            continue
        w = 1.0 / a.nh
        sizes[j][len(a.seq)] += w
        total[j] += w
        if a.strand == "+":
            plus[j] += w
    out = []
    for j in range(len(islands)):
        if not total[j]:
            out.append({"size": None, "frac_top": None})
        else:
            out.append({"size": sizes[j].most_common(1)[0][0],
                        "frac_top": round(plus[j] / total[j], 4)})
    return out
```

Then in `process_chromosome`, replace the region-building line:

```python
    islands = call_islands(depth, mincov, 0)
    profiles = island_profiles(islands, alns)
    regions = [(chrom, s, e) for (s, e) in
               merge_compatible(islands, profiles, pad, chrom_len=chrom_len)]
```

Note `call_islands` is now called with pad 0 — the merging is done by
`merge_compatible`. Add `chrom_len=0` as a parameter of `process_chromosome`
(default 0, meaning no clamp: pass `None` to `merge_compatible` when it is 0),
and have `extract_by_chrom` pass the reference length it already can read via
`pysam`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_clustering_profiles.py tests/test_clustering.py tests/test_features_signals.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add features.py tests/test_clustering_profiles.py
git commit -m "feat(features): compute the per-island profile that gates merging

call_islands now returns raw islands (pad 0) and merge_compatible does the
joining, so the merge can see the dominant size and strand of each side."
```

---

### Task 4: `--mincov` in RPM

**Files:**
- Modify: `srnaprofiler.py:604`, `classify.py`
- Test: `tests/test_classify_mincov_rpm.py`

**Interfaces:**
- Produces: `classify.mincov_to_depth(mincov_rpm: float, lib_size: float) -> int`

- [ ] **Step 1: Write the failing test**

```python
"""--mincov is read in reads per million and converted once per run.

As raw counts it did not transfer between libraries: the 10 calibrated on a
30.8M-alignment run means nothing at another depth.
"""
from classify import mincov_to_depth


def test_converts_rpm_using_library_size():
    # 1 RPM of a 10M-read library is 10 reads
    assert mincov_to_depth(1.0, 10_000_000) == 10
    assert mincov_to_depth(0.5, 10_000_000) == 5


def test_same_rpm_gives_different_depth_per_library():
    assert mincov_to_depth(1.0, 5_000_000) == 5
    assert mincov_to_depth(1.0, 40_000_000) == 40


def test_never_below_one():
    # a tiny library must still require at least one read
    assert mincov_to_depth(0.1, 1_000) == 1
    assert mincov_to_depth(0.0, 10_000_000) == 1


def test_handles_unknown_library_size():
    assert mincov_to_depth(1.0, 0) == 1
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_classify_mincov_rpm.py -q`
Expected: FAIL, `ImportError: cannot import name 'mincov_to_depth'`

- [ ] **Step 3: Write the implementation**

In `classify.py`:

```python
def mincov_to_depth(mincov_rpm, lib_size):
    """Reads-per-million threshold as an absolute depth for this library."""
    if lib_size <= 0:
        return 1
    return max(1, round(mincov_rpm * lib_size / 1e6))
```

In `run_classify`, after `lib_size` is known and before `extract_by_chrom`:

```python
    min_depth = mincov_to_depth(args.mincov, lib_size)
    log(f"--mincov {args.mincov} RPM -> minimum depth {min_depth} "
        f"for a library of {lib_size:.0f} reads", "INFO")
```

and pass `min_depth` where `args.mincov` was passed.

In `srnaprofiler.py:604`, change the argument:

```python
    cl.add_argument("--mincov", type=float, default=1.0,
                    help="Minimum coverage to seed a locus, in READS PER MILLION "
                         "(default: %(default)s). Converted to an absolute depth "
                         "using the library size, so it transfers between runs.")
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_classify_mincov_rpm.py tests/test_classify_integration.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add classify.py srnaprofiler.py tests/test_classify_mincov_rpm.py
git commit -m "feat(classify): --mincov in reads per million

As a raw count the threshold did not transfer: the 10 calibrated on a
30.8M-alignment library means something different at any other depth. In RPM
the same number carries across libraries and organisms, which is what makes a
cross-organism comparison of --pad meaningful at all.

This invalidates the earlier mincov calibration; it has to be redone in RPM."
```

---

### Task 5: `pad` per config

**Files:**
- Modify: `classify.py`, `rules.yaml`, `rules.plant.yaml`
- Test: `tests/test_config_clustering.py`

**Interfaces:**
- Consumes: `decider.load_config`
- Produces: config key `clustering.pad`; CLI `--pad` overrides it when given

- [ ] **Step 1: Write the failing test**

```python
"""pad lives in the config, because the kingdoms disagree on the optimum.

Measured: at pad 500 Drosophila loses 2 piRNA clusters while Arabidopsis is
still gaining. The CLI flag still wins when the user passes one.
"""
import pytest

from decider import load_config
from classify import resolve_pad


@pytest.mark.parametrize("path,expected", [("rules.yaml", 200),
                                           ("rules.plant.yaml", 500)])
def test_each_config_declares_its_pad(path, expected):
    assert load_config(path)["clustering"]["pad"] == expected


def test_cli_overrides_the_config():
    cfg = {"clustering": {"pad": 500}}
    assert resolve_pad(cli_pad=120, config=cfg) == 120


def test_config_used_when_cli_absent():
    cfg = {"clustering": {"pad": 500}}
    assert resolve_pad(cli_pad=None, config=cfg) == 500


def test_falls_back_when_config_is_silent():
    assert resolve_pad(cli_pad=None, config={}) == 200
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_config_clustering.py -q`
Expected: FAIL, `ImportError: cannot import name 'resolve_pad'`

- [ ] **Step 3: Write the implementation**

In `classify.py`:

```python
DEFAULT_PAD = 200


def resolve_pad(cli_pad, config):
    """CLI wins; otherwise the config's clustering.pad; otherwise the default."""
    if cli_pad is not None:
        return cli_pad
    return (config.get("clustering", {}) or {}).get("pad", DEFAULT_PAD)
```

In `run_classify`, replace uses of `args.pad` with `resolve_pad(args.pad, config)`.

In `srnaprofiler.py`, make the flag default to `None`:

```python
    cl.add_argument("--pad", type=int, default=None,
                    help="Extend each island by this many bp on both sides before "
                         "merging (default: from the config; 200 for arthropods, "
                         "500 for plants)")
```

In `rules.yaml`, at the top level:

```yaml
clustering:
  pad: 200
```

In `rules.plant.yaml`, at the top level:

```yaml
clustering:
  pad: 500
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_config_clustering.py -q`
Expected: PASS, 5 passed

- [ ] **Step 5: Commit**

```bash
git add classify.py srnaprofiler.py rules.yaml rules.plant.yaml tests/test_config_clustering.py
git commit -m "feat(config): pad per organism, since the kingdoms disagree

Measured on both: at pad 500 Drosophila drops from 85/90 to 83/90 known piRNA
clusters while Arabidopsis still gains (TAS 6/6, reads per locus 34 -> 131).
Drosophila's piCs are large regions with internal structure that aggressive
merging erases; plant loci are compact. So pad belongs with the classes and the
weights, not in a single global default."
```

---

### Task 6: miRNA must dominate a locus to claim it

**Files:**
- Modify: `classify.py` (`mark_mirna`)
- Test: `tests/test_classify_mirna_dominance.py`

**Interfaces:**
- Produces: `classify.mark_mirna(records, mirna_intervals, min_fraction=0.5)`

- [ ] **Step 1: Write the failing test**

```python
"""A miRNA only takes a locus's label if it covers most of it.

With bigger loci a 21-nt miRNA can sit inside a cluster holding much more, and
labelling the whole thing miRNA hides everything else. At pad=75 one locus over
2kb was affected; at pad=200 already six.
"""
from records import ClusterRecord
from classify import mark_mirna


def _rec(start, end):
    return ClusterRecord(chrom="c", start=start, end=end,
                         name=f"c:{start}-{end}", reads=100, cls="hcsiRNA")


def test_mirna_covering_the_locus_claims_it():
    r = _rec(1000, 1100)
    mark_mirna([r], {"c": [(1000, 1090)]})
    assert r.cls == "miRNA"


def test_mirna_inside_a_large_locus_does_not_claim_it():
    r = _rec(1000, 11000)                      # 10 kb locus
    mark_mirna([r], {"c": [(5000, 5021)]})     # 21 nt inside it
    assert r.cls == "hcsiRNA"


def test_partial_overlap_below_threshold_does_not_claim():
    r = _rec(1000, 2000)
    mark_mirna([r], {"c": [(1000, 1300)]})     # 30% of the locus
    assert r.cls == "hcsiRNA"


def test_threshold_is_configurable():
    r = _rec(1000, 2000)
    mark_mirna([r], {"c": [(1000, 1300)]}, min_fraction=0.25)
    assert r.cls == "miRNA"


def test_claimed_locus_keeps_the_old_confidence_behaviour():
    r = _rec(1000, 1100)
    mark_mirna([r], {"c": [(1000, 1090)]})
    assert r.confidence == 1.0
    assert r.review is False
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_classify_mirna_dominance.py -q`
Expected: FAIL on `test_mirna_inside_a_large_locus_does_not_claim_it` — the current code claims on any overlap

- [ ] **Step 3: Write the implementation**

Replace `mark_mirna` in `classify.py`:

```python
def mark_mirna(records, mirna_intervals, min_fraction=0.5):
    """Label a locus miRNA only when the miRNA covers most of its span."""
    for r in records:
        span = r.end - r.start
        if span <= 0:
            continue
        covered = 0
        for (s, e) in mirna_intervals.get(r.chrom, []):
            covered = max(covered, min(r.end, e) - max(r.start, s))
        if covered / span >= min_fraction:
            r.cls = "miRNA"
            r.confidence = 1.0
            r.review = False
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_classify_mirna_dominance.py tests/test_classify_integration.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add classify.py tests/test_classify_mirna_dominance.py
git commit -m "fix(classify): a miRNA must dominate a locus to take its label

mark_mirna claimed a locus on any overlap, so a 21-nt miRNA inside a 10 kb
cluster relabelled the whole cluster and hid everything else in it. Harmless
while loci were tiny -- one locus over 2kb at pad=75 -- but it grows with pad,
reaching six at pad=200."
```

---

### Task 7: report the size composition

**Files:**
- Modify: `features.py` (`build_record`), `records.py`, `reporter.py`
- Test: `tests/test_features_size_purity.py`

**Interfaces:**
- Produces: `ClusterRecord.size_purity: float`, `ClusterRecord.dicer_call_2: str`; TSV columns of the same names

- [ ] **Step 1: Write the failing test**

```python
"""A locus can hold two size populations; report it instead of hiding it.

18% of loci have dominant-size purity below 0.5, already at pad=75, so this is
biology we were discarding by reporting only the dominant DicerCall.
"""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _reads(length, n, start=100):
    return [Aln("c", start + i, "+", "A" * length, 1) for i in range(n)]


def test_pure_locus_reports_purity_one_and_no_second_size():
    r = build_record(("c", 100, 300), _reads(24, 50), 1_000_000, CFG)
    assert r.size_purity == 1.0
    assert r.dicer_call_2 == ""


def test_bimodal_locus_reports_both():
    r = build_record(("c", 100, 300), _reads(24, 60) + _reads(21, 40, start=200),
                     1_000_000, CFG)
    assert 0.55 < r.size_purity < 0.65
    assert r.dicer_call == "24"
    assert r.dicer_call_2 == "21"


def test_second_size_omitted_when_purity_is_high():
    r = build_record(("c", 100, 300), _reads(24, 95) + _reads(21, 5, start=200),
                     1_000_000, CFG)
    assert r.size_purity > 0.9
    assert r.dicer_call_2 == ""


def test_empty_locus_is_safe():
    r = build_record(("c", 100, 300), [], 1_000_000, CFG)
    assert r.size_purity == 0.0
    assert r.dicer_call_2 == ""
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `python -m pytest tests/test_features_size_purity.py -q`
Expected: FAIL, `AttributeError: 'ClusterRecord' object has no attribute 'size_purity'`

- [ ] **Step 3: Write the implementation**

In `records.py`, after `dicer_call`:

```python
    size_purity: float = 0.0
    dicer_call_2: str = ""
```

In `reporter.py`, add `"size_purity", "dicer_call_2"` right after `"dicer_call"` in `_COLS`.

In `features.py`, inside `build_record` after `r.dicer_call` is set:

```python
    in_window = {int(k): v for k, v in seq_lens.items() if k.isdigit()}
    win_total = sum(in_window.values())
    if win_total:
        ranked = sorted(in_window.items(), key=lambda kv: kv[1], reverse=True)
        r.size_purity = round(ranked[0][1] / win_total, 4)
        # A second size is only worth reporting when the locus is not dominated
        # by one population.
        if len(ranked) > 1 and r.size_purity < 0.8:
            r.dicer_call_2 = str(ranked[1][0])
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_features_size_purity.py tests/test_reporter.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add features.py records.py reporter.py tests/test_features_size_purity.py
git commit -m "feat(features): report size purity and the second size population

18% of loci hold two size populations, measured at pad=75, so this is not an
artefact of merging -- reporting only the dominant DicerCall was discarding a
real second population. The class is unchanged: the locus still gets one label,
from the dominant population."
```

---

### Task 8: validate on both kingdoms and recalibrate

**Files:**
- Modify: `README.md`, `docs/ROADMAP.md`, possibly `rules.yaml`/`rules.plant.yaml` (pad values)
- Test: the existing validators

- [ ] **Step 1: Sweep mincov in RPM on Drosophila**

```bash
for rpm in 0.5 1 2 5; do
  python srnaprofiler.py classify --sam ~/dros.sorted.bam -o ~/calib_out \
    -p dros_rpm$rpm --mincov $rpm --threads 2 --config rules.yaml
done
```

Record for each: locus count, median span, median reads per locus, review rate,
and piC recall from `python extras/validate_piclusters.py <tsv>`.

- [ ] **Step 2: Pick the RPM value and confirm the pad**

Choose the RPM that keeps piC recall at 85/90 with the highest reads per locus.
Then, with that RPM fixed, confirm `pad` 200 is still the best for Drosophila by
running 200 and 500 and comparing recall — the bilateral extension changed what
`pad` means, so the earlier sweep does not carry over.

- [ ] **Step 3: Repeat on Arabidopsis**

```bash
python srnaprofiler.py classify --sam /home/lucas/plant_data/at_flower.sorted.bam \
  -o /mnt/d/srnaprofiler-teste/plant/classify -p at_new --mincov <chosen> \
  --threads 1 --config rules.plant.yaml \
  --mirna-loci /mnt/d/srnaprofiler-teste/plant/mirdeep_pred.bed
python extras/validate_plant_tas.py \
  /mnt/d/srnaprofiler-teste/plant/classify/at_new.clusters.tsv \
  /mnt/d/srnaprofiler-teste/plant/GCF_000001735.4_TAIR10.1_genomic.gff
```

TAS must stay 6/6. If it drops, stop and report — the plant `pad` needs
adjusting before anything is documented.

- [ ] **Step 4: Check the five validation criteria from the spec**

```bash
python - <<'PY'
import pandas as pd
for lbl, path in (("Drosophila", "/home/lucas/calib_out/dros_new.clusters.tsv"),
                  ("Arabidopsis", "/mnt/d/srnaprofiler-teste/plant/classify/at_new.clusters.tsv")):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    k = d[d.review.astype(str).str.lower() == "false"]
    si = k[k.cls == "siRNA"]
    print(f"{lbl}: {len(d):,} loci | span {(d.end-d.start).median():.0f} | "
          f"reads {k.reads.median():.0f} | review {(d.review.astype(str).str.lower()=='true').mean():.3f}")
    if len(si):
        print(f"   siRNA |frac_top-0.5| median {(si.frac_top-0.5).abs().median():.3f} "
              f"(was 0.262 before this work)")
    print(f"   size_purity median {k.size_purity.median():.3f}")
PY
```

Expected: reads per locus up, `|frac_top-0.5|` of siRNA down, `size_purity`
not lower than the 0.64 measured before.

- [ ] **Step 5: Update the pad values if the sweep says so, and document**

Write the chosen RPM and pad into `rules.yaml` and `rules.plant.yaml`. In
`README.md`, replace the `--mincov` calibration table with the RPM one and state
that `--pad` now extends both sides. In `docs/ROADMAP.md`, mark the clustering
item done with the before/after numbers, and note that the mincov calibration of
2026-08-01 was superseded.

- [ ] **Step 6: Commit**

```bash
git add README.md docs/ROADMAP.md rules.yaml rules.plant.yaml
git commit -m "docs: recalibrate mincov in RPM and confirm pad on both kingdoms"
```

---

## Self-Review

**Spec coverage:** mincov in RPM — Task 4; bilateral extension — Task 1; pad per config — Task 5 (values reconfirmed in Task 8); profile-conditioned merge — Tasks 2 and 3; miRNA dominance — Task 6; size heterogeneity — Task 7; validation on both kingdoms — Task 8.

**Placeholders:** none; every code step carries its code. Task 8 deliberately leaves the *chosen* RPM and pad values open, because they are the output of the sweep — the criteria for choosing them are stated exactly.

**Type consistency:** `extend_and_merge` (Task 1) is superseded as the code path by `merge_compatible` (Task 2) but kept as a tested building block; `island_profiles` (Task 3) returns the `{"size", "frac_top"}` dicts that `merge_compatible` consumes; `mincov_to_depth` and `resolve_pad` (Tasks 4, 5) are both called from `run_classify`; `size_purity`/`dicer_call_2` (Task 7) match the `_COLS` entries.

**One risk to flag for the executor:** Task 3 changes `call_islands` to be called with pad 0, moving all merging into `merge_compatible`. If any other caller of `call_islands` relies on it merging (check `clustering.cluster_regions`, used by `features.extract`), it must be updated in the same task or its behaviour changes silently.
