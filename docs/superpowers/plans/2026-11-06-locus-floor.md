# Locus floor, read-bounded edges, fragments file and phased minimum span — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop reporting tiny leftover pieces as loci (they go to a separate `fragments.tsv`), bound every locus by its reads instead of padding, and require 4 phase cycles before a locus may be called phasiRNA/phasiRNA24.

**Architecture:** Clustering keeps its current island → profile-conditioned merge → miRNA carve-out path. Two steps are added at the end of `features.process_chromosome`, after reads are assigned: (3) each non-carved region is trimmed to the extent of its assigned reads; (1+2) a region whose read mass or read span is below the locus floor is not a locus — it is returned as a fragment and written to `<prefix>.fragments.tsv`. Seed-and-extend (option 1) reduces to this floor because `merge_compatible` already extends every island in both directions across compatible neighbours: a group with no island strong enough is exactly a group whose total evidence is below the floor. (A) is configuration only: the generic per-class `cluster_min_length` downgrade (decider, Task 13) gets values for the phased classes.

**Tech Stack:** Python 3, pysam, pytest; runs on litrp3.

**Spec:** user decisions in this conversation (2026-11-06): "planeja 1 3 2 e A" — 1 = seed/floor so small islands never start a locus; 3 = locus edges from reads, not from the midpoint of the gap; 2 = rejected islands written to a separate file, never silently dropped; A = phasiRNA/phasiRNA24 need ≥ 4 cycles of read span (PhaseTank: ≥ 4 phased reads within 84 nt). No separate spec document; ledger `.superpowers/sdd/2026-11-03-benchmark/progress.md` holds the rulings.

## Global Constraints

- Test species (Bombyx DRR377498, maize SRR15530410/11/18/19) are NEVER used to choose any value in this plan. Calibration uses only fly SRR30544370, Arabidopsis SRR24236124, rice DRR076318.
- Fixes go into the pipeline root `.py` files and `rules*.yaml`, not into post-processing.
- Every new config key is optional; absent or `null` reproduces today's behaviour exactly (so a replay with an old config is unchanged).
- A carved miRNA interval (`--mirna-loci`) is never trimmed and never becomes a fragment.
- Tests run in the throwaway venv in the scratchpad (suite at e3a91da: 277 passed, 1 skipped). Real-data runs on litrp3.
- Commits on `feat/benchmark`, conventional messages, trailer `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`. No push, no merge.
- Scientific text added to BASELINE/results docs: inline citations + `## Referências` with DOIs.

## Review Focus

1. A region whose assigned reads all start before its window (reads straddling the left edge, assigned by overlap) — trimmed bounds must be clipped to the original window so neighbouring loci never overlap.
2. A region with zero assigned reads (possible after carve-out) — it is a fragment (mass 0), not a crash on `min()` of an empty list.
3. Replay of an old `.clusters.tsv` (no `read_span`) through the decider with the new phased `cluster_min_length` — falls back to `end - start` (already implemented; pin with a test using a phasiRNA spec).
4. A chromosome where every region is a fragment — `process_chromosome` returns `([], fragments)` and `classify` still writes an empty-bodied `.clusters.tsv` with header.
5. Benchmark tiny counts after trimming — miRNA loci are ~60–120 bp by biology and must not inflate "tiny".

---

### Task 1: Locus floor, read-bounded edges and fragments in `features.process_chromosome`

**Files:**
- Modify: `features.py` (`process_chromosome`, `_chrom_job`, `extract_by_chrom`)
- Modify: `classify.py` (`run_classify`: pass the keys, write fragments, log)
- Modify: `reporter.py` (add `write_fragments_tsv`)
- Test: `benchmark/tests/test_locus_floor.py` (new); update the two existing `process_chromosome` callers in `benchmark/tests/test_clustering_phase_gate.py` (lines ~229–243 and ~379–390) to unpack the new return value.

**Interfaces:**
- Consumes: `_assign_reads(regions, alns) -> list[list[Alignment]]`; `carve_out(regions, intervals)`; config dict `ext_cfg`.
- Produces:
  - `process_chromosome(...) -> tuple[list[ClusterRecord], list[dict]]` (records, fragments). Fragment dict keys, exactly: `chrom, start, end, reads, read_span, reason` (`reason` is `"min_reads"`, `"min_span"` or `"min_reads,min_span"`; `start/end` are the trimmed bounds when trimming is on).
  - `extract_by_chrom(...) -> tuple[list[ClusterRecord], list[dict]]`.
  - `reporter.write_fragments_tsv(fragments: list[dict], path: str) -> None` — header `chrom start end reads read_span reason`, reads to 2 decimals.
  - Config keys read from `ext_cfg` (copied in `run_classify` from the `clustering:` YAML section): `locus_min_reads` (float|None), `locus_min_span` (int|None), `trim_to_reads` (bool, default `False` when absent — today's behaviour).
  - Output file `<prefix>.fragments.tsv` always written (header only when none).

Rules (decided here; the implementer does not choose):
- Order inside `process_chromosome`: islands → merge → carve → `_assign_reads` → trim → floor → `build_record` on survivors → density p-value (unchanged, uses the new `end - start`).
- A region is *carved* if it overlaps any interval in `carve`. Carved regions skip trim and floor.
- Trim: `new_s = max(s, min(a.pos for a in bucket))`, `new_e = min(e, max(a.pos + len(a.seq) for a in bucket))`; empty bucket → leave bounds as they are (it will fail the floor).
- Read mass = `sum(1/a.nh)` over the bucket (same scale as `record.reads`). Read span = `max end - min pos` over the bucket, 0 when empty (same definition as `build_record`).
- Floor: fragment if (`locus_min_reads` is not None and mass < `locus_min_reads`) or (`locus_min_span` is not None and span < `locus_min_span`).
- Fragments' reads are not reassigned to any locus (they are reported, not absorbed).
- `run_classify` logs: `"{n} fragments below the locus floor ({mass:.0f} reads, {pct:.2f}% of the assigned read mass) -> {base}.fragments.tsv"` at INFO.

- [ ] **Step 1: Write the failing tests** in `benchmark/tests/test_locus_floor.py`, using the monkeypatch pattern of `test_process_chromosome_stores_density_pvalue` (fake `read_region`, `chrom_len=1000`):
  - `test_floor_off_by_default`: no new keys → same regions/bounds as before, `fragments == []`.
  - `test_trim_to_reads_bounds_locus`: one island of reads at 200–260, `pad=50`, `trim_to_reads=True` → record `start == 200`, `end == 260`; with `trim_to_reads` absent → `start == 150`, `end == 310`.
  - `test_trim_clipped_to_window`: two strand-blocked neighbours whose split boundary cuts through a read → after trimming, `left.end <= right.start`.
  - `test_min_span_makes_fragment`: a 30-bp stack of 25 reads, `locus_min_span=50` → `records == []`, one fragment with `reason == "min_span"`, `reads == 25.0`, `read_span == 30`.
  - `test_min_reads_makes_fragment`: 5 reads over 80 bp, `locus_min_reads=10` → fragment `reason == "min_reads"`; both keys failing → `reason == "min_reads,min_span"`.
  - `test_carved_mirna_never_fragment_nor_trimmed`: carve interval 400–480 holding 3 reads at 420–442, floor 10 reads / 50 bp, trim on → a record with `start == 400`, `end == 480`.
  - `test_all_fragments_returns_empty_records`: every region below floor → `([], [..])`.
  - `test_write_fragments_tsv_header_only`: `write_fragments_tsv([], p)` → file content is exactly the header line.
- [ ] **Step 2: Run** `pytest benchmark/tests/test_locus_floor.py -v` → FAIL (tuple unpack / missing function).
- [ ] **Step 3: Implement** per the rules above; update `classify.run_classify` to unpack `(records, fragments)`, copy the three keys from the `clustering:` section into `ext_cfg` (same pattern as `strand_min_reads`), write `<base>.fragments.tsv`, log the line above. Update the two existing callers in `test_clustering_phase_gate.py`.
- [ ] **Step 4: Run** the full suite in the scratchpad venv → all pass (277 + new, 1 skipped).
- [ ] **Step 5: Commit** `feat(clustering): locus floor, read-bounded edges and fragments.tsv`.

### Task 2: Phased-class minimum span (A) and config placeholders

**Files:**
- Modify: `rules.plant.yaml` (classes `phasiRNA`, `phasiRNA24`; `clustering:` section)
- Modify: `rules.yaml` (`clustering:` section only)
- Test: `benchmark/tests/test_clustering_phase_gate.py` (append)

**Interfaces:**
- Consumes: `RuleScorer._cluster_failures(spec, record)` and the post-decision downgrade (Task 13), unchanged.
- Produces: config values the calibration (Task 4) overwrites only in the `clustering:` section.

- [ ] **Step 1: Failing tests** (load the real `rules.plant.yaml`):
  - `test_plant_phasirna_min_span_is_four_cycles`: `cfg["classes"]["phasiRNA"]["cluster_min_length"] == 84` and `cfg["classes"]["phasiRNA24"]["cluster_min_length"] == 96`; neither class has `cluster_max_density_pvalue`.
  - `test_phasirna_short_span_downgraded`: a record that wins phasiRNA with `read_span=60` → `cls == "ambiguous"`, `review is True`, evidence starts `nearest=phasiRNA; below_cluster_criterion(length); `.
  - `test_phasirna_replay_without_read_span_uses_window`: same record with `read_span=None`, `end-start=120` → stays phasiRNA.
  - `test_floor_keys_present_and_off`: both YAMLs have `clustering.locus_min_reads`, `locus_min_span`, `trim_to_reads` present with values `null`, `null`, `false` (placeholders until Task 4).
- [ ] **Step 2: Run** → FAIL.
- [ ] **Step 3: Edit YAML.** Phased classes: `cluster_min_length: 84` / `96` with the comment `# 4 phase cycles of read span (21/24 nt x 4): fewer cycles cannot show phasing; PhaseTank requires >= 4 phased reads within 84 nt (Guo et al., 2015). Below it the locus is downgraded to ambiguous/review (user decision 2026-11-06).` `clustering:` in both files: the three keys with the comment `# Locus floor (user decision 2026-11-06): a region below it is not a locus; it goes to <prefix>.fragments.tsv. trim_to_reads bounds each locus by its reads instead of padding. Values set by tuning-only calibration (fly, Arabidopsis, rice).`
- [ ] **Step 4: Run** full suite → pass.
- [ ] **Step 5: Commit** `feat(rules): phased classes need four cycles of read span`.

### Task 3: Benchmark — tiny counts on read span, fragments summary

**Files:**
- Modify: `benchmark/metrics_delim.py` (tiny counts), `benchmark/report.py` (Loci counts table), the module that loads the work directory's `.clusters.tsv` for metrics (find where `metrics_delim` gets `loci`).
- Test: `benchmark/tests/test_metrics_delim.py`, `benchmark/tests/test_report.py`

**Interfaces:**
- Produces in `metrics["delimitation"]["global"]`: `tiny_counts` now computed on `read_span` (fallback `end - start` when the locus has no `read_span`) and **excluding loci with `cls == "miRNA"`**; new keys `n_fragments` (int|None) and `fragment_read_frac` (float|None) = fragment read mass / (fragment mass + sum of locus `reads`), both None when `<prefix>.fragments.tsv` is absent (older runs).
- Report "Loci counts" table gains columns `fragments` and `fragment_reads` (fraction, 3 decimals; "—" when None).

- [ ] **Step 1: Failing tests:** tiny counts ignore a 70-bp miRNA locus and count a 70-bp read_span siRNA locus inside a 1,000-bp window; `fragment_read_frac == 0.2` for loci reads [40, 40] and fragments reads [10, 10]; absent file → both None; report renders the two columns and "—".
- [ ] **Step 2: Run** → FAIL. **Step 3:** implement; add one sentence to BASELINE.md "How to read" (section 2): tiny counts use the read span and exclude miRNA loci from this commit on; earlier numbers used the window. Do not change dated numbers.
- [ ] **Step 4: Run** full suite → pass. **Step 5: Commit** `feat(benchmark): tiny counts on read span; report fragments` (benchmark/ is outside the pipeline hash).

### Task 4: Calibration on tuning libraries (controller, litrp3)

Not a subagent task: run `classify` only, on fly SRR30544370, Arabidopsis SRR24236124, rice DRR076318, from a checkout at Task 3's commit, in `~/calib14/<config>/`, then `python -m benchmark metrics` per library (reuse the calib12 symlink layout and `cmp12.py`).

Grid (both YAMLs edited per config; phased minimum from Task 2 on in all):
| config | trim_to_reads | locus_min_reads | locus_min_span |
|---|---|---|---|
| t14_off | false | null | null |
| t14_trim | true | null | null |
| t14_r10s30 | true | 10 | 30 |
| t14_r10s50 | true | 10 | 50 |
| t14_r20s50 | true | 20 | 50 |

Selection rule (fixed before looking): take the strictest config such that, versus `t14_trim`, no class loses more than 0.01 absolute recall or 0.01 precision in any of the three libraries, and `fragment_read_frac <= 0.02` in all three. Ties → fewer tiny loci. If `t14_trim` itself costs > 0.01 recall versus `t14_off` in any class, rule on trim separately and record it.

- [ ] Record the table and the ruling in the ledger; set the chosen values in both YAMLs; full suite; commit `chore(rules): locus floor from tuning calibration` with the table summary in the body.

### Task 5: Full run #5 (all 14 libraries) and results document

- [ ] Ask the user before launching (cost ~10 h on litrp3). Run on the Task 4 commit with the run #4 driver pattern (`~/bench_all5.sh`, status `~/bench_all5.status`), then `benchmark report`.
- [ ] Write `benchmark/RESULTS-run5.md` comparing run #4 and run #5 for tuning and test (test is reported, never tuned on), including the note that per-locus precision depends on fragmentation (Bombyx piRNA: tp 1670 → 321, fp 71 → 60 between BASELINE and run #4), with `## Referências` (ShortStack, proTRAC, PICB, PhaseTank, PHASIS with DOIs from `lit-cluster-thresholds.md`).
