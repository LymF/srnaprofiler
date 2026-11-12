# HTML report v2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the cluster page into a figure-first page in the lab's R-plot style (sense up, antisense mirrored down), drop the glossary and the sentence checklist, add an About tab, and make the build fast and the file small.

**Architecture:** `report_signals` is rewritten to read each cluster's alignments once from the BAM and aggregate into counters (no per-read list), producing all figure data; `signals_for` can fan out over clusters. `report_html` ships the data in a compact quantised format, with the library-context percentiles from `report_data`. `app.js` draws the figures as inline SVG and the About tab; a JavaScriptCore smoke test (macOS) exercises `app.js` with a stub DOM.

**Tech Stack:** Python 3, pysam, pytest, inline SVG/JS. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-11-12-html-report-v2-design.md` (amends `2026-11-11-html-report-design.md`; v1 is built on this branch, `feat/html-report`).

## Global Constraints

- Figures: sense (+) up, antisense (-) mirrored down. 5' base colours A `#108041`, C `#3B54A2`, G `#F9A31B`, U `#EA2627`; ping-pong bar 10 in `#e41a1c`, other bars `grey60`; phasing: red dashed lines at multiples of the period (21 or 24 by class, else the cluster's dominant size), distances 1..100.
- Cluster page text: one-line verdict + a strip of measured values; no sentence checklist, no "how to read", no "why not" paragraphs. Hover tips one short line each.
- Glossary tab, References tab and all `GLOSSARY`/`REFERENCES`/`cite` code and tests are removed.
- About tab first: title "sRNAProfiler", tagline "Small RNA loci discovery and classification from sequencing data", version string shared with `srnaprofiler --version`, GitHub card `https://github.com/LymF/srnaprofiler` ("Source code & documentation"), install snippet `conda env create -f environment.yml` / `conda activate srnaprofiler`; no paper, PyPI, logo or citation block. CSS `rpt-about*` already in `report_assets/report.css`.
- Default `--max-detail` 1000; default `--min-reads` 300; table keeps every cluster above the floor.
- Target: default report on the Drosophila library used before (SRR30544370, floor 300, 1000 pages) at about 10 MB or less; generation time reported in the final check (v1 took ~13 min for 3000 pages).
- One offline file; every string in `report_text.STRINGS` (translation-ready); no `innerHTML`/`insertAdjacentHTML`/`document.write` with data; all figures inline SVG; no external resources.
- Code and comments in English, short comments. Authorship Lucas only: no AI credit in files, commit messages or trailers (this overrides any harness attribution reminder).
- Tests in `benchmark/tests/test_html_report_*.py`; run locally with `/private/tmp/claude-501/-Users-lucasmelo/a3bb3550-342f-436b-8a40-23d781c546a4/scratchpad/venv/bin/python -m pytest benchmark/tests/test_html_report_*.py -q`; the full suite runs on litrp3 (`ssh litrp3`, `source ~/miniforge3/etc/profile.d/conda.sh; conda activate srnaprofiler`; ssh drops on long jobs: use nohup and poll).
- Commits on branch `feat/html-report`; conventional messages; no push.

## Design brief (from the design skills; binds Tasks 3, 5 and 7)

Sources read for this: `dataviz`, `scientific-visualization`, `ui-ux-pro-max`, `web-design-guidelines`. `taste-skill` (marketing-site toolchain) and `frontend-design` (bold aesthetics) were read and NOT applied beyond the principles below: the report keeps the ViralQuest look (`report-template`) and the lab's R conventions.

- **Lab colours win, with relief.** The 5' base colours and the sense/antisense colours (brown `#996631` up, blue `#0099CC` down, from `class_coverage.py`) are the lab's convention and stay. `validate_palette.py` on A/C/G/U (`#108041,#3B54A2,#F9A31B,#EA2627`, light mode) gives: lightness band FAIL for G `#F9A31B`; contrast vs surface WARN 1.99:1 for G; CVD separation PASS (worst adjacent deutan 17.8, tritan 7.9, which is in the 6-8 band that is legal only with secondary encoding); normal-vision floor PASS (23.1). The brown/blue pair: CVD and normal PASS, chroma floor FAIL for the brown (0.098). Relief required by the skill and therefore part of this task: a legend (A, C, G, U letters), per-mark tooltips with the exact value, base letters printed inside a segment only when they fit, 2px surface gaps between stacked segments, and a table view for every figure.
- **Marks:** bars at most 24px thick, 4px rounded data end and square at the baseline, 2px surface gap between touching segments (no stroke borders), lines 2px with round join/cap, markers >= 8px with a 2px surface ring, hairline solid recessive grid/axes. The red dashed period lines in the phasing panels are reference markers (the lab's `plotPhasing.R` convention) and stay dashed; gridlines never dashed.
- **Text:** never in the series colour (labels, values, legends use text tokens); sentence case; units in axis labels ("small RNA length (nt)", "5'-5' distance (nt)", "Reads"); on-screen text at least 12 px; one short panel title per figure with the panel letter (A-H) and the read count n it is drawn from. Direct labels only for the dominant size/peak (selective), never a value on every bar.
- **Interaction (dataviz + ui-ux):** every mark is its own tooltip target (`pointermove` and `focus`), the hovered mark lifts; hit areas larger than the mark (at least 24 px); tooltips enhance, never gate: every value is also in the table view; keyboard focus shows the same as hover; visible focus ring; reduced motion respected (no animation, nothing hidden); touch targets 44 px for controls; no information by colour alone.
- **Table view and export for every figure** (the "writing" use case): a toggle that shows the plotted values as a table (header states the precision: values reproduce the drawn ones within 1% of the panel maximum), and a "Download SVG" action per figure (serialise the inline SVG; the template's `charts.js` `exportSVG` is a reference). The sequences table stays exact.
- **Form/axes:** one axis per panel (no dual axes); sense up / antisense down share one axis with a zero line; panels of the same kind share the y scale across strands so the halves are comparable; zero baseline for bars; a strand with no reads draws an empty half with a "no reads on this strand" note, not a blank.
- **Page chrome:** one filter/control row above the content it scopes (class filter, sort, in-page floor); loading never shifts layout; reserve figure space so nothing jumps; contrast of text tokens at least 4.5:1 (check light and the dark-mode steps if the stylesheet has them); no emoji icons, SVG only.
- **Anti-patterns to check against at review:** `dataviz/references/anti-patterns.md` (thick saturated blocks, number on every point, border around marks, clipped labels, fixed height without the axis band, tooltip as the only way to read a value).

## Review Focus

1. A deep locus (millions of alignments, mostly multimappers): memory bounded by distinct (position, strand, length, 5' base) keys, not by read count; no per-read list.
2. No reads / failed region / chromosome absent from the BAM: the page still draws the measured-value strip from the TSV and shows a "no reads found" message in place of the figures.
3. A sense-only or antisense-only cluster (and a strand with zero reads in the logo): mirrored panels render with an empty half; no division by zero, no NaN in the JSON.
4. Very short cluster (fewer bases than bins), megabase-long cluster, and fewer than 10 distinct sequences: coverage/5'-end bins stay bounded and the sequences table just shows fewer rows.
5. Realistic payload size: the size-budget test must use non-zero, realistic magnitudes (v1's budget test under-measured with zero-filled data); a cluster name with `<`, `&`, quotes or `</script>` stays escaped.

---

### Task 1: Single-pass figure data (`report_signals.py`)

**Files:**
- Modify: `report_signals.py` (replace `cluster_signals` internals; keep module small)
- Test: `benchmark/tests/test_html_report_signals.py` (rewrite for the new shape; keep the failed-region and missing-chromosome tests)

**Interfaces:**
- Consumes: `samio` (alignment filtering `_keep`, NH handling via `prefers_xm_as_nh`, see `read_region_handle`/`_to_aln` for conventions), `class_signals` (reference implementations used by the tests only).
- Produces: `cluster_signals(fh, chrom, start, end, *, xm_as_nh, dominant=None, size_lo=18, size_hi=35, max_bins=400, top_seqs=10) -> dict` reading the region once from an open `pysam.AlignmentFile` (no list of alignments). Keys, all floats rounded to 3 decimals, never NaN:
  - `n_alns` (int, records kept), `dominant` (int|None; `dominant` argument or the most abundant length inside the window).
  - `sizes`: `{"lo": size_lo, "A"|"C"|"G"|"U": {"+": [per-size weighted counts], "-": [...]}}` (RNA-oriented 5' base; 1/NH weights), same numbers as `class_signals.base_composition`.
  - `phasing`: `{"+": [100 values, distances 1..100], "-": [...]}` equal to `class_signals.phasing_series`; `pingpong`: 30 values equal to `class_signals.pingpong_offsets`.
  - `cov`: `{"start", "bin", "+": [...], "-": [...]}` mean raw depth per base in each bin (difference array, not weighted by 1/NH, as `coverage_track`), at most `max_bins` bins; `ends`: `{"bin", "+": [...], "-": [...]}` 1/NH-weighted counts of 5' end positions per bin (same binning).
  - `logo`: `{"+": [26 x {"A","C","G","U"} fractions], "-": [...]}` over reads of the dominant size, per strand, RNA-oriented, T as U; all zero where a strand has no reads.
  - `seqs`: up to `top_seqs` rows `[sequence (RNA, oriented), length, strand, weighted reads, 5' base]`, most abundant first, ties by sequence.
  - `signals_for(bam_path, clusters, xm_as_nh=None, threads=1)` keeps its contract (one entry per cluster name; failed or coordinate-less clusters get `{"failed": True}`; chromosome absent from the BAM gets `n_alns == 0`; one WARN with count and first error), now calling the new function.

- [ ] **Step 1: Write the failing tests** using small BAMs from `benchmark/tests/bamutil.py` (extend `make_bam` or add a helper in the test for reads with chosen sequence, strand and NH; keep the old helper working): `test_matches_class_signals_reference` (build the equivalent `Aln` list and assert `sizes` sums per 5' base equal `class_signals.base_composition`, `phasing` equals `phasing_series` for both strands, `pingpong` equals `pingpong_offsets`, and the sum of `cov` x bin equals the sum of `coverage_track`); `test_nh_weighting`; `test_minus_strand_is_reverse_complemented_in_logo_and_seqs`; `test_sense_only_cluster_has_zero_antisense` (all `-` arrays zero, no NaN); `test_bins_bounded` (short region: bins <= length; 1 Mb region: <= 400 bins); `test_seqs_top_and_fewer_than_ten`; `test_zero_reads` (`n_alns == 0`, all arrays zero-filled); `test_memory_is_bounded_by_distinct_keys` (200,000 identical records produce one counter key: assert via a counting wrapper or by timing-free structural check that no list of alignments is built; describe the check in the test); keep `test_missing_chromosome` and the WARN tests.
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_signals.py -q` — expect FAIL.
- [ ] **Step 3: Implement** `cluster_signals` per the interface: one pass over `fh.fetch(chrom, start, end)` applying the same record filter as `samio._keep`, accumulating `Counter` keys `(strand, pos, length, nh, first_base)`; derive sizes, 5'-end positions, coverage (difference array over the cluster window, clipped like `coverage_track`), ping-pong and phasing from those counters; get the sequence only as needed (5' base for every kept record; full sequence only for dominant-size records, collapsed by sequence). Reuse `class_signals.offset_spectrum`-style logic by calling `signals.offset_spectrum` on the end Counters as `class_signals` does.
- [ ] **Step 4: Run** the same command — expect PASS; also run `pytest benchmark/tests/test_html_report_cli.py -q` (it calls `signals_for`).
- [ ] **Step 5: Commit** `feat(report): single-pass figure data per cluster`.

### Task 2: Parallel figure data and the version single source

**Files:**
- Modify: `report_signals.py` (`signals_for(..., threads=1)`), `report_cli.py`, `srnaprofiler.py` (parser `report`: `--threads`; `--version` string), Create: `version.py`
- Test: `benchmark/tests/test_html_report_cli.py`, `benchmark/tests/test_html_report_signals.py`

**Interfaces:**
- Consumes: Task 1's `cluster_signals`.
- Produces: `signals_for(bam_path, clusters, xm_as_nh=None, threads=1)` — with `threads > 1`, clusters are split in chunks over a `concurrent.futures.ProcessPoolExecutor`, each worker opens its own `AlignmentFile`; results merged into the same dict, same WARN semantics; `report --threads N` (int >= 1, default 1; invalid values give a clear error, no traceback); `version.py` with `__version__ = "1.0"` used by `srnaprofiler --version` and the report meta (`meta.version`).

- [ ] **Step 1: Write the failing tests:** `test_threads_give_identical_signals` (same BAM, threads 1 vs 2: equal dicts); `test_failed_cluster_still_zero_filled_with_threads`; `test_threads_must_be_positive` (CLI exits with a clear message for 0 and -1); `test_version_single_source` (the string printed by `srnaprofiler --version` contains `version.__version__`).
- [ ] **Step 2: Run** the affected test files — expect FAIL.
- [ ] **Step 3: Implement** per the interface (small chunks of ~50 clusters; deterministic merge order).
- [ ] **Step 4: Run** the affected test files — expect PASS.
- [ ] **Step 5: Commit** `feat(report): --threads for figure data and a shared version string`.

### Task 3: Compact payload, library context, About data; remove the glossary (Python side)

**Files:**
- Modify: `report_html.py`, `report_data.py`, `report_text.py`
- Test: `benchmark/tests/test_html_report_html.py`, `benchmark/tests/test_html_report_data.py`, `benchmark/tests/test_html_report_text.py`

**Interfaces:**
- Consumes: Task 1 signal dict; Task 2 `version.__version__`.
- Produces:
  - `report_data.library_context(clusters) -> dict[str, dict]` keyed by cluster name: `{"reads_pct": 0-100 int (percentile of reads among all shown clusters), "conf_pct": 0-100 int | None (percentile of confidence within the cluster's class; None for ambiguous or a class with one cluster)}`; ties share the mean rank.
  - `report_html._signals(sig)` emits the quantised format used by Task 5: for each array a `[max, [ints 0..99]]` pair (`max` the float maximum, ints = round(99 * v / max), trailing zeros dropped); `logo` as ints 0..99 per base; `seqs` as rows; `phasing`, `pingpong`, `sizes` (four 5' bases x two strands), `cov`, `ends` as above; `{"err":1}`/`{"none":1}` markers unchanged. A round-trip helper `dequantise(pair)` in `report_html` (test-only use) bounds the error by `max/198`.
  - Payload: no `glossary` or `references` keys; `checks`/`note`/`why`/`cautions` sentences replaced by a `strip`: list of `[key, level, value]` (value as number or null) for the measured values and the `v` verdict key; `meta` gains `version`, `github`, `tagline`, `install` (strings from `STRINGS`/`version.py`, not hard-coded in JS) and `ctx` (per-row `[reads_pct, conf_pct]` inside `details`).
  - `report_text`: delete `GLOSSARY`, `REFERENCES`, the glossary/reference tests and the sentence builders no longer used (`explain` keeps only what the verdict line and the strip need); `ui.*` keys for About, figure titles/axes and hover tips are added (one short line each).

- [ ] **Step 1: Write the failing tests:** `test_payload_has_no_glossary_or_references`; `test_quantise_roundtrip_error_bound`; `test_quantise_all_zero_array` (max 0, no division by zero, no NaN); `test_library_context_percentiles` (ties, a class with one cluster -> `conf_pct` None, ambiguous -> None); `test_strip_has_measured_values_for_old_tsv` (old TSV: values None, levels "unavailable"); `test_size_budget_realistic` (1000 synthetic clusters whose signal arrays have realistic random magnitudes and full lengths: file under 12 MB; assert the per-cluster JSON stays under a stated bound, e.g. 6 KB); `test_name_with_script_tag_still_escaped` (keep from v1); `test_about_strings_in_table` (about title/tagline/github/install keys exist, no `<`).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_html.py benchmark/tests/test_html_report_data.py benchmark/tests/test_html_report_text.py -q` — expect FAIL.
- [ ] **Step 3: Implement** the interfaces; one place for the quantiser (`_quant(values) -> [max, ints]`).
- [ ] **Step 4: Run** the same command — expect PASS.
- [ ] **Step 5: Commit** `feat(report): compact figure payload, library context, no glossary`.

### Task 4: JavaScriptCore smoke harness for `app.js`

**Files:**
- Create: `benchmark/tests/jsc_harness.js` (stub DOM + runner), `benchmark/tests/test_html_report_js.py`

**Interfaces:**
- Consumes: `report_html.build_report` output (the page's `report-data` JSON and `app.js`).
- Produces: a pytest module that, on macOS with `osascript` available, runs `app.js` under `osascript -l JavaScript` with a stub DOM on these datasets and asserts no exception and the expected section/SVG/row counts: full (many clusters, with and without detail), empty selection, capped, ambiguous-only, plant-like (hcsiRNA/phasiRNA), a cluster named `</script><b>"x"&`, an old-TSV dataset (values null), a no-reads and a failed-region cluster. `pytest.mark.skipif` when `osascript` is missing. Opening a cluster (calling the open function through the stub's click) must be exercised for each dataset.

- [ ] **Step 1: Write the harness and tests** so they fail against a deliberately broken copy of app.js (add one test that injects a syntax error / thrown exception into a temp copy and expects the harness to fail, proving the harness detects errors).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_js.py -q` — expect FAIL for the broken-copy detection to be proven, PASS for the rest on the current app.js (v1 front end).
- [ ] **Step 3: Implement** the stub DOM (createElement, createElementNS, createTextNode, appendChild, setAttribute, addEventListener, classList, style, querySelector(All), getElementById, textContent, children traversal sufficient for `app.js`), plus a way to read rendered text and count `svg`/`rect`/`path` elements.
- [ ] **Step 4: Run** `pytest benchmark/tests/test_html_report_js.py -q` — expect PASS.
- [ ] **Step 5: Commit** `test(report): JavaScriptCore smoke harness for the front end`.

### Task 5: Figure-first cluster page, About tab, removal of glossary and checklist (front end)

**Files:**
- Modify: `report_assets/app.js`, `report_assets/template.html`, `report_assets/report.css` (appended `rpt-srna-` block only), `report_text.py` (strings)
- Test: `benchmark/tests/test_html_report_js.py`, `benchmark/tests/test_html_report_html.py`

**Interfaces:**
- Consumes: Task 3 payload (quantised signals, strip, `ctx`, About meta), Task 4 harness.
- Produces: the cluster page with, in order: header (location 1-based, class tag, confidence or "n/a (nearest: X)"), one-line verdict, measured-value strip with levels and context markers and the library-context percentiles, then the figures: size x strand x 5' base (stacked bars, mirrored), phasing (two panels, red dashed period lines), whole-cluster coverage (mirrored), ping-pong (bar 10 red), 5' end track (mirrored), nucleotide logo (two strand panels), top-sequences table; "no reads found"/error message in place of the figures when `none`/`err`. About tab first, the Glossary and References tabs gone, no glossary/reference strings left in `STRINGS`.

- [ ] **Step 1: Write the failing tests:** extend the JS harness test to assert per cluster page the presence of the expected figure containers (counts of mirrored bar groups, two phasing panels, 8 figure blocks in total, a table with at most 10 rows), the About tab as first tab with the GitHub link, and absence of Glossary/References tabs; assert per figure a table-view toggle and a Download SVG control exist, that the legend lists A/C/G/U, that every bar/segment carries a tooltip handler and a focus target, and that the antisense-only/sense-only datasets render a \"no reads on this strand\" note; add `test_app_js_has_no_glossary_code` and keep the strings-table tests (every key used exists; no hard-coded English phrase in app.js; no `<` in strings; no `innerHTML`).
- [ ] **Step 2: Run** the JS and HTML test files — expect FAIL.
- [ ] **Step 3: Implement** the figures as small SVG helpers (`mirrorBars`, `lineTwoPanel`, `trackMirror`, `logoPanel`, `seqTable`) reading the quantised arrays (multiply back by `max/99`); keep text to the verdict and strip; follow the Design brief above for marks, colours with their relief, text rules, per-mark tooltips (pointer and keyboard focus), panel letters and n, the table-view toggle and the per-figure \"Download SVG\" action; hover tips via `title` attributes (strings table).
- [ ] **Step 4: Run** the JS and HTML test files — expect PASS.
- [ ] **Step 5: Commit** `feat(report): figure-first cluster page and About tab`.

### Task 6: CLI defaults and README

**Files:**
- Modify: `srnaprofiler.py` (`--max-detail` default 1000), `report_cli.py`, `README.md`
- Test: `benchmark/tests/test_html_report_cli.py`

**Interfaces:**
- Produces: `--max-detail` default 1000 (still validated); README "Report" section updated: default 1000 pages, `--threads`, figures list, About tab, no glossary, the display-only thresholds footer statement kept.

- [ ] **Step 1: Write failing tests:** `test_max_detail_default_is_1000`; `test_readme_mentions_threads_and_no_glossary` (simple text asserts on README).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_cli.py -q` — expect FAIL.
- [ ] **Step 3: Implement** the default and README edits; run the full local report test set.
- [ ] **Step 4: Run** `pytest benchmark/tests/test_html_report_*.py -q` — expect PASS.
- [ ] **Step 5: Commit** `feat(report): 1000 detail pages by default, README`.

### Task 7 (controller, litrp3): real-data check

- [ ] Deploy HEAD to `~/srnaprofiler-report3` (new directory; never write elsewhere except `~/report-v2`), run `report` on fly SRR30544370 and Arabidopsis SRR24236124 (existing `clusters.tsv` + `bowtie.bam` under `~/srnaprofiler/data/benchmark/results/<acc>/work/`) with default flags and `--threads 4`, in the background with nohup, polling.
- [ ] Record wall time and file size against the targets (about 10 MB, v1 was 13 min); run the full suite on litrp3 (`pytest benchmark/tests -q`).
- [ ] Copy both reports to the scratchpad and open them; Lucas checks the look (the author has no browser access). Stop and ask before any further wording/figure change.

### Task 8 (controller): design review of the finished page

- [ ] Re-run `validate_palette.py` on the colours actually used (light, and dark if the stylesheet defines it) and record the result in the final report next to the relief that covers each WARN/FAIL.
- [ ] Run the `web-design-guidelines` review (fetch the guidelines from its source URL with WebFetch, as the skill says) on `report_assets/app.js`, `report_assets/template.html` and the appended block of `report_assets/report.css`; read the findings; fix the real ones in one dispatch, with a scoped re-review.
- [ ] Check `dataviz/references/anti-patterns.md` against the two real reports; list anything that matches.
- [ ] Stop and ask Lucas before any change of the lab colours or wording.
