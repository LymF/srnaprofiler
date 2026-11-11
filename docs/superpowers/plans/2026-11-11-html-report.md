# HTML report (`srnaprofiler report`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new subcommand that turns a `classify` run (`.clusters.tsv` + BAM) into one self-contained HTML file a layperson can read: library overview, cluster list, cluster pages, glossary.

**Architecture:** Four small root modules (flat layout, as the rest of the repo): `report_data` (read the TSV, apply the depth floor, pick which clusters get pages), `report_text` (rule-based plain-language explanations, string table, glossary), `report_signals` (per-cluster figure data from the BAM via `class_signals`), `report_html` (assemble JSON + template into one offline file). The CLI wiring lives in `srnaprofiler.py` as subcommand `report`. Nothing in classification changes.

**Tech Stack:** Python 3, pysam, PyYAML (already used), pytest. Front end: plain HTML/CSS/JS with inline SVG, styled from `~/.claude/skills/report-template/assets/report.css` (copied into `report_assets/`).

**Spec:** `docs/superpowers/specs/2026-11-11-html-report-design.md`

## Global Constraints

- One offline HTML file: CSS, JS and data inline; no external scripts, fonts or images.
- Depth floor default **300 reads**, in **raw reads** (`reads` column), user-adjustable with `--min-reads`; applied only in the report. The pipeline's `locus_min_reads`/`fragments.tsv` are untouched.
- Cluster pages: every cluster above the floor up to **3,000**; above that the 3,000 with highest confidence and a visible notice (`--max-detail`).
- No "validate first" badge. No LLM. Order by confidence.
- Always state how many clusters the floor hid and how many passed.
- Subcommand `report` is opt-in; `classify` and every existing output stay byte-identical.
- Report text in English; every user-facing string lives in the table in `report_text.py`.
- Evidence-dot thresholds are display-only and documented as such in the module and in the report footer.
- Code and comments in English, short comments. Authorship: Lucas only (no AI credit, also not in commit messages or trailers).
- Glossary statements of biological fact carry a citation; the report ends with a References section with DOIs. Each reference must be read from its source (paper lookup), not written from memory, before it goes in.
- Tests live in `benchmark/tests/` (that is where the repo's pytest suite sits), run on litrp3 in the `srnaprofiler` conda env (local pytest is missing).
- Commits on a new branch `feat/html-report`; conventional messages; no push.

## Review Focus

1. Floor above every cluster (zero rows pass): all tables empty with a clear message, no crash, no division by zero.
2. A `.clusters.tsv` from a run older than the `e_*`/`gate_*` columns (like `results-teste/dros_ref.clusters.tsv`): the report still builds, dots show "unavailable", explanations say so instead of guessing.
3. Fields with `NA`, `N` (dicer_call) or empty strings, and class names the report has no text for (for example `phasiRNA24`, `hcsiRNA`, `other`): parsed without error, shown with a generic explanation.
4. Cluster names, chromosome names and sequences containing `<`, `&`, quotes or `</script>`: the embedded JSON and the rendered text must not break the page or inject markup.
5. A cluster whose region has zero reads in the BAM, or a chromosome absent from the BAM: its page builds with "no reads found", and the rest of the report is unaffected.

---

### Task 1: Read the TSV and apply the floor (`report_data.py`)

**Files:**
- Create: `report_data.py`
- Test: `benchmark/tests/test_html_report_data.py`

**Interfaces:**
- Produces:
  - `@dataclass Cluster`: `name, chrom, start, end, reads, reads_rpm, cls, confidence, dicer_call (int|None), size_purity, frac_top, pct_1u, pct_10a, pingpong_z, phase_pvalue, multimap_frac, unique_frac, review (bool), major_rna, evidence (str), e (dict[str, float|None]` keyed `1u 5a fractop_bid pingpong phase`), `gates (dict[str, bool])`, `has_evidence (bool)`; unparsable numbers become `None`.
  - `load_clusters(path: str) -> list[Cluster]`: tolerant of missing optional columns (`e_*`, `gate_*`, `p_real`).
  - `@dataclass Selection`: `shown: list[Cluster]` (≥ floor, sorted by confidence desc then reads desc), `detail: list[Cluster]` (first `max_detail` of `shown`), `n_total`, `n_hidden`, `floor`, `capped (bool)`.
  - `select(clusters, floor: float, max_detail: int = 3000) -> Selection`.
  - `library_summary(clusters, fragments_path: str | None) -> dict`: keys `library_reads` (median of `reads/reads_rpm*1e6` over rows with `reads_rpm > 0`, else `None`), `cluster_reads`, `unique_share` (read-weighted mean `unique_frac`), `fragments` (count, or `None`), `fragment_reads`.

- [ ] **Step 1: Write failing tests** — `test_load_tolerates_old_columns` (the TSV header of `results-teste/dros_ref.clusters.tsv` without `e_*`: `has_evidence` is `False`, `e` values `None`); `test_na_and_N_become_none` (`dicer_call="N"`, `score="NA"`); `test_select_floor` (reads 299/300/301 with floor 300 → shown excludes 299, `n_hidden` counts it); `test_select_sorted_by_confidence_then_reads`; `test_select_cap` (3,001 clusters, `max_detail=3000` → `capped` true, `len(detail)==3000`, `shown` keeps all); `test_floor_above_all` (empty `shown`, `n_hidden == n_total`); `test_library_summary_without_fragments` (`fragments is None`).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_data.py -v` — expect FAIL (module missing).
- [ ] **Step 3: Implement** the interfaces above in `report_data.py`; read with `csv.DictReader(delimiter="\t")`, one helper `_num(x)` for float/None.
- [ ] **Step 4: Run** the same command — expect all PASS.
- [ ] **Step 5: Commit** `feat(report): read clusters.tsv and apply the depth floor`.

### Task 2: Plain-language explanations (`report_text.py`)

**Files:**
- Create: `report_text.py`
- Test: `benchmark/tests/test_html_report_text.py`

**Interfaces:**
- Consumes: `Cluster`, `Selection` from Task 1.
- Produces:
  - `STRINGS: dict[str, str]` (every user-facing sentence, keyed), `GLOSSARY: list[dict]` (`term`, `text`, `cite` list of reference keys), `REFERENCES: dict[str, str]` (key → full reference with DOI).
  - `evidence_status(c: Cluster) -> dict[str, str]`: keys `size, 1u, 10a, pingpong, phasing`; values `"present" | "weak" | "absent" | "unavailable"`. Size from `size_purity` (≥0.5 present, ≥0.3 weak); the others from `c.e` (≥0.5 present, ≥0.25 weak); `"unavailable"` when the value is `None`. Thresholds are module constants `DISPLAY_PRESENT = 0.5`, `DISPLAY_WEAK = 0.25`, documented as display-only.
  - `explain(c: Cluster) -> dict`: `verdict` (one sentence for the class), `checks` (list of `{key, status, title, detail}` for size, 1U, 10A, strand, ping-pong, phasing, multimapping), `why_not` (list of `{cls, reason}` for the other classes that have a gate column set to false, or a generic line when gates are unavailable), `cautions` (list of strings), `missing` (for `ambiguous`: which evidence was absent).
  - `summary_sentence(counts: dict[str, int], floor: float) -> str`.
  - `warnings(sel: Selection, summary: dict) -> list[str]`: banner when fewer than 20 clusters pass the floor (with the floor that would keep 20, when computable), when `capped`, and when `summary["unique_share"] < 0.3`.

- [ ] **Step 1: Write failing tests** — `test_status_thresholds` (purity 0.5 → present, 0.3 → weak, 0.29 → absent; `e` of 0.25 → weak); `test_status_unavailable_on_old_tsv`; `test_explain_piRNA_names_pingpong` (a piRNA with `e_pingpong=0.9` has a `checks` entry `pingpong` with status `present`); `test_explain_ambiguous_lists_missing`; `test_explain_unknown_class_is_generic` (`cls="other"` and `"phasiRNA24"`: no exception, non-empty `verdict`); `test_summary_sentence_names_top_class`; `test_warning_few_clusters`; `test_warning_cap`; `test_every_cited_key_has_reference` (every `cite` key in `GLOSSARY` exists in `REFERENCES`); `test_strings_have_no_markup` (no `<` in `STRINGS`).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_text.py -v` — expect FAIL.
- [ ] **Step 3: Implement** the interfaces. Glossary terms: small RNA, cluster, piRNA, siRNA, miRNA, phasiRNA, 1U/10A, ping-pong, phasing, ambiguous, confidence, multimapping. Look each cited paper up before adding it to `REFERENCES` (read the abstract/page; take authors, title, journal, year, DOI from it).
- [ ] **Step 4: Run** the same command — expect PASS.
- [ ] **Step 5: Commit** `feat(report): rule-based explanations, strings and glossary`.

### Task 3: Per-cluster figure data from the BAM (`report_signals.py`)

**Files:**
- Create: `report_signals.py`
- Test: `benchmark/tests/test_html_report_signals.py`

**Interfaces:**
- Consumes: `class_signals` (`base_composition`, `pingpong_offsets`, `coverage_track`, `_oriented`) and `samio.read_region_handle`.
- Produces:
  - `cluster_signals(alns, start: int, end: int, size_lo: int = 18, size_hi: int = 35, max_bins: int = 120) -> dict` with keys: `sizes` (`{"lo", "plus": [...], "minus": [...]}`, 1/NH-weighted counts per size), `nt` (list of 26 dicts `A,C,G,U` fractions for RNA positions 1..26, using the dominant-size reads), `coverage` (`{"start", "bin", "plus": [...], "minus": [...]}` binned to at most `max_bins`), `pingpong` (list of 30 floats, overlap 1..30), `n_alns` (int). Zero alignments → all lists present and zero-filled, `n_alns == 0`.
  - `signals_for(bam_path: str, clusters, xm_as_nh: bool | None = None) -> dict[str, dict]`: opens the BAM once, one entry per cluster name; a chromosome missing from the BAM or an exception for one region gives the zero-filled entry, never aborts the rest. `xm_as_nh=None` calls `samio.prefers_xm_as_nh(bam_path)`.

- [ ] **Step 1: Write failing tests** with `collections.namedtuple("Aln", "chrom pos seq nh strand")` as in `test_locus_floor.py`: `test_sizes_weighted_by_nh` (one read with `nh=2` counts 0.5); `test_nt_position1_is_U_for_T_start` (T-first reads on `+` give `U` fraction 1.0 at position 1); `test_pingpong_peak_at_10` (a plus read and a minus read whose 5′ ends overlap by 10 → `pingpong[9]` is the maximum); `test_coverage_binned` (a 10,000 bp region gives ≤ 120 bins whose sum equals the unbinned sum); `test_zero_alignments_zero_filled`; `test_missing_chromosome` (build a BAM with `bamutil.make_bam`, ask for `chrZ` → zero-filled entry, other cluster fine).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_signals.py -v` — expect FAIL.
- [ ] **Step 3: Implement.** Round floats to 3 decimals to keep the JSON small.
- [ ] **Step 4: Run** — expect PASS.
- [ ] **Step 5: Commit** `feat(report): per-cluster figure data from the BAM`.

### Task 4: Assemble the HTML (`report_html.py`, `report_assets/`)

**Files:**
- Create: `report_html.py`, `report_assets/report.css` (copy of the skill's `report.css`, prefix and tokens kept), `report_assets/app.js`, `report_assets/template.html`
- Test: `benchmark/tests/test_html_report_html.py`

**Interfaces:**
- Consumes: Tasks 1–3 outputs.
- Produces: `build_report(sel: Selection, summary: dict, signals: dict[str, dict], title: str, out_path: str) -> None`. Embeds one JSON object (`library`, `counts`, `warnings`, `rows`, `details`, `strings`, `glossary`, `references`) in a `<script type="application/json" id="report-data">` after replacing `<` with `<` so no data can close the tag. `app.js` renders the four tabs from that object.

Front end content (from the spec, mockup `scratchpad/mockup-report-v1.html` and spike `spike-report-v2.html` are the visual reference; they are throwaway and not copied):
- Overview: summary sentence, four number cards, class donut, dominant-size histogram stacked by class, library cards (library size, unique share, cluster reads, fragments below the pipeline floor), warning banners, a "how to read" line under each chart.
- Cluster list: class filters, sort by confidence/reads/name, 50 per page, class tag, confidence bar, five evidence dots, reads; for ambiguous rows the `missing` text; row click opens the cluster page; a floor input that can only lower the view to rows already embedded (rows below the CLI floor are not embedded) with a note saying to rerun with `--min-reads` to go lower.
- Cluster page: verdict, checklist, four SVG figures from `signals`, "why not another class", cautions, dominant sequence.
- Glossary and References tabs. Footer states the display-only thresholds.

- [ ] **Step 1: Write failing tests** — `test_output_is_self_contained` (no `src="http`, no `href="http`, no `<script src`, no `@import`); `test_embedded_json_roundtrips` (parse the `report-data` block; `rows` length equals `len(sel.shown)`); `test_script_tag_cannot_be_closed` (a cluster named `</script><b>x` produces output with exactly one `</script>` per script element and the name recoverable from the JSON); `test_empty_selection_builds` (Review Focus 1); `test_old_tsv_builds` (Review Focus 2, using a Cluster list with `e` all `None`); `test_size_budget` (1,625 synthetic clusters with zero-filled signals produce a file under 6 MB).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_html.py -v` — expect FAIL.
- [ ] **Step 3: Implement** `build_report`: read the three asset files, inline CSS and JS into the template, write the file. `app.js` escapes every inserted string (`textContent` or an `esc()` helper), never `innerHTML` with raw data.
- [ ] **Step 4: Run** — expect PASS.
- [ ] **Step 5: Commit** `feat(report): self-contained HTML assembly and front end`.

### Task 5: CLI subcommand (`report`)

**Files:**
- Modify: `srnaprofiler.py` (add `report` parser next to `classplots`, around line 473)
- Create: `report_cli.py` with `run_report(args) -> None`
- Test: `benchmark/tests/test_html_report_cli.py`

**Interfaces:**
- Consumes: Tasks 1–4.
- Produces: `srnaprofiler report --clusters run.clusters.tsv --bam sample.sorted.bam -o OUT [-p PREFIX] [--min-reads 300] [--max-detail 3000] [--fragments run.fragments.tsv] [--title TEXT]` writing `OUT/<prefix>.report.html` (default prefix `report`). `--fragments` defaults to the `.fragments.tsv` beside `--clusters` when it exists. Logs with `utils.log`: counts shown/hidden, and a WARN when the cap applies.

- [ ] **Step 1: Write failing tests** — `test_report_end_to_end` (a hand-written 3-row `.clusters.tsv` plus a BAM from `bamutil.make_bam`; run `run_report` with a `SimpleNamespace`; output exists, JSON has 3 rows when `--min-reads 0`); `test_parser_accepts_flags` (build the real parser and parse the flags above); `test_classify_outputs_untouched` (running `run_report` does not modify the input TSV: compare sha256 before/after).
- [ ] **Step 2: Run** `pytest benchmark/tests/test_html_report_cli.py -v` — expect FAIL.
- [ ] **Step 3: Implement** `run_report` and the parser entry (`rp.set_defaults(func=run_report)`; import it where `run_classplots` is imported).
- [ ] **Step 4: Run** the full suite on litrp3: `pytest benchmark/tests -q` — expect the previous 475 passed, 1 skipped plus the new tests.
- [ ] **Step 5: Commit** `feat(report): srnaprofiler report subcommand`.

### Task 6 (controller, litrp3): real-data check

- [ ] Deploy HEAD to `~/srnaprofiler-ml`; run `classify` (rules mode) on the fly library used before, so `clusters.tsv` carries `e_*`/`gate_*`; run `report` with default floor.
- [ ] Check: file size, the numbers on the overview match `awk` counts on the TSV, one piRNA and one ambiguous cluster read correctly, the page works in a browser (the author has no browser access: ask Lucas to open it and report anything that looks wrong or does not load).
- [ ] Run the same on one plant library with `rules.plant.yaml` to see the plant classes' texts.
- [ ] Add a short "Report" section to `README.md` (usage, floor, what the evidence dots mean) with the references it cites; commit.
- [ ] Stop and ask Lucas before any further change to wording or thresholds.
