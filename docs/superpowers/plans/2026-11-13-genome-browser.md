# Genome-browser view per locus Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Each detailed locus in the HTML report gets a "Genomic context" figure: locus ± flank on one coordinate axis, with sRNA coverage by strand and size class, gene models, transposons, structural RNA and PHAS tracks, zoom and pan, tooltips and a table alternative. A last, labelled-experimental task adds predicted ORFs.

**Architecture:** `report` takes the same repeatable `--annotation label:path` flag as `classify` plus `--flank` (default 2000). A new `report_genome.py` slices features from the annotation layers per locus window and builds a compact payload under a short key; `report_signals.py` gains a window-coverage function; `app.js` draws the figure with the existing `figure()` helper in inline SVG. The figure exists only when annotation was given, so plain reports keep their data keys.

**Tech Stack:** Python 3 (pysam, existing `annotation.py`), vanilla JS SVG in `report_assets/app.js`, pytest with the stub-DOM approach already used in `benchmark/tests/test_html_report_js.py`, headless Chrome for real-browser screenshots.

**Spec:** `docs/superpowers/specs/2026-11-13-genome-browser-design.md`

## Global Constraints

- Self-contained single HTML, offline, no external library or CDN.
- Every string that comes from an annotation file (names, ids, families) is inserted as a text node only; no `innerHTML`.
- A report built without `--annotation` has no `gm` keys and no genomic-context figure; its data payload is unchanged (the HTML bytes change only through the shipped JS/strings).
- Flank: `--flank N` bp on each side, default `2000`, `N >= 0`, window clamped to the contig `[0, length)`.
- Report size budget 50 MB (existing); features per track per window are capped and the cap is shown ("N more not shown").
- The figure never changes a class, confidence or review flag.
- The report has no dark theme today (found 2026-11-13): colours come from the existing tokens and size palette, light only.
- Tests: pytest in `benchmark/tests/`, synthetic fixtures only; run `/private/tmp/claude-501/-Users-lucasmelo/98c1215c-4a0e-4faa-b4de-b553a40b3bd2/scratchpad/venv/bin/python -m pytest benchmark/tests -q` (842 passed, 11 skipped at the start of this plan).
- Real-browser rule: every JS task states in its report what was checked in headless Chrome (`/Applications/Google Chrome.app/Contents/MacOS/Google Chrome --headless=new --disable-gpu --hide-scrollbars --window-size=1400,H --virtual-time-budget=20000 --screenshot=<png> file://<html>`); the report has no URL routing, so use a temporary copy in the scratchpad dir with an appended `<script>` that opens the wanted tab/locus. Nothing temporary stays in the repo.
- Commit messages end with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.

## Review Focus

1. Locus near a contig start or end: window clamps to `[0, length)`, ruler and coverage still align.
2. A gene much longer than the window (e.g. 200 kb gene, 5 kb window): drawn clipped with a continuation arrow, never stretched or dropped.
3. Window with no feature in any layer, and a layer file that was given but has no feature on this contig: the figure says "nothing annotated in this window" per track, not silence.
4. Gene-dense window (dozens of transcripts of one gene, overlapping genes): one row per gene, no overlapping boxes, cap respected, "N more not shown".
5. Unstranded features (`.`) and BED rows without name: drawn without arrow, tooltip says "strand unknown".

---

## File Structure

- Modify `annotation.py`: `Feature` gains `parent` and `label`; parsers fill them.
- Create `report_genome.py`: window computation, feature slicing, gene-model grouping, compact payload.
- Modify `report_signals.py`: `window_tracks` (coverage by strand and size class over a window).
- Modify `report_cli.py`, `srnaprofiler.py`: `--annotation`, `--flank` on `report`.
- Modify `report_html.py`, `report_text.py`: `gm` in the detail dict, contract comment, strings.
- Modify `report_assets/app.js`, `report_assets/report.css`: the figure, interaction.
- Create `report_orfs.py` (Task 7, experimental).
- Tests: `benchmark/tests/test_genome_*.py`.

---

### Task 1: Feature structure (parent and label)

**Files:**
- Modify: `annotation.py`
- Test: `benchmark/tests/test_genome_features.py`

**Interfaces:**
- Produces: `Feature = namedtuple("Feature", "start end strand ftype name fid parent label", defaults=("", "", ""))` (new trailing fields `parent`, `label`, both `""` when absent; existing positional uses of the first six fields keep working).
  - GFF3: `parent` = first value of `Parent`; GTF: `parent` = `transcript_id` for exon-like rows, `gene_id` for transcript rows, `""` for gene rows. `label` = first non-empty of `Name`, `gene_name`, `gene`, `locus_tag`, `ID`/`gene_id` (GTF), after the existing percent-decoding.
  - BED: `label` = column 4 when present and not `.`; RepeatMasker `.out`: `label` = repeat name (column 10), `name` stays `class/family`.
- Consumes: existing parsers in `annotation.py`.

- [ ] **Step 1: Write failing tests** (`test_genome_features.py`): GFF3 Ensembl-style `gene`/`mRNA`/`exon` lines give `parent` chains (`exon.parent == "transcript:T1"`, `mRNA.parent == "gene:G1"`, `gene.parent == ""`) and `label` from `Name` (fall back to ID); NCBI-style `Parent=rna-XM_1` and `gene-LOC1`; GTF exon `parent == transcript_id`; BED 4-column `label`; RM `.out` `label == "Gypsy-1_DM-I"` while `name == "LTR/Gypsy"`; every existing parser test still passes; `Feature(1, 2, "+", "exon", "")` (old six-field construction) still works with `parent == ""`.
- [ ] **Step 2: Run** `python -m pytest benchmark/tests/test_genome_features.py -v` → FAIL.
- [ ] **Step 3: Implement** the extra fields in `Feature` and the four parsers; interning of the new strings like `name`.
- [ ] **Step 4: Run** the new file and the full suite → PASS.
- [ ] **Step 5: Commit** `git add annotation.py benchmark/tests/test_genome_features.py && git commit -m "feat(annotation): features keep parent and label for gene models"`.

---

### Task 2: Window features and gene models (`report_genome.py`)

**Files:**
- Create: `report_genome.py`
- Test: `benchmark/tests/test_genome_windows.py`

**Interfaces:**
- Consumes: `annotation.Layer` (`.kind`, `.index[chrom].overlaps(start, end)`), `annotation.Feature` (Task 1).
- Produces:
  - `FEATURE_CAP = 60`
  - `window_for(start: int, end: int, flank: int, contig_len: int) -> tuple[int, int]` → `(max(0, start - flank), min(contig_len, end + flank))`.
  - `window_features(chrom: str, ws: int, we: int, locus: tuple[int, int], layers: list) -> dict` returning
    `{"g": [[label, strand, s, e, hl, [[t, s, e], ...]], ...], "t": [[label, family, strand, s, e, hl], ...], "r": [[label, ftype, strand, s, e, hl], ...], "p": [[label, ftype, strand, s, e, hl], ...], "mo": {"g": n, "t": n, "r": n, "p": n}, "ly": ["gene", "te", ...]}`
    - Coordinates are 0-based half-open genomic positions, features clipped never: clipping is the renderer's job, so a feature longer than the window keeps its real `s`/`e`.
    - `hl` is 1 when the feature overlaps the locus `[locus[0], locus[1])`, else 0.
    - `g` (gene layer): one entry per GENE: group exon-like rows (`exon`, `CDS`, `five_prime_UTR`, `three_prime_UTR`, `UTR`, case-insensitive) under their transcript by `parent`, transcripts under their gene by `parent`; keep the LONGEST transcript (by exon span) per gene; blocks `[t, s, e]` with `t` in `"E"` (exon/noncoding exon), `"C"` (CDS), `"U"` (UTR); label = gene `label`; if the file has no gene/transcript rows (BED), each feature is its own entry with a single `E` block.
    - `t` (te layer): `family` = feature `name`; `r` and `p` analogously from `ncrna` and `phas` layers.
    - Each list is sorted by start, truncated to `FEATURE_CAP` keeping highlighted ones first, and the dropped count is in `mo`.
    - `ly` lists the layer kinds that were loaded (so the renderer can say "no te layer given" versus "nothing in this window").
  - `build_windows(clusters, layers, flank: int, contig_lengths: dict) -> dict` → `{cluster name: {"w": [ws, we], "L": [start, end], **window_features(...)}}` for the clusters passed.

- [ ] **Step 1: Write failing tests:** `window_for` clamps at 0 and at the contig end, flank 0 returns the locus; an Ensembl-style fixture (gene G1, mRNA T1 and T2 of different length, exons/CDS/UTR) gives ONE gene entry using the longest transcript with `E`/`C`/`U` blocks; a BED gene layer gives one `E` block per row; a TE layer gives `family` and `hl` 1/0 correctly at the locus boundary (half-open: a TE ending exactly at the locus start is `hl == 0`); a feature longer than the window keeps its real coordinates; 100 TEs in a window give 60 entries, highlighted first, `mo["t"] == 40`; layer kinds listed in `ly`; a contig with no features in a loaded layer gives empty lists (not an error); `build_windows` keys by cluster name and uses the clamped window.
- [ ] **Step 2: Run** `python -m pytest benchmark/tests/test_genome_windows.py -v` → FAIL.
- [ ] **Step 3: Implement** `window_for`, `window_features`, `build_windows` in `report_genome.py`; pure functions, no I/O.
- [ ] **Step 4: Run** the file and the full suite → PASS.
- [ ] **Step 5: Commit** `git add report_genome.py benchmark/tests/test_genome_windows.py && git commit -m "feat(report): slice annotation features per locus window into gene models"`.

---

### Task 3: Window coverage by strand and size class

**Files:**
- Modify: `report_signals.py`
- Test: `benchmark/tests/test_genome_coverage.py`

**Interfaces:**
- Consumes: the BAM weighting helpers already in `report_signals.py`/`samio.py` (1/NH weighting, `xm_as_nh`), `_quant` quantisation (Pairs `[max, ints 0..99]`).
- Produces:
  - `SIZE_CLASSES = ("18-20", "21", "22", "23", "24", "25-35")` (reads outside 18-35 are ignored).
  - `window_tracks(fh, chrom: int|str, ws: int, we: int, *, xm_as_nh: bool, max_bins: int = 400) -> dict` → `{"s": ws, "b": bin_size, "nb": n_bins, "k": list(SIZE_CLASSES), "+": {cls: Pairs}, "-": {cls: Pairs}, "e": {"+": Pairs, "-": Pairs}}` where each entry is the mean 1/NH-weighted depth per bin of reads of that size class and strand, `e` the 5' end counts per bin, `bin_size = max(1, ceil((we - ws) / max_bins))`, `nb = ceil((we - ws) / bin_size)`.
  - `window_tracks_for(bam_path: str, windows: dict, xm_as_nh=None, threads: int = 1) -> dict` → `{cluster name: window_tracks(...)}` for `windows = {name: (chrom, ws, we)}`; same chunking/ProcessPool pattern and the same "failing region is zero-filled and counted in one WARN" behaviour as `signals_for`.

- [ ] **Step 1: Write failing tests** on a synthetic pysam BAM (build with pysam in tmp_path): two reads of 21 nt `+` and one of 24 nt `-` in known positions give the expected per-class per-strand per-bin depths and 5' end bins (for `-` reads the 5' end is the right end); NH=2 reads count 0.5; reads outside the window and reads of 16 nt are excluded; `bin_size`/`nb` for a 2,000 bp window and a 20,000 bp window (cap 400); `window_tracks_for` with `threads=2` equals `threads=1`; a window on an unknown contig is zero-filled with exactly one WARN.
- [ ] **Step 2: Run** `python -m pytest benchmark/tests/test_genome_coverage.py -v` → FAIL.
- [ ] **Step 3: Implement** `window_tracks` and `window_tracks_for` in `report_signals.py` (reuse the accumulation helpers; do not duplicate the 1/NH logic).
- [ ] **Step 4: Run** the file and the full suite → PASS.
- [ ] **Step 5: Commit** `git add report_signals.py benchmark/tests/test_genome_coverage.py && git commit -m "feat(report): window coverage by strand and size class"`.

---

### Task 4: CLI, payload and contract

**Files:**
- Modify: `srnaprofiler.py`, `report_cli.py`, `report_html.py`, `report_text.py`
- Test: `benchmark/tests/test_genome_report.py`

**Interfaces:**
- Consumes: `annotation.load_layers`, `annotation.parse_layer_spec`, `report_genome.build_windows`, `report_signals.window_tracks_for`, `samio.bam_references`.
- Produces:
  - `report --annotation LABEL:PATH` (repeatable, same labels and parsing as `classify`, validated before any BAM work, clean ERROR and non-zero exit on a bad label, a missing file, or a bad/truncated `.gz`) and `report --flank N` (`type=int`, default `2000`, must be `>= 0`, ERROR otherwise).
  - In `run_report`: when annotation layers were given, build the windows for `sel.detail` (contig lengths from the BAM header), compute `window_tracks_for` over the same windows (honouring `--threads`), and pass both to `report_html.build_report(..., genome=...)`; log the number of windows and the payload bytes added (INFO) and WARN when the added payload exceeds 20 MB.
  - `build_report(sel, summary, signals, title, out, mirna=None, genome=None)`; `genome` is `{cluster name: {"w", "L", "g", "t", "r", "p", "mo", "ly", "cv": window_tracks dict}}`. `_detail` adds key `gm` (that dict) to a cluster's detail only when present. The data-contract comment at the top of `report_html.py` documents `gm` field by field. `meta["gm"] = 1` and `meta["flank"] = N` are set only when `genome` is given.
  - New STRINGS keys under `ui.gm.*` for the figure title, track names, "nothing annotated in this window", "no <kind> layer given", "N more not shown", zoom control labels and tooltip rows (the JS tasks consume them; add them all here so the table is edited once).

- [ ] **Step 1: Write failing tests:** parser accepts repeatable `--annotation` and `--flank`, default 2000, negative flank → error; a plain `run_report` (no annotation) payload has no `gm` and no `meta.gm`; with a tiny synthetic BAM + GFF the payload detail for a locus carries `gm` with the window `w == [max(0, s-2000), min(len, e+2000)]`, `L`, gene entry and `cv`; `--flank 500` changes `w`; a bad label / missing file / fake `.gz` exits non-zero before BAM work; the contract comment mentions `gm`; byte-identity of the no-annotation payload against a stored baseline built in the test from the same fixture with the flag absent.
- [ ] **Step 2: Run** `python -m pytest benchmark/tests/test_genome_report.py -v` → FAIL.
- [ ] **Step 3: Implement** the flags, wiring, payload key and strings; validation order mirrors `classify` (labels, files and gz before the BAM is opened).
- [ ] **Step 4: Run** the file and the full suite → PASS.
- [ ] **Step 5: Commit** `git add srnaprofiler.py report_cli.py report_html.py report_text.py benchmark/tests/test_genome_report.py && git commit -m "feat(report): --annotation and --flank, genome-view payload"`.

---

### Task 5: The figure (tracks, tooltips, table)

**Files:**
- Modify: `report_assets/app.js`, `report_assets/report.css`
- Test: `benchmark/tests/test_genome_figure_js.py` (stub-DOM, same approach as `test_html_report_js.py`)

**Interfaces:**
- Consumes: `d.gm` (Task 4), the existing `figure({...})` helper, `ticks`, `niceStep`, `xTicks`, `hits`, `showTip`, `STRAND_FILL`, size colours of the size histogram.
- Produces: `figGenome(d, r)` registered as the FIRST figure on the detail page when `d.gm` exists (above the coverage figure; take the next free figure letter and keep letters consistent with the existing figures), drawing in one SVG, top to bottom:
  1. ruler with 1-based ticks (`niceStep`), the locus drawn as a shaded band across all tracks with its edges marked;
  2. coverage: plus strand above the axis, minus strand below, stacked by size class with the existing size colours, one shared y scale, y ticks;
  3. 5' ends (thin track);
  4. one track per loaded layer kind in this order: genes, TEs, structural RNA, PHAS. Genes: one lane per gene (greedy lane packing so boxes never overlap), exon blocks as boxes, CDS taller than UTR, introns as a thin line with strand chevrons, label to the right or inside when it fits; a gene or feature longer than the window is clipped at the window edge with a continuation arrow. TEs coloured by class (the part of `family` before `/`: LTR, LINE, SINE, DNA, RC, other) with a legend; highlighted features (`hl`) drawn with a stronger outline; unstranded features drawn without chevrons;
  5. for a miRNA locus with `d.mi`, a marker for the precursor (`d.mi.lo`, `d.mi.L`) on the coverage track.
  - Hover tooltip per feature (label, type/family, strand or "strand unknown", 1-based coordinates, "overlaps the locus" when `hl`) via `showTip`; per-track empty states from `ly` and the lists ("nothing annotated in this window" vs "no te layer given"); "N more not shown" from `mo`.
  - Table alternative (the figure helper's `table`): feature rows with track, label, type, strand, start, end, overlaps-locus.
  - Legend: strand fill, size classes, TE classes present.

- [ ] **Step 1: Write failing tests** (stub DOM, built from a hand-made `gm`): the figure is the first on the page when `d.gm` exists and absent otherwise; one `rect`/`path` per exon of each gene lane, no two features of one lane overlap (assert on the packed lanes); the gene longer than the window is clipped (its drawn x range equals the window range) and has a continuation marker; unstranded feature has no chevron; names containing `<script>` and `"` appear only as text nodes; "nothing annotated" vs "no te layer given" strings appear for the right tracks; `mo` produces the "more not shown" text; the table rows count equals the feature count; precursor marker only for a miRNA with `d.mi`.
- [ ] **Step 2: Run** `python -m pytest benchmark/tests/test_genome_figure_js.py -v` → FAIL.
- [ ] **Step 3: Implement** `figGenome` and its helpers in `app.js` and the few CSS rules (classes `rpt-gm-*`) using existing tokens; lane packing and clipping as small pure functions so the test can address them.
- [ ] **Step 4: Run** the file and the full suite → PASS. Then build a real report from a small synthetic BAM+GFF fixture (reuse the Task 4 fixture), and screenshot in headless Chrome: a gene-dense locus, a TE locus, an empty-window locus; read the PNGs and fix what looks wrong; list in the report what you saw.
- [ ] **Step 5: Commit** `git add report_assets/app.js report_assets/report.css benchmark/tests/test_genome_figure_js.py && git commit -m "feat(report): genomic context figure with coverage and annotation tracks"`.

---

### Task 6: Zoom, pan and gene expansion

**Files:**
- Modify: `report_assets/app.js`, `report_assets/report.css`, `report_text.py` (labels already added in Task 4; add any missing)
- Test: `benchmark/tests/test_genome_zoom_js.py`

**Interfaces:**
- Consumes: `figGenome` (Task 5).
- Produces: inside the figure, a view state `{vs, ve}` (initially the whole window), controls "+", "−", "reset", "◀", "▶" as real buttons, mouse-wheel zoom around the pointer, drag to pan, keyboard (arrow keys pan, `+`/`-` zoom, `0` reset) when the figure has focus; zoom is limited to `[40 bp, whole window]`; everything redraws from the stored data (no new fetch); the visible range is shown as text ("chr:start-end, N bp"); a "Show all transcripts" toggle is NOT built (one transcript per gene is the decided rule); a "Locus" button zooms to the locus ± 10%.

- [ ] **Step 1: Write failing tests:** zoom in halves the range around the centre and is clamped to 40 bp; zoom out never exceeds the window; pan clamps at the window edges; reset restores the window; keyboard handlers; wheel zoom keeps the pointer position fixed; the visible-range text updates; the buttons are `button` elements with aria-labels from STRINGS; redraw after zoom keeps the lane assignment stable (same gene same lane).
- [ ] **Step 2: Run** `python -m pytest benchmark/tests/test_genome_zoom_js.py -v` → FAIL.
- [ ] **Step 3: Implement** the view state, controls and handlers; throttle redraw with `requestAnimationFrame`.
- [ ] **Step 4: Run** the file and the full suite → PASS. Then in headless Chrome verify the initial view and a zoomed view (use the temporary-script technique to dispatch the zoom click) and report what you saw.
- [ ] **Step 5: Commit** `git add report_assets/app.js report_assets/report.css report_text.py benchmark/tests/test_genome_zoom_js.py && git commit -m "feat(report): zoom and pan in the genomic context figure"`.

---

### Task 7: Predicted ORFs (experimental spike)

**Files:**
- Create: `report_orfs.py`
- Modify: `srnaprofiler.py`, `report_cli.py`, `report_genome.py`, `report_assets/app.js`
- Test: `benchmark/tests/test_genome_orfs.py`

**Interfaces:**
- Produces:
  - `find_orfs(seq: str, offset: int, min_aa: int = 100) -> list[tuple[int, int, str, int]]` → `(start, end, strand, frame)` genomic 0-based half-open positions of ATG-to-stop ORFs on the six frames of `seq` (`offset` is the genomic start of `seq`); ORFs touching the sequence edge are not reported; `min_aa` counts amino acids without the stop.
  - `report --predict-orfs` (EXPERIMENTAL flag, off by default, requires `--genome FASTA` indexed with `.fai`, and `--annotation`), `--orf-min-aa N` (default 100). The ORFs are added to `gm` under `"o"` for each window, EXCLUDING any ORF that overlaps an annotated CDS in the window, capped like the other tracks.
  - A track "Predicted ORFs (experimental)" in the figure with its own legend and a visible "experimental, not annotation" note; never part of the annotation verdict, never touches class/confidence/review.
- This task is a labelled throwaway spike: after the tests pass, the controller runs it on real windows and the user decides whether to keep it. Keep the code small and removable (one module, one flag, one track).

- [ ] **Step 1: Write failing tests:** a synthetic sequence with a known ORF on `+` and one on `-` returns the right coordinates and frames; an ORF shorter than `min_aa` is dropped; edge-touching ORFs are dropped; ORFs overlapping a CDS are excluded from the payload; `--predict-orfs` without `--genome` or without `--annotation` is a clean error; default report has no `"o"` key.
- [ ] **Step 2: Run** `python -m pytest benchmark/tests/test_genome_orfs.py -v` → FAIL.
- [ ] **Step 3: Implement** `find_orfs` (plain Python, standard code table for ATG/stops only), the flags, the payload key and the track (reuse the TE track drawing with its own colour and legend).
- [ ] **Step 4: Run** the file and the full suite → PASS. Screenshot one synthetic window with the ORF track in headless Chrome and report what you saw.
- [ ] **Step 5: Commit** `git add report_orfs.py srnaprofiler.py report_cli.py report_genome.py report_assets/app.js benchmark/tests/test_genome_orfs.py && git commit -m "feat(report): experimental predicted-ORF track (--predict-orfs)"`.

---

## Decisions (mine, not discussed)

- The figure exists only when `--annotation` is given, so a plain report is untouched.
- One transcript (the longest) per gene; no transcript expansion toggle.
- Coverage is coloured by six size classes; phasing-register ticks from the spec are dropped (the existing phasing figure stays; genomic phasing ticks are not derivable from the stored data).
- Features capped at 60 per track per window.
- ORFs: ATG-to-stop, six frames, minimum 100 aa, only where no CDS is annotated; experimental flag.
- Data prerequisites (not code tasks, done by the controller on the server): rebuild `rmsk_te.bed` for Drosophila with `repClass/repName` in column 4; find and fetch a TE GFF for Arabidopsis (Araport11 or TAIR10); then regenerate the two real reports for review.
