# HTML report v2: figure-first cluster pages and an About page

Date: 2026-11-12. Status: draft for review. Amends `2026-11-11-html-report-design.md` (v1, built on
branch `feat/html-report`); where this document and v1 disagree, this one wins.

## Intent

v1 explains every check in sentences, which the user finds too much hand-holding and not informative
cluster by cluster. v2 turns the cluster page into a figure page that shows the data the way the
lab's own R scripts do, keeps text to a one-line verdict and a strip of measured values, drops the
glossary, and adds an About page like ViralQuest's.

## Goals

- **Cluster page figures** (all with sense (+) up and antisense (-) mirrored down, as in the R plots
  in the repo):
  1. Read size x strand x 5' base: stacked bars A green `#108041`, C blue `#3B54A2`, G orange
     `#F9A31B`, U red `#EA2627`, x axis "small RNA length (nt)" (`calcDistributionPerBase_publication.R`;
     numbers from `class_signals.base_composition`).
  2. Phasing: 5'-5' distance (1..100 nt) as a line, one panel for sense and one for antisense, red
     dashed lines at multiples of the period (class period: 21 or 24; otherwise the cluster's dominant
     size) (`plotPhasing.R`; numbers from `class_signals.phasing_series`).
  3. Coverage over the whole cluster, sense up and antisense down, full extent of the cluster
     (`class_coverage.py` look).
  4. Ping-pong 5'-5' overlap 1..30 nt with the 10 bar in red `#e41a1c`, others grey `grey60`
     (`plotPingPong.R`).
  5. 5' end track: positions of read 5' ends along the cluster, per strand.
  6. Nucleotide logo per position for the dominant size, split by strand.
  7. Table of the most abundant distinct sequences (sequence, length, strand, reads, 5' base).
  8. Library context: where the cluster sits in the library's distribution of reads and, within its
     class, of confidence (percentile).
- Text on the cluster page: one-line verdict and a compact strip of measured values (size, 1U, 10A,
  ping-pong z, phasing p, strand, multimapping) with values and the colours/context markers v1 added;
  no "how to read" sentences, no per-check paragraphs. Short hover tips are allowed.
- **About page**, first tab, like ViralQuest's `About`: hero with title "sRNAProfiler", the tagline
  taken from the README ("Small RNA loci discovery and classification from sequencing data"), the
  version (the string behind `srnaprofiler --version`), a GitHub link card
  (`https://github.com/LymF/srnaprofiler`, "Source code & documentation") and an install snippet
  (`conda env create -f environment.yml` / `conda activate srnaprofiler`). The About CSS already
  exists in `report_assets/report.css` (`rpt-about*`). Everything the project does not have yet is
  left out, not faked: no paper card, no PyPI card, no logo, no citation block. They are added when
  they exist.
- Detail pages limited: default `--max-detail` **3000** (revised 2026-11-12: the user raised the size budget to about 50 MB, so the v1 default stays); the table keeps every cluster above the floor (300 reads by default). Larger values stay available; the report warns above 5000 pages (browser memory).
- The report stays one offline file, opt-in `srnaprofiler report`, classify untouched.

## Non-goals

- No change to clustering, classification or the ML classifier.
- No LLM / AI text (still phase 2).
- No multi-library page.

## Removed

- Glossary tab and References tab, and every `GLOSSARY` / `REFERENCES` entry and `cite` plumbing in
  `report_text.py`, `report_html.py`, `app.js` and tests. (User: "acho desnecessario".) The README
  section is updated. The per-fact citation rule applied to the glossary only; the About page makes
  no biological claims.
- The sentence checklist, the "how to read" lines and the "why not another class" paragraphs on the
  cluster page. The verdict line and the measured-value strip stay.

## Constraints

- Data size and speed are part of this work (previously deferred as findings 7 and 8):
  - Richer figures make each page heavier. Target: a default report (1000 detail pages) of about
    50 MB or less (user: "ate uns 50 da joia"; about 21 MB measured for 3000 dense synthetic pages) on the Drosophila library used before, with the default 300-read floor.
  - Quantise per-cluster signals compactly (for example integer 0-99 scaled to a stored per-cluster
    maximum, trailing zeros dropped) and send keys and values rather than rendered sentences; the
    page formats text from the strings table.
  - Generation time: the first v1 run took ~13 minutes for 3,000 pages. v2 must read each cluster's
    reads once, aggregate in a single pass (counters keyed by position, strand, length), build
    coverage with a difference array instead of a per-base loop, and not hold every alignment of a
    deep locus in memory. Optional parallelism over clusters (each worker opens its own BAM handle)
    with a `--threads` flag, default 1.
- Strings stay in the strings table (translation-ready); no `innerHTML` with data; every figure is
  inline SVG; no external resources.
- Sense/antisense separation in every figure where strand is defined; the 5' base colours are the
  ones above.
- Code and comments in English, short comments. Authorship Lucas only, no AI credit.

## Parts left to the builder

- Resolution of the whole-cluster coverage (bin count vs cluster length), bin rule, and the
  quantisation scheme; the figure sizes and layout of the cluster page.
- How the "library context" percentiles are computed and stored (one pass over the clusters in
  `report_data`, shared across pages).
- Whether the sequences table is capped at 10 or 20 rows (builder's call, default 10).
- Where the version string comes from so that `srnaprofiler --version` and the About page cannot
  drift.
- How the overview and list tabs adjust once the glossary is gone (no behaviour change beyond that).

## Decisions I made (the user can overrule)

- Default `--max-detail` 3000 (revised from 1000 when the user raised the size budget to about 50 MB).
- The About page shows only what exists; no placeholders for the paper, PyPI or a logo.
- The 5' end track and the nucleotide logo use the strand-split convention of the other figures.
- Hover tips stay, one short line each.
- Library context is percentile-based and shown as a small strip, not a full chart.

## Open questions

None. (Hover tips: kept short; decided by the builder-controller, the user said to proceed.)
