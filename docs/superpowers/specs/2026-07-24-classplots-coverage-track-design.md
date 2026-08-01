# classplots — per-locus coverage track (addendum)

**Author:** Lucas
**Date:** 2026-07-24
**Status:** design (approved)

Addendum to `docs/superpowers/specs/2026-07-23-classplots-design.md`. Replaces the
per-class normalized coverage metaplot with a per-locus absolute coverage track in the
style of `extras/genome-coverage.py`.

---

## 1. Goal

Give each **chosen locus** (the `--top-n` most abundant per class, and each `--region`) a
genome-browser-style coverage track: sense coverage drawn upward, antisense downward, over
the locus's own coordinates. This is the aesthetic of `extras/genome-coverage.py`, reused
from the reads classplots already fetches — no new dependency.

This **replaces** the previous per-class normalized coverage metaplot (`coverage_profile`
+ `plotCoverageProfile.R`), which is removed from the flow.

## 2. Scope decisions (from brainstorming)

- **Per chosen locus only.** The track is produced only for single-locus units
  (`len(members) == 1`): the top-N loci and each `--region`. The pooled-per-class unit no
  longer emits a coverage plot (absolute genomic position summed across loci on different
  chromosomes is meaningless); it keeps sizedist, pattern, ping-pong, phasing, and logo.
- **Raw depth (IGV style), NOT 1/NH.** Each alignment contributes +1 to every position it
  covers, split by strand — matching `genome-coverage.py` and genome browsers. This is the
  one intentional exception to classplots' 1/NH weighting (like the WebLogo distinct-seq
  exception); it must be documented in code.
- **Matplotlib, no bedtools/samtools subprocess.** Computed in Python from the already-
  fetched `samio.Aln` list; rendered with matplotlib (PDF backend is native — no Ghostscript).

## 3. Components

- `class_signals.coverage_track(alns, region_start, region_end) -> (list[int], list[int])`
  — per-base raw depth over `[region_start, region_end)` as `(sense, antisense)`, each of
  length `region_end - region_start`. For each alignment, `+1` to every covered position in
  the region on its strand lane (`+` → sense, `-` → antisense). No weighting.
- A matplotlib plotter (new small module `class_coverage.py`, function
  `plot_coverage_track(sense, antisense, chrom, start, out_path)`), reproducing
  `genome-coverage.py`'s look:
  - `fill_between(x, sense, color="#996631")` (label "Sense (+)"),
  - `fill_between(x, [-a for a in antisense], color="#0099CC")` (label "Antisense (-)"),
  - symmetric `ylim(-m, m)` with `m = 1.1 * max(max(sense), max(antisense), 1)`,
  - `axhline(0, color="black")`, dashed grid, legend,
  - x-axis = genomic position (`start + i`), title `Coverage — {chrom}:{start}-{end}`,
  - `matplotlib.use("Agg")`, save to PDF, `plt.close(fig)`.
- Orchestrator wiring in `_signal_plots_for`: **replace** the current coverage block. The
  function takes an explicit `coverage=False` flag; the per-locus (top-N) and `--region`
  calls pass `coverage=True`, the pooled-per-class call leaves it `False`. Coverage is
  emitted only when `coverage and len(members) == 1` — the explicit flag (not just the
  member count) is needed because a class with a single locus would otherwise get a coverage
  plot on its pooled call too. It computes `coverage_track` and calls `plot_coverage_track`,
  writing `{prefix}.{name}.coverage.pdf`.
  - `_signal_plots_for` needs the locus `chrom` for the title/track; pass the chrom for
    single-member units (thread it through, or derive it — the members already carry
    `(alns, start, end)`; add `chrom` so members become `(alns, chrom, start, end)` or pass
    a separate `chrom` argument used only for the single-member case).

## 4. Removals

- `class_signals.coverage_profile` and its test file are removed.
- The `plotCoverageProfile.R` invocation is removed from `_signal_plots_for`; the `.R` file
  stays in the repo (unused) to avoid touching unrelated tooling.
- The `_write_coverage_tsv` helper is removed (the new track writes no intermediate TSV, or
  writes a simple `pos<TAB>sense<TAB>antisense` table under `tables/` for reproducibility —
  keep a table for parity with the other plots).

## 5. Testing (TDD)

- `coverage_track`: synthetic `Aln` list — a `+` read over part of the region raises the
  sense lane at the covered positions only (antisense stays 0); a `-` read raises antisense;
  reads clipped to the region; two overlapping reads stack (depth 2). Raw counts (no 1/NH):
  a read with `nh=2` still contributes +1.
- `plot_coverage_track`: renders headless (Agg) to a temp PDF that exists and is non-empty.
- Orchestrator: a single-member unit writes `*.coverage.pdf`; a multi-member (pooled) unit
  writes none. (Extend the existing skipif-guarded signals smoke test.)

## 6. Risks

- Very wide loci (tens of kb) produce long per-base arrays; fine for typical sRNA clusters
  (hundreds–few kb). No downsampling in v1.
- Raw depth inflates repetitive piRNA loci (96% multimappers) visually — intended (IGV-like);
  documented as the deliberate non-1/NH exception.
