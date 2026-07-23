# classplots — per-class plotting for the classification engine

**Author:** Lucas
**Date:** 2026-07-23
**Status:** design (pending review)

Related: `docs/superpowers/specs/2026-07-22-classification-engine-design.md` (classify v1),
memory `srnaprofiler-pipeline-reference` (all signal formulas), ROADMAP item #8.

---

## 1. Goal

Give the classification engine a per-class visualization layer. After `classify`
produces a `.clusters.tsv` (one locus per row, with its assigned `cls`), `classplots`
re-aggregates the alignments **grouped by class** and produces publication-grade signal
plots plus a QC summary. This turns the scalar feature matrix and the raw signal shapes
into figures that let a human confirm the biology of each class (e.g. 1U + ping-pong in
piRNA loci, balanced strands + 21nt in siRNA loci).

### Non-goals (YAGNI)
- No re-classification and no changes to `classify`, `records`, `decider`, or `rules.yaml`.
- No per-locus small multiples beyond the `--top-n` most abundant loci per class.
- No new signals beyond what the engine already computes; scalars stay in the QC summary,
  not in bespoke R scripts.

---

## 2. CLI

New standalone subcommand (decoupled from `classify`, re-runnable without re-classifying):

```
srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam -o OUTDIR
    [--prefix NAME] [--top-n 5] [--classes siRNA,piRNA,miRNA]
    [--min-length 15 --max-length 35] [--no-signals] [--no-summary] [--threads N]
```

- `--clusters` / `--bam` required. BAM must be sorted+indexed (same input `classify` uses).
- `xm_as_nh` is auto-detected from the BAM `@PG` (`samio.prefers_xm_as_nh`), not a flag.
- Classes plotted = the distinct values of the `cls` column (siRNA / piRNA / miRNA).
  `--classes` only **restricts** that set; it never invents classes.
- `--top-n 0` disables the per-locus plots (pooled-only).
- `--no-signals` / `--no-summary` toggle the two layers independently.

---

## 3. Fidelity to the engine (core principle)

The legacy `plots.py` reads raw SAM, dedups by `read_id`, and does **no 1/NH weighting**.
`classplots` must instead read through `samio` so the plots reflect exactly what the
classifier saw:

- `samio.read_region` (indexed `fetch`), `-F 256` (primary only),
- weight every read by `w = 1/NH` (NH via `XM` for bowtie1, per `samio._nh`),
- reverse-complement reads on the `-` strand before reading first base / building logos,
  so composition is measured in the RNA 5'→3' orientation.

All per-class accumulators consume `samio.Aln` tuples and carry these conventions.

---

## 4. Architecture

New module `classplots.py` with `run_classplots(args)` as the orchestrator, plus small
**pure, testable** units. It reuses the R-invocation half of `plots.py`
(`run_r_script`, `check_r_dependencies`, `get_r_scripts_path`, `generate_*`) and the
matplotlib logo core of `weblogo_new.py`.

### Units (pure functions)
- `load_clusters(tsv) -> list[Locus]` — parse the `.clusters.tsv`; keep coords, `name`,
  `cls`, `reads`, and the scalar columns needed by the QC summary.
- `group_by_class(loci, classes) -> dict[str, list[Locus]]`.
- `top_n_loci(loci, n) -> list[Locus]` — by `reads` descending.
- **Signal accumulators** (input: iterable of `samio.Aln`; output: matrices in the exact
  format each downstream plotter expects; all 1/NH-weighted, `-` revcomp'd):
  - `base_composition(alns, minlen, maxlen)` → `{size: {+/-: {A,C,G,U}}}` (float counts).
  - `pingpong_offsets(alns)` → `{overlap: weighted_pairs}` for 1..30.
  - `phasing_series(alns, phase_length)` → per-strand periodicity series.
  - `pattern_matrix(loci_alns)` → `#`-separated matrix, **one row per locus**, cols 0–20 =
    `+` strand size counts (15–35nt) z-scored within the row, 21–41 = `-` strand. Ported
    from `clusters.py::_generate_zscore_matrix` (originally one row per chromosome; the
    per-class pooling makes the row = locus, so `plotPattern`'s error bars become
    locus-to-locus variability within the class).
  - `positional_composition(alns, length)` → position→{A,C,G,U} frequencies for the logo.
  - `coverage_profile(loci_alns)` → per-strand normalized coverage across loci (each locus
    rescaled to [0,1] and averaged).

### Targeted refactor (serves reuse, not scope creep)
Extract the `plus5/minus5 → overlap_counts` accumulation currently inlined in
`features.build_record` into `signals.py` (e.g. `signals.overlap_counts(plus5, minus5)`),
so `classify` and `classplots` share **one** source of truth for the ping-pong signal.
`features` imports it; behavior unchanged (covered by existing tests).

---

## 5. Plot catalog

### Signal plots — per class (pooled) and per top-N locus

| Plot | Backend | Script | Status |
|------|---------|--------|--------|
| Size distribution + first base (+ sum-strand, + piRNA zoom) | R | `calcDistributionPerBase*.R` | reuse |
| Size×strand z-score profile (error bars across loci) | R | `plotPattern_publication.R` | reuse (generator ported from `clusters.py`) |
| Ping-pong signature (overlap 1–30, highlight 10, annotate `pingpong_frac`/z) | R | `plotPingPong.R` | **new** |
| Phasing periodicity (per strand, period = `DicerCall`, overhang −2) | R | `plotPhasing.R` | **new** |
| Sequence logo (positional A/C/G/U) | Python/matplotlib | reuse `weblogo_new.py` core (self-contained TextPath; **not** `weblogo.py`, which shells out to the external `weblogo` CLI that is unavailable) | **new wiring** |
| Coverage / strand profile (normalized metaplot per class) | R | `plotCoverageProfile.R` | **new** |

The legacy `plotOffsetPhasingSeries.R` is **retired** and replaced by the two correct,
separated scripts (`plotPingPong.R`, `plotPhasing.R`). It bundled three panels and used an
ad-hoc `scale(Count + 1e-9)` z-score that does not match the engine's conventions.

### QC summary — from `.clusters.tsv` scalars (Python/matplotlib)
One multi-panel figure set, works even when `Rscript` is absent:
- loci count per class;
- `dicer_call` distribution per class;
- box/violin of `pingpong_frac`, `phase_score`, `pct_1u`, `gc`, `multimap_frac` per class;
- `confidence` distribution and `review` rate per class.

---

## 6. Phasing / ping-pong correctness (explicit verification task)

The old phasing/offset R was written before the engine's conventions existed, so it is
**not** assumed correct. Before either new script is accepted:

1. Define the per-class series to match `signals.py` / `phasing.py` exactly:
   - ping-pong: 5'–5' cross-strand overlap, canonical peak at **offset 10**;
   - phasing: register `(pos5 − region_start) mod L` with `L = int(DicerCall)`, and the
     **−2 overhang correction on the `-` strand** (see `srnaprofiler-pipeline-reference §5`).
2. Prove alignment with synthetic `Aln` fixtures: a constructed pair whose 5' ends overlap
   by exactly 10 must land in the offset-10 bin; a constructed phased ladder at period `L`
   must peak at its known register. Tests assert the accumulator output before any R runs.
3. If the R rendering (axis, z definition) disagrees with the engine, fix the R — do not
   bend the engine to the old plot.

---

## 7. Output layout

```
OUTDIR/tables/   {prefix}.{class}.sizedist.tsv        (size/±/ACGU, weighted)
                 {prefix}.{class}.pattern.matrix       (# 43-col z-score)
                 {prefix}.{class}.pingpong.tsv         (overlap, pairs)
                 {prefix}.{class}.phasing.{plus,minus}.tsv
                 {prefix}.{class}.coverage.tsv
OUTDIR/plots/    {prefix}.{class}.sizedist.pdf (+ _sumstrand, + _zoom for piRNA)
                 {prefix}.{class}.pattern.pdf
                 {prefix}.{class}.pingpong.pdf
                 {prefix}.{class}.phasing.pdf
                 {prefix}.{class}.logo.pdf
                 {prefix}.{class}.coverage.pdf
                 {prefix}.{class}.locus.{name}.*.pdf   (top-N)
                 {prefix}.summary.{counts,dicer,signals,review}.pdf
```

Intermediate `tables/` are kept for reproducibility and re-plotting outside the tool.

---

## 8. Error handling
- Class with no loci → warn and skip.
- `Rscript` absent (`check_r_dependencies` false) → warn, skip **signal** plots, still emit
  the matplotlib QC summary and the sequence logo (both pure Python).
- Empty region / locus with zero kept reads → guard, skip that accumulation.
- Malformed `.clusters.tsv` row → skip row with a warning; missing required columns → error.

---

## 9. Testing (TDD)
- Accumulators on synthetic `Aln` lists: 1/NH weighting sums correctly; `-` strand
  revcomp for first base / logo; ping-pong pair at offset 10; phased ladder register;
  coverage normalization to [0,1].
- `load_clusters` parses a fixture TSV (and rejects a missing-column one).
- `group_by_class` / `top_n_loci` ordering.
- R invocation: smoke test that **skips** when `Rscript` is unavailable (mirrors
  `check_r_dependencies`), otherwise checks a PDF is produced from a tiny fixture.
- QC summary and logo render headless (matplotlib `Agg`) without R.

---

## 10. Risks / open items
- **Resolved:** the `pattern` matrix generator is `clusters.py::_generate_zscore_matrix`
  (Python, from Perl `calcPatternInSamFile.pl`) — port it, swapping the dedup'd SAM read
  for weighted per-locus `Aln`. The logo is `weblogo_new.py` (matplotlib); the old
  `weblogo.py` (external `weblogo` CLI) is the one that errored and is avoided.
- Coverage profile could partially echo `frac_top` in the QC summary; kept because the
  **positional shape** (edge vs center, strand asymmetry) is not captured by the scalar.
- Weighted (float) counts flow into R scripts written for integer counts; verify ggplot
  handles floats (expected yes) during the port.
- Both legacy generators dedup by read and are unweighted; the port must instead use
  1/NH-weighted `Aln` (Section 3) so the plots stay consistent with the classifier.
