# HTML report for srnaprofiler: design

Date: 2026-11-11. Status: draft for review.

## Intent

A person who is new to small RNAs runs srnaprofiler and gets a single HTML file that tells them,
in plain language, what is in their library, which clusters deserve attention, and why each
cluster got its class. The same file must serve three uses: understanding the library, picking
clusters to validate at the bench, and writing (figures and text a user can cite).

The user explicitly said the reader is a layperson, so the report explains the biology; it does not
only show numbers.

## Goals

- One self-contained, offline HTML file per library (CSS, JS and data inline), in the visual style of
  the `report-template` skill (`~/.claude/skills/report-template`, the ViralQuest look).
- Three layers: library overview, cluster list, single-cluster page. Plus a glossary tab.
- Every class call is explained with the evidence behind it and the reason the other classes were
  rejected. An `ambiguous` call says what was missing.
- A depth floor of **300 reads** by default (user decision), user-adjustable, applied only in the
  report.
- Texts come from rules in this first version. They are written so that a later model can replace or
  extend them (phase 2).

## Non-goals

- No change to clustering, classification, the ML classifier or the benchmark. The pipeline's own
  locus floor (`mincov`, `fragments.tsv`) stays as it is; the report floor is a view filter.
- No LLM in this version. No "validate first" badge (user decision): ordering by confidence covers it.
- No multi-library comparison page.

## Anti-goals (failure even if it technically works)

- Reporting a confidence number without saying what it means.
- Silently dropping clusters: the report always states how many clusters were hidden by the floor and
  how many passed it.
- Presenting evidence thresholds as biological truth when they are only display thresholds.

## Constraints

- Inputs are what `classify` already writes: `<prefix>.clusters.tsv`, `.clusters.gff3`, run JSON.
  Per-cluster figures need reads, so they are computed from the BAM with the existing
  `class_signals` helpers (base composition, ping-pong offsets, phasing, coverage track) and stored as
  small JSON arrays, not PNGs.
- Size: with the floor at 300 the example library keeps 1,625 of 25,738 clusters (1,267 piRNA, 190
  ambiguous, 100 siRNA, 68 miRNA). The file must stay in the low single-digit MB. Libraries where very
  few clusters pass the floor must show a warning and a one-click way to lower it.
- Offline: external scripts or fonts are not allowed in the output file.
- Code and comments in English; report text in English first, with strings kept in one table so a
  Portuguese version is a translation, not a rewrite. Authorship: Lucas only.
- A new subcommand `report` (own code path), so existing subcommands and outputs are unchanged.
  Opt-in: nothing runs unless the user asks.

## What the user decided

- Audience: layperson; uses = understand, choose clusters, write.
- Three layers, all important; glossary included.
- Depth floor 300 reads (default), no "validate first" flag.
- Spike first, then expand to a model (phase 2).
- Cluster pages for every cluster above the floor up to 3,000; above that, the 3,000 with highest
  confidence, with a visible notice (user agreed 2026-11-11).
- Separate subcommand `srnaprofiler report`, run after `classify` (user agreed 2026-11-11).
- Depth floor unit: raw reads by default (builder's reading of "agreed with all three"; the user
  was shown the RPM trade-off and did not pick it). RPM stays a possible option, to be confirmed.

## Content

**Overview.** One sentence summary built from the class shares; key numbers (clusters analysed,
hidden by floor, ambiguous, reads in analysed clusters); class donut; dominant-size histogram
stacked by class; library cards (library size and unique-mapping share derived from the clusters TSV, reads in clusters,
fragments below the pipeline floor; the mapped percentage is not in the classify outputs, so v1 omits it) with automatic warnings; each chart
has a one-line "how to read".

**Cluster list.** Filter by class; sort by confidence, reads, name; pagination; per row: location,
class tag, confidence bar, five evidence dots (size, 1U, 10A, ping-pong, phasing), reads; for
ambiguous rows, what was missing. Click opens the cluster page.

**Cluster page.** Verdict (class, confidence, one sentence); checklist of evidence with a plain
explanation each (size, 1U, 10A, strand, ping-pong, phasing, multimapping); four figures (read-size
histogram, nucleotide per position, coverage by strand, ping-pong overlap histogram); table "why not
another class"; cautions (for example high multimapping); dominant sequence.

**Glossary.** About a dozen terms. Each definition that states a biological fact carries a citation,
and the report ends with a References section with DOIs.

## Parts left to the builder

- Mapping from rule/gate columns to the explanation sentences. The `e_*` and `gate_*` columns exist in
  `--classifier ml` runs and in current rules runs; the sample `results-teste` predates them and
  cannot be used to test the explanations.
- Colour thresholds for the evidence dots: either reuse the rules' own thresholds from
  `rules.yaml` / `rules.plant.yaml` or document them as display-only. They must not be invented
  silently (the spike's thresholds were illustrative guesses).
- Plant classes (phasiRNA, hcsiRNA) need explanations written, not only the animal ones.
- File layout, test strategy, and the generator's internal structure.

## Phase 2 (not part of this build)

Add an AI-written paragraph per library and per cluster, always labelled as a suggestion.
ViralQuest does this with a local Ollama model selected by `--model-type ollama --model-name <model>`;
the showcase page shows an off-the-shelf model there, and it does not say the model was trained.
The user's idea (2026-11-11) is to base the model on the literature (papers) rather than on
benchmark labels; the approach (for example retrieval over papers versus training) is deferred. The report's
text layer should expose one function per paragraph so it can be swapped later.

## Open questions

1. Confirm raw reads (default) versus RPM for the depth floor.

## Spike evidence

A throwaway spike (not committed) rendered the real `dros_ref` clusters with the floor adjustable.
At 100 reads 7,355 clusters remain; at 300, 1,625. The spike has no per-cluster figures because the
TSV alone cannot provide them. The spike's JavaScript was not exercised in a browser by the author.
