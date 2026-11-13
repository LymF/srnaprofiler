# Genome-browser view per locus — design (draft for approval)

## Intent

The HTML report shows each locus as a stack of separate figures (coverage,
5' ends, logo, sizes, phasing, hairpin). What a reader of an sRNA paper wants to
see is the locus *in its genomic context*, the way a genome browser shows it:
reads piled on the genome, with the genes, transposons and other annotated
features sitting underneath on the same coordinates. That is how one judges
"is this 24-nt cluster a TE-derived siRNA, an exon-derived endo-siRNA, or a
piRNA cluster antisense to its TE?" at a glance, and it is the natural home for
the annotation evidence (supports / conflicts / neutral) built in
`feat/annotation-evidence`.

## Goals

- On each locus's detail page, a "Genomic context" figure: locus ± flank, one
  coordinate axis, stacked tracks, zoom and pan inside the window, hover and
  click for details. Beautiful, consistent with the report's design (same
  typography, same colour tokens; the report has no dark theme today, so none is
  added here).
- Tracks (top to bottom): ruler; sRNA coverage, plus strand up and minus strand
  down, coloured by dominant read length class (reusing the report's size
  colours); 5' ends; gene models (exon, CDS, UTR, introns drawn as thin lines
  with strand arrows); transposons coloured by class (LTR, LINE, DNA, other)
  with the family in the tooltip; structural RNA and PHAS loci when those layers
  were given; for miRNA loci the existing hairpin stays a separate figure but its
  precursor is marked on the coverage track. (Phasing-register ticks were dropped: genomic phasing positions are not derivable from the stored data; the existing phasing figure stays.)
- Annotation overlap already summarised in the locus header (chip + text) is
  drawn here: the features that produced it are the ones highlighted.
- Self-contained single HTML, offline, no external libraries or CDN (the
  report works offline today and must keep doing so).
- Works for any organism: it draws whatever annotation the user supplies.

## Non-goals

- Not IGV: no whole-genome navigation, no loading tracks in the browser, no
  editing, no BAM streaming. The window is fixed per locus at build time.
- No new analysis: the figure draws data, it does not change any class,
  confidence or review flag.
- Predicted ORFs are NOT in the first version (see "Later, spike").

## Anti-goals (would be a failure even if it technically works)

- Visual clutter: a gene-dense window (dozens of transcripts) must stay
  readable. One row per gene by default (longest transcript), expandable.
- Slow or huge reports: 3000 detail pages today, budget 50 MB. Feature data is
  sliced per window and capped; a report without annotation is unchanged.
- Drawing absence as evidence: an empty feature track means "nothing annotated
  here", and the figure says so in words when a layer was not provided at all.

## Data flow (my decisions)

- `report` gains the same repeatable `--annotation label:path` flag as
  `classify` (same parsers, same layer kinds gene|te|ncrna|phas, same contig
  rules). Without it the figure shows coverage and ends only, and the report
  is as today.
- The parser keeps feature structure (`Parent`/transcript id, type, strand,
  name, family/class, id) so exons and CDS group under their transcript; the
  current `Feature` tuple gets the extra fields.
- For each detailed locus the report slices the features overlapping the window
  from the interval index and embeds them compactly under a short key in the
  detail dict, with a hard cap per window and a "N more not shown" row.
- Coverage and ends reuse the existing `cv`/`en` bins; the window is widened to
  locus ± flank (default 2 kb, `--flank`, clamped to the contig), so
  those signals are computed over the window, not only the locus. The extra BAM
  work is the main runtime cost and must be measured.
- Rendering is inline SVG built in `app.js` with the existing text-node-only
  rule (annotation names come from user files). A "Table" toggle gives the
  accessible alternative, as the other figures have.

## Verification

The report's JS has so far been tested only in a stub DOM. Acceptance for this
feature includes real-browser checks: headless Chrome screenshots,
of at least: a TE-dense Drosophila piRNA locus, a gene-dense Arabidopsis
locus, a miRNA locus, and a locus with no annotation at all, reviewed by a
person.

## Prerequisites on the data side (found on the servers, 2026-11-13)

- Drosophila: `rmsk_te.bed` has 3 columns (no repeat family). Rebuild it with
  `repClass/repName` in column 4 (`benchmark/scripts/build_te.py`).
- Arabidopsis: no TE annotation on the server (Ensembl TAIR10.63 has none). A
  TE GFF (Araport11 or TAIR10 TE) has to be found and fetched.

## Later, spike (not in the first version)

Predicted ORFs: a six-frame scan on the window with a minimum length, shown as a
track only where no CDS is annotated. The user wants to see it in practice
before deciding; it is built as a labelled throwaway spike on real windows
after the annotated version works, to judge the noise, and kept only if the
user says so.

## Left to the builder

Exact tracks heights, colours (from the existing tokens), tooltip wording,
pan/zoom interaction details, feature cap numbers, compact encoding of the
feature slices.

## Settled with the user (2026-11-13)

1. Flank: default 2 kb on each side, adjustable at build time with `report --flank N`
   (bp, clamped to the contig). Replaces the earlier `max(1 kb, locus length)` guess.
2. Scope: annotated features first, then the ORF spike; the user said to do everything.
