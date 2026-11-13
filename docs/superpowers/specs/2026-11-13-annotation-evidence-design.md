# Annotation as evidence — design

## Intent

`classify` calls loci from small-RNA signatures only (size, 5' base, ping-pong,
phasing, strand bias). It never looks at what the locus overlaps in the genome.
A reader (and a reviewer) wants to know: is this 24-nt cluster on an LTR
retrotransposon? Is this piRNA cluster antisense to its TE? Is this "siRNA" a
tRNA fragment? That context is both a biological sanity check and the natural
way to describe the results.

Add an optional, organism-agnostic annotation layer that **describes** every
locus and, where biology supports it, **supports or contradicts** the class
call. It never decides a class on its own.

## Goals

- Accept optional annotation (GFF3, GTF, BED; plain or gzip), several layers at
  once, each with a label: genes, TEs/repeats, structural RNA, PHAS/other.
- For each locus: overlap per layer (feature, class/family, fraction of the
  locus covered, strand relation sense/antisense/both, feature ID).
- Write it to `clusters.tsv` (new columns), `clusters.gff3` (attributes) and the
  HTML report (a per-locus annotation line, and a per-class summary).
- An `annotation_support` flag per locus: `supports`, `conflicts`, `neutral`,
  `none` (no annotation given), derived from a small biology table in the
  config (so plant and arthropod can differ).
- `conflicts` and `supports` feed a **review flag and a bounded confidence
  nudge**, never a class change.
- Optional helper to annotate only the loci (not the genome) with RepeatMasker
  or a TE library, producing a GFF the layer above can read.

## Non-goals

- Running RepeatModeler/EDTA/de novo TE discovery. Documented, not wrapped.
- Annotating the whole genome with RepeatMasker inside the pipeline.
- Changing the class from annotation. Not in the ML features either: the ML
  model stays trained without it.
- Gene-level analysis (target prediction, GO).

## Anti-goals

- Annotation must not leak into the benchmark. Benchmark runs without
  `--annotation` must be byte-identical to today.
- Absence of an annotated TE must never count as evidence against a TE class
  (annotations are incomplete).
- No hard dependency: no annotation given means no new columns beyond
  `annotation_support = none`, and no new required tools.

## Biology table (default, editable in the rules YAML)

| Class | Supports | Conflicts |
|---|---|---|
| piRNA | TE/repeat overlap; sense/antisense to TE consistent with ping-pong | exon/CDS of a protein-coding gene with no TE |
| siRNA 24-nt / hc-siRNA (plant) | TE/repeat | CDS only |
| siRNA 21/22-nt | TE (easiRNA), gene (endo-siRNA) | none |
| phasiRNA | PHAS layer, or gene of NLR/PPR/TAS family | none |
| miRNA | intergenic or intron | TE (flag only, many MIRNAs derive from MITEs) |
| tRNA/rRNA/snoRNA fragments | structural RNA layer | TE |

A locus overlapping nothing is `neutral`. The table maps (class, feature kind,
strand relation) to supports/conflicts/neutral; the nudge size is a config
value, default small and capped.

## Constraints

- Reuse `pybedtools` (already a dependency) or a pure interval index; no new
  heavy dependency.
- Large genomes: index by chromosome, sorted/binary search, memory bounded.
- Contig names must match the BAM; mismatches reported once, with how many
  features were skipped.
- GFF `Parent`/`Name`/`Class`/`family` attributes read generically. RepeatMasker
  `.out` and EDTA `.TEanno.gff3` supported as named formats because they are
  the common TE inputs.

## Interface (proposal, builder may refine)

```
srnaprofiler classify ... \
  --annotation genes:genome.gff3 \
  --annotation TE:repeatmasker.out \
  --annotation ncRNA:rfam.gff
```

- `label:path` repeatable; label becomes the layer name and column suffix.
- A separate subcommand `annotate-loci` takes `clusters.tsv` + genome FASTA and
  runs RepeatMasker (or a user TE FASTA via minimap2/BLAST) on the locus
  sequences only, writing a GFF in genome coordinates. Minutes, not hours.
- `report` reads the new columns when present.

## Decided by me (not discussed, please check)

- Evidence-only, bounded nudge, no class change.
- Annotation excluded from ML features and from the benchmark.
- Locus-only RepeatMasker as the answer to the runtime concern.

## Left to the builder

- Exact column names and nudge formula.
- Interval index implementation.
- How the report shows it (follow the existing report patterns).
- Test fixtures: tiny GFF/RepeatMasker .out/EDTA files and a toy BAM.

## Settled with the user (2026-11-13)

1. Small, capped confidence nudge: approved.
2. Which real annotations exist for the paper's organisms is unknown; it is a
   separate search task. Implementation uses synthetic fixtures and supports all
   listed formats generically, so it does not wait on that search.
