# sRNAProfiler

Small RNA loci discovery and classification from sequencing data. It calls loci from
coverage, measures the signatures each biogenesis pathway leaves behind (size, strand
balance, 5′ nucleotide, ping-pong overlap, phasing), and assigns a class with a transparent,
editable rule set. Plants and animals share one engine, and the biology lives in a config file.

Classes come out of `rules.yaml` (animal: siRNA, piRNA) or `rules.plant.yaml` (plant:
phasiRNA, hcsiRNA, siRNA, and optional 24-nt phasiRNA). miRNA is not called by signature; it
is taken from annotation or prediction and carved out of the surrounding locus.

Two optional layers sit on top of the signature call and never change a class:

- **Annotation as evidence** (`classify --annotation`): genes, transposons, structural RNA and
  PHAS loci overlapping each locus support or contradict the call, nudge its confidence by a small
  capped amount and flag conflicts for review (see [Annotation](#annotation-optional)).
- **A self-contained HTML report** (`report`): one offline page with a plain-language reading of every
  locus, the miRNA hairpin, and, when an annotation is given, a genome-browser-style figure of the
  locus in its context (see [Report](#report)).

### Status

| part | state |
|---|---|
| signature classifier (rules), plant and animal | stable; the benchmark runs on it |
| `--classifier ml` | optional, needs a trained model |
| `mirna` (miRDeep2 / miRDeep-P2) | stable |
| `report` | stable; checked in headless Chrome on real Drosophila and Arabidopsis libraries |
| `classify --annotation` | stable; without it the outputs are byte-identical to a run from before the layer existed, so benchmark numbers stay valid |
| genome figure in `report` | stable; zoom, pan, table and SVG export |
| `report --predict-orfs` | **experimental**: ORFs from the genome sequence, never used by the classifier |
| `annotate-loci` (RepeatMasker on the loci only) | implemented and unit-tested; see the note in [Annotation](#annotation-optional) before relying on it |

Real-data check of the annotation layer (class calls identical with and without it):
a *Drosophila melanogaster* library with the Ensembl dm6 genes and the RepeatMasker track
(10750 loci: 3488 piRNA `supports`, 282 `conflicts`, 698 `neutral`), and an *Arabidopsis thaliana*
library with Ensembl TAIR10.63 genes and the TAIR10 TE annotation (hcsiRNA: 12601 `supports`, 1468
`conflicts`, 6878 `neutral`; without a TE layer every hcsiRNA is `neutral`).

---

## Install

```bash
conda env create -f environment.yml
conda activate srnaprofiler
```

Two tools conda cannot install are handled by the `setup` module:

```bash
python3 srnaprofiler.py setup --check              # what is present, and which modules can run
python3 srnaprofiler.py setup --tools              # miRDeep-P2 (SourceForge tarball + a one-line patch)
python3 srnaprofiler.py setup --data DIR           # Rfam index for the ncRNA filter, ~2 h to build
```

`setup` records the paths it prepared in `~/.config/srnaprofiler/paths.yaml`, so `mapping` and
`mirna` find them without flags.

---

## Order of the steps

```
qc -g genome.fa  →  mapping  →  mirna  →  classify --mirna-loci <BED>  →  classplots
```

`qc -g` is **not optional**: it standardizes FASTA headers and maps IUPAC codes to N. The
miRNA engines reject a genome with either, and NCBI genomes have both.

`classify` has no miRNA class of its own: miRNA loci come from the BED passed to
`--mirna-loci`. A precursor that is not in that BED is classified by its signature like any
other locus, and since miRNAs share their size range and 5' bias with other pathways, it may
be assigned to a neighbouring class. Passing a curated annotation, or the validated
predictions from `mirna`, is what keeps them labelled correctly.

---

## Modules

| module | use it to |
|---|---|
| `setup` | check the toolchain and prepare reference data, once per machine |
| `qc` | trim reads and standardize the genome |
| `mapping` | align to the genome, optionally dropping structural RNA |
| `mirna` | find and validate miRNA loci, so `classify` can carve them out |
| `classify` | call loci, measure their signatures, assign a class (optionally with `--annotation`) |
| `annotate-loci` | annotate transposons on the classified loci only (RepeatMasker), for organisms without a TE annotation |
| `te-library` | build a TE library for a genome (RepeatModeler2) and drop gene contaminants, to feed `annotate-loci --te-library` |
| `classplots` | look at the result: per-class and per-locus figures |
| `report` | read the result: one offline HTML page with plain-language explanations |

### Main parameters

**`qc`** — `-p/--prefix`, `-o/--outdir`, `--reads1` (`--reads2` for paired), `-g/--genome`,
`--adapter` (auto-detected if omitted), `--min-len` (default 15), `--qtrim` (20), `--threads`.

**`mapping`** — `--reads`, `-g/--genome` or `--index-prefix`, `--mismatches`,
`--k` (max alignments per read, default 20), `--ncrna-db` (Bowtie index of rRNA/tRNA/snoRNA/
snRNA; drops those reads before clustering), `--threads`. A sorted, indexed BAM is the default
output.

**`mirna`** — `--engine {mirdeep2,mirdeep-p2}` (**no default**: mirdeep2 for animals,
mirdeep-p2 for plants; the wrong one halves precision), `-g/--genome`, `--index-prefix`
(required; without it the engine aligns nothing and reports success with zero predictions),
`--reads`, `--species`, `--mature`/`--hairpin` (miRBase FASTA, optional), `--no-validate`.

**`classify`** — `--sam` (SAM or BAM), `-g/--genome`, `--mincov` (seeding depth in **reads per
million**, default 1.0), `--pad` (island extension in bp; sets the merge distance, and with `trim_to_reads` the reported bounds come from the reads), `--config` (`rules.yaml` or
`rules.plant.yaml`), `--mirna-loci` (BED), `--organelle-contigs` (comma-separated contig names
to drop), `--norm`, `--threads`, `--classifier {rules,ml}` (default `rules`) and `--ml-model`.

**Optional ML classifier.** `--classifier ml` replaces the rule scoring with a trained model
that gives each locus a class (or `ambiguous` when its top probability is below the model's
threshold) and `p_real`, the probability that the locus is a real small-RNA locus. Loci carved
from `--mirna-loci` keep `miRNA` and get `p_real = NA`. `--ml-model` points to the `.joblib`;
by default it is the `ml.model` key of `--config` (`models/arthropod.joblib` for `rules.yaml`,
`models/plant.joblib` for `rules.plant.yaml`), read from the repo root. A `<model>.schema.json`
(features, classes, threshold, algorithm, seed, commit) must sit beside it, and its feature list
must equal `mlfeatures.FEATURES`. The default `rules` mode needs none of the packages below:
`pip install scikit-learn joblib` (plus `xgboost` or `lightgbm` if the model uses them). Models
are produced by `python -m benchmark ml-train --profile {arthropod,plant} --algorithm <name>
-o models/<profile>.joblib`.

**`classplots`** — `--clusters` (the `.clusters.tsv`), `--bam`, `--top-n`, `--classes`,
`--region`, `--logo-format`.

---

## The config

One file per kingdom, and it is meant to be edited. Everything below is per class.

**Locus calling**

| key | meaning |
|---|---|
| `clustering.pad` | how far a coverage island grows on each side before neighbours may merge |
| `clustering.locus_min_reads` | read mass (`sum 1/NH`) below which a region is not a locus but a fragment (shipped: 20) |
| `clustering.locus_min_span` | read span in bp (first read start to last read end) below which a region is a fragment (shipped: 50) |
| `clustering.locus_span_exempt_reads` | the span test applies only below this read mass (shipped: 100), so a narrow stack with that many reads, such as an unannotated miRNA or a tRNA/snoRNA fragment, stays a locus. `null`/absent: the span test always applies. `locus_min_reads` is unaffected |
| `clustering.trim_to_reads` | `true`/`false` (strict; other values are an error). When on, each locus is bounded by its reads (clipped to the merged window) instead of the padded window. Known-miRNA loci (`--mirna-loci`) are never trimmed or made fragments |
| `extractor.dicermin` / `dicermax` | read lengths kept (18–35) |
| `extractor.strand_cutoff` | above this share on one strand, the locus is called single-stranded |

**Decision**

| key | meaning |
|---|---|
| `decision.review_confidence` | below this softmax confidence, the locus is flagged `review`. Its rigour depends on how many classes exist, so it is not comparable between configs |
| `decision.min_reads` | weighted read count below which a locus is flagged `review` (a flag, not a filter; use `--mincov` to remove noise). With the locus floor on, loci below `clustering.locus_min_reads` are already fragments, so with the shipped configs (floor 20 > 10) this flag never fires on non-miRNA loci |
| `decision.min_margin` | minimum gap to the runner-up; a no-op with only two classes |

**Per class**

| key | meaning |
|---|---|
| `size: {lo, hi}` | the length window that defines the class |
| `w_size` | weight of falling inside that window |
| `w_1u` / `w_5a` | weight of a 5′ U (piRNA, miRNA) or 5′ A (AGO4-loaded hc-siRNA) |
| `w_fractop_mid` | rewards a balanced strand ratio (double-stranded precursor) |
| `w_pingpong` | rewards 10-nt 5′–5′ overlap between strands (secondary piRNA) |
| `w_pingpong_penalty` | penalizes it, for classes where it should be absent |
| `w_phase` | rewards periodic 5′ spacing (phasiRNA, Zucchini-processed piRNA) |
| `phase_x0`, `phase_scale` | midpoint and width of the sigmoid applied to `phase_score` |
| `phase_on: cycles` | score phasing by the run of cycles in register rather than by depth |
| `enabled: false` | switch a class off entirely. This also removes its window from clustering |
| `defining_evidence` | which signatures must be present for the reported `score`; a nested list is an OR (`[[1u, pingpong], size]` = 5′ U **or** ping-pong, **and** size) |

Every weight is multiplied by its own class's size window, so a signature cannot drag a locus
into a class whose length range it does not occupy.

---

## Report

`report` turns a `classify` run into one offline HTML file (fonts embedded, no network needed) for
readers who are new to small RNAs. It is an opt-in subcommand you run after `classify`; nothing in
`classify` changes.

```
srnaprofiler report --clusters run.clusters.tsv --bam sample.sorted.bam -o out -p run
```

It writes `out/run.report.html`. The miRNA hairpin figure needs ViennaRNA (`RNAfold`, plus `RNAplot` or its Python binding) for the secondary structure; it is in `environment.yml` and `setup --check` lists it under `report`. It only reads the `classify` outputs and the BAM.

The file has three parts:

- **About tab**: the tool name, tagline, version, GitHub link and install commands.
- **Overview and cluster list**: the library summary, the library quality strip (below), then a
  sortable, filterable list of clusters.
- **Cluster page**: the measured values, a "Why this class" checklist (below), then these
  figures: size x strand x 5' base, phasing, whole-cluster coverage, ping-pong, 5' end track,
  nucleotide logo, top sequences and library context (a miRNA cluster matched to a `--mirna` table also gets the hairpin figure, below). Sense is drawn up and antisense down. Every
  figure except the top-sequences one has a table view and a "Download SVG" button.
- **Genomic context** (only with `--annotation`; on a miRNA page it comes before the hairpin): the locus with `--flank` bases each side on its
  contig, one track each for coverage (by read size class), 5' ends (summed per bin, not per base),
  genes (exons, CDS, UTR, strand arrows), transposons, structural RNA and PHAS loci, plus the
  predicted ORFs with `--predict-orfs` and a miRNA precursor window on miRNA clusters. Zoom with the
  buttons, Ctrl/Cmd + wheel (the plain wheel only while the figure has the focus), the keys `+`,
  `-`, `0` and the arrows, and pan by dragging. The Table and Download SVG buttons reflect the
  current view. A layer whose contigs match none of the BAM references says so in its track
  (usually chr1 vs 1 naming).

| flag | default | meaning |
|---|---|---|
| `--min-reads` | 300 | depth floor in raw reads (the `reads` column, 1/NH-weighted): clusters below it are left out of the report and counted as hidden |
| `--max-detail` | 3000 | clusters above the floor get a detail page up to this number (half the budget is split evenly across the classes present, best first; with `--annotation` a fifth is kept for loci whose annotation `conflicts` with the class, since those are what a reviewer opens first; the rest follows the global order, resolved clusters by confidence before ambiguous ones); the rest stay in the list. Each page is about 7 KB (about 13 KB with `--annotation`), so about 7 MB (13 MB) per 1000 pages; above 5000 `report` warns that the file can get large and slow to open in a browser |
| `--threads` | 1 | worker processes for the figure data. Pays off only for large inputs and on Linux (fork) |
| `--fragments` | `<prefix>.fragments.tsv` beside `--clusters` | regions below the `classify` locus floor, shown in the library cards |
| `--mirna` | `<prefix>.mirna.tsv` beside `--clusters` | the `mirna` module's table; draws the miRNA hairpin figure (below) on the page of each miRNA cluster it matches. An explicit file that does not exist is an error; an absent default file just means no figure |
| `--fastp` | `<prefix>.fastp.json` beside `--clusters` | fastp JSON for reads before/after trimming and the adapter-trimmed share in the library quality strip |
| `--qc-summary` | `<prefix>.qc_summary.txt` beside `--clusters` | `qc` summary, the fallback for reads before/after trimming |
| `--no-library-qc` | off | skip the extra BAM pass that builds the library quality strip |
| `--title` | built-in | page title |
| `--annotation` | none | `LABEL:PATH`, repeatable, same labels (`gene`, `te`, `ncrna`, `phas`) and formats as `classify --annotation`. Adds the "Genomic context" figure (below) to each cluster page |
| `--flank` | 2000 | bases shown each side of the locus in that figure, 0 to 100000. Only with `--annotation` (otherwise `report` warns that it is ignored) |
| `--predict-orfs` | off | EXPERIMENTAL. Adds a track of ORFs predicted from the genome sequence (not annotation, never used by the classifier). Needs `--annotation` and `--genome` |
| `--genome` | none | genome FASTA indexed with `samtools faidx` (`.fai`), read by `--predict-orfs` only |
| `--orf-min-aa` | 100 | shortest predicted ORF, in amino acids (only with `--predict-orfs`) |

**Size.** `--annotation` adds about 6 KB to each detail page, so a report with it is about 13 KB
per page. The coverage series (`cv`) dominates that and does not depend on `--flank`; what controls
the file size is `--max-detail` (about 80 MB at `--max-detail 6000` on a Drosophila run). `report`
logs the bytes the genome view added and warns when the finished file is over 50 MiB; lower
`--max-detail` to shrink it.

**Use the same annotation in both commands.** The annotation text and support chip in a cluster's
header come from `classify --annotation`; the Genomic context figure comes from `report
--annotation`. Pass the SAME files to both or they can disagree. `report` reads `<prefix>.run.json`
beside `--clusters` (where `classify` records its annotation) and logs one warning naming the
difference when the labels and file names differ, or when only one of the two had an annotation.

The BAM must be indexed (`samtools index`); without the index `report` stops with an error.
The floor only filters the view. It does not change how loci are called or classified.

The kingdom (animal or plant) is detected from the clusters file: a run counts as plant when the
table has a plant-only class (`phasiRNA`, `phasiRNA24`, `hcsiRNA`) as a `gate_*` column or a call.
The class definitions are read from the shipped `rules.yaml` / `rules.plant.yaml`, so a run with a
custom `--config` may define classes differently. The evidence dots and strip use display-only
cutoffs, stated in the page footer; they are not the classifier's thresholds. Run `classify` with
a version that writes the `e_*` columns, or the dots show "not available".

### Why this class

Each cluster page lists, before the figures, one row per check the class rules apply: the check,
the value observed in the cluster, the rule, and the result as an icon plus text (Met, Not met,
Not available, Not tested). The rows come from the shipped `rules.yaml` / `rules.plant.yaml`,
so they follow what the classifier does:

- **Size window** (for example piRNA 24-31 nt, siRNA 20-22 nt): the dominant size (`dicer_call`)
  must fall inside it. This is a hard limit; outside every window the cluster is ambiguous.
- **Defining evidence** (for example 5' U or ping-pong for piRNA, strand balance for siRNA,
  phasing for phasiRNA, 5' A for hcsiRNA). The classifier weighs these in a score and has no cutoff
  on them, so a row is *Met* when its effect strength (the `e_*` column, a 95% Wilson lower bound
  against the signature's null) is above 0, that is, above chance. The observed percentage or
  statistic from the table is shown next to the strength. `phasiRNA24` also needs the phasing
  test at p of 0.01 or less, as the rules require. Without the `e_*` columns (a `classify` run
  from before they were written) these rows read *Not available* and only the raw value shows.
- **Cluster criteria** (piRNA: reads span at least 1000 bp and denser than chance, p 0.01):
  read from the cluster's evidence trail. The observed length is the `read_span` column and the
  density row shows `density_pvalue`; a table from before these columns existed falls back to
  the locus window (end minus start) and shows no density value.
- An **ambiguous** cluster lists the checks of its nearest class, so you can see which failed;
  criteria are *Not tested* when the size already failed, as in the classifier. An ambiguous
  call from `--classifier ml` has no nearest class and no checklist. For an ML call the rows are
  the class's basic checks and the evidence rows are for reference.
- **miRNA** clusters get one line instead: the class comes from the miRDeep2 hairpin search.

The rows reflect the shipped rules files; a run with a custom `--config` may differ.

### miRNA hairpin

With a `mirna` table (`--mirna`, or `<prefix>.mirna.tsv` beside `--clusters`), the page of a
cluster classified miRNA gets, before the other figures, a figure "miRNA hairpin" and two small
tables. A table row is matched to a cluster by chromosome and overlap, the rule of
`classify --mirna-loci`: the row's interval must cover at least 50% of the cluster's span; with
several candidates the largest overlap wins, a `Pass` before a `Fail`. Only clusters that have a
detail page are drawn.

- **Precursor and fold.** The precursor sequence of the table (5' to 3', on the miRNA's strand, as
  RNA) in monospace, with the mature, loop and star segments tinted red, yellow and purple as
  miRDeep2 does and also named in text. The colours are the report's existing red, amber and
  purple; the letters stay in dark ink on a light tint, so the yellow does not have to carry
  contrast by itself. Below it the dot-bracket structure from RNAfold (ViennaRNA), found as the
  `RNAfold` executable on `PATH`, else the python module `RNA`. Without either the structure line is
  left out and the figure says "Structure not available: RNAfold not installed." Nothing is cached
  on disk. The sequence is the one miRScore was given: the hairpin plus up to 25 nt of genomic
  flank on each side.
- **2D structure.** At the top of the figure, the classic secondary-structure drawing of the hairpin
  (ViennaRNA naview layout via the python module, else `RNAplot -t 1`, 5' end at the top-left),
  with bonds, mature/loop/star tints, every 10th base numbered and the MFE. It scales to the card
  width; letters appear only when the bases are at least 13 px apart. The per-base table view adds
  the paired position. Without a layout tool the panel is omitted.
- **Reads.** Read depth along the precursor (all collapsed reads), and under the sequence the BAM
  reads that lie fully inside the precursor and align without indels or clipping, on the
  precursor's strand (for a minus-strand precursor the reads are reverse-complemented and the
  coordinates flipped). They are collapsed by position and sequence, counted with the usual 1/NH
  weights, and the 25 most abundant are drawn with their read count as a number and a bar. Lower
  case matches the precursor; an upper-case letter on a dark cell is a mismatch. A read is marked
  mature or star when at least half of it overlaps that segment. Reads on the other strand and
  reads running past the precursor ends are counted in a note, not drawn. The precursor is placed on
  the genome from the table's `start`/`end` and the flank rule above; if its length does not fit,
  the sequence is drawn without reads and a note says so. The reads are read in the same pass over the
  BAM as the other figures. The table view lists every position and the drawn read rows; Download SVG
  exports the figure (which scrolls sideways when the precursor is wide).
- **Numbers.** Mature and star reads, precision, mismatches, bulges, precursor length (with
  flanks) and hairpin locus length, miRDeep2 score and probability, randfold p-value, Rfam alert,
  miRBase hit and same-seed hit, as given in the table (`NA` when the run lacked them).
- **Plant criteria, shown for reference.** Three criteria of Axtell and Meyers (2018) that the
  table can answer, each with Met / Not met / Not available as icon plus text, for every kingdom:
  hairpin locus of 300 nt or less, at most 5 mismatches in the duplex, precision of at least 75%.
  They are not used to call the class. The other criteria of that paper are not checked.

Not drawn: a mature/star duplex view with the 2-nt 3' overhang, which would need the base pairs
resolved from the fold; the duplex's mismatch count comes from the table instead. A table without
a precursor sequence (miRScore did not run) gives the numbers and a message instead of the figure.

Reference: Axtell, M. J., & Meyers, B. C. (2018). Revisiting criteria for plant MIRNA annotation
in the era of big data. *The Plant Cell*, 30(2), 272-284. https://doi.org/10.1105/tpc.17.00851

### Library quality

The Overview tab ends with a library quality strip, in the style of miRTrace. A whole-library
read-length histogram (16 to 40 nt, plus one bin for shorter and one for longer reads) is built
from the BAM in one extra pass (by contig, parallel with `--threads`). Reads are the primary
mapped alignments, weighted by 1/NH like the cluster counts, so the total is on the scale of the
library size behind RPM. Unmapped reads are not in the histogram, so shares are of mapped reads;
the number of unmapped records is read from the BAM index and shown when it is not zero. Three flags,
each with its number and a text status:

| flag | status Check when |
|---|---|
| share of reads of 20 to 25 nt | under 25% |
| dominant length | outside 18 to 35 nt, the sizes the classifier keeps |
| share of reads shorter than 18 nt | more than 25% |

The 25% floor on the 20 to 25 nt share is the cutoff the strip is built on; the other two cutoffs
are display-only choices of this report, not miRTrace values. If a `fastp` JSON
(`<prefix>.fastp.json`) or a `qc` summary (`<prefix>.qc_summary.txt`) sits beside `--clusters`,
or is given with `--fastp` / `--qc-summary`, the strip also shows reads before and after trimming
and, from the fastp JSON, the share of reads with an adapter trimmed. Absent or unreadable files
just leave that part out. `--no-library-qc` skips the extra BAM pass and the whole strip.

---

## Annotation (optional)

`classify --annotation LABEL:PATH` (repeatable) reads genome annotation as extra evidence for each
locus. It is off by default. The class is never changed by it: annotation only nudges the
confidence and can flag a locus for review.

```
srnaprofiler classify ... --annotation gene:genes.gff3 --annotation te:repeats.out.gz
```

- **Formats**: GFF3, GTF, BED, RepeatMasker `.out` and EDTA `.TEanno.gff3`, each plain or `.gz`.
- **Labels**: `gene` (exon, intron), `te` (family and class when the file gives them), `ncrna` and
  `phas` (known PHAS loci). A locus is tagged with a layer when the layer covers at least half of it
  (`annotation.min_overlap`).
- **Output**: two columns just before `review`, `annotation` (for example
  `te:LTR/Gypsy:antisense:0.82|gene:HA8:sense:0.10`, empty when nothing overlaps) and
  `annotation_support` (`supports`, `conflicts`, `neutral` or `none`). The GFF3 gets the same two
  attributes, and the `report` shows them on the cluster page.

The support call comes from a biology table in `rules.yaml` / `rules.plant.yaml` (`annotation:`
section, editable in a custom `--config`):

| class | supports | conflicts |
|---|---|---|
| piRNA (animal) | TE | gene exon |
| hcsiRNA (plant) | TE | gene exon |
| siRNA | TE, gene exon | none |
| phasiRNA, phasiRNA24 (plant) | known PHAS locus | none |
| miRNA | intergenic, gene intron | TE |

A locus whose contexts match neither column is `neutral`. For `supports` the confidence goes up by
`annotation.nudge` (0.05, capped at 1); for `conflicts` it goes down by the same amount and
`review` is set. `miRNA` and `ambiguous` calls get the annotation text and the support label, but
no nudge and no review flag. A `conflicts` nudge can take the confidence below the class's usual
floor (1/K), so a flagged locus may read lower than any unannotated call of that class could.
In a gene file only gene, transcript and exon-like rows (exon, CDS, UTR, ...) count; region,
chromosome and match rows are ignored, and so is any type containing `transposable_element`,
together with all the rows below it (the Ensembl dm6 GFF3 carries FlyBase TE insertions as
`transposable_element_gene` and `transposable_element`; Araport11 lists TEs as
`transposable_element_gene` > `mRNA` > `exon`/`CDS`: the TE gene rows and their children are all
ignored by the gene layer, since TEs belong to the `te` layer, not to the genes). A locus on an annotated contig that overlaps none of them is
`intergenic` (written as `gene:intergenic:.:1.00`). Text entries are
`kind:name:relation:fraction` (for a gene the name is the gene's Name when the file gives one, for
a TE the family or class), plus the feature ID as a fifth field when it has no `|` or `:`.
The report adds a per-class table of supports/conflicts/neutral over the loci it lists.

Annotation is incomplete, above all for TEs and for non-model species, so a missing feature is
never evidence against a class. Only a feature that is there can support or conflict.
An absence-type `conflicts` verdict (gene exon, intron or intergenic with NO TE) claims that the
supporting feature is missing, so it is only given when every layer kind that could have supplied
the class's `supports` evidence was loaded; with a gene layer only, piRNA and hcsiRNA loci on exons
stay `neutral` (and a note is logged). A presence-type conflict (a miRNA locus overlapping a TE)
needs only the te layer.

**TE annotation without a whole-genome run.** `annotate-loci` runs RepeatMasker on the classified
loci only, which is far cheaper than masking the genome:

```
srnaprofiler annotate-loci --clusters run.clusters.tsv -g genome.fa -o out -p loci --species insects
srnaprofiler classify ... --annotation te:out/loci.te.gff3
```

RepeatMasker caps sequence ids at 50 characters, so the loci go in as `L1`, `L2`, ... and
`out/loci.loci.map.tsv` (id, chrom, 0-based start, end) maps them back to the genome. A failed run
removes its scratch files and can be rerun; only an existing `loci.te.gff3` blocks it.

Use `--te-library` instead of `--species` for a custom library. RepeatMasker is an optional
install (`environment-te.yml`, below).

**Building a TE library (`te-library`).** Three routes to a TE annotation: (1) a curated
library from Dfam/RepeatMasker (`annotate-loci --species`), enough for model clades; (2) a library built for
the genome with `te-library`; (3) an annotation you already have (`--annotation te:...`). Measured on
*A. thaliana*, Dfam 4.0 curated has 9 families and RepeatMasker with it found 0 TEs in 22584 loci
(recall 0.2% against TAIR10); RepeatModeler2 on the genome (3h49, 32 threads) gave 631 families and,
used in `annotate-loci`, 68.8% recall and 83.2% precision (130 s for the loci). RepeatModeler2 libraries
can also hold multicopy gene families, so `te-library` filters them out. On that library the filter
(24 s) dropped 102 families as genes and 1 as snRNA, and moved recall from 68.8% to 66.9% and precision
from 83.2% to 87.4%; the hcsiRNA loci called TE by the library alone that sit on a gene exon fell from
20.8% to 11.9%:

```
srnaprofiler te-library -g genome.fa -o OUTDIR -p PREFIX [--threads N] [--resume]
    [--organelle-contigs Pt,Mt] [--genes GENES.gff3|gtf[.gz]] [--exclude-fasta FASTA ...]
    [--dfam-species "Latin name"] [--families EXISTING.fa]
srnaprofiler annotate-loci --clusters run.clusters.tsv -g genome.fa -o out --te-library OUTDIR/PREFIX-families.fa
```

Steps: write the nuclear genome (without `--organelle-contigs`); `BuildDatabase` + `RepeatModeler -LTRStruct`
(`--resume` passes `-recoverDir` to an interrupted run; `--families` skips this step); drop each family
that `minimap2 -x asm20` aligns to the exons of the non-TE genes (`--genes`) or to `--exclude-fasta`
sequences with identity >= 0.80 over >= 50% of the family; families classed rRNA, tRNA, snRNA, scRNA
or srpRNA are always dropped; optionally append curated Dfam families (`--dfam-species`, needs `famdb.py`).
It writes `OUTDIR/PREFIX-families.fa` and `OUTDIR/PREFIX.te-library.report.tsv` (kept/dropped, reason,
matched sequence per family). Without `--genes`/`--exclude-fasta` no sequence filtering happens.
RepeatModeler, minimap2 and famdb.py are optional installs; `classify` and `report` never call them.

**Installing the TE tools.** They live in their own environment, so the main one stays small:

```
conda env create -f environment-te.yml        # RepeatMasker, RepeatModeler2, rmblast, LTR tools, minimap2
conda activate srnaprofiler-te
```

Run `annotate-loci` and `te-library` inside it (it also carries the Python packages srnaprofiler needs).
RepeatMasker has no repeat database of its own. For `--species` it needs Dfam (FamDB): download
`dfam40.0.h5.gz` and `dfam40.curated.consensus.0.h5.gz` (about 90 MB, curated consensi for all
taxa) from https://www.dfam.org/releases/current/families/FamDB/, gunzip them into
`<env>/share/RepeatMasker/Libraries/famdb/`, and set `FAMDB_DATA_DIR` in `<env>/share/famdb-*/famdb.conf`
to that folder. `--te-library` needs no database at all. Check the coverage of your clade first
(`famdb.py -i <dir> lineage -a "Genus species"`): *A. thaliana* has 9 curated families, *D. melanogaster* 399.

Without `--annotation`, nothing changes: `classify` writes byte-identical outputs to a run before
the layer existed, which is what keeps the benchmark numbers valid.

---

## Output

`classify` writes `<prefix>.clusters.tsv`, `<prefix>.clusters.gff3`, `<prefix>.run.json` and
`<prefix>.fragments.tsv`.

`<prefix>.fragments.tsv` lists the regions rejected by the locus floor, so nothing is dropped
silently. Columns: `chrom`, `start`, `end` (the region, trimmed when `trim_to_reads` is on),
`reads` (1/NH-weighted read mass), `read_span` (bp from first read start to last read end; 0 for a
region with no reads) and `reason` (`min_reads`, `min_span` or `min_reads,min_span`). Fragment
reads are not reassigned to neighbouring loci. The file is header-only when the floor is off or
nothing falls below it.

**Counts are 1/NH-weighted.** A read mapping to *N* places contributes 1/*N* to each, so
counts are *not* comparable with tools that assign a multi-read to a single position. Those
are larger by the multiplicity factor, 3.19× in one Drosophila library.

### Main columns

| column | meaning |
|---|---|
| `chrom`, `start`, `end`, `name` | locus position (0-based, half-open) |
| `reads`, `reads_rpm` | weighted abundance, raw and per million |
| `distinct_seqs` | distinct sequences. The sample size of the binomial tests, not an abundance |
| `major_rna`, `dicer_call` | dominant read length, and the length class called from it |
| `frac_in_window` | share of the mass inside the winning class's size window |
| `size_purity` | share of the mass within ±1 nt of the peak. Low means two populations share the locus |
| `dicer_call_2` | a second peak ≥2 nt away carrying ≥20% of the mass |
| `frac_top`, `strand_call` | share on the plus strand, and the resulting call |
| `pct_1u`, `pct_5a`, `pct_10a` | 5′ U, 5′ A, and adenine at position 10 |
| `gc` | GC content |
| `pingpong_z`, `pingpong_frac` | strength of the 10-nt 5′–5′ overlap signature |
| `pingpong_pairs10` | pairs available to measure it. **Filter on this before averaging `pingpong_frac`**, which is 0 below 10 pairs, meaning not measurable rather than absent |
| `phase_score`, `phase_number`, `phase_ratio` | phasing strength, cycles in register, and share of the mass in the dominant register |
| `phase_pvalue` | hypergeometric test, Bonferroni-corrected over the scanned windows |
| `multimap_frac`, `unique_frac` | share of alignments that are multi-mapping, and share of the mass from uniquely mapping reads. `unique_frac` near 0 marks a locus whose signature is borrowed from repeats |
| `cls`, `confidence` | assigned class and softmax confidence |
| `review` | true when the classes compete, the locus is too shallow, or evidence is thin. A flag to inspect, not a rejection |
| `p_real` | `--classifier ml` only: probability the locus is real; `NA` in rules mode and for carved miRNA |
| `score` | strength of the weakest defining evidence, in [0, 1]; a normalized Wilson lower bound against each signature's null, so it ranks reliability rather than sample size |
| `s_1u`, `s_5a`, `s_fractop_bid`, `s_pingpong`, `s_phase` | the individual p-values, as −log10 |
| `evidence` | which terms fired for the assigned class |
| `annotation`, `annotation_support` | only with `--annotation`: overlapping features and whether they support or conflict with the class (see Annotation, above) |

`score` and the `s_*` columns answer different questions: `s_*` say how surprising a signature
is, which grows with depth; `score` says how strong the effect is, which does not. Use `score`
to compare loci.

## License

MIT — see [LICENSE](LICENSE).
