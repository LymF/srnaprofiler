# sRNAProfiler

Small RNA loci discovery and classification from sequencing data. It calls loci from
coverage, measures the signatures each biogenesis pathway leaves behind (size, strand
balance, 5′ nucleotide, ping-pong overlap, phasing), and assigns a class with a transparent,
editable rule set. Plants and animals share one engine, and the biology lives in a config file.

Classes come out of `rules.yaml` (animal: siRNA, piRNA) or `rules.plant.yaml` (plant:
phasiRNA, hcsiRNA, siRNA, and optional 24-nt phasiRNA). miRNA is not called by signature; it
is taken from annotation or prediction and carved out of the surrounding locus.

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
| `classify` | call loci, measure their signatures, assign a class |
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

| flag | default | meaning |
|---|---|---|
| `--min-reads` | 300 | depth floor in raw reads (the `reads` column, 1/NH-weighted): clusters below it are left out of the report and counted as hidden |
| `--max-detail` | 3000 | clusters above the floor get a detail page up to this number (resolved clusters by confidence first, ambiguous ones after them); the rest stay in the list. Each page is about 7 KB, so about 7 MB per 1000 pages; above 5000 `report` warns that the file can get large and slow to open in a browser |
| `--threads` | 1 | worker processes for the figure data. Pays off only for large inputs and on Linux (fork) |
| `--fragments` | `<prefix>.fragments.tsv` beside `--clusters` | regions below the `classify` locus floor, shown in the library cards |
| `--mirna` | `<prefix>.mirna.tsv` beside `--clusters` | the `mirna` module's table; draws the miRNA hairpin figure (below) on the page of each miRNA cluster it matches. An explicit file that does not exist is an error; an absent default file just means no figure |
| `--fastp` | `<prefix>.fastp.json` beside `--clusters` | fastp JSON for reads before/after trimming and the adapter-trimmed share in the library quality strip |
| `--qc-summary` | `<prefix>.qc_summary.txt` beside `--clusters` | `qc` summary, the fallback for reads before/after trimming |
| `--no-library-qc` | off | skip the extra BAM pass that builds the library quality strip |
| `--title` | built-in | page title |

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
  read from the cluster's evidence trail. The table does not carry the read span or the density
  p-value, so the observed length is the locus window (equal to the read span when
  `trim_to_reads` is on) and the density row shows no value.
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

`score` and the `s_*` columns answer different questions: `s_*` say how surprising a signature
is, which grows with depth; `score` says how strong the effect is, which does not. Use `score`
to compare loci.

## License

MIT — see [LICENSE](LICENSE).
