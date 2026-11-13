# Annotation as evidence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `classify` optionally takes genome annotation layers (genes, TEs, structural RNA, PHAS), describes every locus against them, and uses a small capped confidence nudge plus a review flag when the overlap supports or contradicts the class call.

**Architecture:** One new module `annotation.py` (parsers, interval index, overlap summary, support rule, nudge) called from `classify.py` after classification and before `write_tsv`. New columns are written only when `--annotation` is given, so runs without it stay byte-identical. A second module `annotate_loci.py` runs RepeatMasker on the locus sequences only and emits a GFF that `annotation.py` reads. The HTML report reads the new columns when present.

**Tech Stack:** Python 3, pysam (FASTA access), pyyaml, pytest. No new dependency: `pybedtools` is not installed in every env, so the interval index is pure Python.

**Spec:** `docs/superpowers/specs/2026-11-13-annotation-evidence-design.md`

## Global Constraints

- Without `--annotation`, `clusters.tsv`, `.gff3` and `run.json` are byte-identical to today (the benchmark must not move). The `annotation` key in `run.json` appears only through `params`, which already dumps all args: a new `annotation=None` arg would change `run.json`, so `vars(args)` must drop it when `None`.
- Annotation never changes `cls`. It changes only `confidence` (by at most `annotation.nudge`, default `0.05`), `review` (set True on `conflicts`, never cleared) and the new columns.
- Absence of an annotated feature is never evidence: no overlap is `neutral`.
- Layer label must be one of `gene`, `te`, `ncrna`, `phas`; anything else is a usage error naming the allowed labels. (My call, see Decisions.)
- Formats: GFF3, GTF, BED, RepeatMasker `.out`, EDTA `.TEanno.gff3` (read as GFF3), each plain or `.gz`.
- Contig names are matched exactly to the BAM's; features on unknown contigs are skipped and reported once with a count.
- Tests live in `benchmark/tests/` (`test_annotation_*.py`), import modules from the repo root like the existing tests, and use synthetic fixtures only.
- Environment: the default `python3` here lacks pytest/pysam. Before Task 1, create a venv in the session scratchpad directory (not `/tmp`) and `pip install pytest pysam pyyaml biopython pandas matplotlib`. All `Run:` lines use `python -m pytest` from that venv at the repo root.
- Commit messages end with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.

## Review Focus

1. Annotation file with `chr1` while the BAM says `1` (or `NC_000001.11`): expect one warning with the skipped count, the run continues, loci come out `none`/`neutral`, not a crash.
2. Locus on a contig absent from the annotation: `neutral`, empty `annotation` cell (`NA`), no error.
3. GFF with `Parent` cycles, features with `start > end`, missing strand, or `#`/`###`/`##FASTA` sections: skipped or tolerated, never an exception, never reading the embedded FASTA as features.
4. Locus where the class has no row in the biology table (`ambiguous`, `miRNA`, unknown plant classes): `neutral`, confidence and review untouched.
5. Nudge at the bounds: confidence 0.98 + support never exceeds 1.0; 0.02 + conflict never below 0.0; `conflicts` on a locus already `review=True` stays True.

---

## File Structure

- Create `annotation.py`: parsers, `IntervalIndex`, `Layer`, `load_layers`, `summarize_locus`, `support_for`, `apply_annotation`. One responsibility: turn annotation files into per-locus evidence.
- Create `annotate_loci.py`: locus FASTA extraction, RepeatMasker call, `.out` to genome-coordinate GFF.
- Modify `records.py`: two fields on `ClusterRecord`.
- Modify `reporter.py`: optional extra columns in `write_tsv`, attributes in `write_gff3`.
- Modify `classify.py`: call the annotation step.
- Modify `srnaprofiler.py`: `--annotation` on `classify`; `annotate-loci` subcommand.
- Modify `rules.yaml`, `rules.plant.yaml`: `annotation:` block (biology table, nudge, min_overlap).
- Modify `report_data.py`, `report_html.py`: show the new columns.
- Modify `README.md`: usage and the biology table.

---

### Task 1: Interval index and annotation file parsers

**Files:**
- Create: `annotation.py`
- Test: `benchmark/tests/test_annotation_parse.py`

**Interfaces:**
- Produces:
  - `Feature = namedtuple("Feature", "start end strand ftype name")` (0-based half-open, `strand` in `"+","-","."`, `ftype` the GFF type or `"feature"` for BED/RM, `name` the family/class string or `""`)
  - `parse_annotation(path: str) -> dict[str, list[Feature]]` keyed by contig; format chosen by name: `.out`/`.out.gz` RepeatMasker, `.bed`, `.gtf`, otherwise GFF3 (also `.gff`, `.gff3`, `.TEanno.gff3`); transparently gunzips.
  - `class IntervalIndex` with `__init__(features: list[Feature])` and `overlaps(start: int, end: int) -> list[Feature]` (features overlapping `[start, end)`).

- [ ] **Step 1: Write failing tests** in `benchmark/tests/test_annotation_parse.py`, tmp-path fixtures written inline:
  - `test_gff3_exon_and_gene_types`: two lines (`gene`, `exon`) on `chr1`; `parse_annotation` returns both with 0-based start (`GFF start 101 -> 100`), strand kept.
  - `test_gff3_te_name_from_classification_attr`: EDTA-style line with `Name=TE_1;classification=LTR/Gypsy` → `name == "LTR/Gypsy"`; with only `Class=DNA/MITE` or `family=Gypsy` falls back in that order (`classification`, `Class`, `class`, `family`, `Family`), else the type column.
  - `test_gff3_ignores_comments_and_fasta_section`: lines after `##FASTA` are not read; `#` lines skipped; a line with `start > end` skipped.
  - `test_gtf_attributes`: `gene_id "g1"; transcript_id "t1";` parses, `exon` type kept.
  - `test_bed_three_and_six_columns`: strand `.` when absent; `track`/`#` lines skipped.
  - `test_repeatmasker_out`: the 3 header lines are skipped; a row with strand `C` gives `-`; `name` is the `class/family` column; coordinates `begin-1, end`.
  - `test_gzip_input`: same GFF gzipped parses identically.
  - `test_index_overlaps`: features `[0,10)`, `[5,100)`, `[200,210)`; `overlaps(8, 12)` returns the first two, `overlaps(100, 150)` returns none (half-open), a very long early feature is still found for a late query.

- [ ] **Step 2: Run:** `python -m pytest benchmark/tests/test_annotation_parse.py -v` → FAIL (module missing).

- [ ] **Step 3: Implement** `Feature`, `parse_annotation`, `IntervalIndex` in `annotation.py`. Index per contig: features sorted by start, a parallel prefix-max of ends; `overlaps` bisects the last start `< end` and walks back while prefix-max `> start`. Parsers are generators over lines; malformed lines are skipped, never raised.

- [ ] **Step 4: Run** the same command → PASS.

- [ ] **Step 5: Commit**
  ```bash
  git add annotation.py benchmark/tests/test_annotation_parse.py
  git commit -m "feat(annotation): parsers for GFF3/GTF/BED/RepeatMasker and an interval index"
  ```

---

### Task 2: Layers and per-locus overlap summary

**Files:**
- Modify: `annotation.py`
- Test: `benchmark/tests/test_annotation_summary.py`

**Interfaces:**
- Consumes: `Feature`, `parse_annotation`, `IntervalIndex` (Task 1).
- Produces:
  - `LAYER_KINDS = ("gene", "te", "ncrna", "phas")`
  - `class Layer` with attributes `kind: str`, `path: str`, `index: dict[str, IntervalIndex]`, `n_features: int`
  - `parse_layer_spec(spec: str) -> tuple[str, str]` splits `"te:/path/x.out"` into `("te", "/path/x.out")` on the first `:`; raises `ValueError` listing `LAYER_KINDS` when the label is not one of them or the path part is empty.
  - `load_layers(specs: list[str], known_contigs: set[str]) -> list[Layer]` parses each layer, drops contigs not in `known_contigs`, and logs one `WARN` per layer with the count of skipped features and a 3-name sample of unmatched contigs (via `utils.log`).
  - `summarize_locus(chrom: str, start: int, end: int, strand: str, layers: list[Layer], min_overlap: float) -> dict` returning `{"contexts": set[str], "text": str}`:
    - contexts drawn from `te`, `gene_exon`, `gene_intron`, `intergenic`, `ncrna`, `phas`.
    - A layer contributes its context when the union of its overlapping features, clipped to the locus, covers `>= min_overlap` of the locus length.
    - `gene` layer: features whose `ftype` is in `{"gene","mRNA","transcript","ncRNA_gene"}` are span-type, everything else (exon, CDS, BED rows) is exon-like. Exon-like union `>= min_overlap` → `gene_exon`; else span union `>= min_overlap` → `gene_intron`; else, if any gene layer is loaded and no layer contributed anything → `intergenic`.
    - `text` is `"|".join` of one entry per contributing layer: `kind:name:relation:fraction` where `name` is the most-covering feature's `name` (or `ftype` if empty), `relation` is `sense`/`antisense`/`both`/`.` of the locus `strand` against that feature's strand, `fraction` the layer's covered fraction to 2 decimals. Empty string when nothing contributes.

- [ ] **Step 1: Write failing tests** in `benchmark/tests/test_annotation_summary.py`:
  - `test_parse_layer_spec_ok_and_errors`: `"te:a.out"` ok; `"foo:a.out"` raises `ValueError` whose message contains `gene, te, ncrna, phas`; `"te:"` raises.
  - `test_te_context_and_text`: TE feature `LTR/Gypsy` on `-` covering a `+` locus fully → contexts `{"te"}`, text `te:LTR/Gypsy:antisense:1.00`.
  - `test_below_min_overlap_is_not_a_context`: feature covers 30% of the locus, `min_overlap=0.5` → empty contexts, empty text.
  - `test_union_not_sum`: two overlapping TE features that together cover 60% (not 120%) → fraction `0.60`.
  - `test_gene_exon_vs_intron_vs_intergenic`: `gene` + `exon` features → exon locus `gene_exon`; locus inside the gene but outside exons → `gene_intron`; locus on nothing with a gene layer loaded → `intergenic`; the same nothing-locus with only a `te` layer loaded → empty contexts (no `intergenic` claim without a gene layer).
  - `test_bed_gene_layer_counts_as_exon_like`: BED rows → `gene_exon`.
  - `test_unknown_contig_locus_is_empty`: locus on `chrX` absent from the layer → empty, no exception.
  - `test_load_layers_skips_unknown_contigs_and_warns`: layer file with `chr1` and `scaf9`, `known_contigs={"chr1"}`, a captured `log` call mentions `1` skipped and `scaf9`; `Layer.index` has only `chr1`.

- [ ] **Step 2: Run:** `python -m pytest benchmark/tests/test_annotation_summary.py -v` → FAIL.

- [ ] **Step 3: Implement** the above in `annotation.py`. Coverage union: merge the clipped intervals, sum lengths.

- [ ] **Step 4: Run** → PASS.

- [ ] **Step 5: Commit**
  ```bash
  git add annotation.py benchmark/tests/test_annotation_summary.py
  git commit -m "feat(annotation): annotation layers and per-locus overlap summary"
  ```

---

### Task 3: Support rule, nudge, config biology table

**Files:**
- Modify: `annotation.py`, `rules.yaml`, `rules.plant.yaml`
- Test: `benchmark/tests/test_annotation_support.py`

**Interfaces:**
- Consumes: `Layer`, `summarize_locus` (Task 2); `records.ClusterRecord` (gains two fields in this task).
- Produces:
  - `ClusterRecord.annotation: str = ""` and `ClusterRecord.annotation_support: str = "none"` in `records.py`.
  - `support_for(cls: str, contexts: set[str], table: dict) -> str` returning `"supports"`, `"conflicts"` or `"neutral"`: a class absent from `table` or with no contexts → `neutral`; contexts ∩ `supports` non-empty → `supports`; else contexts ∩ `conflicts` non-empty → `conflicts`; else `neutral`. Supports wins over conflicts (a piRNA on a TE inside an exon is supported).
  - `apply_annotation(records: list, layers: list[Layer], cfg: dict) -> dict[str, int]` where `cfg` is the config's `annotation` block (`nudge`, `min_overlap`, `classes`). For every record sets `annotation`, `annotation_support`. For classes in the table and `cls not in {"miRNA", "ambiguous"}`: `supports` → `confidence = min(1.0, round(confidence + nudge, 3))`; `conflicts` → `confidence = max(0.0, round(confidence - nudge, 3))` and `review = True`. Returns counts per support value.
  - `rules.yaml` and `rules.plant.yaml` gain at top level:
    ```yaml
    annotation:
      min_overlap: 0.5     # fraction of the locus a layer must cover to count
      nudge: 0.05          # max confidence change from annotation, up or down
      classes:             # context kinds: te, gene_exon, gene_intron, intergenic, ncrna, phas
        piRNA: {supports: [te], conflicts: [gene_exon]}
        siRNA: {supports: [te, gene_exon], conflicts: []}
        miRNA: {supports: [intergenic, gene_intron], conflicts: [te]}
    ```
    plant file: `hcsiRNA: {supports: [te], conflicts: [gene_exon]}`, `siRNA: {supports: [te, gene_exon], conflicts: []}`, `phasiRNA: {supports: [phas], conflicts: []}`, `phasiRNA24: {supports: [phas], conflicts: []}`, `miRNA` as above. Add a comment on each file that these are the defaults from the biology table in the spec and that the nudge never changes `cls`.

- [ ] **Step 1: Write failing tests** in `benchmark/tests/test_annotation_support.py`:
  - `test_support_for_table`: piRNA + `{"te"}` → `supports`; piRNA + `{"gene_exon"}` → `conflicts`; piRNA + `{"te","gene_exon"}` → `supports`; piRNA + `set()` → `neutral`; unknown class `"foo"` + `{"te"}` → `neutral`.
  - `test_apply_nudges_and_flags`: build records with `ClusterRecord` (cls piRNA, confidence 0.80, review False) on a TE layer fixture → confidence `0.85`, support `supports`, review False; a piRNA on an exon-only layer → confidence `0.75`, review True, support `conflicts`; `cls` unchanged in both.
  - `test_nudge_is_bounded`: confidence 0.98 + supports → 1.0; 0.02 + conflicts → 0.0.
  - `test_review_never_cleared`: record `review=True`, supports → still True.
  - `test_miRNA_and_ambiguous_untouched`: both get `annotation`/`annotation_support` filled but `confidence` and `review` identical to before.
  - `test_no_layers_marks_none`: `apply_annotation(records, [], cfg)` is not called by classify, but `support_for` with empty contexts is `neutral`; here assert `apply_annotation` with `layers=[]` leaves `annotation_support == "none"` and `annotation == ""` for every record.
  - `test_shipped_configs_have_annotation_block`: load both yaml files; each has `annotation.min_overlap == 0.5`, `annotation.nudge == 0.05`, and every class named in the block exists in that file's `classes` or is `miRNA`.

- [ ] **Step 2: Run:** `python -m pytest benchmark/tests/test_annotation_support.py -v` → FAIL.

- [ ] **Step 3: Implement** `support_for`, `apply_annotation` and the record fields; edit both YAML files.

- [ ] **Step 4: Run** → PASS. Also run `python -m pytest benchmark/tests -q` → no regressions (the new yaml key must not break config readers; if a reader rejects unknown top-level keys, fix that reader).

- [ ] **Step 5: Commit**
  ```bash
  git add annotation.py records.py rules.yaml rules.plant.yaml benchmark/tests/test_annotation_support.py
  git commit -m "feat(annotation): support rule, capped confidence nudge and biology table"
  ```

---

### Task 4: Writers and classify wiring

**Files:**
- Modify: `reporter.py`, `classify.py`, `srnaprofiler.py`
- Test: `benchmark/tests/test_annotation_classify.py`

**Interfaces:**
- Consumes: `load_layers`, `apply_annotation`, `parse_layer_spec` (Tasks 2-3).
- Produces:
  - `write_tsv(records, path, gate_classes=(), annotated=False)`: when `annotated` is True, two columns `annotation`, `annotation_support` are inserted immediately before `review`; when False the file is byte-identical to today.
  - `write_gff3(records, path, annotated=False)`: when True the attribute string gains `;annotation_support=<v>;annotation=<text with ; = , & escaped as %3B %3D %2C %26>`.
  - `classify --annotation LABEL:PATH` (`action="append"`, dest `annotation`, default `None`).
  - In `run_classify`, after the classifier/`gate_flags` block and before `write_tsv`: if `args.annotation`, build layers with `load_layers(specs, set(chroms))` (all BAM references, not only the non-excluded ones), call `apply_annotation(records, layers, config.get("annotation", {}))`, log the per-support counts at INFO. `params` for `run.json` drops `annotation` when it is `None`.

- [ ] **Step 1: Write failing tests** in `benchmark/tests/test_annotation_classify.py`:
  - `test_tsv_unchanged_without_annotation`: write a record list with `write_tsv(..., annotated=False)` and compare the header to the expected existing header (assert `"annotation" not in header`).
  - `test_tsv_columns_with_annotation`: header contains `annotation` and `annotation_support` directly before `review`, and a record's cells carry its values; an empty annotation text stays an empty cell (not `NA`).
  - `test_gff3_attributes_escaped`: annotation text `te:LTR/Gypsy:antisense:0.82|gene:exon:sense:0.10` appears with `|` kept and no raw `;` or `=` inside the value.
  - `test_cli_accepts_repeatable_annotation`: parse `["classify","-p","x","-o","o","--sam","s.bam","--annotation","te:a.out","--annotation","gene:b.gff3"]` via the project's parser builder → `args.annotation == ["te:a.out","gene:b.gff3"]`; without the flag it is `None`.
  - `test_run_json_has_no_annotation_key_when_unused`: the params-filter helper (extract it as `_run_params(args) -> dict` in `classify.py`) omits `annotation` when `None` and keeps it when set.

- [ ] **Step 2: Run:** `python -m pytest benchmark/tests/test_annotation_classify.py -v` → FAIL.

- [ ] **Step 3: Implement** the writer changes, `_run_params`, the CLI flag, and the wiring. An invalid label surfaces as a clean `ERROR` log and non-zero exit, like the other argument errors in `classify.py`.

- [ ] **Step 4: Run** the file → PASS, then `python -m pytest benchmark/tests -q` → no regressions, specifically the existing report and pipeline-regression tests.

- [ ] **Step 5: Commit**
  ```bash
  git add reporter.py classify.py srnaprofiler.py benchmark/tests/test_annotation_classify.py
  git commit -m "feat(classify): --annotation layers, new TSV/GFF3 columns when used"
  ```

---

### Task 5: Annotate only the loci (RepeatMasker)

**Files:**
- Create: `annotate_loci.py`
- Modify: `srnaprofiler.py`, `annotation.py` (RepeatMasker `.out` parser gains a coordinate offset hook)
- Test: `benchmark/tests/test_annotate_loci.py`

**Interfaces:**
- Consumes: `parse_annotation` RepeatMasker branch (Task 1).
- Produces:
  - `write_locus_fasta(clusters_tsv: str, genome: str, out_fa: str, flank: int = 0) -> int` writes one record per locus named `chrom:start-end` (0-based start, as in the TSV), sequence from `pysam.FastaFile(genome)`, returns the number written. Loci longer than 1 Mb are skipped with a WARN (the count is returned).
  - `rm_out_to_gff(out_path: str, gff_path: str) -> int` converts a RepeatMasker `.out` whose query names are `chrom:start-end` into a GFF3 in genome coordinates (`start0 + begin`), type `repeat_region`, attribute `Name=<repeat>;Class=<class/family>`, strand `C`→`-`. Returns the feature count.
  - `run_annotate_loci(args) -> None` runs `RepeatMasker -pa N [-species S | -lib FASTA] -dir DIR locus.fa`, then `rm_out_to_gff`. Requires exactly one of `--species` / `--te-library`; checks `RepeatMasker` is on PATH via `utils.which` and fails with an `ERROR` log naming how to install it.
  - Subcommand `annotate-loci`: `--clusters`, `-g/--genome` (required), `-o/--outdir` (required), `-p/--prefix` (default `loci`), `--species`, `--te-library`, `--threads` (default 4). Output `<outdir>/<prefix>.te.gff3`, to be passed back as `--annotation te:<that file>`.

- [ ] **Step 1: Write failing tests** in `benchmark/tests/test_annotate_loci.py` (no RepeatMasker needed; subprocess is monkeypatched):
  - `test_write_locus_fasta`: tiny genome FASTA + `.fai` (build with `pysam.faidx`), a 2-row clusters TSV → 2 records named `chr1:100-160` etc. with the right bases.
  - `test_rm_out_to_gff_offsets_and_strand`: `.out` row on query `chr1:1000-1200`, `begin=11 end=60`, strand `C`, class `LTR/Gypsy` → GFF line on `chr1`, start `1011`, end `1060` (1-based), strand `-`, `Class=LTR/Gypsy`.
  - `test_run_requires_one_of_species_or_library`: neither and both → `SystemExit`/error with a message naming `--species` and `--te-library`.
  - `test_run_fails_clearly_without_repeatmasker`: `which` patched to False → error mentions `RepeatMasker`; no subprocess call.
  - `test_run_calls_repeatmasker_and_writes_gff`: subprocess patched to write a canned `.out` into the `-dir`; the final `<prefix>.te.gff3` exists and is readable by `parse_annotation`.

- [ ] **Step 2: Run:** `python -m pytest benchmark/tests/test_annotate_loci.py -v` → FAIL.

- [ ] **Step 3: Implement** `annotate_loci.py` and register the subcommand next to `classplots` in `srnaprofiler.py`. Output name collisions with existing files: refuse to overwrite and say so.

- [ ] **Step 4: Run** the file → PASS, then the full suite.

- [ ] **Step 5: Commit**
  ```bash
  git add annotate_loci.py annotation.py srnaprofiler.py benchmark/tests/test_annotate_loci.py
  git commit -m "feat(annotate-loci): RepeatMasker on locus sequences only"
  ```

---

### Task 6: Report and README

**Files:**
- Modify: `report_data.py`, `report_html.py`, `README.md`
- Test: `benchmark/tests/test_annotation_report.py`

**Interfaces:**
- Consumes: the `annotation`, `annotation_support` TSV columns (Task 4).
- Produces:
  - `report_data.Cluster` gains `annotation: str = ""` and `annotation_support: str = "none"`, filled by `_row_to_cluster` from the TSV when the columns exist (absent column → defaults; old TSVs keep working).
  - In `report_html._row`, add key `an` (the support string) only when it is not `none`, and in the per-locus detail dict key `an` (the `annotation` text) when non-empty. The data contract comment at the top of `report_html.py` documents both keys.
  - The report page shows, per locus with annotation, the text (e.g. `te:LTR/Gypsy:antisense:0.82`) and a coloured support chip (`supports`/`conflicts`/`neutral`); a legend sentence says annotation is evidence and never changes the class. Follow the existing JS/CSS patterns in `report_assets`.
  - README section "Annotation (optional)": the flag, accepted formats, the labels, the biology table, the `annotate-loci` recipe, the benchmark-identity guarantee, and the caveat that annotation is incomplete so absence means nothing.

- [ ] **Step 1: Write failing tests** in `benchmark/tests/test_annotation_report.py`:
  - `test_load_clusters_reads_annotation_columns`: a TSV with the two columns → fields filled; a TSV without them → defaults `""`/`"none"`.
  - `test_row_has_an_only_when_annotated`: `_row` output has no `an` key for `none`, has `an == "conflicts"` otherwise.
  - `test_detail_carries_annotation_text`: detail dict has the text under `an`.
  - `test_report_html_mentions_annotation_legend_only_when_used`: rendering a report from an annotated TSV contains the legend sentence; from a plain TSV it does not.

- [ ] **Step 2: Run:** `python -m pytest benchmark/tests/test_annotation_report.py -v` → FAIL.

- [ ] **Step 3: Implement** the data and HTML changes, and write the README section. The JS has never run in a real browser (see the html-report notes): after the tests pass, open one generated report in a browser and say in the report to the user what was and was not checked.

- [ ] **Step 4: Run** the file → PASS, then `python -m pytest benchmark/tests -q` → all green.

- [ ] **Step 5: Commit**
  ```bash
  git add report_data.py report_html.py README.md benchmark/tests/test_annotation_report.py
  git commit -m "feat(report): show annotation evidence per locus; document the layer"
  ```

---

## Decisions (mine, not discussed)

- Layer labels are the fixed set `gene|te|ncrna|phas` so the biology table can key on them.
- The spec's "sense/antisense to the TE consistent with ping-pong" is descriptive only (in `text`), not part of the support rule: ping-pong loci are two-stranded, so strand is a weak test.
- The spec's tRNA/rRNA/snoRNA row is descriptive only (`ncrna` context shows in `text`): the classifier has no class for them.
- `miRNA` calls come from the miRDeep2 BED, so annotation describes them but does not nudge or flag them.
- `annotate-loci` supports RepeatMasker only (`--species` or `--te-library`); minimap2/BLAST not wrapped (YAGNI).
- pybedtools was offered as installable, but the pure-Python index stays: it needs no `bedtools` binary, tests run anywhere, and the overlap logic is small. Revisit only if profiling on a large genome shows it slow.
- Task 5 also adds `repeatmasker` (bioconda) to `environment.yml` as an optional tool and mentions it in `setup`'s toolchain check as optional (warn, never fail).
- Open for later: real annotations per organism (arthropods, plants, LSV) still have to be found; the tests use synthetic files only.
