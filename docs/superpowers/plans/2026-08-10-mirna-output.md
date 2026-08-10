# miRNA Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `mirna` subcommand deliver one row per prediction with its verdict and evidence, the three validated FASTAs, the structure plots, and a layout that matches the rest of the tool.

**Architecture:** A new `mirna_report.py` owns the deliverables — building records and writing them — mirroring what `reporter.py` already does for `classify`. `mirna.py` keeps orchestration and third-party interop and shrinks by handing the writing off. The join happens in one place: `build_records` takes the engine's predictions, miRScore's rows, and (for `mirdeep2`) the raw `result_*.csv` rows, and returns one record per prediction.

**Tech Stack:** Python 3.13, stdlib `csv`/`json`/`dataclasses`, `pytest`. No new dependency.

## Global Constraints

- **Code, comments, docstrings, logs and help text in English.** Comments short — only what is needed; long context belongs in the commit message.
- **Author is Lucas.** Never add a `Co-Authored-By: Claude` trailer to a commit, and never credit AI in any document.
- **Coordinates in the TSV and BED are 0-based half-open.** The engines report 1-based inclusive; conversion happens once, in `build_records`.
- **Absent value is the string `"NA"`**, matching `reporter._cell`.
- **`result` is never empty.** Its four values are `Pass`, `Fail`, `NotScored`, `NotValidated`.
- **Fail loud, never silent-zero.** A required miRScore column that is missing raises; an optional passthrough column that is missing becomes `"NA"`.
- Run tests with `pytest tests/ -q` from the repo root, inside the `srnaprofiler` conda env.

---

### Task 1: Read miRScore's results as full rows

`parse_mirscore_results` returns only the set of passing names. The table needs every column, and a missing column must fail rather than become `NA`.

**Files:**
- Modify: `mirna.py:961-968` (`parse_mirscore_results`)
- Test: `tests/test_mirna_report_rows.py` (create)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `read_mirscore_rows(path) -> dict[str, dict]`, keyed by the `name` column, each value the raw row with miRScore's own header names. `MIRSCORE_REQUIRED` is the tuple of columns the table depends on. `parse_mirscore_results(path) -> set[str]` keeps its signature and is reimplemented on top.

- [ ] **Step 1: Write the failing test**

```python
"""miRScore's results read as rows, with a missing column failing loud."""
import pytest

from mirna import read_mirscore_rows, parse_mirscore_results

HEADER = ("name,mSeq,mLen,mStart,mStop,msSeq,msLen,msStart,msStop,precSeq,"
          "precLen,result,flags,mismatches,bulges,mReads,msReads,totReads,"
          "numLibraries,precision\n")
ROWS = ("ath-1,TGGACTGAAGGGAGCTCCCTT,21,15,35,GGGGAATGAAGCCTGGTCCGA,21,60,80,"
        "ACGT,120,Pass,NA,1,0,8214,112,8500,1,0.98\n"
        "ath-2,TTGACAGAAGATAGAGAGCAC,21,20,40,NA,NA,NA,NA,ACGT,110,Fail,"
        "Star not detected,3,1,4401,0,4500,1,0.91\n")


def write_csv(tmp_path, text):
    p = tmp_path / "miRScore_results.csv"
    p.write_text(text)
    return p


def test_rows_are_keyed_by_name_and_keep_every_column(tmp_path):
    rows = read_mirscore_rows(write_csv(tmp_path, HEADER + ROWS))
    assert set(rows) == {"ath-1", "ath-2"}
    assert rows["ath-1"]["msSeq"] == "GGGGAATGAAGCCTGGTCCGA"
    assert rows["ath-2"]["flags"] == "Star not detected"


def test_a_missing_required_column_fails_loud(tmp_path):
    """Silently writing NA is how a column rename becomes a wrong table."""
    header = HEADER.replace(",msReads", ",starReads")
    with pytest.raises(RuntimeError, match="msReads"):
        read_mirscore_rows(write_csv(tmp_path, header + ROWS))


def test_passing_names_still_come_out_of_the_same_read(tmp_path):
    assert parse_mirscore_results(write_csv(tmp_path, HEADER + ROWS)) == {"ath-1"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_report_rows.py -q`
Expected: FAIL with `ImportError: cannot import name 'read_mirscore_rows'`

- [ ] **Step 3: Write minimal implementation**

Replace `parse_mirscore_results` in `mirna.py` with:

```python
MIRSCORE_REQUIRED = ("name", "mSeq", "msSeq", "msLen", "precSeq", "precLen",
                     "mStart", "mStop", "msStart", "msStop", "result", "flags",
                     "mismatches", "bulges", "mReads", "msReads", "totReads",
                     "precision")


def read_mirscore_rows(path):
    """miRScore's results keyed by name, with every column kept.

    A missing column raises: writing NA for a renamed column would put a wrong
    table in front of the user with no sign anything went wrong.
    """
    import csv as _csv
    with open(path) as fh:
        reader = _csv.DictReader(fh)
        missing = [c for c in MIRSCORE_REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise RuntimeError(
                f"{path} lacks the column(s) {', '.join(missing)}; miRScore's "
                f"output format changed and the table cannot be built from it")
        return {r["name"]: dict(r) for r in reader}


def parse_mirscore_results(path):
    """Names miRScore passed. Anything not an explicit pass is a fail."""
    return {n for n, r in read_mirscore_rows(path).items()
            if str(r.get("result", "")).strip().lower().startswith("pass")}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mirna_report_rows.py -q`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add mirna.py tests/test_mirna_report_rows.py
git commit -m "feat(mirna): read miRScore's results as rows, failing on a renamed column"
```

---

### Task 2: Keep miRDeep2's raw novel rows

The `mirdeep2`-only columns (score, rfam alert, randfold, miRBase hit) live in the same CSV that `parse_mirdeep2_novel` already walks, and are thrown away.

**Files:**
- Modify: `mirna.py:897-920` (`parse_mirdeep2_novel`)
- Test: `tests/test_mirna_report_rows.py` (append)

**Interfaces:**
- Consumes: nothing.
- Produces: `mirdeep2_novel_rows(path) -> list[dict]` — the novel section's raw rows with miRDeep2's own header names. `parse_mirdeep2_novel(path) -> list[dict]` keeps its signature and its output shape (`chrom, strand, name, start, end, mature`) and is rebuilt on top.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_mirna_report_rows.py`:

```python
from mirna import mirdeep2_novel_rows, parse_mirdeep2_novel

NOVEL = "\n".join([
    "some preamble line",
    "\t".join(["provisional id", "miRDeep2 score", "rfam alert",
               "significant randfold p-value", "miRBase miRNA",
               "example miRBase miRNA with the same seed",
               "consensus mature sequence", "precursor coordinate"]),
    "\t".join(["nov-1", "96.3", "-", "yes", "-", "ath-miR166a",
               "uggacugaagggagcucccuu", "Chr1:1234..1325"]),
    "\t".join(["nov-2", "10.0", "rRNA/tRNA", "no", "-", "-",
               "uugacagaagauagagagcac", "Chr2:88100..88190"]),
    "",
]) + "\n"


def test_raw_novel_rows_keep_the_mirdeep_only_columns(tmp_path):
    p = tmp_path / "result_01_01_2026.csv"
    p.write_text(NOVEL)
    rows = mirdeep2_novel_rows(p)
    assert [r["provisional id"] for r in rows] == ["nov-1", "nov-2"]
    assert rows[0]["example miRBase miRNA with the same seed"] == "ath-miR166a"
    assert rows[1]["rfam alert"] == "rRNA/tRNA"


def test_predictions_are_unchanged_by_the_refactor(tmp_path):
    p = tmp_path / "result_01_01_2026.csv"
    p.write_text(NOVEL)
    preds = parse_mirdeep2_novel(p)
    assert preds[0] == {"chrom": "Chr1", "strand": ".", "name": "nov-1",
                        "start": 1234, "end": 1325,
                        "mature": "TGGACTGAAGGGAGCTCCCTT"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_report_rows.py -q`
Expected: FAIL with `ImportError: cannot import name 'mirdeep2_novel_rows'`

- [ ] **Step 3: Write minimal implementation**

Replace `parse_mirdeep2_novel` in `mirna.py` with:

```python
def mirdeep2_novel_rows(path):
    """The novel section of miRDeep2's result CSV, rows as written."""
    import csv as _csv
    rows = list(_csv.reader(open(path), delimiter="\t"))
    start = None
    for i, r in enumerate(rows):
        if r and r[0] == "provisional id":
            hdr, start = r, i + 1
            break
    if start is None:
        return []
    out = []
    for r in rows[start:]:
        if not r or not r[0].strip():
            break
        out.append(dict(zip(hdr, r)))
    return out


def parse_mirdeep2_novel(path):
    """miRDeep2's novel predictions, in the shape the engines share."""
    out = []
    for d in mirdeep2_novel_rows(path):
        m = re.match(r"(\S+):(\d+)\.\.(\d+)", d.get("precursor coordinate", ""))
        if not m:
            continue
        out.append({"chrom": m.group(1), "strand": ".", "name": d["provisional id"],
                    "start": int(m.group(2)), "end": int(m.group(3)),
                    "mature": d["consensus mature sequence"].strip().upper().replace("U", "T")})
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mirna_report_rows.py -q`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add mirna.py tests/test_mirna_report_rows.py
git commit -m "refactor(mirna): keep miRDeep2's raw novel rows, predictions built on top"
```

---

### Task 3: The record, and the join

The core of the spec: one record per prediction, four `result` values, coordinates converted once, name collisions refused.

**Files:**
- Create: `mirna_report.py`
- Test: `tests/test_mirna_report_build.py` (create)

**Interfaces:**
- Consumes: the prediction dicts (`chrom, strand, name, start, end, mature`, 1-based inclusive) that `parse_mirdeep_p2_predictions` and `parse_mirdeep2_novel` return; the row dicts from Task 1 and Task 2.
- Produces:
  - `MirnaRecord` dataclass with the fields listed in `MIRNA_COLS`.
  - `MIRNA_COLS: list[str]` — the TSV column order.
  - `build_records(predictions, engine, mirscore_rows=None, submitted=None, mirdeep_rows=None) -> list[MirnaRecord]`. `mirscore_rows=None` means miRScore did not run. `submitted` is the set of names handed to miRScore; `None` is treated as "all of them".

- [ ] **Step 1: Write the failing test**

```python
"""One record per prediction: the verdict, the evidence, and the coordinates."""
import pytest

from mirna_report import build_records, MirnaRecord, MIRNA_COLS

PREDS = [
    {"chrom": "Chr1", "strand": "+", "name": "p1", "start": 1234, "end": 1325,
     "mature": "TGGACTGAAGGGAGCTCCCTT"},
    {"chrom": "Chr2", "strand": "-", "name": "p2", "start": 88100, "end": 88190,
     "mature": "TTGACAGAAGATAGAGAGCAC"},
    {"chrom": "Chr3", "strand": "+", "name": "p3", "start": 551020, "end": 551108,
     "mature": "TAAGGGATCTTACAGGTAACG"},
]

SCORE_ROWS = {
    "p1": {"name": "p1", "mSeq": "TGGACTGAAGGGAGCTCCCTT", "msSeq": "GGGGAATGAAGCC",
           "msLen": "13", "precSeq": "ACGTACGT", "precLen": "120",
           "mStart": "15", "mStop": "35", "msStart": "60", "msStop": "80",
           "result": "Pass", "flags": "NA", "mismatches": "1", "bulges": "0",
           "mReads": "8214", "msReads": "112", "totReads": "8500",
           "precision": "0.98"},
    "p2": {"name": "p2", "mSeq": "TTGACAGAAGATAGAGAGCAC", "msSeq": "NA",
           "msLen": "NA", "precSeq": "ACGT", "precLen": "110",
           "mStart": "20", "mStop": "40", "msStart": "NA", "msStop": "NA",
           "result": "Fail", "flags": "Star not detected", "mismatches": "3",
           "bulges": "1", "mReads": "4401", "msReads": "0", "totReads": "4500",
           "precision": "0.91"},
}


def test_one_record_per_prediction():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert [r.name for r in recs] == ["p1", "p2", "p3"]


def test_coordinates_become_zero_based_half_open():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert (recs[0].start, recs[0].end) == (1233, 1325)


def test_the_four_verdicts():
    """p1 passed, p2 failed, p3 was never submitted, and nothing is blank."""
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert [r.result for r in recs] == ["Pass", "Fail", "NotScored"]
    assert recs[2].flags == "precursor could not be re-cut from the genome"

    none_ran = build_records(PREDS, "mirdeep-p2", mirscore_rows=None)
    assert {r.result for r in none_ran} == {"NotValidated"}
    assert all(r.star_seq == "NA" for r in none_ran)


def test_submitted_but_absent_from_the_results_is_not_a_failure():
    """miRScore dropping a locus is not miRScore rejecting it."""
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS,
                         submitted={"p1", "p2", "p3"})
    assert recs[2].result == "NotScored"
    assert recs[2].flags == "absent from miRScore results"


def test_evidence_columns_come_from_mirscore():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert recs[0].star_seq == "GGGGAATGAAGCC"
    assert recs[0].star_len == "13"
    assert recs[0].mature_reads == "8214"
    assert recs[0].precision == "0.98"
    assert recs[0].mature_len == 21          # ours, from the prediction


def test_mirdeep_only_columns_are_na_for_the_plant_engine():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert recs[0].mirdeep_score == "NA"
    assert recs[0].mirbase_hit == "NA"
    assert recs[0].engine == "mirdeep-p2"


def test_mirdeep_only_columns_are_filled_when_the_engine_provides_them():
    rows = {"p1": {"provisional id": "p1", "miRDeep2 score": "96.3",
                   "rfam alert": "-", "significant randfold p-value": "yes",
                   "miRBase miRNA": "-",
                   "example miRBase miRNA with the same seed": "ath-miR166a"}}
    recs = build_records(PREDS, "mirdeep2", SCORE_ROWS, submitted={"p1", "p2"},
                         mirdeep_rows=rows)
    assert recs[0].mirdeep_score == "96.3"
    assert recs[0].mirbase_same_seed == "ath-miR166a"
    assert recs[1].mirdeep_score == "NA"     # absent row, not an error


def test_a_repeated_name_fails_loud():
    """The join is by name; a duplicate would silently multiply rows."""
    dup = PREDS + [dict(PREDS[0])]
    with pytest.raises(RuntimeError, match="p1"):
        build_records(dup, "mirdeep-p2", SCORE_ROWS, submitted={"p1"})


def test_every_column_exists_on_the_record():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    for col in MIRNA_COLS:
        assert hasattr(recs[0], col), col
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_report_build.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'mirna_report'`

- [ ] **Step 3: Write minimal implementation**

Create `mirna_report.py`:

```python
"""The mirna module's deliverables: one record per prediction, and the files.

The engine knows where a locus is and miRScore knows whether it is a miRNA;
neither knows both. This is where the two are joined.
"""
from dataclasses import dataclass, field

NA = "NA"

MIRNA_COLS = [
    "name", "chrom", "start", "end", "strand", "engine",
    "mature_seq", "mature_len", "star_seq", "star_len",
    "precursor_seq", "precursor_len",
    "mature_start", "mature_end", "star_start", "star_end",
    "result", "flags", "mismatches", "bulges",
    "mature_reads", "star_reads", "total_reads", "precision",
    "mirdeep_score", "mirdeep_prob", "rfam_alert", "randfold",
    "mirbase_hit", "mirbase_same_seed",
]

# miRScore column -> our column. Ours are snake_case and say what they hold.
_FROM_MIRSCORE = {
    "msSeq": "star_seq", "msLen": "star_len",
    "precSeq": "precursor_seq", "precLen": "precursor_len",
    "mStart": "mature_start", "mStop": "mature_end",
    "msStart": "star_start", "msStop": "star_end",
    "flags": "flags", "mismatches": "mismatches", "bulges": "bulges",
    "mReads": "mature_reads", "msReads": "star_reads",
    "totReads": "total_reads", "precision": "precision",
}

# miRDeep2 column -> our column. Passthrough: absent means NA, not an error.
_FROM_MIRDEEP = {
    "miRDeep2 score": "mirdeep_score",
    "estimated probability that the miRNA candidate is a true positive": "mirdeep_prob",
    "rfam alert": "rfam_alert",
    "significant randfold p-value": "randfold",
    "miRBase miRNA": "mirbase_hit",
    "example miRBase miRNA with the same seed": "mirbase_same_seed",
}


@dataclass
class MirnaRecord:
    """One prediction, with wherever it got to."""
    name: str
    chrom: str
    start: int                  # 0-based, half-open, like .clusters.tsv
    end: int
    strand: str
    engine: str
    mature_seq: str
    mature_len: int
    result: str                 # Pass | Fail | NotScored | NotValidated
    star_seq: str = NA
    star_len: str = NA
    precursor_seq: str = NA
    precursor_len: str = NA
    mature_start: str = NA
    mature_end: str = NA
    star_start: str = NA
    star_end: str = NA
    flags: str = NA
    mismatches: str = NA
    bulges: str = NA
    mature_reads: str = NA
    star_reads: str = NA
    total_reads: str = NA
    precision: str = NA
    mirdeep_score: str = NA
    mirdeep_prob: str = NA
    rfam_alert: str = NA
    randfold: str = NA
    mirbase_hit: str = NA
    mirbase_same_seed: str = NA


def _verdict(row):
    """miRScore's own result, normalised to Pass/Fail."""
    return "Pass" if str(row.get("result", "")).strip().lower().startswith("pass") else "Fail"


def build_records(predictions, engine, mirscore_rows=None, submitted=None,
                  mirdeep_rows=None):
    """One record per prediction, joined by name.

    `mirscore_rows=None` means miRScore did not run at all -- which is a
    different thing from a locus it rejected, so it gets its own verdict.
    """
    seen = set()
    for p in predictions:
        if p["name"] in seen:
            raise RuntimeError(
                f"prediction name {p['name']} appears more than once; the "
                f"table joins miRScore by name and would multiply rows")
        seen.add(p["name"])

    by_name = dict(mirdeep_rows or {})

    out = []
    for p in predictions:
        rec = MirnaRecord(
            name=p["name"], chrom=p["chrom"],
            start=p["start"] - 1, end=p["end"],       # 1-based inclusive in
            strand=p.get("strand", "."), engine=engine,
            mature_seq=p["mature"], mature_len=len(p["mature"]),
            result="NotValidated",
        )
        if mirscore_rows is not None:
            row = mirscore_rows.get(p["name"])
            if row is not None:
                rec.result = _verdict(row)
                for src, dst in _FROM_MIRSCORE.items():
                    setattr(rec, dst, row.get(src, NA) or NA)
            elif submitted is not None and p["name"] in submitted:
                rec.result = "NotScored"
                rec.flags = "absent from miRScore results"
            else:
                rec.result = "NotScored"
                rec.flags = "precursor could not be re-cut from the genome"
        for src, dst in _FROM_MIRDEEP.items():
            val = by_name.get(p["name"], {}).get(src)
            if val not in (None, ""):
                setattr(rec, dst, val)
        out.append(rec)
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mirna_report_build.py -q`
Expected: 9 passed

- [ ] **Step 5: Commit**

```bash
git add mirna_report.py tests/test_mirna_report_build.py
git commit -m "feat(mirna): one record per prediction, joining the engine and miRScore"
```

---

### Task 4: The writers, and the guarantee that the four files agree

**Files:**
- Modify: `mirna_report.py`
- Test: `tests/test_mirna_report_write.py` (create)

**Interfaces:**
- Consumes: `MirnaRecord`, `MIRNA_COLS`, `build_records` from Task 3.
- Produces:
  - `write_table(records, path) -> None`
  - `write_bed(records, path) -> int` — validated only, returns the count.
  - `write_fasta(records, path, field) -> int` — validated only; `field` is one of `"mature_seq"`, `"star_seq"`, `"precursor_seq"`.
  - `validated(records) -> list[MirnaRecord]`

- [ ] **Step 1: Write the failing test**

```python
"""The table has everything; the BED and the FASTAs have exactly the passes."""
from mirna_report import (build_records, write_table, write_bed, write_fasta,
                          validated, MIRNA_COLS)

PREDS = [
    {"chrom": "Chr1", "strand": "+", "name": "p1", "start": 1234, "end": 1325,
     "mature": "TGGACTGAAGGGAGCTCCCTT"},
    {"chrom": "Chr2", "strand": "-", "name": "p2", "start": 88100, "end": 88190,
     "mature": "TTGACAGAAGATAGAGAGCAC"},
]
ROWS = {
    "p1": {"result": "Pass", "msSeq": "GGGGAATGAAGCC", "precSeq": "ACGTACGT",
           "flags": "NA"},
    "p2": {"result": "Fail", "msSeq": "NA", "precSeq": "ACGT",
           "flags": "Star not detected"},
}


def recs():
    return build_records(PREDS, "mirdeep-p2", ROWS, submitted={"p1", "p2"})


def test_the_table_has_a_row_per_prediction_and_the_declared_header(tmp_path):
    p = tmp_path / "at.mirna.tsv"
    write_table(recs(), p)
    lines = p.read_text().splitlines()
    assert lines[0].split("\t") == MIRNA_COLS
    assert len(lines) == 3
    assert lines[1].split("\t")[MIRNA_COLS.index("result")] == "Pass"


def test_a_tab_in_a_field_cannot_break_a_row(tmp_path):
    r = recs()
    r[1].flags = "star\tnot\tdetected"
    p = tmp_path / "at.mirna.tsv"
    write_table(r, p)
    for line in p.read_text().splitlines():
        assert len(line.split("\t")) == len(MIRNA_COLS)


def test_the_bed_carries_only_the_passes_in_half_open_coordinates(tmp_path):
    p = tmp_path / "at.mirna_loci.bed"
    assert write_bed(recs(), p) == 1
    assert p.read_text() == "Chr1\t1233\t1325\tp1\n"


def test_the_three_fastas_and_the_bed_are_the_same_set(tmp_path):
    r = recs()
    names = {rec.name for rec in validated(r)}
    for field in ("mature_seq", "star_seq", "precursor_seq"):
        path = tmp_path / f"at.{field}.fa"
        assert write_fasta(r, path, field) == len(names)
        got = {ln[1:] for ln in path.read_text().splitlines() if ln.startswith(">")}
        assert got == names


def test_the_fasta_holds_the_sequence_of_the_requested_field(tmp_path):
    p = tmp_path / "at.star.fa"
    write_fasta(recs(), p, "star_seq")
    assert p.read_text() == ">p1\nGGGGAATGAAGCC\n"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_report_write.py -q`
Expected: FAIL with `ImportError: cannot import name 'write_table'`

- [ ] **Step 3: Write minimal implementation**

Append to `mirna_report.py`:

```python
def _cell(value):
    """Serialises a field, neutralising TAB/newline that would break the row."""
    if value is None:
        return NA
    return str(value).replace("\t", " ").replace("\n", " ")


def write_table(records, path):
    with open(path, "w") as fh:
        fh.write("\t".join(MIRNA_COLS) + "\n")
        for r in records:
            fh.write("\t".join(_cell(getattr(r, c)) for c in MIRNA_COLS) + "\n")


def validated(records):
    """The records miRScore passed -- the set every other file agrees on."""
    return [r for r in records if r.result == "Pass"]


def write_bed(records, path):
    """The validated precursors, ready for `classify --mirna-loci`."""
    n = 0
    with open(path, "w") as fh:
        for r in validated(records):
            fh.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.name}\n")
            n += 1
    return n


def write_fasta(records, path, field):
    """One validated sequence per record, named as in the BED."""
    if field not in ("mature_seq", "star_seq", "precursor_seq"):
        raise ValueError(f"not a sequence field: {field}")
    n = 0
    with open(path, "w") as fh:
        for r in validated(records):
            fh.write(f">{r.name}\n{getattr(r, field)}\n")
            n += 1
    return n
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mirna_report_write.py -q`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add mirna_report.py tests/test_mirna_report_write.py
git commit -m "feat(mirna): the table, the BED and the three FASTAs, one set between them"
```

---

### Task 5: Promote the structure plots

They already exist, four levels down, produced by miRScore for both engines.

**Files:**
- Modify: `mirna_report.py`
- Test: `tests/test_mirna_report_structures.py` (create)

**Interfaces:**
- Consumes: `MirnaRecord` from Task 3.
- Produces: `collect_structures(mirscore_out, records, dest) -> int` — copies the plot for each scored record into `dest/passed/` or `dest/failed/`, preferring `strucVis/` (structure with read depth) over `RNAplots/`, and returns how many were copied.

- [ ] **Step 1: Write the failing test**

```python
"""The plots miRScore already makes, moved where someone will find them."""
from pathlib import Path

from mirna_report import build_records, collect_structures

PREDS = [
    {"chrom": "Chr1", "strand": "+", "name": "p1", "start": 10, "end": 100,
     "mature": "TGGACTGAAG"},
    {"chrom": "Chr1", "strand": "+", "name": "p2", "start": 200, "end": 300,
     "mature": "TTGACAGAAG"},
    {"chrom": "Chr1", "strand": "+", "name": "p3", "start": 400, "end": 500,
     "mature": "TAAGGGATCT"},
]
ROWS = {"p1": {"result": "Pass"}, "p2": {"result": "Fail"}}


def recs():
    return build_records(PREDS, "mirdeep-p2", ROWS, submitted={"p1", "p2"})


def make_mirscore_out(tmp_path, subdir, ext):
    out = tmp_path / "mirscore"
    for verdict in ("passed", "failed"):
        (out / subdir / verdict).mkdir(parents=True)
    (out / subdir / "passed" / f"p1{ext}").write_text("plot p1")
    (out / subdir / "failed" / f"p2{ext}").write_text("plot p2")
    return out


def test_strucvis_plots_are_split_by_verdict(tmp_path):
    out = make_mirscore_out(tmp_path, "strucVis", ".ps")
    dest = tmp_path / "at.mirna_structures"
    assert collect_structures(out, recs(), dest) == 2
    assert (dest / "passed" / "p1.ps").read_text() == "plot p1"
    assert (dest / "failed" / "p2.ps").read_text() == "plot p2"


def test_rnaplots_are_the_fallback_when_strucvis_is_absent(tmp_path):
    out = make_mirscore_out(tmp_path, "RNAplots", "_ss.eps")
    dest = tmp_path / "at.mirna_structures"
    assert collect_structures(out, recs(), dest) == 2
    assert (dest / "passed" / "p1_ss.eps").exists()


def test_strucvis_wins_when_both_exist(tmp_path):
    out = make_mirscore_out(tmp_path, "strucVis", ".ps")
    for verdict in ("passed", "failed"):
        (out / "RNAplots" / verdict).mkdir(parents=True)
    (out / "RNAplots" / "passed" / "p1_ss.eps").write_text("worse plot")
    dest = tmp_path / "at.mirna_structures"
    collect_structures(out, recs(), dest)
    assert (dest / "passed" / "p1.ps").exists()
    assert not (dest / "passed" / "p1_ss.eps").exists()


def test_no_plot_source_is_a_warning_not_a_crash(tmp_path):
    dest = tmp_path / "at.mirna_structures"
    assert collect_structures(tmp_path / "mirscore", recs(), dest) == 0


def test_a_never_scored_prediction_has_no_plot(tmp_path):
    """p3 never reached miRScore, so there is nothing to copy for it."""
    out = make_mirscore_out(tmp_path, "strucVis", ".ps")
    dest = tmp_path / "at.mirna_structures"
    collect_structures(out, recs(), dest)
    assert not list(dest.rglob("p3*"))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_report_structures.py -q`
Expected: FAIL with `ImportError: cannot import name 'collect_structures'`

- [ ] **Step 3: Write minimal implementation**

Append to `mirna_report.py` (add `import shutil` and `from pathlib import Path` and `from utils import log` at the top):

```python
def collect_structures(mirscore_out, records, dest):
    """miRScore's plots, split by verdict, where the user will look.

    strucVis draws read depth on the fold and is the one worth looking at;
    RNAplots is the plain structure and stands in when strucVis did not run.
    """
    mirscore_out, dest = Path(mirscore_out), Path(dest)
    source = next((s for s in ("strucVis", "RNAplots")
                   if (mirscore_out / s).is_dir()), None)
    if source is None:
        log(f"no structure plots under {mirscore_out}: neither strucVis/ nor "
            f"RNAplots/. Is strucVis on PATH?", "WARN")
        return 0

    verdicts = {"Pass": "passed", "Fail": "failed"}
    n = 0
    for r in records:
        sub = verdicts.get(r.result)
        if sub is None:                      # NotScored/NotValidated: no plot exists
            continue
        for plot in sorted((mirscore_out / source / sub).glob(f"{r.name}.*")) or \
                    sorted((mirscore_out / source / sub).glob(f"{r.name}_*")):
            (dest / sub).mkdir(parents=True, exist_ok=True)
            shutil.copy2(plot, dest / sub / plot.name)
            n += 1
    return n
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mirna_report_structures.py -q`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add mirna_report.py tests/test_mirna_report_structures.py
git commit -m "feat(mirna): the structure plots become a deliverable, split by verdict"
```

---

### Task 6: miRScore's validation returns rows, not just names

`validate_with_mirscore` returns a set of passing names, which cannot tell a rejection from a locus that was never submitted.

**Files:**
- Modify: `mirna.py:1032-1085` (`validate_with_mirscore`)
- Test: `tests/test_mirna_validation_result.py` (create)

**Interfaces:**
- Consumes: `read_mirscore_rows` from Task 1.
- Produces: `Validation` dataclass with `rows: dict` and `submitted: set`. `validate_with_mirscore(...) -> Validation | None` — `None` still means "did not run".

- [ ] **Step 1: Write the failing test**

```python
"""Validation carries what was submitted, not only what passed."""
from mirna import Validation


def test_validation_holds_rows_and_the_submitted_set():
    v = Validation(rows={"p1": {"result": "Pass"}}, submitted={"p1", "p2"})
    assert v.rows["p1"]["result"] == "Pass"
    assert "p2" in v.submitted


def test_a_submitted_name_missing_from_rows_is_visible():
    """This is the difference between 'rejected' and 'never got there'."""
    v = Validation(rows={"p1": {"result": "Pass"}}, submitted={"p1", "p2"})
    assert v.submitted - set(v.rows) == {"p2"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_validation_result.py -q`
Expected: FAIL with `ImportError: cannot import name 'Validation'`

- [ ] **Step 3: Write minimal implementation**

In `mirna.py`, add near the other dataclasses/imports:

```python
from dataclasses import dataclass


@dataclass
class Validation:
    """What miRScore said, and about which predictions it was asked."""
    rows: dict          # name -> its row in miRScore_results.csv
    submitted: set      # names handed to miRScore at all
```

In `validate_with_mirscore`, collect the submitted names in the loop that writes the FASTAs and return a `Validation`:

```python
    kept = set()
    with open(mature_fa, "w") as mf, open(hairpin_fa, "w") as hf:
        for p in predictions:
            hp = hairpin_with_flank(genome, p["chrom"], p["start"], p["end"],
                                    p["mature"], flank=flank)
            if hp is None:
                continue
            mf.write(f">{p['name']}\n{p['mature']}\n")
            hf.write(f">{p['name']}\n{hp}\n")
            kept.add(p["name"])
    if not kept:
        log("No prediction could be re-cut against the genome; skipping "
            "validation", "WARN")
        return None
    log(f"{len(kept)}/{len(predictions)} predictions prepared for miRScore", "INFO")
```

and at the end, in place of the `passing` return:

```python
    rows = read_mirscore_rows(results)
    passing = {n for n, r in rows.items()
               if str(r.get("result", "")).strip().lower().startswith("pass")}
    if not passing and no_reads_counted(out / "reads.csv"):
        raise RuntimeError(
            f"miRScore counted zero reads at every locus ({out/'reads.csv'}); "
            "that is its silent failure, not a biological result -- check that "
            "the reads it was given are the library the predictions came from")
    log(f"miRScore validated {len(passing)}/{len(kept)}", "INFO")
    return Validation(rows=rows, submitted=kept)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mirna_validation_result.py -q && pytest tests/ -q`
Expected: 2 passed, then the whole suite green

- [ ] **Step 5: Commit**

```bash
git add mirna.py tests/test_mirna_validation_result.py
git commit -m "refactor(mirna): validation returns miRScore's rows and what was submitted"
```

---

### Task 7: `finish_mirna` writes every deliverable

**Files:**
- Modify: `mirna.py:487-513` (`finish_mirna`), and delete `write_mirna_loci_bed` (`mirna.py:979-991`), now superseded by `mirna_report.write_bed`
- Test: `tests/test_mirna_deliverables.py` (create)

**Interfaces:**
- Consumes: `build_records`, `write_table`, `write_bed`, `write_fasta`, `collect_structures`, `validated` (Tasks 3-5); `Validation` (Task 6).
- Produces: `finish_mirna(predictions, args, work_dir, reads, mirdeep_rows=None) -> Path` — writes every deliverable under `args.outdir` and returns the BED path. `deliverable_paths(outdir, prefix) -> dict[str, Path]` names them in one place.

- [ ] **Step 1: Write the failing test**

```python
"""Every deliverable, named the same way whichever engine ran."""
from pathlib import Path
from types import SimpleNamespace

from mirna import finish_mirna, deliverable_paths

PREDS = [
    {"chrom": "Chr1", "strand": "+", "name": "p1", "start": 1234, "end": 1325,
     "mature": "TGGACTGAAGGGAGCTCCCTT"},
    {"chrom": "Chr2", "strand": "-", "name": "p2", "start": 88100, "end": 88190,
     "mature": "TTGACAGAAGATAGAGAGCAC"},
]


def args_for(tmp_path):
    return SimpleNamespace(prefix="at", outdir=str(tmp_path), engine="mirdeep-p2",
                           genome=None, no_validate=True, threads=1)


def test_the_deliverables_are_named_off_the_prefix(tmp_path):
    p = deliverable_paths(tmp_path, "at")
    assert p["table"].name == "at.mirna.tsv"
    assert p["bed"].name == "at.mirna_loci.bed"
    assert p["mature"].name == "at.mirna.mature.fa"
    assert p["star"].name == "at.mirna.star.fa"
    assert p["precursor"].name == "at.mirna.precursor.fa"
    assert p["run_json"].name == "at.mirna.run.json"
    assert p["structures"].name == "at.mirna_structures"


def test_without_validation_the_table_is_written_and_no_fasta_is(tmp_path):
    """--no-validate: every row NotValidated, so there is no validated set."""
    finish_mirna(PREDS, args_for(tmp_path), tmp_path / "work", reads=None)
    p = deliverable_paths(tmp_path, "at")
    lines = p["table"].read_text().splitlines()
    assert len(lines) == 3
    assert all("NotValidated" in ln for ln in lines[1:])
    assert not p["mature"].exists()
    assert not p["star"].exists()
    assert not p["precursor"].exists()


def test_the_bed_is_still_written_unvalidated_so_the_run_can_go_on(tmp_path):
    bed = finish_mirna(PREDS, args_for(tmp_path), tmp_path / "work", reads=None)
    assert bed == deliverable_paths(tmp_path, "at")["bed"]
    assert len(bed.read_text().splitlines()) == 2


def test_the_run_json_records_what_ran(tmp_path):
    import json
    finish_mirna(PREDS, args_for(tmp_path), tmp_path / "work", reads=None)
    meta = json.loads(deliverable_paths(tmp_path, "at")["run_json"].read_text())
    assert meta["engine"] == "mirdeep-p2"
    assert meta["validated"] is False
    assert meta["counts"]["predicted"] == 2
    assert meta["counts"]["passed"] == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_deliverables.py -q`
Expected: FAIL with `ImportError: cannot import name 'deliverable_paths'`

- [ ] **Step 3: Write minimal implementation**

In `mirna.py`, delete `write_mirna_loci_bed` and replace `finish_mirna` with:

```python
def deliverable_paths(outdir, prefix):
    """Every file this module hands over, named off the prefix in one place."""
    outdir = Path(outdir)
    return {
        "table": outdir / f"{prefix}.mirna.tsv",
        "bed": outdir / f"{prefix}.mirna_loci.bed",
        "mature": outdir / f"{prefix}.mirna.mature.fa",
        "star": outdir / f"{prefix}.mirna.star.fa",
        "precursor": outdir / f"{prefix}.mirna.precursor.fa",
        "run_json": outdir / f"{prefix}.mirna.run.json",
        "structures": outdir / f"{prefix}.mirna_structures",
    }


def finish_mirna(predictions, args, work_dir, reads, mirdeep_rows=None):
    """Validate, write everything the module hands over, and say what happened."""
    import json
    from mirna_report import (build_records, write_table, write_bed,
                              write_fasta, collect_structures, validated)

    kingdom = "plant" if args.engine == "mirdeep-p2" else "animal"
    result = None
    if getattr(args, "no_validate", False):
        log("--no-validate: writing predictions without miRScore", "WARN")
    elif getattr(args, "genome", None):
        result = validate_with_mirscore(
            predictions, args.genome, reads, kingdom,
            Path(work_dir) / "validation", threads=getattr(args, "threads", 1))
    else:
        log("no -g/--genome: skipping miRScore validation, which needs the "
            "genome to re-cut precursors with flanks", "WARN")

    records = build_records(
        predictions, args.engine,
        mirscore_rows=result.rows if result else None,
        submitted=result.submitted if result else None,
        mirdeep_rows=mirdeep_rows)

    paths = deliverable_paths(args.outdir, args.prefix)
    write_table(records, paths["table"])
    n_bed = write_bed(records, paths["bed"]) if result else _write_bed_unvalidated(
        records, paths["bed"])
    if validated(records):
        for key, field in (("mature", "mature_seq"), ("star", "star_seq"),
                           ("precursor", "precursor_seq")):
            write_fasta(records, paths[key], field)
    if result:
        collect_structures(Path(work_dir) / "validation" / "mirscore", records,
                           paths["structures"])

    counts = {"predicted": len(records),
              "submitted": len(result.submitted) if result else 0,
              "passed": sum(1 for r in records if r.result == "Pass"),
              "failed": sum(1 for r in records if r.result == "Fail"),
              "not_scored": sum(1 for r in records if r.result == "NotScored")}
    paths["run_json"].write_text(json.dumps(
        {"engine": args.engine, "kingdom": kingdom,
         "validated": result is not None, "counts": counts}, indent=2))

    if result:
        log(f"{n_bed} validated loci written to {paths['bed']} "
            f"(of {len(records)} predicted)", "INFO")
    else:
        log(f"{n_bed} UNVALIDATED loci written to {paths['bed']}", "WARN")
    log(f"table: {paths['table']}", "INFO")
    return paths["bed"]


def _write_bed_unvalidated(records, path):
    """Without miRScore there is no pass set, so every prediction is written.

    A false prediction becomes a carved locus with miRNA asserted over it --
    which is what the WARN in the caller is about.
    """
    with open(path, "w") as fh:
        for r in records:
            fh.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.name}\n")
    return len(records)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mirna_deliverables.py -q`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add mirna.py tests/test_mirna_deliverables.py
git commit -m "feat(mirna): the module hands over a table, three FASTAs and the plots"
```

---

### Task 8: One layout, whichever engine ran

The last divergence: `mirdeep2` builds `mirna/processed/` + `mirna/results/`, `mirdeep-p2` builds `mirna/mirdp2/`.

**Files:**
- Modify: `mirna.py:515-523` and `mirna.py:549-553` (the two `run_mirna_analysis` branches), `mirna.py:697-700` (the call into `finish_mirna`)
- Test: `tests/test_mirna_layout.py` (create)

**Interfaces:**
- Consumes: `deliverable_paths` from Task 7.
- Produces: `work_paths(outdir, prefix) -> dict[str, Path]` with keys `root`, `engine`, `validation`.

- [ ] **Step 1: Write the failing test**

```python
"""The work tree stops depending on --engine."""
from mirna import work_paths


def test_the_work_tree_is_named_off_the_prefix(tmp_path):
    w = work_paths(tmp_path, "at")
    assert w["root"] == tmp_path / "at.mirna_work"
    assert w["engine"] == tmp_path / "at.mirna_work" / "engine"
    assert w["validation"] == tmp_path / "at.mirna_work" / "validation"


def test_the_engine_directory_does_not_carry_the_engine_name(tmp_path):
    """Whoever switches --engine finds the raw output in the same place.

    work_paths takes no engine argument, which is the point: there is no
    branch left that could name the directory after the tool that ran.
    """
    assert work_paths(tmp_path, "at")["engine"].name == "engine"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mirna_layout.py -q`
Expected: FAIL with `ImportError: cannot import name 'work_paths'`

- [ ] **Step 3: Write minimal implementation**

Add to `mirna.py`:

```python
def work_paths(outdir, prefix):
    """Where the raw third-party output lands, same tree for every engine."""
    root = Path(outdir) / f"{prefix}.mirna_work"
    return {"root": root, "engine": root / "engine",
            "validation": root / "validation"}
```

In `run_mirna_analysis`, replace the `mirdeep-p2` branch's directory setup:

```python
    if getattr(args, "engine", "mirdeep2") == "mirdeep-p2":
        outdir = ensure_outdir(Path(args.outdir))
        work = work_paths(outdir, args.prefix)
        preds = run_mirdeep_p2(args, args.reads, ensure_outdir(work["engine"]))
        finish_mirna(preds, args, work["root"], args.reads)
        log("miRNA analysis completed successfully!", "INFO")
        return
```

and, in the `mirdeep2` branch, replace the `mirna_dir`/`processed_dir`/`results_dir` setup:

```python
    outdir = ensure_outdir(Path(args.outdir))
    work = work_paths(outdir, args.prefix)
    processed_dir = ensure_outdir(work["engine"] / "processed")
    results_dir = ensure_outdir(work["engine"] / "results")
```

and the call at the end of that branch:

```python
    if novel:
        log(f"Step 7: validating {len(novel)} novel predictions...", "INFO")
        rows = {r["provisional id"]: r
                for csv_path in sorted(results_dir.glob("result_*.csv"))
                for r in mirdeep2_novel_rows(csv_path)}
        finish_mirna(novel, args, work["root"], input_files["reads"],
                     mirdeep_rows=rows)
    else:
        log("no novel predictions to validate", "INFO")
```

- [ ] **Step 4: Run the whole suite**

Run: `pytest tests/ -q`
Expected: every test passes. Any test that asserted a path under `outdir/mirna/` must be updated to `deliverable_paths`/`work_paths` — grep for it first:

```bash
grep -rn "mirna/" tests/ | grep -v "\.pyc"
```

- [ ] **Step 5: Commit**

```bash
git add mirna.py tests/test_mirna_layout.py
git commit -m "refactor(mirna): one work tree and one set of names, whichever engine ran"
```

---

### Task 9: The run on real data

The acceptance criteria are numbers the ROADMAP already fixes: 87 predictions, 56 validated, 56/56 in miRBase.

**Files:**
- Modify: `docs/ROADMAP.md` (record the result), `README.md` (the new output layout)

- [ ] **Step 1: Run the module end to end**

Use the Arabidopsis library with ncRNA-filtered reads, the same inputs as the clean-install run recorded in the ROADMAP:

```bash
python3 srnaprofiler.py mirna --engine mirdeep-p2 \
    -p at -o <outdir> -g <TAIR10.fa> --reads <at.ncrna_filtered.fastq> \
    --species ath --threads 8
```

- [ ] **Step 2: Check the seven acceptance criteria**

```bash
cd <outdir>
test $(($(wc -l < at.mirna.tsv) - 1)) -eq 87 && echo "rows ok"
awk -F'\t' 'NR>1 && $17=="Pass"' at.mirna.tsv | wc -l          # expect 56
wc -l < at.mirna_loci.bed                                       # expect 56
for f in mature star precursor; do grep -c '^>' at.mirna.$f.fa; done   # 56 each
diff <(cut -f4 at.mirna_loci.bed | sort) \
     <(awk -F'\t' 'NR>1 && $17=="Pass" {print $1}' at.mirna.tsv | sort)
ls at.mirna_structures/passed | wc -l
python3 -c "import json;print(json.load(open('at.mirna.run.json'))['counts'])"
```

Expected: 87 rows, 56 passes, 56 BED lines, 56 records in each FASTA, an empty `diff`, and `counts` agreeing with all of it.

- [ ] **Step 3: Confirm the biology did not move**

```bash
python3 extras/validate_mirna.py --bed at.mirna_loci.bed --mirbase <ath.gff3>
```

Expected: 56/56 in miRBase, as recorded for the clean install.

- [ ] **Step 4: Record the numbers and the layout**

Add to `docs/ROADMAP.md`, under the miRNA section, the measured counts and the new output tree. Update the `mirna` section of `README.md` with the deliverable names.

- [ ] **Step 5: Commit**

```bash
git add docs/ROADMAP.md README.md
git commit -m "docs: what the mirna module now hands over, measured on Arabidopsis"
```

---

## What this plan does not do

Each is out of scope by an explicit decision in the spec, repeated here so nobody adds it while passing through:

- **isomiR profiling** — dedicated tools and a standard format (mirtop) exist.
- **MFE/AMFE/MFEI** — miRScore computes MFE and discards it; its criteria are the duplex geometry.
- **Our own miRBase matching** — miRDeep2's is a seed 2-8 lookup, an animal-centric notion. The miRBase columns are passthrough.
- **Cleaning the work tree** and converting `.ps` to PDF.
