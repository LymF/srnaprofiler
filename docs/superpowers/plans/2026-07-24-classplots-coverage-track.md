# classplots coverage-track Implementation Plan

> Executed inline (tightly-coupled change across class_signals.py / class_coverage.py / classplots.py). TDD, frequent commits.

**Goal:** Replace the per-class normalized coverage metaplot with a per-locus IGV-style mirror coverage track (sense up, antisense down), raw depth, matplotlib, for single-locus units (top-N + --region).

**Spec:** `docs/superpowers/specs/2026-07-24-classplots-coverage-track-design.md`

## Global Constraints
- English code/comments/logs; author Lucas; no AI attribution.
- Coverage track is **raw** depth (each alignment +1 per covered position, per strand) — the one documented exception to classplots' 1/NH weighting.
- No new external dependency (matplotlib only). `plotCoverageProfile.R` stays in the repo but is no longer called.

---

### Task 1: `class_signals.coverage_track`; remove `coverage_profile`

- Modify: `class_signals.py` (add `coverage_track`, remove `coverage_profile`)
- Remove: `tests/test_class_signals_coverage.py`
- Create: `tests/test_class_signals_coverage_track.py`

Add:
```python
def coverage_track(alns, region_start, region_end):
    """Per-base RAW depth over [region_start, region_end): (sense, antisense).

    IGV-style: each alignment adds +1 to every position it covers on its strand
    lane, NOT weighted by 1/NH. This is the one intentional exception to
    classplots' 1/NH weighting (a genome-browser coverage view).
    """
    length = max(region_end - region_start, 0)
    sense = [0] * length
    antisense = [0] * length
    for a in alns:
        lane = sense if a.strand == "+" else antisense
        s = max(a.pos, region_start)
        e = min(a.pos + len(a.seq), region_end)
        for p in range(s, e):
            lane[p - region_start] += 1
    return sense, antisense
```
Delete the whole `coverage_profile` function.

Test (`tests/test_class_signals_coverage_track.py`):
```python
from samio import Aln
from class_signals import coverage_track


def test_coverage_track_strand_lanes_and_clip():
    alns = [Aln("chr", 0, "+", "A" * 50, 1), Aln("chr", 60, "-", "A" * 20, 1)]
    sense, anti = coverage_track(alns, 0, 100)
    assert len(sense) == 100 and len(anti) == 100
    assert sense[0] == 1 and sense[49] == 1 and sense[50] == 0
    assert anti[60] == 1 and anti[79] == 1 and anti[80] == 0
    assert all(v == 0 for v in anti[:60])


def test_coverage_track_is_raw_and_stacks():
    # nh=2 still counts +1 (raw, not 1/NH); two overlapping + reads stack to 2
    alns = [Aln("chr", 0, "+", "A" * 10, 2), Aln("chr", 0, "+", "A" * 10, 1)]
    sense, _ = coverage_track(alns, 0, 10)
    assert sense[0] == 2


def test_coverage_track_empty_region():
    assert coverage_track([], 5, 5) == ([], [])
```

Steps: write test → run (fail import) → implement + delete old fn → delete old test file → `pytest tests/test_class_signals_coverage_track.py -q` → `pytest -q` → commit `feat(class_signals): raw per-base coverage_track; drop normalized coverage_profile`.

---

### Task 2: `class_coverage.plot_coverage_track` (matplotlib)

- Create: `class_coverage.py`
- Create: `tests/test_class_coverage.py`

```python
"""IGV-style mirror coverage track (matplotlib, native PDF, no Ghostscript)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_coverage_track(sense, antisense, chrom, start, out_path):
    """Sense up (#996631), antisense down (#0099CC), over genomic position."""
    x = [start + i for i in range(len(sense))]
    neg = [-a for a in antisense]
    m = max(max(sense, default=0), max(antisense, default=0), 1) * 1.1
    fig, ax = plt.subplots(figsize=(15, 7))
    ax.fill_between(x, sense, color="#996631", label="Sense (+)")
    ax.fill_between(x, neg, color="#0099CC", label="Antisense (-)")
    ax.set_ylim(-m, m)
    ax.set_xlabel(f"Position on {chrom}", fontsize=14)
    ax.set_ylabel("sRNA density (coverage)", fontsize=14)
    ax.set_title(f"Coverage - {chrom}:{start}-{start + len(sense)}",
                 fontsize=16, fontweight="bold")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
```

Test:
```python
from pathlib import Path
from class_coverage import plot_coverage_track


def test_plot_coverage_track_writes_pdf(tmp_path):
    out = tmp_path / "cov.pdf"
    plot_coverage_track([0, 1, 5, 2, 0], [0, 2, 3, 1, 0], "chr1", 100, str(out))
    assert Path(out).exists() and Path(out).stat().st_size > 0
```

Steps: write test → run (fail) → implement → `pytest tests/test_class_coverage.py -q` → `pytest -q` → commit `feat(class_coverage): IGV-style mirror coverage track plot`.

---

### Task 3: wire into the orchestrator; members carry chrom; drop the metaplot call

- Modify: `classplots.py`
- Modify: `tests/test_classplots_run.py` (empty-unit member becomes a 4-tuple; add a coverage assertion)

Changes in `classplots.py`:
1. Add import: `import class_coverage`.
2. Rewrite `_write_coverage_tsv`:
```python
def _write_coverage_tsv(sense, antisense, start, path):
    with open(path, "w") as fh:
        fh.write("pos\tsense\tantisense\n")
        for i, (s, a) in enumerate(zip(sense, antisense)):
            fh.write(f"{start + i}\t{s}\t{a}\n")
```
3. `_signal_plots_for` docstring: `members: list of (alns, chrom, region_start, region_end).`
4. Replace the whole coverage block (the `acc = {...}` loop + `_write_coverage_tsv` + `plotCoverageProfile.R` run) with:
```python
    # coverage track: only for a single-locus unit (top-N locus or --region),
    # where an absolute genomic-position view is meaningful. Raw depth (IGV-style).
    if len(members) == 1:
        _, chrom, rs, re = members[0]
        sense, antisense = cs.coverage_track(alns_list[0], rs, re)
        _write_coverage_tsv(sense, antisense, rs, tables / f"{prefix}.{name}.coverage.tsv")
        class_coverage.plot_coverage_track(
            sense, antisense, chrom, rs, str(plotsdir / f"{prefix}.{name}.coverage.pdf"))
```
5. In `run_classplots`, add `l.chrom`/`reg.chrom` to the three member constructions:
```python
        members = [(_fetch(fh, l, xm), l.chrom, l.start, l.end) for l in cls_loci]
        ...
            m = [(_fetch(fh, l, xm), l.chrom, l.start, l.end)]
        ...
        m = [(_fetch(fh, reg, xm), reg.chrom, reg.start, reg.end)]
```

Changes in `tests/test_classplots_run.py`:
- Empty-unit test member → `[([], "chr", 0, 100)]`.
- In the signals smoke test, additionally assert a per-locus coverage PDF is written and the pooled class has none:
```python
    assert list((out / "plots").glob("*piRNA.locus.*coverage.pdf")), "per-locus coverage expected"
    assert not (out / "plots" / "mini.piRNA.coverage.pdf").exists(), "pooled class has no coverage"
```

Steps: apply edits → `pytest tests/test_classplots_run.py -q` (signals smoke runs R; weblogo absent OK) → `pytest -q` → commit `feat(classplots): per-locus IGV coverage track; members carry chrom; drop metaplot call`.

---

### Task 4: docs
- Modify `README.md` classplots section: note coverage is a per-locus IGV-style track (top-N + --region).
- Commit `docs: note per-locus coverage track in classplots`.
