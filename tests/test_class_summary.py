from pathlib import Path
from class_summary import write_summary

TSV = "tests/data/mini.clusters.tsv"


def test_write_summary_creates_panels(tmp_path):
    paths = write_summary(TSV, str(tmp_path), "mini")
    assert paths, "expected at least one panel"
    for p in paths:
        assert Path(p).exists() and Path(p).stat().st_size > 0
    names = {Path(p).name for p in paths}
    assert "mini.summary.counts.pdf" in names
