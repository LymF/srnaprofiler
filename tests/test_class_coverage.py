from pathlib import Path
from class_coverage import plot_coverage_track


def test_plot_coverage_track_writes_pdf(tmp_path):
    out = tmp_path / "cov.pdf"
    plot_coverage_track([0, 1, 5, 2, 0], [0, 2, 3, 1, 0], "chr1", 100, str(out))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_plot_coverage_track_all_zero_no_crash(tmp_path):
    out = tmp_path / "z.pdf"
    plot_coverage_track([0, 0, 0], [0, 0, 0], "chr1", 0, str(out))
    assert Path(out).exists()
