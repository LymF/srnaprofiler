"""Step 7 of the animal run: one results file, and a hand-over either way.

miRDeep2 names its results `result_DD_MM_YYYY_t_HH_MM_SS.csv`, so sorting them
by name orders by day of month, not by run. Coordinates and extra columns must
therefore come from the same file, or a row shows run A's locus with run B's
score -- the serials in `provisional id` are per run.
"""
import json
from pathlib import Path
from types import SimpleNamespace

from mirna import (deliverable_paths, finish_animal_run, novel_with_rows)

HEADER = ("provisional id\tmiRDeep2 score\testimated probability that the "
          "miRNA candidate is a true positive\trfam alert\texample miRBase "
          "miRNA with the same seed\tmiRBase miRNA\tsignificant randfold "
          "p-value\tconsensus mature sequence\tconsensus star sequence\t"
          "consensus precursor sequence\tprecursor coordinate")


def results_csv(path, score, coord):
    path.write_text(
        "novel miRNAs predicted by miRDeep2\n\n"
        + HEADER + "\n"
        + f"1_1\t{score}\t0.9\t-\t-\t-\tyes\tugagguaguagguuguauagu\t"
          f"cuauacaaucuacugucuuuc\tuggaaua\t{coord}\n\n")


def args_for(tmp_path):
    return SimpleNamespace(prefix="at", outdir=str(tmp_path), engine="mirdeep2",
                           genome=None, no_validate=True, threads=1)


def test_the_row_and_the_coordinates_come_from_the_same_results_file(tmp_path):
    """Two runs in the same directory: never one's locus with the other's score."""
    results_csv(tmp_path / "result_02_01_2026_t_10_00_00.csv", "5.5", "chr1:100..200:+")
    results_csv(tmp_path / "result_01_02_2026_t_10_00_00.csv", "99.9", "chr9:900..1000:-")

    novel, rows = novel_with_rows(tmp_path)
    assert len(novel) == 1
    score = rows["1_1"]["miRDeep2 score"]
    assert (novel[0]["chrom"], score) in (("chr1", "5.5"), ("chr9", "99.9"))


def test_no_results_file_at_all_yields_nothing_rather_than_a_stale_pick(tmp_path):
    assert novel_with_rows(tmp_path) == ([], {})


def test_zero_novel_predictions_still_hands_over_the_deliverables(tmp_path):
    """Success with no files is a lie: the empty table is the honest answer."""
    out = tmp_path / "out"
    out.mkdir()
    finish_animal_run(tmp_path / "results", args_for(out), tmp_path / "work",
                      reads=None)

    p = deliverable_paths(out, "at")
    assert p["table"].read_text().splitlines()[1:] == []
    assert p["bed"].read_text() == ""
    meta = json.loads(p["run_json"].read_text())
    assert meta["counts"]["predicted"] == 0
    assert meta["validated"] is False
