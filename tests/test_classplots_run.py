import shutil
import types

import pysam
import pytest

from classplots import run_classplots


def _make_bam(path):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 2000}],
              "PG": [{"ID": "bowtie", "PN": "bowtie"}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for pos in range(100, 380, 21):
            a = pysam.AlignedSegment()
            a.query_name = f"r{pos}"; a.query_sequence = "T" + "A" * 23
            a.flag = 0; a.reference_id = 0; a.reference_start = pos
            a.mapping_quality = 255
            a.cigar = [(0, 24)]
            a.set_tag("XM", 1)
            bam.write(a)
    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


def _args(**kw):
    base = dict(clusters="tests/data/mini.clusters.tsv", bam=None, outdir=None,
               prefix="mini", top_n=1, classes=None, region=None, min_length=15,
               max_length=35, logo_trim=None, no_signals=True, no_summary=False)
    base.update(kw)
    return types.SimpleNamespace(**base)


def test_run_classplots_summary_only(tmp_path):
    bam = tmp_path / "reads.bam"; _make_bam(bam)
    out = tmp_path / "out"
    run_classplots(_args(bam=str(bam), outdir=str(out)))
    pdfs = list((out / "plots").glob("mini.summary.*.pdf"))
    assert pdfs, "summary panels should be written even with --no-signals"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript absent")
def test_run_classplots_signals_smoke(tmp_path):
    bam = tmp_path / "reads.bam"; _make_bam(bam)
    out = tmp_path / "out"
    run_classplots(_args(bam=str(bam), outdir=str(out), no_signals=False, no_summary=True, top_n=1))
    # the piRNA locus (chr1:100-400) has reads -> its signal plots should be produced
    assert list((out / "plots").glob("*piRNA*.pingpong.pdf")), "expected piRNA ping-pong plot"
    assert list((out / "tables").glob("*piRNA*.sizedist.tsv")), "expected piRNA sizedist table"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript absent")
def test_run_classplots_topn_filenames_sanitized(tmp_path):
    bam = tmp_path / "reads.bam"; _make_bam(bam)
    out = tmp_path / "out"
    run_classplots(_args(bam=str(bam), outdir=str(out), no_signals=False, no_summary=True, top_n=1))
    produced = list((out / "plots").glob("*locus*")) + list((out / "tables").glob("*locus*"))
    assert produced, "expected per-locus outputs"
    assert all(":" not in p.name for p in produced), "locus filenames must not contain ':'"
