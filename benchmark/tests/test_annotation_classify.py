"""Task 4: writers, CLI flag and classify wiring for annotation evidence."""
import argparse
import sys
from pathlib import Path

import pytest

from classify import _run_params, run_classify
from records import ClusterRecord
from reporter import write_gff3, write_tsv

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _recs():
    a = ClusterRecord(chrom="chr1", start=10, end=30, strand_call="+",
                      cls="piRNA", confidence=0.8, review=False,
                      annotation="te:LTR/Gypsy:antisense:0.82|gene:exon:sense:0.10",
                      annotation_support="supports")
    b = ClusterRecord(chrom="chr1", start=50, end=70, strand_call="-",
                      cls="siRNA", confidence=0.6, review=True)
    return [a, b]


def _header(path):
    return Path(path).read_text().splitlines()[0].split("\t")


def test_tsv_unchanged_without_annotation(tmp_path):
    old, new = tmp_path / "old.tsv", tmp_path / "new.tsv"
    write_tsv(_recs(), old)
    write_tsv(_recs(), new, annotated=False)
    assert old.read_bytes() == new.read_bytes()
    assert "annotation" not in _header(new)
    assert "annotation_support" not in _header(new)


def test_gff3_unchanged_without_annotation(tmp_path):
    old, new = tmp_path / "old.gff3", tmp_path / "new.gff3"
    write_gff3(_recs(), old)
    write_gff3(_recs(), new, annotated=False)
    assert old.read_bytes() == new.read_bytes()
    assert "annotation" not in new.read_text()


def test_tsv_columns_with_annotation(tmp_path):
    p = tmp_path / "a.tsv"
    write_tsv(_recs(), p, annotated=True)
    lines = Path(p).read_text().splitlines()
    head = lines[0].split("\t")
    i = head.index("annotation")
    assert head[i + 1] == "annotation_support" and head[i + 2] == "review"
    r1 = lines[1].split("\t")
    assert r1[i] == "te:LTR/Gypsy:antisense:0.82|gene:exon:sense:0.10"
    assert r1[i + 1] == "supports"
    r2 = lines[2].split("\t")
    assert r2[i] == "" and r2[i + 1] == "none"


def test_gff3_attributes_escaped(tmp_path):
    p = tmp_path / "a.gff3"
    recs = _recs()
    recs[0].annotation = "te:a;b=c,d&e|gene:x"
    write_gff3(recs, p, annotated=True)
    attrs = Path(p).read_text().splitlines()[1].split("\t")[8]
    assert "annotation_support=supports" in attrs
    val = attrs.split(";annotation=")[1]
    assert val == "te:a%3Bb%3Dc%2Cd%26e|gene:x"
    assert ";" not in val and "=" not in val

    write_gff3(_recs(), p, annotated=True)
    val = Path(p).read_text().splitlines()[1].split("\t")[8].split(";annotation=")[1]
    assert val == "te:LTR/Gypsy:antisense:0.82|gene:exon:sense:0.10"


def _parser():
    from srnaprofiler import build_parser
    return build_parser()


def test_cli_accepts_repeatable_annotation():
    base = ["classify", "-p", "x", "-o", "o", "--sam", "s.bam"]
    args = _parser().parse_args(base + ["--annotation", "te:a.out",
                                        "--annotation", "gene:b.gff3"])
    assert args.annotation == ["te:a.out", "gene:b.gff3"]
    assert _parser().parse_args(base).annotation is None


def test_run_json_has_no_annotation_key_when_unused():
    ns = argparse.Namespace(prefix="x", annotation=None, func=lambda a: 0)
    assert "annotation" not in _run_params(ns)
    assert _run_params(ns) == {"prefix": "x"}
    ns.annotation = ["te:a.out"]
    assert _run_params(ns)["annotation"] == ["te:a.out"]


def test_invalid_annotation_label_fails_cleanly(tmp_path):
    ns = argparse.Namespace(outdir=str(tmp_path), prefix="x", config="rules.yaml",
                            sam="s.bam", annotation=["bogus:a.out"])
    with pytest.raises(ValueError, match="invalid annotation layer"):
        run_classify(ns)
