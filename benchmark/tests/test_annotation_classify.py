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


def test_gff3_escape_percent_first_and_whitespace(tmp_path):
    from reporter import _gff3_escape
    assert _gff3_escape("a%b;c") == "a%25b%3Bc"
    assert _gff3_escape("a\tb\nc") == "a%09b%0Ac"
    p = tmp_path / "a.gff3"
    recs = _recs()
    recs[0].annotation = "x\ty\n%z"
    write_gff3(recs, p, annotated=True)
    lines = p.read_text().splitlines()
    assert len(lines) == 3 and len(lines[1].split("\t")) == 9
    assert lines[1].endswith(";annotation=x%09y%0A%25z")


def _run_main(monkeypatch, argv):
    import srnaprofiler
    monkeypatch.setattr(srnaprofiler, "check_dependencies", lambda: True)
    monkeypatch.setattr(sys, "argv", ["srnaprofiler"] + argv)
    srnaprofiler.main()


def test_invalid_annotation_label_exits_nonzero_with_error(tmp_path, monkeypatch, capsys):
    with pytest.raises(SystemExit) as e:
        _run_main(monkeypatch, ["classify", "-p", "x", "-o", str(tmp_path),
                                "--sam", "s.bam", "--annotation", "bogus:a.out"])
    assert e.value.code == 1
    err = capsys.readouterr().err
    assert "[ERROR]" in err and "invalid annotation layer" in err


def test_missing_annotation_file_exits_nonzero_with_error(tmp_path, monkeypatch, capsys):
    with pytest.raises(SystemExit) as e:
        _run_main(monkeypatch, ["classify", "-p", "x", "-o", str(tmp_path),
                                "--sam", "s.bam",
                                "--annotation", f"te:{tmp_path}/nope.out"])
    assert e.value.code == 1
    err = capsys.readouterr().err
    assert "[ERROR]" in err and "annotation file not found" in err


# ---- end to end: real run_classify on a tiny synthetic BAM -------------------

def _build_inputs(tmp_path, chroms, loci):
    """loci: list of (chrom, start). 30 reads of 24 nt per locus, spread over ~64 bp."""
    import pysam
    sys.path.insert(0, str(REPO / "benchmark" / "tests"))
    import bamutil
    reads = []
    for chrom, start in loci:
        for k in range(30):
            reads.append((chrom, start + (k % 40), "ACGT" * 6, "+", 1))
    bam = tmp_path / "in.bam"
    bamutil.make_seq_bam(bam, reads, chroms=chroms)
    fa = tmp_path / "g.fa"
    import random
    rng = random.Random(1)
    fa.write_text("".join(
        f">{n}\n{''.join(rng.choice('ACGT') for _ in range(ln))}\n" for n, ln in chroms))
    pysam.faidx(str(fa))
    return bam, fa


def _args(tmp_path, bam, fa, outdir, **kw):
    argv = ["classify", "-p", "run", "-o", str(outdir), "--sam", str(bam),
            "-g", str(fa), "--norm", "1000000", "--mincov", "1", "--threads", "1",
            "--config", str(REPO / "rules.yaml")]
    for k, v in kw.items():
        if k == "annotation":
            for spec in v:
                argv += ["--annotation", spec]
        else:
            argv += [f"--{k.replace('_', '-')}", v]
    return _parser().parse_args(argv)


def _table(path):
    lines = Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    return head, [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


@pytest.fixture
def e2e(tmp_path):
    chroms = (("chr1", 5000), ("chrM", 5000))
    bam, fa = _build_inputs(tmp_path, chroms, [("chr1", 1000), ("chrM", 2000)])
    te = tmp_path / "te.out"
    te.write_text(
        "   SW  perc perc perc  query  position in query  matching  repeat  position in repeat\n"
        "score  div. del. ins.  sequence  begin end (left)  repeat class/family begin end (left) ID\n\n"
        "  500  1.0  0.0  0.0  chr1  900  1100  (3900)  +  rep1  LTR/Gypsy  1  200  (0)  1\n"
        "  500  1.0  0.0  0.0  chrM  1900  2100  (2900)  +  rep2  LTR/Gypsy  1  200  (0)  2\n")
    return tmp_path, bam, fa, te


def test_e2e_without_annotation_has_no_annotation_output(e2e):
    tmp_path, bam, fa, te = e2e
    out = tmp_path / "o1"
    run_classify(_args(tmp_path, bam, fa, out))
    head, rows = _table(out / "run.clusters.tsv")
    assert rows, "synthetic BAM must yield loci"
    assert "annotation" not in head and "annotation_support" not in head
    assert "annotation" not in (out / "run.clusters.gff3").read_text()
    import json
    assert "annotation" not in json.loads((out / "run.run.json").read_text())


def test_e2e_with_annotation_columns_and_unchanged_cls(e2e):
    tmp_path, bam, fa, te = e2e
    base, ann = tmp_path / "base", tmp_path / "ann"
    run_classify(_args(tmp_path, bam, fa, base))
    run_classify(_args(tmp_path, bam, fa, ann, annotation=[f"te:{te}"]))
    h0, r0 = _table(base / "run.clusters.tsv")
    h1, r1 = _table(ann / "run.clusters.tsv")
    i = h1.index("annotation")
    assert h1[i + 1] == "annotation_support" and h1[i + 2] == "review"
    assert [r["cls"] for r in r1] == [r["cls"] for r in r0]
    assert any(r["annotation"] for r in r1)
    assert "annotation_support=" in (ann / "run.clusters.gff3").read_text()
    import json
    assert json.loads((ann / "run.run.json").read_text())["annotation"] == [f"te:{te}"]


def test_e2e_organelle_contigs_still_pass_all_bam_references(e2e, monkeypatch):
    import classify
    tmp_path, bam, fa, te = e2e
    seen = {}
    real = classify.load_layers

    def spy(specs, known):
        seen["known"] = set(known)
        return real(specs, known)

    monkeypatch.setattr(classify, "load_layers", spy)
    out = tmp_path / "o"
    run_classify(_args(tmp_path, bam, fa, out, annotation=[f"te:{te}"],
                       organelle_contigs="chrM"))
    assert seen["known"] == {"chr1", "chrM"}
    _, rows = _table(out / "run.clusters.tsv")
    assert {r["chrom"] for r in rows} == {"chr1"}


def test_e2e_contig_name_mismatch_warns_once_and_is_neutral(tmp_path, monkeypatch, capsys):
    bam, fa = _build_inputs(tmp_path, (("1", 5000), ("2", 5000)), [("1", 1000)])
    te = tmp_path / "te.out"
    te.write_text(
        "score  div. del. ins.  sequence  begin end (left)  repeat class/family begin end (left) ID\n\n"
        "  500  1.0  0.0  0.0  chr1  900  1100  (3900)  +  rep1  LTR/Gypsy  1  200  (0)  1\n")
    out = tmp_path / "o"
    run_classify(_args(tmp_path, bam, fa, out, annotation=[f"te:{te}"]))
    err = capsys.readouterr().err
    assert err.count("[WARN]") == 1 and "not in the reference" in err
    _, rows = _table(out / "run.clusters.tsv")
    assert rows
    # No contig of the layer matches the BAM: the code yields "neutral" (absent
    # features are never evidence) with an empty summary, not "none".
    assert {r["annotation_support"] for r in rows} == {"neutral"}
    assert {r["annotation"] for r in rows} == {""}
