"""Task 5: annotate only the loci with RepeatMasker."""
import argparse
import os
import sys

import pysam
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import annotate_loci
from annotation import parse_annotation

HEADER = "chrom\tstart\tend\tname\tcls\treads\tdicer_call\n"


def _genome(tmp_path):
    seq = "".join("ACGT"[i % 4] for i in range(400))
    fa = tmp_path / "g.fa"
    fa.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(fa))
    return str(fa), seq


def _clusters(tmp_path, rows):
    p = tmp_path / "c.tsv"
    p.write_text(HEADER + "".join(
        f"{c}\t{s}\t{e}\tn{i}\tpiRNA\t10\t24\n" for i, (c, s, e) in enumerate(rows)))
    return str(p)


def _args(tmp_path, **kw):
    fa, _ = _genome(tmp_path)
    d = dict(clusters=_clusters(tmp_path, [("chr1", 100, 160)]), genome=fa,
             outdir=str(tmp_path / "out"), prefix="loci", species="insects",
             te_library=None, threads=2)
    d.update(kw)
    return argparse.Namespace(**d)


def test_write_locus_fasta(tmp_path):
    fa, seq = _genome(tmp_path)
    tsv = _clusters(tmp_path, [("chr1", 100, 160), ("chr1", 200, 230)])
    out = str(tmp_path / "l.fa")
    assert annotate_loci.write_locus_fasta(tsv, fa, out) == 2
    recs = {}
    name = None
    for line in open(out):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            recs[name] = ""
        else:
            recs[name] += line
    assert recs == {"chr1:100-160": seq[100:160], "chr1:200-230": seq[200:230]}


def test_write_locus_fasta_skips_huge(tmp_path):
    fa, _ = _genome(tmp_path)
    tsv = _clusters(tmp_path, [("chr1", 0, 50), ("chr1", 0, 1_000_001)])
    assert annotate_loci.write_locus_fasta(tsv, fa, str(tmp_path / "l.fa")) == 1


def test_rm_out_to_gff_offsets_and_strand(tmp_path):
    out = tmp_path / "x.out"
    out.write_text(
        "   SW  perc perc perc  query  position in query  matching  repeat  class/family  position in repeat\n"
        "score  div. del. ins.  sequence  begin  end  (left)  repeat  class/family  begin end (left)  ID\n\n"
        "  500  1.0  0.0  0.0  chr1:1000-1200  11  60  (140)  C  rnd-1  LTR/Gypsy  (0)  50  1  1\n")
    gff = str(tmp_path / "x.gff3")
    assert annotate_loci.rm_out_to_gff(str(out), gff) == 1
    cols = [l for l in open(gff) if not l.startswith("#")][0].rstrip("\n").split("\t")
    assert cols[0] == "chr1" and cols[2] == "repeat_region"
    assert (cols[3], cols[4], cols[6]) == ("1011", "1060", "-")
    assert "Class=LTR/Gypsy" in cols[8] and "Name=rnd-1" in cols[8]


def test_run_requires_one_of_species_or_library(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(annotate_loci, "which", lambda n: True)
    for kw in (dict(species=None, te_library=None),
               dict(species="x", te_library="lib.fa")):
        with pytest.raises(SystemExit):
            annotate_loci.run_annotate_loci(_args(tmp_path, **kw))


def test_run_fails_clearly_without_repeatmasker(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(annotate_loci, "which", lambda n: False)
    calls = []
    monkeypatch.setattr(annotate_loci, "_run_repeatmasker", lambda *a, **k: calls.append(a))
    with pytest.raises(SystemExit):
        annotate_loci.run_annotate_loci(_args(tmp_path))
    assert "RepeatMasker" in capsys.readouterr().err
    assert not calls


def _fake_rm(cmd):
    d = cmd[cmd.index("-dir") + 1]
    fa = cmd[-1]
    with open(os.path.join(d, os.path.basename(fa) + ".out"), "w") as fh:
        fh.write("score div del ins seq begin end left rep class\n\n"
                 "  500  1.0  0.0  0.0  chr1:100-160  11  40  (20)  +  rnd-1  LINE/R1  (0)  30  1  1\n")


def test_run_calls_repeatmasker_and_writes_gff(tmp_path, monkeypatch):
    monkeypatch.setattr(annotate_loci, "which", lambda n: True)
    seen = []
    monkeypatch.setattr(annotate_loci, "_run_repeatmasker",
                        lambda cmd: (seen.append(cmd), _fake_rm(cmd)))
    annotate_loci.run_annotate_loci(_args(tmp_path))
    cmd = seen[0]
    assert cmd[0] == "RepeatMasker" and "-species" in cmd and "-pa" in cmd
    gff = tmp_path / "out" / "loci.te.gff3"
    feats = parse_annotation(str(gff))
    f = feats["chr1"][0]
    assert (f.start, f.end, f.name) == (110, 140, "LINE/R1")


def test_run_uses_library_flag(tmp_path, monkeypatch):
    monkeypatch.setattr(annotate_loci, "which", lambda n: True)
    seen = []
    monkeypatch.setattr(annotate_loci, "_run_repeatmasker",
                        lambda cmd: (seen.append(cmd), _fake_rm(cmd)))
    annotate_loci.run_annotate_loci(_args(tmp_path, species=None, te_library="lib.fa"))
    assert "-lib" in seen[0] and "-species" not in seen[0]


def test_run_refuses_to_overwrite(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(annotate_loci, "which", lambda n: True)
    monkeypatch.setattr(annotate_loci, "_run_repeatmasker", _fake_rm)
    a = _args(tmp_path)
    os.makedirs(a.outdir)
    open(os.path.join(a.outdir, "loci.te.gff3"), "w").write("keep")
    with pytest.raises(SystemExit):
        annotate_loci.run_annotate_loci(a)
    assert "exists" in capsys.readouterr().err
    assert open(os.path.join(a.outdir, "loci.te.gff3")).read() == "keep"


def test_setup_reports_repeatmasker_optional():
    import setup_module
    reqs = setup_module.MODULE_TOOLS["annotate-loci"]
    rm = [r for r in reqs if r["name"] == "RepeatMasker"][0]
    assert "optional" in rm["note"]
    assert setup_module.summarise({"annotate-loci": [
        {"ok": False, "note": rm["note"]}]})["annotate-loci"] is True
