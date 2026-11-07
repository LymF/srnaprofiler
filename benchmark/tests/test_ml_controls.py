import gzip
import math

import pandas as pd
import pytest

from benchmark.ml import dataset
from benchmark.ml.controls import is_control
from benchmark.ml.replicates import repro_features
from benchmark.scripts import build_controls as bc


def _read_bed(path):
    return [tuple(l.rstrip("\n").split("\t")) for l in open(path) if l.strip()]


GFF = (
    "##gff-version 3\n"
    "c1\tx\tgene\t1\t5000\t.\t+\t.\tID=g1\n"
    "c1\tx\tmRNA\t1\t5000\t.\t+\t.\tID=t1;Parent=g1\n"
    "c1\tx\texon\t101\t300\t.\t+\t.\tID=e1;Parent=t1\n"
    "c1\tx\texon\t1001\t1200\t.\t+\t.\tID=e2;Parent=t1\n"
    "c1\tx\texon\t2001\t2100\t.\t+\t.\tID=e3;Parent=t2,t1\n"
    "c1\tx\tlnc_RNA\t3001\t4000\t.\t+\t.\tID=t3\n"
    "c1\tx\texon\t3001\t3500\t.\t+\t.\tID=e4;Parent=t3\n"
    "c1\tx\trRNA\t4501\t4800\t.\t-\t.\tID=r1\n"
    "ctgA\tx\texon\t11\t60\t.\t+\t.\tID=e5;Parent=t1\n"
)


def test_controls_drop_truth_overlap(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(GFF)
    truth = tmp_path / "t.bed"
    truth.write_text("c1\t1100\t1150\tmiR1\n")  # inside exon e2 (0-based 1000-1200)
    out = tmp_path / "controls.bed"
    bc.main([str(gff), str(truth), "-o", str(out)])
    rows = _read_bed(out)
    coords = {(r[0], int(r[1]), int(r[2]), r[3]) for r in rows}
    assert ("c1", 1000, 1200, "exon") not in coords
    assert ("c1", 100, 300, "exon") in coords
    assert ("c1", 2000, 2100, "exon") in coords
    assert ("c1", 4500, 4800, "rRNA") in coords
    assert not any(r[1:3] == ("3000", "3500") for r in rows)  # non-mRNA exon excluded


def test_controls_gz_aliases_and_gff_truth(tmp_path):
    gff = tmp_path / "a.gff3.gz"
    with gzip.open(gff, "wt") as fh:
        fh.write(GFF)
    truth = tmp_path / "t.gff3"
    truth.write_text("c1\ts\tmiRNA_primary_transcript\t4501\t4600\t.\t+\t.\tID=m\n")
    alias = tmp_path / "al.tsv"
    alias.write_text("ctgA\tCM0001.1\n")
    out = tmp_path / "controls.bed"
    bc.main([str(gff), str(truth), "--aliases", str(alias), "-o", str(out)])
    rows = _read_bed(out)
    assert ("CM0001.1", "10", "60", "exon") in rows
    assert not any(r[3] == "rRNA" for r in rows)  # rRNA overlaps GFF truth


def _table(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def test_is_control_rule(tmp_path):
    bed = tmp_path / "c.bed"
    bed.write_text("chr1\t900\t1050\texon\n")
    t = _table([("chr1", 1000, 1100),   # 50 of 100 in the exon -> control
                ("chr1", 1010, 1110),   # 40 of 100 -> not
                ("chr2", 1000, 1100)])
    assert list(is_control(t, bed)) == [True, False, False]


def _rep(rows):
    return pd.DataFrame(rows, columns=["acc", "organism", "tissue", "chrom", "start", "end",
                                       "reads_rpm"])


def test_repro_three_replicates():
    df = _rep([("A", "o", "t", "c", 1000, 1100, 10.0),
               ("B", "o", "t", "c", 1000, 1100, 20.0),
               ("C", "o", "t", "c", 5000, 5100, 5.0),
               ("A", "o", "t", "c", 9000, 9100, 4.0)])
    out = repro_features(df)
    assert out.loc[0, "repro_frac"] == 0.5
    assert out.loc[0, "repro_rpm_ratio"] == 2.0
    assert out.loc[3, "repro_frac"] == 0.0 and math.isnan(out.loc[3, "repro_rpm_ratio"])
    assert out.loc[1, "repro_frac"] == 0.5


def test_repro_single_replicate_nan():
    out = repro_features(_rep([("A", "o", "t", "c", 1, 100, 1.0)]))
    assert out.isna().all().all() and list(out.columns) == ["repro_frac", "repro_rpm_ratio"]


def test_repro_groups_do_not_mix():
    df = _rep([("A", "o", "t1", "c", 1000, 1100, 10.0),
               ("B", "o", "t2", "c", 1000, 1100, 10.0)])
    assert repro_features(df)["repro_frac"].isna().all()


def test_build_table_fills_controls_and_repro(tmp_path, monkeypatch):
    from benchmark.tests.test_ml_dataset import MANIFEST_HEADER, _write_clusters
    from benchmark.truth import Region
    manifest = tmp_path / "m.tsv"
    with open(manifest, "w") as fh:
        fh.write("\t".join(MANIFEST_HEADER) + "\n")
        for acc, rep in (("T1", "1"), ("T2", "2")):
            fh.write("\t".join([acc, "Drosophila melanogaster", "tuning", "arthropod", "ovary",
                                rep, "u", "x.gff"]) + "\n")
    res = tmp_path / "res"
    _write_clusters(res / "T1" / "work" / "T1.clusters.tsv",
                    [("chr1", 1000, 1100, "piRNA"), ("chr1", 5000, 5100, "siRNA")])
    _write_clusters(res / "T2" / "work" / "T2.clusters.tsv", [("chr1", 1000, 1100, "piRNA")])
    monkeypatch.setattr(dataset, "load_truth", lambda lib: [Region("chr1", 0, 10, "piRNA", "s", "r")])
    tdir = tmp_path / "truth"
    (tdir / "dmel").mkdir(parents=True)
    (tdir / "dmel" / "controls.bed").write_text("chr1\t900\t1200\texon\n")
    monkeypatch.setattr(dataset, "TRUTH_DIR", tdir)
    df = dataset.build_table(["T1", "T2"], res, manifest)
    assert list(df.is_control) == [True, False, True]
    assert list(df.repro_frac) == [1.0, 0.0, 1.0]
    (tdir / "dmel" / "controls.bed").unlink()
    df = dataset.build_table(["T1", "T2"], res, manifest)
    assert not df.is_control.any()
