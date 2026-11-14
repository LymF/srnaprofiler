"""Task 5: annotate only the loci with RepeatMasker."""
import argparse
import os
import subprocess
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


def _read_fa(path):
    recs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            recs[name] = ""
        else:
            recs[name] += line
    return recs


def _read_map(path):
    return [l.rstrip("\n").split("\t") for l in open(path)]


def test_write_locus_fasta(tmp_path):
    fa, seq = _genome(tmp_path)
    tsv = _clusters(tmp_path, [("chr1", 100, 160), ("chr1", 200, 230)])
    out, mp = str(tmp_path / "l.fa"), str(tmp_path / "l.map.tsv")
    assert annotate_loci.write_locus_fasta(tsv, fa, out, mp) == 2
    assert _read_fa(out) == {"L1": seq[100:160], "L2": seq[200:230]}
    assert _read_map(mp) == [["L1", "chr1", "100", "160"], ["L2", "chr1", "200", "230"]]


def test_write_locus_fasta_short_ids_for_long_contig_names(tmp_path):
    seq = "ACGT" * 50
    long_name = "scaffold_" + "x" * 80
    fa = tmp_path / "g.fa"
    fa.write_text(f">{long_name}\n{seq}\n")
    pysam.faidx(str(fa))
    tsv = _clusters(tmp_path, [(long_name, 10, 60)])
    out, mp = str(tmp_path / "l.fa"), str(tmp_path / "l.map.tsv")
    annotate_loci.write_locus_fasta(tsv, str(fa), out, mp)
    assert all(len(k) <= 50 for k in _read_fa(out))
    assert _read_map(mp) == [["L1", long_name, "10", "60"]]


def test_write_locus_fasta_skips_huge(tmp_path):
    fa, _ = _genome(tmp_path)
    tsv = _clusters(tmp_path, [("chr1", 0, 50), ("chr1", 0, 1_000_001)])
    out, mp = str(tmp_path / "l.fa"), str(tmp_path / "l.map.tsv")
    assert annotate_loci.write_locus_fasta(tsv, fa, out, mp) == 1
    assert _read_map(mp) == [["L1", "chr1", "0", "50"]]


def test_rm_out_to_gff_offsets_and_strand(tmp_path):
    out = tmp_path / "x.out"
    out.write_text(
        "   SW  perc perc perc  query  position in query  matching  repeat  class/family  position in repeat\n"
        "score  div. del. ins.  sequence  begin  end  (left)  repeat  class/family  begin end (left)  ID\n\n"
        "  500  1.0  0.0  0.0  L2  11  60  (140)  C  rnd-1  LTR/Gypsy  (0)  50  1  1\n"
        "  500  1.0  0.0  0.0  L9  11  60  (140)  C  rnd-1  LTR/Gypsy  (0)  50  1  1\n")
    mp = tmp_path / "x.map.tsv"
    mp.write_text("L1\tchrA\t5\t90\nL2\tchr1\t1000\t1200\n")
    gff = str(tmp_path / "x.gff3")
    assert annotate_loci.rm_out_to_gff(str(out), str(mp), gff) == 1
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
                 "  500  1.0  0.0  0.0  L1  11  40  (20)  +  rnd-1  LINE/R1  (0)  30  1  1\n")


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


def _failing_rm(cmd):
    raise subprocess.CalledProcessError(1, cmd)


def _silent_rm(cmd):
    pass


def _leftovers(a):
    out = os.path.join(a.outdir)
    return sorted(n for n in os.listdir(out))


@pytest.mark.parametrize("rm", [_failing_rm, _silent_rm])
def test_failed_run_cleans_up_and_can_be_rerun(tmp_path, monkeypatch, rm):
    monkeypatch.setattr(annotate_loci, "which", lambda n: True)
    a = _args(tmp_path)
    monkeypatch.setattr(annotate_loci, "_run_repeatmasker", rm)
    with pytest.raises(SystemExit):
        annotate_loci.run_annotate_loci(a)
    assert _leftovers(a) == []
    monkeypatch.setattr(annotate_loci, "_run_repeatmasker", _fake_rm)
    annotate_loci.run_annotate_loci(a)
    assert "loci.te.gff3" in _leftovers(a)


def test_failure_keeps_preexisting_files(tmp_path, monkeypatch):
    monkeypatch.setattr(annotate_loci, "which", lambda n: True)
    monkeypatch.setattr(annotate_loci, "_run_repeatmasker", _failing_rm)
    a = _args(tmp_path)
    os.makedirs(os.path.join(a.outdir, "loci.repeatmasker"))
    with pytest.raises(SystemExit):
        annotate_loci.run_annotate_loci(a)
    assert _leftovers(a) == ["loci.repeatmasker"]


# ---- structural RNA (--ncrna-cm) with fake cmscan/cmpress on PATH ----
import stat

import annotate_loci as al


def _tool(path, body):
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def _tbl_row(fam, acc, q, f, t, strand, ev, inc="!", clan="-"):
    return (f"1 {fam} {acc} {clan} {q} - cm 1 100 {f} {t} {strand} no 1 0.5 0.0 "
            f"50.0 {ev} {inc} * 0 0 0 0 0 0 desc text\n")


class Nc:
    def __init__(self, tmp_path, monkeypatch, tools=("cmscan", "cmpress"), rows=""):
        self.bin = tmp_path / "bin"
        self.bin.mkdir()
        self.log = tmp_path / "calls.log"
        self.log.write_text("")
        self.tbl = tmp_path / "canned.tbl"
        self.tbl.write_text("# header\n" + rows)
        self.scan_copy = tmp_path / "scan_copy.fa"
        monkeypatch.setenv("FAKE_LOG", str(self.log))
        monkeypatch.setenv("FAKE_TBL", str(self.tbl))
        monkeypatch.setenv("FAKE_SCAN", str(self.scan_copy))
        monkeypatch.setenv("PATH", f"{self.bin}:/bin:/usr/bin")
        if "cmscan" in tools:
            _tool(self.bin / "cmscan",
                  'echo "cmscan $*" >> "$FAKE_LOG"\n'
                  'while [ $# -gt 1 ]; do [ "$1" = "--tblout" ] && out="$2"; last="$1"; shift; done\n'
                  'cp "$FAKE_TBL" "$out"\ncp "$1" "$FAKE_SCAN"\n')
        if "cmpress" in tools:
            _tool(self.bin / "cmpress", 'echo "cmpress $*" >> "$FAKE_LOG"\ntouch "$2.i1m"\n')

    def calls(self):
        return self.log.read_text().splitlines()

    def scanned(self):
        return [l[1:].strip() for l in self.scan_copy.read_text().splitlines() if l.startswith(">")]


def _nc_args(tmp_path, rows=None, cm_index=True, **kw):
    fa, _ = _genome(tmp_path)
    rows = rows or [("chr1", 0, 60, "hcsiRNA"), ("chr1", 100, 160, "miRNA"),
                    ("chr1", 200, 260, "ambiguous")]
    p = tmp_path / "c.tsv"
    p.write_text(HEADER + "".join(
        f"{c}\t{s}\t{e}\tn{i}\t{cl}\t10\t24\n" for i, (c, s, e, cl) in enumerate(rows)))
    cm = tmp_path / "Rfam.cm"
    cm.write_text("cm")
    if cm_index:
        (tmp_path / "Rfam.cm.i1m").write_text("x")
    d = dict(clusters=str(p), genome=fa, outdir=str(tmp_path / "out"), prefix="loci",
             species=None, te_library=None, threads=2, ncrna_cm=str(cm),
             ncrna_clanin=None, ncrna_classes=al.DEFAULT_NCRNA_CLASSES)
    d.update(kw)
    return argparse.Namespace(**d)


def test_ncrna_default_classes_scan_only_selected(tmp_path, monkeypatch):
    nc = Nc(tmp_path, monkeypatch)
    al.run_annotate_loci(_nc_args(tmp_path))
    assert nc.scanned() == ["L2", "L3"]


def test_ncrna_all_scans_every_locus(tmp_path, monkeypatch):
    nc = Nc(tmp_path, monkeypatch)
    al.run_annotate_loci(_nc_args(tmp_path, ncrna_classes="all"))
    assert nc.scanned() == ["L1", "L2", "L3"]


def test_cmpress_only_when_index_missing(tmp_path, monkeypatch):
    nc = Nc(tmp_path, monkeypatch)
    al.run_annotate_loci(_nc_args(tmp_path, cm_index=True))
    assert not any(c.startswith("cmpress") for c in nc.calls())
    import shutil
    shutil.rmtree(tmp_path / "out")
    os.remove(tmp_path / "Rfam.cm.i1m")
    al.run_annotate_loci(_nc_args(tmp_path, cm_index=False))
    assert any(c.startswith("cmpress -F") for c in nc.calls())


def test_cmscan_flags_and_z_value(tmp_path, monkeypatch):
    nc = Nc(tmp_path, monkeypatch)
    al.run_annotate_loci(_nc_args(tmp_path))
    cmd = [c for c in nc.calls() if c.startswith("cmscan")][0]
    # two loci of 60 bp: 2 * 0.00012 Mb < 1, so Z = 1
    assert " -Z 1 " in cmd and "--cut_ga --rfam --nohmmonly --fmt 2 --cpu 2" in cmd
    assert "--oclan" not in cmd and "-o /dev/null" in cmd and "--tblout" in cmd


def test_z_value_scales(tmp_path):
    f = tmp_path / "s.fa"
    f.write_text(">a\n" + "A" * 1_250_000 + "\n")
    assert al._z_value(str(f)) == "2.5"


def test_clanin_adds_oclan_oskip(tmp_path, monkeypatch):
    nc = Nc(tmp_path, monkeypatch)
    cl = tmp_path / "Rfam.clanin"
    cl.write_text("x")
    al.run_annotate_loci(_nc_args(tmp_path, ncrna_clanin=str(cl)))
    cmd = [c for c in nc.calls() if c.startswith("cmscan")][0]
    assert f"--clanin {cl} --oclan --oskip" in cmd


def test_coordinates_classes_and_inc_filter(tmp_path, monkeypatch):
    rows = (_tbl_row("MIR166", "RF00644", "L2", 5, 25, "+", "1e-5") +
            _tbl_row("tRNA", "RF00005", "L3", 30, 10, "-", "2e-9") +
            _tbl_row("snoR14", "RF00001", "L3", 1, 20, "+", "1e-3") +
            _tbl_row("5S_rRNA", "RF00001", "L3", 1, 20, "+", "1e-3", inc="?") +
            _tbl_row("U2", "RF00004", "L2", 40, 55, "+", "1e-4") +
            _tbl_row("5_8S_rRNA", "RF00002", "L2", 1, 5, "+", "1e-4") +
            _tbl_row("RNaseP_nuc", "RF00009", "L2", 6, 9, "+", "1e-4") +
            _tbl_row("SSU_rRNA_plant", "RF01960", "L2", 10, 12, "+", "1e-4"))
    Nc(tmp_path, monkeypatch, rows=rows)
    a = _nc_args(tmp_path)
    al.run_annotate_loci(a)
    gff = tmp_path / "out" / "loci.ncrna.gff3"
    lines = [l.rstrip("\n").split("\t") for l in open(gff) if not l.startswith("#")]
    assert len(lines) == 7                       # the inc "?" row is gone
    assert [(l[0], int(l[3])) for l in lines] == sorted((l[0], int(l[3])) for l in lines)
    by = {l[8].split(";")[0]: l for l in lines}
    mir = by["Name=MIR166"]                       # L2 starts at 100 (0-based)
    assert (mir[3], mir[4], mir[6], mir[1], mir[2]) == ("105", "125", "+", "Infernal", "ncRNA")
    assert "Accession=RF00644" in mir[8] and "Class=miRNA" in mir[8] and "evalue=1e-5" in mir[8]
    trna = by["Name=tRNA"]                        # reverse hit 30..10 on L3 (start 200)
    assert (trna[3], trna[4], trna[6]) == ("210", "230", "-") and "Class=tRNA" in trna[8]
    cls = {k: v[8].split("Class=")[1].split(";")[0] for k, v in by.items()}
    assert cls["Name=snoR14"] == "snoRNA" and cls["Name=U2"] == "snRNA"
    assert cls["Name=5_8S_rRNA"] == "rRNA" and cls["Name=SSU_rRNA_plant"] == "rRNA"
    assert cls["Name=RNaseP_nuc"] == "ncRNA"


def test_output_loads_as_ncrna_layer(tmp_path, monkeypatch):
    Nc(tmp_path, monkeypatch, rows=_tbl_row("MIR166", "RF00644", "L2", 5, 25, "+", "1e-5"))
    al.run_annotate_loci(_nc_args(tmp_path))
    from annotation import load_layers
    L = load_layers([f"ncrna:{tmp_path / 'out' / 'loci.ncrna.gff3'}"], {"chr1"})[0]
    assert L.kind == "ncrna" and L.n_features == 1
    f = parse_annotation(str(tmp_path / "out" / "loci.ncrna.gff3"))["chr1"][0]
    assert (f.start, f.end) == (104, 125)


def test_both_modes_in_one_call(tmp_path, monkeypatch):
    Nc(tmp_path, monkeypatch, rows=_tbl_row("MIR166", "RF00644", "L2", 5, 25, "+", "1e-5"))
    monkeypatch.setattr(al, "which", lambda n: True)
    monkeypatch.setattr(al, "_run_repeatmasker", _fake_rm)
    al.run_annotate_loci(_nc_args(tmp_path, species="insects"))
    out = tmp_path / "out"
    assert (out / "loci.te.gff3").exists() and (out / "loci.ncrna.gff3").exists()


def test_neither_mode_is_error(tmp_path, monkeypatch):
    Nc(tmp_path, monkeypatch)
    with pytest.raises(SystemExit):
        al.run_annotate_loci(_nc_args(tmp_path, ncrna_cm=None))


def test_both_te_options_still_error_with_ncrna(tmp_path, monkeypatch):
    Nc(tmp_path, monkeypatch)
    monkeypatch.setattr(al, "which", lambda n: True)
    with pytest.raises(SystemExit):
        al.run_annotate_loci(_nc_args(tmp_path, species="x", te_library="l.fa"))


def test_existing_ncrna_output_blocks_rerun(tmp_path, monkeypatch, capsys):
    Nc(tmp_path, monkeypatch)
    a = _nc_args(tmp_path)
    os.makedirs(a.outdir)
    open(os.path.join(a.outdir, "loci.ncrna.gff3"), "w").write("keep")
    with pytest.raises(SystemExit):
        al.run_annotate_loci(a)
    assert "exists" in capsys.readouterr().err
    assert open(os.path.join(a.outdir, "loci.ncrna.gff3")).read() == "keep"


@pytest.mark.parametrize("tools,name", [(("cmpress",), "cmscan"), (("cmscan",), "cmpress")])
def test_missing_infernal_tool_is_error(tmp_path, monkeypatch, capsys, tools, name):
    Nc(tmp_path, monkeypatch, tools=tools)
    with pytest.raises(SystemExit):
        al.run_annotate_loci(_nc_args(tmp_path, cm_index=(name == "cmscan") and False))
    err = capsys.readouterr().err
    assert name in err and "infernal" in err


def test_failed_cmscan_cleans_scratch(tmp_path, monkeypatch):
    nc = Nc(tmp_path, monkeypatch)
    _tool(nc.bin / "cmscan", "exit 1\n")
    a = _nc_args(tmp_path)
    with pytest.raises(SystemExit):
        al.run_annotate_loci(a)
    assert _leftovers(a) == []


def test_no_locus_of_selected_classes_is_error(tmp_path, monkeypatch):
    Nc(tmp_path, monkeypatch)
    a = _nc_args(tmp_path, rows=[("chr1", 0, 60, "hcsiRNA")])
    with pytest.raises(SystemExit):
        al.run_annotate_loci(a)
    assert _leftovers(a) == []


def test_setup_lists_infernal_tools():
    import setup_module
    names = {r["name"]: r for r in setup_module.MODULE_TOOLS["annotate-loci"]}
    assert "optional" in names["cmscan"]["note"] and "optional" in names["cmpress"]["note"]
