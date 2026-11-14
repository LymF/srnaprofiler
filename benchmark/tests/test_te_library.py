"""te-library: RepeatModeler2 + contaminant filter, with fake tools on PATH."""
import argparse
import os
import stat
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import te_library


def _script(path, body):
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


class Env:
    def __init__(self, tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler",
                                                     "minimap2", "famdb.py"),
                 sys_path=True):
        self.tmp = tmp_path
        self.bin = tmp_path / "bin"
        self.bin.mkdir()
        self.log = tmp_path / "calls.log"
        self.log.write_text("")
        self.paf = tmp_path / "fake.paf"
        self.paf.write_text("")
        self.fams = tmp_path / "fake_fams.fa"
        self.dfam = tmp_path / "fake_dfam.fa"
        self.dfam.write_text("")
        env = {"FAKE_LOG": str(self.log), "FAKE_PAF": str(self.paf),
               "FAKE_FAMS": str(self.fams), "FAKE_DFAM": str(self.dfam)}
        for k, v in env.items():
            monkeypatch.setenv(k, v)
        monkeypatch.setenv("PATH", str(self.bin) + (":/bin:/usr/bin" if sys_path else ""))
        if "BuildDatabase" in tools:
            _script(self.bin / "BuildDatabase", 'echo "BuildDatabase $*" >> "$FAKE_LOG"\n')
        if "RepeatModeler" in tools:
            _script(self.bin / "RepeatModeler",
                    'echo "RepeatModeler $*" >> "$FAKE_LOG"\n'
                    'db=$2\nmkdir -p RM_1\ncp "$FAKE_FAMS" "$db-families.fa"\n')
        if "minimap2" in tools:
            _script(self.bin / "minimap2", 'echo "minimap2 $*" >> "$FAKE_LOG"\ncat "$FAKE_PAF"\n')
        if "famdb.py" in tools:
            _script(self.bin / "famdb.py", 'echo "famdb.py $*" >> "$FAKE_LOG"\ncat "$FAKE_DFAM"\n')

    def calls(self):
        return self.log.read_text().splitlines()


def _genome(tmp_path):
    p = tmp_path / "g.fa"
    p.write_text(">chr1\n" + "ACGT" * 50 + "\n>Pt\n" + "TTTT" * 20 + "\n")
    return str(p)


def _fa(path, recs):
    path.write_text("".join(f">{h}\n{s}\n" for h, s in recs))


def _args(tmp_path, **kw):
    d = dict(genome=_genome(tmp_path), outdir=str(tmp_path / "out"), prefix="sp",
             threads=2, resume=False, organelle_contigs=None, genes=None,
             exclude_fasta=None, dfam_species=None, families=None)
    d.update(kw)
    return argparse.Namespace(**d)


def _names(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def _report(a):
    rows = [l.rstrip("\n").split("\t")
            for l in open(os.path.join(a.outdir, "sp.te-library.report.tsv"))
            if not l.startswith("#")]
    return {r[0]: r for r in rows[1:]}


def _paf(env, q, qlen, qs, qe, match, blk, t="exon1"):
    env.paf.write_text(f"{q}\t{qlen}\t{qs}\t{qe}\t+\t{t}\t1000\t0\t{qe - qs}\t{match}\t{blk}\t60\n")


def test_nuclear_fasta_excludes_organelle(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    a = _args(tmp_path, organelle_contigs="Pt")
    te_library.run_te_library(a)
    nuc = os.path.join(a.outdir, "sp.nuclear.fa")
    assert _names(nuc) == ["chr1"]
    assert any(c.startswith("BuildDatabase") and "sp.nuclear.fa" in c for c in env.calls())


def test_families_from_repeatmodeler_are_picked_up(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25), ("b#DNA/hAT", "TTGA" * 25)])
    a = _args(tmp_path)
    te_library.run_te_library(a)
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy", "b#DNA/hAT"]
    rm = [c for c in env.calls() if c.startswith("RepeatModeler")][0]
    assert "-LTRStruct" in rm and "-database sp" in rm and "-threads 2" in rm
    assert "-recoverDir" not in rm


def test_resume_passes_recoverdir(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    a = _args(tmp_path, resume=True)
    os.makedirs(os.path.join(a.outdir, "RM_123.1", "round-1"))
    os.makedirs(os.path.join(a.outdir, "RM_999.9"))  # RepeatClassifier temp dir: not a run
    te_library.run_te_library(a)
    rm = [c for c in env.calls() if c.startswith("RepeatModeler")][0]
    assert "-recoverDir" in rm and "RM_123.1" in rm and "RM_999.9" not in rm


def _gene_run(tmp_path, monkeypatch, match, blk, qs=0, qe=100, qlen=100):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("fam1#LTR/Gypsy", "ACGT" * (qlen // 4))])
    _paf(env, "fam1#LTR/Gypsy", qlen, qs, qe, match, blk)
    ex = tmp_path / "ex.fa"
    _fa(ex, [("exon1", "ACGT" * 50)])
    a = _args(tmp_path, exclude_fasta=[str(ex)])
    te_library.run_te_library(a)
    return a, env


def test_high_identity_high_coverage_dropped_as_gene(tmp_path, monkeypatch):
    a, env = _gene_run(tmp_path, monkeypatch, match=90, blk=100, qe=80)
    r = _report(a)["fam1#LTR/Gypsy"]
    assert r[2] == "dropped" and r[3] == "gene" and r[4] == "exon1"
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == []
    assert any("minimap2 -x asm20" in c for c in env.calls())


def test_low_identity_kept(tmp_path, monkeypatch):
    a, _ = _gene_run(tmp_path, monkeypatch, match=70, blk=100, qe=80)
    assert _report(a)["fam1#LTR/Gypsy"][2] == "kept"


def test_low_coverage_kept(tmp_path, monkeypatch):
    a, _ = _gene_run(tmp_path, monkeypatch, match=90, blk=100, qe=30)
    assert _report(a)["fam1#LTR/Gypsy"][2] == "kept"


def test_genes_option_extracts_exons(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("fam1#LTR/Gypsy", "ACGT" * 25)])
    _paf(env, "fam1#LTR/Gypsy", 100, 0, 100, 95, 100, t="exon1|chr1:10-60")
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n"
                   "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                   "chr1\tx\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n"
                   "chr1\tx\texon\t11\t60\t.\t+\t.\tID=e1;Parent=t1\n"
                   "chr1\tx\ttransposable_element_gene\t101\t200\t.\t+\t.\tID=te1\n"
                   "chr1\tx\texon\t101\t200\t.\t+\t.\tID=e2;Parent=te1\n")
    a = _args(tmp_path, genes=str(gff))
    te_library.run_te_library(a)
    assert _report(a)["fam1#LTR/Gypsy"][3] == "gene"
    ex = os.path.join(a.outdir, "sp.gene-exons.fa")
    assert _names(ex) == ["exon1|chr1:10-60"]


@pytest.mark.parametrize("cls", ["rRNA", "tRNA", "snRNA", "scRNA", "srpRNA"])
def test_rna_classes_dropped_without_filter(tmp_path, monkeypatch, cls):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    _fa(env.fams, [(f"x#{cls}", "ACGT" * 25), ("y#LTR/Gypsy", "ACGT" * 25)])
    a = _args(tmp_path)
    te_library.run_te_library(a)
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["y#LTR/Gypsy"]
    assert _report(a)[f"x#{cls}"][3] == cls


def test_satellite_kept(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    _fa(env.fams, [("s#Satellite", "ACGT" * 25)])
    a = _args(tmp_path)
    te_library.run_te_library(a)
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["s#Satellite"]


def test_no_filter_options_warns_and_keeps_non_rna(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25), ("b#rRNA", "ACGT" * 25),
                   ("c", "ACGT" * 25)])
    a = _args(tmp_path)
    te_library.run_te_library(a)
    assert "no --genes or --exclude-fasta" in capsys.readouterr().err
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy", "c"]
    assert "no sequence filtering" in open(
        os.path.join(a.outdir, "sp.te-library.report.tsv")).read()


def test_families_option_skips_repeatmodeler(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch, tools=("minimap2",))
    lib = tmp_path / "lib.fa"
    _fa(lib, [("a#LTR/Gypsy", "ACGT" * 25)])
    ex = tmp_path / "ex.fa"
    _fa(ex, [("exon1", "ACGT" * 50)])
    a = _args(tmp_path, families=str(lib), exclude_fasta=[str(ex)])
    te_library.run_te_library(a)
    assert not any(c.startswith(("BuildDatabase", "RepeatModeler")) for c in env.calls())
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy"]


def test_dfam_appended_with_prefix(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler", "famdb.py"))
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    _fa(env.dfam, [("AtREP1#DNA/Helitron", "GGCC" * 10), ("Sat1#Satellite", "AATT" * 10)])
    a = _args(tmp_path, dfam_species="Arabidopsis thaliana")
    te_library.run_te_library(a)
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == [
        "a#LTR/Gypsy", "dfam|AtREP1#DNA/Helitron", "dfam|Sat1#Satellite"]
    fd = [c for c in env.calls() if c.startswith("famdb.py")][0]
    assert "--curated" in fd and "-a" in fd and "-ad" in fd and "fasta_name" in fd
    assert "Arabidopsis thaliana" in fd


def test_missing_famdb_is_a_warning(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    a = _args(tmp_path, dfam_species="Arabidopsis thaliana")
    te_library.run_te_library(a)
    err = capsys.readouterr().err
    assert "famdb.py" in err and "WARN" in err
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy"]


def test_missing_minimap2_with_filter_is_error(tmp_path, monkeypatch, capsys):
    Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"), sys_path=False)
    ex = tmp_path / "ex.fa"
    _fa(ex, [("exon1", "ACGT" * 50)])
    a = _args(tmp_path, exclude_fasta=[str(ex)])
    with pytest.raises(SystemExit):
        te_library.run_te_library(a)
    err = capsys.readouterr().err
    assert "ERROR" in err and "minimap2" in err


def test_missing_repeatmodeler_is_error(tmp_path, monkeypatch, capsys):
    Env(tmp_path, monkeypatch, tools=(), sys_path=False)
    with pytest.raises(SystemExit):
        te_library.run_te_library(_args(tmp_path))
    assert "BuildDatabase" in capsys.readouterr().err


def test_existing_output_blocks_rerun(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    a = _args(tmp_path)
    te_library.run_te_library(a)
    n = len(env.calls())
    with pytest.raises(SystemExit):
        te_library.run_te_library(a)
    assert "exists" in capsys.readouterr().err
    assert len(env.calls()) == n


def test_setup_lists_te_library_tools_as_optional():
    import setup_module
    reqs = setup_module.MODULE_TOOLS["te-library"]
    assert {r["name"] for r in reqs} == {"BuildDatabase", "RepeatModeler", "minimap2", "famdb.py"}
    assert all("optional" in r["note"] for r in reqs)


def test_minimap2_argv_has_base_level_flags(tmp_path, monkeypatch):
    a, env = _gene_run(tmp_path, monkeypatch, match=90, blk=100, qe=80)
    mm = [c for c in env.calls() if c.startswith("minimap2")][0]
    assert "-x asm20" in mm and " -c " in mm and "-s 60" in mm and "-N 5" in mm


def test_repeatmodeler_output_kept_as_raw(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25), ("b#rRNA", "ACGT" * 25)])
    a = _args(tmp_path)
    te_library.run_te_library(a)
    raw = os.path.join(a.outdir, "sp-families.raw.fa")
    assert _names(raw) == ["a#LTR/Gypsy", "b#rRNA"]
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy"]
    assert "sp-families.raw.fa" in open(
        os.path.join(a.outdir, "sp.te-library.report.tsv")).read()


def test_finished_repeatmodeler_without_report_is_not_rerun(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    a = _args(tmp_path)
    os.makedirs(a.outdir)
    _fa(__import__("pathlib").Path(a.outdir) / "sp-families.fa", [("a#LTR/Gypsy", "ACGT" * 25)])
    te_library.run_te_library(a)
    assert env.calls() == []
    assert "skipped" in capsys.readouterr().err
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy"]


def test_missing_inputs_fail_before_repeatmodeler(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    for kw in (dict(genes=str(tmp_path / "nope.gff3")),
               dict(exclude_fasta=[str(tmp_path / "nope.fa")]),
               dict(families=str(tmp_path / "nope.fa"))):
        with pytest.raises(SystemExit):
            te_library.run_te_library(_args(tmp_path, **kw))
    assert env.calls() == []


def test_zero_exons_names_contig_mismatch(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n1\tx\texon\t11\t60\t.\t+\t.\tID=e1\n")
    with pytest.raises(SystemExit):
        te_library.run_te_library(_args(tmp_path, genes=str(gff)))
    err = capsys.readouterr().err
    assert "naming" in err and "1 features" in err


def test_empty_exclude_fasta_is_error(tmp_path, monkeypatch):
    Env(tmp_path, monkeypatch)
    ex = tmp_path / "e.fa"
    ex.write_text("")
    with pytest.raises(SystemExit):
        te_library.run_te_library(_args(tmp_path, exclude_fasta=[str(ex)]))


def test_genes_work_on_genome_without_fai_in_readonly_dir(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("fam1#LTR/Gypsy", "ACGT" * 25)])
    gff = tmp_path / "g.gff3"
    gff.write_text("chr1\tx\texon\t11\t60\t.\t+\t.\tID=e1\n")
    a = _args(tmp_path, genes=str(gff))
    te_library.run_te_library(a)
    assert not os.path.exists(a.genome + ".fai")


def test_dfam_rna_classes_dropped(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler", "famdb.py"))
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    _fa(env.dfam, [("U1#snRNA", "GGCC" * 10), ("X#DNA/hAT", "AATT" * 10)])
    a = _args(tmp_path, dfam_species="x")
    te_library.run_te_library(a)
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy", "dfam|X#DNA/hAT"]
    assert "1 Dfam families of RNA classes" in capsys.readouterr().err


def test_resume_without_run_dir_starts_fresh(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    te_library.run_te_library(_args(tmp_path, resume=True))
    assert "starting fresh" in capsys.readouterr().err
    assert "-recoverDir" not in [c for c in env.calls() if c.startswith("RepeatModeler")][0]


def test_resume_prefers_dir_mentioning_prefix(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch)
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    a = _args(tmp_path, resume=True)
    for n, txt in (("RM_1.a", "sp run"), ("RM_2.b", "other")):
        os.makedirs(os.path.join(a.outdir, n))
        open(os.path.join(a.outdir, n, "rmod.log"), "w").write(txt)
    te_library.run_te_library(a)
    rm = [c for c in env.calls() if c.startswith("RepeatModeler")][0]
    assert "RM_1.a" in rm


def test_duplicates_empty_and_zero_families(tmp_path, monkeypatch, capsys):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25), ("a#LTR/Gypsy extra", "TTTT" * 25),
                   ("e#DNA/hAT", "")])
    a = _args(tmp_path)
    te_library.run_te_library(a)
    err = capsys.readouterr().err
    assert "1 duplicate" in err and "zero-length" in err
    assert _names(os.path.join(a.outdir, "sp-families.fa")) == ["a#LTR/Gypsy"]
    assert _report(a)["e#DNA/hAT"][3] == "empty"


def test_zero_families_is_error(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    env.fams.write_text("")
    with pytest.raises(SystemExit):
        te_library.run_te_library(_args(tmp_path))


def test_crlf_fasta(tmp_path):
    p = tmp_path / "c.fa"
    p.write_bytes(b">a#LTR/Gypsy desc\r\nACGT\r\nAC\r\n")
    assert list(te_library.read_fasta(str(p))) == [("a#LTR/Gypsy desc", "ACGTAC")]


def test_partial_dfam_file_removed_on_failure(tmp_path, monkeypatch):
    env = Env(tmp_path, monkeypatch, tools=("BuildDatabase", "RepeatModeler"))
    _script(env.bin / "famdb.py", 'echo partial\nexit 1\n')
    _fa(env.fams, [("a#LTR/Gypsy", "ACGT" * 25)])
    a = _args(tmp_path, dfam_species="x")
    te_library.run_te_library(a)
    assert not os.path.exists(os.path.join(a.outdir, "sp.dfam.fa"))
