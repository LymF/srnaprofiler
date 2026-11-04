import hashlib
import json
from pathlib import Path

import pytest

from benchmark import run as run_mod
from benchmark.manifest import Library


def make_lib(group="arthropod", organism="Bombyx mori"):
    return Library("SRR9", organism, "test", group, "ovary", 1, "http://x", ["t"])


@pytest.fixture
def fake(monkeypatch):
    calls = []
    state = {"bed": None, "clusters": True}

    def fake_run(cmd, **kwargs):
        if cmd[0] == "git":
            raise run_mod.subprocess.CalledProcessError(1, cmd)
        calls.append(cmd)
        sub = cmd[2]
        out = Path(cmd[cmd.index("-o") + 1])
        acc = cmd[cmd.index("-p") + 1]
        if sub == "mirna" and state["bed"] is not None:
            (out / f"{acc}.mirna_loci.bed").write_text(state["bed"])
        if sub == "classify" and state["clusters"]:
            (out / f"{acc}.clusters.tsv").write_text("x\n")

        class R:
            returncode = 0

        return R()

    monkeypatch.setattr(run_mod.subprocess, "run", fake_run)
    return calls, state


def go(tmp_path, lib=None, **kw):
    rules = tmp_path / "rules.yaml"
    rules.write_bytes(b"a: 1\n")
    return run_mod.run_library(lib or make_lib(), tmp_path / "r.fq.gz",
                               tmp_path / "g.fa", tmp_path / "w", rules, **kw)


def test_command_sequence(fake, tmp_path):
    calls, _ = fake
    res = go(tmp_path)
    assert [c[2] for c in calls] == ["qc", "mapping", "mirna", "classify"]
    w = str(tmp_path / "w")
    qc, mp, mi, cl = calls
    for c in calls:
        assert c[c.index("-p") + 1] == "SRR9" and c[c.index("-o") + 1] == w
        assert c[c.index("--threads") + 1] == "8"
    assert qc[qc.index("--reads1") + 1] == str(tmp_path / "r.fq.gz")
    assert qc[qc.index("-g") + 1] == str(tmp_path / "g.fa")
    assert mp[mp.index("--reads") + 1].endswith("SRR9.trimmed.R1.fastq.gz")
    assert mp[mp.index("-g") + 1].endswith("SRR9.genome.std.fa")
    assert "--ncrna-db" not in mp
    assert mi[mi.index("--index-prefix") + 1] == mp[mp.index("-o") + 1] + "/SRR9.bowtie_index"
    assert cl[cl.index("--sam") + 1].endswith("SRR9.bowtie.bam")
    assert cl[cl.index("--config") + 1] == str(tmp_path / "rules.yaml")
    assert res.clusters_tsv.name == "SRR9.clusters.tsv" and res.bam.name == "SRR9.bowtie.bam"


def test_mirna_runs_for_every_library_without_mature(fake, tmp_path):
    calls, _ = fake
    go(tmp_path, make_lib(organism="Bombyx mori"))
    mi = calls[2]
    assert "--mature" not in mi and "--hairpin" not in mi
    assert mi[mi.index("--species") + 1] == "Bombyx mori"


@pytest.mark.parametrize("group,engine", [("arthropod", "mirdeep2"), ("plant", "mirdeep-p2")])
def test_engine_by_taxon_group(fake, tmp_path, group, engine):
    calls, _ = fake
    go(tmp_path, make_lib(group))
    mi = calls[2]
    assert mi[mi.index("--engine") + 1] == engine


def test_empty_mirna_bed_does_not_fail(fake, tmp_path):
    calls, state = fake
    state["bed"] = ""
    res = go(tmp_path)
    assert "--mirna-loci" not in calls[3]
    assert json.loads(res.meta.read_text())["mirna_loci_used"] is False


def test_nonempty_mirna_bed_is_passed(fake, tmp_path):
    calls, state = fake
    state["bed"] = "chr1\t1\t60\tm1\n"
    res = go(tmp_path)
    cl = calls[3]
    assert cl[cl.index("--mirna-loci") + 1].endswith("SRR9.mirna_loci.bed")
    assert json.loads(res.meta.read_text())["mirna_loci_used"] is True


def test_meta_records_rules_sha256(fake, tmp_path):
    res = go(tmp_path, threads=2)
    meta = json.loads(res.meta.read_text())
    assert meta["rules_sha256"] == hashlib.sha256(b"a: 1\n").hexdigest()
    assert meta["accession"] == "SRR9" and meta["mirna_engine"] == "mirdeep2"
    assert meta["tool_git_sha"] == "unknown"
    assert meta["role"] == "test" and len(meta["commands"]) == 4
    assert meta["commands"][0][meta["commands"][0].index("--threads") + 1] == "2"


def test_missing_clusters_raises(fake, tmp_path):
    _, state = fake
    state["clusters"] = False
    with pytest.raises(FileNotFoundError):
        go(tmp_path)


def test_interpreter_from_env(fake, tmp_path, monkeypatch):
    calls, _ = fake
    monkeypatch.setenv("SRNAPROFILER_PYTHON", "/opt/env/bin/python")
    go(tmp_path)
    assert all(c[0] == "/opt/env/bin/python" for c in calls)
