import hashlib
import json
import os
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
        if sub == "qc":
            (out / f"{acc}.trimmed.R1.fastq.gz").write_text("t")
            (out / f"{acc}.genome.std.fa").write_text(">c\nA\n")
        if sub == "mapping":
            (out / f"{acc}.bowtie.bam").write_text("b")
            (out / f"{acc}.bowtie.bam.bai").write_text("i")
            if "-g" in cmd:  # mapping -g builds the index in the workdir
                for suffix in ("1", "2", "3", "4", "rev.1", "rev.2"):
                    (out / f"{acc}.bowtie_index.{suffix}.ebwt").write_text("x")
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
    if not rules.exists():
        rules.write_bytes(b"a: 1\n")
    for name, data in (("r.fq.gz", b"reads"), ("g.fa", b">c\nA\n")):
        if not (tmp_path / name).exists():
            (tmp_path / name).write_bytes(data)
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


def set_mtimes(directory, t):
    for f in Path(directory).rglob("*"):
        if f.is_file():
            os.utime(f, (t, t))


def age_inputs(tmp_path):
    for name in ("r.fq.gz", "g.fa", "rules.yaml"):
        os.utime(tmp_path / name, (500, 500))


def test_second_run_skips_every_step(fake, tmp_path):
    calls, state = fake
    go(tmp_path)
    assert len(calls) == 4
    set_mtimes(tmp_path / "w", 1000)
    age_inputs(tmp_path)
    calls.clear()
    res = go(tmp_path)
    assert calls == []
    meta = json.loads(res.meta.read_text())
    assert meta["skipped_steps"] == ["qc", "mapping", "mirna", "classify"]
    assert len(meta["commands"]) == 4  # the commands of the full pipeline stay recorded


def test_newer_reads_rerun_everything(fake, tmp_path):
    calls, _ = fake
    go(tmp_path)
    set_mtimes(tmp_path / "w", 1000)
    age_inputs(tmp_path)
    os.utime(tmp_path / "r.fq.gz", (2000, 2000))
    calls.clear()
    res = go(tmp_path)
    assert [c[2] for c in calls] == ["qc", "mapping", "mirna", "classify"]
    assert json.loads(res.meta.read_text())["skipped_steps"] == []


def test_only_missing_output_step_and_downstream_rerun(fake, tmp_path):
    calls, _ = fake
    go(tmp_path)
    set_mtimes(tmp_path / "w", 1000)
    age_inputs(tmp_path)
    (tmp_path / "w" / "SRR9.clusters.tsv").unlink()
    calls.clear()
    go(tmp_path)
    assert [c[2] for c in calls] == ["classify"]
    # a BAM without its index counts as incomplete: mapping reruns, and classify follows
    set_mtimes(tmp_path / "w", 1000)
    (tmp_path / "w" / "SRR9.bowtie.bam.bai").unlink()
    calls.clear()
    go(tmp_path)
    assert [c[2] for c in calls] == ["mapping", "classify"]


def test_mirna_with_no_loci_is_not_rerun(fake, tmp_path):
    calls, state = fake
    state["bed"] = None  # mirna produced no bed at all
    go(tmp_path)
    assert (tmp_path / "w" / "SRR9.mirna.done").exists()
    set_mtimes(tmp_path / "w", 1000)
    age_inputs(tmp_path)
    calls.clear()
    go(tmp_path)
    assert calls == []


def test_failed_mirna_leaves_no_done_marker(fake, tmp_path, monkeypatch):
    orig = run_mod.subprocess.run

    def boom(cmd, **kw):
        if len(cmd) > 2 and cmd[2] == "mirna":
            raise run_mod.subprocess.CalledProcessError(1, cmd)
        return orig(cmd, **kw)

    monkeypatch.setattr(run_mod.subprocess, "run", boom)
    with pytest.raises(run_mod.subprocess.CalledProcessError):
        go(tmp_path)
    assert not (tmp_path / "w" / "SRR9.mirna.done").exists()


def test_shared_index_is_built_once_per_species(fake, tmp_path):
    calls, _ = fake
    shared = tmp_path / "idx" / "Bombyx_mori"
    prefix = shared / "Bombyx_mori"
    lib1 = make_lib()
    go(tmp_path, lib1, index_prefix=prefix)
    mp, mi = calls[1], calls[2]
    assert "-g" in mp and "--index-prefix" not in mp
    assert mi[mi.index("--index-prefix") + 1] == str(prefix)
    assert (shared / "Bombyx_mori.1.ebwt").exists() and (shared / "Bombyx_mori.rev.2.ebwt").exists()

    calls.clear()
    lib2 = Library("SRR10", "Bombyx mori", "test", "arthropod", "ovary", 2, "http://x", ["t"])
    run_mod.run_library(lib2, tmp_path / "r.fq.gz", tmp_path / "g.fa", tmp_path / "w2",
                        tmp_path / "rules.yaml", index_prefix=prefix)
    mp, mi = calls[1], calls[2]
    assert mp[mp.index("--index-prefix") + 1] == str(prefix) and "-g" not in mp
    assert mi[mi.index("--index-prefix") + 1] == str(prefix)


def test_pipeline_sha256_covers_root_py_and_rules_only(tmp_path):
    (tmp_path / "a.py").write_text("x = 1\n")
    (tmp_path / "rules.yaml").write_text("r: 1\n")
    (tmp_path / "rules.plant.yaml").write_text("p: 1\n")
    (tmp_path / "benchmark").mkdir()
    (tmp_path / "benchmark" / "b.py").write_text("y = 1\n")
    (tmp_path / "other.yaml").write_text("z: 1\n")
    base = run_mod.pipeline_sha256(tmp_path)
    assert len(base) == 64 and base == run_mod.pipeline_sha256(tmp_path)
    (tmp_path / "benchmark" / "b.py").write_text("y = 2\n")
    (tmp_path / "other.yaml").write_text("z: 2\n")
    assert run_mod.pipeline_sha256(tmp_path) == base  # benchmark/ and other files do not count
    for name in ("a.py", "rules.yaml", "rules.plant.yaml"):
        f = tmp_path / name
        old = f.read_text()
        f.write_text(old + "# edit\n")
        assert run_mod.pipeline_sha256(tmp_path) != base, name
        f.write_text(old)
    (tmp_path / "new.py").write_text("")
    assert run_mod.pipeline_sha256(tmp_path) != base


def test_meta_records_pipeline_sha_and_git_state(fake, tmp_path):
    res = go(tmp_path)
    meta = json.loads(res.meta.read_text())
    assert meta["pipeline_sha256"] == run_mod.pipeline_sha256(run_mod.REPO)
    assert meta["tool_git_sha"] == "unknown" and meta["git_dirty"] is None  # no git available


def test_git_dirty_reflects_porcelain_output(monkeypatch):
    class R:
        def __init__(self, out):
            self.stdout = out

    monkeypatch.setattr(run_mod.subprocess, "run", lambda cmd, **kw: R(" M srnaprofiler.py\n"))
    assert run_mod._git_dirty() is True
    monkeypatch.setattr(run_mod.subprocess, "run", lambda cmd, **kw: R(""))
    assert run_mod._git_dirty() is False


def test_meta_records_md5_of_reads_and_genome(fake, tmp_path):
    res = go(tmp_path)
    meta = json.loads(res.meta.read_text())
    assert meta["reads_md5"] == hashlib.md5(b"reads").hexdigest()
    assert meta["genome_md5"] == hashlib.md5(b">c\nA\n").hexdigest()


def _prepare_fresh_run(tmp_path):
    go(tmp_path)
    set_mtimes(tmp_path / "w", 1000)
    age_inputs(tmp_path)


def test_resume_skips_steps_under_same_pipeline_hash(fake, tmp_path):
    calls, _ = fake
    _prepare_fresh_run(tmp_path)
    calls.clear()
    res = go(tmp_path)
    assert calls == []
    meta = json.loads(res.meta.read_text())
    assert meta["resumed_from_pipeline_sha256"] == run_mod.pipeline_sha256()


def test_resume_reruns_everything_under_different_pipeline_hash(fake, tmp_path):
    calls, _ = fake
    _prepare_fresh_run(tmp_path)
    meta_path = tmp_path / "w" / "benchmark_run.json"
    stored = json.loads(meta_path.read_text())
    stored["pipeline_sha256"] = "0" * 64
    meta_path.write_text(json.dumps(stored))
    set_mtimes(tmp_path / "w", 1000)
    calls.clear()
    res = go(tmp_path)
    assert [c[2] for c in calls] == ["qc", "mapping", "mirna", "classify"]
    meta = json.loads(res.meta.read_text())
    assert meta["skipped_steps"] == []
    assert meta["resumed_from_pipeline_sha256"] == "0" * 64
    assert meta["pipeline_sha256"] == run_mod.pipeline_sha256()


def test_no_previous_meta_runs_everything_even_with_fresh_outputs(fake, tmp_path):
    calls, _ = fake
    _prepare_fresh_run(tmp_path)
    (tmp_path / "w" / "benchmark_run.json").unlink()
    set_mtimes(tmp_path / "w", 1000)
    calls.clear()
    res = go(tmp_path)
    assert [c[2] for c in calls] == ["qc", "mapping", "mirna", "classify"]
    assert json.loads(res.meta.read_text())["resumed_from_pipeline_sha256"] is None


def test_index_ready_needs_the_six_exact_nonempty_files(tmp_path):
    prefix = tmp_path / "idx"
    for suffix in ("1", "2", "3", "4", "rev.1"):
        (tmp_path / f"idx.{suffix}.ebwt").write_text("x")
    (tmp_path / "idx.rev.2.ebwt.part").write_text("x")
    assert not run_mod._index_ready(prefix)
    (tmp_path / "idx.rev.2.ebwt").write_text("")
    assert not run_mod._index_ready(prefix)  # empty file does not count
    (tmp_path / "idx.rev.2.ebwt").write_text("x")
    assert run_mod._index_ready(prefix)


def _classify(calls):
    return next(c for c in calls if c[2] == "classify")


def test_default_classify_argv_has_no_classifier(fake, tmp_path):
    go(tmp_path)
    argv = _classify(fake[0])
    assert "--classifier" not in argv and "--ml-model" not in argv


def test_ml_classify_argv_has_classifier_and_model(fake, tmp_path):
    model = tmp_path / "m.joblib"
    model.write_bytes(b"model")
    go(tmp_path, classifier="ml", ml_model=model)
    argv = _classify(fake[0])
    assert argv[argv.index("--classifier") + 1] == "ml"
    assert argv[argv.index("--ml-model") + 1] == str(model)
    meta = json.loads((tmp_path / "w" / "benchmark_run.json").read_text())
    assert meta["classifier"] == "ml"
    assert meta["ml_model_sha256"] == hashlib.sha256(b"model").hexdigest()


def test_ml_model_is_a_classify_freshness_input(fake, tmp_path, monkeypatch):
    model = tmp_path / "m.joblib"
    model.write_bytes(b"model")
    seen = []
    orig = run_mod._fresh
    monkeypatch.setattr(run_mod, "_fresh", lambda o, i: seen.append(list(i)) or orig(o, i))
    go(tmp_path, classifier="ml", ml_model=model)
    go(tmp_path, classifier="ml", ml_model=model)
    assert any(model in inputs for inputs in seen)


def test_ml_model_defaults_to_rules_yaml_path(fake, tmp_path, monkeypatch):
    monkeypatch.setattr(run_mod, "REPO", tmp_path)
    (tmp_path / "models").mkdir()
    (tmp_path / "models" / "x.joblib").write_bytes(b"m")
    (tmp_path / "rules.yaml").write_text("ml:\n  model: models/x.joblib\n")
    go(tmp_path, classifier="ml")
    assert "--ml-model" not in _classify(fake[0])
    meta = json.loads((tmp_path / "w" / "benchmark_run.json").read_text())
    assert meta["ml_model_path"] == str(tmp_path / "models" / "x.joblib")


def test_switching_classifier_reruns_classify_once(fake, tmp_path):
    calls, _ = fake
    model = tmp_path / "m.joblib"
    model.write_bytes(b"model")
    _prepare_fresh_run(tmp_path)  # fresh rules-mode outputs and meta
    calls.clear()
    go(tmp_path, classifier="ml", ml_model=model)
    assert [c[2] for c in calls] == ["classify"]
    set_mtimes(tmp_path / "w", 2000)
    age_inputs(tmp_path)
    os.utime(model, (500, 500))
    calls.clear()
    go(tmp_path, classifier="ml", ml_model=model)
    assert calls == []
    # a different model forces classify again
    model.write_bytes(b"other model")
    calls.clear()
    go(tmp_path, classifier="ml", ml_model=model)
    assert [c[2] for c in calls] == ["classify"]
