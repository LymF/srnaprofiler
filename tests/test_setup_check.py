"""`setup --check` answers "will a run work here", not "is something installed".

Every module has a different set of tools, and the failures that cost this
project the most time were silent: a tool present but unusable. So the check is
per module, and for the pieces with a known silent failure it verifies the thing
that actually breaks -- ViennaRNA's Python binding rather than the package name,
the miRDeep-P2 patch rather than the file's existence.
"""
from pathlib import Path

import pytest

from setup_module import (check_python_package, check_tool, check_module,
                          MODULE_TOOLS, summarise)


def test_a_tool_on_the_path_is_found():
    assert check_tool("python3")["ok"] is True


def test_a_missing_tool_is_reported_not_raised():
    r = check_tool("definitely-not-a-real-tool-xyz")
    assert r["ok"] is False
    assert "not-a-real-tool" in r["name"]


def test_an_importable_package_is_found():
    assert check_python_package("json")["ok"] is True


def test_a_missing_package_is_reported_with_its_import_name():
    r = check_python_package("nonexistent_module_xyz")
    assert r["ok"] is False
    assert "nonexistent_module_xyz" in r["detail"]


def test_viennarna_is_checked_by_its_binding_not_its_package_name():
    """`conda install viennarna` succeeds on python 3.14 and `import RNA` fails.

    The binary and the binding are different interfaces: miRDeep calls RNAfold,
    miRScore imports RNA. Checking the package name would pass while miRScore
    could not start.
    """
    names = [c["name"] for c in MODULE_TOOLS["mirna"]]
    assert "RNAfold" in names
    assert any("RNA" == c.get("import_name") for c in MODULE_TOOLS["mirna"])


def test_classify_needs_almost_nothing_external():
    """A run starting from a BAM should not require the whole toolchain."""
    external = [c for c in MODULE_TOOLS["classify"] if c.get("kind") == "tool"]
    assert [c["name"] for c in external] == ["samtools"]


def test_every_module_is_covered():
    assert set(MODULE_TOOLS) >= {"qc", "mapping", "classify", "mirna", "plots"}


def test_a_module_check_reports_one_result_per_requirement():
    results = check_module("classify")
    assert len(results) == len(MODULE_TOOLS["classify"])
    assert all("ok" in r and "name" in r for r in results)


def test_the_summary_says_which_modules_are_usable():
    ok = summarise({"qc": [{"name": "fastp", "ok": True, "detail": ""}],
                    "mirna": [{"name": "RNAfold", "ok": False, "detail": "x"}]})
    assert ok["qc"] is True
    assert ok["mirna"] is False


# --- setup --data: prepare the references and record where they are ------

from setup_module import write_config, link_rfam_for_mirdeep_p2, rfam_index_ready


def test_recorded_paths_survive_a_round_trip(tmp_path, monkeypatch):
    cfg = tmp_path / "sub" / "paths.yaml"
    monkeypatch.setenv("SRNAPROFILER_CONFIG", str(cfg))
    write_config({"ncrna_db": "/data/ncrna"})
    from utils import user_config
    assert user_config()["ncrna_db"] == "/data/ncrna"


def test_recording_a_path_keeps_the_ones_already_there(tmp_path, monkeypatch):
    """Running `--data` twice for different things must not erase the first."""
    cfg = tmp_path / "paths.yaml"
    monkeypatch.setenv("SRNAPROFILER_CONFIG", str(cfg))
    write_config({"ncrna_db": "/data/ncrna"})
    write_config({"mirdp2_pipeline": "/opt/p.bash"})
    from utils import user_config
    got = user_config()
    assert got["ncrna_db"] == "/data/ncrna"
    assert got["mirdp2_pipeline"] == "/opt/p.bash"


def _fake_index(dirpath, prefix):
    dirpath.mkdir(parents=True, exist_ok=True)
    for e in ("1", "2", "3", "4", "rev.1", "rev.2"):
        (dirpath / f"{prefix}.{e}.ebwt").write_text("x")
    return dirpath / prefix


def test_an_index_is_ready_only_with_all_six_files(tmp_path):
    p = _fake_index(tmp_path / "idx", "ncrna")
    assert rfam_index_ready(p) is True
    (tmp_path / "idx" / "ncrna.rev.2.ebwt").unlink()
    assert rfam_index_ready(p) is False


def test_the_rfam_index_is_linked_into_mirdeep_p2(tmp_path):
    """One index serves both: miRDeep-P2's rfam_index is a separate download,
    and without it the pipeline runs to completion and predicts nothing."""
    rfam = _fake_index(tmp_path / "rfam", "ncrna")
    pipeline = tmp_path / "1.1.4" / "miRDP2.bash"
    (pipeline.parent / "scripts" / "index").mkdir(parents=True)
    pipeline.write_text("#!/bin/bash\n")

    n = link_rfam_for_mirdeep_p2(rfam, pipeline)

    assert n == 6
    linked = pipeline.parent / "scripts" / "index" / "rfam_index.1.ebwt"
    assert linked.is_symlink()
    assert linked.resolve() == (tmp_path / "rfam" / "ncrna.1.ebwt").resolve()


def test_linking_refuses_an_incomplete_index(tmp_path):
    rfam = _fake_index(tmp_path / "rfam", "ncrna")
    (tmp_path / "rfam" / "ncrna.3.ebwt").unlink()
    pipeline = tmp_path / "1.1.4" / "miRDP2.bash"
    (pipeline.parent / "scripts" / "index").mkdir(parents=True)
    pipeline.write_text("#!/bin/bash\n")
    with pytest.raises(RuntimeError, match="incomplete"):
        link_rfam_for_mirdeep_p2(rfam, pipeline)


from setup_module import patch_mirdeep_p2

BUGGY_LINE = '\tseq="$results_folder/${filename}.ori.fa"\n'


def _p2(tmp_path, line):
    p = tmp_path / "1.1.4" / "miRDP2.bash"
    (p.parent / "scripts" / "index").mkdir(parents=True)
    p.write_text("#!/bin/bash\necho hi\n" + line)
    return p


def test_the_path_bug_is_patched_and_the_original_kept(tmp_path):
    p = _p2(tmp_path, BUGGY_LINE)
    assert patch_mirdeep_p2(p) is True
    assert "${filename}/${filename}.ori.fa" in p.read_text()
    assert p.with_suffix(p.suffix + ".orig").exists()


def test_patching_twice_changes_nothing(tmp_path):
    """setup must be idempotent: running it again cannot corrupt an install."""
    p = _p2(tmp_path, BUGGY_LINE)
    patch_mirdeep_p2(p)
    first = p.read_text()
    assert patch_mirdeep_p2(p) is False
    assert p.read_text() == first


def test_an_unrecognised_script_is_refused_rather_than_edited(tmp_path):
    """Their next version may fix it differently; do not rewrite blindly."""
    p = _p2(tmp_path, 'seq="something else entirely"\n')
    with pytest.raises(RuntimeError, match="not find"):
        patch_mirdeep_p2(p)


# --- setup --tools: install into the active env, keep scripts in the repo ---

from setup_module import conda_install_command, tools_dir, CONDA_PACKAGES


def test_the_install_targets_the_active_environment(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/opt/conda/envs/mine")
    cmd = conda_install_command(["bowtie"])
    assert "-p" in cmd and "/opt/conda/envs/mine" in cmd
    assert "-y" in cmd


def test_installing_outside_a_conda_environment_is_refused(monkeypatch):
    """Without a target env, conda would write to base -- someone else's."""
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    with pytest.raises(RuntimeError, match="conda"):
        conda_install_command(["bowtie"])


def test_the_channels_are_the_ones_the_packages_live_in(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/opt/env")
    cmd = " ".join(conda_install_command(["bowtie"]))
    assert "bioconda" in cmd and "conda-forge" in cmd


def test_the_package_list_covers_every_module(monkeypatch):
    assert {"bowtie", "samtools", "fastp", "viennarna"} <= set(CONDA_PACKAGES)


def test_downloaded_scripts_live_beside_the_code():
    """Not in the user's home: an installation stays self-contained."""
    import setup_module
    assert tools_dir().name == "tools"
    assert tools_dir().parent == Path(setup_module.__file__).resolve().parent
