"""The work tree stops depending on --engine."""
from mirna import work_paths


def test_the_work_tree_is_named_off_the_prefix(tmp_path):
    w = work_paths(tmp_path, "at")
    assert w["root"] == tmp_path / "at.mirna_work"
    assert w["engine"] == tmp_path / "at.mirna_work" / "engine"
    assert w["validation"] == tmp_path / "at.mirna_work" / "validation"


def test_the_engine_directory_does_not_carry_the_engine_name(tmp_path):
    """Whoever switches --engine finds the raw output in the same place.

    work_paths takes no engine argument, which is the point: there is no
    branch left that could name the directory after the tool that ran.
    """
    assert work_paths(tmp_path, "at")["engine"].name == "engine"
