"""`setup --check` answers "will a run work here", not "is something installed".

Every module has a different set of tools, and the failures that cost this
project the most time were silent: a tool present but unusable. So the check is
per module, and for the pieces with a known silent failure it verifies the thing
that actually breaks -- ViennaRNA's Python binding rather than the package name,
the miRDeep-P2 patch rather than the file's existence.
"""
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
