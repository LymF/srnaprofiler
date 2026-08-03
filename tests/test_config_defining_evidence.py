"""Every class must declare what defines it, and two classes must not be
defined by the same set -- that would make them indistinguishable."""
import pytest

from decider import load_config

CONFIGS = ["rules.yaml", "rules.plant.yaml"]
KNOWN = {"1u", "5a", "fractop", "pingpong", "phase", "density"}


@pytest.mark.parametrize("path", CONFIGS)
def test_every_class_declares_defining_evidence(path):
    for name, spec in load_config(path)["classes"].items():
        assert spec.get("defining_evidence"), f"{path}:{name}"


@pytest.mark.parametrize("path", CONFIGS)
def test_defining_evidence_names_are_known(path):
    for name, spec in load_config(path)["classes"].items():
        assert set(spec["defining_evidence"]) <= KNOWN, f"{path}:{name}"


@pytest.mark.parametrize("path", CONFIGS)
def test_no_two_classes_share_the_same_defining_set(path):
    seen = {}
    for name, spec in load_config(path)["classes"].items():
        key = frozenset(spec["defining_evidence"])
        assert key not in seen, f"{path}: {name} and {seen[key]} are indistinguishable"
        seen[key] = name
