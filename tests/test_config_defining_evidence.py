"""Every class must declare what defines it, and two classes must not be
defined by the same set -- that would make them indistinguishable."""
import pytest

from decider import load_config, RuleScorer

CONFIGS = ["rules.yaml", "rules.plant.yaml"]
KNOWN = set(RuleScorer.KNOWN_EVIDENCE)


def _flat_names(defining_evidence):
    """Flattens OR slots (nested lists) into a single set of evidence names."""
    names = set()
    for item in defining_evidence:
        names.update(item) if isinstance(item, list) else names.add(item)
    return names


def _key(spec):
    """Hashable key for the defining set: OR slots become inner frozensets.

    The size window is part of the key when `size` is one of the defining
    evidences, because e_size is scored against the class's own window: two
    classes both defined by [phase, size] are still telling different things
    apart if one asks for 21-22 nt and the other for 24.
    """
    names = frozenset(
        frozenset(item) if isinstance(item, list) else item
        for item in spec["defining_evidence"]
    )
    if "size" in _flat_names(spec["defining_evidence"]):
        window = spec.get("size") or {}
        return names, (window.get("lo"), window.get("hi"))
    return names, None


@pytest.mark.parametrize("path", CONFIGS)
def test_every_class_declares_defining_evidence(path):
    for name, spec in load_config(path)["classes"].items():
        assert spec.get("defining_evidence"), f"{path}:{name}"


@pytest.mark.parametrize("path", CONFIGS)
def test_defining_evidence_names_are_known(path):
    for name, spec in load_config(path)["classes"].items():
        assert _flat_names(spec["defining_evidence"]) <= KNOWN, f"{path}:{name}"


@pytest.mark.parametrize("path", CONFIGS)
def test_no_two_classes_share_the_same_defining_set(path):
    seen = {}
    for name, spec in load_config(path)["classes"].items():
        key = _key(spec)
        assert key not in seen, f"{path}: {name} and {seen[key]} are indistinguishable"
        seen[key] = name
