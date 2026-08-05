"""pad lives in the config: it is calibrated per organism.

Measured at 1 RPM: 500 keeps the piRNA-cluster recall of 200 in Drosophila with
twice the reads per locus, and 1000 erodes it. The CLI flag still wins when the
user passes one.
"""
import pytest

from decider import load_config
from classify import resolve_pad


@pytest.mark.parametrize("path,expected", [("rules.yaml", 500),
                                           ("rules.plant.yaml", 500)])
def test_each_config_declares_its_pad(path, expected):
    assert load_config(path)["clustering"]["pad"] == expected


def test_cli_overrides_the_config():
    cfg = {"clustering": {"pad": 500}}
    assert resolve_pad(cli_pad=120, config=cfg) == 120


def test_config_used_when_cli_absent():
    cfg = {"clustering": {"pad": 300}}
    assert resolve_pad(cli_pad=None, config=cfg) == 300


def test_falls_back_when_config_is_silent():
    assert resolve_pad(cli_pad=None, config={}) == 500


# The merge groups neighbouring sizes by the class windows of the config, so a
# 25-nt island joins a 26-nt one (both piRNA) but not a 21-nt one.
from classify import class_size_windows


def test_windows_come_from_the_class_definitions():
    cfg = {"classes": {"siRNA": {"size": {"lo": 20, "hi": 22}},
                       "piRNA": {"size": {"lo": 24, "hi": 31}}}}
    assert sorted(class_size_windows(cfg)) == [(20, 22), (24, 31)]


def test_classes_without_a_size_window_are_skipped():
    cfg = {"classes": {"miRNA": {}, "piRNA": {"size": {"lo": 24, "hi": 31}}}}
    assert class_size_windows(cfg) == [(24, 31)]


def test_empty_config_gives_no_windows():
    assert class_size_windows({}) == []


def test_real_configs_declare_windows():
    assert len(class_size_windows(load_config("rules.yaml"))) >= 2
    assert len(class_size_windows(load_config("rules.plant.yaml"))) >= 2
