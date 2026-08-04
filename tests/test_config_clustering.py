"""pad lives in the config, because the kingdoms disagree on the optimum.

Measured: at pad 500 Drosophila loses 2 piRNA clusters while Arabidopsis is
still gaining. The CLI flag still wins when the user passes one.
"""
import pytest

from decider import load_config
from classify import resolve_pad


@pytest.mark.parametrize("path,expected", [("rules.yaml", 200),
                                           ("rules.plant.yaml", 500)])
def test_each_config_declares_its_pad(path, expected):
    assert load_config(path)["clustering"]["pad"] == expected


def test_cli_overrides_the_config():
    cfg = {"clustering": {"pad": 500}}
    assert resolve_pad(cli_pad=120, config=cfg) == 120


def test_config_used_when_cli_absent():
    cfg = {"clustering": {"pad": 500}}
    assert resolve_pad(cli_pad=None, config=cfg) == 500


def test_falls_back_when_config_is_silent():
    assert resolve_pad(cli_pad=None, config={}) == 200
