import shutil
from pathlib import Path
import pytest
import class_logo
from class_logo import plot_logo


def test_plot_logo_empty_returns_false(tmp_path):
    assert plot_logo([], "sense", str(tmp_path / "o.pdf")) is False


def test_plot_logo_missing_logomaker_returns_false(tmp_path, monkeypatch):
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **kw):
        if name == "logomaker":
            raise ImportError("simulated missing logomaker")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    assert plot_logo(["UACG", "UACG"], "sense", str(tmp_path / "o.pdf")) is False


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("logomaker") is None,
    reason="logomaker not installed")
def test_plot_logo_writes_file(tmp_path):
    out = tmp_path / "logo.pdf"
    ok = plot_logo(["UAGCU", "UAGCU", "UAGGU", "UACCU"], "sense", str(out))
    assert ok is True
    assert Path(out).exists() and Path(out).stat().st_size > 0
