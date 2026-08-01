# tests/test_r_scripts_smoke.py
import shutil, subprocess
from pathlib import Path
import pytest

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript absent")
ROOT = Path(__file__).resolve().parent.parent


def test_plot_pingpong_runs(tmp_path):
    tsv = tmp_path / "pp.tsv"
    lines = ["overlap\tpairs"] + [f"{o}\t{5 if o != 10 else 90}" for o in range(1, 31)]
    tsv.write_text("\n".join(lines) + "\n")
    out = tmp_path / "x"
    r = subprocess.run(["Rscript", str(ROOT / "plotPingPong.R"), str(tsv), str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "x.pingpong.pdf").exists()
