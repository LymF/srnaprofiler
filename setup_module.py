#!/usr/bin/env python3
"""Prepare and verify an installation.

The question this answers is "will a run work here", which is not the same as
"is something installed". The failures that cost this project the most time were
tools that were present and unusable: `conda install viennarna` succeeds on
Python 3.14 while `import RNA` still fails, and miRDeep-P2 runs to completion
and predicts nothing when its path bug is unpatched or its Rfam index is absent.

So each requirement is checked through the interface that actually breaks.
"""
import importlib
import os
import shutil
import subprocess
from pathlib import Path

from utils import log, config_path


def _tool(name, note=""):
    return {"kind": "tool", "name": name, "note": note}


def _pkg(name, import_name=None, note=""):
    return {"kind": "package", "name": name,
            "import_name": import_name or name, "note": note}


# What each module needs, and nothing more: someone who only runs `classify`
# from an existing BAM should not be told to install fastp.
MODULE_TOOLS = {
    "qc": [
        _tool("fastp", "or cutadapt"),
        _tool("fastqc", "optional: the report is skipped without it"),
    ],
    "mapping": [
        _tool("bowtie"),
        _tool("bowtie-build"),
        _tool("samtools"),
    ],
    "classify": [
        _tool("samtools"),
        _pkg("pysam"),
        _pkg("pyyaml", "yaml"),
        _pkg("pandas"),
        _pkg("numpy"),
        _pkg("scipy"),
    ],
    "mirna": [
        _tool("bowtie"),
        _tool("RNAfold", "ViennaRNA binary, used by both engines"),
        _pkg("viennarna", "RNA",
             "ViennaRNA's Python binding, used by miRScore -- a different "
             "interface from the RNAfold binary, and missing on Python 3.14"),
        _tool("miRScore", "optional: validation is skipped without it"),
    ],
    "plots": [
        _tool("Rscript"),
        _pkg("matplotlib"),
    ],
}


def check_tool(name):
    """Is the executable on PATH, and which one would run."""
    found = shutil.which(name)
    return {"name": name, "ok": found is not None,
            "detail": found or f"{name} not found in PATH"}


def check_python_package(import_name):
    """Can it actually be imported -- not: is a package of that name installed."""
    try:
        mod = importlib.import_module(import_name)
    except Exception as exc:
        return {"name": import_name, "ok": False,
                "detail": f"cannot import {import_name}: {exc}"}
    return {"name": import_name, "ok": True,
            "detail": getattr(mod, "__file__", "") or "built-in"}


def check_module(module):
    """One result per requirement of `module`, in declaration order."""
    out = []
    for req in MODULE_TOOLS[module]:
        if req["kind"] == "tool":
            r = check_tool(req["name"])
        else:
            r = check_python_package(req["import_name"])
            r["name"] = req["name"]
        r["note"] = req.get("note", "")
        out.append(r)
    return out


def summarise(results_by_module):
    """Whether each module can run: every non-optional requirement satisfied."""
    usable = {}
    for module, results in results_by_module.items():
        usable[module] = all(r["ok"] or "optional" in r.get("note", "")
                             for r in results)
    return usable


def run_check(args):
    """Print what works, what does not, and which modules that leaves usable."""
    results = {m: check_module(m) for m in MODULE_TOOLS}
    for module, rows in results.items():
        print(f"\n{module}")
        for r in rows:
            mark = "ok  " if r["ok"] else ("warn" if "optional" in r["note"]
                                           else "MISS")
            detail = r["detail"] if r["ok"] else r["detail"]
            print(f"  [{mark}] {r['name']:<12} {detail[:60]}")
            if not r["ok"] and r["note"]:
                print(f"         {r['note']}")

    usable = summarise(results)
    print("\nmodules that can run here:")
    for module, ok in usable.items():
        print(f"  {module:<10} {'yes' if ok else 'no'}")

    cfg = config_path()
    print(f"\nrecorded paths: {cfg}" if cfg.exists()
          else f"\nno recorded paths yet ({cfg}); run `setup --data`")
    return 0 if all(usable.values()) else 1
