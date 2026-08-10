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


_EBWT = ("1", "2", "3", "4", "rev.1", "rev.2")


def write_config(paths):
    """Record prepared paths, merging with whatever is already there.

    `--data` is run once per kind of reference, so overwriting would make the
    second run erase the first.
    """
    import yaml
    from utils import user_config

    current = user_config()
    current.update({k: str(v) for k, v in paths.items()})
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        yaml.safe_dump(current, fh, sort_keys=True)
    log(f"recorded in {path}: {', '.join(sorted(paths))}", "INFO")
    return path


def rfam_index_ready(prefix):
    """All six bowtie files, not just some: a partial index aligns nothing."""
    prefix = Path(prefix)
    return all((prefix.parent / f"{prefix.name}.{e}.ebwt").exists()
               for e in _EBWT)


def link_rfam_for_mirdeep_p2(rfam_prefix, pipeline):
    """Point miRDeep-P2's rfam_index at the index built for --ncrna-db.

    It ships `mature_index` but not `rfam_index`, which is a separate download,
    and without it the pipeline runs to completion and predicts nothing. The
    paths are derived from the script's own directory and have no flag, so a
    symlink is the way in that survives their updates -- and it avoids holding
    the same 2.9 GB twice.
    """
    rfam_prefix = Path(rfam_prefix)
    if not rfam_index_ready(rfam_prefix):
        raise RuntimeError(f"the Rfam index at {rfam_prefix} is incomplete; "
                           f"build it with extras/fetch_ncrna_db.sh")
    index_dir = Path(pipeline).parent / "scripts" / "index"
    if not index_dir.is_dir():
        raise RuntimeError(f"{index_dir} does not exist -- is {pipeline} really "
                           f"a miRDeep-P2 install?")
    n = 0
    for e in _EBWT:
        target = rfam_prefix.parent / f"{rfam_prefix.name}.{e}.ebwt"
        link = index_dir / f"rfam_index.{e}.ebwt"
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(target.resolve())
        n += 1
    log(f"linked {n} Rfam index files into {index_dir}", "INFO")
    return n


_P2_BUG = 'seq="$results_folder/${filename}.ori.fa"'
_P2_FIX = 'seq="$results_folder/${filename}/${filename}.ori.fa"'


def patch_mirdeep_p2(pipeline):
    """Fix the 1.1.4 path bug, keeping the original beside it.

    The script writes the collapsed reads to `$out/$name/$name.ori.fa` and reads
    them back from `$out/$name.ori.fa`, so every later step sees an empty file
    and the run ends with `end...` and zero predictions.

    Returns True if it patched, False if it was already patched. Refuses a
    script it does not recognise instead of rewriting it: a later version may
    fix this differently, and editing third-party code on a guess is worse than
    stopping.
    """
    pipeline = Path(pipeline)
    text = pipeline.read_text()
    if _P2_FIX in text:
        log(f"{pipeline.name} is already patched", "INFO")
        return False
    if _P2_BUG not in text:
        raise RuntimeError(
            f"could not find the known path bug in {pipeline}, and it is not "
            f"patched either. Refusing to edit a script I do not recognise; "
            f"check the version.")
    backup = Path(str(pipeline) + ".orig")
    if not backup.exists():
        backup.write_text(text)
    pipeline.write_text(text.replace(_P2_BUG, _P2_FIX))
    log(f"patched {pipeline.name} (original kept at {backup.name})", "INFO")
    return True


def run_data(args):
    """Prepare the reference data and record where it went."""
    recorded = {}

    if args.rfam_dir:
        prefix = Path(args.rfam_dir) / "ncrna"
        if rfam_index_ready(prefix):
            log(f"Rfam index already built at {prefix}", "INFO")
        else:
            script = Path(__file__).parent / "extras" / "fetch_ncrna_db.sh"
            log(f"building the Rfam index in {args.rfam_dir} "
                f"(downloads ~1.8 GB, then bowtie-build)", "INFO")
            subprocess.run(["bash", str(script), str(args.rfam_dir)], check=True)
            if not rfam_index_ready(prefix):
                raise RuntimeError(
                    f"{script.name} finished but {prefix}.*.ebwt is incomplete")
        recorded["ncrna_db"] = prefix

    if args.mirdp2:
        pipeline = Path(args.mirdp2)
        patch_mirdeep_p2(pipeline)
        rfam = recorded.get("ncrna_db")
        if rfam is None:
            from utils import user_config
            rfam = user_config().get("ncrna_db")
        if rfam:
            link_rfam_for_mirdeep_p2(rfam, pipeline)
        else:
            log("no Rfam index recorded yet, so miRDeep-P2's rfam_index was not "
                "linked -- run `setup --data --rfam-dir DIR` first, or it will "
                "predict nothing", "WARN")
        recorded["mirdp2_pipeline"] = pipeline

    if not recorded:
        log("nothing to do: pass --rfam-dir and/or --mirdp2", "WARN")
        return 1
    write_config(recorded)
    return 0


# Everything conda can install, in one transaction: a partial environment is
# how you get a tool that is present and unusable.
CONDA_PACKAGES = [
    "bowtie", "samtools", "fastp", "fastqc", "bedtools",
    "viennarna", "mirdeep2", "mirscore",
    "pysam", "pyyaml", "biopython", "pandas", "numpy", "scipy",
    "matplotlib", "pybedtools", "tqdm", "pytest",
]

MIRDP2_URL = ("https://cfhcable.dl.sourceforge.net/project/mirdp2/"
              "version%201.1.4/miRDP2-v1.1.4.tar.gz")


def tools_dir():
    """Where downloaded third-party scripts go: beside the code, not in $HOME.

    An installation that keeps its own tools can be moved, copied or deleted as
    one thing, and two checkouts do not fight over the same patched script.
    """
    return Path(__file__).resolve().parent / "tools"


def conda_install_command(packages):
    """Install into the environment that is active, never into base.

    Without CONDA_PREFIX conda would resolve to base, which is usually shared
    and not what anyone means by "install this for my project".
    """
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        raise RuntimeError(
            "no active conda environment (CONDA_PREFIX unset). Create and "
            "activate one first: conda env create -f environment.yml && "
            "conda activate srnaprofiler")
    return ["conda", "install", "-y", "-p", prefix,
            "-c", "bioconda", "-c", "conda-forge", *packages]


def install_mirdeep_p2(dest=None):
    """Fetch miRDeep-P2 into the tools directory and patch it.

    Downloaded from a named mirror: SourceForge's "latest/download" serves the
    *test data* rather than the program, and its default mirror often refuses
    connections outright.
    """
    dest = Path(dest or tools_dir())
    dest.mkdir(parents=True, exist_ok=True)
    pipeline = dest / "1.1.4" / "miRDP2-v1.1.4_pipeline.bash"
    if pipeline.exists():
        log(f"miRDeep-P2 already in {dest}", "INFO")
    else:
        tarball = dest / "miRDP2-v1.1.4.tar.gz"
        log(f"downloading miRDeep-P2 to {tarball}", "INFO")
        subprocess.run(["curl", "-sSL", "--fail", "-m", "1800",
                        "-o", str(tarball), MIRDP2_URL], check=True)
        subprocess.run(["tar", "-xzf", str(tarball), "-C", str(dest)], check=True)
        tarball.unlink()
        if not pipeline.exists():
            raise RuntimeError(f"the archive did not contain {pipeline.name}; "
                               f"the download may have served the test data")
    patch_mirdeep_p2(pipeline)
    return pipeline


def run_tools(args):
    """Install what conda can, fetch what it cannot, and verify the result."""
    log(f"installing into {os.environ.get('CONDA_PREFIX', '(no env)')}", "INFO")
    subprocess.run(conda_install_command(CONDA_PACKAGES), check=True)

    pipeline = install_mirdeep_p2()
    rfam = None
    from utils import user_config
    rfam = user_config().get("ncrna_db")
    if rfam and rfam_index_ready(rfam):
        link_rfam_for_mirdeep_p2(rfam, pipeline)
    else:
        log("miRDeep-P2's rfam_index is not linked yet: run "
            "`setup --data --rfam-dir DIR` or it will predict nothing", "WARN")
    write_config({"mirdp2_pipeline": pipeline})

    log("verifying...", "INFO")
    return run_check(args)
