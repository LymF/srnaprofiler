#!/usr/bin/env python3
"""
setup_env.py
=============
Automatic dependency installer and environment checker for sRNA Profiler.

This script detects missing dependencies (Python, Perl, R, and external tools)
and installs them automatically when possible.
"""

import os
import sys
import subprocess
from utils import log

# -------------------------------------------------------
# Helper functions
# -------------------------------------------------------

def safe_run(cmd, desc=None, check=False):
    """Run a shell command with logging."""
    desc = desc or " ".join(cmd)
    try:
        subprocess.run(cmd, check=check)
        log(f"✅ {desc} completed successfully", "INFO")
    except subprocess.CalledProcessError as e:
        log(f"⚠️  Failed: {desc}", "WARN")
        log(f"{e}", "DEBUG")


# -------------------------------------------------------
# Python dependencies
# -------------------------------------------------------

def install_python_deps():
    required = ["pandas", "numpy", "scipy", "matplotlib", "seaborn", "biopython"]
    log("🐍 Checking Python packages...", "INFO")
    for pkg in required:
        try:
            __import__(pkg)
            log(f"✅ Python package found: {pkg}", "DEBUG")
        except ImportError:
            log(f"Installing missing Python package: {pkg}", "INFO")
            safe_run([sys.executable, "-m", "pip", "install", pkg])


# -------------------------------------------------------
# Perl dependencies
# -------------------------------------------------------

def install_perl_deps():
    required = [
        "Bio::SeqIO",
        "Bio::Seq::Quality",
        "Getopt::Long",
        "Statistics::Basic",
        "Statistics::RankCorrelation",
        "Number::Format",
    ]
    log("🧬 Checking Perl modules...", "INFO")

    perl_exec = shutil.which("perl")
    if not perl_exec:
        log("❌ Perl not found! Please install Perl or activate your Conda environment.", "ERROR")
        return

    for mod in required:
        test_cmd = ["perl", "-M" + mod, "-e", "1"]
        result = subprocess.run(test_cmd, capture_output=True)
        if result.returncode != 0:
            log(f"Installing missing Perl module: {mod}", "INFO")
            safe_run(["cpanm", "--notest", mod])
        else:
            log(f"✅ Perl module found: {mod}", "DEBUG")

    # Configure @INC for Conda Perl
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        site_path = Path(conda_prefix) / "lib/perl5/site_perl"
        site_path.mkdir(parents=True, exist_ok=True)
        sitecustom = site_path / "sitecustomize.pl"
        content = (
            f'use lib "{conda_prefix}/lib/perl5/site_perl";\n'
            f'use lib "{conda_prefix}/lib/perl5/site_perl/x86_64-linux-thread-multi";\n'
            f'use lib "{conda_prefix}/lib/site_perl";\n'
            f'use lib "{conda_prefix}/lib/site_perl/x86_64-linux-thread-multi";\n1;\n'
        )
        sitecustom.write_text(content)
        log("✅ sitecustomize.pl configured for Conda Perl", "INFO")


# -------------------------------------------------------
# R dependencies
# -------------------------------------------------------

def install_r_deps():
    required = ["ggplot2", "reshape2", "plyr", "scales"]
    r_exec = shutil.which("Rscript")
    if not r_exec:
        log("❌ Rscript not found! Install R or add it to PATH.", "ERROR")
        return

    log("📊 Checking R packages...", "INFO")
    for pkg in required:
        cmd = [
            "Rscript",
            "-e",
            f"if (!requireNamespace('{pkg}', quietly=TRUE)) install.packages('{pkg}', repos='https://cloud.r-project.org')",
        ]
        result = subprocess.run(cmd)
        if result.returncode == 0:
            log(f"✅ R package ready: {pkg}", "DEBUG")
        else:
            log(f"⚠️  Could not install or verify R package: {pkg}", "WARN")


# -------------------------------------------------------
# External tools
# -------------------------------------------------------

import shutil
from pathlib import Path

def check_tools():
    required = ["bowtie", "bowtie-build", "samtools", "fastqc", "cutadapt", "fastp"]
    log("🔧 Checking external executables...", "INFO")
    for tool in required:
        path = shutil.which(tool)
        if path:
            log(f"✅ {tool} found at {path}", "DEBUG")
        else:
            log(f"⚠️  {tool} not found! Install via: conda install -c bioconda {tool}", "WARN")


# -------------------------------------------------------
# Main setup workflow
# -------------------------------------------------------

def run_setup():
    log("=" * 60, "INFO")
    log("🚀 Starting full environment setup for sRNA Profiler", "INFO")
    log("=" * 60, "INFO")

    install_python_deps()
    install_perl_deps()
    install_r_deps()
    check_tools()

    log("=" * 60, "INFO")
    log("✅ All dependencies checked/installed successfully.", "INFO")
    log("You are ready to run sRNA Profiler!", "INFO")
    log("=" * 60, "INFO")
