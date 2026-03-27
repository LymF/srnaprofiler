#!/usr/bin/env python3
"""
Utility functions for sRNA Profiler
Shared functions across all modules
"""
import sys
import gzip
import shlex
import subprocess
import datetime
from pathlib import Path
from typing import Optional, Tuple, Union
from shutil import which as _which

# Configuration constants
DEFAULT_KMER_SIZE = 2000
DEFAULT_MIN_LEN = 15
DEFAULT_MAX_LEN = 36
DEFAULT_QTRIM = 20
DEFAULT_SI_THRESHOLD = 100
DEFAULT_PI_THRESHOLD = 300
DEFAULT_THREADS = 8
DEFAULT_MISMATCHES = 1
DEFAULT_MAX_ALIGNMENTS = 20

FASTQ_ADAPTERS = {
    'illumina_srna': 'TGGAATTCTCGGGTGCCAAGG',
    'truseq': 'AGATCGGAAGAGC',
    'nextera': 'CTGTCTCTTATACACATCT'
}

# Color codes for logging
LOG_COLORS = {
    "INFO": "\033[92m",    # Green
    "WARN": "\033[93m",    # Yellow
    "ERROR": "\033[91m",   # Red
    "DEBUG": "\033[94m",   # Blue
    "RESET": "\033[0m"
}


def log(msg: str, level: str = "INFO") -> None:
    """Enhanced logging with timestamps and colors."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    color = LOG_COLORS.get(level, "")
    reset = LOG_COLORS["RESET"]
    sys.stderr.write(f"{color}[{timestamp}] [{level}] {msg}{reset}\n")
    sys.stderr.flush()


def which(name: str) -> bool:
    """Check if a program is available in PATH."""
    return _which(name) is not None


def ensure_outdir(d: Path) -> Path:
    """Create output directory if it doesn't exist."""
    d.mkdir(parents=True, exist_ok=True)
    log(f"Output directory ready: {d}", "DEBUG")
    return d


def validate_file_exists(filepath: Path, description: str = "File") -> None:
    """Validate that a file exists."""
    if not filepath.exists():
        raise FileNotFoundError(f"{description} not found: {filepath}")
    log(f"{description} validated: {filepath}", "DEBUG")


def run_command(cmd: list[str], env: Optional[dict] = None) -> None:
    """Run a command with proper error handling."""
    log("Running: " + " ".join(shlex.quote(c) for c in cmd), "DEBUG")
    try:
        result = subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
        if result.stdout:
            log(f"Command output: {result.stdout.strip()}", "DEBUG")
    except subprocess.CalledProcessError as e:
        log(f"Command failed: {' '.join(cmd)}", "ERROR")
        if e.stderr:
            log(f"Error details: {e.stderr}", "ERROR")
        raise


def open_auto(filepath: Path, mode: str = "rt"):
    """Automatically handle gzipped and regular files."""
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, mode)
    return open(filepath, mode)


def detect_file_format(filepath: Path) -> str:
    """Detect if file is FASTQ or FASTA format."""
    with open_auto(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == '@':
                return "fastq"
            elif line[0] == '>':
                return "fasta"
            break
    return "unknown"


def detect_phred_encoding(filepath: Path, sample_size: int = 10000) -> int:
    """Detect Phred encoding (33 or 64) from FASTQ file."""
    if detect_file_format(filepath) != "fastq":
        return 33  # Default for non-FASTQ files
    
    min_qual = 126
    max_qual = 0
    
    with open_auto(filepath) as f:
        line_no = 0
        for line in f:
            line_no += 1
            if line_no % 4 == 0:  # Quality line
                for char in line.strip():
                    qual = ord(char)
                    min_qual = min(min_qual, qual)
                    max_qual = max(max_qual, qual)
            if line_no >= sample_size:
                break
    
    # Phred+33: ASCII 33-73 typically
    # Phred+64: ASCII 64-104 typically
    if min_qual < 59:
        return 33
    return 64


def convert_uracil_to_thymine_optimized(file_in: Path, file_out: Path) -> bool:
    """
    Optimized single-pass U->T conversion for RNA sequences.
    Detects format and converts in one pass.
    Returns True if conversion happened, False otherwise.
    """
    log(f"Checking for uracil bases in {file_in.name}...", "DEBUG")
    
    # Detect format
    file_format = detect_file_format(file_in)
    if file_format == "unknown":
        log(f"Unknown file format for {file_in.name}, skipping U->T conversion", "WARN")
        return False
    
    # Single-pass conversion with detection
    found_u = False
    temp_out = file_out.with_suffix('.tmp' + file_out.suffix)
    
    try:
        with open_auto(file_in) as fin, open_auto(temp_out, "wt") as fout:
            if file_format == "fastq":
                line_no = 0
                for line in fin:
                    line_no += 1
                    if line_no % 4 == 2:  # Sequence line
                        if not found_u and ('U' in line or 'u' in line):
                            found_u = True
                        fout.write(line.replace('U', 'T').replace('u', 't'))
                    else:
                        fout.write(line)
            
            elif file_format == "fasta":
                for line in fin:
                    if line and line[0] == '>':
                        fout.write(line)
                    else:
                        if not found_u and ('U' in line or 'u' in line):
                            found_u = True
                        fout.write(line.replace('U', 'T').replace('u', 't'))
        
        # Only keep converted file if U was found
        if found_u:
            temp_out.rename(file_out)
            log(f"Converted uracil (U) to thymine (T) in {file_in.name}", "INFO")
            return True
        else:
            temp_out.unlink()  # Remove temp file
            log(f"No uracil bases found in {file_in.name}", "DEBUG")
            return False
            
    except Exception as e:
        if temp_out.exists():
            temp_out.unlink()
        log(f"Error during U->T conversion: {e}", "ERROR")
        raise


def infer_format(path: str, user_choice: str = "auto") -> str:
    """Infer file format from extension or user choice."""
    if user_choice != "auto":
        return user_choice
    
    p = path.lower()
    if p.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
        return "fasta"
    elif p.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
        return "fastq"
    
    # Try to detect from content
    format_detected = detect_file_format(Path(path))
    if format_detected != "unknown":
        return format_detected
    
    # Default fallback
    return "fastq"


def get_read_stats(filepath: Path, max_reads: int = 100000) -> dict:
    """Get basic statistics from read file."""
    stats = {
        'total_reads': 0,
        'total_bases': 0,
        'min_length': float('inf'),
        'max_length': 0,
        'has_uracil': False,
        'format': detect_file_format(filepath)
    }
    
    with open_auto(filepath) as f:
        if stats['format'] == 'fastq':
            line_no = 0
            for line in f:
                line_no += 1
                if line_no % 4 == 2:  # Sequence line
                    seq = line.strip()
                    stats['total_reads'] += 1
                    stats['total_bases'] += len(seq)
                    stats['min_length'] = min(stats['min_length'], len(seq))
                    stats['max_length'] = max(stats['max_length'], len(seq))
                    if not stats['has_uracil'] and ('U' in seq or 'u' in seq):
                        stats['has_uracil'] = True
                if stats['total_reads'] >= max_reads:
                    break
        
        elif stats['format'] == 'fasta':
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        seq = ''.join(current_seq)
                        stats['total_reads'] += 1
                        stats['total_bases'] += len(seq)
                        stats['min_length'] = min(stats['min_length'], len(seq))
                        stats['max_length'] = max(stats['max_length'], len(seq))
                        if not stats['has_uracil'] and ('U' in seq or 'u' in seq):
                            stats['has_uracil'] = True
                        current_seq = []
                        if stats['total_reads'] >= max_reads:
                            break
                else:
                    current_seq.append(line)
            
            # Don't forget last sequence
            if current_seq and stats['total_reads'] < max_reads:
                seq = ''.join(current_seq)
                stats['total_reads'] += 1
                stats['total_bases'] += len(seq)
                stats['min_length'] = min(stats['min_length'], len(seq))
                stats['max_length'] = max(stats['max_length'], len(seq))
                if not stats['has_uracil'] and ('U' in seq or 'u' in seq):
                    stats['has_uracil'] = True
    
    if stats['total_reads'] == 0:
        stats['min_length'] = 0
    
    stats['avg_length'] = stats['total_bases'] / stats['total_reads'] if stats['total_reads'] > 0 else 0
    
    return stats


def check_disk_space(path: Path, required_mb: int = 1000) -> bool:
    """Check if there's enough disk space."""
    import shutil
    stat = shutil.disk_usage(path.parent if path.exists() else Path.cwd())
    available_mb = stat.free / (1024 * 1024)
    if available_mb < required_mb:
        log(f"Low disk space warning: {available_mb:.1f}MB available, {required_mb}MB recommended", "WARN")
        return False
    return True


def validate_arguments(args) -> None:
    """Common argument validation."""
    # Check output directory
    outdir = Path(args.outdir)
    if not check_disk_space(outdir, 1000):
        log("Consider freeing up disk space before proceeding", "WARN")
    
    # Validate input files
    if hasattr(args, 'reads1') and args.reads1:
        validate_file_exists(Path(args.reads1), "Input reads file")
    
    if hasattr(args, 'reads2') and args.reads2:
        validate_file_exists(Path(args.reads2), "Paired reads file")
    
    if hasattr(args, 'reads') and args.reads:
        validate_file_exists(Path(args.reads), "Input reads file")
    
    if hasattr(args, 'genome') and args.genome:
        validate_file_exists(Path(args.genome), "Genome file")
    
    if hasattr(args, 'sam') and args.sam:
        validate_file_exists(Path(args.sam), "SAM file")
    
    # Check for required tools
    if hasattr(args, 'subcmd'):
        if args.subcmd == 'qc':
            if not which("fastqc"):
                log("FastQC not found in PATH - QC plots will be skipped", "WARN")
            if not which("fastp") and not which("cutadapt"):
                raise RuntimeError("Neither fastp nor cutadapt found in PATH")
        
        elif args.subcmd == 'mapping':
            if not which("bowtie"):
                raise RuntimeError("bowtie not found in PATH")
            if not which("samtools"):
                raise RuntimeError("samtools not found in PATH")
        
        elif args.subcmd == 'clusters':
            # Check Python dependencies for clusters
            missing_py_deps = []
            try:
                import pandas
            except ImportError:
                missing_py_deps.append("pandas")
            try:
                import numpy
            except ImportError:
                missing_py_deps.append("numpy")
            try:
                import scipy
            except ImportError:
                missing_py_deps.append("scipy")
            
            if missing_py_deps:
                log(f"Missing Python packages for clusters: {', '.join(missing_py_deps)}", "ERROR")
                log("Install with: pip install " + " ".join(missing_py_deps), "INFO")
                raise RuntimeError("Missing Python dependencies for cluster analysis")
            
            # Check for R if plots are requested
            if not hasattr(args, 'no_plots') or not args.no_plots:
                if not which("Rscript"):
                    log("Rscript not found - plots will be skipped", "WARN")