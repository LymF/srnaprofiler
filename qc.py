#!/usr/bin/env python3
"""
Quality Control module for sRNA Profiler
Handles read QC, trimming, and genome preparation
"""
import sys
import subprocess
from pathlib import Path
from typing import Optional

from utils import (
    log, which, ensure_outdir, run_command, validate_arguments,
    convert_uracil_to_thymine_optimized, get_read_stats,
    detect_phred_encoding, validate_file_exists,
    DEFAULT_MIN_LEN, DEFAULT_QTRIM,
    DEFAULT_THREADS, FASTQ_ADAPTERS
)


def fastqc(reads: str, outdir: Path, threads: int) -> None:
    """Run FastQC on reads."""
    if not which("fastqc"):
        log("FastQC not found in PATH; skipping FastQC analysis", "WARN")
        return
    
    try:
        cmd = ["fastqc", "-o", str(outdir), "--quiet", "--threads", str(threads), reads]
        run_command(cmd)
        log(f"FastQC completed for {Path(reads).name}", "INFO")
    except Exception as e:
        log(f"FastQC failed: {e}", "WARN")


def trim_with_fastp(
    in1: str,
    out1: Path,
    in2: Optional[str],
    out2: Optional[Path],
    min_len: int,
    qtrim: int,
    adapter: Optional[str],
    threads: int,
    phred: int,
    prefix: str
) -> dict:
    """Trim reads using fastp and return statistics."""
    if not which("fastp"):
        raise RuntimeError("fastp not found in PATH")
    
    html_report = out1.parent / f"{prefix}.fastp.html"
    json_report = out1.parent / f"{prefix}.fastp.json"
    
    cmd = [
        "fastp",
        "-i", in1, "-o", str(out1),
        "--length_required", str(min_len),
        "--qualified_quality_phred", str(qtrim),
        "-w", str(threads),
        "-h", str(html_report),
        "-j", str(json_report),
    ]
    
    if in2 and out2:
        cmd += ["-I", in2, "-O", str(out2)]
    
    if adapter:
        cmd += ["--adapter_sequence", adapter]
    else:
        # Auto-detect adapters
        cmd += ["--detect_adapter_for_pe" if in2 else "--adapter_sequence", FASTQ_ADAPTERS['illumina_srna']]
    
    # Handle Phred encoding
    if phred == 64:
        cmd += ["--phred64"]
    
    run_command(cmd)
    log(f"Trimming completed with fastp", "INFO")
    
    # Parse statistics from JSON if available
    stats = {}
    try:
        import json
        if json_report.exists():
            with open(json_report) as f:
                data = json.load(f)
                stats = {
                    'total_reads_before': data.get('summary', {}).get('before_filtering', {}).get('total_reads', 0),
                    'total_reads_after': data.get('summary', {}).get('after_filtering', {}).get('total_reads', 0),
                    'total_bases_after': data.get('summary', {}).get('after_filtering', {}).get('total_bases', 0)
                }
                log(f"Reads retained: {stats['total_reads_after']:,}/{stats['total_reads_before']:,} "
                    f"({100*stats['total_reads_after']/stats['total_reads_before']:.1f}%)", "INFO")
    except Exception as e:
        log(f"Could not parse fastp statistics: {e}", "DEBUG")
    
    return stats


def trim_with_cutadapt(
    in1: str,
    out1: Path,
    in2: Optional[str],
    out2: Optional[Path],
    min_len: int,
    qtrim: int,
    adapter: Optional[str],
    threads: int,
    phred: int,
    prefix: str
) -> dict:
    """Trim reads using cutadapt and return statistics."""
    if not which("cutadapt"):
        raise RuntimeError("cutadapt not found in PATH")
    
    report_file = out1.parent / f"{prefix}.cutadapt.report.txt"
    
    base_cmd = [
        "cutadapt",
        "-m", str(min_len),
        "-q", str(qtrim),
        "-j", str(threads),
        "--report", "minimal",
    ]
    
    # Add adapter sequences
    if adapter:
        if in2:
            base_cmd += ["-a", adapter, "-A", adapter]
        else:
            base_cmd += ["-a", adapter]
    else:
        default_adapter = FASTQ_ADAPTERS['illumina_srna']
        if in2:
            base_cmd += ["-a", default_adapter, "-A", default_adapter]
        else:
            base_cmd += ["-a", default_adapter]
    
    # Handle Phred encoding
    if phred == 64:
        base_cmd += ["--quality-base=64"]
    
    # Add input/output files
    cmd = base_cmd + ["-o", str(out1), in1]
    if in2 and out2:
        cmd += ["-p", str(out2), in2]
    
    # Run with report capture
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(report_file, 'w') as f:
            f.write(result.stdout)
        
        # Parse basic statistics from output
        stats = {}
        for line in result.stdout.split('\n'):
            if 'Total reads processed:' in line:
                stats['total_reads_before'] = int(line.split()[-1].replace(',', ''))
            elif 'Reads written (passing filters):' in line:
                stats['total_reads_after'] = int(line.split()[-2].replace(',', ''))
        
        if 'total_reads_before' in stats and 'total_reads_after' in stats:
            log(f"Reads retained: {stats['total_reads_after']:,}/{stats['total_reads_before']:,} "
                f"({100*stats['total_reads_after']/stats['total_reads_before']:.1f}%)", "INFO")
        
        log(f"Trimming completed with cutadapt", "INFO")
        
    except subprocess.CalledProcessError as e:
        log(f"Cutadapt failed: {e.stderr}", "ERROR")
        raise
    
    return stats


def standardize_fasta_headers(in_fa: Path, out_fa: Path) -> int:
    """
    Standardize FASTA headers to keep only the first word.
    Returns number of sequences processed.
    """
    import re
    log(f"Standardizing headers in {in_fa.name}", "INFO")
    
    seq_count = 0
    with open(in_fa, "rt") as i, open(out_fa, "wt") as o:
        for line in i:
            if line.startswith(">"):
                # Extract only the first word (sequence ID)
                match = re.match(r"^>(\S+)", line)
                if match:
                    o.write(f">{match.group(1)}\n")
                else:
                    o.write(line)
                seq_count += 1
            else:
                o.write(line)
    
    log(f"Standardized {seq_count} sequence headers", "INFO")
    return seq_count


def prepare_genome(genome_path, outdir, prefix):
    """Standardizes the FASTA (uppercase) without tiling. Returns the output path."""
    import os
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{prefix}.genome.fa")
    with open(genome_path) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                fout.write(line.strip().upper() + "\n")
    return out


def generate_qc_report(stats: dict, outdir: Path, prefix: str) -> None:
    """Generate a summary QC report."""
    report_file = outdir / f"{prefix}.qc_summary.txt"
    
    with open(report_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write(f"QC Report for {prefix}\n")
        f.write("="*60 + "\n\n")
        
        if 'read_stats_before' in stats:
            f.write("Input Reads Statistics:\n")
            f.write("-"*30 + "\n")
            rs = stats['read_stats_before']
            f.write(f"  Format:       {rs.get('format', 'unknown')}\n")
            f.write(f"  Total reads:  {rs.get('total_reads', 0):,}\n")
            f.write(f"  Total bases:  {rs.get('total_bases', 0):,}\n")
            f.write(f"  Read length:  {rs.get('min_length', 0)}-{rs.get('max_length', 0)} "
                   f"(avg: {rs.get('avg_length', 0):.1f})\n")
            f.write(f"  Has uracil:   {'Yes' if rs.get('has_uracil', False) else 'No'}\n\n")
        
        if 'trim_stats' in stats:
            f.write("Trimming Statistics:\n")
            f.write("-"*30 + "\n")
            ts = stats['trim_stats']
            if 'total_reads_before' in ts and 'total_reads_after' in ts:
                retained = 100 * ts['total_reads_after'] / ts['total_reads_before']
                f.write(f"  Reads before: {ts['total_reads_before']:,}\n")
                f.write(f"  Reads after:  {ts['total_reads_after']:,} ({retained:.1f}%)\n\n")
        
        if 'genome_stats' in stats:
            f.write("Genome Processing:\n")
            f.write("-"*30 + "\n")
            gs = stats['genome_stats']
            f.write(f"  Sequences:    {gs.get('sequences', 0)}\n\n")
        
        f.write("="*60 + "\n")
    
    log(f"QC report saved to {report_file}", "INFO")


def run_qc(args) -> None:
    """Main QC workflow."""
    # Validate arguments
    validate_arguments(args)
    
    # Setup output directory
    outdir = ensure_outdir(Path(args.outdir))
    
    # Collect statistics
    all_stats = {}
    
    # Process reads if provided
    if args.reads1:
        log("="*60, "INFO")
        log(f"Starting QC analysis for {args.prefix}", "INFO")
        log("="*60, "INFO")
        
        # Get initial read statistics
        log("Analyzing input reads...", "INFO")
        all_stats['read_stats_before'] = get_read_stats(Path(args.reads1))
        
        # Detect Phred encoding if auto
        if not hasattr(args, 'phred') or args.phred == 33:
            detected_phred = detect_phred_encoding(Path(args.reads1))
            if detected_phred != 33:
                log(f"Detected Phred+{detected_phred} encoding", "INFO")
                args.phred = detected_phred
        
        # Convert U->T if needed
        reads1_processed = Path(args.reads1)
        if all_stats['read_stats_before'].get('has_uracil', False):
            converted1 = outdir / f"{args.prefix}.reads1.DNA.fastq.gz"
            if convert_uracil_to_thymine_optimized(Path(args.reads1), converted1):
                reads1_processed = converted1
        
        reads2_processed = None
        if args.reads2:
            reads2_stats = get_read_stats(Path(args.reads2))
            if reads2_stats.get('has_uracil', False):
                converted2 = outdir / f"{args.prefix}.reads2.DNA.fastq.gz"
                if convert_uracil_to_thymine_optimized(Path(args.reads2), converted2):
                    reads2_processed = converted2
            else:
                reads2_processed = Path(args.reads2)
        
        # Run FastQC
        log("Running FastQC...", "INFO")
        fastqc(str(reads1_processed), outdir, args.threads)
        if reads2_processed:
            fastqc(str(reads2_processed), outdir, args.threads)
        
        # Trimming
        log("Trimming reads...", "INFO")
        trimmed1 = outdir / f"{args.prefix}.trimmed.R1.fastq.gz"
        trimmed2 = outdir / f"{args.prefix}.trimmed.R2.fastq.gz" if reads2_processed else None
        
        # Auto-detect tool if needed
        tool = args.tool
        if tool == "auto":
            if which("fastp"):
                tool = "fastp"
                log("Using fastp for trimming", "INFO")
            elif which("cutadapt"):
                tool = "cutadapt"
                log("Using cutadapt for trimming", "INFO")
            else:
                raise RuntimeError("Neither fastp nor cutadapt found in PATH")
        
        # Run trimming
        if tool == "fastp":
            trim_stats = trim_with_fastp(
                in1=str(reads1_processed), out1=trimmed1,
                in2=str(reads2_processed) if reads2_processed else None,
                out2=trimmed2,
                min_len=args.min_len, qtrim=args.qtrim,
                adapter=args.adapter, threads=args.threads,
                phred=args.phred, prefix=args.prefix
            )
        else:
            trim_stats = trim_with_cutadapt(
                in1=str(reads1_processed), out1=trimmed1,
                in2=str(reads2_processed) if reads2_processed else None,
                out2=trimmed2,
                min_len=args.min_len, qtrim=args.qtrim,
                adapter=args.adapter, threads=args.threads,
                phred=args.phred, prefix=args.prefix
            )
        
        all_stats['trim_stats'] = trim_stats
        
        # Run FastQC on trimmed reads
        log("Running FastQC on trimmed reads...", "INFO")
        fastqc(str(trimmed1), outdir, args.threads)
        if trimmed2:
            fastqc(str(trimmed2), outdir, args.threads)
        
        log(f"Trimmed reads saved:", "INFO")
        log(f"  R1: {trimmed1}", "INFO")
        if trimmed2:
            log(f"  R2: {trimmed2}", "INFO")
    
    # Process genome if provided
    if args.genome:
        log("-"*60, "INFO")
        log("Processing genome...", "INFO")
        
        validate_file_exists(Path(args.genome), "Genome FASTA")

        # Standardize headers
        std_genome = outdir / f"{args.prefix}.genome.std.fa"
        seq_count = standardize_fasta_headers(Path(args.genome), std_genome)

        # Prepare genome (uppercase, no tiling — map to full genome)
        prepared_genome = prepare_genome(str(std_genome), str(outdir), args.prefix)

        all_stats['genome_stats'] = {
            'sequences': seq_count,
        }

        log(f"Genome files created:", "INFO")
        log(f"  Standardized: {std_genome}", "INFO")
        log(f"  Prepared (no tiling): {prepared_genome}", "INFO")
    
    # Generate summary report
    generate_qc_report(all_stats, outdir, args.prefix)
    
    log("="*60, "INFO")
    log(f"QC analysis completed successfully!", "INFO")
    log(f"Results saved in: {outdir}", "INFO")
    log("="*60, "INFO")
