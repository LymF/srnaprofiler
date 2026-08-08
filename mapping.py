#!/usr/bin/env python3
"""
Mapping module for sRNA Profiler
Handles read alignment using Bowtie1
"""
import os
import sys
import gzip
import subprocess
from pathlib import Path
from typing import Optional, Tuple

from utils import (
    log, which, ensure_outdir, run_command, validate_arguments,
    convert_uracil_to_thymine_optimized, infer_format,
    validate_file_exists, get_read_stats, check_disk_space,
    DEFAULT_THREADS, DEFAULT_MISMATCHES, DEFAULT_MAX_ALIGNMENTS
)


def count_reads(path: str) -> int:
    """Number of reads in a FASTA/FASTQ, gzipped or not."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = fh.readline()
        if not first:
            return 0
        if first.startswith("@"):
            return 1 + sum(1 for i, _ in enumerate(fh) if i % 4 == 3)
        return 1 + sum(1 for line in fh if line.startswith(">"))


def build_ncrna_filter_cmd(db_prefix: str, reads_in: str, reads_out: str,
                           threads: int = 1) -> list:
    """bowtie command that keeps only the reads NOT matching the ncRNA database.

    tRNA, rRNA, snoRNA and snRNA fragments are abundant in a small RNA library
    and are not siRNA or piRNA. Removing them here, before clustering, means no
    locus is built from them -- filtering afterwards leaves the locus in place
    with its abundance already wrong.
    """
    return ["bowtie", "-v", "1", "-f", "-p", str(threads),
            db_prefix, reads_in, "--un", reads_out, "-S", "/dev/null"]


def bowtie_build(genome_fa: str, out_prefix: str, threads: int = 1) -> None:
    """Build Bowtie1 index from genome FASTA."""
    if not which("bowtie-build"):
        raise RuntimeError("bowtie-build not found in PATH")
    
    log(f"Building Bowtie index for {Path(genome_fa).name}", "INFO")
    
    cmd = ["bowtie-build", "--threads", str(threads), genome_fa, out_prefix]
    run_command(cmd)
    
    log(f"Bowtie index built: {out_prefix}", "INFO")


def check_bowtie_index(index_prefix: str) -> bool:
    """Check if Bowtie index files exist."""
    index_path = Path(index_prefix)
    required_extensions = ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']
    
    for ext in required_extensions:
        if not (Path(str(index_path) + ext).exists() or 
                Path(str(index_path) + ext.replace('.ebwt', '.ebwtl')).exists()):
            return False
    return True


def run_bowtie_alignment(
    index_prefix: str,
    reads1: str,
    reads2: Optional[str],
    sam_out: Path,
    threads: int,
    mismatches: int,
    max_alignments: int,
    fmt: str,
    keep_unmapped: bool,
    prefix: str
) -> dict:
    """Run Bowtie alignment and return statistics."""
    if not which("bowtie"):
        raise RuntimeError("bowtie not found in PATH")
    if not which("samtools"):
        raise RuntimeError("samtools not found in PATH")
    
    # Determine format flag
    fmt_flag = "-f" if fmt == "fasta" else "-q"
    
    # Build bowtie command
    bowtie_cmd = [
        "bowtie",
        fmt_flag,
        "-S",  # SAM output
        "-p", str(threads),
        "-k", str(max_alignments),  # Report up to k alignments
        "-v", str(mismatches),  # Allow v mismatches
        "--best",  # Report best alignments
        "--strata",  # Report only from best stratum
        "-x", index_prefix
    ]
    
    # Add read files
    if reads2:
        bowtie_cmd += ["-1", reads1, "-2", reads2]
    else:
        bowtie_cmd += [reads1]
    
    # Build samtools command for filtering
    samtools_cmd = ["samtools", "view", "-h", "-S"]
    if not keep_unmapped:
        samtools_cmd += ["-F", "4"]  # Exclude unmapped reads
    
    # Add threading for samtools if available
    samtools_cmd += ["-@", str(max(1, threads // 2))]
    
    # Run pipeline
    log(f"Starting alignment with Bowtie1...", "INFO")
    log(f"  Index: {index_prefix}", "DEBUG")
    log(f"  Reads: {reads1}" + (f" and {reads2}" if reads2 else ""), "DEBUG")
    log(f"  Output: {sam_out} (keep unmapped: {keep_unmapped})", "DEBUG")
    
    stats = {}
    stderr_file = sam_out.parent / f"{prefix}.bowtie.stderr.txt"
    
    try:
        # Run bowtie and capture stderr for statistics
        with open(sam_out, "wb") as out_file, open(stderr_file, "w") as err_file:
            # Run bowtie
            p1 = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE, stderr=err_file)
            # Pipe through samtools
            p2 = subprocess.Popen(samtools_cmd, stdin=p1.stdout, stdout=out_file, stderr=subprocess.PIPE)
            p1.stdout.close()
            
            # Wait for completion
            stderr2 = p2.communicate()[1]
            rc2 = p2.returncode
            rc1 = p1.wait()
            
            if rc1 != 0 or rc2 != 0:
                log(f"Bowtie return code: {rc1}, Samtools return code: {rc2}", "ERROR")
                if stderr2:
                    log(f"Samtools error: {stderr2.decode()}", "ERROR")
                raise RuntimeError("Bowtie alignment pipeline failed")
        
        # Parse bowtie statistics from stderr
        if stderr_file.exists():
            with open(stderr_file) as f:
                content = f.read()
                for line in content.split('\n'):
                    if 'reads processed:' in line:
                        stats['total_reads'] = int(line.split(':')[1].strip())
                    elif 'reads with at least one' in line and 'alignment:' in line:
                        # Format: "# reads with at least one alignment: 12345 (67.89%)"
                        parts = line.split(':')[1].strip().split()
                        stats['aligned_reads'] = int(parts[0])
                        if len(parts) > 1 and parts[1].startswith('('):
                            stats['alignment_rate'] = parts[1].strip('()%')
                    elif 'reads that failed to align:' in line:
                        parts = line.split(':')[1].strip().split()
                        stats['unaligned_reads'] = int(parts[0])
        
        # Log statistics
        if stats:
            log(f"Alignment completed:", "INFO")
            if 'total_reads' in stats:
                log(f"  Total reads: {stats['total_reads']:,}", "INFO")
            if 'aligned_reads' in stats:
                log(f"  Aligned reads: {stats['aligned_reads']:,} ({stats.get('alignment_rate', 'N/A')})", "INFO")
            if 'unaligned_reads' in stats:
                log(f"  Unaligned reads: {stats['unaligned_reads']:,}", "INFO")
        
        log(f"SAM file created: {sam_out}", "INFO")
        
    except Exception as e:
        log(f"Alignment failed: {e}", "ERROR")
        if stderr_file.exists():
            with open(stderr_file) as f:
                log(f"Bowtie stderr: {f.read()}", "ERROR")
        raise
    
    return stats


def available_memory_mb() -> int:
    """Physical memory currently free, in MB (0 if it cannot be probed)."""
    try:
        return (os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")) // (1024 ** 2)
    except (ValueError, OSError, AttributeError):
        return 0


def sort_memory_per_thread_mb(threads: int, budget_mb: int = None) -> int:
    """Per-thread share for `samtools sort -m`, so the total stays bounded.

    samtools defaults to 768 MB PER THREAD, which is not a budget at all: it
    scales with --threads and asks for ~9 GB at 12 threads, killing the sort on
    a smaller machine. We spend at most half of the free memory (capped at
    4 GB, since more buys little) and split that across the threads.
    """
    if budget_mb is None:
        free = available_memory_mb()
        budget_mb = min(4096, max(768, free // 2)) if free else 768
    return max(64, budget_mb // max(1, threads))


def convert_sam_to_bam(sam_file: Path, bam_file: Path, threads: int, sort: bool = True) -> None:
    """Convert SAM to BAM format with optional sorting."""
    if not which("samtools"):
        raise RuntimeError("samtools not found in PATH")

    log(f"Converting SAM to BAM format...", "INFO")

    if sort:
        per_thread = sort_memory_per_thread_mb(threads)
        log(f"Sort memory: {per_thread} MB x {threads} thread(s)", "DEBUG")
        # Convert and sort in one step. -T keeps the temporary files beside the
        # output under a known prefix instead of a guessed one, so two runs
        # writing to the same directory cannot collide.
        cmd = [
            "samtools", "sort",
            "-@", str(threads),
            "-m", f"{per_thread}M",
            "-T", str(bam_file) + ".tmp",
            "-o", str(bam_file),
            str(sam_file)
        ]
    else:
        # Just convert
        cmd = [
            "samtools", "view",
            "-bS",
            "-@", str(threads),
            "-o", str(bam_file),
            str(sam_file)
        ]
    
    run_command(cmd)
    
    # Create index
    if sort:
        log(f"Creating BAM index...", "INFO")
        run_command(["samtools", "index", str(bam_file)])
    
    log(f"BAM file created: {bam_file}", "INFO")


def get_alignment_stats(sam_file: Path, threads: int = 1) -> dict:
    """Extract detailed alignment statistics from SAM file."""
    if not which("samtools"):
        log("samtools not found, skipping detailed statistics", "WARN")
        return {}
    
    stats = {}
    
    try:
        # flagstat with -@ threads. No short timeout: on a large (multi-GB) SAM
        # a 120 s cap silently killed this and returned empty stats. flagstat is
        # a single linear pass, so let it run to completion.
        log(f"Running samtools flagstat with {threads} thread(s)...", "INFO")
        result = subprocess.run(
            ["samtools", "flagstat", "-@", str(threads), str(sam_file)],
            capture_output=True, text=True, check=True
        )

        for line in result.stdout.split('\n'):
            if 'in total' in line:
                stats['total_alignments'] = int(line.split()[0])
            elif 'mapped (' in line and 'primary' not in line:
                stats['mapped_alignments'] = int(line.split()[0])
            elif 'properly paired' in line:
                stats['proper_pairs'] = int(line.split()[0])

    except Exception as e:
        log(f"Could not get detailed statistics: {e}", "WARN")

    return stats


def generate_mapping_report(stats: dict, outdir: Path, prefix: str) -> None:
    """Generate mapping statistics report."""
    report_file = outdir / f"{prefix}.mapping_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write(f"Mapping Report for {prefix}\n")
        f.write("="*60 + "\n\n")
        
        if 'read_stats' in stats:
            f.write("Input Read Statistics:\n")
            f.write("-"*30 + "\n")
            rs = stats['read_stats']
            f.write(f"  Format:       {rs.get('format', 'unknown')}\n")
            f.write(f"  Total reads:  {rs.get('total_reads', 0):,}\n")
            f.write(f"  Has uracil:   {'Yes' if rs.get('has_uracil', False) else 'No'}\n\n")
        
        if 'alignment_stats' in stats:
            f.write("Alignment Statistics:\n")
            f.write("-"*30 + "\n")
            alns = stats['alignment_stats']
            if 'total_reads' in alns:
                f.write(f"  Total reads:     {alns['total_reads']:,}\n")
            if 'aligned_reads' in alns:
                f.write(f"  Aligned reads:   {alns['aligned_reads']:,}")
                if 'alignment_rate' in alns:
                    f.write(f" ({alns['alignment_rate']})")
                f.write("\n")
            if 'unaligned_reads' in alns:
                f.write(f"  Unaligned reads: {alns['unaligned_reads']:,}\n")
            if 'total_alignments' in alns:
                f.write(f"  Total alignments: {alns['total_alignments']:,}\n")
            f.write("\n")
        
        if 'parameters' in stats:
            f.write("Parameters Used:\n")
            f.write("-"*30 + "\n")
            params = stats['parameters']
            f.write(f"  Mismatches:      {params.get('mismatches', 'N/A')}\n")
            f.write(f"  Max alignments:  {params.get('max_alignments', 'N/A')}\n")
            f.write(f"  Threads:         {params.get('threads', 'N/A')}\n")
            f.write(f"  Keep unmapped:   {params.get('keep_unmapped', False)}\n\n")
        
        f.write("="*60 + "\n")
    
    log(f"Mapping report saved to {report_file}", "INFO")


def run_mapping(args) -> None:
    """Main mapping workflow."""
    # Validate arguments
    validate_arguments(args)
    
    # Setup output directory
    outdir = ensure_outdir(Path(args.outdir))
    
    # Check disk space
    if not check_disk_space(outdir, 5000):  # Need more space for SAM files
        log("Consider freeing up disk space before proceeding", "WARN")
    
    # Determine output file
    sam_out = Path(args.sam_out) if args.sam_out else outdir / f"{args.prefix}.bowtie.sam"
    
    log("="*60, "INFO")
    log(f"Starting mapping for {args.prefix}", "INFO")
    log("="*60, "INFO")
    
    # Collect statistics
    all_stats = {'parameters': {
        'mismatches': args.mismatches,
        'max_alignments': args.k,
        'threads': args.threads,
        'keep_unmapped': args.keep_unmapped
    }}
    
    # Process reads - convert U->T if needed
    reads_to_map = None
    reads2_to_map = None
    
    if args.reads:
        log("Preparing single-end reads...", "INFO")
        reads_path = Path(args.reads)
        validate_file_exists(reads_path, "Input reads")
        
        # Get read statistics
        read_stats = get_read_stats(reads_path)
        all_stats['read_stats'] = read_stats
        
        # Convert U->T if needed
        if read_stats.get('has_uracil', False):
            converted = outdir / f"{args.prefix}.reads.DNA.fastq.gz"
            if convert_uracil_to_thymine_optimized(reads_path, converted):
                reads_to_map = str(converted)
        else:
            reads_to_map = args.reads
            
    elif args.reads_paired:
        log("Preparing paired-end reads...", "INFO")
        r1, r2 = args.reads_paired
        r1_path, r2_path = Path(r1), Path(r2)
        validate_file_exists(r1_path, "Read 1 file")
        validate_file_exists(r2_path, "Read 2 file")
        
        # Get statistics for R1
        read_stats = get_read_stats(r1_path)
        all_stats['read_stats'] = read_stats
        
        # Convert R1
        if read_stats.get('has_uracil', False):
            c1 = outdir / f"{args.prefix}.reads1.DNA.fastq.gz"
            if convert_uracil_to_thymine_optimized(r1_path, c1):
                reads_to_map = str(c1)
        else:
            reads_to_map = r1
        
        # Convert R2
        r2_stats = get_read_stats(r2_path)
        if r2_stats.get('has_uracil', False):
            c2 = outdir / f"{args.prefix}.reads2.DNA.fastq.gz"
            if convert_uracil_to_thymine_optimized(r2_path, c2):
                reads2_to_map = str(c2)
        else:
            reads2_to_map = r2
    
    # Build or verify index
    index_prefix = args.index_prefix
    
    if args.genome:
        log("Building Bowtie index...", "INFO")
        validate_file_exists(Path(args.genome), "Genome FASTA")
        index_prefix = str(outdir / f"{args.prefix}.bowtie_index")
        bowtie_build(args.genome, index_prefix, args.threads)
    elif args.index_prefix:
        log(f"Using existing index: {args.index_prefix}", "INFO")
        if not check_bowtie_index(args.index_prefix):
            raise RuntimeError(f"Bowtie index files not found for prefix: {args.index_prefix}")
    else:
        raise RuntimeError("Either --genome or --index-prefix must be provided")
    
    # Infer format
    fmt = infer_format(reads_to_map, args.fmt)
    log(f"Input format: {fmt}", "INFO")
    
    # Run alignment
    alignment_stats = run_bowtie_alignment(
        index_prefix=index_prefix,
        reads1=reads_to_map,
        reads2=reads2_to_map,
        sam_out=sam_out,
        threads=args.threads,
        mismatches=args.mismatches,
        max_alignments=args.k,
        fmt=fmt,
        keep_unmapped=args.keep_unmapped,
        prefix=args.prefix
    )
    
    all_stats['alignment_stats'] = alignment_stats
    
    # Get additional statistics from SAM
    if sam_out.exists():
        detailed_stats = get_alignment_stats(sam_out, args.threads)
        all_stats['alignment_stats'].update(detailed_stats)
        
        # Optionally convert to BAM
        if hasattr(args, 'create_bam') and args.create_bam:
            bam_out = sam_out.with_suffix('.bam')
            convert_sam_to_bam(sam_out, bam_out, args.threads)
    
    # Generate report
    generate_mapping_report(all_stats, outdir, args.prefix)
    
    log("="*60, "INFO")
    log(f"Mapping completed successfully!", "INFO")
    log(f"Results saved in: {outdir}", "INFO")
    log("="*60, "INFO")