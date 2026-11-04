#!/usr/bin/env python3
"""
Quality Control module for sRNA Profiler
Handles read QC, trimming, and genome preparation
"""
import sys
import subprocess
import gzip
import os
from pathlib import Path
from typing import Optional

from utils import (
    log, which, ensure_outdir, run_command, validate_arguments,
    convert_uracil_to_thymine_optimized, get_read_stats,
    detect_phred_encoding, validate_file_exists,
    DEFAULT_MIN_LEN, DEFAULT_QTRIM,
    DEFAULT_THREADS, FASTQ_ADAPTERS
)


def adapter_fractions(reads: str, candidates: dict[str, str], n_reads: int = 100_000,
                      probe_len: int = 12) -> dict[str, float]:
    """Fraction of the first `n_reads` reads containing each candidate's probe."""
    probes = {name: seq[:probe_len] for name, seq in candidates.items()}
    counts = {name: 0 for name in candidates}
    total = 0
    opener = gzip.open if str(reads).endswith(".gz") else open
    with opener(reads, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            total += 1
            seq = line.strip().upper()
            for name, probe in probes.items():
                if probe in seq:
                    counts[name] += 1
            if total >= n_reads:
                break
    return {name: (c / total if total else 0.0) for name, c in counts.items()}


def detect_adapter(reads: str, candidates: dict[str, str], n_reads: int = 100_000,
                   min_frac: float = 0.05, probe_len: int = 12) -> Optional[str]:
    """Pick the 3' adapter carried by the reads, or None.

    Returns the full sequence of the candidate whose first `probe_len` nt occur
    in the highest fraction of the first `n_reads` reads, if that fraction
    reaches `min_frac`; ties go to the first candidate in dict order.
    """
    fractions = adapter_fractions(reads, candidates, n_reads, probe_len)
    if not fractions:
        return None
    best = max(fractions, key=lambda n: fractions[n])  # first wins on ties
    return candidates[best] if fractions[best] >= min_frac else None


def _detect_from_reads(reads: str) -> Optional[str]:
    """detect_adapter over FASTQ_ADAPTERS, logging the outcome."""
    found = detect_adapter(reads, FASTQ_ADAPTERS)
    if found is None:
        log("No known 3' adapter detected in the reads", "INFO")
        return None
    name = next(n for n, s in FASTQ_ADAPTERS.items() if s == found)
    frac = adapter_fractions(reads, {name: found})[name]
    log(f"Detected 3' adapter '{name}' ({found}) in {100 * frac:.1f}% of reads", "INFO")
    return found


def fastqc(reads: str, outdir: Path, threads: int) -> None:
    """Run FastQC on reads."""
    if not which("fastqc"):
        log("FastQC not found in PATH; skipping FastQC analysis", "WARN")
        return
    
    try:
        cmd = ["fastqc", "--outdir", str(outdir), "--quiet", "--threads", str(threads), reads]
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
    
    if not adapter:
        adapter = _detect_from_reads(in1)
    if adapter:
        cmd += ["--adapter_sequence", adapter]
    elif in2:
        cmd += ["--detect_adapter_for_pe"]
    # else single-end with no known adapter: fastp auto-detects by itself
    
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
    if not adapter:
        adapter = _detect_from_reads(in1)
        if not adapter:
            adapter = FASTQ_ADAPTERS['illumina_srna']
            log("No adapter detected; cutadapt cannot auto-detect, so the "
                f"illumina_srna default ({adapter}) is a guess", "WARN")
    if in2:
        base_cmd += ["-a", adapter, "-A", adapter]
    else:
        base_cmd += ["-a", adapter]
    
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


_IUPAC_TO_N = str.maketrans("BDHKMRSVWY", "NNNNNNNNNN")


def prepare_genome(genome_path, outdir, prefix):
    """Standardizes the FASTA (uppercase, IUPAC -> N) without tiling.

    miRDeep2 rejects any base outside [ACGTNacgtn], and reference genomes do
    carry ambiguity codes. Returns the output path.
    """
    import os
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{prefix}.genome.fa")
    with open(genome_path) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                fout.write(line.strip().upper().translate(_IUPAC_TO_N) + "\n")
    return out


RANDOM_ENDS_SAMPLE = 20_000
RANDOM_ENDS_LEN = 4


def choose_random_ends(rate0: float, rate4: float) -> int:
    """4 when trimming 4 nt from both ends at least doubles the alignment rate and reaches 30%."""
    return RANDOM_ENDS_LEN if rate4 >= 2 * rate0 and rate4 >= 0.30 else 0


def _open_text(path, mode: str):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _bowtie_rate(index_prefix: str, reads: Path, threads: int, trim: int) -> float:
    """Fraction of reads with an alignment (bowtie1 -v 1 -k 1), trimming `trim` nt from both ends."""
    cmd = ["bowtie", "-v", "1", "-k", "1", "-p", str(threads), "-x", index_prefix,
           "--trim5", str(trim), "--trim3", str(trim), str(reads)]
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                          text=True, check=True)
    processed = aligned = None
    for line in proc.stderr.splitlines():
        if line.startswith("# reads processed:"):
            processed = int(line.split(":")[1])
        elif line.startswith("# reads with at least one"):
            aligned = int(line.split(":")[1].split()[0])
    if not processed or aligned is None:
        raise RuntimeError(f"could not parse bowtie summary: {proc.stderr[-300:]}")
    return aligned / processed


def measure_random_end_rates(trimmed: Path, index_prefix: str, threads: int,
                             n_reads: int = RANDOM_ENDS_SAMPLE) -> tuple[float, float]:
    """Align the first n_reads trimmed reads untrimmed and with 4 nt cut from each end.

    Returns (rate0, rate4).
    """
    if not which("bowtie"):
        raise RuntimeError("bowtie not found in PATH")
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        sample = Path(tmp) / "sample.fastq"
        with _open_text(trimmed, "rt") as src, open(sample, "w") as dst:
            for i, line in enumerate(src):
                if i >= 4 * n_reads:
                    break
                dst.write(line)
        return (_bowtie_rate(index_prefix, sample, threads, 0),
                _bowtie_rate(index_prefix, sample, threads, RANDOM_ENDS_LEN))


def trim_random_ends(paths: list[Path], min_len: int, n: int = RANDOM_ENDS_LEN) -> tuple[int, int]:
    """Remove n nt from both ends of every read, in place (temp file + rename).

    Reads (pairs, when two files are given) shorter than min_len afterwards are dropped from
    all files. Returns (kept, dropped).
    """
    paths = [Path(p) for p in paths]
    tmps = [p.with_name("tmp." + p.name) for p in paths]
    kept = dropped = 0
    srcs = [_open_text(p, "rt") for p in paths]
    dsts = [_open_text(t, "wt") for t in tmps]
    try:
        while True:
            recs = []
            for src in srcs:
                rec = [src.readline() for _ in range(4)]
                recs.append(rec if rec[0] else None)
            if all(r is None for r in recs):
                break
            if any(r is None for r in recs):
                raise ValueError("paired FASTQ files have different numbers of reads")
            cut = [(r[0], r[1].rstrip("\n")[n:-n], r[2], r[3].rstrip("\n")[n:-n]) for r in recs]
            if any(len(c[1]) < min_len for c in cut):
                dropped += 1
                continue
            kept += 1
            for dst, (h, seq, plus, qual) in zip(dsts, cut):
                dst.write(f"{h}{seq}\n{plus}{qual}\n")
    finally:
        for f in srcs + dsts:
            f.close()
    for p, t in zip(paths, tmps):
        os.replace(t, p)
    return kept, dropped


def _detection_index(args, std_genome: Optional[Path], outdir: Path) -> Optional[str]:
    """Bowtie1 index for detection: --index-prefix, else the one `mapping -g` would build."""
    from mapping import bowtie_build, bowtie_index_ready
    if getattr(args, "index_prefix", None):
        return args.index_prefix
    if std_genome is None:
        return None
    prefix = str(outdir / f"{args.prefix}.bowtie_index")
    if bowtie_index_ready(prefix):
        log(f"reusing index {prefix}", "INFO")
    else:
        bowtie_build(str(std_genome), prefix, args.threads)
    return prefix


def resolve_random_ends(args, trimmed: list[Path], std_genome: Optional[Path],
                        outdir: Path) -> dict:
    """Decide (auto) or apply (4) random-end trimming on the trimmed reads; returns stats."""
    mode = str(getattr(args, "random_ends", "auto"))
    stats = {"random_ends": 0, "random_ends_rate0": None, "random_ends_rate4": None}
    if mode == "auto":
        index = _detection_index(args, std_genome, outdir) \
            if (std_genome or getattr(args, "index_prefix", None)) else None
        if index is None:
            log("Random-end detection needs a genome (-g) or --index-prefix; assuming 0", "INFO")
            return stats
        rate0, rate4 = measure_random_end_rates(trimmed[0], index, args.threads)
        stats["random_ends"] = choose_random_ends(rate0, rate4)
        stats["random_ends_rate0"], stats["random_ends_rate4"] = rate0, rate4
        log(f"Random ends: untrimmed {100 * rate0:.2f}% aligned, 4 nt trimmed per end "
            f"{100 * rate4:.2f}% aligned; choice {stats['random_ends']}", "INFO")
    else:
        stats["random_ends"] = int(mode)
    if stats["random_ends"]:
        kept, dropped = trim_random_ends(trimmed, args.min_len, stats["random_ends"])
        stats["random_ends_kept"], stats["random_ends_dropped"] = kept, dropped
        log(f"Removed {stats['random_ends']} nt from both ends: kept {kept:,}, "
            f"dropped {dropped:,} reads shorter than {args.min_len}", "INFO")
    return stats


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
        
        if 'random_ends' in stats:
            f.write("Randomized insert ends:\n")
            f.write("-"*30 + "\n")
            f.write(f"  Trimmed per end: {stats['random_ends']} nt\n")
            if stats.get('random_ends_rate0') is not None:
                f.write(f"  Aligned untrimmed: {100 * stats['random_ends_rate0']:.2f}%\n")
                f.write(f"  Aligned trimmed:   {100 * stats['random_ends_rate4']:.2f}%\n")
            f.write("\n")

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
    
    # Genome first: random-end detection aligns trimmed reads to it
    # Process genome if provided
    std_genome = None
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

        all_stats.update(resolve_random_ends(
            args, [trimmed1] + ([trimmed2] if trimmed2 else []), std_genome, outdir))
        
        # Run FastQC on trimmed reads
        log("Running FastQC on trimmed reads...", "INFO")
        fastqc(str(trimmed1), outdir, args.threads)
        if trimmed2:
            fastqc(str(trimmed2), outdir, args.threads)
        
        log(f"Trimmed reads saved:", "INFO")
        log(f"  R1: {trimmed1}", "INFO")
        if trimmed2:
            log(f"  R2: {trimmed2}", "INFO")
    
    # Generate summary report
    generate_qc_report(all_stats, outdir, args.prefix)
    
    log("="*60, "INFO")
    log(f"QC analysis completed successfully!", "INFO")
    log(f"Results saved in: {outdir}", "INFO")
    log("="*60, "INFO")
