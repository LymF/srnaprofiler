#!/usr/bin/env python3
"""
Plots module for sRNA Profiler
Handles all plot generation using R scripts
"""
import sys
import subprocess
from pathlib import Path
from typing import Optional, List, Dict

from utils import log, which, ensure_outdir, run_command, validate_file_exists
import gzip


def open_sam_file(filepath):
    """Open SAM file (handles gzipped files)."""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def calculate_base_composition_from_sam(
    sam_file: Path,
    chromosome: str = None,
    min_length: int = 15,
    max_length: int = 35
) -> dict:
    """
    Calculate first base composition by size from SAM file.
    This is what the R scripts expect for base distribution plots.
    
    Returns:
        Dictionary: {size: {'+': {'A': count, 'C': count, 'G': count, 'T': count}, 
                            '-': {'A': count, ...}}}
    """
    base_comp = {}
    for size in range(min_length, max_length + 1):
        base_comp[size] = {
            '+': {'A': 0, 'C': 0, 'G': 0, 'T': 0},
            '-': {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        }
    
    try:
        seen_reads = set()
        
        with open_sam_file(sam_file) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                
                read_id = parts[0]
                chrom = parts[2]
                flag = int(parts[1])
                seq = parts[9].upper()
                
                # Skip duplicates
                if read_id in seen_reads:
                    continue
                seen_reads.add(read_id)
                
                # Filter by chromosome if specified
                if chromosome and chrom != chromosome:
                    continue
                
                # Skip unmapped
                if chrom == '*':
                    continue
                
                # Get read length
                read_len = len(seq)
                if not (min_length <= read_len <= max_length):
                    continue
                
                # Get first base
                # If reverse strand, get reverse complement first
                if flag & 0x10:  # Reverse strand
                    first_base = reverse_complement(seq)[0] if seq else 'N'
                    strand = '-'
                else:  # Forward strand
                    first_base = seq[0] if seq else 'N'
                    strand = '+'
                
                # Count if valid base
                if first_base in ['A', 'C', 'G', 'T']:
                    base_comp[read_len][strand][first_base] += 1
        
        return base_comp
        
    except Exception as e:
        log(f"Error calculating base composition: {e}", "ERROR")
        return base_comp


def write_base_composition_tsv(
    base_comp: dict,
    output_file: Path
) -> bool:
    """
    Write base composition to TSV format expected by R scripts.
    
    Format (each size appears twice - once for + strand, once for - strand):
    size	A	C	G	U
    15	100	50	30	20
    15	-50	-25	-15	-10
    16	120	60	35	25
    16	-60	-30	-18	-12
    """
    try:
        ensure_outdir(output_file.parent)
        
        with open(output_file, 'w') as f:
            # Header (note: U not T, for R script compatibility)
            f.write("size\tA\tC\tG\tU\n")
            
            for size in sorted(base_comp.keys()):
                # Positive strand (forward)
                a_pos = base_comp[size]['+']['A']
                c_pos = base_comp[size]['+']['C']
                g_pos = base_comp[size]['+']['G']
                t_pos = base_comp[size]['+']['T']
                f.write(f"{size}\t{a_pos}\t{c_pos}\t{g_pos}\t{t_pos}\n")
                
                # Negative strand (reverse) - negative values
                a_neg = base_comp[size]['-']['A'] * -1
                c_neg = base_comp[size]['-']['C'] * -1
                g_neg = base_comp[size]['-']['G'] * -1
                t_neg = base_comp[size]['-']['T'] * -1
                f.write(f"{size}\t{a_neg}\t{c_neg}\t{g_neg}\t{t_neg}\n")
        
        if output_file.exists():
            log(f"Created base composition TSV: {output_file}", "DEBUG")
            return True
        return False
        
    except Exception as e:
        log(f"Error writing base composition TSV: {e}", "ERROR")
        return False


def get_r_scripts_path() -> Path:
    """Return the path to R scripts directory."""
    script_dir = Path(__file__).parent.absolute()
    
    # Check common locations for R scripts
    possible_paths = [
        script_dir / "R",  # R subdirectory (MOST COMMON)
        script_dir,        # Same directory
        script_dir / "scripts",
        script_dir / "r_scripts",
        script_dir / "R_scripts",
        Path("/mnt/scripts"),
        Path("/opt/scripts"),
    ]
    
    log(f"Searching for R scripts starting from: {script_dir}", "DEBUG")
    
    for path in possible_paths:
        if path.exists():
            log(f"  Checking: {path}", "DEBUG")
            # Check for at least one required R script
            required_scripts = [
                "calcDistributionPerBase_publication.R",
                "plotPattern_publication.R"
            ]
            for script in required_scripts:
                if (path / script).exists():
                    log(f"Found R scripts directory: {path}", "INFO")
                    return path
        else:
            log(f"  Path does not exist: {path}", "DEBUG")
    
    # If not found, try R subdirectory as default
    default_path = script_dir / "R"
    log(f"R scripts not found in common locations. Using default: {default_path}", "WARN")
    log(f"Please ensure R scripts are in one of these locations:", "WARN")
    for p in possible_paths:
        log(f"  - {p}", "WARN")
    
    return default_path


def check_r_dependencies() -> bool:
    """Check if R and required packages are available."""
    if not which("Rscript"):
        log("Rscript not found in PATH", "WARN")
        return False
    
    # Check for required R packages
    required_packages = ["ggplot2", "reshape2", "plyr"]
    
    test_script = """
    packages <- c('ggplot2', 'reshape2', 'plyr')
    missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) {
        cat(paste(missing, collapse=','))
        quit(status=1)
    }
    """
    
    try:
        result = subprocess.run(
            ["Rscript", "-e", test_script],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode != 0:
            missing = result.stdout.strip()
            log(f"Missing R packages: {missing}", "WARN")
            log(f"Install with: install.packages(c({missing}))", "INFO")
            return False
        
        return True
        
    except Exception as e:
        log(f"Error checking R dependencies: {e}", "DEBUG")
        return False


def run_r_script(
    script_name: str,
    args: List[str],
    r_scripts_path: Optional[Path] = None,
    timeout: int = 300,
    cwd: Optional[Path] = None
) -> bool:
    """
    Run an R script with given arguments.
    
    Args:
        script_name: Name of the R script (e.g., 'plotPattern_publication.R')
        args: List of arguments to pass to the script
        r_scripts_path: Path to R scripts directory
        timeout: Timeout in seconds
        cwd: Working directory for R script execution
    
    Returns:
        True if successful, False otherwise
    """
    if not which("Rscript"):
        log("Rscript not found in PATH", "ERROR")
        return False
    
    # Determine R scripts path
    if r_scripts_path is None:
        r_scripts_path = get_r_scripts_path()
    
    script_path = r_scripts_path / script_name
    
    if not script_path.exists():
        log(f"R script not found: {script_path}", "ERROR")
        return False
    
    log(f"Running R script: {script_name}", "INFO")
    log(f"  Script path: {script_path}", "DEBUG")
    log(f"  Arguments: {args}", "DEBUG")
    if cwd:
        log(f"  Working directory: {cwd}", "DEBUG")
    
    cmd = ["Rscript", str(script_path)] + [str(arg) for arg in args]
    
    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False
        )
        
        if result.returncode != 0:
            log(f"R script failed with return code {result.returncode}", "ERROR")
            if result.stderr:
                log(f"R stderr: {result.stderr}", "ERROR")
            if result.stdout:
                log(f"R stdout: {result.stdout}", "DEBUG")
            return False
        
        if result.stdout:
            log(f"R output: {result.stdout.strip()}", "DEBUG")
        
        return True
        
    except subprocess.TimeoutExpired:
        log(f"R script timed out after {timeout}s", "ERROR")
        return False
    except Exception as e:
        log(f"Error running R script: {e}", "ERROR")
        return False


def generate_base_distribution_plot(
    input_tsv: Path,
    output_prefix: Path,
    r_scripts_path: Optional[Path] = None,
    plot_type: str = "per_base"
) -> bool:
    """
    Generate base distribution plots from TSV file.
    
    Args:
        input_tsv: Path to input TSV file with distribution data
        output_prefix: Output prefix (without extension)
        r_scripts_path: Path to R scripts directory
        plot_type: Type of plot - 'per_base', 'per_base_zoom', or 'sum_strand'
    
    Returns:
        True if successful, False otherwise
    """
    validate_file_exists(input_tsv, "Distribution TSV file")
    
    # Select appropriate R script
    script_map = {
        "per_base": "calcDistributionPerBase_publication.R",
        "per_base_zoom": "calcDistributionPerBase_publication_zoom_piRNAs.R",
        "sum_strand": "calcDistributionPerBaseSumStrand_publication.R"
    }
    
    script_name = script_map.get(plot_type)
    if not script_name:
        log(f"Unknown plot type: {plot_type}", "ERROR")
        return False
    
    # Ensure output directory exists
    ensure_outdir(output_prefix.parent)
    
    # Run R script
    # These scripts expect: input_file, output_prefix
    args = [str(input_tsv), str(output_prefix)]
    
    success = run_r_script(script_name, args, r_scripts_path)
    
    if success:
        log(f"Generated base distribution plot: {output_prefix}.pdf", "INFO")
    
    return success


def generate_pattern_plot(
    input_matrix: Path,
    output_prefix: Path,
    r_scripts_path: Optional[Path] = None,
    both_strands: bool = False
) -> bool:
    """
    Generate pattern plots (Z-score) from matrix file.
    
    Args:
        input_matrix: Path to input matrix CSV file
        output_prefix: Output prefix (without extension)
        r_scripts_path: Path to R scripts directory
        both_strands: If True, plot both strands separately
    
    Returns:
        True if successful, False otherwise
    """
    validate_file_exists(input_matrix, "Pattern matrix file")
    
    # Select appropriate R script
    script_name = ("plotPattern_publication_both_strands.R" if both_strands 
                   else "plotPattern_publication.R")
    
    # Ensure output directory exists
    ensure_outdir(output_prefix.parent)
    
    # Run R script
    # These scripts expect: input_matrix, output_prefix
    args = [str(input_matrix), str(output_prefix)]
    
    success = run_r_script(script_name, args, r_scripts_path)
    
    if success:
        log(f"Generated pattern plot: {output_prefix}.pdf", "INFO")
    
    return success


def generate_all_base_distribution_plots(
    tables_dir: Path,
    plots_dir: Path,
    prefix: str,
    categories: List[str],
    r_scripts_path: Optional[Path] = None
) -> Dict[str, bool]:
    """
    Generate all base distribution plots for given categories.
    
    Args:
        tables_dir: Directory containing TSV files
        plots_dir: Output directory for plots
        prefix: Output file prefix
        categories: List of categories (e.g., ['putative_siRNA', 'putative_piRNA'])
        r_scripts_path: Path to R scripts directory
    
    Returns:
        Dictionary mapping plot names to success status
    """
    log("Generating base distribution plots...", "INFO")
    
    # Check R availability
    if not check_r_dependencies():
        log("R dependencies not available, skipping plots", "WARN")
        return {}
    
    ensure_outdir(plots_dir)
    results = {}
    
    for category in categories:
        tsv_file = tables_dir / f"{prefix}.{category}.distribution.tsv"
        
        if not tsv_file.exists():
            log(f"TSV file not found: {tsv_file.name}", "WARN")
            results[f"{category}_per_base"] = False
            continue
        
        # Generate per-base plot
        output_prefix = plots_dir / f"{prefix}.{category}.distribution"
        success = generate_base_distribution_plot(
            input_tsv=tsv_file,
            output_prefix=output_prefix,
            r_scripts_path=r_scripts_path,
            plot_type="per_base"
        )
        results[f"{category}_per_base"] = success
        
        # Generate zoomed plot for piRNAs if applicable
        if "piRNA" in category and success:
            output_prefix_zoom = plots_dir / f"{prefix}.{category}.distribution_zoom"
            success_zoom = generate_base_distribution_plot(
                input_tsv=tsv_file,
                output_prefix=output_prefix_zoom,
                r_scripts_path=r_scripts_path,
                plot_type="per_base_zoom"
            )
            results[f"{category}_per_base_zoom"] = success_zoom
        
        # Generate sum strand plot
        if success:
            output_prefix_sum = plots_dir / f"{prefix}.{category}.distribution_sumstrand"
            success_sum = generate_base_distribution_plot(
                input_tsv=tsv_file,
                output_prefix=output_prefix_sum,
                r_scripts_path=r_scripts_path,
                plot_type="sum_strand"
            )
            results[f"{category}_sum_strand"] = success_sum
    
    success_count = sum(1 for v in results.values() if v)
    log(f"Generated {success_count}/{len(results)} base distribution plots", "INFO")
    
    return results


def generate_all_pattern_plots(
    tables_dir: Path,
    plots_dir: Path,
    prefix: str,
    categories: List[str],
    r_scripts_path: Optional[Path] = None
) -> Dict[str, bool]:
    """
    Generate all pattern plots for given categories.
    
    Args:
        tables_dir: Directory containing pattern matrix files
        plots_dir: Output directory for plots
        prefix: Output file prefix
        categories: List of categories (e.g., ['putative_siRNA', 'putative_piRNA'])
        r_scripts_path: Path to R scripts directory
    
    Returns:
        Dictionary mapping plot names to success status
    """
    log("Generating pattern plots...", "INFO")
    
    # Check R availability
    if not check_r_dependencies():
        log("R dependencies not available, skipping plots", "WARN")
        return {}
    
    ensure_outdir(plots_dir)
    results = {}
    
    for category in categories:
        # Look for zscore matrix file
        matrix_file = tables_dir / f"{prefix}.{category}.zscore_matrix.csv"
        
        # Also check for .pattern extension (legacy)
        if not matrix_file.exists():
            matrix_file = tables_dir / f"{prefix}.{category}.pattern"
        
        if not matrix_file.exists():
            log(f"Pattern matrix not found for {category}", "WARN")
            results[f"{category}_pattern"] = False
            continue
        
        # Generate combined strand plot
        output_prefix = plots_dir / f"{prefix}.{category}.pattern"
        success = generate_pattern_plot(
            input_matrix=matrix_file,
            output_prefix=output_prefix,
            r_scripts_path=r_scripts_path,
            both_strands=False
        )
        results[f"{category}_pattern"] = success
        
        # Generate both strands plot
        if success:
            output_prefix_both = plots_dir / f"{prefix}.{category}.pattern_both_strands"
            success_both = generate_pattern_plot(
                input_matrix=matrix_file,
                output_prefix=output_prefix_both,
                r_scripts_path=r_scripts_path,
                both_strands=True
            )
            results[f"{category}_pattern_both"] = success_both
    
    success_count = sum(1 for v in results.values() if v)
    log(f"Generated {success_count}/{len(results)} pattern plots", "INFO")
    
    return results


def generate_all_plots(
    tables_dir: Path,
    plots_dir: Path,
    prefix: str,
    categories: List[str],
    r_scripts_path: Optional[Path] = None,
    include_base_distribution: bool = True,
    include_patterns: bool = True
) -> Dict[str, Dict[str, bool]]:
    """
    Generate all plots (base distribution and patterns).
    
    Args:
        tables_dir: Directory containing data files
        plots_dir: Output directory for plots
        prefix: Output file prefix
        categories: List of categories to plot
        r_scripts_path: Path to R scripts directory
        include_base_distribution: Generate base distribution plots
        include_patterns: Generate pattern plots
    
    Returns:
        Dictionary with 'base_distribution' and 'patterns' keys containing results
    """
    log("="*60, "INFO")
    log("Starting plot generation", "INFO")
    log("="*60, "INFO")
    
    results = {}
    
    if include_base_distribution:
        results['base_distribution'] = generate_all_base_distribution_plots(
            tables_dir=tables_dir,
            plots_dir=plots_dir,
            prefix=prefix,
            categories=categories,
            r_scripts_path=r_scripts_path
        )
    
    if include_patterns:
        results['patterns'] = generate_all_pattern_plots(
            tables_dir=tables_dir,
            plots_dir=plots_dir,
            prefix=prefix,
            categories=categories,
            r_scripts_path=r_scripts_path
        )
    
    # Summary
    total_success = sum(
        sum(1 for v in cat_results.values() if v)
        for cat_results in results.values()
    )
    total_attempted = sum(
        len(cat_results)
        for cat_results in results.values()
    )
    
    log("="*60, "INFO")
    log(f"Plot generation completed: {total_success}/{total_attempted} successful", "INFO")
    log("="*60, "INFO")
    
    return results


def get_chromosomes_from_sam(sam_file: Path) -> list:
    """Extract chromosome names from SAM header."""
    chromosomes = []
    
    try:
        with open_sam_file(sam_file) as f:
            for line in f:
                if line.startswith('@SQ'):
                    parts = line.strip().split('\t')
                    for part in parts:
                        if part.startswith('SN:'):
                            chrom = part[3:]
                            chromosomes.append(chrom)
                elif not line.startswith('@'):
                    break
        
        log(f"Found {len(chromosomes)} chromosomes in SAM", "INFO")
        return chromosomes
        
    except Exception as e:
        log(f"Error reading SAM: {e}", "ERROR")
        return []


def generate_per_chromosome_plots(
    sam_file: Path,
    output_dir: Path,
    prefix: str,
    chromosomes: list = None,
    min_length: int = 15,
    max_length: int = 35,
    r_scripts_path: Path = None
) -> dict:
    """
    Generate plots for each chromosome separately.
    
    Args:
        sam_file: Path to SAM file
        output_dir: Base output directory
        prefix: Output file prefix
        chromosomes: List of chromosomes (auto-detect if None)
        min_length: Minimum read length
        max_length: Maximum read length
        r_scripts_path: Path to R scripts
    
    Returns:
        Dictionary with success status for each chromosome
    """
    log("="*60, "INFO")
    log("Per-Chromosome Plot Generation", "INFO")
    log("="*60, "INFO")
    
    # Check dependencies
    if not check_r_dependencies():
        log("R dependencies not met, cannot generate plots", "ERROR")
        return {}
    
    # Get chromosomes
    if chromosomes is None:
        chromosomes = get_chromosomes_from_sam(sam_file)
        if not chromosomes:
            log("No chromosomes found in SAM file", "ERROR")
            return {}
    
    # Create output structure
    # Check if output_dir already ends with 'plots'
    if output_dir.name == "plots":
        plots_dir = output_dir / "per_chromosome"
    else:
        plots_dir = output_dir / "plots" / "per_chromosome"
    
    temp_dir = plots_dir / ".temp"
    ensure_outdir(temp_dir)
    
    log(f"Output directory: {plots_dir}", "DEBUG")
    
    results = {}
    total = len(chromosomes)
    
    log(f"Processing {total} chromosomes...\n", "INFO")
    
    for idx, chrom in enumerate(chromosomes, 1):
        log(f"[{idx}/{total}] Processing {chrom}...", "INFO")
        
        # Calculate base composition
        base_comp = calculate_base_composition_from_sam(
            sam_file=sam_file,
            chromosome=chrom,
            min_length=min_length,
            max_length=max_length
        )
        
        # Check for data
        total_reads = sum(
            base_comp[size][strand][base]
            for size in base_comp
            for strand in ['+', '-']
            for base in ['A', 'C', 'G', 'T']
        )
        
        if total_reads == 0:
            log(f"  No reads found for {chrom}, skipping", "WARN")
            results[chrom] = False
            continue
        
        log(f"  Found {total_reads:,} reads", "INFO")
        
        # Write TSV
        tsv_file = temp_dir / f"{prefix}.{chrom}.tsv"
        if not write_base_composition_tsv(base_comp, tsv_file):
            results[chrom] = False
            continue
        
        # Create chromosome output directory
        chrom_dir = plots_dir / chrom
        ensure_outdir(chrom_dir)
        
        # Track success for this chromosome
        chrom_results = {}
        
        # Generate 3 types of plots
        
        # 1. Per-base plot (both strands)
        log(f"  Generating per-base plot...", "INFO")
        output_prefix_base = chrom_dir / f"{prefix}.{chrom}"
        success_base = generate_base_distribution_plot(
            input_tsv=tsv_file,
            output_prefix=output_prefix_base,
            r_scripts_path=r_scripts_path,
            plot_type="per_base"
        )
        chrom_results['per_base'] = success_base
        
        # 2. Zoomed plot (piRNA range)
        log(f"  Generating zoomed plot...", "INFO")
        output_prefix_zoom = chrom_dir / f"{prefix}.{chrom}_zoom"
        success_zoom = generate_base_distribution_plot(
            input_tsv=tsv_file,
            output_prefix=output_prefix_zoom,
            r_scripts_path=r_scripts_path,
            plot_type="per_base_zoom"
        )
        chrom_results['zoom'] = success_zoom
        
        # 3. Sum strands plot
        log(f"  Generating sum strands plot...", "INFO")
        output_prefix_sum = chrom_dir / f"{prefix}.{chrom}_sumstrand"
        success_sum = generate_base_distribution_plot(
            input_tsv=tsv_file,
            output_prefix=output_prefix_sum,
            r_scripts_path=r_scripts_path,
            plot_type="sum_strand"
        )
        chrom_results['sum_strand'] = success_sum
        
        # Mark as successful if at least one plot succeeded
        results[chrom] = any(chrom_results.values())
        
        # Summary for this chromosome
        successful_plots = sum(1 for v in chrom_results.values() if v)
        if successful_plots == 3:
            log(f"  ✓ All 3 plots generated", "INFO")
        elif successful_plots > 0:
            log(f"  ⚠ {successful_plots}/3 plots generated", "WARN")
        else:
            log(f"  ✗ All plots failed", "WARN")
    
    # Cleanup
    try:
        import shutil
        shutil.rmtree(temp_dir)
    except:
        pass
    
    # Summary
    log("\n" + "="*60, "INFO")
    successful = sum(1 for v in results.values() if v)
    log(f"Summary: {successful}/{total} chromosomes successful", "INFO")
    log("="*60, "INFO")
    
    return results


def run_from_sam(args) -> None:
    """Generate plots directly from SAM file without requiring pre-computed tables."""
    sam_file = Path(args.sam)
    if not sam_file.exists():
        log(f"SAM file not found: {sam_file}", "ERROR")
        raise FileNotFoundError(f"SAM file not found: {sam_file}")
    
    output_dir = Path(args.outdir)
    ensure_outdir(output_dir)
    
    # Determine mode
    per_chrom = getattr(args, 'per_chrom', False)
    
    log("="*60, "INFO")
    log("Generating plots directly from SAM file", "INFO")
    log(f"SAM: {sam_file}", "INFO")
    log(f"Output: {output_dir}", "INFO")
    log(f"Mode: {'PER-CHROMOSOME' if per_chrom else 'COMBINED'}", "INFO")
    log("="*60, "INFO")
    
    # Check R dependencies
    if not check_r_dependencies():
        log("R dependencies not met", "ERROR")
        raise RuntimeError("R dependencies not met")
    
    # Determine R scripts path
    r_scripts_path = Path(args.r_scripts_path) if hasattr(args, 'r_scripts_path') and args.r_scripts_path else None
    
    min_len = getattr(args, 'min_len', 15)
    max_len = getattr(args, 'max_len', 35)
    
    if per_chrom:
        # Per-chromosome mode (using internal function)
        chromosomes = getattr(args, 'chromosomes', None)
        
        results = generate_per_chromosome_plots(
            sam_file=sam_file,
            output_dir=output_dir,
            prefix=args.prefix,
            chromosomes=chromosomes,
            min_length=min_len,
            max_length=max_len,
            r_scripts_path=r_scripts_path
        )
        
        # Summary
        successful = sum(1 for v in results.values() if v)
        total = len(results)
        log(f"\nResults: {successful}/{total} chromosomes successful", "INFO")
        
    else:
        # Combined mode
        log("Calculating base composition from SAM...", "INFO")
        base_comp = calculate_base_composition_from_sam(
            sam_file=sam_file,
            chromosome=None,
            min_length=min_len,
            max_length=max_len
        )
        
        # Check for data
        total_reads = sum(
            base_comp[size][strand][base]
            for size in base_comp
            for strand in ['+', '-']
            for base in ['A', 'C', 'G', 'T']
        )
        
        if total_reads == 0:
            log("No reads found in SAM file", "ERROR")
            raise ValueError("No reads found in SAM file")
        
        log(f"Found {total_reads:,} total reads", "INFO")
        
        # Setup output directories
        plots_dir = output_dir / "plots"
        tables_dir = output_dir / "tables"
        ensure_outdir(plots_dir)
        ensure_outdir(tables_dir)
        
        # Write TSV
        tsv_file = tables_dir / f"{args.prefix}.Geral_base_distribution"
        if not write_base_composition_tsv(base_comp, tsv_file):
            raise RuntimeError("Failed to write base composition TSV")
        
        log(f"Wrote TSV: {tsv_file}", "INFO")
        
        # Generate multiple plot types
        output_prefix_base = plots_dir / args.prefix
        plot_results = {}
        
        # 1. Per-base plot (both strands)
        log("Generating per-base plot...", "INFO")
        success_base = generate_base_distribution_plot(
            input_tsv=tsv_file,
            output_prefix=output_prefix_base,
            r_scripts_path=r_scripts_path,
            plot_type="per_base"
        )
        plot_results['per_base'] = success_base
        
        # 2. Zoomed plot (piRNA range)
        log("Generating zoomed plot (piRNA range)...", "INFO")
        output_prefix_zoom = plots_dir / f"{args.prefix}_zoom"
        success_zoom = generate_base_distribution_plot(
            input_tsv=tsv_file,
            output_prefix=output_prefix_zoom,
            r_scripts_path=r_scripts_path,
            plot_type="per_base_zoom"
        )
        plot_results['per_base_zoom'] = success_zoom
        
        # 3. Sum strands plot
        log("Generating sum strands plot...", "INFO")
        output_prefix_sum = plots_dir / f"{args.prefix}_sumstrand"
        success_sum = generate_base_distribution_plot(
            input_tsv=tsv_file,
            output_prefix=output_prefix_sum,
            r_scripts_path=r_scripts_path,
            plot_type="sum_strand"
        )
        plot_results['sum_strand'] = success_sum
        
        # Summary
        successful = sum(1 for v in plot_results.values() if v)
        total = len(plot_results)
        
        log("="*60, "INFO")
        log(f"Plot generation summary: {successful}/{total} successful", "INFO")
        for plot_type, success in plot_results.items():
            status = "✓" if success else "✗"
            log(f"  {status} {plot_type}", "INFO")
        log("="*60, "INFO")
        
        if successful == 0:
            log("✗ All plots failed", "ERROR")
            raise RuntimeError("Plot generation failed")


def run_plots(args) -> None:
    """Main function for plots module - can be called from srnaprofiler.py"""
    # Setup paths
    tables_dir = Path(args.tables_dir) if hasattr(args, 'tables_dir') and args.tables_dir else Path(args.outdir) / "tables"
    plots_dir = Path(args.plots_dir) if hasattr(args, 'plots_dir') and args.plots_dir else Path(args.outdir) / "plots"
    
    # Validate tables directory
    if not tables_dir.exists():
        log(f"Tables directory not found: {tables_dir}", "ERROR")
        raise FileNotFoundError(f"Tables directory not found: {tables_dir}")
    
    # Ensure plots directory
    ensure_outdir(plots_dir)
    
    # Determine categories to plot
    categories = []
    if hasattr(args, 'categories') and args.categories:
        categories = args.categories
    else:
        # Auto-detect from available files
        log("Auto-detecting categories from available files...", "INFO")
        for category in ['putative_siRNA', 'putative_piRNA', 'ambiguous']:
            # Check for distribution TSV
            tsv_file = tables_dir / f"{args.prefix}.{category}.distribution.tsv"
            if tsv_file.exists():
                categories.append(category)
                log(f"Found data for: {category}", "INFO")
    
    if not categories:
        log("No data files found for plotting", "WARN")
        log(f"Looked in: {tables_dir}", "INFO")
        return
    
    log(f"Will generate plots for: {', '.join(categories)}", "INFO")
    
    # Determine R scripts path
    r_scripts_path = None
    if hasattr(args, 'r_scripts_path') and args.r_scripts_path:
        r_scripts_path = Path(args.r_scripts_path)
    
    # Generate plots
    results = generate_all_plots(
        tables_dir=tables_dir,
        plots_dir=plots_dir,
        prefix=args.prefix,
        categories=categories,
        r_scripts_path=r_scripts_path,
        include_base_distribution=getattr(args, 'base_distribution', True),
        include_patterns=getattr(args, 'patterns', True)
    )
    
    # Report
    log("Plot generation results:", "INFO")
    for plot_type, type_results in results.items():
        success = sum(1 for v in type_results.values() if v)
        total = len(type_results)
        log(f"  {plot_type}: {success}/{total} successful", "INFO")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate plots for sRNA Profiler analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes of operation:
  1. From tables (default): Use pre-computed TSV files from clusters analysis
  2. From SAM: Generate plots directly from SAM file

Examples:
  # Mode 1: From pre-computed tables (fast)
  %(prog)s -p sample1 -o results/
  
  # Mode 2: From SAM file (generates distribution plots)
  %(prog)s --sam mapped.sam -p sample1 -o results/
  
  # Mode 2 + per-chromosome
  %(prog)s --sam mapped.sam -p sample1 -o results/ --per-chrom
  
  # Mode 2 + specific chromosomes
  %(prog)s --sam mapped.sam -p sample1 -o results/ --per-chrom --chromosomes chr1 chr2
"""
    )
    
    # Required arguments
    parser.add_argument(
        "-p", "--prefix",
        required=True,
        help="Output file prefix"
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    
    # SAM mode arguments
    parser.add_argument(
        "--sam",
        help="SAM file (enables direct SAM mode)"
    )
    parser.add_argument(
        "--per-chrom",
        action="store_true",
        help="Generate separate plots for each chromosome (SAM mode only)"
    )
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        help="Specific chromosomes to plot (with --per-chrom)"
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=15,
        help="Minimum read length (default: 15)"
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=35,
        help="Maximum read length (default: 35)"
    )
    
    # Table mode arguments
    parser.add_argument(
        "--tables-dir",
        help="Directory containing data tables (default: <outdir>/tables)"
    )
    parser.add_argument(
        "--plots-dir",
        help="Output directory for plots (default: <outdir>/plots)"
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        help="Categories to plot (default: auto-detect)"
    )
    
    # Common arguments
    parser.add_argument(
        "--r-scripts-path",
        help="Path to R scripts directory (auto-detect if not specified)"
    )
    parser.add_argument(
        "--no-base-distribution",
        action="store_false",
        dest="base_distribution",
        help="Skip base distribution plots"
    )
    parser.add_argument(
        "--no-patterns",
        action="store_false",
        dest="patterns",
        help="Skip pattern plots"
    )
    
    args = parser.parse_args()
    
    try:
        # Determine mode based on presence of --sam
        if args.sam:
            # SAM mode
            run_from_sam(args)
        else:
            # Table mode (original behavior)
            run_plots(args)
            
    except Exception as e:
        log(f"Plot generation failed: {e}", "ERROR")
        import traceback
        traceback.print_exc()
        sys.exit(1)