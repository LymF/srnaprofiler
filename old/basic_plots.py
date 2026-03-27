#!/usr/bin/env python3
"""
basic-plots.py - Plot Generation Module for sRNA Profiler
Separate module for generating plots from processed data files
Enables easy reproduction and editing without re-running analysis
"""
import sys
import argparse
import subprocess
from pathlib import Path
from typing import Optional, List, Tuple

from utils import log, ensure_outdir, validate_file_exists


def get_script_directory():
    """Return the directory where the current Python script is located."""
    return Path(__file__).parent.absolute()


def find_r_script(script_name: str, r_scripts_path: str = "") -> Optional[Path]:
    """
    Find R script in the specified path or in the same directory as this script.
    
    Args:
        script_name: Name of the R script file
        r_scripts_path: Custom path to R scripts (optional)
    
    Returns:
        Path to R script or None if not found
    """
    # Try custom path first
    if r_scripts_path:
        custom_path = Path(r_scripts_path) / script_name
        if custom_path.exists():
            return custom_path
    
    # Try same directory as this script
    script_dir = get_script_directory()
    local_path = script_dir / script_name
    if local_path.exists():
        return local_path
    
    # Try parent directory (if running from subdirectory)
    parent_path = script_dir.parent / script_name
    if parent_path.exists():
        return parent_path
    
    log(f"R script not found: {script_name}", "WARN")
    return None


def run_r_script(script_path: Path, args: List[str], output_dir: Path) -> Tuple[bool, str]:
    """
    Run an R script with arguments.
    
    Args:
        script_path: Path to R script
        args: List of arguments for the script
        output_dir: Directory where output will be created
    
    Returns:
        Tuple of (success: bool, message: str)
    """
    if not script_path or not script_path.exists():
        return False, f"Script not found: {script_path}"
    
    cmd = ["Rscript", str(script_path)] + [str(a) for a in args]
    
    try:
        log(f"Running: {' '.join(cmd)}", "DEBUG")
        result = subprocess.run(
            cmd,
            cwd=str(output_dir),
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode == 0:
            return True, f"Success: {script_path.name}"
        else:
            error_msg = result.stderr if result.stderr else "Unknown error"
            return False, f"Failed: {script_path.name} - {error_msg}"
            
    except subprocess.TimeoutExpired:
        return False, f"Timeout: {script_path.name} (>5 min)"
    except Exception as e:
        return False, f"Error: {script_path.name} - {str(e)}"


def generate_base_distribution_plot(
    tsv_file: Path,
    output_prefix: str,
    output_dir: Path,
    ylabel: str = "RPM",
    script_name: str = "calcDistributionPerBase_publication.R",
    r_scripts_path: str = ""
) -> bool:
    """
    Generate base distribution plot.
    
    Args:
        tsv_file: Input TSV file with base distribution data
        output_prefix: Prefix for output files (without extension)
        output_dir: Directory for output files
        ylabel: Label for Y axis
        script_name: Name of R script to use
        r_scripts_path: Custom path to R scripts
    
    Returns:
        True if successful, False otherwise
    """
    if not tsv_file.exists():
        log(f"Input file not found: {tsv_file}", "ERROR")
        return False
    
    script_path = find_r_script(script_name, r_scripts_path)
    if not script_path:
        return False
    
    args = [str(tsv_file), output_prefix, ylabel]
    success, msg = run_r_script(script_path, args, output_dir)
    
    if success:
        log(f"Generated plot: {output_prefix}.pdf", "INFO")
    else:
        log(msg, "ERROR")
    
    return success


def generate_pattern_plot(
    pattern_file: Path,
    output_prefix: str,
    output_dir: Path,
    both_strands: bool = False,
    r_scripts_path: str = ""
) -> bool:
    """
    Generate pattern analysis plot.
    
    Args:
        pattern_file: Input pattern file with Z-scores
        output_prefix: Prefix for output files
        output_dir: Directory for output files
        both_strands: If True, plot both strands; if False, single strand
        r_scripts_path: Custom path to R scripts
    
    Returns:
        True if successful, False otherwise
    """
    if not pattern_file.exists():
        log(f"Input file not found: {pattern_file}", "ERROR")
        return False
    
    script_name = (
        "plotPattern_publication_both_strands.R" if both_strands 
        else "plotPattern_publication.R"
    )
    
    script_path = find_r_script(script_name, r_scripts_path)
    if not script_path:
        return False
    
    args = [str(pattern_file), output_prefix]
    success, msg = run_r_script(script_path, args, output_dir)
    
    if success:
        suffix = ".both_strands" if both_strands else ""
        log(f"Generated plot: {output_prefix}{suffix}-patterns.pdf", "INFO")
    else:
        log(msg, "ERROR")
    
    return success


def generate_all_base_distribution_plots(
    tables_dir: Path,
    plots_dir: Path,
    prefix: str,
    categories: List[str] = None,
    r_scripts_path: str = ""
) -> dict:
    """
    Generate all base distribution plots for specified categories.
    
    Args:
        tables_dir: Directory containing TSV files
        plots_dir: Directory for output plots
        prefix: Project prefix
        categories: List of categories to plot (e.g., ['putative_siRNA', 'putative_piRNA'])
        r_scripts_path: Custom path to R scripts
    
    Returns:
        Dictionary with results for each plot
    """
    if categories is None:
        categories = ['putative_siRNA', 'putative_piRNA', 'ambiguous']
    
    # Define plot configurations
    plot_configs = []
    
    for category in categories:
        # Standard base distribution plot
        tsv_file = tables_dir / f"{prefix}.{category}.distribution.tsv"
        if tsv_file.exists():
            plot_configs.append({
                'type': 'base_distribution',
                'category': category,
                'tsv_file': tsv_file,
                'output_prefix': f"{prefix}.{category}",
                'ylabel': 'RPM',
                'script': 'calcDistributionPerBase_publication.R'
            })
        
        # Sum strand plot
        tsv_sum = tables_dir / f"{prefix}.{category}.distribution_sum.tsv"
        if tsv_sum.exists():
            plot_configs.append({
                'type': 'base_distribution',
                'category': category,
                'tsv_file': tsv_sum,
                'output_prefix': f"{prefix}.{category}.sum",
                'ylabel': 'RPM',
                'script': 'calcDistributionPerBaseSumStrand_publication.R'
            })
        
        # piRNA zoom plot (only for piRNA)
        if 'piRNA' in category:
            if tsv_file.exists():
                plot_configs.append({
                    'type': 'base_distribution',
                    'category': category,
                    'tsv_file': tsv_file,
                    'output_prefix': f"{prefix}.{category}.zoom",
                    'ylabel': 'RPM',
                    'script': 'calcDistributionPerBase_publication_zoom_piRNAs.R'
                })
    
    results = {}
    
    # Sequential processing
    log(f"Generating {len(plot_configs)} plots sequentially...", "INFO")
    for config in plot_configs:
        key = f"{config['category']}_{config.get('script', 'unknown')}"
        success = generate_base_distribution_plot(
            tsv_file=config['tsv_file'],
            output_prefix=config['output_prefix'],
            output_dir=plots_dir,
            ylabel=config.get('ylabel', 'RPM'),
            script_name=config['script'],
            r_scripts_path=r_scripts_path
        )
        results[key] = success
    
    # Summary
    success_count = sum(1 for v in results.values() if v)
    log(f"Plot generation complete: {success_count}/{len(results)} successful", "INFO")
    
    return results


def generate_all_pattern_plots(
    tables_dir: Path,
    plots_dir: Path,
    prefix: str,
    categories: List[str] = None,
    r_scripts_path: str = ""
) -> dict:
    """
    Generate all pattern analysis plots for specified categories.
    
    Args:
        tables_dir: Directory containing pattern files
        plots_dir: Directory for output plots
        prefix: Project prefix
        categories: List of categories to plot
        r_scripts_path: Custom path to R scripts
    
    Returns:
        Dictionary with results for each plot
    """
    if categories is None:
        categories = ['putative_siRNA', 'putative_piRNA', 'ambiguous']
    
    plot_configs = []
    
    for category in categories:
        pattern_file = tables_dir / f"{prefix}.{category}.pattern"
        if pattern_file.exists():
            # Single strand plot
            plot_configs.append({
                'type': 'pattern',
                'category': category,
                'pattern_file': pattern_file,
                'output_prefix': f"{prefix}.{category}",
                'both_strands': False
            })
            
            # Both strands plot
            plot_configs.append({
                'type': 'pattern',
                'category': category,
                'pattern_file': pattern_file,
                'output_prefix': f"{prefix}.{category}",
                'both_strands': True
            })
    
    results = {}
    
    # Sequential processing
    log(f"Generating {len(plot_configs)} pattern plots sequentially...", "INFO")
    for config in plot_configs:
        strand_type = "both" if config['both_strands'] else "single"
        key = f"{config['category']}_pattern_{strand_type}"
        success = generate_pattern_plot(
            pattern_file=config['pattern_file'],
            output_prefix=config['output_prefix'],
            output_dir=plots_dir,
            both_strands=config['both_strands'],
            r_scripts_path=r_scripts_path
        )
        results[key] = success
    
    success_count = sum(1 for v in results.values() if v)
    log(f"Pattern plot generation complete: {success_count}/{len(results)} successful", "INFO")
    
    return results


def discover_plot_files(tables_dir: Path, prefix: str) -> dict:
    """
    Discover available files for plotting.
    
    Args:
        tables_dir: Directory to search for files
        prefix: Project prefix
    
    Returns:
        Dictionary with available files by category
    """
    available = {
        'base_distribution': [],
        'pattern': []
    }
    
    # Look for base distribution files
    for pattern in ['*.distribution.tsv', '*.distribution_sum.tsv']:
        for file in tables_dir.glob(f"{prefix}.{pattern}"):
            available['base_distribution'].append(file)
    
    # Look for pattern files
    for file in tables_dir.glob(f"{prefix}.*.pattern"):
        available['pattern'].append(file)
    
    return available


def main():
    """Standalone execution for basic-plots module."""
    parser = argparse.ArgumentParser(
        description="Generate plots from processed sRNA data files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate all plots for a project
  %(prog)s -p myproject -o results/

  # Generate only base distribution plots
  %(prog)s -p myproject -o results/ --base-only

  # Generate only pattern plots
  %(prog)s -p myproject -o results/ --pattern-only

  # Use custom R scripts directory
  %(prog)s -p myproject -o results/ --r-scripts /path/to/scripts/

  # Specify categories
  %(prog)s -p myproject -o results/ --categories putative_siRNA putative_piRNA
        """
    )
    
    # Required arguments
    parser.add_argument(
        "-p", "--prefix",
        required=True,
        help="Project prefix (same as used in analysis)"
    )
    parser.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory (should contain tables/ subdirectory)"
    )
    
    # Optional arguments
    parser.add_argument(
        "--tables-dir",
        help="Custom path to tables directory (default: outdir/tables)"
    )
    parser.add_argument(
        "--plots-dir",
        help="Custom path to plots directory (default: outdir/plots)"
    )
    parser.add_argument(
        "--r-scripts",
        default="",
        help="Custom path to R scripts directory"
    )
    
    # Plot selection
    parser.add_argument(
        "--base-only",
        action="store_true",
        help="Generate only base distribution plots"
    )
    parser.add_argument(
        "--pattern-only",
        action="store_true",
        help="Generate only pattern plots"
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        help="Specific categories to plot (default: all found)"
    )
    
    # Discovery mode
    parser.add_argument(
        "--discover",
        action="store_true",
        help="Only discover available files, don't generate plots"
    )
    
    args = parser.parse_args()
    
    # Setup directories
    outdir = Path(args.outdir)
    tables_dir = Path(args.tables_dir) if args.tables_dir else outdir / "tables"
    plots_dir = Path(args.plots_dir) if args.plots_dir else outdir / "plots"
    
    # Validate
    if not outdir.exists():
        log(f"Output directory does not exist: {outdir}", "ERROR")
        sys.exit(1)
    
    if not tables_dir.exists():
        log(f"Tables directory does not exist: {tables_dir}", "ERROR")
        sys.exit(1)
    
    # Create plots directory
    ensure_outdir(plots_dir)
    
    log("="*60, "INFO")
    log(f"sRNA Profiler - Plot Generation Module", "INFO")
    log("="*60, "INFO")
    log(f"Project: {args.prefix}", "INFO")
    log(f"Tables:  {tables_dir}", "INFO")
    log(f"Plots:   {plots_dir}", "INFO")
    
    # Discovery mode
    if args.discover:
        log("\nDiscovering available files...", "INFO")
        available = discover_plot_files(tables_dir, args.prefix)
        
        log(f"\nBase distribution files ({len(available['base_distribution'])}):", "INFO")
        for f in sorted(available['base_distribution']):
            log(f"  - {f.name}", "INFO")
        
        log(f"\nPattern files ({len(available['pattern'])}):", "INFO")
        for f in sorted(available['pattern']):
            log(f"  - {f.name}", "INFO")
        
        sys.exit(0)
    
    # Generate plots
    all_results = {}
    
    if not args.pattern_only:
        log("\n--- Base Distribution Plots ---", "INFO")
        base_results = generate_all_base_distribution_plots(
            tables_dir=tables_dir,
            plots_dir=plots_dir,
            prefix=args.prefix,
            categories=args.categories,
            r_scripts_path=args.r_scripts
        )
        all_results.update(base_results)
    
    if not args.base_only:
        log("\n--- Pattern Analysis Plots ---", "INFO")
        pattern_results = generate_all_pattern_plots(
            tables_dir=tables_dir,
            plots_dir=plots_dir,
            prefix=args.prefix,
            categories=args.categories,
            r_scripts_path=args.r_scripts
        )
        all_results.update(pattern_results)
    
    # Final summary
    log("="*60, "INFO")
    total = len(all_results)
    success = sum(1 for v in all_results.values() if v)
    failed = total - success
    
    log(f"Plot generation completed:", "INFO")
    log(f"  Total:      {total}", "INFO")
    log(f"  Successful: {success}", "INFO")
    log(f"  Failed:     {failed}", "INFO")
    log(f"  Output:     {plots_dir}", "INFO")
    log("="*60, "INFO")
    
    # Exit with error if any failed
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()