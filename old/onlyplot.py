#!/usr/bin/env python3
"""
onlyplot.py - Plot generation module for sRNA Profiler
Generate plots from existing files without re-running analysis
"""
import sys
import os
import argparse
from pathlib import Path
from typing import Optional

from utils import log, ensure_outdir, validate_file_exists
from clusters import PlotGeralDistribution


def get_sam_basename(sam_path: str) -> str:
    """Extract base name from SAM file path for use in plot names."""
    sam_file = Path(sam_path)
    # Remove extensions and get the base name
    basename = sam_file.stem
    if sam_file.suffix == '.gz':
        # Handle .sam.gz files
        basename = Path(sam_file.stem).stem
    return basename


def run_onlyplot(args):
    """Generate plots from existing SAM file without re-running analysis."""
    
    # Validate input files
    validate_file_exists(Path(args.sam), "SAM file")
    validate_file_exists(Path(args.genome), "Genome file")
    if args.reads:
        validate_file_exists(Path(args.reads), "Reads file")
    
    # Setup output directory
    outdir = ensure_outdir(Path(args.outdir))
    
    log("="*60, "INFO")
    log(f"Starting plot generation for {args.prefix}", "INFO")
    log("="*60, "INFO")
    
    log(f"Input files:", "INFO")
    log(f"  SAM:     {args.sam}", "INFO")
    log(f"  Genome:  {args.genome}", "INFO")
    log(f"  Reads:   {args.reads or 'Not provided'}", "INFO")
    
    # Create directory structure
    tables_dir = ensure_outdir(outdir / "tables")
    plots_dir = ensure_outdir(outdir / "plots")
    profiles_dir = ensure_outdir(outdir / "profiles")
    
    # Generate base distribution plots
    log("Generating base distribution plots...", "INFO")
    
    # Use SAM file basename as suffix
    sam_basename = get_sam_basename(args.sam)
    analyses = [(args.sam, sam_basename, f"SAM: {sam_basename}")]
    
    # Check if separated SAM files exist (optional)
    sam_dir = outdir / "sam"
    if sam_dir.exists():
        si_sam = sam_dir / f"{args.prefix}.putative_siRNA.sam"
        pi_sam = sam_dir / f"{args.prefix}.putative_piRNA.sam"
        amb_sam = sam_dir / f"{args.prefix}.ambiguous.sam"
        
        if si_sam.exists() and si_sam.stat().st_size > 0:
            analyses.append((str(si_sam), "putative_siRNA", "putative_siRNA reads"))
            log(f"Found siRNA SAM: {si_sam}", "INFO")
        
        if pi_sam.exists() and pi_sam.stat().st_size > 0:
            analyses.append((str(pi_sam), "putative_piRNA", "putative_piRNA reads"))
            log(f"Found piRNA SAM: {pi_sam}", "INFO")
        
        if amb_sam.exists() and amb_sam.stat().st_size > 0:
            analyses.append((str(amb_sam), "ambiguous", "Ambiguous reads"))
            log(f"Found ambiguous SAM: {amb_sam}", "INFO")
    
    # Generate plots for each analysis
    for sam_file, suffix, desc in analyses:
        log(f"Generating plots for {desc}...", "INFO")
        basename = f"{args.prefix}.{suffix}"
        
        # Use default values for plot generation
        plotter = PlotGeralDistribution(
            sam_file=sam_file,
            basename=basename,
            tables_dir=tables_dir,
            profiles_dir=profiles_dir,
            plots_dir=plots_dir,
            start=15,  # Default min length
            end=35,    # Default max length
            norm=getattr(args, 'norm', None),
            minimo=30,  # Default minimum
            plot=True,
            antisense=False,
            r_scripts_path=""
        )
        plotter.analyze()
    
    # Generate pattern plots if requested
    if getattr(args, 'run_patterns', False):
        log("Generating pattern analysis plots...", "INFO")
        # Add pattern analysis plot generation here if needed
        # This would use the existing classification.tsv and other results
        pass
    
    log("="*60, "INFO")
    log(f"Plot generation completed successfully!", "INFO")
    log(f"Plots saved in: {plots_dir}", "INFO")
    log("="*60, "INFO")


def main():
    """Standalone execution for onlyplot module."""
    parser = argparse.ArgumentParser(
        description="Generate plots from existing SAM file without re-running analysis"
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
    parser.add_argument(
        "--sam",
        required=True,
        help="SAM/BAM file with mapped reads (can be gzipped)"
    )
    parser.add_argument(
        "-g", "--genome",
        required=True,
        help="Reference genome FASTA"
    )
    
    # Optional inputs
    parser.add_argument(
        "--reads",
        help="Original reads file for additional statistics"
    )
    
    # Normalization options
    parser.add_argument(
        "-n", "--norm",
        type=int,
        help="Library size for normalization (RPM)"
    )
    
    # Analysis options
    parser.add_argument(
        "--run-patterns",
        action="store_true",
        help="Generate pattern analysis plots"
    )
    
    args = parser.parse_args()
    run_onlyplot(args)


if __name__ == "__main__":
    main()