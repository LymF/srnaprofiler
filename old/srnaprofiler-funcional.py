#!/usr/bin/env python3
"""
sRNA Profiler - Small RNA Analysis Pipeline
Main entry point for the analysis pipeline
"""
import argparse
import sys
from pathlib import Path

# Import the analysis modules
from qc import run_qc
from mapping import run_mapping
from clusters import run_clusters
from onlyplot import run_onlyplot
from mirna import run_mirna_analysis
from utils import (
    log, which, 
    DEFAULT_KMER_SIZE, DEFAULT_MIN_LEN, DEFAULT_QTRIM,
    DEFAULT_THREADS, DEFAULT_MISMATCHES, DEFAULT_MAX_ALIGNMENTS,
    DEFAULT_SI_THRESHOLD, DEFAULT_PI_THRESHOLD,
    DEFAULT_MAX_LEN, FASTQ_ADAPTERS
)


def print_header():
    """Version 5: Clean and professional with DNA helix symbol"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    BG_BLUE = '\033[104m'  # Light blue background
    
    header = f"""
{BOLD}{BLUE}                    ╔══════════════════════════════╗
                    ║  {WHITE}sRNA{GREEN}Profiler{BLUE}  │  v1.0       ║
                    ║                              ║
                    ╚══════════════════════════════╝{RESET}

                {WHITE}🧬  Small RNA Sequencing Analysis Pipeline  🧬{RESET}

{BLUE}    ┌─────────────────────────────────────────────────────────┐
    │ {WHITE}• Quality Control & Read Trimming                       {BLUE}│
    │ {WHITE}• Genome Preparation & Tiling                           {BLUE}│
    │ {WHITE}• Read Alignment with Bowtie                            {BLUE}│
    │ {WHITE}• Cluster Identification & Classification               {BLUE}│
    │ {WHITE}• miRNA Prediction & Quantification                     {BLUE}│
    │ {WHITE}• Visualization & Reporting                             {BLUE}│
    └─────────────────────────────────────────────────────────┘{RESET}
"""
    print(header)


def build_parser():
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="srnaprofiler",
        description="Small RNA Profiler - Comprehensive small RNA analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run QC on reads and prepare genome
  %(prog)s qc -p myproject -o results/ --reads1 sample.fastq.gz -g genome.fa

  # Map trimmed reads to tiled genome  
  %(prog)s mapping -p myproject -o results/ --reads trimmed.fastq.gz -g genome_tiles.fa

  # Analyze clusters with miRNA filtering
  %(prog)s clusters -p myproject -o results/ --sam mapped.sam -g genome_tiles.fa --mirna-fasta mirnas.fa

  # Full pipeline with pattern analysis
  %(prog)s clusters -p myproject -o results/ --sam mapped.sam -g genome.fa --run-patterns

  # Generate plots from existing files
  %(prog)s onlyplot -p myproject -o results/ --sam mapped.sam -g genome.fa --reads reads.fastq.gz

  # Predict miRNAs with miRDeep2 (build new index)
  %(prog)s mirna -p myproject -o results/ --reads reads.fastq.gz -g genome.fa --species "My_Species"

  # Predict miRNAs using existing index from mapping
  %(prog)s mirna -p myproject -o results/ --reads reads.fastq.gz --index-prefix results/myproject.bowtie_index --species "My_Species"

For more information on each module:
  %(prog)s qc --help
  %(prog)s mapping --help
  %(prog)s clusters --help
  %(prog)s onlyplot --help
  %(prog)s mirna --help
        """
    )
    
    parser.add_argument(
        "--version", 
        action="version", 
        version="%(prog)s 1.0"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    subparsers = parser.add_subparsers(
        dest="subcmd",
        required=True,
        help="Analysis modules"
    )
    
    # ==================== QC Module ====================
    qc_parser = subparsers.add_parser(
        "qc",
        help="Quality control, trimming, and genome preparation",
        description="Perform quality control on reads and prepare genome for mapping"
    )
    
    # Required arguments
    qc_required = qc_parser.add_argument_group("required arguments")
    qc_required.add_argument(
        "-p", "--prefix",
        required=True,
        help="Output file prefix for all generated files"
    )
    qc_required.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory for results"
    )
    
    # Input files
    qc_input = qc_parser.add_argument_group("input files")
    qc_input.add_argument(
        "--reads1",
        help="Input reads file (FASTQ/FASTA, can be gzipped)"
    )
    qc_input.add_argument(
        "--reads2",
        help="Paired-end mate file (FASTQ/FASTA, can be gzipped)"
    )
    qc_input.add_argument(
        "-g", "--genome",
        help="Genome FASTA file to standardize and tile"
    )
    
    # Trimming options
    qc_trim = qc_parser.add_argument_group("trimming options")
    qc_trim.add_argument(
        "--adapter",
        help=f"Adapter sequence (default: auto-detect or {FASTQ_ADAPTERS['illumina_srna']})"
    )
    qc_trim.add_argument(
        "--tool",
        choices=["auto", "fastp", "cutadapt"],
        default="auto",
        help="Trimming tool to use (default: %(default)s)"
    )
    qc_trim.add_argument(
        "--phred",
        type=int,
        choices=[33, 64],
        default=33,
        help="Phred quality encoding (default: %(default)s)"
    )
    qc_trim.add_argument(
        "--qtrim",
        type=int,
        default=DEFAULT_QTRIM,
        help="Quality trimming threshold (default: %(default)s)"
    )
    qc_trim.add_argument(
        "--min-len",
        type=int,
        default=DEFAULT_MIN_LEN,
        help="Minimum read length after trimming (default: %(default)s). For miRNA analyses --min-len 18 is optimal"
    )
    
    # Genome processing
    qc_genome = qc_parser.add_argument_group("genome processing")
    qc_genome.add_argument(
        "--kmer",
        type=int,
        default=DEFAULT_KMER_SIZE,
        help="Tile size for genome splitting in bp (default: %(default)s)"
    )
    qc_genome.add_argument(
        "--no-tiling",
        action="store_true",
        help="Disable genome tiling (use for miRNA analysis)"
    )
    
    # Performance options
    qc_perf = qc_parser.add_argument_group("performance options")
    qc_perf.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Number of threads to use (default: %(default)s)"
    )
    qc_perf.add_argument(
        "--max-mem-mb",
        type=int,
        help="Maximum memory in MB (optional)"
    )
    
    qc_parser.set_defaults(func=run_qc)
    
    # ==================== Mapping Module ====================
    map_parser = subparsers.add_parser(
        "mapping",
        help="Align reads using Bowtie1",
        description="Map small RNA reads to reference genome or transcriptome"
    )
    
    # Required arguments
    map_required = map_parser.add_argument_group("required arguments")
    map_required.add_argument(
        "-p", "--prefix",
        required=True,
        help="Output file prefix"
    )
    map_required.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    
    # Input reads (mutually exclusive)
    map_reads = map_parser.add_argument_group("input reads (choose one)")
    reads_group = map_reads.add_mutually_exclusive_group(required=True)
    reads_group.add_argument(
        "--reads",
        help="Single-end reads (FASTQ/FASTA, can be gzipped)"
    )
    reads_group.add_argument(
        "--reads-paired",
        nargs=2,
        metavar=("R1", "R2"),
        help="Paired-end reads"
    )
    
    # Reference (mutually exclusive)
    map_ref = map_parser.add_argument_group("reference (choose one)")
    ref_group = map_ref.add_mutually_exclusive_group(required=True)
    ref_group.add_argument(
        "-g", "--genome",
        help="Reference genome FASTA (will build index)"
    )
    ref_group.add_argument(
        "--index-prefix",
        help="Pre-built Bowtie1 index prefix"
    )
    
    # Alignment options
    map_align = map_parser.add_argument_group("alignment options")
    map_align.add_argument(
        "--fmt",
        choices=["auto", "fasta", "fastq"],
        default="auto",
        help="Input format (default: %(default)s)"
    )
    map_align.add_argument(
        "--mismatches",
        type=int,
        default=DEFAULT_MISMATCHES,
        help="Maximum mismatches allowed (default: %(default)s)"
    )
    map_align.add_argument(
        "--k",
        type=int,
        default=DEFAULT_MAX_ALIGNMENTS,
        help="Maximum alignments to report per read (default: %(default)s)"
    )
    map_align.add_argument(
        "--keep-unmapped",
        action="store_true",
        help="Keep unmapped reads in output SAM file"
    )
    
    # Output options
    map_output = map_parser.add_argument_group("output options")
    map_output.add_argument(
        "--sam-out",
        help="Custom SAM output path (default: <outdir>/<prefix>.bowtie.sam)"
    )
    map_output.add_argument(
        "--create-bam",
        action="store_true",
        help="Also create sorted BAM file with index"
    )
    
    # Performance options
    map_perf = map_parser.add_argument_group("performance options")
    map_perf.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Number of threads (default: %(default)s)"
    )
    
    map_parser.set_defaults(func=run_mapping)
    
    # ==================== Clusters Module ====================
    cluster_parser = subparsers.add_parser(
        "clusters",
        help="Analyze small RNA clusters and patterns",
        description="Classify contigs, filter miRNAs, and analyze small RNA patterns"
    )
    
    # Required arguments
    cluster_required = cluster_parser.add_argument_group("required arguments")
    cluster_required.add_argument(
        "-p", "--prefix",
        required=True,
        help="Output file prefix"
    )
    cluster_required.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    cluster_required.add_argument(
        "--sam",
        required=True,
        help="SAM/BAM file with mapped reads (can be gzipped)"
    )
    cluster_required.add_argument(
        "-g", "--genome",
        required=True,
        help="Reference genome (tiled) FASTA"
    )
    
    # Optional inputs
    cluster_input = cluster_parser.add_argument_group("optional inputs")
    cluster_input.add_argument(
        "--reads",
        help="Original reads file for additional statistics"
    )
    
    # Classification thresholds
    cluster_class = cluster_parser.add_argument_group("classification thresholds")
    cluster_class.add_argument(
        "--si",
        type=int,
        default=DEFAULT_SI_THRESHOLD,
        help="Minimum reads for siRNA classification (default: %(default)s)"
    )
    cluster_class.add_argument(
        "--pi",
        type=int,
        default=DEFAULT_PI_THRESHOLD,
        help="Minimum reads for piRNA classification (default: %(default)s)"
    )
    cluster_class.add_argument(
        "--len-min",
        type=int,
        default=DEFAULT_MIN_LEN,
        help="Minimum small RNA length (default: %(default)s)"
    )
    cluster_class.add_argument(
        "--len-max",
        type=int,
        default=DEFAULT_MAX_LEN,
        help="Maximum small RNA length (default: %(default)s)"
    )
    
    # Normalization options
    cluster_norm = cluster_parser.add_argument_group("normalization options")
    cluster_norm.add_argument(
        "-n", "--norm",
        type=int,
        help="Library size for normalization (RPM)"
    )
    
    # miRNA filtering options
    cluster_mirna = cluster_parser.add_argument_group("miRNA filtering")
    mirna_ref_group = cluster_mirna.add_mutually_exclusive_group()
    mirna_ref_group.add_argument(
        "--mirna-fasta",
        help="miRNA FASTA reference (headers with coordinates accepted)"
    )
    mirna_ref_group.add_argument(
        "--mirna-gff",
        help="miRNA GFF reference (coordinates)"
    )
    mirna_ref_group.add_argument(
        "--mirna-bed",
        help="miRNA BED reference (coordinates)"
    )
    cluster_mirna.add_argument(
        "--mirna-keep",
        action="store_true",
        help="Keep miRNA clusters in downstream analysis (label only)"
    )
    
    # Analysis options
    cluster_analysis = cluster_parser.add_argument_group("analysis options")
    cluster_analysis.add_argument(
        "--run-patterns",
        action="store_true",
        help="Analyze nucleotide patterns and ping-pong signatures"
    )
    cluster_analysis.add_argument(
        "--force",
        action="store_true",
        help="Force re-run of all steps"
    )
    
    # Alignment filtering options
    cluster_align = cluster_parser.add_argument_group("alignment filtering")
    cluster_align.add_argument(
        "--mapq",
        type=int,
        default=0,
        help="Minimum MAPQ quality score (default: %(default)s)"
    )
    cluster_align.add_argument(
        "--nm",
        type=int,
        default=1_000_000_000,
        help="Maximum number of mismatches (default: %(default)s)"
    )
    
    cluster_parser.set_defaults(func=run_clusters)
    
    # ==================== OnlyPlot Module ====================
    onlyplot_parser = subparsers.add_parser(
        "onlyplot",
        help="Generate plots from existing files without re-running analysis",
        description="Create plots from SAM, genome, and reads files without re-analyzing"
    )
    
    # Required arguments
    onlyplot_required = onlyplot_parser.add_argument_group("required arguments")
    onlyplot_required.add_argument(
        "-p", "--prefix",
        required=True,
        help="Output file prefix"
    )
    onlyplot_required.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    onlyplot_required.add_argument(
        "--sam",
        required=True,
        help="SAM/BAM file with mapped reads (can be gzipped)"
    )
    onlyplot_required.add_argument(
        "-g", "--genome",
        required=True,
        help="Reference genome FASTA"
    )
    
    # Optional inputs
    onlyplot_input = onlyplot_parser.add_argument_group("optional inputs")
    onlyplot_input.add_argument(
        "--reads",
        help="Original reads file for additional statistics"
    )
    
    # Normalization options
    onlyplot_norm = onlyplot_parser.add_argument_group("normalization options")
    onlyplot_norm.add_argument(
        "-n", "--norm",
        type=int,
        help="Library size for normalization (RPM)"
    )
    
    # Analysis options
    onlyplot_analysis = onlyplot_parser.add_argument_group("analysis options")
    onlyplot_analysis.add_argument(
        "--run-patterns",
        action="store_true",
        help="Generate pattern analysis plots"
    )
    
    onlyplot_parser.set_defaults(func=run_onlyplot)
    
    # ==================== miRNA Module ====================
    mirna_parser = subparsers.add_parser(
        "mirna",
        help="miRNA analysis with miRDeep2",
        description="Predict and quantify miRNAs using miRDeep2"
    )
    
    # Required arguments
    mirna_required = mirna_parser.add_argument_group("required arguments")
    mirna_required.add_argument(
        "-p", "--prefix",
        required=True,
        help="Output file prefix"
    )
    mirna_required.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    mirna_required.add_argument(
        "--reads",
        required=True,
        help="Input reads file (FASTA/FASTQ, can be gzipped)"
    )
    mirna_required.add_argument(
        "--species",
        required=True,
        help="Species name for miRNA annotation"
    )
    
    # Reference options (mutually exclusive group)
    mirna_ref = mirna_parser.add_argument_group("reference (choose one)")
    ref_group = mirna_ref.add_mutually_exclusive_group(required=True)
    ref_group.add_argument(
        "-g", "--genome",
        help="Reference genome FASTA (will build index)"
    )
    ref_group.add_argument(
        "--index-prefix",
        help="Pre-built Bowtie index prefix from mapping module"
    )
    
    # Reference files (optional)
    mirna_ref_files = mirna_parser.add_argument_group("reference files (optional)")
    mirna_ref_files.add_argument(
        "--mature",
        help="Mature miRNA reference FASTA (from miRBase)"
    )
    mirna_ref_files.add_argument(
        "--hairpin",
        help="Hairpin/precursor miRNA reference FASTA (from miRBase)"
    )
    
    # Analysis options
    mirna_analysis = mirna_parser.add_argument_group("analysis options")
    mirna_analysis.add_argument(
        "--lib-size",
        type=int,
        help="Library size for normalization (RPM)"
    )
    mirna_analysis.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads (default: %(default)s)"
    )
    
    mirna_parser.set_defaults(func=run_mirna_analysis)
    
    return parser


def check_dependencies():
    """Check for required dependencies."""
    log("Checking dependencies...", "INFO")
    
    required = {
        'qc': ['fastp or cutadapt'],  # fastqc is optional
        'mapping': ['bowtie', 'bowtie-build', 'samtools'],
        'clusters': ['python packages: pandas, numpy, matplotlib, biopython, pybedtools'],
        'onlyplot': ['python packages: pandas, numpy, matplotlib, biopython, pybedtools'],
        'mirna': ['bowtie-build', 'bowtie', 'fastq2fasta.pl', 'collapse_reads_md.pl', 
                 'mapper.pl', 'miRDeep2.pl', 'quantifier.pl']
    }
    
    missing = []
    warnings = []
    
    # Check QC tools
    if not which('fastp') and not which('cutadapt'):
        missing.append('fastp or cutadapt (required for trimming)')
    if not which('fastqc'):
        warnings.append('fastqc (optional for QC reports)')
    
    # Check mapping tools
    if not which('bowtie'):
        missing.append('bowtie (required for mapping)')
    if not which('bowtie-build'):
        missing.append('bowtie-build (required for index building)')
    if not which('samtools'):
        missing.append('samtools (required for SAM processing)')
    
    # Check Python packages for clusters and onlyplot
    py_packages = ['pandas', 'numpy', 'Bio', 'pybedtools']
    for pkg in py_packages:
        try:
            if pkg == 'Bio':
                __import__('Bio.SeqIO')
            elif pkg == 'pybedtools':
                __import__('pybedtools')
            else:
                __import__(pkg)
        except ImportError:
            missing.append(f'{pkg} (Python package for clusters/onlyplot)')
    
    # Check miRDeep2 dependencies (only for mirna module)
    mirna_deps = ['fastq2fasta.pl', 'collapse_reads_md.pl', 'mapper.pl', 'miRDeep2.pl', 'quantifier.pl']
    for dep in mirna_deps:
        if not which(dep):
            missing.append(f'{dep} (required for miRNA analysis)')
    
    # Report results
    if missing:
        log("Missing required tools:", "ERROR")
        for tool in missing:
            log(f"  - {tool}", "ERROR")
    
    if warnings:
        log("Missing optional tools:", "WARN")
        for tool in warnings:
            log(f"  - {tool}", "WARN")
    
    if not missing and warnings:
        log("All required dependencies found, but some optional tools are missing", "INFO")
        log("Some features may have limited functionality", "INFO")
    elif not missing:
        log("All dependencies found", "INFO")
    
    return len(missing) == 0


def main():
    """Main entry point."""
    # Print header
    print_header()
    
    # Parse arguments
    parser = build_parser()
    args = parser.parse_args()
    
    # Set up logging verbosity
    if hasattr(args, 'verbose') and args.verbose:
        import utils
        utils.LOG_LEVEL = "DEBUG"
    
    # Check dependencies
    all_deps_available = check_dependencies()
    
    if not all_deps_available and hasattr(args, 'subcmd'):
        if args.subcmd in ['clusters', 'onlyplot']:
            log("Warning: Some dependencies for cluster/plot analysis are missing", "WARN")
            log("Analysis may fail or have limited functionality", "WARN")
        elif args.subcmd == 'mirna':
            log("Error: Required dependencies for miRNA analysis are missing", "ERROR")
            log("miRNA analysis cannot proceed", "ERROR")
            sys.exit(1)
        elif args.subcmd in ['qc', 'mapping']:
            log("Error: Required dependencies are missing", "ERROR")
            log("Analysis cannot proceed", "ERROR")
            sys.exit(1)
    
    # Run the selected module
    try:
        log(f"Starting {args.subcmd} module", "INFO")
        log(f"Output directory: {args.outdir}", "INFO")
        log(f"Project prefix: {args.prefix}", "INFO")
        
        args.func(args)
        log(f"{args.subcmd} module completed successfully", "INFO")
        
    except KeyboardInterrupt:
        log("Pipeline interrupted by user", "WARN")
        sys.exit(130)
    
    except Exception as e:
        log(f"Pipeline failed: {e}", "ERROR")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()