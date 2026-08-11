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
from mirna import run_mirna_analysis
from setup_module import run_check, run_data, run_tools
from utils import (
    log, which, validate_prefix,
    DEFAULT_MIN_LEN, DEFAULT_QTRIM,
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
        prog="python3 srnaprofiler.py",
        description="Small RNA Profiler - Comprehensive small RNA analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
The order, and what each step is for:

  setup       once per machine: check the toolchain, fetch what conda cannot
  qc          trim reads and standardize the genome (needed before mirna)
  mapping     align to the genome, optionally dropping structural RNA
  mirna       find and validate miRNA loci, so classify can carve them out
  classify    call loci, measure their signatures, assign a class
  classplots  look at the result: per-class and per-locus figures

A typical run:

  %(prog)s setup --check
  %(prog)s qc       -p run -o out/ --reads1 reads.fastq.gz -g genome.fa
  %(prog)s mapping  -p run -o out/ --reads out/run.trimmed.R1.fastq.gz -g out/run.genome.std.fa
  %(prog)s mirna    -p run -o out/ --engine mirdeep-p2 -g out/run.genome.std.fa \\
                    --index-prefix out/run.bowtie_index --reads out/run.trimmed.R1.fastq.gz --species osa
  %(prog)s classify -p run -o out/ --sam out/run.bowtie.bam -g out/run.genome.std.fa \\
                    --mirna-loci out/run.mirna_loci.bed --config rules.plant.yaml
  %(prog)s classplots --clusters out/run.clusters.tsv --bam out/run.bowtie.bam -p run -o out/

Use --config rules.yaml for animals and rules.plant.yaml for plants; the classes
themselves are defined there. Every module takes --help.
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
        help="Trim reads and standardize the genome",
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
        help="Genome FASTA file to standardize (mapped as full genome, not tiled)"
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
    
    # ==================== Setup Module ====================
    setup_parser = subparsers.add_parser(
        "setup",
        help="Check the toolchain and prepare reference data",
        description="Check that a run will work here, and prepare the "
                    "reference data the pipeline needs")
    setup_parser.add_argument(
        "--check",
        action="store_true",
        help="Report what each module needs and whether it can run here")
    setup_parser.add_argument(
        "--tools",
        action="store_true",
        help="Install the toolchain into the ACTIVE conda environment, fetch "
             "miRDeep-P2 into ./tools, patch it, and verify the result")
    setup_parser.add_argument(
        "--data",
        metavar="DIR",
        help="Build the Rfam structural RNA index in DIR and record where it "
             "went (downloads ~1.8 GB and runs bowtie-build; skipped if DIR "
             "already holds it)")
    setup_parser.add_argument(
        "--mirdp2",
        metavar="PIPELINE",
        help="Record a miRDeep-P2 installed elsewhere instead of fetching one: "
             "patches its path bug and links the Rfam index it needs")
    def _setup(a):
        if a.tools:
            sys.exit(run_tools(a))
        sys.exit(run_data(a) if (a.data or a.mirdp2) else run_check(a))

    setup_parser.set_defaults(func=_setup)

    # ==================== Mapping Module ====================
    map_parser = subparsers.add_parser(
        "mapping",
        help="Align reads to the genome, with an optional ncRNA filter",
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
        "--ncrna-db",
        help="Bowtie1 index of structural RNA (build with "
             "extras/fetch_ncrna_db.sh); matching reads are dropped before the "
             "genome alignment"
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
        help="Accepted for compatibility: the sorted BAM is built by default"
    )
    map_output.add_argument(
        "--no-bam",
        action="store_true",
        help="Stop at the SAM and build no BAM (classify would have to sort it "
             "itself)"
    )
    map_output.add_argument(
        "--keep-sam",
        action="store_true",
        help="Keep the SAM after the BAM is written (default: remove it once "
             "the sorted BAM verifies; SAMs run tens of GB)"
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
    
    
    
    # ==================== miRNA Module ====================
    mirna_parser = subparsers.add_parser(
        "mirna",
        help="Predict and validate miRNA loci, for classify --mirna-loci",
        description="Predict miRNA loci with miRDeep2 or miRDeep-P2 and validate them with miRScore"
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
    mirna_parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip miRScore validation. On by default: an unvalidated false "
             "prediction becomes a carved locus with miRNA asserted over it, "
             "while a missed one merely falls into the ordinary classes"
    )
    mirna_required.add_argument(
        "--engine",
        required=True,
        choices=["mirdeep2", "mirdeep-p2"],
        help="Prediction engine: mirdeep2 for animals, mirdeep-p2 for plants "
             "(no default -- the wrong engine halves precision, see README)"
    )
    
    # Reference options (mutually exclusive group)
    # miRDeep2 needs both: the index to map, the FASTA to fold precursors.
    mirna_ref = mirna_parser.add_argument_group("reference")
    mirna_ref.add_argument(
        "-g", "--genome",
        required=True,
        help="Reference genome FASTA. Required: the engines fold precursors "
             "from it, and validation re-cuts them from it with flanks"
    )
    mirna_ref.add_argument(
        "--index-prefix",
        help="Pre-built Bowtie index prefix (e.g. from the mapping module). "
             "Without -g, the genome FASTA is looked up next to it"
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

    # ==================== Classify Module ====================
    from classify import run_classify
    cl = subparsers.add_parser("classify", help="Call loci and classify them by signature")
    cl.add_argument("-p", "--prefix", required=True)
    cl.add_argument("-o", "--outdir", required=True)
    cl.add_argument("--sam", required=True, help="SAM/BAM mapped to the whole genome")
    cl.add_argument("-g", "--genome",
                    help="Genome FASTA (indexed, .fai). Used as the null for the "
                         "5' base scores: a random cut in a locus yields a 5'U as "
                         "often as that locus carries T. Without it the null falls "
                         "back to a uniform 0.25 and s_1u/s_5a are inflated.")
    cl.add_argument("--mincov", type=float, default=1.0,
                    help="Minimum coverage to seed a locus, in READS PER MILLION "
                         "(default: %(default)s). Converted to an absolute depth "
                         "using the library size, so it transfers between runs.")
    cl.add_argument("--pad", type=int, default=None,
                    help="Extend each island by this many bp on both sides before "
                         "merging (default: from the config; 500 in both shipped "
                         "configs)")
    cl.add_argument("--mirna-loci", dest="mirna_loci", help="High-confidence miRNA BED")
    cl.add_argument("--organelle-contigs", dest="organelle_contigs",
                    help="Comma-separated contigs to leave out, e.g. chloroplast "
                         "and mitochondrion (NC_000932.1,NC_037304.1 in TAIR10.1)")
    cl.add_argument("--config", default="rules.yaml",
                    help="Scoring config (default: %(default)s, calibrated for "
                         "arthropods). Use rules.plant.yaml for plants.")
    cl.add_argument("--norm", type=int, help="Library size for RPM")
    cl.add_argument("--threads", type=int, default=4,
                    help="Threads for BAM sorting and per-chromosome processing (default: %(default)s)")
    cl.set_defaults(func=run_classify)

    from classplots import run_classplots
    cp = subparsers.add_parser("classplots",
        help="Per-class plots from a classify .clusters.tsv + BAM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  # all classes, 10 top loci each\n"
            "  srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \\\n"
            "      -o out --top-n 10\n\n"
            "  # only piRNA and siRNA (class names are case-insensitive)\n"
            "  srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \\\n"
            "      -o out --classes piRNA,siRNA\n\n"
            "  # one explicit locus, regardless of class\n"
            "  srnaprofiler classplots --clusters run.clusters.tsv --bam sample.sorted.bam \\\n"
            "      -o out --region chr2:100-4000\n"
        ))
    cp.add_argument("--clusters", required=True,
                    help="classify .clusters.tsv (e.g. run.clusters.tsv)")
    cp.add_argument("--bam", required=True,
                    help="Sorted+indexed BAM used by classify (e.g. sample.sorted.bam)")
    cp.add_argument("-o", "--outdir", required=True)
    cp.add_argument("-p", "--prefix", default="classplots")
    cp.add_argument("--top-n", dest="top_n", type=int, default=5,
                    help="Per-locus plots for the N most abundant loci per class (0 disables)")
    cp.add_argument("--max-loci", dest="max_loci", type=int, default=500,
                    help="Cap the loci aggregated into the per-class plots, taking "
                         "the most abundant (default: %(default)s, 0 = no cap). "
                         "Without it a large class holds most of the BAM in memory: "
                         "the 1000 densest hc-siRNA loci of an Arabidopsis run are "
                         "10.7M alignments. The log reports how much of the class "
                         "read mass the plot covers.")
    cp.add_argument("--classes",
                    help="Comma-separated subset of classes to plot, case-insensitive "
                         "(e.g. 'piRNA,siRNA' or 'pirna'); default: all classes in the TSV")
    cp.add_argument("--region", action="append", metavar="CHROM:START-END",
                    help="Plot one explicit coordinate span, e.g. chr2:100-4000 "
                         "(repeatable), regardless of class")
    cp.add_argument("--min-length", dest="min_length", type=int, default=15,
                    help="Lower bound of the logo size window (size-dist/pattern plots use "
                         "the fixed 15-35 publication window of their R scripts)")
    cp.add_argument("--max-length", dest="max_length", type=int, default=35,
                    help="Upper bound of the logo size window (see --min-length)")
    cp.add_argument("--logo-trim", dest="logo_trim", type=int, default=None,
                    help="Trim reads to N nt for the logo (default: class dominant DicerCall)")
    cp.add_argument("--logo-format", dest="logo_format", default="pdf",
                    choices=["pdf", "svg", "png", "eps"],
                    help="Sequence-logo output format (default: pdf). Rendered by logomaker "
                         "via matplotlib (no Ghostscript needed)")
    cp.add_argument("--no-signals", dest="no_signals", action="store_true",
                    help="Skip R/WebLogo signal plots (QC summary only)")
    cp.add_argument("--no-summary", dest="no_summary", action="store_true",
                    help="Skip the matplotlib QC summary")
    cp.set_defaults(func=run_classplots)

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
        if args.subcmd in ['clusters', 'plots']:
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
        # setup has neither: it inspects the installation, it produces no run
        if hasattr(args, "outdir"):
            log(f"Output directory: {args.outdir}", "INFO")
        if hasattr(args, "prefix"):
            # One seam for every subcommand: they all name files off it.
            validate_prefix(args.prefix)
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