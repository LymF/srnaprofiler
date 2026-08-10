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
from mirna import run_mirna_analysis
from setup_module import run_check, run_data, run_tools
from utils import (
    log, which, 
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

  # Generate plots directly from SAM (quick QC)
  %(prog)s plots --sam mapped.sam -p myproject -o results/

  # Generate plots per chromosome
  %(prog)s plots --sam mapped.sam -p myproject -o results/ --per-chrom

  # Re-generate plots from existing analysis
  %(prog)s plots -p myproject -o results/

  # Predict miRNAs with miRDeep2 (build new index)
  %(prog)s mirna -p myproject -o results/ --reads reads.fastq.gz -g genome.fa --species "My_Species"

  # Predict miRNAs using existing index from mapping
  %(prog)s mirna -p myproject -o results/ --reads reads.fastq.gz --index-prefix results/myproject.bowtie_index --species "My_Species"

For more information on each module:
  %(prog)s qc --help
  %(prog)s mapping --help
  %(prog)s clusters --help
  %(prog)s plots --help
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
        help="Verify and prepare an installation",
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
        action="store_true",
        help="Prepare reference data and record where it went")
    setup_parser.add_argument(
        "--rfam-dir",
        help="Where to build the Rfam structural RNA index (with --data). "
             "Downloads ~1.8 GB and runs bowtie-build; skipped if already built")
    setup_parser.add_argument(
        "--mirdp2",
        help="Path to miRDP2-v1.1.4_pipeline.bash (with --data): patches its "
             "path bug, links the Rfam index it needs, and records it")
    def _setup(a):
        if a.tools:
            sys.exit(run_tools(a))
        sys.exit(run_data(a) if a.data else run_check(a))

    setup_parser.set_defaults(func=_setup)

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
        help="Also create sorted BAM file with index"
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
    
    # ==================== Plots Module ====================
    plots_parser = subparsers.add_parser(
        "plots",
        help="Generate plots from SAM file or existing tables",
        description="""Generate sRNA distribution and pattern plots in multiple modes:
  
  Mode 1: From SAM file directly 
  Mode 2: From existing tables (re-generate plots after analysis)
  Mode 3: Per-chromosome mode (detailed analysis per chromosome/contig)

Examples:
  # Quick QC from SAM (combined plot)
  %(prog)s --sam data.sam -p sample1 -o results/
  
  # Per-chromosome analysis
  %(prog)s --sam data.sam -p sample1 -o results/ --per-chrom
  
  # Specific chromosomes only
  %(prog)s --sam data.sam -p sample1 -o results/ --per-chrom --chromosomes chr1 chr2
  
  # Re-generate from existing tables
  %(prog)s -p sample1 -o results/
  
  # Custom size range 
  %(prog)s --sam data.sam -p sample1 -o results/ --min-len 21 --max-len 24
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    plots_required = plots_parser.add_argument_group("required arguments")
    plots_required.add_argument(
        "-p", "--prefix",
        required=True,
        help="Output file prefix"
    )
    plots_required.add_argument(
        "-o", "--outdir",
        required=True,
        help="Output directory"
    )
    
    # SAM mode arguments (optional)
    plots_sam = plots_parser.add_argument_group("SAM mode arguments")
    plots_sam.add_argument(
        "--sam",
        help="SAM/BAM file (enables direct SAM processing mode)"
    )
    plots_sam.add_argument(
        "--per-chrom",
        action="store_true",
        help="Generate separate plots for each chromosome (requires --sam)"
    )
    plots_sam.add_argument(
        "--chromosomes",
        nargs="+",
        metavar="CHR",
        help="Specific chromosomes to plot (use with --per-chrom)"
    )
    plots_sam.add_argument(
        "--min-len",
        type=int,
        default=15,
        help="Minimum read length (default: %(default)s)"
    )
    plots_sam.add_argument(
        "--max-len",
        type=int,
        default=35,
        help="Maximum read length (default: %(default)s)"
    )
    
    # Table mode arguments
    plots_tables = plots_parser.add_argument_group("table mode arguments")
    plots_tables.add_argument(
        "--tables-dir",
        help="Directory containing data tables (default: <outdir>/tables)"
    )
    plots_tables.add_argument(
        "--categories",
        nargs="+",
        metavar="CAT",
        help="Specific categories to plot (e.g., putative_siRNA putative_piRNA)"
    )
    
    # Plot options
    plots_options = plots_parser.add_argument_group("plot options")
    plots_options.add_argument(
        "--no-base-distribution",
        action="store_false",
        dest="base_distribution",
        help="Skip base distribution plots"
    )
    plots_options.add_argument(
        "--no-patterns",
        action="store_false",
        dest="patterns",
        help="Skip pattern plots"
    )
    
    # Common options
    plots_common = plots_parser.add_argument_group("common options")
    plots_common.add_argument(
        "--r-scripts-path",
        help="Path to R scripts directory (auto-detect if not specified)"
    )
    
    plots_parser.set_defaults(func=run_plots)
    
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
    mirna_parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip miRScore validation. On by default: an unvalidated false "
             "prediction becomes a carved locus with miRNA asserted over it, "
             "while a missed one merely falls into the ordinary classes"
    )
    mirna_parser.add_argument(
        "--mirdp2",
        help="Path to miRDP2-v1.1.4_pipeline.bash. Only needed if "
             "`srnaprofiler setup` has not recorded it"
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
    cl = subparsers.add_parser("classify", help="Classify sRNA clusters (measure->decide engine)")
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


def run_plots(args):
    """
    Generate plots - delegates to plots.py module.
    
    This function is a wrapper that calls the plots module,
    which handles both SAM mode and table mode automatically.
    """
    try:
        from plots import run_plots as plots_run_plots, run_from_sam
        
        # Determine mode and delegate
        if hasattr(args, 'sam') and args.sam:
            # SAM mode
            log("Running plots in SAM mode", "INFO")
            run_from_sam(args)
        else:
            # Table mode
            log("Running plots in table mode", "INFO")
            plots_run_plots(args)
            
    except ImportError as e:
        log(f"Failed to import plots module: {e}", "ERROR")
        log("Please ensure plots.py is in the same directory", "ERROR")
        sys.exit(1)
    except Exception as e:
        log(f"Plot generation failed: {e}", "ERROR")
        import traceback
        traceback.print_exc()
        sys.exit(1)


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