#!/usr/bin/env python3
"""
mirna.py - miRNA Analysis Module for sRNA Profiler
Integrated miRDeep2 analysis for miRNA prediction and quantification
"""
import sys
import os
import argparse
import subprocess
import shutil
import gzip
import re
from pathlib import Path
from typing import Optional, Dict, List
from collections import defaultdict

from utils import log, ensure_outdir, validate_file_exists, which, run_command


def check_mirdeep2_dependencies():
    """Check if miRDeep2 and dependencies are available."""
    dependencies = ['bowtie-build', 'bowtie', 'fastq2fasta.pl', 'collapse_reads_md.pl', 
                   'mapper.pl', 'miRDeep2.pl', 'quantifier.pl']
    
    missing = []
    for dep in dependencies:
        if not which(dep):
            missing.append(dep)
    
    if missing:
        log(f"Missing miRDeep2 dependencies: {', '.join(missing)}", "ERROR")
        log("Please install miRDeep2 and ensure it's in PATH", "ERROR")
        return False
    
    log("All miRDeep2 dependencies found", "INFO")
    return True


def count_fasta_sequences(fasta_path: Path) -> int:
    """Count the number of sequences in a FASTA file."""
    try:
        count = 0
        with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    except Exception as e:
        log(f"Error counting sequences in {fasta_path}: {e}", "ERROR")
        return 0


def is_fastq_file(file_path: Path) -> bool:
    """Check if file is FASTQ format."""
    fastq_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
    return any(str(file_path).endswith(ext) for ext in fastq_extensions)


def is_fasta_file(file_path: Path) -> bool:
    """Check if file is FASTA format."""
    fasta_extensions = ['.fa', '.fasta', '.fa.gz', '.fasta.gz']
    return any(str(file_path).endswith(ext) for ext in fasta_extensions)


def validate_fasta_file(fasta_path: Path) -> bool:
    """Validate FASTA file format."""
    try:
        # Skip validation for ARF files
        if fasta_path.suffix == '.arf' or 'arf' in fasta_path.name:
            log(f"Skipping FASTA validation for ARF file: {fasta_path}", "DEBUG")
            return True
            
        with open(fasta_path, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                log(f"Invalid FASTA file: first line doesn't start with '>': {first_line}", "ERROR")
                return False
            
            # Check if file has at least one sequence
            seq_count = count_fasta_sequences(fasta_path)
            if seq_count == 0:
                log(f"FASTA file has no sequences: {fasta_path}", "ERROR")
                return False
                
        return True
    except Exception as e:
        log(f"Error validating FASTA file {fasta_path}: {e}", "ERROR")
        return False


READ_FORMATS = (".fastq", ".fq", ".fasta", ".fa", ".fna")


def reads_suffix(path: Path) -> str:
    """Extension to keep when copying reads: .suffix drops .fastq from .fastq.gz."""
    sufs = [s.lower() for s in path.suffixes]
    if sufs and sufs[-1] in (".gz", ".bz2"):
        if len(sufs) >= 2 and sufs[-2] in READ_FORMATS:
            return sufs[-2] + sufs[-1]
        return sufs[-1]
    return path.suffix


def copy_input_file(source_path: Path, dest_path: Path, file_type: str) -> Path:
    """Copy input file to pipeline directory structure if needed."""
    if source_path == dest_path or dest_path.exists():
        return dest_path
    
    log(f"Copying {file_type} to pipeline directory: {source_path} -> {dest_path}", "INFO")
    shutil.copy2(source_path, dest_path)
    return dest_path


def convert_fastq_to_fasta(fastq_path: Path, fasta_path: Path) -> bool:
    """Convert FASTQ to FASTA format using miRDeep2 script."""
    try:
        # Handle gzipped files by decompressing first
        if str(fastq_path).endswith('.gz'):
            log(f"Decompressing {fastq_path}...", "INFO")
            decompressed_fastq = fasta_path.parent / f"{fastq_path.stem}.fastq"
            with gzip.open(fastq_path, 'rb') as f_in:
                with open(decompressed_fastq, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            fastq_to_convert = decompressed_fastq
        else:
            fastq_to_convert = fastq_path
        
        cmd = ['fastq2fasta.pl', str(fastq_to_convert)]
        with open(fasta_path, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        
        # Clean up decompressed file if created
        if str(fastq_path).endswith('.gz') and 'decompressed_fastq' in locals() and decompressed_fastq.exists():
            decompressed_fastq.unlink()
        
        if result.returncode != 0:
            log(f"fastq2fasta.pl failed: {result.stderr}", "ERROR")
            return False
        
        # Verify the output file has sequences
        seq_count = count_fasta_sequences(fasta_path)
        if seq_count == 0:
            log(f"Converted FASTA file is empty: {fasta_path}", "ERROR")
            return False
        
        log(f"Converted FASTQ to FASTA: {fasta_path} ({seq_count} sequences)", "INFO")
        return True
        
    except Exception as e:
        log(f"Error converting FASTQ to FASTA: {e}", "ERROR")
        return False


def collapse_reads_mirdeep(fasta_path: Path, collapsed_path: Path, prefix: str = "lcg") -> bool:
    """Collapse identical reads using miRDeep2 collapse_reads_md.pl."""
    try:
        # Remove existing file to avoid conflicts
        if collapsed_path.exists():
            collapsed_path.unlink()
            
        cmd = ['collapse_reads_md.pl', str(fasta_path), prefix]
        with open(collapsed_path, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            log(f"collapse_reads_md.pl failed: {result.stderr}", "ERROR")
            return False
        
        # Verify the collapsed file
        collapsed_count = count_fasta_sequences(collapsed_path)
        if collapsed_count == 0:
            log(f"Collapsed FASTA file is empty: {collapsed_path}", "ERROR")
            return False
        
        log(f"Collapsed reads using miRDeep2: {collapsed_path} ({collapsed_count} sequences)", "INFO")
        return True
        
    except Exception as e:
        log(f"Error collapsing reads with miRDeep2: {e}", "ERROR")
        return False


def check_bowtie_index(index_prefix: str) -> bool:
    """Check if Bowtie index files exist."""
    index_path = Path(index_prefix)
    required_extensions = ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']
    
    for ext in required_extensions:
        if not (Path(str(index_path) + ext).exists() or 
                Path(str(index_path) + ext.replace('.ebwt', '.ebwtl')).exists()):
            return False
    return True


def build_bowtie_index(genome_path: Path, index_prefix: Path) -> bool:
    """Build Bowtie index for genome."""
    try:
        cmd = ['bowtie-build', str(genome_path), str(index_prefix)]
        run_command(cmd)
        log(f"Built Bowtie index: {index_prefix}", "INFO")
        return True
    except Exception as e:
        log(f"Error building Bowtie index: {e}", "ERROR")
        return False


def run_mirdeep2_mapper(reads_path: Path, index_prefix: Path, output_dir: Path, 
                       prefix: str, threads: int = 4) -> Optional[Dict[str, Path]]:
    """Run miRDeep2 mapper.pl to map reads to genome."""
    try:
        mapped_reads = output_dir / f"{prefix}.mapped.fa"
        arf_file = output_dir / f"{prefix}.genome.arf"
        
        # Remove existing files to avoid conflicts
        if mapped_reads.exists():
            mapped_reads.unlink()
        if arf_file.exists():
            arf_file.unlink()
        
        # First, verify the input file has sequences
        input_count = count_fasta_sequences(reads_path)
        if input_count == 0:
            log(f"Input FASTA file is empty: {reads_path}", "ERROR")
            return None
        
        log(f"Mapping {input_count} sequences to genome...", "INFO")
        
        cmd = [
            'mapper.pl',
            str(reads_path),
            '-c', '-i', '-j', '-q',  # Collapse, ignore unmapped, remove duplicates, quiet
            '-l', '18',              # Minimum length
            '-p', str(index_prefix),
            '-s', str(mapped_reads),
            '-t', str(arf_file),
            '-o', str(threads)
        ]
        
        run_command(cmd)
        
        # Verify mapping produced output
        if not mapped_reads.exists() or count_fasta_sequences(mapped_reads) == 0:
            log(f"Mapping produced no aligned reads: {mapped_reads}", "ERROR")
            return None
        
        mapped_count = count_fasta_sequences(mapped_reads)
        log(f"Reads mapped successfully: {mapped_count}/{input_count} sequences aligned", "INFO")
        log(f"  Mapped reads: {mapped_reads}", "INFO")
        log(f"  ARF file: {arf_file}", "INFO")
        
        return {
            'mapped_reads': mapped_reads,
            'arf_file': arf_file
        }
        
    except Exception as e:
        log(f"Error in miRDeep2 mapping: {e}", "ERROR")
        return None


def run_mirdeep2_prediction(mapped_reads: Path, genome: Path, arf_file: Path,
                           output_dir: Path, prefix: str, species: str,
                           lib_size: Optional[int] = None,
                           mature_ref: Optional[Path] = None,
                           hairpin_ref: Optional[Path] = None) -> bool:
    """Run miRDeep2.pl for miRNA prediction."""
    try:
        species_safe = species.replace(' ', '_')
        
        # Usa os.path.relpath para calcular caminhos relativos confiáveis
        def get_relpath(target_path: Path) -> str:
            return os.path.relpath(str(target_path), str(output_dir))
        
        cmd = [
            'miRDeep2.pl',
            get_relpath(mapped_reads),
            get_relpath(genome),
            get_relpath(arf_file),
            get_relpath(mature_ref) if mature_ref else 'none',
            'none',
            get_relpath(hairpin_ref) if hairpin_ref else 'none',
            '-t', species_safe
        ]
        
        if lib_size:
            cmd.extend(['-P', str(lib_size)])
        
        log(f"Running miRDeep2 from directory: {output_dir}", "DEBUG")
        log(f"Using os.path.relpath calculated paths:", "DEBUG")
        log(f"  reads: {get_relpath(mapped_reads)}", "DEBUG")
        log(f"  genome: {get_relpath(genome)}", "DEBUG")
        log(f"  arf: {get_relpath(arf_file)}", "DEBUG")
        
        # Execute from output directory
        result = subprocess.run(cmd, cwd=str(output_dir), check=True)
        
        log(f"miRDeep2 prediction completed for {species}", "INFO")
        return True
        
    except Exception as e:
        log(f"Error in miRDeep2 prediction: {e}", "ERROR")
        return False


def run_quantifier(mapped_reads: Path, output_dir: Path, species: str,
                  hairpin_ref: Optional[Path] = None,
                  mature_ref: Optional[Path] = None) -> bool:
    """Run quantifier.pl for miRNA expression quantification."""
    try:
        species_safe = species.replace(' ', '_')
        
        def get_relpath(target_path: Path) -> str:
            return os.path.relpath(str(target_path), str(output_dir))
        
        cmd = ['quantifier.pl', '-r', get_relpath(mapped_reads), '-t', species_safe]
        
        if hairpin_ref:
            cmd.extend(['-p', get_relpath(hairpin_ref)])
        if mature_ref:
            cmd.extend(['-m', get_relpath(mature_ref)])
        
        log(f"Running quantifier from directory: {output_dir}", "DEBUG")
        log(f"Using os.path.relpath for reads: {get_relpath(mapped_reads)}", "DEBUG")
        
        # Execute from output directory
        result = subprocess.run(cmd, cwd=str(output_dir), check=True)
        
        log("miRNA quantification completed", "INFO")
        return True
        
    except Exception as e:
        log(f"Error in miRNA quantification: {e}", "ERROR")
        return False


def parse_mirdeep2_results(output_dir: Path, prefix: str) -> Dict:
    """Parse miRDeep2 results and generate summary."""
    results = {
        'predicted_mirnas': 0,
        'novel_mirnas': 0,
        'known_mirnas': 0,
        'output_files': []
    }
    
    # Look for result files
    result_files = [
        'result_*.html',
        'miRNAs_expressed_all_samples_*.csv',
        'expression_analyses/*.csv'
    ]
    
    for pattern in result_files:
        for file_path in output_dir.glob(pattern):
            results['output_files'].append(str(file_path))
    
    # Try to parse the main result file
    result_csv = None
    for file_path in output_dir.glob('miRNAs_expressed_*.csv'):
        result_csv = file_path
        break
    
    if result_csv and result_csv.exists():
        try:
            import pandas as pd
            df = pd.read_csv(result_csv)
            results['predicted_mirnas'] = len(df)
            results['novel_mirnas'] = len(df[df['miRNA'].str.contains('novel', na=False)])
            results['known_mirnas'] = len(df) - results['novel_mirnas']
        except Exception as e:
            log(f"Could not parse result file: {e}", "WARN")
    
    return results


def generate_mirna_report(results: Dict, output_dir: Path, prefix: str):
    """Generate a summary report for miRNA analysis."""
    report_file = output_dir / f"{prefix}.mirna_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write(f"miRNA Analysis Report - {prefix}\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total miRNAs predicted: {results['predicted_mirnas']}\n")
        f.write(f"Novel miRNAs: {results['novel_mirnas']}\n")
        f.write(f"Known miRNAs: {results['known_mirnas']}\n\n")
        
        f.write("OUTPUT FILES:\n")
        f.write("-" * 30 + "\n")
        for file_path in results['output_files']:
            f.write(f"{Path(file_path).name}\n")
    
    log(f"miRNA report generated: {report_file}", "INFO")


def clean_mirna_directories(processed_dir: Path, results_dir: Path):
    """Clean miRNA analysis directories before starting."""
    try:
        # Clean processed files
        patterns = ['*.mapped.fa', '*.genome.arf']
        
        for pattern in patterns:
            for file_path in processed_dir.glob(pattern):
                file_path.unlink()
                log(f"Removed existing file: {file_path}", "INFO")
        
        # Clean result files
        result_patterns = ['result_*', 'miRNAs_*', 'expression_*', '*.csv', '*.html']
        for pattern in result_patterns:
            for file_path in results_dir.glob(pattern):
                if file_path.is_file():
                    file_path.unlink()
                    log(f"Removed existing result: {file_path}", "INFO")
        
    except Exception as e:
        log(f"Warning: Could not clean all directories: {e}", "WARN")


GENOME_SUFFIXES = (".fa", ".fasta", ".fna", ".fas")


def find_genome_beside_index(index_prefix: str):
    """Locate the genome FASTA sitting next to a bowtie index, or None.

    Covers .fna as well as .fa/.fasta: genomes downloaded from NCBI (the common
    case) use .fna, and missing it aborted the module with the file in plain
    sight. Compressed FASTA is ignored on purpose — miRDeep2 needs it
    uncompressed and would fail later, less clearly.

    With several FASTAs in the directory, the one whose name matches the index
    prefix wins; picking arbitrarily would silently analyse the wrong genome.
    """
    index_path = Path(index_prefix)
    candidates = [p for p in index_path.parent.iterdir()
                  if p.is_file() and p.suffix.lower() in GENOME_SUFFIXES]
    if not candidates:
        return None
    for p in candidates:
        if p.stem == index_path.name:
            return p
    return sorted(candidates)[0]


def setup_input_files(args, processed_dir: Path) -> Dict[str, Path]:
    """Copy all input files to pipeline directory structure for portability."""
    input_files = {}
    
    # 1. Handle genome file
    if hasattr(args, 'genome') and args.genome:
        genome_source = Path(args.genome)
        genome_local = processed_dir / f"{args.prefix}.genome.fa"
        input_files['genome'] = copy_input_file(genome_source, genome_local, "genome")
    elif hasattr(args, 'index_prefix') and args.index_prefix:
        # Try to find genome from index directory
        genome_source = find_genome_beside_index(args.index_prefix)
        if genome_source:
            log(f"Using genome found beside the index: {genome_source.name}", "INFO")
            genome_local = processed_dir / f"{args.prefix}.genome.fa"
            input_files['genome'] = copy_input_file(genome_source, genome_local, "genome")
        else:
            log(f"No genome FASTA ({'/'.join(GENOME_SUFFIXES)}) next to the index "
                f"at {Path(args.index_prefix).parent} — pass -g/--genome "
                f"explicitly", "ERROR")
            sys.exit(1)
    
    # 2. Handle reads file
    reads_source = Path(args.reads)
    reads_local = processed_dir / f"{args.prefix}.input_reads{reads_suffix(reads_source)}"
    input_files['reads'] = copy_input_file(reads_source, reads_local, "reads")
    
    # 3. Handle miRNA reference files
    if args.mature:
        mature_source = Path(args.mature)
        mature_local = processed_dir / f"{args.prefix}.mature.fa"
        input_files['mature'] = copy_input_file(mature_source, mature_local, "mature reference")
    
    if args.hairpin:
        hairpin_source = Path(args.hairpin)
        hairpin_local = processed_dir / f"{args.prefix}.hairpin.fa"
        input_files['hairpin'] = copy_input_file(hairpin_source, hairpin_local, "hairpin reference")
    
    return input_files


def run_mirna_analysis(args):
    """Main function for miRNA analysis module."""
    
    # Check dependencies
    if not check_mirdeep2_dependencies():
        sys.exit(1)
    
    # Validate input files exist
    validate_file_exists(Path(args.reads), "Reads file")
    
    if hasattr(args, 'genome') and args.genome:
        validate_file_exists(Path(args.genome), "Genome file")
    
    if args.mature:
        validate_file_exists(Path(args.mature), "Mature reference file")
    
    if args.hairpin:
        validate_file_exists(Path(args.hairpin), "Hairpin reference file")
    
    # Check if using existing index
    if hasattr(args, 'index_prefix') and args.index_prefix:
        index_files = [f"{args.index_prefix}.1.ebwt", f"{args.index_prefix}.2.ebwt"]
        if not any(Path(f).exists() for f in index_files):
            log(f"Bowtie index files not found for prefix: {args.index_prefix}", "ERROR")
            sys.exit(1)
    
    # Setup output directory structure
    outdir = ensure_outdir(Path(args.outdir))
    mirna_dir = ensure_outdir(outdir / "mirna")
    processed_dir = ensure_outdir(mirna_dir / "processed")
    results_dir = ensure_outdir(mirna_dir / "results")
    
    # Clean existing files to avoid conflicts
    clean_mirna_directories(processed_dir, results_dir)
    
    # Copy all input files to pipeline directory for portability
    log("Setting up input files for portable analysis...", "INFO")
    input_files = setup_input_files(args, processed_dir)
    
    log("=" * 60, "INFO")
    log(f"Starting miRNA analysis for {args.prefix}", "INFO")
    log("=" * 60, "INFO")
    
    log(f"Species: {args.species}", "INFO")
    log(f"Library size: {args.lib_size or 'Not specified'}", "INFO")
    log(f"Mature reference: {args.mature or 'None'}", "INFO")
    log(f"Hairpin reference: {args.hairpin or 'None'}", "INFO")
    
    # Step 1: Prepare reads
    log("Step 1: Preparing reads...", "INFO")
    
    reads_path = input_files['reads']
    
    # Handle input files - always convert to FASTA first
    if is_fastq_file(reads_path):
        # Convert FASTQ to FASTA
        fasta_reads = processed_dir / f"{args.prefix}.reads.fa"
        log(f"Converting FASTQ to FASTA: {reads_path} -> {fasta_reads}", "INFO")
        if not convert_fastq_to_fasta(reads_path, fasta_reads):
            sys.exit(1)
    elif is_fasta_file(reads_path):
        # Use FASTA file directly
        fasta_reads = reads_path
        log(f"Using existing FASTA file: {fasta_reads}", "INFO")
    else:
        log(f"Unsupported file format: {reads_path}. Expected FASTQ or FASTA.", "ERROR")
        sys.exit(1)
    
    # Verify FASTA file has sequences
    initial_count = count_fasta_sequences(fasta_reads)
    if initial_count == 0:
        log(f"Input FASTA file is empty: {fasta_reads}", "ERROR")
        sys.exit(1)
    
    log(f"Input contains {initial_count} sequences", "INFO")
    
    # Collapse reads using miRDeep2
    collapsed_reads = processed_dir / f"{args.prefix}.collapsed.fa"
    log(f"Collapsing reads with miRDeep2: {fasta_reads} -> {collapsed_reads}", "INFO")
    if not collapse_reads_mirdeep(fasta_reads, collapsed_reads, "lcg"):
        sys.exit(1)
    
    # Step 2: Prepare genome index
    log("Step 2: Preparing genome index...", "INFO")
    
    if hasattr(args, 'index_prefix') and args.index_prefix:
        # Use existing index from mapping module
        index_prefix = args.index_prefix
        log(f"Using existing Bowtie index: {index_prefix}", "INFO")
        
        if not check_bowtie_index(index_prefix):
            log(f"Bowtie index files not found for prefix: {index_prefix}", "ERROR")
            sys.exit(1)
    else:
        # Build new index
        index_prefix = str(processed_dir / f"{args.prefix}.genome_index")
        log(f"Building new Bowtie index: {index_prefix}", "INFO")
        if not build_bowtie_index(input_files['genome'], Path(index_prefix)):
            sys.exit(1)
    
    # Step 3: Map reads to genome
    log("Step 3: Mapping reads to genome...", "INFO")
    mapping_results = run_mirdeep2_mapper(
        collapsed_reads, Path(index_prefix), processed_dir, args.prefix, args.threads
    )
    
    if not mapping_results:
        sys.exit(1)
    
    # Step 4: Run miRNA prediction
    log("Step 4: Predicting miRNAs...", "INFO")
    mature_ref = input_files.get('mature')
    hairpin_ref = input_files.get('hairpin')
    
    # Use the local genome copy
    genome_path = input_files['genome']
    
    # Verify input files before prediction
    log("Verifying input files for miRDeep2...", "INFO")
    required_files = {
        'collapsed_reads': collapsed_reads,
        'genome': genome_path, 
        'arf_file': mapping_results['arf_file']
    }

    for file_type, file_path in required_files.items():
        if not file_path.exists():
            log(f"Required file not found: {file_path} ({file_type})", "ERROR")
            sys.exit(1)
        elif file_type == 'arf_file':
            # ARF files don't need FASTA validation
            log(f"Validated {file_type}: {file_path}", "DEBUG")
        elif not validate_fasta_file(file_path):
            log(f"Invalid FASTA file: {file_path}", "ERROR")
            sys.exit(1)
        else:
            log(f"Validated {file_type}: {file_path}", "DEBUG")
    
    if not run_mirdeep2_prediction(
        collapsed_reads,  # Use collapsed reads for prediction
        genome_path,
        mapping_results['arf_file'],
        results_dir,
        args.prefix,
        args.species,
        args.lib_size,
        mature_ref,
        hairpin_ref
    ):
        sys.exit(1)
    
    # Step 5: Quantify miRNA expression
    log("Step 5: Quantifying miRNA expression...", "INFO")
    if not run_quantifier(
        collapsed_reads,  # Use collapsed reads for quantification
        results_dir,
        args.species,
        hairpin_ref,
        mature_ref
    ):
        log("Quantification failed, but prediction completed", "WARN")
    
    # Step 6: Parse results and generate report
    log("Step 6: Generating results summary...", "INFO")
    results = parse_mirdeep2_results(results_dir, args.prefix)
    generate_mirna_report(results, mirna_dir, args.prefix)
    
    log("=" * 60, "INFO")
    log(f"miRNA analysis completed successfully!", "INFO")
    log(f"Results saved in: {mirna_dir}", "INFO")
    log("=" * 60, "INFO")


def main():
    """Standalone execution for miRNA module."""
    parser = argparse.ArgumentParser(
        description="miRNA analysis with miRDeep2"
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
        "--reads",
        required=True,
        help="Input reads file (FASTA/FASTQ)"
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Species name for annotation"
    )
    
    # Reference options (mutually exclusive)
    ref_group = parser.add_mutually_exclusive_group(required=True)
    ref_group.add_argument(
        "-g", "--genome",
        help="Reference genome FASTA (will build index)"
    )
    ref_group.add_argument(
        "--index-prefix",
        help="Pre-built Bowtie index prefix from mapping module"
    )
    
    # Optional reference files
    parser.add_argument(
        "--mature",
        help="Mature miRNA reference FASTA (from miRBase)"
    )
    parser.add_argument(
        "--hairpin", 
        help="Hairpin miRNA reference FASTA (from miRBase)"
    )
    
    # Analysis options
    parser.add_argument(
        "--lib-size",
        type=int,
        help="Library size for normalization"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    
    args = parser.parse_args()
    run_mirna_analysis(args)


if __name__ == "__main__":
    main()

# --- interop with the external miRNA tools -------------------------------
#
# Each function here exists because a third-party tool fails silently without
# it: it runs to completion, reports success, and produces nothing or garbage.

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def hairpin_with_flank(genome, chrom, start, end, mature, flank=25):
    """Precursor re-cut from the genome so the mature is not at either edge.

    miRDeep2 reports the precursor flush against the mature, and miRScore
    refuses that -- it needs context on both sides to judge the fold. The
    coordinates are miRDeep2's own (1-based, inclusive); only the flanks are
    added. Returns None when the mature is on neither strand of the region,
    which means the two do not describe the same locus and no hairpin should be
    invented for them.
    """
    lo = max(0, start - 1 - flank)
    hi = min(genome.get_reference_length(chrom), end + flank)
    seq = genome.fetch(chrom, lo, hi).upper()
    mature = mature.upper().replace("U", "T")
    for candidate in (seq, _revcomp(seq)):
        i = candidate.find(mature)
        if 0 < i and i + len(mature) < len(candidate):
            return candidate
    return None


def relative_to(path, workdir):
    """Path as miRScore needs it: relative to the directory it runs in.

    miRScore chdirs into its output directory and prefixes `../` to the reads
    path. An absolute path becomes `..//abs/path`, bowtie finds nothing, the BAM
    comes out empty, and the run only fails later, in samtools sort.
    """
    return Path(os.path.relpath(Path(path).resolve(), Path(workdir).resolve()))


def mirdeep_p2_command(pipeline, genome, index, reads, outdir, threads=1):
    """The miRDeep-P2 invocation, with the two inputs it silently mishandles.

    Gzipped reads are refused here rather than passed on: given a `.gz` the
    pipeline runs to completion, writes `end...` in its progress log and emits
    zero predictions. The output directory has to exist before it starts -- it
    does not create its own -- and that is the caller's job.
    """
    if str(reads).endswith(".gz"):
        raise ValueError(
            f"miRDeep-P2 cannot read gzipped input ({reads}): it would run to "
            f"completion and predict nothing. Decompress it first.")
    return ["bash", str(pipeline), "-g", str(genome), "-x", str(index),
            "-q", "-i", str(reads), "-p", str(threads), "-o", str(outdir)]
