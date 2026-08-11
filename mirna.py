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
from dataclasses import dataclass

from utils import (log, ensure_outdir, validate_file_exists, which, run_command,
                   ensure_within, remove_path, validate_prefix)


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


def require_bowtie_index(index_prefix):
    """The index miRDeep-P2 needs, or a loud refusal.

    Its `-x` took whatever it was given: omitted, the argument became the string
    "None", bowtie aligned nothing, every intermediate came out empty, and the
    pipeline still finished with `end...` and zero predictions. That is its
    fourth silent failure, and the only one that looked like a real run.
    """
    if not index_prefix:
        raise RuntimeError(
            "miRDeep-P2 needs --index-prefix. Without it bowtie aligns nothing "
            "and the pipeline reports success with zero predictions.")
    if not check_bowtie_index(index_prefix):
        raise RuntimeError(
            f"no bowtie index at {index_prefix} -- looked for .1.ebwt and its "
            f"siblings (or .ebwtl). miRDeep-P2 would align nothing and predict "
            f"nothing without failing. Build one with `bowtie-build`, or point "
            f"--index-prefix at the index the mapping module used.")
    return index_prefix


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


def work_paths(outdir, prefix):
    """Where the raw third-party output lands, same tree for every engine."""
    root = Path(outdir) / f"{prefix}.mirna_work"
    return {"root": root, "engine": root / "engine",
            "validation": root / "validation"}


def deliverable_paths(outdir, prefix):
    """Every file this module hands over, named off the prefix in one place."""
    outdir = Path(outdir)
    return {
        "table": outdir / f"{prefix}.mirna.tsv",
        "bed": outdir / f"{prefix}.mirna_loci.bed",
        "mature": outdir / f"{prefix}.mirna.mature.fa",
        "star": outdir / f"{prefix}.mirna.star.fa",
        "precursor": outdir / f"{prefix}.mirna.precursor.fa",
        "run_json": outdir / f"{prefix}.mirna.run.json",
        "structures": outdir / f"{prefix}.mirna_structures",
    }


def drop_validated_deliverables(paths, outdir):
    """Remove what only a validated run may leave behind, and say so.

    An unvalidated run writes no FASTA and no plots; leaving the previous
    run's would make the outdir claim a validated set this BED contradicts.
    These are deletes aimed by the prefix, so each one is confined to
    `outdir` before it happens.
    """
    for key in ("mature", "star", "precursor", "structures"):
        ensure_within(outdir, paths[key])
    for key in ("mature", "star", "precursor", "structures"):
        if paths[key].exists() or paths[key].is_symlink():
            remove_path(paths[key])
            log(f"removed stale {paths[key]}: this run did not validate", "WARN")


NA = "NA"

# A tool's own version flag. miRDeep-P2 ships none, so it answers NA rather
# than a number nobody printed.
_ENGINE_VERSION_CMD = {"mirdeep2": ["miRDeep2.pl", "--version"]}


def tool_version(argv):
    """What the tool says its version is, or NA -- never a guess."""
    if not shutil.which(argv[0]):
        return NA
    try:
        proc = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return NA
    m = re.search(r"\d+\.\d+(?:\.\d+)*", (proc.stdout or "") + (proc.stderr or ""))
    return m.group(0) if m else NA


def engine_version(engine):
    cmd = _ENGINE_VERSION_CMD.get(engine)
    return tool_version(cmd) if cmd else NA


def finish_mirna(predictions, args, work_dir, reads, mirdeep_rows=None):
    """Validate, write everything the module hands over, and say what happened."""
    import json
    from mirna_report import (build_records, write_table, write_bed,
                              write_bed_unvalidated, write_fasta,
                              collect_structures)

    kingdom = "plant" if args.engine == "mirdeep-p2" else "animal"
    result = None
    if getattr(args, "no_validate", False):
        log("--no-validate: writing predictions without miRScore", "WARN")
    elif getattr(args, "genome", None):
        result = validate_with_mirscore(
            predictions, args.genome, reads, kingdom,
            Path(work_dir) / "validation", threads=getattr(args, "threads", 1))
    else:
        log("no -g/--genome: skipping miRScore validation, which needs the "
            "genome to re-cut precursors with flanks", "WARN")

    records = build_records(
        predictions, args.engine,
        mirscore_rows=result.rows if result else None,
        submitted=result.submitted if result else None,
        mirdeep_rows=mirdeep_rows)

    paths = deliverable_paths(args.outdir, args.prefix)
    write_table(records, paths["table"])
    n_bed = write_bed(records, paths["bed"]) if result else write_bed_unvalidated(
        records, paths["bed"])
    if result:
        # Written whenever validation ran: an empty FASTA is a truthful answer,
        # and one left from a previous run would contradict this BED.
        for key, field in (("mature", "mature_seq"), ("star", "star_seq"),
                           ("precursor", "precursor_seq")):
            write_fasta(records, paths[key], field)
        n_plots = collect_structures(
            Path(work_dir) / "validation" / "mirscore", records,
            paths["structures"], args.outdir)
        log(f"{n_plots} structure plots in {paths['structures']}", "INFO")
    else:
        drop_validated_deliverables(paths, args.outdir)

    counts = {"predicted": len(records),
              "submitted": len(result.submitted) if result else "NA",
              "passed": sum(1 for r in records if r.result == "Pass"),
              "failed": sum(1 for r in records if r.result == "Fail"),
              "not_scored": sum(1 for r in records if r.result == "NotScored")}
    paths["run_json"].write_text(json.dumps(
        {"engine": args.engine, "engine_version": engine_version(args.engine),
         "kingdom": kingdom, "validated": result is not None,
         "mirscore_version": tool_version(["miRScore", "--version"]) if result else NA,
         "counts": counts}, indent=2))

    if result:
        log(f"{n_bed} validated loci written to {paths['bed']} "
            f"(of {len(records)} predicted)", "INFO")
    else:
        log(f"{n_bed} UNVALIDATED loci written to {paths['bed']}", "WARN")
    log(f"table: {paths['table']}", "INFO")
    return paths["bed"]


def novel_with_rows(results_dir):
    """The novel predictions and their raw row, from ONE result_*.csv.

    miRDeep2 names the file `result_DD_MM_YYYY_t_HH_MM_SS.csv`, so sorting by
    name orders by day of the month; and `provisional id` serials restart every
    run. Mixing two files puts one run's coordinates beside another's score.
    """
    for csv_path in sorted(Path(results_dir).glob("result_*.csv")):
        novel = parse_mirdeep2_novel(csv_path)
        if novel:
            return novel, {r["provisional id"]: r
                           for r in mirdeep2_novel_rows(csv_path)}
    return [], {}


def finish_animal_run(results_dir, args, work_root, reads):
    """Step 7, run whether or not there were predictions.

    Zero predictions is a result, and it is reported in the table and the
    run.json -- not by leaving the output directory empty.
    """
    novel, rows = novel_with_rows(results_dir)
    if novel:
        log(f"Step 7: validating {len(novel)} novel predictions...", "INFO")
    else:
        log("no novel predictions: writing an empty table, not nothing", "WARN")
    return finish_mirna(novel, args, work_root, reads, mirdeep_rows=rows)


def run_mirna_analysis(args):
    """Main function for miRNA analysis module."""
    
    if getattr(args, "engine", "mirdeep2") == "mirdeep-p2":
        outdir = ensure_outdir(Path(args.outdir))
        work = work_paths(outdir, args.prefix)
        preds = run_mirdeep_p2(args, args.reads, ensure_outdir(work["engine"]))
        finish_mirna(preds, args, work["root"], args.reads)
        log("miRNA analysis completed successfully!", "INFO")
        return

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
    work = work_paths(outdir, args.prefix)
    processed_dir = ensure_outdir(work["engine"] / "processed")
    results_dir = ensure_outdir(work["engine"] / "results")
    
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
    generate_mirna_report(results, work["root"], args.prefix)

    # Step 7: validate the novel predictions and emit the BED for --mirna-loci
    finish_animal_run(results_dir, args, work["root"], input_files["reads"])

    log("=" * 60, "INFO")
    log(f"miRNA analysis completed successfully!", "INFO")
    log(f"Results saved in: {Path(args.outdir)}", "INFO")
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
    # This entry point does not go through srnaprofiler.main(), where every
    # other prefix is checked, and the run below deletes paths named off it.
    try:
        validate_prefix(args.prefix)
    except ValueError as e:
        log(str(e), "ERROR")
        sys.exit(1)
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

    Absolute, not resolved: the reads are handed over as a symlink named so
    miRScore can key on it, and resolving it would restore the name we staged
    away from.
    """
    return Path(os.path.relpath(os.path.abspath(path), os.path.abspath(workdir)))


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


# The line miRDeep-P2 1.1.4 ships broken: it writes the collapsed reads to
# $out/$name/$name.ori.fa and then reads them from $out/$name.ori.fa.
_MIRDP2_PATH_BUG = 'seq="$results_folder/${filename}.ori.fa"'


def check_mirdeep_p2_install(pipeline):
    """Refuse a miRDeep-P2 install that would predict nothing and call it success.

    Three of its failure modes are silent -- the pipeline writes `end...` in its
    progress log and emits zero predictions -- so they are checked here, before
    a run, instead of being discovered in an empty result. Fixing them is
    `srnaprofiler setup`'s job: an analysis should not patch third-party code.
    """
    pipeline = Path(pipeline)
    if not pipeline.exists():
        raise RuntimeError(f"miRDeep-P2 not found at {pipeline}. Install it and "
                           f"prepare it with `srnaprofiler setup`.")

    if _MIRDP2_PATH_BUG in pipeline.read_text():
        raise RuntimeError(
            f"{pipeline.name} still carries the 1.1.4 path bug: it writes "
            f"$out/$name/$name.ori.fa and reads $out/$name.ori.fa, so it runs to "
            f"completion and predicts nothing. Run `srnaprofiler setup` to patch it.")

    index = pipeline.parent / "scripts" / "index"
    missing = [e for e in ("1", "2", "3", "4", "rev.1", "rev.2")
               if not (index / f"rfam_index.{e}.ebwt").exists()]
    if missing:
        raise RuntimeError(
            f"miRDeep-P2's rfam_index is missing from {index} (it ships only "
            f"mature_index). Without it the run predicts nothing. `srnaprofiler "
            f"setup` links it to the Rfam index built for --ncrna-db.")


def parse_mirdeep_p2_predictions(path):
    """miRDeep-P2's `*_filter_P_prediction` as a list of predictions.

    Columns: chrom, strand, read id, name, mature coords, precursor coords,
    mature sequence, precursor sequence. Coordinates are 1-based inclusive.
    A line that does not parse is skipped rather than guessed at -- a wrong
    coordinate here carves a locus in the wrong place and asserts miRNA on it.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            m = re.match(r"(\d+)\.\.(\d+)$", f[5].strip())
            if not m:
                continue
            out.append({"chrom": f[0], "strand": f[1], "name": f[3],
                        "start": int(m.group(1)), "end": int(m.group(2)),
                        "mature": f[6].strip().upper().replace("U", "T")})
    return out


def mirdeep2_novel_rows(path):
    """The novel section of miRDeep2's result CSV, rows as written."""
    import csv as _csv
    rows = list(_csv.reader(open(path), delimiter="\t"))
    start = None
    for i, r in enumerate(rows):
        if r and r[0] == "provisional id":
            hdr, start = r, i + 1
            break
    if start is None:
        return []
    out = []
    for r in rows[start:]:
        if not r or not r[0].strip():
            break
        out.append(dict(zip(hdr, r)))
    return out


def parse_mirdeep2_novel(path):
    """miRDeep2's novel predictions, in the shape the engines share."""
    out = []
    for d in mirdeep2_novel_rows(path):
        m = re.match(r"(\S+):(\d+)\.\.(\d+)(?::([+-]))?", d.get("precursor coordinate", ""))
        if not m:
            continue
        out.append({"chrom": m.group(1), "strand": m.group(4) or ".", "name": d["provisional id"],
                    "start": int(m.group(2)), "end": int(m.group(3)),
                    "mature": d["consensus mature sequence"].strip().upper().replace("U", "T")})
    return out


def mirscore_command(mature, hairpin, reads, kingdom, outdir, workdir, threads=1):
    """miRScore's invocation, with every path relative to `workdir`.

    It chdirs into its output directory and prefixes `../` to the reads path, so
    an absolute path becomes `..//abs/path`: bowtie aligns nothing, the BAM comes
    out empty, and the run dies later in samtools sort with an unrelated message.
    """
    rel = lambda p: str(relative_to(p, workdir))
    return ["miRScore", "-mature", rel(mature), "-hairpin", rel(hairpin),
            "-fastq", rel(reads), "-kingdom", kingdom,
            "-threads", str(threads), "-out", rel(outdir)]


def staged_reads_name(name):
    """A read file name miRScore can key on: no dots before the extension.

    It cuts the library name at the first dot but names the BAM with the whole
    stem, so `x.ncrna_filtered.fastq` looks up `x` against `x.ncrna_filtered_s`,
    finds nothing, and reports zero reads for every locus -- while the reads sit
    aligned in the BAM. Our own mapping module emits exactly that name.
    """
    stem, dot, ext = str(name).rpartition(".")
    if not dot:
        return str(name)
    return stem.replace(".", "_") + dot + ext


def stage_reads_for_mirscore(reads, workdir):
    """The reads under a name miRScore counts correctly, linked not copied."""
    reads = Path(reads).resolve()
    staged = Path(workdir) / staged_reads_name(reads.name)
    if staged.resolve() == reads:
        return staged
    if staged.is_symlink() or staged.exists():
        staged.unlink()
    staged.symlink_to(reads)
    return staged


MIRSCORE_REQUIRED = ("name", "mSeq", "msSeq", "msLen", "precSeq", "precLen",
                     "mStart", "mStop", "msStart", "msStop", "result", "flags",
                     "mismatches", "bulges", "mReads", "msReads", "totReads",
                     "precision")


def read_mirscore_rows(path):
    """miRScore's results keyed by name, with every column kept.

    A missing column raises: writing NA for a renamed column would put a wrong
    table in front of the user with no sign anything went wrong.
    """
    import csv as _csv
    with open(path) as fh:
        reader = _csv.DictReader(fh)
        missing = [c for c in MIRSCORE_REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise RuntimeError(
                f"{path} lacks the column(s) {', '.join(missing)}; miRScore's "
                f"output format changed and the table cannot be built from it")
        return {r["name"]: dict(r) for r in reader}


def parse_mirscore_results(path):
    """Names miRScore passed. Anything not an explicit pass is a fail."""
    return {n for n, r in read_mirscore_rows(path).items()
            if str(r.get("result", "")).strip().lower().startswith("pass")}


def no_reads_counted(reads_csv):
    """Every locus at zero total reads: impossible for a library that aligned."""
    import csv as _csv
    if not Path(reads_csv).exists():
        return False
    with open(reads_csv) as fh:
        rows = list(_csv.DictReader(fh))
    return bool(rows) and all(float(r.get("allReads") or 0) == 0 for r in rows)


def run_mirdeep_p2(args, reads, outdir):
    """The plant engine, with the three things it needs to not fail silently."""
    from utils import user_config
    recorded = user_config().get("mirdp2_pipeline")
    if not recorded:
        raise RuntimeError(
            "no miRDeep-P2 recorded: run `srnaprofiler setup --tools` to fetch "
            "and patch one, or `setup --mirdp2 PATH` to record an existing "
            "install")
    pipeline = Path(recorded)
    check_mirdeep_p2_install(pipeline)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)   # it does not create its own

    reads = Path(reads)
    if str(reads).endswith(".gz"):              # it would predict zero, silently
        plain = outdir / reads.name[:-3]
        log(f"miRDeep-P2 cannot read gzip; decompressing to {plain}", "INFO")
        with gzip.open(reads, "rb") as src, open(plain, "wb") as dst:
            shutil.copyfileobj(src, dst)
        reads = plain

    index = require_bowtie_index(getattr(args, "index_prefix", None))

    cmd = mirdeep_p2_command(pipeline, args.genome, index, reads,
                             outdir, threads=getattr(args, "threads", 1))
    log("Running miRDeep-P2...", "INFO")
    subprocess.run(cmd, check=True)

    hits = sorted(outdir.glob("*/*_filter_P_prediction"))
    if not hits:
        raise RuntimeError(
            f"miRDeep-P2 produced no prediction file under {outdir}. It reports "
            f"success even when it predicts nothing, so this is checked here.")
    preds = parse_mirdeep_p2_predictions(hits[0])
    log(f"miRDeep-P2: {len(preds)} predictions", "INFO")
    return preds


@dataclass
class Validation:
    """What miRScore said, and about which predictions it was asked."""
    rows: dict          # name -> its row in miRScore_results.csv
    submitted: set      # names handed to miRScore at all


def open_genome(path):
    """Seam for tests: the only place pysam.FastaFile is called."""
    import pysam
    return pysam.FastaFile(str(path))


def validate_with_mirscore(predictions, genome_path, reads, kingdom, workdir,
                           threads=1, flank=25):
    """miRScore's verdict per prediction, or None when it cannot be run.

    Takes the precursor from the genome by the engine's own coordinates, with
    flanks: the engines report it flush against the mature, and miRScore refuses
    a mature at the edge of the hairpin.
    """
    if not shutil.which("miRScore"):
        log("miRScore not found in PATH; skipping validation. Install it with "
            "`conda install -c bioconda miRScore`", "WARN")
        return None

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    mature_fa, hairpin_fa = workdir / "mature.fa", workdir / "hairpin.fa"
    genome = open_genome(genome_path)

    kept = set()
    with open(mature_fa, "w") as mf, open(hairpin_fa, "w") as hf:
        for p in predictions:
            hp = hairpin_with_flank(genome, p["chrom"], p["start"], p["end"],
                                    p["mature"], flank=flank)
            if hp is None:
                continue
            mf.write(f">{p['name']}\n{p['mature']}\n")
            hf.write(f">{p['name']}\n{hp}\n")
            kept.add(p["name"])
    if not kept:
        log("No prediction could be re-cut against the genome; skipping "
            "validation", "WARN")
        return None
    log(f"{len(kept)}/{len(predictions)} predictions prepared for miRScore", "INFO")

    out = workdir / "mirscore"
    if out.exists():
        shutil.rmtree(out)          # it refuses an existing output directory
    cmd = mirscore_command(mature_fa, hairpin_fa,
                           stage_reads_for_mirscore(reads, workdir),
                           kingdom, out, workdir, threads=threads)
    subprocess.run(cmd, cwd=workdir, check=True)

    results = out / "miRScore_results.csv"
    if not results.exists():
        raise RuntimeError(f"miRScore wrote no results to {results}")
    rows = read_mirscore_rows(results)
    passing = {n for n, r in rows.items()
               if str(r.get("result", "")).strip().lower().startswith("pass")}
    if not passing and no_reads_counted(out / "reads.csv"):
        raise RuntimeError(
            f"miRScore counted zero reads at every locus ({out/'reads.csv'}); "
            "that is its silent failure, not a biological result -- check that "
            "the reads it was given are the library the predictions came from")
    log(f"miRScore validated {len(passing)}/{len(kept)}", "INFO")
    return Validation(rows=rows, submitted=kept)
