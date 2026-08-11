#!/usr/bin/env python3
"""
Viral sRNA Assembly Pipeline

Steps:
1. Trimming with cutadapt (automatic adapter detection)
2. Mapping against reference genomes (host, bacterial RefSeq, transposons)
3. Separation of unmapped reads by size (20-23nt and 24-30nt)
4. Assembly with metaSPAdes (k-mer 15) for each size range
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import logging
from typing import List, Tuple

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ViralSRNAPipeline:
    """Pipeline for viral sRNA analysis"""
    
    def __init__(self, output_dir, host_genome, bacterial_refseq=None, 
                 transposon_db=None, threads=4, viralquest_params=None):
        """
        Initialize the pipeline
        
        Args:
            output_dir: Output directory
            host_genome: Host reference genome
            bacterial_refseq: Bacterial RefSeq (optional)
            transposon_db: Transposon database (optional)
            threads: Number of threads
            viralquest_params: Dictionary with ViralQuest parameters (optional)
        """
        self.output_dir = Path(output_dir)
        self.host_genome = Path(host_genome)
        self.bacterial_refseq = Path(bacterial_refseq) if bacterial_refseq else None
        self.transposon_db = Path(transposon_db) if transposon_db else None
        self.threads = threads
        self.viralquest_params = viralquest_params or {}
        
        # Create main output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def setup_sample_directories(self, sample_name):
        """
        Create directory structure for a sample
        
        Args:
            sample_name: Name of the sample
        """
        sample_dir = self.output_dir / sample_name
        sample_dir.mkdir(exist_ok=True)
        
        trimmed_dir = sample_dir / "01_trimmed"
        mapping_dir = sample_dir / "02_mapping"
        filtered_dir = sample_dir / "03_size_filtered"
        assembly_dir = sample_dir / "04_assembly"
        viralquest_dir = sample_dir / "05_viralquest"
        
        for d in [trimmed_dir, mapping_dir, filtered_dir, assembly_dir, viralquest_dir]:
            d.mkdir(exist_ok=True)
        
        return {
            'sample_dir': sample_dir,
            'trimmed_dir': trimmed_dir,
            'mapping_dir': mapping_dir,
            'filtered_dir': filtered_dir,
            'assembly_dir': assembly_dir,
            'viralquest_dir': viralquest_dir
        }
    
    def run_cutadapt(self, input_fastq, dirs):
        """
        Step 1: Trimming with cutadapt (automatic adapter detection)
        
        Args:
            input_fastq: Input FASTQ file
            dirs: Dictionary with sample directories
        
        Returns:
            Path to trimmed FASTQ file
        """
        logger.info("=" * 60)
        logger.info("STEP 1: Trimming with cutadapt")
        logger.info("=" * 60)
        
        input_fastq = Path(input_fastq)
        output_fastq = dirs['trimmed_dir'] / f"{input_fastq.stem}_trimmed.fastq"
        report_file = dirs['trimmed_dir'] / "cutadapt_report.txt"
        
        cmd = [
            "cutadapt",
            "-a", "auto",  # Automatic adapter detection
            "-j", str(self.threads),
            "-o", str(output_fastq),
            "--minimum-length", "18",  # Minimum for sRNAs
            "--maximum-length", "35",  # Maximum for sRNAs
            "-q", "20",  # Minimum quality
            str(input_fastq)
        ]
        
        logger.info(f"Command: {' '.join(cmd)}")
        
        with open(report_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, 
                          text=True, check=True)
        
        logger.info(f"✓ Trimming completed: {output_fastq}")
        logger.info(f"✓ Report saved: {report_file}")
        
        return output_fastq
    
    def map_to_reference_bowtie1(self, input_fastq, reference, ref_name, dirs):
        """
        Map reads against a reference genome with 0 mismatches using Bowtie1
        
        Args:
            input_fastq: Input FASTQ file
            reference: Reference genome
            ref_name: Reference name (for output files)
            dirs: Dictionary with sample directories
            
        Returns:
            FASTQ file with unmapped reads
        """
        logger.info(f"  Mapping against: {ref_name}")
        
        # Index reference if necessary
        if not Path(f"{reference}.1.ebwt").exists():
            logger.info(f"  Indexing {ref_name}...")
            index_cmd = ["bowtie-build", str(reference), str(reference)]
            subprocess.run(index_cmd, check=True, capture_output=True)
        
        # Output files
        sam_file = dirs['mapping_dir'] / f"{ref_name}.sam"
        unmapped_fastq = dirs['mapping_dir'] / f"unmapped_{ref_name}.fastq"
        stats_file = dirs['mapping_dir'] / f"{ref_name}_stats.txt"
        
        # Map with bowtie1 (0 mismatches, -v 0)
        # -v 0: report alignments with at most 0 mismatches
        # --un: write unmapped reads to file
        # -p: number of threads
        # -S: output in SAM format
        map_cmd = [
            "bowtie",
            "-v", "0",  # 0 mismatches allowed
            "-p", str(self.threads),
            "-S",  # SAM output
            "--un", str(unmapped_fastq),  # Unmapped reads
            str(reference),
            str(input_fastq),
            str(sam_file)
        ]
        
        logger.info(f"  Command: {' '.join(map_cmd)}")
        
        with open(stats_file, 'w') as f:
            subprocess.run(map_cmd, stderr=f, check=True)
        
        logger.info(f"  ✓ Unmapped reads: {unmapped_fastq}")
        
        # Remove SAM to save space
        sam_file.unlink()
        
        return unmapped_fastq
    
    def map_against_all_references(self, input_fastq, dirs):
        """
        Step 2: Map against all reference genomes
        
        Args:
            input_fastq: Input FASTQ file
            dirs: Dictionary with sample directories
            
        Returns:
            FASTQ file with final unmapped reads
        """
        logger.info("=" * 60)
        logger.info("STEP 2: Mapping against reference genomes")
        logger.info("=" * 60)
        
        current_fastq = input_fastq
        
        # 1. Map against host
        logger.info("→ Removing host reads...")
        current_fastq = self.map_to_reference_bowtie1(
            current_fastq, self.host_genome, "host", dirs
        )
        
        # 2. Map against bacterial RefSeq (optional)
        if self.bacterial_refseq:
            logger.info("→ Removing bacterial reads...")
            current_fastq = self.map_to_reference_bowtie1(
                current_fastq, self.bacterial_refseq, "bacterial", dirs
            )
        else:
            logger.info("→ Skipping bacterial mapping (not provided)")
        
        # 3. Map against transposons (optional)
        if self.transposon_db:
            logger.info("→ Removing transposon reads...")
            current_fastq = self.map_to_reference_bowtie1(
                current_fastq, self.transposon_db, "transposon", dirs
            )
        
        logger.info(f"✓ Final unmapped reads: {current_fastq}")
        return current_fastq
    
    def filter_by_size(self, input_fastq, dirs):
        """
        Step 3: Separate reads by size ranges
        
        Args:
            input_fastq: FASTQ file with unmapped reads
            dirs: Dictionary with sample directories
            
        Returns:
            Tuple with (fastq_20_23, fastq_24_30)
        """
        logger.info("=" * 60)
        logger.info("STEP 3: Size-based separation")
        logger.info("=" * 60)
        
        fastq_20_23 = dirs['filtered_dir'] / "srna_20-23nt.fastq"
        fastq_24_30 = dirs['filtered_dir'] / "srna_24-30nt.fastq"
        stats_file = dirs['filtered_dir'] / "size_stats.txt"
        
        count_20_23 = 0
        count_24_30 = 0
        count_other = 0
        
        with open(input_fastq, 'r') as f_in, \
             open(fastq_20_23, 'w') as f_20_23, \
             open(fastq_24_30, 'w') as f_24_30:
            
            while True:
                # Read 4 lines of FASTQ
                header = f_in.readline()
                if not header:
                    break
                seq = f_in.readline().strip()
                plus = f_in.readline()
                qual = f_in.readline()
                
                seq_len = len(seq)
                
                # Classify by size
                if 20 <= seq_len <= 23:
                    f_20_23.write(header)
                    f_20_23.write(seq + '\n')
                    f_20_23.write(plus)
                    f_20_23.write(qual)
                    count_20_23 += 1
                elif 24 <= seq_len <= 30:
                    f_24_30.write(header)
                    f_24_30.write(seq + '\n')
                    f_24_30.write(plus)
                    f_24_30.write(qual)
                    count_24_30 += 1
                else:
                    count_other += 1
        
        # Save statistics
        with open(stats_file, 'w') as f:
            f.write(f"Reads 20-23nt: {count_20_23}\n")
            f.write(f"Reads 24-30nt: {count_24_30}\n")
            f.write(f"Reads outside ranges: {count_other}\n")
            f.write(f"Total: {count_20_23 + count_24_30 + count_other}\n")
        
        logger.info(f"  20-23nt: {count_20_23} reads → {fastq_20_23}")
        logger.info(f"  24-30nt: {count_24_30} reads → {fastq_24_30}")
        logger.info(f"  Other: {count_other} reads (discarded)")
        logger.info(f"✓ Statistics saved: {stats_file}")
        
        return fastq_20_23, fastq_24_30
    
    def run_metaspades(self, input_fastq, size_range, dirs):
        """
        Step 4: Assembly with metaSPAdes
        
        Args:
            input_fastq: Input FASTQ file
            size_range: String indicating size range (e.g., "20-23")
            dirs: Dictionary with sample directories
            
        Returns:
            Path to assembled contigs
        """
        logger.info(f"  Assembling {size_range}nt reads with metaSPAdes...")
        
        output_subdir = dirs['assembly_dir'] / f"assembly_{size_range}nt"
        
        cmd = [
            "metaspades.py",
            "-s", str(input_fastq),
            "-o", str(output_subdir),
            "-k", "15",  # k-mer 15
            "-t", str(self.threads),
            "--only-assembler"  # Skip error correction for sRNAs
        ]
        
        logger.info(f"  Command: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        
        contigs = output_subdir / "contigs.fasta"
        scaffolds = output_subdir / "scaffolds.fasta"
        
        logger.info(f"  ✓ Contigs: {contigs}")
        logger.info(f"  ✓ Scaffolds: {scaffolds}")
        
        return contigs
    
    def run_viralquest(self, contigs_fasta, size_range, dirs):
        """
        Step 5: Viral identification with ViralQuest
        
        Args:
            contigs_fasta: Path to assembled contigs
            size_range: String indicating size range (e.g., "20-23")
            dirs: Dictionary with sample directories
            
        Returns:
            Path to ViralQuest output directory
        """
        logger.info(f"  Running ViralQuest for {size_range}nt contigs...")
        
        output_subdir = dirs['viralquest_dir'] / f"viralquest_{size_range}nt"
        
        # Get ViralQuest parameters from config
        vq = self.viralquest_params
        
        cmd = [
            "viralquest",
            "-in", str(contigs_fasta),
            "-out", str(output_subdir),
            "-ref", vq.get('refseq_virus', ''),
            "--blastn_online", vq.get('email', ''),
            "-x", vq.get('nr', ''),
            "-pfam", vq.get('pfam', ''),
            "-rvdb", vq.get('rvdb', ''),
            "-eggnog", vq.get('eggnog', ''),
            "-vfam", vq.get('vfam', ''),
            "-orf", str(vq.get('orf', '5')),
            "-cpu", str(self.threads),
            "--model-type", vq.get('model_type', 'google'),
            "--model-name", vq.get('model_name', 'gemini-3-flash-preview'),
            "--api-key", vq.get('api_key', '')
        ]
        
        # Remove empty parameters
        cmd = [x for x in cmd if x]
        
        logger.info(f"  Command: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        
        logger.info(f"  ✓ ViralQuest output: {output_subdir}")
        
        return output_subdir
    
    def assemble_all(self, fastq_20_23, fastq_24_30, dirs):
        """
        Step 4: Assemble both size ranges
        
        Args:
            fastq_20_23: FASTQ file with 20-23nt reads
            fastq_24_30: FASTQ file with 24-30nt reads
            dirs: Dictionary with sample directories
            
        Returns:
            Dictionary with assembly and ViralQuest results
        """
        logger.info("=" * 60)
        logger.info("STEP 4: Assembly with metaSPAdes (k-mer 15)")
        logger.info("=" * 60)
        
        # Assemble 20-23nt
        logger.info("→ Assembly 20-23nt:")
        contigs_20_23 = self.run_metaspades(fastq_20_23, "20-23", dirs)
        
        # Assemble 24-30nt
        logger.info("→ Assembly 24-30nt:")
        contigs_24_30 = self.run_metaspades(fastq_24_30, "24-30", dirs)
        
        logger.info("✓ Assembly completed!")
        
        results = {
            'contigs_20_23': contigs_20_23,
            'contigs_24_30': contigs_24_30
        }
        
        # Step 5: ViralQuest (if parameters provided)
        if self.viralquest_params:
            logger.info("\n" + "=" * 60)
            logger.info("STEP 5: Viral identification with ViralQuest")
            logger.info("=" * 60)
            
            logger.info("→ ViralQuest 20-23nt:")
            vq_20_23 = self.run_viralquest(contigs_20_23, "20-23", dirs)
            results['viralquest_20_23'] = vq_20_23
            
            logger.info("→ ViralQuest 24-30nt:")
            vq_24_30 = self.run_viralquest(contigs_24_30, "24-30", dirs)
            results['viralquest_24_30'] = vq_24_30
            
            logger.info("✓ ViralQuest completed!")
        
        return results
    
    def process_single_sample(self, input_fastq):
        """
        Process a single FASTQ file
        
        Args:
            input_fastq: Input FASTQ file path
            
        Returns:
            Dictionary with results
        """
        input_fastq = Path(input_fastq)
        sample_name = input_fastq.stem.replace('.fastq', '').replace('.fq', '')
        
        logger.info("\n" + "=" * 60)
        logger.info(f"PROCESSING SAMPLE: {sample_name}")
        logger.info("=" * 60)
        logger.info(f"Input: {input_fastq}")
        logger.info(f"Threads: {self.threads}")
        logger.info("=" * 60 + "\n")
        
        try:
            # Setup directories
            dirs = self.setup_sample_directories(sample_name)
            
            # Step 1: Trimming
            trimmed_fastq = self.run_cutadapt(input_fastq, dirs)
            
            # Step 2: Mapping
            unmapped_fastq = self.map_against_all_references(trimmed_fastq, dirs)
            
            # Step 3: Size filtering
            fastq_20_23, fastq_24_30 = self.filter_by_size(unmapped_fastq, dirs)
            
            # Step 4: Assembly (and optional Step 5: ViralQuest)
            assembly_results = self.assemble_all(
                fastq_20_23, fastq_24_30, dirs
            )
            
            # Final summary
            logger.info("\n" + "=" * 60)
            logger.info(f"SAMPLE {sample_name} COMPLETED SUCCESSFULLY!")
            logger.info("=" * 60)
            logger.info(f"Contigs 20-23nt: {assembly_results['contigs_20_23']}")
            logger.info(f"Contigs 24-30nt: {assembly_results['contigs_24_30']}")
            if 'viralquest_20_23' in assembly_results:
                logger.info(f"ViralQuest 20-23nt: {assembly_results['viralquest_20_23']}")
                logger.info(f"ViralQuest 24-30nt: {assembly_results['viralquest_24_30']}")
            logger.info(f"Output directory: {dirs['sample_dir']}")
            logger.info("=" * 60 + "\n")
            
            return {
                'sample_name': sample_name,
                'status': 'success',
                'output_dir': dirs['sample_dir'],
                **assembly_results
            }
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error executing command for {sample_name}: {e}")
            return {
                'sample_name': sample_name,
                'status': 'failed',
                'error': str(e)
            }
        except Exception as e:
            logger.error(f"Error processing {sample_name}: {e}")
            return {
                'sample_name': sample_name,
                'status': 'failed',
                'error': str(e)
            }
    
    def process_batch(self, input_files):
        """
        Process multiple FASTQ files
        
        Args:
            input_files: List of FASTQ file paths
            
        Returns:
            List of result dictionaries
        """
        logger.info("\n" + "=" * 60)
        logger.info("STARTING BATCH PROCESSING")
        logger.info("=" * 60)
        logger.info(f"Total samples: {len(input_files)}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Threads per sample: {self.threads}")
        logger.info("=" * 60 + "\n")
        
        results = []
        
        for i, input_file in enumerate(input_files, 1):
            logger.info(f"\n>>> Processing sample {i}/{len(input_files)}")
            result = self.process_single_sample(input_file)
            results.append(result)
        
        # Final batch summary
        logger.info("\n" + "=" * 60)
        logger.info("BATCH PROCESSING COMPLETED")
        logger.info("=" * 60)
        
        successful = sum(1 for r in results if r['status'] == 'success')
        failed = sum(1 for r in results if r['status'] == 'failed')
        
        logger.info(f"Total samples: {len(results)}")
        logger.info(f"Successful: {successful}")
        logger.info(f"Failed: {failed}")
        logger.info("=" * 60)
        
        # List successful samples
        if successful > 0:
            logger.info("\nSuccessful samples:")
            for r in results:
                if r['status'] == 'success':
                    logger.info(f"  ✓ {r['sample_name']}")
        
        # List failed samples
        if failed > 0:
            logger.info("\nFailed samples:")
            for r in results:
                if r['status'] == 'failed':
                    logger.info(f"  ✗ {r['sample_name']}: {r['error']}")
        
        logger.info("")
        
        return results


def find_fastq_files(input_path):
    """
    Find all FASTQ files in a directory or return single file
    
    Args:
        input_path: Path to a file or directory
        
    Returns:
        List of FASTQ file paths
    """
    input_path = Path(input_path)
    
    if input_path.is_file():
        # Single file
        return [input_path]
    elif input_path.is_dir():
        # Find all FASTQ files in directory
        fastq_extensions = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']
        fastq_files = []
        
        for ext in fastq_extensions:
            fastq_files.extend(input_path.glob(ext))
        
        fastq_files = sorted(fastq_files)
        
        if not fastq_files:
            logger.error(f"No FASTQ files found in directory: {input_path}")
            sys.exit(1)
        
        logger.info(f"Found {len(fastq_files)} FASTQ files in {input_path}")
        return fastq_files
    else:
        logger.error(f"Input path does not exist: {input_path}")
        sys.exit(1)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Viral sRNA assembly pipeline with optional ViralQuest analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (single file)
  python viral_srna_pipeline.py \\
    -i input.fastq \\
    -o output_dir \\
    -H host_genome.fa \\
    -t 8
  
  # With all optional features
  python viral_srna_pipeline.py \\
    -i fastq_directory/ \\
    -o output_dir \\
    -H host_genome.fa \\
    -B bacterial_refseq.fa \\
    -T transposon_db.fa \\
    --vq-refseq /path/to/refseq.dmnd \\
    --vq-nr /path/to/nr.dmnd \\
    --vq-pfam /path/to/Pfam-A.hmm \\
    --vq-email your@email.com \\
    --vq-api-key YOUR_API_KEY \\
    -t 8
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTQ file or directory with FASTQ files')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory')
    parser.add_argument('-H', '--host-genome', required=True,
                       help='Host reference genome')
    
    # Optional mapping databases
    parser.add_argument('-B', '--bacterial-refseq', default=None,
                       help='Bacterial RefSeq (optional)')
    parser.add_argument('-T', '--transposon-db', default=None,
                       help='Transposon database (optional)')
    
    # ViralQuest parameters
    vq_group = parser.add_argument_group('ViralQuest options (all optional)')
    vq_group.add_argument('--vq-refseq', dest='vq_refseq',
                         help='RefSeq virus database (.dmnd)')
    vq_group.add_argument('--vq-email', dest='vq_email',
                         help='Email for BLASTN online')
    vq_group.add_argument('--vq-nr', dest='vq_nr',
                         help='NR database (.dmnd)')
    vq_group.add_argument('--vq-pfam', dest='vq_pfam',
                         help='Pfam-A HMM database')
    vq_group.add_argument('--vq-rvdb', dest='vq_rvdb',
                         help='RVDB HMM database')
    vq_group.add_argument('--vq-eggnog', dest='vq_eggnog',
                         help='EggNOG HMM database')
    vq_group.add_argument('--vq-vfam', dest='vq_vfam',
                         help='Vfam HMM database')
    vq_group.add_argument('--vq-orf', dest='vq_orf', type=int, default=5,
                         help='Minimum ORF length (default: 5)')
    vq_group.add_argument('--vq-model-type', dest='vq_model_type', default='google',
                         help='Model type (default: google)')
    vq_group.add_argument('--vq-model-name', dest='vq_model_name', 
                         default='gemini-3-flash-preview',
                         help='Model name (default: gemini-3-flash-preview)')
    vq_group.add_argument('--vq-api-key', dest='vq_api_key',
                         help='API key for model')
    
    # General options
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help='Number of threads (default: 4)')
    
    args = parser.parse_args()
    
    # Find input files
    input_files = find_fastq_files(args.input)
    
    # Prepare ViralQuest parameters (only if at least refseq is provided)
    viralquest_params = None
    if args.vq_refseq:
        viralquest_params = {
            'refseq_virus': args.vq_refseq,
            'email': args.vq_email,
            'nr': args.vq_nr,
            'pfam': args.vq_pfam,
            'rvdb': args.vq_rvdb,
            'eggnog': args.vq_eggnog,
            'vfam': args.vq_vfam,
            'orf': args.vq_orf,
            'model_type': args.vq_model_type,
            'model_name': args.vq_model_name,
            'api_key': args.vq_api_key
        }
        logger.info("ViralQuest analysis will be performed after assembly")
    
    # Create and run pipeline
    pipeline = ViralSRNAPipeline(
        output_dir=args.output,
        host_genome=args.host_genome,
        bacterial_refseq=args.bacterial_refseq,
        transposon_db=args.transposon_db,
        threads=args.threads,
        viralquest_params=viralquest_params
    )
    
    # Process samples
    if len(input_files) == 1:
        pipeline.process_single_sample(input_files[0])
    else:
        pipeline.process_batch(input_files)


if __name__ == '__main__':
    main()