#!/usr/bin/env python3
"""
Phasing and Offset Analysis Module
Implements piRNA ping-pong signature detection

Based on:
  Muerdter et al., 2012 - Production of artificial piRNAs in flies and mice
  RNA (2012) doi: 10.1261/rna.029769.111
"""

import subprocess
from pathlib import Path
from collections import defaultdict
from typing import List, Dict, Optional, Tuple
import gzip

from utils import log, ensure_outdir, which


class PhasingAnalysis:
    """
    Phasing and Offset analysis for piRNA ping-pong signature detection.
    
    This class implements the core algorithm from bed_calculate_offsets.pl
    but works directly with SAM files, eliminating intermediate BED conversions.
    """
    
    def __init__(
        self,
        sam_file: Path,
        region_chrom: Optional[str] = None,
        region_start: Optional[int] = None,
        region_end: Optional[int] = None
    ):
        """
        Initialize phasing analysis for a specific genomic region.
        
        Args:
            sam_file: Path to SAM file
            region_chrom: Chromosome/contig name (None = auto-detect first chromosome)
            region_start: Start position (0-based, None = start of chromosome)
            region_end: End position (exclusive, None = end of chromosome)
        """
        self.sam_file = Path(sam_file)
        self.region_chrom = region_chrom
        self.region_start = region_start
        self.region_end = region_end
        
        # Will be set during coverage calculation if not provided
        self.auto_detected_region = False
        self.region_size = None
        
        # Coverage arrays (will be initialized once region is known)
        self.window_positive_sequences = None
        self.window_negative_sequences = None
        self.window_positive_reads = None
        self.window_negative_reads = None
        
        # Store intervals for offset calculation
        self.intervals = []
        
        # Statistics
        self.stats_total_lines = 0
        self.stats_in_region = 0
        self.stats_total_reads = 0
    
    def _open_sam(self):
        """Open SAM file (handles gzipped files)."""
        if str(self.sam_file).endswith('.gz'):
            return gzip.open(self.sam_file, 'rt')
        return open(self.sam_file, 'r')
    
    @staticmethod
    def _extract_multiplicity(read_name: str) -> int:
        """
        Extract multiplicity from read name.
        Format: N-M where M is the multiplicity count.
        
        Args:
            read_name: Read name from SAM
        
        Returns:
            Multiplicity value (default: 1)
        """
        if '-' in read_name:
            try:
                parts = read_name.rsplit('-', 1)
                multiplicity = int(parts[-1])
                return multiplicity
            except (ValueError, IndexError):
                pass
        return 1
    
    def _detect_region_from_sam(
        self,
        size_min: Optional[int] = None,
        size_max: Optional[int] = None
    ) -> Tuple[str, int, int]:
        """
        Auto-detect region from SAM file.
        
        Scans the SAM to find:
        - First chromosome with reads
        - Min and max positions of reads
        
        Args:
            size_min: Minimum read length filter
            size_max: Maximum read length filter
        
        Returns:
            Tuple of (chrom, start, end)
        """
        log("Auto-detecting region from SAM file...", "INFO")
        
        chrom_positions = {}  # {chrom: [list of positions]}
        
        with self._open_sam() as f:
            for line in f:
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                
                flag = int(parts[1])
                chrom = parts[2]
                start = int(parts[3]) - 1  # Convert to 0-based
                seq = parts[9]
                
                # Skip unmapped
                if chrom == '*':
                    continue
                
                # Filter by size if requested
                read_len = len(seq)
                if size_min and read_len < size_min:
                    continue
                if size_max and read_len > size_max:
                    continue
                
                # Calculate end
                end = start + read_len
                
                # Store positions for this chromosome
                if chrom not in chrom_positions:
                    chrom_positions[chrom] = []
                
                chrom_positions[chrom].append(start)
                chrom_positions[chrom].append(end)
        
        if not chrom_positions:
            raise ValueError("No reads found in SAM file")
        
        # Use the chromosome with most reads
        best_chrom = max(chrom_positions.keys(), key=lambda k: len(chrom_positions[k]))
        positions = chrom_positions[best_chrom]
        
        region_start = min(positions)
        region_end = max(positions)
        
        log(f"Auto-detected region: {best_chrom}:{region_start}-{region_end}", "INFO")
        log(f"Region size: {region_end - region_start:,} bp", "INFO")
        
        return best_chrom, region_start, region_end
    
    @staticmethod
    def get_all_contigs_from_sam(
        sam_file: Path,
        size_min: Optional[int] = None,
        size_max: Optional[int] = None,
        min_reads: int = 10
    ) -> Dict[str, Tuple[int, int]]:
        """
        Get all contigs/chromosomes from SAM file with their regions.
        
        Args:
            sam_file: Path to SAM file
            size_min: Minimum read length filter
            size_max: Maximum read length filter
            min_reads: Minimum number of reads to include contig
        
        Returns:
            Dictionary {contig: (start, end)}
        """
        log("Scanning SAM for all contigs...", "INFO")
        
        chrom_positions = {}  # {chrom: [list of positions]}
        
        # Open file
        if str(sam_file).endswith('.gz'):
            f = gzip.open(sam_file, 'rt')
        else:
            f = open(sam_file, 'r')
        
        try:
            for line in f:
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                
                flag = int(parts[1])
                chrom = parts[2]
                start = int(parts[3]) - 1  # Convert to 0-based
                seq = parts[9]
                
                # Skip unmapped
                if chrom == '*':
                    continue
                
                # Filter by size if requested
                read_len = len(seq)
                if size_min and read_len < size_min:
                    continue
                if size_max and read_len > size_max:
                    continue
                
                # Calculate end
                end = start + read_len
                
                # Store positions for this chromosome
                if chrom not in chrom_positions:
                    chrom_positions[chrom] = []
                
                chrom_positions[chrom].append(start)
                chrom_positions[chrom].append(end)
        
        finally:
            f.close()
        
        if not chrom_positions:
            raise ValueError("No reads found in SAM file")
        
        # Calculate regions for each contig
        contigs = {}
        for chrom, positions in chrom_positions.items():
            # Filter by minimum reads (each read adds 2 positions)
            num_reads = len(positions) // 2
            if num_reads < min_reads:
                log(f"  Skipping {chrom}: only {num_reads} reads (< {min_reads})", "DEBUG")
                continue
            
            region_start = min(positions)
            region_end = max(positions)
            region_size = region_end - region_start
            
            contigs[chrom] = (region_start, region_end)
            log(f"  Found {chrom}: {region_start}-{region_end} ({region_size:,} bp, {num_reads} reads)", "INFO")
        
        log(f"Found {len(contigs)} contigs with >= {min_reads} reads", "INFO")
        
        return contigs
    
    def _initialize_coverage_arrays(self):
        """Initialize coverage arrays based on region size."""
        if self.region_size is None:
            raise ValueError("Region size not set")
        
        self.window_positive_sequences = [0] * self.region_size
        self.window_negative_sequences = [0] * self.region_size
        self.window_positive_reads = [0] * self.region_size
        self.window_negative_reads = [0] * self.region_size
    
    def calculate_coverage_from_sam(
        self,
        size_min: Optional[int] = None,
        size_max: Optional[int] = None,
        use_multiplicity: bool = False
    ) -> None:
        """
        Calculate coverage by counting 5' ends at each position.
        
        This is equivalent to the coverage calculation in bed_calculate_offsets.pl
        but works directly with SAM files.
        
        If region is not specified, auto-detects from SAM file.
        
        Args:
            size_min: Minimum read length to include (None = no filter)
            size_max: Maximum read length to include (None = no filter)
            use_multiplicity: Extract multiplicity from read name
        """
        # Auto-detect region if not provided
        if self.region_chrom is None or self.region_start is None or self.region_end is None:
            self.region_chrom, self.region_start, self.region_end = self._detect_region_from_sam(
                size_min=size_min,
                size_max=size_max
            )
            self.auto_detected_region = True
        
        # Calculate region size and initialize arrays
        self.region_size = self.region_end - self.region_start
        self._initialize_coverage_arrays()
        
        log(f"Calculating coverage for region {self.region_chrom}:{self.region_start}-{self.region_end}", "INFO")
        log(f"Size filter: {size_min}-{size_max} nt" if size_min or size_max else "No size filter", "DEBUG")
        
        seen_reads = set()
        
        with self._open_sam() as f:
            for line in f:
                self.stats_total_lines += 1
                
                if line.startswith('@'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                
                read_name = parts[0]
                flag = int(parts[1])
                chrom = parts[2]
                start = int(parts[3]) - 1  # SAM is 1-based, convert to 0-based
                seq = parts[9]
                
                # Skip duplicates (by name)
                if read_name in seen_reads:
                    continue
                seen_reads.add(read_name)
                
                # Filter by chromosome
                if chrom != self.region_chrom:
                    continue
                
                # Skip unmapped
                if chrom == '*':
                    continue
                
                # Calculate read length and end position
                read_len = len(seq)
                end = start + read_len
                
                # Filter by size
                if size_min and read_len < size_min:
                    continue
                if size_max and read_len > size_max:
                    continue
                
                # Filter by region
                if start < self.region_start or end > self.region_end:
                    continue
                
                self.stats_in_region += 1
                
                # Determine strand
                is_reverse = bool(flag & 0x10)
                strand = '-' if is_reverse else '+'
                
                # Extract multiplicity if requested
                multiplicity = 1
                if use_multiplicity:
                    multiplicity = self._extract_multiplicity(read_name)
                
                self.stats_total_reads += multiplicity
                
                # Calculate 5' position
                # For forward strand: 5' is at start
                # For reverse strand: 5' is at end-1
                if strand == '+':
                    pos_5prime = start
                else:
                    pos_5prime = end - 1
                
                # Update coverage arrays
                index = pos_5prime - self.region_start
                if 0 <= index < self.region_size:
                    if strand == '+':
                        self.window_positive_sequences[index] += 1
                        self.window_positive_reads[index] += multiplicity
                    else:
                        self.window_negative_sequences[index] += 1
                        self.window_negative_reads[index] += multiplicity
                    
                    # Store interval for offset calculation
                    self.intervals.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'name': read_name,
                        'strand': strand,
                        'multiplicity': multiplicity
                    })
        
        log(f"Processed {self.stats_total_lines} lines", "DEBUG")
        log(f"Found {self.stats_in_region} reads in region ({self.stats_total_reads} with multiplicity)", "INFO")
    
    def calculate_offsets(
        self,
        upstream: int = 80,
        downstream: int = 80
    ) -> List[Dict]:
        """
        Calculate offset data for all intervals.
        
        For each read, scans a window around its 5' position and counts
        how many other reads have their 5' ends at each offset position.
        
        Args:
            upstream: Window size upstream of 5' end
            downstream: Window size downstream of 5' end
        
        Returns:
            List of offset data dictionaries
        """
        log(f"Calculating offsets (window: -{upstream} to +{downstream})", "INFO")
        
        results = []
        
        # Sort intervals by strand, then start, then end
        self.intervals.sort(key=lambda x: (x['strand'], x['start'], x['end']))
        
        for interval in self.intervals:
            chrom = interval['chrom']
            start = interval['start']
            end = interval['end']
            strand = interval['strand']
            multiplicity = interval['multiplicity']
            
            # Calculate 5' position relative to region start
            if strand == '+':
                nuc_position = start - self.region_start
            else:
                nuc_position = end - 1 - self.region_start
            
            # Calculate genomic position
            genomic_nuc_position = self.region_start + nuc_position
            
            # Define window around the 5' position
            # Note: upstream/downstream are relative to strand direction
            if strand == '+':
                extended_start = nuc_position - upstream
                extended_end = nuc_position + downstream
            else:
                # For reverse strand, swap upstream/downstream
                extended_start = nuc_position - downstream
                extended_end = nuc_position + upstream
            
            # Clip to region boundaries
            extended_start = max(0, extended_start)
            extended_end = min(self.region_size - 1, extended_end)
            
            # Scan the window
            for window_index in range(extended_start, extended_end + 1):
                # Calculate offset relative to the interval's 5' position
                if strand == '+':
                    offset = window_index - nuc_position
                else:
                    offset = nuc_position - window_index
                
                # Count reads at this position
                if strand == '+':
                    seq_same = self.window_positive_sequences[window_index]
                    read_same = self.window_positive_reads[window_index]
                    seq_opposite = self.window_negative_sequences[window_index]
                    read_opposite = self.window_negative_reads[window_index]
                else:
                    seq_same = self.window_negative_sequences[window_index]
                    read_same = self.window_negative_reads[window_index]
                    seq_opposite = self.window_positive_sequences[window_index]
                    read_opposite = self.window_positive_reads[window_index]
                
                # At offset 0, subtract the current interval from same-strand counts
                if offset == 0:
                    seq_same -= 1
                    read_same -= multiplicity
                
                # Only record if there are reads to report
                if seq_same > 0 or read_same > 0 or seq_opposite > 0 or read_opposite > 0:
                    results.append({
                        'chrom': chrom,
                        'genomic_start': genomic_nuc_position,
                        'strand': strand,
                        'window_position': nuc_position,
                        'offset': offset,
                        'seq_count_same_strand': seq_same,
                        'read_count_same_strand': read_same,
                        'seq_count_opposite_strand': seq_opposite,
                        'read_count_opposite_strand': read_opposite
                    })
        
        log(f"Generated {len(results)} offset data points", "INFO")
        return results
    
    def aggregate_phasing(
        self,
        offset_data: List[Dict],
        by_strand: str = '+'
    ) -> List[Dict]:
        """
        Aggregate phasing data by offset for a specific strand.
        
        This groups by offset and sums read_count_same_strand.
        Equivalent to: sort -k5,5 | groupBy -g 5 -c 7 -o sum
        
        Args:
            offset_data: List of offset dictionaries
            by_strand: Strand to filter by ('+' or '-')
        
        Returns:
            List of {'offset': int, 'count': int} dictionaries
        """
        aggregated = defaultdict(int)
        
        for row in offset_data:
            if row['strand'] == by_strand:
                offset = row['offset']
                count = row['read_count_same_strand']
                aggregated[offset] += count
        
        # Convert to sorted list
        result = [
            {'offset': offset, 'count': aggregated[offset]}
            for offset in sorted(aggregated.keys())
        ]
        
        log(f"Aggregated phasing for strand '{by_strand}': {len(result)} offsets", "DEBUG")
        return result
    
    def aggregate_offset(
        self,
        offset_data: List[Dict]
    ) -> List[Dict]:
        """
        Aggregate overlap between sense and antisense strands.
        
        This groups by offset and sums read_count_opposite_strand.
        Equivalent to: sort -k5,5 | groupBy -g 5 -c 9 -o sum
        
        Args:
            offset_data: List of offset dictionaries
        
        Returns:
            List of {'offset': int, 'count': int} dictionaries
        """
        aggregated = defaultdict(int)
        
        for row in offset_data:
            offset = row['offset']
            count = row['read_count_opposite_strand']
            aggregated[offset] += count
        
        # Convert to sorted list
        result = [
            {'offset': offset, 'count': aggregated[offset]}
            for offset in sorted(aggregated.keys())
        ]
        
        log(f"Aggregated offset: {len(result)} offsets", "DEBUG")
        return result
    
    def write_phasing_file(
        self,
        phasing_data: List[Dict],
        output_file: Path
    ) -> bool:
        """
        Write phasing data to file.
        
        Format: Two tab-separated columns (no header)
        offset    count
        
        Args:
            phasing_data: List of phasing dictionaries
            output_file: Output file path
        
        Returns:
            True if successful
        """
        try:
            ensure_outdir(output_file.parent)
            
            with open(output_file, 'w') as f:
                for row in phasing_data:
                    f.write(f"{row['offset']}\t{row['count']}\n")
            
            log(f"Wrote phasing file: {output_file}", "DEBUG")
            return True
            
        except Exception as e:
            log(f"Error writing phasing file: {e}", "ERROR")
            return False
    
    def write_offset_file(
        self,
        offset_data: List[Dict],
        output_file: Path
    ) -> bool:
        """
        Write offset data to file.
        
        Format: Tab-separated WITHOUT header (to match BEDtools groupBy output)
        offset    count
        
        Args:
            offset_data: List of offset dictionaries
            output_file: Output file path
        
        Returns:
            True if successful
        """
        try:
            ensure_outdir(output_file.parent)
            
            with open(output_file, 'w') as f:
                # NO HEADER - matches BEDtools groupBy output format
                for row in offset_data:
                    f.write(f"{row['offset']}\t{row['count']}\n")
            
            log(f"Wrote offset file: {output_file}", "DEBUG")
            return True
            
        except Exception as e:
            log(f"Error writing offset file: {e}", "ERROR")
            return False
    
    def plot_with_r(
        self,
        plus_phasing_file: Optional[Path],
        minus_phasing_file: Optional[Path],
        offset_file: Path,
        output_prefix: Path,
        r_script_path: Path
    ) -> bool:
        """
        Generate plot using R script.
        
        Supports two R script formats:
        1. 2 args: Rscript script.R <offset_file> <output_prefix>
        2. 4 args: Rscript script.R <plus_phasing> <minus_phasing> <offset_file> <output_prefix>
        
        Args:
            plus_phasing_file: Path to plus phasing file (optional for 2-arg scripts)
            minus_phasing_file: Path to minus phasing file (optional for 2-arg scripts)
            offset_file: Path to offset data file
            output_prefix: Output prefix for PDF
            r_script_path: Path to R plotting script
        
        Returns:
            True if successful
        """
        if not which("Rscript"):
            log("Rscript not found in PATH", "ERROR")
            return False
        
        if not r_script_path.exists():
            log(f"R script not found: {r_script_path}", "ERROR")
            return False
        
        if not offset_file.exists():
            log(f"Offset file not found: {offset_file}", "ERROR")
            return False
        
        # Try 4-arg format first (if phasing files provided)
        if plus_phasing_file and minus_phasing_file and plus_phasing_file.exists() and minus_phasing_file.exists():
            cmd = [
                "Rscript",
                str(r_script_path),
                str(plus_phasing_file),
                str(minus_phasing_file),
                str(offset_file),
                str(output_prefix)
            ]
            
            log(f"Running R plot (4-arg): {' '.join(cmd)}", "INFO")
            
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=300
                )
                
                if result.returncode == 0:
                    pdf_file = Path(str(output_prefix) + ".pdf")
                    if pdf_file.exists():
                        log(f"Generated plot: {pdf_file}", "INFO")
                        return True
                    else:
                        log("R script succeeded but PDF not found", "WARN")
                        log(f"R stdout: {result.stdout}", "DEBUG")
                else:
                    # Check if it's a data-related warning (expected) or real error
                    stderr_lower = result.stderr.lower()
                    
                    # Expected data warnings (biologicall normal, don't show as errors)
                    expected_warnings = [
                        "no non-missing arguments",
                        "need finite 'ylim' values",
                        "no data",
                        "empty",
                        "returning inf",
                        "returning -inf"
                    ]
                    
                    is_data_warning = any(warn in stderr_lower for warn in expected_warnings)
                    
                    if is_data_warning:
                        # Data issue (e.g., empty minus strand), try 2-arg
                        log("Note: Some data files empty (normal for unidirectional clusters), trying 2-arg format", "INFO")
                    elif ("argument" in stderr_lower and "expected" in stderr_lower) or \
                       ("usage:" in stderr_lower and "2" in stderr_lower):
                        # Wrong argument count, try 2-arg
                        log("4-arg format not supported, trying 2-arg", "DEBUG")
                    else:
                        # Real error
                        log(f"R script failed (4-arg): {result.stderr}", "ERROR")
                        if result.stdout:
                            log(f"R stdout: {result.stdout}", "DEBUG")
                        return False
                    
            except subprocess.TimeoutExpired:
                log("R script timeout", "ERROR")
                return False
            except Exception as e:
                log(f"Error running R script (4-arg): {e}", "DEBUG")
        
        # Try 2-arg format
        cmd = [
            "Rscript",
            str(r_script_path),
            str(offset_file),
            str(output_prefix)
        ]
        
        log(f"Running R plot (2-arg): {' '.join(cmd)}", "INFO")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if result.returncode == 0:
                pdf_file = Path(str(output_prefix) + ".pdf")
                if pdf_file.exists():
                    log(f"Generated plot: {pdf_file}", "INFO")
                    return True
                else:
                    log("R script succeeded but PDF not found", "WARN")
                    return False
            else:
                # Check if it's expected data warning or real error
                stderr_lower = result.stderr.lower()
                
                expected_warnings = [
                    "no non-missing arguments",
                    "need finite 'ylim' values",
                    "no data",
                    "empty",
                    "returning inf",
                    "returning -inf"
                ]
                
                is_data_warning = any(warn in stderr_lower for warn in expected_warnings)
                
                if is_data_warning:
                    log("Note: Insufficient data in some files (normal for unidirectional data)", "WARN")
                    # Still return False as plot wasn't generated
                    return False
                else:
                    log(f"R script failed: {result.stderr}", "ERROR")
                    return False
                
        except subprocess.TimeoutExpired:
            log("R script timeout", "ERROR")
            return False
        except Exception as e:
            log(f"Error running R script (2-arg): {e}", "ERROR")
            return False


def run_phasing_analysis_all_contigs(
    sam_file: Path,
    output_dir: Path,
    prefix: str,
    size_min: int = 24,
    size_max: int = 30,
    upstream: int = 80,
    downstream: int = 80,
    use_multiplicity: bool = False,
    r_script_path: Optional[Path] = None,
    min_reads: int = 10
) -> Dict[str, Dict]:
    """
    Run phasing analysis on all contigs/chromosomes in SAM file.
    
    Args:
        sam_file: Input SAM file
        output_dir: Output directory
        prefix: Output file prefix
        size_min: Minimum read size
        size_max: Maximum read size
        upstream: Window upstream
        downstream: Window downstream
        use_multiplicity: Use multiplicity from read names
        r_script_path: Path to R plotting script
        min_reads: Minimum reads per contig to analyze
    
    Returns:
        Dictionary {contig: result_dict}
    """
    log("="*60, "INFO")
    log(f"Phasing Analysis - ALL CONTIGS MODE", "INFO")
    log(f"Prefix: {prefix}", "INFO")
    log(f"Size range: {size_min}-{size_max} nt", "INFO")
    log(f"Min reads per contig: {min_reads}", "INFO")
    log("="*60, "INFO")
    
    # Get all contigs
    contigs = PhasingAnalysis.get_all_contigs_from_sam(
        sam_file=sam_file,
        size_min=size_min,
        size_max=size_max,
        min_reads=min_reads
    )
    
    if not contigs:
        log("No contigs found with sufficient reads", "WARN")
        return {}
    
    log(f"\nWill analyze {len(contigs)} contigs", "INFO")
    
    # Process each contig
    all_results = {}
    successful = 0
    
    for idx, (contig, (contig_start, contig_end)) in enumerate(contigs.items(), 1):
        log(f"\n[{idx}/{len(contigs)}] Processing {contig}...", "INFO")
        
        try:
            result = run_phasing_analysis(
                sam_file=sam_file,
                output_dir=output_dir,
                prefix=f"{prefix}.{contig}",
                region_chrom=contig,
                region_start=contig_start,
                region_end=contig_end,
                size_min=size_min,
                size_max=size_max,
                upstream=upstream,
                downstream=downstream,
                use_multiplicity=use_multiplicity,
                r_script_path=r_script_path
            )
            
            all_results[contig] = result
            
            if result and result.get('offset'):
                successful += 1
                log(f"  ✓ {contig} completed successfully", "INFO")
            else:
                log(f"  ✗ {contig} failed or no data", "WARN")
        
        except Exception as e:
            log(f"  ✗ {contig} failed: {e}", "ERROR")
            all_results[contig] = None
    
    # Summary
    log("\n" + "="*60, "INFO")
    log(f"All Contigs Analysis Summary", "INFO")
    log(f"Total contigs: {len(contigs)}", "INFO")
    log(f"Successful: {successful}", "INFO")
    log(f"Failed: {len(contigs) - successful}", "INFO")
    log("="*60, "INFO")
    
    return all_results


def run_phasing_analysis(
    sam_file: Path,
    output_dir: Path,
    prefix: str,
    region_chrom: Optional[str] = None,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
    size_min: int = 24,
    size_max: int = 30,
    upstream: int = 80,
    downstream: int = 80,
    use_multiplicity: bool = False,
    r_script_path: Optional[Path] = None
) -> Dict[str, Path]:
    """
    Run complete phasing analysis pipeline.
    
    If region is not specified, automatically detects region from SAM file.
    
    Args:
        sam_file: Input SAM file
        output_dir: Output directory
        prefix: Output file prefix
        region_chrom: Chromosome name (None = auto-detect)
        region_start: Region start (None = auto-detect)
        region_end: Region end (None = auto-detect)
        size_min: Minimum read size
        size_max: Maximum read size
        upstream: Window upstream
        downstream: Window downstream
        use_multiplicity: Use multiplicity from read names
        r_script_path: Path to R plotting script
    
    Returns:
        Dictionary with paths to generated files
    """
    log("="*60, "INFO")
    log(f"Phasing Analysis: {prefix}", "INFO")
    
    if region_chrom and region_start is not None and region_end is not None:
        log(f"Region: {region_chrom}:{region_start}-{region_end}", "INFO")
    else:
        log("Region: AUTO-DETECT from SAM", "INFO")
    
    log(f"Size range: {size_min}-{size_max} nt", "INFO")
    log("="*60, "INFO")
    
    # Initialize analysis
    analysis = PhasingAnalysis(
        sam_file=sam_file,
        region_chrom=region_chrom,
        region_start=region_start,
        region_end=region_end
    )
    
    # Calculate coverage (will auto-detect region if needed)
    analysis.calculate_coverage_from_sam(
        size_min=size_min,
        size_max=size_max,
        use_multiplicity=use_multiplicity
    )
    
    # Calculate offsets
    offset_data = analysis.calculate_offsets(
        upstream=upstream,
        downstream=downstream
    )
    
    if not offset_data:
        log("No offset data generated", "WARN")
        return {}
    
    # Aggregate data
    phasing_plus = analysis.aggregate_phasing(offset_data, by_strand='+')
    phasing_minus = analysis.aggregate_phasing(offset_data, by_strand='-')
    offset_agg = analysis.aggregate_offset(offset_data)
    
    # Write output files
    ensure_outdir(output_dir)
    
    size_tag = f"{size_min}-{size_max}" if size_min != size_max else str(size_min)
    
    # Add chromosome to prefix if auto-detected
    if analysis.auto_detected_region:
        file_prefix = f"{prefix}.{analysis.region_chrom}"
    else:
        file_prefix = prefix
    
    plus_phasing_file = output_dir / f"{file_prefix}.plus.{size_tag}.phasing"
    minus_phasing_file = output_dir / f"{file_prefix}.minus.{size_tag}.phasing"
    offset_file = output_dir / f"{file_prefix}.{size_tag}.offset"
    
    analysis.write_phasing_file(phasing_plus, plus_phasing_file)
    analysis.write_phasing_file(phasing_minus, minus_phasing_file)
    analysis.write_offset_file(offset_agg, offset_file)
    
    # Generate plot
    plot_success = False
    if r_script_path and r_script_path.exists():
        output_prefix = output_dir / f"{file_prefix}.{size_tag}"
        plot_success = analysis.plot_with_r(
            plus_phasing_file=plus_phasing_file,
            minus_phasing_file=minus_phasing_file,
            offset_file=offset_file,
            output_prefix=output_prefix,
            r_script_path=r_script_path
        )
    else:
        log("R script not provided or not found, skipping plot", "WARN")
    
    log("="*60, "INFO")
    log("Phasing analysis completed", "INFO")
    if analysis.auto_detected_region:
        log(f"Analyzed region: {analysis.region_chrom}:{analysis.region_start}-{analysis.region_end}", "INFO")
    log("="*60, "INFO")
    
    return {
        'plus_phasing': plus_phasing_file,
        'minus_phasing': minus_phasing_file,
        'offset': offset_file,
        'plot': output_dir / f"{file_prefix}.{size_tag}.pdf" if plot_success else None,
        'region': {
            'chrom': analysis.region_chrom,
            'start': analysis.region_start,
            'end': analysis.region_end
        }
    }


if __name__ == "__main__":
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(
        description="Phasing and offset analysis for piRNA ping-pong signature detection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  1. Specific region:  Provide --chrom and --region
  2. Auto-detect:      Omit --chrom and --region (uses contig with most reads)
  3. All contigs:      Use --all-contigs flag

Examples:
  # Mode 1: Analyze specific region
  %(prog)s --sam reads.sam -p sample1 -o results/ \\
      --chrom chr1 --region 1000-10000

  # Mode 2: Auto-detect best contig
  %(prog)s --sam reads.sam -p sample1 -o results/

  # Mode 3: Analyze ALL contigs (good for small genomes/scaffolds)
  %(prog)s --sam reads.sam -p sample1 -o results/ --all-contigs

  # With size filter and R plotting
  %(prog)s --sam reads.sam -p sample1 -o results/ --all-contigs \\
      --size-min 24 --size-max 30 \\
      --r-script R/plotOffsetPhasingSeries.R

  # Set minimum reads per contig
  %(prog)s --sam reads.sam -p sample1 -o results/ --all-contigs \\
      --min-reads 50
"""
    )
    
    parser.add_argument(
        "--sam", "-s",
        required=True,
        help="Input SAM file"
    )
    parser.add_argument(
        "--prefix", "-p",
        required=True,
        help="Output prefix"
    )
    parser.add_argument(
        "--outdir", "-o",
        required=True,
        help="Output directory"
    )
    
    # Mode selection
    parser.add_argument(
        "--all-contigs",
        action="store_true",
        help="Analyze all contigs/chromosomes (good for fragmented genomes)"
    )
    
    # Region specification (mode 1)
    parser.add_argument(
        "--chrom", "-c",
        help="Chromosome/contig name (optional, auto-detect if not provided)"
    )
    parser.add_argument(
        "--region", "-r",
        help="Region in format 'start-end' (optional, auto-detect if not provided)"
    )
    
    # Analysis parameters
    parser.add_argument(
        "--size-min",
        type=int,
        default=24,
        help="Minimum read size (default: %(default)s)"
    )
    parser.add_argument(
        "--size-max",
        type=int,
        default=30,
        help="Maximum read size (default: %(default)s)"
    )
    parser.add_argument(
        "--upstream",
        type=int,
        default=80,
        help="Window upstream (default: %(default)s)"
    )
    parser.add_argument(
        "--downstream",
        type=int,
        default=80,
        help="Window downstream (default: %(default)s)"
    )
    parser.add_argument(
        "--min-reads",
        type=int,
        default=10,
        help="Minimum reads per contig (for --all-contigs mode, default: %(default)s)"
    )
    parser.add_argument(
        "--multiplicity",
        action="store_true",
        help="Use multiplicity from read names (format: N-M)"
    )
    parser.add_argument(
        "--r-script",
        help="Path to R plotting script"
    )
    
    args = parser.parse_args()
    
    # Validate mode selection
    if args.all_contigs and (args.chrom or args.region):
        log("ERROR: Cannot use --all-contigs with --chrom or --region", "ERROR")
        log("Choose either --all-contigs OR specify region with --chrom/--region", "INFO")
        sys.exit(1)
    
    # Parse region if provided
    region_chrom = None
    region_start = None
    region_end = None
    
    if args.chrom and args.region:
        region_chrom = args.chrom
        try:
            region_start, region_end = map(int, args.region.split('-'))
        except:
            log(f"Invalid region format: {args.region}", "ERROR")
            sys.exit(1)
    elif args.chrom or args.region:
        log("ERROR: Both --chrom and --region must be provided together, or neither", "ERROR")
        log("Use neither for auto-detection, or use --all-contigs for all contigs", "INFO")
        sys.exit(1)
    
    # Run analysis
    try:
        if args.all_contigs:
            # Mode 3: All contigs
            log("Running in ALL CONTIGS mode", "INFO")
            results = run_phasing_analysis_all_contigs(
                sam_file=Path(args.sam),
                output_dir=Path(args.outdir),
                prefix=args.prefix,
                size_min=args.size_min,
                size_max=args.size_max,
                upstream=args.upstream,
                downstream=args.downstream,
                use_multiplicity=args.multiplicity,
                r_script_path=Path(args.r_script) if args.r_script else None,
                min_reads=args.min_reads
            )
            
            if results:
                log("\nResults summary:", "INFO")
                successful = sum(1 for r in results.values() if r and r.get('offset'))
                log(f"  Analyzed: {len(results)} contigs", "INFO")
                log(f"  Successful: {successful}", "INFO")
                log(f"  Failed: {len(results) - successful}", "INFO")
        
        else:
            # Mode 1 or 2: Specific region or auto-detect
            if region_chrom:
                log("Running in SPECIFIC REGION mode", "INFO")
            else:
                log("Running in AUTO-DETECT mode", "INFO")
            
            results = run_phasing_analysis(
                sam_file=Path(args.sam),
                output_dir=Path(args.outdir),
                prefix=args.prefix,
                region_chrom=region_chrom,
                region_start=region_start,
                region_end=region_end,
                size_min=args.size_min,
                size_max=args.size_max,
                upstream=args.upstream,
                downstream=args.downstream,
                use_multiplicity=args.multiplicity,
                r_script_path=Path(args.r_script) if args.r_script else None
            )
            
            if results:
                log("\nGenerated files:", "INFO")
                for key, path in results.items():
                    if key != 'region' and path:
                        log(f"  {key}: {path}", "INFO")
        
    except Exception as e:
        log(f"Analysis failed: {e}", "ERROR")
        import traceback
        traceback.print_exc()
        sys.exit(1)