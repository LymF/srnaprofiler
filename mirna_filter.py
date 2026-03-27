#!/usr/bin/env python3
"""
microRNA Filtering Module for Small RNA Profiling Pipeline

This module implements miRNA detection and filtering to replicate the behavior
of the original Perl pipeline. It detects miRNA clusters based on:
1. Overlap with known miRNA loci (from GFF, BED or FASTA)
2. Size distribution pattern (22-24 nt, strand-biased)
3. Abundance thresholds

Author: Adapted from original Perl pipeline
"""

import os
import logging
from pathlib import Path
from typing import Dict, Set, Tuple, Optional, List
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class MiRNAFilter:
    """
    Handles microRNA detection and filtering for small RNA clusters.
    
    Attributes:
        mirna_fasta: Path to FASTA file with known miRNA sequences
        mirna_gff: Path to GFF file with miRNA loci coordinates
        mirna_bed: Path to BED file with miRNA loci coordinates
        keep_mirnas: If True, keep miRNAs in output but label them; if False, exclude them
        logger: Logging instance
    """
    
    def __init__(
        self,
        mirna_fasta: Optional[str] = None,
        mirna_gff: Optional[str] = None,
        mirna_bed: Optional[str] = None,
        keep_mirnas: bool = False,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the miRNA filter.
        
        Args:
            mirna_fasta: Path to miRNA FASTA reference (optional)
            mirna_gff: Path to miRNA GFF reference (optional)
            mirna_bed: Path to miRNA BED reference (optional)
            keep_mirnas: Whether to keep miRNAs in analysis (default: False)
            logger: Logger instance (creates new if None)
        """
        self.mirna_fasta = mirna_fasta
        self.mirna_gff = mirna_gff
        self.mirna_bed = mirna_bed
        self.keep_mirnas = keep_mirnas
        self.logger = logger or self._setup_logger()
        
        # Storage for miRNA references
        self.mirna_ids: Set[str] = set()
        self.mirna_regions: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        
        # Load reference data
        self._load_references()
    
    def _setup_logger(self) -> logging.Logger:
        """Create a logger if none provided."""
        logger = logging.getLogger('MiRNAFilter')
        logger.setLevel(logging.INFO)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        return logger
    
    def _load_references(self):
        """Load miRNA reference data from FASTA, GFF, or BED."""
        if self.mirna_fasta:
            self._load_mirna_fasta()
        
        if self.mirna_gff:
            self._load_mirna_gff()
        
        if self.mirna_bed:
            self._load_mirna_bed()
        
        if not self.mirna_fasta and not self.mirna_gff and not self.mirna_bed:
            self.logger.info("No miRNA reference provided - skipping miRNA filtering")
    
    def _load_mirna_fasta(self):
        """Load miRNA IDs from FASTA file."""
        self.logger.info(f"Loading miRNA sequences from {self.mirna_fasta}")
        
        try:
            for record in SeqIO.parse(self.mirna_fasta, "fasta"):
                # Store the sequence ID (chromosome/contig name)
                self.mirna_ids.add(record.id)
            
            self.logger.info(f"Loaded {len(self.mirna_ids)} miRNA IDs from FASTA")
        
        except Exception as e:
            self.logger.error(f"Error loading miRNA FASTA: {e}")
            raise
    
    def _load_mirna_gff(self):
        """Load miRNA genomic coordinates from GFF file."""
        self.logger.info(f"Loading miRNA loci from {self.mirna_gff}")
        
        try:
            with open(self.mirna_gff, 'r') as f:
                for line in f:
                    # Skip comments and headers
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    
                    # GFF format: seqid, source, type, start, end, score, strand, phase, attributes
                    seqid = fields[0]
                    feature_type = fields[2]
                    start = int(fields[3])
                    end = int(fields[4])
                    
                    # Filter for miRNA features (adjust as needed for your GFF)
                    if 'miRNA' in feature_type or 'miR' in feature_type:
                        self.mirna_regions[seqid].append((start, end))
            
            total_regions = sum(len(regions) for regions in self.mirna_regions.values())
            self.logger.info(
                f"Loaded {total_regions} miRNA regions across "
                f"{len(self.mirna_regions)} chromosomes from GFF"
            )
        
        except Exception as e:
            self.logger.error(f"Error loading miRNA GFF: {e}")
            raise

    def _load_mirna_bed(self):
        """Load miRNA genomic coordinates from BED file."""
        self.logger.info(f"Loading miRNA loci from {self.mirna_bed}")
        
        try:
            with open(self.mirna_bed, 'r') as f:
                for line in f:
                    # Skip comments and headers
                    if line.startswith('#') or line.startswith('track') or line.startswith('browser') or not line.strip():
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 3:
                        continue
                    
                    # BED format: chrom, start, end, name, score, strand, ...
                    seqid = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    
                    # Store the region
                    self.mirna_regions[seqid].append((start, end))
            
            total_regions = sum(len(regions) for regions in self.mirna_regions.values())
            self.logger.info(
                f"Loaded {total_regions} miRNA regions across "
                f"{len(self.mirna_regions)} chromosomes from BED"
            )
        
        except Exception as e:
            self.logger.error(f"Error loading miRNA BED: {e}")
            raise
    
    def check_mirna_overlap(
        self,
        cluster_id: str,
        chromosome: str,
        start: int,
        end: int
    ) -> bool:
        """
        Check if a cluster overlaps with known miRNA loci.
        
        Args:
            cluster_id: Cluster identifier
            chromosome: Chromosome/contig name
            start: Cluster start position
            end: Cluster end position
        
        Returns:
            True if cluster overlaps a miRNA locus, False otherwise
        """
        # Check by ID match (FASTA-based)
        if cluster_id in self.mirna_ids or chromosome in self.mirna_ids:
            return True
        
        # Check by genomic overlap (GFF/BED-based)
        if chromosome in self.mirna_regions:
            for mirna_start, mirna_end in self.mirna_regions[chromosome]:
                # Check for any overlap
                if not (end < mirna_start or start > mirna_end):
                    return True
        
        return False
    
    def classify_mirna_by_pattern(
        self,
        size_distribution: Dict[int, int],
        size_distribution_pos: Dict[int, int],
        size_distribution_neg: Dict[int, int],
        total_reads: int
    ) -> bool:
        """
        Classify a cluster as miRNA based on size distribution pattern.
        
        This replicates the Perl logic:
        - miRNAs: predominantly 22-24 nt
        - Strong strand bias (ratio > 1000)
        - Sufficient abundance (> 50 reads)
        - Dominates over degradation products (15-19 nt)
        
        Args:
            size_distribution: Total reads per size
            size_distribution_pos: Positive strand reads per size
            size_distribution_neg: Negative strand reads per size
            total_reads: Total cluster reads
        
        Returns:
            True if pattern matches miRNA characteristics
        """
        # Count reads in miRNA size range (22-24 nt)
        mirna_size_total = 0
        mirna_size_pos = 0
        mirna_size_neg = 0
        
        for size in range(22, 25):  # 22, 23, 24
            mirna_size_total += size_distribution.get(size, 0)
            mirna_size_pos += size_distribution_pos.get(size, 0)
            mirna_size_neg += size_distribution_neg.get(size, 0)
        
        # Count degradation reads (15-19 nt)
        degradation_total = sum(
            size_distribution.get(size, 0) for size in range(15, 20)
        )
        
        # Count siRNA reads (20-22 nt)
        sirna_total = sum(
            size_distribution.get(size, 0) for size in range(20, 23)
        )
        
        # Count piRNA reads (24-29 nt)
        pirna_total = sum(
            size_distribution.get(size, 0) for size in range(24, 30)
        )
        
        # Check miRNA criteria (from Perl code)
        if mirna_size_total > 50:  # Minimum abundance
            # Must dominate over degradation, siRNA, and piRNA
            if (mirna_size_total > degradation_total and 
                mirna_size_total > sirna_total and
                mirna_size_total > pirna_total):
                
                # Check for strong strand bias
                pos = max(mirna_size_pos, 1)  # Avoid division by zero
                neg = max(mirna_size_neg, 1)
                
                strand_ratio = pos / neg if pos > neg else neg / pos
                
                # Very strong strand bias indicates miRNA (>1000:1 ratio)
                if strand_ratio > 1000:
                    return True
        
        return False
    
    def filter_clusters(
        self,
        clusters_fasta: str,
        classification_df: pd.DataFrame,
        size_stats: Dict[str, Dict],
        output_prefix: str
    ) -> Tuple[pd.DataFrame, Set[str]]:
        """
        Main filtering function to identify and handle miRNA clusters.
        
        Args:
            clusters_fasta: Path to input clusters FASTA
            classification_df: DataFrame with cluster classifications
            size_stats: Dictionary with size distribution per cluster
            output_prefix: Output file prefix
        
        Returns:
            Tuple of (updated_classification_df, set_of_mirna_ids)
        """
        if not self.mirna_fasta and not self.mirna_gff and not self.mirna_bed:
            self.logger.info("Skipping miRNA filtering (no reference provided)")
            return classification_df, set()
        
        self.logger.info("Starting miRNA cluster identification...")
        
        mirna_clusters = set()
        mirna_records = []
        
        # Parse cluster coordinates and sequences
        cluster_info = {}
        for record in SeqIO.parse(clusters_fasta, "fasta"):
            cluster_info[record.id] = record
        
        # Check each cluster
        for cluster_id, record in cluster_info.items():
            is_mirna = False
            
            # Extract chromosome and coordinates from cluster ID
            # Expected format: chr:start-end or similar
            # Adjust parsing based on your actual format
            chrom, coords = self._parse_cluster_id(cluster_id)
            
            if chrom and coords:
                start, end = coords
                
                # Check overlap with known miRNA loci
                if self.check_mirna_overlap(cluster_id, chrom, start, end):
                    is_mirna = True
                    self.logger.debug(f"Cluster {cluster_id} overlaps known miRNA locus")
            
            # Also check by size distribution pattern
            if cluster_id in size_stats:
                stats = size_stats[cluster_id]
                
                if self.classify_mirna_by_pattern(
                    stats.get('size_dist', {}),
                    stats.get('size_dist_pos', {}),
                    stats.get('size_dist_neg', {}),
                    stats.get('total_reads', 0)
                ):
                    is_mirna = True
                    self.logger.debug(f"Cluster {cluster_id} matches miRNA size pattern")
            
            if is_mirna:
                mirna_clusters.add(cluster_id)
                mirna_records.append(record)
        
        self.logger.info(f"Identified {len(mirna_clusters)} miRNA clusters")
        
        # Write miRNA FASTA
        if mirna_records:
            mirna_fasta_path = f"{output_prefix}.clusters.miRNA.fasta"
            SeqIO.write(mirna_records, mirna_fasta_path, "fasta")
            self.logger.info(f"Wrote miRNA clusters to {mirna_fasta_path}")
        
        # Update classification DataFrame
        classification_df = self._update_classifications(
            classification_df,
            mirna_clusters
        )
        
        return classification_df, mirna_clusters
    
    def _parse_cluster_id(self, cluster_id: str) -> Tuple[Optional[str], Optional[Tuple[int, int]]]:
        """
        Parse cluster ID to extract chromosome and coordinates.
        
        Args:
            cluster_id: Cluster identifier (e.g., "chr1:1000-2000")
        
        Returns:
            Tuple of (chromosome, (start, end)) or (None, None) if parsing fails
        """
        try:
            # Try common formats
            if ':' in cluster_id and '-' in cluster_id:
                # Format: chr:start-end
                chrom, coords = cluster_id.split(':')
                start, end = coords.split('-')
                return chrom, (int(start), int(end))
            
            elif '_' in cluster_id:
                # Format: chr_start_end
                parts = cluster_id.split('_')
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[-2])
                    end = int(parts[-1])
                    return chrom, (start, end)
        
        except (ValueError, IndexError):
            pass
        
        return None, None
    
    def _update_classifications(
        self,
        classification_df: pd.DataFrame,
        mirna_clusters: Set[str]
    ) -> pd.DataFrame:
        """
        Update classification DataFrame with miRNA labels.
        
        Args:
            classification_df: Original classification DataFrame
            mirna_clusters: Set of miRNA cluster IDs
        
        Returns:
            Updated DataFrame
        """
        # Create a copy to avoid modifying original
        df = classification_df.copy()
        
        # Add miRNA classification
        df.loc[df['cluster_id'].isin(mirna_clusters), 'final_classification'] = 'miRNA'
        
        # If not keeping miRNAs, mark them for exclusion
        if not self.keep_mirnas:
            df.loc[df['cluster_id'].isin(mirna_clusters), 'exclude_from_analysis'] = True
            self.logger.info(
                f"Marked {len(mirna_clusters)} miRNA clusters for exclusion "
                "(--mirna-keep not specified)"
            )
        else:
            self.logger.info(
                f"Keeping {len(mirna_clusters)} miRNA clusters in analysis "
                "(--mirna-keep specified)"
            )
        
        return df
    
    def write_filtered_fastas(
        self,
        input_fasta: str,
        output_prefix: str,
        mirna_clusters: Set[str],
        classification_df: pd.DataFrame
    ):
        """
        Write separate FASTA files for each RNA class, excluding miRNAs if needed.
        
        Args:
            input_fasta: Original clusters FASTA
            output_prefix: Output file prefix
            mirna_clusters: Set of miRNA cluster IDs
            classification_df: Classification DataFrame
        """
        # Read all sequences
        all_records = {rec.id: rec for rec in SeqIO.parse(input_fasta, "fasta")}
        
        # Separate by classification
        sirna_records = []
        pirna_records = []
        ambiguous_records = []
        mirna_records = []
        
        for _, row in classification_df.iterrows():
            cluster_id = row['cluster_id']
            classification = row['final_classification']
            
            if cluster_id not in all_records:
                continue
            
            record = all_records[cluster_id]
            
            if classification == 'miRNA':
                mirna_records.append(record)
            elif classification == 'putative_siRNA':
                # Exclude if miRNA and not keeping
                if not (cluster_id in mirna_clusters and not self.keep_mirnas):
                    sirna_records.append(record)
            elif classification == 'putative_piRNA':
                if not (cluster_id in mirna_clusters and not self.keep_mirnas):
                    pirna_records.append(record)
            elif classification == 'ambiguous':
                if not (cluster_id in mirna_clusters and not self.keep_mirnas):
                    ambiguous_records.append(record)
        
        # Write files
        if sirna_records:
            SeqIO.write(sirna_records, f"{output_prefix}.putative_siRNA.fasta", "fasta")
            self.logger.info(f"Wrote {len(sirna_records)} siRNA clusters")
        
        if pirna_records:
            SeqIO.write(pirna_records, f"{output_prefix}.putative_piRNA.fasta", "fasta")
            self.logger.info(f"Wrote {len(pirna_records)} piRNA clusters")
        
        if ambiguous_records:
            SeqIO.write(ambiguous_records, f"{output_prefix}.ambiguous.fasta", "fasta")
            self.logger.info(f"Wrote {len(ambiguous_records)} ambiguous clusters")
        
        if mirna_records:
            SeqIO.write(mirna_records, f"{output_prefix}.miRNA.fasta", "fasta")
            self.logger.info(f"Wrote {len(mirna_records)} miRNA clusters")


def step_filter_mirna(
    clusters_fasta: str,
    classification_table: str,
    size_stats_dict: Dict,
    output_prefix: str,
    mirna_fasta: Optional[str] = None,
    mirna_gff: Optional[str] = None,
    mirna_bed: Optional[str] = None,
    mirna_keep: bool = False,
    logger: Optional[logging.Logger] = None
) -> Tuple[str, Set[str]]:
    """
    Pipeline integration function for miRNA filtering.
    
    This function should be called in the main pipeline after cluster identification
    but before siRNA/piRNA classification and pattern analysis.
    
    Args:
        clusters_fasta: Path to clusters FASTA file
        classification_table: Path to classification TSV file
        size_stats_dict: Dictionary with size statistics per cluster
        output_prefix: Output file prefix
        mirna_fasta: Path to miRNA FASTA reference (optional)
        mirna_gff: Path to miRNA GFF reference (optional)
        mirna_bed: Path to miRNA BED reference (optional)
        mirna_keep: Whether to keep miRNAs in analysis
        logger: Logger instance
    
    Returns:
        Tuple of (updated_classification_table_path, set_of_mirna_cluster_ids)
    
    Example integration in main pipeline:
        ```python
        # After cluster identification, before siRNA/piRNA classification
        updated_classification, mirna_ids = step_filter_mirna(
            clusters_fasta=f"{prefix}.clusters.fasta",
            classification_table=f"{prefix}.classification.tsv",
            size_stats_dict=size_statistics,
            output_prefix=prefix,
            mirna_fasta=args.mirna_fasta,
            mirna_gff=args.mirna_gff,
            mirna_bed=args.mirna_bed,
            mirna_keep=args.mirna_keep,
            logger=logger
        )
        ```
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Initialize filter
    mirna_filter = MiRNAFilter(
        mirna_fasta=mirna_fasta,
        mirna_gff=mirna_gff,
        mirna_bed=mirna_bed,
        keep_mirnas=mirna_keep,
        logger=logger
    )
    
    # Load classification table
    try:
        classification_df = pd.read_csv(classification_table, sep='\t')
    except Exception as e:
        logger.error(f"Error loading classification table: {e}")
        raise
    
    # Filter clusters
    updated_df, mirna_clusters = mirna_filter.filter_clusters(
        clusters_fasta=clusters_fasta,
        classification_df=classification_df,
        size_stats=size_stats_dict,
        output_prefix=output_prefix
    )
    
    # Write updated classification table
    updated_table_path = f"{output_prefix}.classification.miRNA_filtered.tsv"
    updated_df.to_csv(updated_table_path, sep='\t', index=False)
    logger.info(f"Wrote updated classification table to {updated_table_path}")
    
    # Write filtered FASTA files
    mirna_filter.write_filtered_fastas(
        input_fasta=clusters_fasta,
        output_prefix=output_prefix,
        mirna_clusters=mirna_clusters,
        classification_df=updated_df
    )
    
    return updated_table_path, mirna_clusters


if __name__ == "__main__":
    # Example usage for testing
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Filter miRNA clusters from small RNA profiling results"
    )
    parser.add_argument(
        "--clusters-fasta",
        required=True,
        help="Input clusters FASTA file"
    )
    parser.add_argument(
        "--classification-table",
        required=True,
        help="Input classification table (TSV)"
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output file prefix"
    )
    parser.add_argument(
        "--mirna-fasta",
        help="miRNA reference FASTA file"
    )
    parser.add_argument(
        "--mirna-gff",
        help="miRNA reference GFF file"
    )
    parser.add_argument(
        "--mirna-bed",
        help="miRNA reference BED file"
    )
    parser.add_argument(
        "--mirna-keep",
        action="store_true",
        help="Keep miRNAs in analysis (label only, don't exclude)"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Run filtering (with empty size_stats_dict for testing)
    step_filter_mirna(
        clusters_fasta=args.clusters_fasta,
        classification_table=args.classification_table,
        size_stats_dict={},
        output_prefix=args.output_prefix,
        mirna_fasta=args.mirna_fasta,
        mirna_gff=args.mirna_gff,
        mirna_bed=args.mirna_bed,
        mirna_keep=args.mirna_keep
    )