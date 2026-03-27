#!/usr/bin/perl5.32.1 

use strict;
use warnings;
use Getopt::Long;

my $usage = "
Usage:
$0 -s <sam file> -p <prefix> -f <output folder> -c <chromosome> -r <region> -si <start size> -se <end size> [--combined]
$0 -h

Description:
This script performs phasing and offset analysis on sRNA SAM files.
It can run in two modes:

1. INDIVIDUAL (default): Iterates through each size from <start size> to <end size>,
   performing a complete, separate analysis for each size.

2. COMBINED (--combined): Performs a single analysis on all reads within the
   entire size range specified by -si and -se.

Arguments:
  -s          : Input SAM file. (Required)
  -p          : Prefix for output files. (Required)
  -f          : Folder to save the results. (Required)
  -c          : Chromosome or contig name. (Required)
  -r          : Region in 'start-end' format (e.g., 1-2000). (Required)
  -si         : Start of the size range to analyze (inclusive). (Required)
  -se         : End of the size range to analyze (inclusive). (Required)
  --combined  : Optional flag to run a single analysis for the entire size range.
  -h          : Display this help message.

Examples:
# Individual analysis for 21, 22, and 23 nt sRNAs
$0 -s r.sam -p l1 -f res -c chr1 -r 1-2000 -si 21 -se 23

# Combined analysis for all sRNAs from 21 to 23 nt
$0 -s r.sam -p l1 -f res -c chr1 -r 1-2000 -si 21 -se 23 --combined
";

$| = 1;     # Force immediate flushing of output, no buffering.

my $sam;
my $prefix;
my $folder;
my $chr;
my $region;
my $si;
my $se;
my $help;
# --- NEW ARGUMENT VARIABLE ---
my $combined_mode;

GetOptions (
    "s=s"  => \$sam,
    "p=s"  => \$prefix,
    "f=s"  => \$folder,
    "r=s"  => \$region,
    "c=s"  => \$chr,
    "si=i" => \$si,
    "se=i" => \$se, 
    "h!"   => \$help,
    # --- NEW ARGUMENT ---
    "combined!" => \$combined_mode
);

if ($help) {
    die $usage;
}

# --- Argument Validation ---
# (Omitted for brevity, same as before)
die "\nERROR: Please provide a SAM input file with -s.", $usage if not defined $sam;
die "\nERROR: Please provide an output prefix with -p.", $usage if not defined $prefix;
die "\nERROR: Please provide an output folder with -f.", $usage if not defined $folder;
die "\nERROR: Please provide a chromosome name with -c.", $usage if not defined $chr;
die "\nERROR: Please provide a region in 'start-end' format with -r.", $usage if not defined $region;
die "\nERROR: Invalid value for start size (-si). Must be >= 0.", $usage if not defined $si or $si < 0;
die "\nERROR: Invalid value for end size (-se). Must be >= start size.", $usage if not defined $se or $se < 0 or $se < $si;


# Check if the output directory exists. If not, create it.
if (not -d $folder) {
    print "Creating output directory: $folder\n";
    mkdir $folder or die "ERROR: Could not create directory '$folder': $!";
}

# --- Path to external scripts ---
my $split_script = "/home/ericgdp/bin/splitSamByStrand.pl";
my $filter_script = "/home/ericgdp/bin/filterSamBowtieBySize.pl";
my $offset_script = "/home/ericgdp/bin/bed_calculate_offsets.pl";
my $plot_script = "/home/ericgdp/bin/R/phasing-plot-teste.R";

# Check if the R script exists before proceeding
if (not -f $plot_script) {
    die "ERROR: R plotting script not found at: $plot_script\n";
}

# --- Pre-computation Step (run only once for both modes) ---
print "1. Splitting SAM by strand (this is done only once)...\n";
my $plus_sam = "$folder/$prefix.plus.sam";
my $minus_sam = "$folder/$prefix.minus.sam";
system("perl5.34.0 $split_script -in $sam -out $folder/$prefix") == 0 or die "Failed to execute split: $?";


# --- NEW: CHOOSE ANALYSIS MODE BASED ON --combined FLAG ---
if ($combined_mode) {
    # --- MODE 1: COMBINED ANALYSIS (for a range of sizes) ---
    my $size_tag = "$si-$se";
    print "\n" . "="x60 . "\n";
    print "=== Running in COMBINED mode for range: $size_tag nt\n";
    print "="x60 . "\n\n";

    run_analysis($size_tag, $si, $se);

} else {
    # --- MODE 2: INDIVIDUAL ANALYSIS (loop for each size) ---
    print "\n" . "="x60 . "\n";
    print "=== Running in INDIVIDUAL mode for each size from $si to $se nt\n";
    print "="x60 . "\n\n";

    for my $current_size ($si .. $se) {
        print "--- Processing size: $current_size nt ---\n";
        run_analysis($current_size, $current_size, $current_size);
        print "--- Finished processing size $current_size nt ---\n\n";
    }
}

print "\n" . "="x60 . "\n";
print "=== All analyses completed successfully! ===\n";
print "="x60 . "\n\n";


# --- NEW: SUBROUTINE TO AVOID CODE DUPLICATION ---
# This subroutine contains the core analysis pipeline.
# It's called differently depending on the chosen mode.
sub run_analysis {
    my ($tag, $start_size, $end_size) = @_;

    # Define filenames using the tag (e.g., "21" or "21-24")
    my $plus_sized_sam = "$folder/$prefix.plus.$tag.sam";
    my $minus_sized_sam = "$folder/$prefix.minus.$tag.sam";
    my $total_sized_sam = "$folder/$prefix.$tag.sam";
    my $plus_bed_file = "$folder/$prefix.plus.$tag.bed";
    my $minus_bed_file = "$folder/$prefix.minus.$tag.bed";
    my $total_bed_file = "$folder/$prefix.$tag.bed";
    my $plus_phasing = "$folder/$prefix.plus.$tag.phasing";
    my $minus_phasing = "$folder/$prefix.minus.$tag.phasing";
    my $offset_file = "$folder/$prefix.$tag.offset";

    # Step 2: Filter SAM files
    print "Step 2: Filtering SAM for size(s) $tag nt...\n";
    system("perl5.34.0 $filter_script -i $plus_sam -si $start_size -se $end_size -o $plus_sized_sam");
    system("perl5.34.0 $filter_script -i $minus_sam -si $start_size -se $end_size -o $minus_sized_sam");
    system("perl5.34.0 $filter_script -i $sam -si $start_size -se $end_size -o $total_sized_sam");

    # Step 3: Convert SAM to BED
    print "Step 3: Converting SAM to BED...\n";
    system("samtools view -bS $plus_sized_sam | bedtools bamtobed -i - > $plus_bed_file");
    system("samtools view -bS $minus_sized_sam | bedtools bamtobed -i - > $minus_bed_file");
    system("samtools view -bS $total_sized_sam | bedtools bamtobed -i - > $total_bed_file");

    # Step 4: Calculate Phasing and Offset
    print "Step 4: Calculating phasing and offset...\n";
    system("perl5.34.0 $offset_script --region $chr:$region --up 80 --down 80 --header $plus_bed_file | grep -v \"intervals\" | grep -v \"chrom\"  | sed 1d | sort -k5,5 | groupBy -g 5 -c 7 -o sum > $plus_phasing"); 
    system("perl5.34.0 $offset_script --region $chr:$region --up 80 --down 80 --header $minus_bed_file | grep -v \"intervals\" | grep -v \"chrom\"  | sed 1d | sort -k5,5 | groupBy -g 5 -c 7 -o sum > $minus_phasing"); 
    system("perl5.34.0 $offset_script --region $chr:$region --up 80 --down 80 --header $total_bed_file | grep -v \"intervals\" | grep -v \"chrom\"  | sed 1d | sort -k5,5 | groupBy -g 5 -c 9 -o sum > $offset_file");

    # Step 5: Generate plots with R
    print "Step 5: Generating plots with R...\n";
    my $output_prefix_size = "$folder/$prefix.$tag"; 
    
    if (-s $plus_phasing && -s $minus_phasing && -s $offset_file) {
        my $cmd_plot = "Rscript $plot_script $plus_phasing $minus_phasing $offset_file $output_prefix_size";
        print "R command: $cmd_plot\n";
        system($cmd_plot);
        
        if (-e "$output_prefix_size.phasing_offset_line.pdf") {
            print "Successfully created PDF for size(s) $tag nt.\n";
        } else {
            print "ERROR: R script ran but failed to create PDF for size(s) $tag nt. Check R errors above.\n";
        }
    } else {
        print "WARNING: Skipping plot for size(s) $tag nt because one or more input files are empty or missing.\n";
    }
}
