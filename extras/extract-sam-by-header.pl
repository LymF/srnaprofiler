#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
Usage:
$0 -i <input.sam> -c <chromosome_name> -o <output.sam>
$0 -h

Description:
This script efficiently extracts all header lines and all alignment records
for a specific chromosome/contig from a large SAM file into a new SAM file.
The original input file is not modified.

Arguments:
  -i, --input    : Path to the large input SAM file. (Required)
  -c, --chrom    : Name of the chromosome/contig to extract. (Required)
  -o, --output   : Path for the new, smaller output SAM file. (Required)
  -h, --help     : Display this help message.

Example:
# Extracts all data for 'Ac2012_v1_629_kmer_6001-8000' from 'finalSIRNA.sam'
# and saves it to 'kmer629_only.sam'.

$0 -i finalSIRNA.sam -c Ac2012_v1_629_kmer_6001-8000 -o kmer629_only.sam
";

my ($input_sam, $output_sam, $target_chrom, $help);

GetOptions(
    "i|input=s"  => \$input_sam,
    "o|output=s" => \$output_sam,
    "c|chrom=s"  => \$target_chrom,
    "h|help"     => \$help
);

if ($help || !$input_sam || !$output_sam || !$target_chrom) {
    die $usage;
}

# Open input and output files, with error handling
open my $in_fh, '<', $input_sam or die "Error: Cannot open input file '$input_sam': $!\n";
open my $out_fh, '>', $output_sam or die "Error: Cannot open output file '$output_sam': $!\n";

print "Extracting chromosome '$target_chrom' from '$input_sam'...\n";

my $found_alignments = 0;

# Read the input SAM file line by line
while (my $line = <$in_fh>) {
    # If it's a header line, print it directly to the output
    if ($line =~ /^@/) {
        print $out_fh $line;
        next; # Move to the next line
    }

    # If it's an alignment line, split it by tabs
    my @fields = split /\t/, $line;

    # The chromosome/RNAME is in the 3rd column (index 2)
    # Check if this line belongs to our target chromosome
    if (defined $fields[2] && $fields[2] eq $target_chrom) {
        print $out_fh $line;
        $found_alignments++;
    }
}

# Close the filehandles
close $in_fh;
close $out_fh;

print "Done.\n";
print "Found and wrote $found_alignments alignment records to '$output_sam'.\n";

