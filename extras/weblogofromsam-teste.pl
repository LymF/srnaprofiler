#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

# --- Subroutine to check for required executables ---
sub check_dependency {
    my ($cmd) = @_;
    # A simple way to check if a command exists in the system's PATH
    return 1 if `which $cmd 2>/dev/null`;
    warn "WARNING: Command '$cmd' not found in your PATH. This script will likely fail.\n";
    return 0;
}

my $usage = "
Usage:
$0 -sam <sam_file> -p <prefix> -s <start_size> -e <end_size> -t <trim_length>
$0 -h

Description:
This script generates a WebLogo from sRNA sequences in a SAM file.
It splits reads by strand, filters by size, trims, and generates a logo.

Arguments:
  -sam <file>    : Input SAM file. (Required)
  -p <prefix>    : Prefix for all output files. (Required)
  -s <int>       : Minimum sRNA size to include (e.g., 21). (Required)
  -e <int>       : Maximum sRNA size to include (e.g., 24). (Required)
  -t <int>       : Length to trim all sequences to. (Required)
  -h             : Display this help message.

Example:
$0 -sam my_reads.sam -p my_locus -s 21 -e 24 -t 21
";

# --- Argument Handling ---
my ($sam, $prefix, $start, $end, $trim, $help);

GetOptions (
    "sam=s"    => \$sam,
    "prefix=s" => \$prefix, 
    "s=i"      => \$start, 
    "e=i"      => \$end, 
    "t=i"      => \$trim, 
    "h!"       => \$help
);

if ($help) {
    die $usage;
}

# --- Argument Validation ---
if (not defined $sam) {
    die "\nERROR: Please provide a SAM input file with -sam\n", $usage;
}
if (not defined $prefix) {
    die "\nERROR: Please provide an output prefix with -p\n", $usage;
}
if (not defined $start or $start <= 0) {
    die "\nERROR: Please provide a valid start size with -s\n", $usage;
}
if (not defined $end or $end < $start) {
    die "\nERROR: Please provide a valid end size with -e (must be >= start size)\n", $usage;
}
if (not defined $trim or $trim <= 0) {
    die "\nERROR: Please provide a valid trim length with -t\n", $usage;
}

# --- Dependency Checks ---
print "--- Checking for required programs ---\n";
check_dependency("samtools");
check_dependency("awk");
check_dependency("sort");
check_dependency("uniq");
check_dependency("fastx_trimmer");
check_dependency("weblogo");
print "NOTE: This script also depends on 'filter_fasta_by_size.py' being in the same directory or in your PATH.\n";
print "-------------------------------------\n\n";

# --- File Naming ---
my $size_range_tag = "${start}to${end}";
my $trim_tag = "${size_range_tag}.trim${trim}";

my $plus_sam_file = "$prefix.plus.sam";
my $minus_sam_file = "$prefix.minus.sam";
my $plus_fasta_file = "$prefix.byStrand.plus.fasta";
my $minus_fasta_file = "$prefix.byStrand.minus.fasta";
my $plus_filtered_fasta = "$prefix.byStrand.plus.$size_range_tag.fasta";
my $minus_filtered_fasta = "$prefix.byStrand.minus.$size_range_tag.fasta";
my $plus_trimmed_fasta = "$prefix.byStrand.plus.$trim_tag.fasta";
my $minus_trimmed_fasta = "$prefix.byStrand.minus.$trim_tag.fasta";
my $plus_uniq_fasta = "$prefix.byStrand.plus.$trim_tag.uniq.fasta";
my $minus_uniq_fasta = "$prefix.byStrand.minus.$trim_tag.uniq.fasta";
my $plus_pdf = "$prefix.byStrand.plus.$trim_tag.bits.weblogo.pdf";
my $minus_pdf = "$prefix.byStrand.minus.$trim_tag.bits.weblogo.pdf";

# --- Pipeline Execution ---

print "Step 1 of 7: Splitting SAM by strand...\n";
system("samtools view -F 0x10 $sam > $plus_sam_file") == 0 or die "ERROR: samtools failed to extract forward strand reads. $?";
system("samtools view -f 0x10 $sam > $minus_sam_file") == 0 or die "ERROR: samtools failed to extract reverse strand reads. $?";

print "Step 2 of 7: Converting SAM to FASTA...\n";
system("awk '/^@/{ next; } { print \">\"\$1; print \$10 }' $plus_sam_file > $plus_fasta_file") == 0 or die "ERROR: awk failed to convert plus strand SAM to FASTA. $?";
system("awk '/^@/{ next; } { print \">\"\$1; print \$10 }' $minus_sam_file > $minus_fasta_file") == 0 or die "ERROR: awk failed to convert minus strand SAM to FASTA. $?";

print "Step 3 of 7: Filtering FASTA by size (>= $start, <= $end)...\n";
system("python3 filter_fasta_by_size.py $minus_fasta_file $start $end $minus_filtered_fasta") == 0 or die "ERROR: Python script 'filter_fasta_by_size.py' failed for minus strand. $?";
system("python3 filter_fasta_by_size.py $plus_fasta_file $start $end $plus_filtered_fasta") == 0 or die "ERROR: Python script 'filter_fasta_by_size.py' failed for plus strand. $?";

print "Step 4 of 7: Trimming sequences to $trim nt...\n";
system("fastx_trimmer -l $trim -i $minus_filtered_fasta -o $minus_trimmed_fasta") == 0 or die "ERROR: fastx_trimmer failed for minus strand. $?";
system("fastx_trimmer -l $trim -i $plus_filtered_fasta -o $plus_trimmed_fasta") == 0 or die "ERROR: fastx_trimmer failed for plus strand. $?";

print "Step 5 of 7: Converting T to U for RNA WebLogo...\n";
system("perl -pi -e 's/T/U/gi' $minus_trimmed_fasta") == 0 or die "ERROR: Failed to convert T to U in minus strand file. $?";
system("perl -pi -e 's/T/U/gi' $plus_trimmed_fasta") == 0 or die "ERROR: Failed to convert T to U in plus strand file. $?";

print "Step 6 of 7: Getting unique reads...\n";
# --- CHANGE APPLIED HERE: Replaced non-standard 'nsort -u' with standard 'sort | uniq' ---
system("cat $minus_trimmed_fasta | grep -v '>' | sort  | uniq | awk '{ print \">seq_\" NR; print \$1 }' > $minus_uniq_fasta") == 0 or die "ERROR: Failed to get unique reads for minus strand. $?";
system("cat $plus_trimmed_fasta | grep -v '>' | sort  | uniq | awk '{ print \">seq_\" NR; print \$1 }' > $plus_uniq_fasta") == 0 or die "ERROR: Failed to get unique reads for plus strand. $?";

print "Step 7 of 7: Plotting WebLogo...\n";
my $weblogo_opts = "-F pdf -A rna -U bits -s large --resolution 300 --errorbars NO -x 'Read coordinates (nt)' -S 1 --ticmarks 0.5";
system("weblogo $weblogo_opts -f $minus_uniq_fasta -c classic -P 'antisense' -o $minus_pdf") == 0 or die "ERROR: weblogo command failed for minus strand. $?";
system("weblogo $weblogo_opts -f $plus_uniq_fasta -c classic -P 'sense' -o $plus_pdf") == 0 or die "ERROR: weblogo command failed for plus strand. $?";

print "\n--- Pipeline finished successfully! ---\n";
print "Cleaning up intermediate files...\n";
`rm -f $plus_sam_file $minus_sam_file $plus_fasta_file $minus_fasta_file $plus_filtered_fasta $minus_filtered_fasta $plus_trimmed_fasta $minus_trimmed_fasta $plus_uniq_fasta $minus_uniq_fasta`;
print "Cleanup complete.\n";
