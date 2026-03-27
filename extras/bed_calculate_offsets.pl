#!/usr/bin/env perl
=pod
  piRNA 5" starting position counter

  Copyright (C) 2011 A. Gordon (gordon at cshl dot edu)
                 and F. Muerdter (fmuerdte at cshl dot edu)

  License: Free for personal / academic / non-commercial use.

  See related publication:
     [ Production of artificial piRNAs in flies and mice ]
     [ Muerdter, F.*, Olovnikov, I.*, Molaro, A.*, Rozhkov, N.V.*, Czech, B., Gordon, A., Hannon, G.J., Aravin, A.A. ]
     [ RNA (2012) ]
     [ doi: 10.1261/rna.029769.111 ]

  For questions, please contact the corresponding authors:
    A. A. Aravin, E-mail aaa@caltech.edu.
	G. J. Hannon, E-mail hannon@cshl.edu.

  When using this script (or derivatives there-of),
  please cite: [ Muerdter et al., 2012 ]
=cut

##
## Version 0.2
##
use strict;
use warnings;
use Data::Dumper;
use Carp;
use Getopt::Long;



##
## Forward function declarations
##
sub parse_command_line();
sub show_help();
sub set_filter_region($);
sub load_external_intervals();
sub get_next_interval();
sub interval_in_region($$$);
sub initialize_empty_windows();
sub calculate_region_coverage();
sub scan_region_windows();
sub sort_intervals_to_process();

##
## Global Parameters
##   Set by user on command line
my $region_chrom;	# Process reads from this region only, ignore all other reads
my $region_start;
my $region_end;
my $use_multiplity_count; # Use multiplicity count from read name, column 4 (N-M) where M is the multiplicity count.
my $extend_upstream=10;
my $extend_downstream=10;
my $header; # print header line
my $debug; # print a lot of debugging information to STDERR
my $external_intervals; #if set, we calculate the coverage from one BED file, but scan the intervals from another BED file

## Global variables
my $region_size;
my $stats_lines_count=0;
my $stats_sequences_count=0;
my $stats_reads_count=0;

# Sequences = multiplitiy of 1 (one per line/interval)
# Reads = with multiplicity values
my @window_positive_sequences;
my @window_negative_sequences;
my @window_positive_reads;
my @window_negative_reads;
my @intervals_to_process;

##
## Program start
##

parse_command_line();
initialize_empty_windows();
print STDERR "Loading intervals..." if $debug;
calculate_region_coverage();
if ($external_intervals) {
	print STDERR "Loading external intervals..." if $debug;
	load_external_intervals();
}
print STDERR "Sorting intervals..." if $debug;
sort_intervals_to_process();
if ($header) {
	print join("\t", qw/chrom genomic_start strand window_position offset seq_count_same_strand read_count_same_strand seq_count_opposite_strand read_count_opposite_strand/), "\n";
}
scan_region_windows();

if ($debug) {
	print STDERR "Read $stats_lines_count lines\n";
	print STDERR "Processed $stats_sequences_count intervals ($stats_reads_count reads) in region $region_chrom:$region_start-$region_end\n";
}



#chrX:21,539,852-21,540,175

##
## Program End
##

sub show_help()
{
	print<<EOF;
piRNA 5" starting position counter
Copyright (C) 2011 A. Gordon   (gordonatcshl.edu)
               and F. Muerdter (fmuerdteatcshl.edu)

  License: Free for personal / academic / non-commercial use.

  See related publication:
     [ Production of artificial piRNAs in flies and mice ]
     [ Muerdter, F.*, Olovnikov, I.*, Molaro, A.*, Rozhkov, N.V.*, Czech, B., Gordon, A., Hannon, G.J., Aravin, A.A. ]
     [ RNA (2012) ]
     [ doi: 10.1261/rna.029769.111 ]

  For questions, please contact the corresponding authors:
    A. A. Aravin, E-mail aaaat.edu.
	G. J. Hannon, E-mail hannonat.edu.

  When using this script (or derivatives there-of),
  please cite: [ Muerdter et al., 2012 ]

Usage: $0 [OPTIONS] INPUT.BED

Options:
 -h   --help    =    This help screen.

 -r chrX:start-end
 --region chrX:start-end  =  Specify region of interest. Reads outside this region will be ignored.

 --up N         = Start window N upstream of the intervals 5" position.
 --down N       = End window N downstream of the intervals 5" position.

 -m
 --mult         = Take multiplicity count from name column (4th column).
                  Multiplicity should be in "N-M" format, where M is the multiplicity value.

 --intervals FILE.BED = Compare the intervals in FILE.BED against the coverage of INPUT.BED.
                        By default, the intervals in INPUT.BED are compared against their own coverage.

 -h
 --header       = Print header line.

 --debug        = Print debug information

Example:

\$ cat test.bed
chr1	95	125	1-7	0	+
chr1	100	150	1-13	0	+
chr1	100	140	1-4	0	+
chr1	80	111	1-20	0	-

\$ $0 --region chr1:0-2000 --mult --up 10 --down 10 --header test.bed
chrom	genomic_start	strand	window_position	offset	seq_count_same_strand	read_count_same_strand	seq_count_opposite_strand	read_count_opposite_strand
chr1	95	+	95	5	2	17	0	0
chr1	100	+	100	-5	1	7	0	0
chr1	100	+	100	0	1	13	0	0
chr1	100	+	100	10	0	0	1	20
chr1	100	+	100	-5	1	7	0	0
chr1	100	+	100	0	1	4	0	0
chr1	100	+	100	10	0	0	1	20
chr1	110	+	110	-10	0	0	2	17

EOF
	exit(1);
}

sub parse_command_line()
{
	my $region;

	my $rc = GetOptions(
			"help|h" => \&show_help,
			"region|r=s" => \$region,
			"mult|m" => \$use_multiplity_count,
			"debug" => \$debug,
			"up|upstream=i" => \$extend_upstream,
			"down|downstream=i" => \$extend_downstream,
			"intervals=s" => \$external_intervals,
			"h|header" => \$header,
		);

	die "Error: invalid command line options.\n" unless $rc;
	die "Error: missing region of interest (example: --region chrX:1000-2000).\n" unless $region;
	die "Error: external intervals file ($external_intervals) not found\n" if defined $external_intervals && ! ( -e $external_intervals ) ;

	set_filter_region($region);
}

sub set_filter_region($)
{
	my ($region) = shift or carp;

	($region_chrom, $region_start, $region_end) = $region =~ /^([\w-]+):(\d+)-(\d+)$/
		or die "Error: invalid region ($region), expecting 'chrX:1000-2000'.\n";

	$region_end--; #we assume the BED END position is 1-based

	die "Error: invalid region ($region), start is bigger than end position.\n" if ($region_start > $region_end);

	die "Error: invalid region ($region), size (end-start) is zero.\n" if ($region_start == $region_end);

	$region_size = $region_end - $region_start  + 1;

	print STDERR "Filter region: $region_chrom:$region_start-$region_end\n" if $debug;
}

=pod
  piRNA 5" starting position counter

  Copyright (C) 2011 A. Gordon (gordon at cshl dot edu)
                 and F. Muerdter (fmuerdte at cshl dot edu)

  License: Free for personal / academic / non-commercial use.

  See related publication:
     [ Production of artificial piRNAs in flies and mice ]
     [ Muerdter, F.*, Olovnikov, I.*, Molaro, A.*, Rozhkov, N.V.*, Czech, B., Gordon, A., Hannon, G.J., Aravin, A.A. ]
     [ RNA (2012) ]
     [ doi: 10.1261/rna.029769.111 ]

  For questions, please contact the corresponding authors:
    A. A. Aravin, E-mail aaa@caltech.edu.
	G. J. Hannon, E-mail hannon@cshl.edu.

  When using this script (or derivatives there-of),
  please cite: [ Muerdter et al., 2012 ]
=cut

sub get_next_interval()
{
	my $line = <>;
	return () unless defined $line; #End of file ?
	chomp $line;

	my ($chrom, $start, $end, $name, $quality, $strand, @other) = split /\t/, $line;

	# Do some input validation
	die "Input error: not enough columns in line $. (expecting at least 6 BED columns)\n" unless defined $strand;
	die "Input error: bad start value ($start) on line $. (expecting a number)\n" unless $start =~ /^\d+$/;
	die "Input error: bad end value ($end) on line $. (expecting a number)\n" unless $end =~ /^\d+$/;
	die "Input error: bad strand value ($strand) on line $. (expecting + or -)\n" unless $strand eq "+" || $strand eq "-";

	my $multiplicity = 1 ;

	#extract multiplicity value from read name
	if ($use_multiplity_count) {
		(undef, $multiplicity) = $name =~ /^(\d+)-([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)$/
			or die "Input error: invalid multiplicity value on column 4 ($name) on line 4. (expecting N-M)\n";
	}

	return $chrom, $start, $end-1, $name, $strand, $multiplicity;
}

sub interval_in_region($$$)
{
	my ($chrom, $start, $end) = @_ or carp;

	return undef unless $chrom eq $region_chrom;
	return undef unless $start >= $region_start;
	return undef unless $end <= $region_end;

	return 1;
}

sub initialize_empty_windows()
{
	$#window_positive_sequences = $region_size;
	$#window_negative_sequences = $region_size;
	$#window_positive_reads = $region_size;
	$#window_negative_reads = $region_size;

	foreach my $index ( 0 .. $region_size ) {
		$window_positive_sequences[$index] = 0 ;
		$window_negative_sequences[$index] = 0 ;
		$window_positive_reads[$index] = 0 ;
		$window_negative_reads[$index] = 0 ;
	}
}

sub load_external_intervals()
{
	open FILE,"<",$external_intervals
		or die "Error: failed to open external intervals file ($external_intervals): $!\n";
	while ( my $line = <FILE> ) {
		chomp $line;
		my ($chrom, $start, $end, $name, $quality, $strand, @other) = split /\t/, $line;

		# Do some input validation
		die "Input error: not enough columns in file $external_intervals line $. (expecting at least 6 BED columns)\n" unless defined $strand;
		die "Input error: bad start value ($start) in file $external_intervals line $. (expecting a number)\n" unless $start =~ /^\d+$/;
		die "Input error: bad end value ($end) in file $external_intervals line $. (expecting a number)\n" unless $end =~ /^\d+$/;
		die "Input error: bad strand value ($strand) in file $external_intervals line $. (expecting + or -)\n" unless $strand eq "+" || $strand eq "-";

		my $multiplicity = 1 ;

		push @intervals_to_process, [ $chrom, $start, $end, $name, $strand, $multiplicity ] ;
	}
	close FILE;
}

sub calculate_region_coverage()
{
	while (my ($chrom, $start, $end, $name, $strand, $multiplicity) = get_next_interval()) {
		$stats_lines_count++;
		next unless interval_in_region($chrom, $start, $end);
		$stats_sequences_count++;
		$stats_reads_count+=$multiplicity;

		if ( $strand eq "+" ) {
			$window_positive_sequences[$start - $region_start] += 1;
			$window_positive_reads[$start - $region_start] += $multiplicity;
		} else {
			$window_negative_sequences[$end - $region_start] += 1;
			$window_negative_reads[$end - $region_start] += $multiplicity;
		}

		# load the intervals from the same input file used for coverage calculations
		if (! defined $external_intervals ) {
			push @intervals_to_process, [ $chrom, $start, $end, $name, $strand, $multiplicity ] ;
		}
	}
}

sub sort_intervals_to_process()
{
	## this copmlicated sort function sorts by: Strand [4], then Start[1], then end [2]
	@intervals_to_process = sort {
			 ( $a->[4] ne $b->[4] ) ? ( $a->[4] cmp $b->[4] ) :
				 (
					 ( $a->[1] != $b->[1] ) ? ( $a->[1] <=> $b->[1] ) :
					 (
						$a->[2] <=> $b->[2]
					 )
				 );
			 } @intervals_to_process ;
}

sub scan_region_windows()
{
	print "intervals ".scalar @intervals_to_process ,"\n";
	## Iterave over every nucleotide in the region
	foreach my $interval_ref ( @intervals_to_process ) {

		my ($interval_chrom, $interval_start, $interval_end,
			$interval_name, $interval_strand, $interval_multiplicity ) = @{$interval_ref};

		print STDERR "processing interval $interval_chrom:$interval_start-$interval_end strand $interval_strand name '$interval_name' multiplicity $interval_multiplicity\n" if $debug;

		my $nuc_position = (( $interval_strand eq "+") ? $interval_start : $interval_end ) -  $region_start;
		my $genomic_nuc_position = $region_start + $nuc_position;


		##
		## Calculate the window around this read (up/down stream is relative to the interval's strand)
		##
		my ($extended_start, $extended_end);
		if ( $interval_strand eq "+" ) {
			$extended_start = $nuc_position - $extend_upstream;
			$extended_end   = $nuc_position + $extend_downstream;
		} else {
			$extended_start = $nuc_position - $extend_downstream;
			$extended_end   = $nuc_position + $extend_upstream;
		}
		$extended_start = 0 if $extended_start < 0 ;
		$extended_end = $region_size-1 if $extended_end > $region_size-1;
		
		 #print "extended region (start= $extended_start end= $extended_end nuc_pos=$nuc_position ) \n";
		## Scan the window around (before/after) the nucleotide position
		foreach my $window_index ( $extended_start .. $extended_end ) {
			my ($sequences_count_same, $sequences_count_opposite,
				$reads_count_same, $reads_count_opposite, $offset);

			##
			## Find the counts of reads/intervals in this window position
			##
			if ( $interval_strand eq "+" ) {
				$sequences_count_same = $window_positive_sequences[$window_index];
				$reads_count_same = $window_positive_reads[$window_index];
				$sequences_count_opposite = $window_negative_sequences[$window_index];
				$reads_count_opposite = $window_negative_reads[$window_index];

				$offset = $window_index - $nuc_position ;
			} else {
				$sequences_count_same = $window_negative_sequences[$window_index];
				$reads_count_same = $window_negative_reads[$window_index];
				$sequences_count_opposite = $window_positive_sequences[$window_index];
				$reads_count_opposite = $window_positive_reads[$window_index];

				$offset = $nuc_position - $window_index;
			}
			#print "indice = $window_index  ---off =  $offset\n";
			# At offset 0, subtract the current interval from the counts
			# (only if NOT using external intervals )
			
					
			
			if ( $offset == 0 && (!defined $external_intervals) ) {
				$sequences_count_same -= 1 ;
				$reads_count_same -= $interval_multiplicity;
			}

			if ( $sequences_count_same != 0 || $reads_count_same != 0 ||
				$sequences_count_opposite != 0 || $reads_count_opposite != 0) {
				print $region_chrom, "\t",
					$genomic_nuc_position, "\t",
					"+", "\t",
					$nuc_position, "\t",
					$offset, "\t",  ## This is the window offset
					$sequences_count_same, "\t",
					$reads_count_same, "\t",
					$sequences_count_opposite, "\t",
					$reads_count_opposite,
					"\n";
			}
		}
	}
}
