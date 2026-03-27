#!/usr/bin/perl5.34.0

use strict; # Added for better error checking
use warnings; # Added for better error checking
use Getopt::Long;

my $usage = "

#####################################################
# Script to calculate offset between phasing series #
#####################################################

$0 -sam <sam file>  -d <min density> -size <size genome> -p <prefix> [--plot]
$0 -h

-sam    : Library to be analyzed
-d      : Minimum density of siRNA
-size   : Size genome
-p      : Prefix for output files
-ss     : Minimum size of phasing series (number of siRNAs)
-h      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $sam;
our $gsize;
our $pre;
my $plot; # This variable is not used in the script, but kept for consistency with usage
my $density;
my $image; # This variable is not used in the script, but kept for consistency with usage
my $ss;
my $help;

GetOptions (
    "sam=s"    => \$sam,
    "size=i"   => \$gsize,
    "p=s"      => \$pre,
    "ss=i"     => \$ss,
    "d=i"      => \$density,
    "h!"       => \$help
);

if ($help) {
    die $usage;
}

if (not(defined($sam))) {
    die "\nERROR: [[ Give an input SAM file   ]]\n", $usage;
}
if (not(defined($gsize)) or $gsize < 0) {
    die "\nERROR: [[ Give the size of the genome used in sam files! ]]\n", $usage;
}

if (not(defined($density)) or $density < 0) {
    die "\nERROR: [[ Give the minimum density of siRNAs! ]]\n", $usage;
}
if (not(defined($ss)) or $ss < 0) {
    die "\nERROR: [[ Give the minimum  size of phasing serie! ]]\n", $usage;
}

#
#Loading SAM file and generating FPKM
#

our %pos;
my $max=0; # This variable is calculated but not used later. Kept for consistency.

sub loadLib{
    my ($file, $array_ref, $prefix) = @_; # Use my for lexical variables and destructure array

    open(my $L, "<", $file) or die "Cannot open $file: $!"; # Use lexical filehandle and proper error handling

    ${$array_ref}{"header"} = ""; # Initialize header

    while(<$L>){
        chomp;
        if(/^\@/){ # SAM header lines start with '@'
             ${$array_ref}{"header"} .= $_ ."\n";
             next; # Skip to next line if it's a header
        }

        # Process alignment lines
        my @f           = split(/\t/,$_);
        # Ensure there are enough fields before accessing them
        if (scalar @f < 11) { # A minimal SAM line has 11 fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL)
            warn "Skipping malformed SAM line: $_ \n";
            next;
        }

        my $strand_flag = $f[1]; # SAM flag for strand
        my $chromossome = $f[2];
        my $seq_size    = length($f[9]);
        my $start       = $f[3];
        my $end         = $start + $seq_size - 1;
        my $ps          = "sense"; # Default to sense

        # Check SAM flag for strand (16 indicates reverse strand)
        if (($strand_flag & 16) == 16){ # Bitwise AND to check the 16th bit
            $ps="antisense";
        }

        my $id="$prefix"."_$ps"."_$start"."_$end";

        if($start > $max){ $max=$start; } # Update max position

        # Use the dereferenced hash for updates
        if(not exists(${$array_ref}{$id})){
            ${$array_ref}{$id}{"strand"}    = $strand_flag;
            ${$array_ref}{$id}{"sense"}     = $ps;
            ${$array_ref}{$id}{"start"}     = $start;
            ${$array_ref}{$id}{"end"}       = $end;
            ${$array_ref}{$id}{"ref"}       = $f[2];
            ${$array_ref}{$id}{"seq"}       = $f[9];
            ${$array_ref}{$id}{"cigar"}     = $f[5]; # f[5] is CIGAR string
            ${$array_ref}{$id}{"qual"}      = $f[10]; # f[10] is QUAL
            ${$array_ref}{$id}{"score"}     = 1;
        }else{
            ${$array_ref}{$id}{"score"} += 1;
        }
    }
    close($L); # Corrected filehandle to $L
}

our %lib; # Declare %lib as our (global) for consistency with original script

loadLib($sam, \%lib, "lib1");


# Populate %pos based on density filter
# Corrected logic for %pos population
foreach my $key (keys %lib){
    # Skip the "header" key, as it's not a read entry
    next if $key eq "header"; # <--- THIS IS THE NEWLY ADDED FIX

    my $sc    = $lib{$key}{"score"};
    my $sense = $lib{$key}{"sense"};
    my $start = $lib{$key}{"start"};

    if ($sc >= $density){
        # Store the score directly. $pos{$sense}{$start} will hold the count/score.
        # This is what's expected by the phasing series calculation.
        $pos{$sense}{$start} = $sc;
    }
}

my $cont=1; # This variable is not used after this point, but kept for consistency.

open(my $T_log, ">", "$pre.T.log") or die "Cannot open $pre.T.log: $!"; # Use lexical filehandle
print $T_log "Processing log for $pre\n"; # Added a simple log message
close($T_log);


my %e_p; # Stores the count of siRNAs in a phasing series for 'sense'
my %in_serie_p; # To mark positions that have already been part of a series for 'sense'

# Calculate phasing series for sense strand
foreach my $initial_start (sort {$a<=>$b} keys %{$pos{"sense"}}){
    if (exists $in_serie_p{$initial_start}) {
        next; # Already processed as part of an earlier series
    }

    my $current_start = $initial_start;
    my $count = 0;
    while (exists $pos{"sense"}{$current_start}) {
        $in_serie_p{$current_start} = 1; # Mark as processed
        $count++;
        $current_start += 21; # Check the next expected phased position (21nt interval)
    }

    if ($count > 1) { # A series must have at least 2 siRNAs (original + 1 offset)
        $e_p{$initial_start} = $count;
        # print "Sense series found: Start $initial_start, Count $count\n"; # Debugging line
    }
}

my %e_n; # Stores the count of siRNAs in a phasing series for 'antisense'
my %in_serie_n; # To mark positions that have already been part of a series for 'antisense'

# Calculate phasing series for antisense strand
foreach my $initial_start (sort {$a<=>$b} keys %{$pos{"antisense"}}){
    if (exists $in_serie_n{$initial_start}) {
        next; # Already processed as part of an earlier series
    }

    my $current_start = $initial_start;
    my $count = 0;
    while (exists $pos{"antisense"}{$current_start}) {
        $in_serie_n{$current_start} = 1; # Mark as processed
        $count++;
        $current_start += 21; # Check the next expected phased position (21nt interval)
    }

    if ($count > 1) { # A series must have at least 2 siRNAs
        $e_n{$initial_start} = $count;
        # print "Antisense series found: Start $initial_start, Count $count\n"; # Debugging line
    }
}

###creating vector with starts in sense
my @start_s;
foreach my $id (keys %e_p){
    if($e_p{$id} >= $ss){
        # Add the initial start of the series
        push(@start_s, $id);
        # Add subsequent positions in the series
        for (my $i = 1; $i < $e_p{$id}; $i++){
            my $next = $id + ($i * 21);
            push(@start_s, $next);
        }
    }
}
# print "Sense starts in series: " . join(", ", sort {$a<=>$b} @start_s) . "\n"; # Debugging line

###creating vector with starts in antisense
my @start_n;
foreach my $id (keys %e_n){
    if($e_n{$id} >= $ss){
        # Add the initial start of the series
        push(@start_n, $id);
        # Add subsequent positions in the series
        for (my $i = 1; $i < $e_n{$id}; $i++){
            my $next = $id + ($i * 21);
            push(@start_n, $next);
        }
    }
}
# print "Antisense starts in series: " . join(", ", sort {$a<=>$b} @start_n) . "\n"; # Debugging line


my %offset;

# Initialize offset hash with zeros for range -30 to 30
for (my $i=-30;$i<31;$i++){
    $offset{$i}=0;
}

# Calculating offset
# Re-evaluate the biological meaning of $startN2 and $overlap if the results are not as expected.
# The current calculation is: $startP - ($startN - 20)
# This means: (Sense 5' end) - (Antisense 5' end - 20)
# If you expect a 10bp overlap between 5' ends, you might want $startP - $startN.
foreach my $startP (@start_s){
    foreach my $startN (@start_n){
        my $startN2 = $startN - 20; # This specific shift might be the source of unexpected offsets
        my $overlap = $startP - $startN2;

        if(exists $offset{$overlap}){
            $offset{$overlap}++; # Increment count for this specific overlap
        }
    }
}

open(my $O_file, ">", "$pre.offset") or die "Cannot open $pre.offset: $!"; # Use lexical filehandle

print $O_file "Overlap\tCount\n";
for (my $i=-30;$i<31;$i++){
    print $O_file "$i\t". $offset{$i}."\n";
}
close($O_file);


`Rscript ~/bin/R/plotOffsetPhasingSeries-teste.R "$pre.offset" "$pre.offset"`;
# Or better, just pass the filename once if the R script only needs one.
# If args[4] is meant to be the prefix, then:
# `Rscript ~/bin/R/plotOffsetPhasingSeries.R "$pre.offset" "$pre"`;