#!/usr/bin/perl5.32.1 

use Getopt::Long;

my $usage = "

$0 -s <sam file> -p <prefix> -f <output folder> -c <chromosome> -r <region> 
$0 -h

-l1     : Library 1
-s1     : Size library l1
-l1     : Library 2
-s1     : Size library l2
-size   : Size genome
-h      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $sam;
my $prefix;
my $folder;
my $chr;
my $region;
my $si;
my $se;

GetOptions ("s=s" => \$sam,"p=s" => \$prefix,"f=s" => \$folder,"r=s" => \$region,"c=s" => \$chr,"si=i" => \$si,"se=i" => \$se, "h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($sam))) {
        die "\nGive an sam input file",$usage;
}
if (not(defined($prefix))) {
        die "\nGive an prefix file",$usage;
}
if (not(defined($folder))) {
        die "\nGive an folder name",$usage;
}
if (not(defined($chr))) {
        die "\nGive an chromosome name",$usage;
}

if (not(defined($region))) {
        die "\nGive an region in format 1-1000",$usage;
}

if(not defined($se) or $se < 0 or $se < $si){
    die("\n\n\t\t[Bad value for end]\n");

}

if(not defined($si) or $si < 0){
    die("\n\n\t\t[Bad value for start]\n");
}


my $aux1 ="perl5.34.0 /home/ericgdp/bin/splitSamByStrand.pl -in $sam -out $folder/$prefix";

my $c2_1 = "perl5.34.0 /home/ericgdp/bin/filterSamBowtieBySize.pl -i $folder/$prefix.plus.sam -si $si -se $se -o $folder/$prefix.plus.$si-$se.sam";
my $c2_2 = "perl5.34.0 /home/ericgdp/bin/filterSamBowtieBySize.pl -i $folder/$prefix.minus.sam -si $si -se $se -o $folder/$prefix.minus.$si-$se.sam";
my $c2_3 = "perl5.34.0 /home/ericgdp/bin/filterSamBowtieBySize.pl -i $sam -si $si -se $se -o $folder/$prefix.$si-$se.sam";

my $c4_2 = "bash /home/ericgdp/bin/samToBed.sh $folder/$prefix.minus.$si-$se.sam $folder/$prefix.minus.$si-$se";
my $c4_1 = "bash /home/ericgdp/bin/samToBed.sh $folder/$prefix.plus.$si-$se.sam $folder/$prefix.plus.$si-$se";
my $c4_3 = "bash /home/ericgdp/bin/samToBed.sh $folder/$prefix.$si-$se.sam $folder/$prefix.$si-$se";


my $pha_pos = "perl5.34.0 /home/ericgdp/bin/bed_calculate_offsets.pl --region $chr:$region --up 80 --down 80 --header $folder/$prefix.plus.$si-$se.bed | grep -v \"intervals\" | grep -v \"chrom\"  | sed 1d | sort -k5,5 | groupBy -g 5 -c 7 -o sum > $folder/$prefix.plus.$si-$se.phasing "; 

my $pha_neg = "perl5.34.0 /home/ericgdp/bin/bed_calculate_offsets.pl --region $chr:$region --up 80 --down 80 --header $folder/$prefix.minus.$si-$se.bed | grep -v \"intervals\" | grep -v \"chrom\"  | sed 1d | sort -k5,5 | groupBy -g 5 -c 7 -o sum > $folder/$prefix.minus.$si-$se.phasing "; 

my $off = "perl5.34.0 /home/ericgdp/bin/bed_calculate_offsets.pl --region $chr:$region --up 80 --down 80 --header $folder/$prefix.$si-$se.bed | grep -v \"intervals\" | grep -v \"chrom\"  | sed 1d | sort -k5,5 | groupBy -g 5 -c 9 -o sum > $folder/$prefix.$si-$se.offset ";


print "split sam...\n";
`$aux1`;

print "filter sam by interval size...\n";
`$c2_1`;
`$c2_2`;
`$c2_3`;

`$c4_1`;
`$c4_2`;
`$c4_3`;

`$pha_pos`;
`$pha_neg`;
`$off`;

#`usr/bin/R --no-save $folder/$prefix.plus.$si-$se.phasing $folder/$prefix.minus.$si-$se.phasing $folder/$prefix.$si-$se.offset $folder/$prefix < /home/ericgdp/bin/R/plotPhasingOffset_2018.R 2>/dev/null`;
print "R --no-save $folder/$prefix.plus.$si-$se.phasing $folder/$prefix.minus.$si-$se.phasing $folder/$prefix.$si-$se.offset $folder/$prefix < /home/ericgdp/bin/phasing/R/plotPhasingOffset_2018.R 2>/dev/null\n";
`usr/bin/R --no-save $folder/$prefix.plus.$si-$se.phasing $folder/$prefix.minus.$si-$se.phasing $folder/$prefix.$si-$se.offset $folder/$prefix < /home/ericgdp/bin/R/phasing-plot-teste.R 2>/dev/null`;


print  "cleaning up...\n";
#`rm -rf $folder/$prefix.$si-$se.sam $folder/$prefix.minus.$si-$se.sam $folder/$prefix.plus.$si-$se.sam`;
#`rm -rf $folder/$prefix.plus.sam $folder/$prefix.minus.sam`;
#`rm -rf $folder/$prefix.plus.$si-$se.phasing $folder/$prefix.minus.$si-$se.phasing $folder/$prefix.$si-$se.offset`;
#`rm -rf $folder/$prefix.plus.$si-$se.bed $folder/$prefix.minus.$si-$se.bed $folder/$prefix.$si-$se.bed `;

