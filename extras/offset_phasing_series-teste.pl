#!/usr/bin/perl5.34.0

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
-ss	: Minimum size of phasing series (number of siRNAs)
-h      : Help message

";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $l1;
our $gsize;
our $pre;
my $plot;
my $density;
my $image;
my $ss;

GetOptions ("sam=s" => \$sam,"size=i" => \$gsize,"p=s" => \$pre,"ss=i" => \$ss, "d=i" => \$density,
	"h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($sam))) {
        die "\nERROR: [[ Give an input SAM file   ]]",$usage;
}
if (not(defined($gsize)) or $gsize < 0) {
        die "\nERROR: [[ Give the size of the genome used in sam files! ]]\n",$usage;
}

if (not(defined($density)) or $density < 0) {
        die "\nERROR: [[ Give the minimum density of siRNAs! ]]\n",$usage;
}
if (not(defined($ss)) or $ss < 0) {
        die "\nERROR: [[ Give the minimum  size of phasing serie! ]]\n",$usage;
}
my %dic;

if(exists $dic{$l1}){$s1 =$dic{$l1};}
if(exists $dic{$l2}){$s2 =$dic{$l2};}

#
#Loading SAM file and generating FPKM
#

our %pos;
my $max=0;
sub loadLib{
    $file=shift;
    $array=shift;
    $prefix = shift;
    
    open(L ,"<$file");
    
    ${$array}{"header"} = "";
     
    while(<L>){
        chomp;
        if(/^\w+/){
            
            @f           = split(/\t/,$_);
            $strand      = $f[1];
            $chromossome = $f[2];
            $seq_size    = length($f[9]);
            $start       = $f[3];
            $end         = $start + $seq_size -1;
            $ps          = "sense";
            
            #negative strand
            if($strand eq 16){
               #print "$start \t $end\n";
                $ps="antisense";
            }
            $id="$prefix"."_$ps"."_$start"."_$end";
                  
            if($start>$max){ $max=$start; }
                    
            if(not exists(${$array}{$id})){
                ${$array}{$id}{"strand"} = $strand;
                ${$array}{$id}{"sense"} = $ps;
                ${$array}{$id}{"start"} = $start;
                ${$array}{$id}{"end"} = $end;
               # print $f[9]." $seq_size start ($start) end ($end) \n";
                ${$array}{$id}{"ref"} = $f[2];
                ${$array}{$id}{"seq"} = $f[9];
                ${$array}{$id}{"si"} = $f[5];
                 ${$array}{$id}{"qual"} = $f[8];
                ${$array}{$id}{"score"} = 1;
            }else{
                ${$array}{$id}{"score"} += 1;
            }
            
        }else{
             ${$array}{"header"}.=$_ ."\n";;
        }
    }
    close(L1);

}

 
our %lib;

loadLib($sam,\%lib,"lib1");


our %Greads;

foreach $key (keys %lib){
        $st1 = ${$lib1}{"strand"};
        
        $sc= $lib{$key}{"score"};
        $sense= $lib{$key}{"sense"};
        $start=$lib{$key}{"start"};
        
       # print "$key \t $sc \n";
        if($sc >=$density){
           #print "$key \t $sc \n";
            $pos{$sense}{$start} =$sc;
            $pos{$sense}{$start}{"s"} = $sc;
        }
        
}



my $cont=1;


open(T,">$pre.T.log");



my %e_p;
my %e_n;

my %in_serie_p;
my %in_serie_n;
foreach $start1 (sort {$a<=>$b} keys %{$pos{"sense"}}){       
       
           # print "\n######################## $start1 \n";
            $cont=1;
            $id=$start1;
            $next = $start1 + ($cont*21);
            if(not exists $in_serie_p{$start1}){
	          	$start1=$next;
	            while(exists $pos{"sense"}{$next}){
	            	$in_serie_p{$next}=1;
	            	$next = $start1 + ($cont*21);
	            	$cont++;
	            	$start1= $next;
	            }
	            if($cont>1){
		            $e_p{$id}=$cont;
		          #  print "$id\t$cont\n";
	            }
	         }
}

#print "Antisense\n";
foreach $start1 (sort {$a<=>$b} keys %{$pos{"antisense"}}){       
       
           # print "\n######################## $start1 \n";
            $cont=1;
            $id=$start1;
            $next = $start1 + ($cont*21);
          	if(not exists $in_serie_n{$start1}){
	          	$start1=$next;
	            while(exists $pos{"antisense"}{$next}){
	            	$in_serie_n{$next}=1;
	            	$next = $start1 + ($cont*21);
	            	$cont++;
	            	$start1= $next;
	            }
	            if($cont > 1){
		            $e_n{$id}=$cont;
		          #  print "$id\t$cont\n";
	            }
			}
}

###creating vector with starts in sense    
my @start_s;        
foreach $id (keys %e_p){
	
	my $i=1;
	if($e_p{$id} >= $ss){
		push(@start_s,$id);
		while ($i < $e_p{$id}){
			$next = $id + ($i*21);
			push(@start_s,$next);
			$i++;
		}
	}

}
###creating vector with starts in antisense
my @start_n;        
foreach $id (keys %e_n){
	
	my $i=1;
	if($e_n{$id} >= $ss){
		push(@start_n,$id);
		while ($i < $e_n{$id}){
			$next = $id + ($i*21);
			push(@start_n,$next);
			$i++;
		}
	}

}


my %offset;

#calculando offset;
for (my $i=-30;$i<31;$i++){
	$offset{$i}=0;
}


foreach $startP (@start_s){
	
	foreach $startN (@start_n){
		$startN2 = $startN -20;
		$overlap = $startP - $startN2;
		 if(exists $offset{$overlap}){
		 	$offset{$overlap}= $offset{$overlap}+1;
		 }
	}
	
}
open(O,">$pre.offset");

print O "Overlap\tCount\n";
for (my $i=-30;$i<31;$i++){
	print O "$i\t". $offset{$i}."\n";
}
close(O);


`R --no-save $pre.offset $pre.offset < ~/bin/R/plotOffsetPhasingSeries-teste.R`;
