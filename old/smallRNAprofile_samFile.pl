#!/usr/bin/perl



use Bio::SeqIO;  
use Bio::Seq::Quality;
use Getopt::Long;

my $usage = "

$0 -sam <sam file> -fa <fasta file> -p <prefix name> -r <reads file> -o <output image folder> -n <normalizing value> -si <minimum #reads in cluster siRNAs> -pi <minimum #reads in cluster piRNAs> [--profile] [--fastq]
$0 -h

-i <input file>         : Input file with columns delimited by a specific delimiter
-h                      : Help message

#
# calc stattistics by reads
#
";

$| = 1;     # forces immediate prints into files rather than the buffer.


my $sam;
my $delimiter;
my $index;
my $profile;
my $fasta;
my $counts;
my $prefix;
my $mp;
my $mS;
my $out_folder;
my $reads;
my $norm;
my $fastq;

GetOptions ("sam=s" => \$sam,"fa=s" => \$fasta,"o=s" => \$out_folder,"si=i" => \$ms,"pi=i" => \$mp,"p=s" => \$prefix,"r=s" => \$reads,"n=s" => \$norm,
"fastq!" => \$fastq, "profile!" => \$profile,"h!" => \$help);

if ($help) {
        die $usage;
}

if (not(defined($sam))) {
        die "\nGive an input file name",$usage;
}
if (not(defined($prefix))) {
        die "\nGive an prefix file name",$usage;
}
if ( defined($profile) and not(defined($fasta))) {
        die "\nGive an fasta file ",$usage;
}

if ( not defined($mp) or $mp <0 ) {
        die "\nGive an minimum of reads for piRNA clusters ",$usage;
}
if ( not defined($ms) or $ms <0 ) {
        die "\nGive an minimum number of reads for siRNA clusters ",$usage;
}

if ( not defined($norm) or $norm <0 ) {
        die "\nGive an minimum normalizing value ",$usage;
}

if ( not(defined($reads))) {
        die "\nGive an reads input file ",$usage;
}


open(IN,"<$sam");

our %tamanhos;
our %tamanhos_p;
our %tamanhos_n;
our $cp=0;
our $cn=0;
our %ch;
our %ch2;

my $min=1000;
my $max=0;


if(defined($fasta)){
	our %contig;

	warn("Loading FASTA file...\n");
	my $in = Bio::SeqIO->new(-format    => 'fasta',
							   -file      => $fasta);
	while (my $data = $in->next_seq) {
		$size = $data->length;
		$contig{$data->primary_id}=$data->seq;
	}



}

my %count_piRNAs;

my %reads;
warn("Loading SAM file...\n");
while (<IN>){
	
  if(/^\w+/){
	@campos = split(/\t/,$_);
	$size = length($campos[9]);
	$c = $campos[2];	
	if(not exists($reads{$campos[2]}{$campos[0]})){
	
	$reads{$campos[2]}{$campos[0]}=1;
	 
	 if (not exists $tamanhos_p{$size} ){
                $tamanhos_p{$size} = 0;
        }
	   if (not exists $ch{$c}{"pos"} ){
                $ch{$c}{"pos"} = 0;
        }
        
        if (not exists $ch{$c}{"neg"} ){
                $ch{$c}{"neg"} = 0;
        }

	 if (not exists $tamanhos_n{$size} ){
                $tamanhos_n{$size} = 0;
        }
	
	if(not exists $ch2{$c}{"pos"}{$size}){
	    $ch2{$c}{"pos"}{$size}=0;
	}
	
	if(not exists $ch2{$c}{"neg"}{$size}){
	    $ch2{$c}{"neg"}{$size}=0;
	}
	
	if($size < $min) { $min=$size;}
	if($size > $max) { $max=$size;}
	
	if ($campos[1] eq "0"){   
	
		$cp++;
		$ch{$c}{"pos"} = $ch{$c}{"pos"} + 1;
		$tamanhos_p{$size}= $tamanhos_p{$size} + 1;
	
		$ch2{$c}{"pos"}{"$size"}= $ch2{$c}{"pos"}{"$size"} +1;
		#print "POSITIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"pos"}{$size} ."\n";
		
	}elsif($campos[1] eq "16"){
	
		$cn++;
		$ch{$c}{"neg"} = $ch{$c}{"neg"} + 1; 
		$tamanhos_n{$size}= $tamanhos_n{$size} + 1;
		$ch2{$c}{"neg"}{$size}= $ch2{$c}{"neg"}{$size} + 1;
		#print "NEGATIVE in [$c] fita [pos]  size [$size] amount: ".$ch{$c}{"neg"}{$size} ."\n";
		
	}

  	
  }else{
  	#read exists
  }
  
  }
}
close(IN);
$t = $cp+$cn;


#
#filtering by size
#
warn("Filtering by minimum size ($minSize)...\n");
my %ch_after_filter_size;
foreach $c (keys %ch){
        $pa = $ch{$c}{"pos"};
        $na=  $ch{$c}{"neg"};
        $ta=$pa+$na;
        if($ta >= $ms or $ta >=$mp){
        	$ch_after_filter_size{$c}=1;
        }
}



open(S,">$prefix.geral_count.stats");
print S "Size\tSENSE\tANTISENSE\tTOTAL\n";
foreach (sort { $a <=> $b } keys(%tamanhos_p) )
{
	$s = $tamanhos_p{$_}  + $tamanhos_n{$_};
        if ($s ne "0"){
                print S "$_\t".$tamanhos_p{$_}."\t".$tamanhos_n{$_}."\t".$s. "\n";
        }

}
close(S);

open(S,">$prefix.count_by_chromosome.stats");
print S "Chr\tSize\tSENSE\tANTISENSE\tTOTAL\n";
foreach $c (keys %ch){
        $pa = $ch{$c}{"pos"};
        $na=  $ch{$c}{"neg"};
        $ta=$pa+$na;
        
        if($pa <=0){ $pa=0;} 
        if($na <=0){ $na=0;} 
            
        print S $c."\t". $ta . "\t" . $na .  "\t". $ta ."\n";
    }
close(S);
    


##
warn("Classifying sequences according to small RNA profile...\n");
##

my %signal;

foreach $c (keys %ch_after_filter_size){
     $signal{$c}{"si" } = "no" ;
     $signal{$c}{"pi" } = "no" ;
     $signal{$c}{"mi" } = "no" ;
     $signal{$c}{"final" } = "undefined" ;
     
	$sumT1P=0;
	$sumT1N=0;
	$allT1=0;
	$sumT2P=0;
	$sumT2N=0;
	$allT2=0;
	$allT4=0;
	$allT1=0;
	$allT4_pos=0;
	$allT4_neg=0;
	$scorePI=0;
	$allT3=0;
	
	for (my $s=$min;$s<=$max;$s++){
        $pa = $ch2{$c}{"pos"}{$s};
        $na=  $ch2{$c}{"neg"}{$s};
        
        if($pa <=0){ $pa=0;} 
        if($na <=0){ $na=0;} 
       
        $ta=$pa+$na;
        
        #degradation 
        if($s >= 15 and $s <= 19){
    
           	$allT1= $allT1 +$ta;
        }
        
        #siRNAs
        if($s>=20 and $s <=22){
        	
        	$allT4= $allT4 + $ta;
        	$allT4_pos=$allT4_pos + $pa;
        	$allT4_neg=$allT4_neg + $na;
        	
        	if($s==22){
        		
        		if($allT4_pos ==0){$allT4_pos=1;}
        		if($allT4_neg ==0){$allT4_neg=1;}
        		
        		if(($allT4_pos > 50 or $allT4_neg > 50) and ($allT4 > $allT1) ){
        			$r=$allT4_pos / $allT4_neg;
        			$r=$allT4_neg/$allT4_pos if $allT4_neg>$allT4_pos;
        		
        			$signal{$c}{"si" } = "yes" ;
        		 
        		}
        	}
                 
        }
        #miRNAs
        if($s>=22 and $s<=24){
        	$sumT2P+=$pa;
        	$sumT2N+=$na;
        	$ta=$pa +$na;
        	$allT2= $allT2 +$ta;
        	
        	if($s==24){
        		if(($allT2 > 50) and ($allT2 > $allT1) and ($allT2 > $allT3) ){
        			if($sumT2P ==0){$sumT2P=1;}
        			if($sumT2N ==0){$sumT2N=1;}
        			
        			$r=$sumT2P / $sumT2N;
        			$r=$sumT2N/$sumT2P if $sumT2N>$sumT2P;
        			
        			if($r > 1000){
        				$signal{$c}{"mi" } = "yes" ;
        			} 
        		}
        	}
        
        }
        #piRNAs
        if($s>=24 and $s<=29){
        	
        	 	if(($pa > 10 or $na > 10) ){
					$scorePI +=1;
				}
		
			    $allT3 = $allT3 + $pa + $na;  
		      #	print "$s->$allT3\n";
        }
        if($s==30){
        	$a=$allT2 +$allT1;
        	#print "$s -> piRNAs($c)\t$allT3\t$scorePI \t sumt1t2($a)\t classification:";
        	if( ($allT3 > $allT1 ) and ($scorePI >= 4) ) {
        		 $signal{$c}{"pi" } ="yes" ;
        	}
        
        }
        
        
        #print $s."\t$pa\t$na\t$ta\n";
    }
    #$allT4: siRNAs
    #$allT3: piRNAs
    #$allT2: miRNAs
    #$allT1: 15-19nt degradation
    
    if( ( $allT4 > $allT3) && ($allT4>$allT1) && $signal{$c}{"si" } eq "yes"  &&  ($r > 0.3 and $r < 10) && ($allT4>=$ms) ) {
    	$signal{$c}{"final" } ="siRNA";
    }
    if( ($allT3/3 > $allT4) &&  ($allT3>$allT1) && ($signal{$c}{"pi" } eq "yes")  && ($allT3>=$mp)){
    	$signal{$c}{"final" } ="piRNA";
    }
     if( ($allT2 > $allT1) &&  ($allT2>$allT3) &&  ($allT2>$allT4) &&  $signal{$c}{"mi" } eq "yes" ){
    	$signal{$c}{"final" } ="miRNA";
    }
    
    if( ($allT1 > $allT3) &&  ($allT1>$allT4)  &&  ($allT1>$allT2) ){
    	$signal{$c}{"final" } ="degradation";
    }
    print   $signal{$c}{"final" } ."\n";
   # warn("chromosome($c), $allT1,$allT4,$allT3,$allT2,[".$signal{$c}{"final" }."]\n");
    
}

open(S,">$prefix.classification.stats");
print S "Contig\tsiRNA\tpiRNA\tmiRNA\tFINAL\n" ;
foreach $c (keys %signal){

	print S "$c\t".$signal{$c}{"si" } ."\t" . $signal{$c}{"pi" }."\t"  . $signal{$c}{"mi" }."\t" . $signal{$c}{"final" }." \n" ;
}
close(S);

if (defined($profile)){
	open(S,">$prefix.profile");
	open(SI,">$prefix.withSiRNA.fasta");
	open(PI,">$prefix.withPiRNA.fasta");
	open(SIPI,">$prefix.withSiRNA_and_PiRNA.fasta");
	
	foreach $c (keys %signal){
		print S " $c\t".$signal{$c}{"si" } ." \t" . $signal{$c}{"pi" }." \n" ;
		if($signal{$c}{"si"} eq "yes" and $signal{$c}{"pi"} eq "no" ){
			print SI ">$c\n".$contig{$c}."\n";
		}
		if($signal{$c}{"pi"} eq "yes" and $signal{$c}{"si"} eq "no"){
			print PI ">$c\n".$contig{$c}."\n";
		}
		if($signal{$c}{"pi"} eq "yes" and $signal{$c}{"si"} eq "yes"){
			print SIPI ">$c\n".$contig{$c}."\n";
		}
	}
	
	close(S);
	close(SI);
	close(PI);
	close(SIPI);

	foreach $c (keys %signal){
	
		$profile = $signal{$c}{"final"};
		open(S,">>$prefix.FINAL.$profile.fasta");
		print S ">$c\n".$contig{$c}."\n";
		close(S);

	}
}

if($min > 5) { $min=5;} 
if($max < 30) { $max=35;} 

open(S,">$prefix.distribution_stats");
foreach $c (keys %ch2){
	print S "\nChromossome [$c]\n";
	for (my $s=$min;$s<=$max;$s++){
        $pa = $ch2{$c}{"pos"}{$s};
        $na=  $ch2{$c}{"neg"}{$s};
        
        if($pa <=0){ $pa=0;} 
        if($na <=0){ $na=0;} 
       
        $ta=$pa+$na;
        print S $s."\t$pa\t$na\t$ta\n";
    }
}
close(S);


if(1>2){
warn("generating sam files...\n");
open(IN,"<$sam");
while (<IN>){
	
  if(/^\w+/){
	@campos = split(/\t/,$_);
	$size = length($campos[9]);
	$c = $campos[2];	
	$profile = $signal{$c}{"final"};
	open(S,">>$prefix.$profile.sam");
  	print S $_;
  	close(S);
		
	  
  }else{
  	foreach $chr (keys %signal){
  		$profile = $signal{$chr}{"final"};
  		open(S,">>$prefix.$profile.sam");
  		print S $_;
  		close(S);
  	}
  }
}
close(IN);

}

`mkdir $out_folder`;

warn("Formatting databases...\n");
`bowtie-build $prefix.FINAL.siRNA.fasta $prefix.FINAL.siRNA.fasta`;
`bowtie-build $prefix.FINAL.piRNA.fasta $prefix.FINAL.piRNA.fasta`;


if(defined($fastq)){
	$type="-q";
}else{
		$type="-f";
}

warn("Mapping reads in databases...\n");
`bowtie $type -S  -p 40 -k 10  -v 1  $prefix.FINAL.siRNA.fasta $reads   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.siRNAs.aux.sam`;
`bowtie $type -S  -p 40 -k 10  -v 1  $prefix.FINAL.piRNA.fasta $reads   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.piRNAs.aux.sam`;


warn("Strading clusters by high abundant strand...\n");
`samStatistics_v2.pl -sam $prefix.FINAL.siRNAs.aux.sam --counts > $out_folder/$prefix.stats`;
`samStatistics_v2.pl -sam $prefix.FINAL.piRNAs.aux.sam --counts > $out_folder/$prefix.stats`;

`defineStrandBasedOnCounts.pl -c $prefix.FINAL.siRNAs.sam.aux.counts -f $prefix.FINAL.siRNA.fasta -p $out_folder/$prefix.FINAL.siRNA.stranded.fasta`;
`defineStrandBasedOnCounts.pl -c $prefix.FINAL.piRNAs.sam.aux.counts -f $prefix.FINAL.piRNA.fasta -p $out_folder/$prefix.FINAL.piRNA.stranded.fasta`;

`rm -rf $prefix.FINAL.siRNAs.aux.sam $prefix.FINAL.piRNAs.aux.sam`;

`bowtie-build $out_folder/$prefix.FINAL.siRNA.stranded.fasta $out_folder/$prefix.FINAL.siRNA.stranded.fasta`;
`bowtie-build $out_folder/$prefix.FINAL.piRNA.stranded.fasta $out_folder/$prefix.FINAL.piRNA.stranded.fasta`;



warn("Separating single and multiple mappings...\n");
`sepSinMult.sh $sam $reads $prefix.host`;

`bowtie $type -S  -p 40 -k 10  -v 1  $out_folder/$prefix.FINAL.siRNA.stranded.fasta $reads   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.siRNAs.sam`;

`bowtie $type -S  -p 40 -k 10  -v 1  $out_folder/$prefix.FINAL.siRNA.stranded.fasta $prefix.host.single_ids.fasta   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.siRNAs.SINGLE.sam`;
`bowtie $type -S  -p 40 -k 10  -v 1  $out_folder/$prefix.FINAL.siRNA.stranded.fasta $prefix.host.multiple_ids.fasta   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.siRNAs.MULTIPLE.sam`;


`bowtie $type -S  -p 40 -k 10  -v 1  $out_folder/$prefix.FINAL.piRNA.stranded.fasta $reads   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.piRNAs.sam`;

`bowtie $type -S  -p 40 -k 10  -v 1  $out_folder/$prefix.FINAL.piRNA.stranded.fasta $prefix.host.single_ids.fasta   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.piRNAs.SINGLE.sam`;
`bowtie $type -S  -p 40 -k 10  -v 1  $out_folder/$prefix.FINAL.piRNA.stranded.fasta $prefix.host.multiple_ids.fasta   | awk -F '\t' '{if( \$2 != 4) print \$0}' > $prefix.FINAL.piRNAs.MULTIPLE.sam`;


warn("Calculatting pattern...\n");
`calcPatternInSamFile.pl -s $prefix.FINAL.siRNAs.sam -o $out_folder/$prefix.SIRNAS.pattern`;
`calcPatternInSamFile.pl -s $prefix.FINAL.piRNAs.sam -o $out_folder/$prefix.PIRNAS.pattern`;
`calcPatternInSamFile.pl -s $prefix.FINAL.siRNAs.SINGLE.sam -o $out_folder/$prefix.SIRNAS.unique.pattern`;
`calcPatternInSamFile.pl -s $prefix.FINAL.piRNAs.SINGLE.sam -o $out_folder/$prefix.PIRNAS.unique.pattern`;

warn("Mapping reads in databases...\n");
`plotGeralDistributionPerBaseByReads.pl -sam $prefix.FINAL.siRNAs.sam -s 15 -e 35 -p $out_folder/$prefix.SIRNAS --plot --norm $norm`;
`plotGeralDistributionPerBaseByReads.pl -sam $prefix.FINAL.siRNAs.SINGLE.sam -s 15 -e 35 -p $out_folder/$prefix.SIRNAS.unique --plot --norm $norm`;
`plotGeralDistributionPerBaseByReads.pl -sam $prefix.FINAL.piRNAs.sam -s 15 -e 35 -p $out_folder/$prefix.PIRNAS --plot --norm $norm`;
`plotGeralDistributionPerBaseByReads.pl -sam $prefix.FINAL.piRNAs.SINGLE.sam -s 15 -e 35 -p $out_folder/$prefix.PIRNAS.unique --plot --norm $norm`;

