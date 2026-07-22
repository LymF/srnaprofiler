#!/usr/bin/bash




grep -v "@" $1  | cut -f 1 | nsort -T /media/temp_data/ | uniq -c | awk '{  if($1==1) print $2} ' > $3.single.ids
grep -v "@" $1  | cut -f 1 | nsort -T /media/temp_data/ | uniq -c | awk '{  if($1>1) print $2} ' > $3.multiple.ids 


seqtk subseq $2 $3.single.ids > $3.single_ids.fasta
seqtk subseq $2 $3.multiple.ids > $3.multiple_ids.fasta