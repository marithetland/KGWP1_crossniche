for fasta in $(ls *.fasta); do number_of_replicons=$(grep -c ">" $fasta) ; for i in $(seq 1 $number_of_replicons); do 
sed -i.bak "s/__contig_1 /__chr /g" $fasta ; plasmids=$((i+1)) ; sed -i.bak "s/__contig_${plasmids} /__pl_${i} /g " 
$fasta ; sed -i.bak 's/_flye//g' $fasta ; sed -i.bak 's/_unicycler//g' $fasta  ;  done ; done 
