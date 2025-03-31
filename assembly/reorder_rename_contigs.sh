#!/usr/bin/env bash
################################################
# KGWP1 cross-niche paper
# 2023-06-29, Hetland MAK
# Commands to:
# - rotate circular contigs to avoid genes being split at the start/end of the sequence
# - reorder contigs by size
# - rename contigs to chr and pl_[0-9]
################################################

#Rotate contigs
for f in $(ls *.fasta | sed 's/.fasta//g') ; do circlator fixstart ${f}.fasta ${f}_rotated.fasta ; done ;

#Make it so there is only one line per fasta record
for f in $(ls *.fasta | sed 's/.fasta//g') ; do bash ~/GitHub/KGWP1_crossniche_paper/assembly/one_line_per_fasta.sh ${f}.fasta ; done ;

#Reorder contigs by size
conda activate seqkit_env ; 
for f in $(ls *.fasta | sed 's/.fasta//g') ; do seqkit sort --by-length --reverse ${f}.fasta >> ${f}_sorted.fasta ; done ;
conda deactivate ;

#Remove contig name from headers
sed -i.bak "s/_contig_[0-9][0-9]/_/g" *_sorted.fasta
sed -i.bak "s/_contig_[0-9]/_/g" *_sorted.fasta

#Remove the unicycler depth decimal points
sed -i.bak "s/\.[0-9][0-9]x//g" *_sorted.fasta

#Add replicon numbers from 1->
conda activate seqkit_env ; 
for f in $(ls *_rotated_sorted.fasta | sed 's/.fasta//g') ; do seqkit replace --pattern "__" --replacement '__contig_{nr} ' ${f}.fasta >> ${f}_renamed.fasta ; done
conda deactivate ;

#Rename contigs so that the first (largest) is chr, and the following pl_[0-9]:
for fasta in $(ls *_rotated_sorted_renamed.fasta); do number_of_replicons=$(grep -c ">" $fasta) ; for i in $(seq 1 $number_of_replicons); do sed -i.bak "s/__contig_1 /__chr /g" $fasta ; plasmids=$((i+1)) ; sed -i.bak "s/__contig_${plasmids} /__pl_${i} /g " $fasta ; sed -i.bak 's/_flye//g' $fasta ; sed -i.bak 's/_unicycler//g' $fasta  ;  done ; done ;

#Lastly, edit the headers so they have the structure: isolate_contig length circular (i.e. remove depth and any second lengths)
for fasta in $(ls *_rotated_sorted_renamed.fasta); do sed -E -i.bak 's/( depth=[0-9]+)//' $fasta ; sed -E -i.bak 's/( length=[0-9]+)$//' $fasta ; done