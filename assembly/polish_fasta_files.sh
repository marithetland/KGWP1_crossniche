#!/usr/bin/env bash

# Script to perform long (medaka) and short-read (polypolish polca) polishing of fasta files from flye or unicycler and annotate with bakta
# Hetland M, 2023-04-21

#Give as input the isolate id and the medaka model, e.g.:
#bash polish_fasta_files.sh T7-394 r941_min_sup_g507 16

#Activate conda environments
source /opt/anaconda/anaconda3/etc/profile.d/conda.sh ;

#Input parameters
isolate_id=$1 #Read in isolate ID from arguments
illumina1=$(echo ${isolate_id}"_1.fastq.gz")
illumina2=$(echo ${isolate_id}"_2.fastq.gz")
medaka_model=$2 #"r941_min_sup_g507"
threads=$3



echo "Running long and short read polishing of: "${isolate_id}" using medaka model: " ${medaka_model}

conda activate clinopore_v2 ;

##########################################
#Run Medaka onfasta files
##########################################
echo "Running long-read polishing with the medaka model: " ${medaka_model}

#Run medaka consensus
medaka_consensus -d ${isolate_id}.fasta -o . -i ${isolate_id}_filtered.fastq.gz -t ${threads} -m ${medaka_model} ;

#Rename contigs
mv consensus.fasta ${isolate_id}_medaka_inter1.fasta ;
python ~/Programs/clinopore-nf/bin/contig_renaming.py ${isolate_id}.fasta  ${isolate_id}_medaka_inter1.fasta ${isolate_id}_medaka_inter2.fasta ${isolate_id}_medaka_inter3.fasta ;

#Sort the contig order from largest to smalles (i.e. place chromosome as first contig)
seqkit sort --by-length --reverse ${isolate_id}_medaka_inter3.fasta > ${isolate_id}_medaka.fasta ;

medaka_polished_assembly=${isolate_id}_medaka.fasta ;


##########################################
#Run Polypolish on fasta files
##########################################
echo "Running short-read polishing with polypolish"

bwa index ${medaka_polished_assembly} ;
bwa mem -t ${threads} -a ${medaka_polished_assembly} ${illumina1} > ${isolate_id}_r1.sam ;
bwa mem -t ${threads} -a ${medaka_polished_assembly} ${illumina2} > ${isolate_id}_r2.sam ;

#Run polypolish
python ~/Programs/clinopore-nf/bin/polypolish_insert_filter.py --in1 ${isolate_id}_r1.sam --in2 ${isolate_id}_r2.sam --out1 ${isolate_id}_filtered_r1.sam --out2 ${isolate_id}_filtered_r2.sam ;

~/Programs/clinopore-nf/bin/polypolish ${medaka_polished_assembly} ${isolate_id}_filtered_r1.sam ${isolate_id}_filtered_r2.sam| sed 's/_polypolish//' > ${isolate_id}_medaka_polypolish1.fasta ;

python ~/Programs/clinopore-nf/bin/contig_renaming.py ${isolate_id}.fasta ${isolate_id}_medaka_polypolish1.fasta ${isolate_id}_inter.fasta ${isolate_id}_medaka_polypolish2.fasta ;

seqkit sort --by-length --reverse ${isolate_id}_medaka_polypolish2.fasta > ${isolate_id}_medaka_polypolish.fasta ;
rm *sam ;

polypolish_polished_assembly=${isolate_id}_medaka_polypolish.fasta ;

conda deactivate ;

##########################################
#Run POLCA on fasta files
##########################################
echo "Running short-read polishing with polca"

conda activate polca_v2 ;

polca.sh -a ${polypolish_polished_assembly} -r "${illumina1} ${illumina2}" -t ${threads} -m 4G ;

mv ${isolate_id}_medaka_polypolish.fasta.PolcaCorrected.fa ${isolate_id}_medaka_polypolish_polca_intermediate.fasta ;

seqkit sort --by-length --reverse ${isolate_id}_medaka_polypolish_polca_intermediate.fasta > ${isolate_id}_medaka_polypolish_polca.fasta

python ~/Programs/clinopore-nf/bin/contig_renaming.py ${isolate_id}.fasta ${isolate_id}_medaka_polypolish_polca.fasta ${isolate_id}_inter.fasta ${isolate_id}_medaka_polypolish_polca.fasta ;

conda deactivate ;

##########################################
#Run BAKTA annotation of polished fastas
##########################################
#echo "Running bakta annotation"

#conda activate bakta_env ;
#bakta --db /opt/anaconda/anaconda3/envs/bakta_env/lib/python3.10/site-packages/bakta/db/db --complete --threads ${threads} --prefix ${isolate_id} --output ${isolate_id}_bakta -v ${isolate_id}_medaka_polypolish_polca.fasta ; 
#conda deactivate ;
