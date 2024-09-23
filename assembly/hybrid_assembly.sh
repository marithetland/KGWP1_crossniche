#!/usr/bin/env bash
################################################
# KGWP1 cross-niche paper
# Commands to:
# - run hybrid assembly of 536 genomes
################################################

################################################
#Clinopore assembly
################################################
cd ./ONT/FASTQ/ ; 
#Set up clinopore folder
mkdir clinopore-nf && cd $_
cp ~/Programs/clinopore-nf/* .
cp -r ~/Programs/clinopore-nf/bin .
cd ..

#Symlink Illumina FASTQs for use in pipeline: 
ls *.fastq.gz | sed 's/.fastq.gz//g' >> list_of_genomes.txt
bash ~/Scripts/symlink_trimmedPEfastqfiles.sh list_of_genomes.txt ./Illumina/trimmed_reads/ .
rename 's/_1_val_1.fq.gz/_1.fastq.gz/g' *gz #Rename so it will work with clinopore
rename 's/_2_val_2.fq.gz/_2.fastq.gz/g' *gz #Rename so it will work with clinopore

#Change these four parameters/paths in the nextflow.config file to look like:
#params.max_retries = 60
#params.queue_size = 60
#params.processors = 60
#params.reads='./ONT/FASTQ/*fastq.gz'

#Run clinopore
conda activate clinopore_v2 ;
nextflow run clinopore.nf --clinopore_env /opt/anaconda/anaconda3/envs/clinopore_v2 --polca_env /opt/anaconda/anaconda3/envs/polca_v2 | tee clinopore_log.txt ; 
conda deactivate ;

################################################
#Unicycler assembly
################################################
#Run filtlong first (assumes that the trimmed short-reads are also present in the folders)
conda activate clinopore_v2 ;
parallel --jobs 10 ' echo {} ; filtlong --min_length 1000 --keep_percent 95 {}.fastq.gz | gzip > {}_filtered.fastq.gz ' ::: $(ls *_1.fastq.gz | sed 's/_1.fastq.gz//g') ;
conda deactivate ;
#Then unicycler
conda activate unicycler_for_klebgapont ;
for f in $(ls *_1.fastq.gz | sed 's/_1.fastq.gz//g') ; do unicycler -1 ${f}_1.fastq.gz -2 ${f}_2.fastq.gz -l ${f}_filtered.fastq.gz -o ${f}_unicycler_hybrid_assembly -t 40 --keep 2 --verbosity 2  ; done ;
conda deactivate ;

#Rename the final unicycler fasta and graph files
for f in $(ls -d *assembly) ; do name=$(ls -d $f) ; echo $name ; mv ${f}/assembly.fasta ${f}/${name}.fasta ; done
for f in $(ls -d *assembly) ; do name=$(ls -d $f) ; echo $name ; mv ${f}/assembly.gfa ${f}/${name}.gfa ; done
