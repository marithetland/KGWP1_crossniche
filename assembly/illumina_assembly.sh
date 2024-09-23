#!/usr/bin/env bash
#################################################
# KGWP1 cross-niche paper
# Commands to:
# - run assembly of 3255 Illumina genomes
#################################################

#################################################
#Environment and versions
################################################# 
cd ./Illumina/FASTQ/ ;
conda activate klebgap_assembly_2022 ; 
#Versions 
fastqc -version #0.11.9
multiqc --version #1.14
trim_galore -v #v0.6.7
cutadapt --version #v3.7 
spades.py -v #v.3.15.4
quast.py -v #v5.2.0
python -—version #v3.9.0

#################################################
#Read trimming
#################################################
#Trim reads with trim-galore, paired-end mode (otherwise default)
parallel --jobs 20 'echo {} ; trim_galore --paired --cores 4 {}_1.fastq.gz {}_2.fastq.gz ' ::: $(ls *_1.fastq.gz | sed 's/_1.fastq.gz//g') ;
mkdir ../trimmed_reads ; mv *val*gz ../trimmed_reads ;

#################################################
#SPAdes assembly
################################################# 
cd ./Illumina/trimmed_reads/ ;
#Assemble trimmed reads with spades using isolat option (otherwise default)
parallel --jobs 5 'echo {} ; spades.py -o {}_spades_assembly --isolate -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz --threads 16 ' ::: $(ls *_1_val_1.fq.gz | sed 's/_1_val_1.fq.gz//g') ; 

#Set output assembly files prefix to parent dir (isolate name)  
parallel --jobs 80 'echo {} ; cd {} ; bash ~/Scripts/prefix_files_with_parentName_recursively.sh ; cd .. ' ::: $(ls -d *y ) ;

#################################################
#QC
#################################################
#FastQC
cd ./Illumina/ ;
parallel --jobs 80 'echo {} ; fastqc {} -o ./QC/fastQC_trimmed/ >> ./logs/{}_fastqc_trimmed.log 2>&1' ::: $(ls ./trimmed_reads/*gz )
multiqc ./QC/fastQC_trimmed/ -f -o ./QC/multiqc_trimmed/

#Quast 
quast.py *.fasta

#Read Count
cd ~/Programs/ ; git clone https://github.com/rrwick/MinION-desktop ; make ; #Install/Setup
cd ./Illumina/FASTQ ;
~/Programs/MinION-desktop/fast_count >> KLEBGAP_fastcount.txt #Print header
~/Programs/MinION-desktop/fast_count *gz >> KLEBGAP_fastcount.txt #Get fast count stats for each file

#To estimate average read depth per genome, take total number of bases from fast_count and divide by the total length of the assembly from quast.