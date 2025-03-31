#!/usr/bin/env bash
################################################
# KGWP1 cross-niche paper
# Pipeline to:
# - run alignment
# - run snp-dists to get pairwise distances
# For each of 107 SLs that are niche-overlapping
################################################

################################################
#Set up folders
################################################
#Set up one folder per SL (from list of the 107 SLs) 
while read l ; do mkdir ${l} ; done < list_of_107_overlapping_SLs.txt ;

#For each SL folder, create subfolders for illumina reads, reference and analyses
for f in $(ls -d SL*) ; do mkdir ${f}/illumina_fastq ${f}/reference_gbk ${f}/prep_gubbins ${f}/run_gubbins ${f}/run_raxml ${f}/final_raxml_tree ${f}/figures ${f}/transmission_events ; done ;

#Get the list of genomes per SL
while read l ; do grep -w ${l} list_of_genomes_with_SL.txt  >> ${l}/${l}_list_of_genomes_with_SL.txt ; done < list_of_107_overlapping_SLs.txt ;

#Get only the list of genomes without the SL
for f in $(ls -d SL*) ; do cut -d"," -f1 ${f}/${f}_list_of_genomes_with_SL.txt >> ${f}/list_of_genomes.txt ; done ;

#Symlink all illumina FASTQ files to these folders:
for f in $(ls -d SL*) ; do cd ${f}/illumina_fastq/ ; bash ~/Scripts/symlink_PEfastq_files.sh ../list_of_genomes.txt ./Illumina/FASTQ/symlink_all_FASTQ/ . ; cd ./transmission_analyses/SL_RedDog_alignments/ ; done ;

#Set up reference genomes (from list of reference genome, SL):
mkdir references ; cd references ;
#Copy the fasta file for the chromosome only
while read -r file sl ; do cp ./hybrid_assemblies/Hybrid_fastas__split_by_replicon/${file}__chr.fasta . ; done < <(cat ../list_reference_genomes.txt | sed 's/,/\t/g') 

#Copy the GenBank file 
while read -r file sl ; do cp ./hybrid_assemblies/bakta_annotated_hybrid_assemblies/${file}*/${file}*gbff . ; done < <(cat ../list_reference_genomes.txt | sed 's/,/\t/g') 
while read -r file sl ; do echo $file ; done < <(cat ../list_reference_genomes.txt | sed 's/,/\t/g')  

#Split the genbank file to align only the chromosome locus
for f in *.gbff; do python ~/GitHub/KGWP1_crossniche_paper/transmission_analyses/split_multi_gbk.py < ${f} ; done

while read -r file sl ; do mv ${file}__chr.fasta ${file}__chr.gbk ../${sl}/reference_gbk/ ; done < <(cat ../list_reference_genomes.txt | sed 's/,/\t/g' ) 


################################################
#RedDog alignment
################################################
#Activate environment 
conda activate /home/markus/Programs/reddog-nf/conda_env ;

#Set up reddog dir in each SL folder
git clone https://github.com/scwatts/reddog-nf.git ;
for d in $(ls -d SL*) ; do cp -r reddog-nf ./${d}/ ; done ;

#Edit reddog nextflow config file for each SL
for SL in $(ls -d SL*) ; do echo ${SL} ; cd ${SL}/reference_gbk/ ;  path_gbk=$(ls *gbk ) ; echo ${path_gbk} ; cd ../.. ; cd ${SL}/ ; path=$(pwd) ; sed -i "s|reads = ''|reads = '${path}/illumina_fastq/*fastq.gz'|g" reddog-nf/nextflow.config ; sed -i "s|reference = ''|reference = '${path}/reference_gbk/${path_gbk}'|g" reddog-nf/nextflow.config ; sed -i "s|output_dir = ''|output_dir = '${path}/reddog_output/'|g" reddog-nf/nextflow.config ; cd  ..  ; done ;

#Run RedDog
for SL in $(ls -d SL*) ; do echo ${SL} ; cd ${SL}/reddog-nf/ ; nextflow ./reddog.nf  | tee reddog_log.txt ; cd ./transmission_analyses/SL_RedDog_alignments/merged_lineages_by_verticall/ ; done ;

conda deactivate;

################################################
#Run snp-dists
################################################
#On RedDog output
for SL in $(ls -d SL*); do snp-dists ${SL}/reddog_output/*__chr_cons0.95.mfasta -m -c >> ${SL}/transmission_events/${SL}_reddog_snpdists.csv ; done

conda deactivate ;

#Note: We originally ran gubbins to detect and filter recombinations. We later swapped this with verticall (for downstream dated phylogenies). We did not use the recombination-filtered files for the snp-dists calling, as was more likely to introduce errors if recombinations were incorrectly called. We did run both with/without recombination filtering, and the difference in results between the two was very small. This is also the same method used in Thorpe et al 2022 Nat Micro. 
