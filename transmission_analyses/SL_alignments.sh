#!/usr/bin/env bash
################################################
# KGWP1 cross-niche paper
# Pipeline to:
# - run alignment
# - run gubbins and filter recombinations
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

#TEMP
#Set up lists for symlinking
while read l ; do grep -w ${l} list_with_SL.txt  >> ${l}/${l}_list_of_genomes_with_SL.txt ; done < list_of_overlapping_SLs.txt 
for f in $(ls -d SL*) ; do cut -f1 ${f}/${f}_list_of_genomes_with_SL.txt >> ${f}/list_of_genomes.txt ; done ;
#Symlink local files
for f in $(ls -d SL*) ; do cd ${f}/illumina_fastq/ ; bash ~/Scripts/symlink_PEfastq_files.sh ../list_of_genomes.txt ./Illumina/FASTQ/symlink_all_FASTQ/ . ; cd ./transmission_analyses/global_BactDate/setup/ ; done ;
#Symlink global files
for f in $(ls -d SL*) ; do cd ${f}/illumina_fastq/ ; bash ~/Scripts/symlink_PEfastq_files.sh ../list_of_genomes.txt ./transmission_analyses/global_BactDate/downloaded_genomes/ . ; cd ./transmission_analyses/global_BactDate/setup/ ; done ;
#Set up reference (first copy reference files into reference_gbk folder)
for SL in $(ls -d SL37 SL45) ; do echo ${SL} ; cd ${SL}/reference_gbk/ ;  path_gbk=$(ls *gbk ) ; echo ${path_gbk} ; cd ../.. ; cd ${SL}/ ; path=$(pwd) ; sed -i "s|reads = ''|reads = '${path}/illumina_fastq/*fastq.gz'|g" reddog-nf/nextflow.config ; sed -i "s|reference = ''|reference = '${path}/reference_gbk/${path_gbk}'|g" reddog-nf/nextflow.config ; sed -i "s|output_dir = ''|output_dir = '${path}/reddog_output/'|g" reddog-nf/nextflow.config ; cd  ..  ; done ;
#Run RedDog and gubbins
conda activate /home/markus/Programs/reddog-nf/conda_env ; cd ./transmission_analyses/global_BactDate/SL17/reddog-nf/ ; nextflow ./reddog.nf -resume --force | tee reddog_resume_log.txt ; rm -rf work ; cd ./transmission_analyses/global_BactDate/ ; for SL in $(ls -d SL* | grep -wv "SL10" | grep -wv "SL107" | grep -wv "SL111" | grep -wv "SL17" ) ; do echo ${SL} ; cd ${SL}/reddog-nf/ ; nextflow ./reddog.nf  | tee reddog_log.txt ; rm -rf work ; cd ./transmission_analyses/global_BactDate/ ; done ; conda deactivate ; cd ./transmission_analyses/global_BactDate/ ; for SL in $(ls -d SL*); do cp ${SL}/reddog_output/*__chr_alleles_cons0.95.csv ${SL}/prep_gubbins/ ; done ; for SL in $(ls -d SL*); do head -1 ${SL}/prep_gubbins/*__chr_alleles_cons0.95.csv  | tr ',' '\n' | sort | grep -v Pos | grep -v Ref >> ${SL}/prep_gubbins/sequence_list.txt ; done ; for SL in $(ls -d SL*); do python /home/markus/Programs/RedDog-1beta.11/parseSNPtable3.py -s ${SL}/prep_gubbins/*__chr_alleles_cons0.95.csv  -l ${SL}/prep_gubbins/sequence_list.txt -m aln -p ${SL}/prep_gubbins/${SL}_crossniche_ | tee ${SL}/prep_gubbins/${SL}_parseSNPtable_log.txt ; done ; for SL in $(ls -d SL*); do python /home/markus/Programs/RedDog-1beta.11/snpTable2GenomeAlignment.py --ref ${SL}/reference_gbk/*__chr.fasta --snps ${SL}/prep_gubbins/*_alleles_cons0.95_*strains_var.csv > ${SL}/run_gubbins/WholeGenomeAlignment.fasta ; done ; conda activate gubbins_env ; for SL in $(ls -d SL*); do cd ${SL}/run_gubbins/ ; run_gubbins.py WholeGenomeAlignment.fasta --iterations 100 --converge-method weighted_robinson_foulds -p ${SL}_crossniche_gubbins --tree-builder raxmlng --model GTRGAMMA --threads 20 | tee ${SL}_gubbins_log.txt ; cd ../.. ; done ; conda deactivate ;

################################################
#Run snp-dists
################################################
#On RedDog output
for SL in $(ls -d SL*); do snp-dists ${SL}/reddog_output/*__chr_cons0.95.mfasta -m -c >> ${SL}/transmission_events/${SL}_reddog_snpdists.csv ; done


#On Gubbins output
for SL in $(ls -d SL*); do snp-dists ${SL}/run_gubbins/${SL}_crossniche_gubbins.filtered_polymorphic_sites.fasta -m -c >> ${SL}/transmission_events/${SL}_recombfree_snpdists.csv ; done

conda deactivate ;
