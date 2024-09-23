#!/usr/bin/env bash
#################################################
# KGWP1 cross-niche paper
# Commands to:
# - run plasmidfinder via abricate
#################################################

#Installation and setup
mamba create -—name abricate -c conda-forge -c bioconda -c defaults abricate
abricate --check
abricate --list

###Download and set up the plasmidfinder database:
#PlasmidFinder database (v2023-01-18, last commit 2023-03-17, downloaded 2023-08-23)

#1. First, download the latest release of the PlasmidFinder database from: https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src:
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git

#2. Concatenate all the fasta files in the directory to create one plasmidfinder db:
cd plasmidfinder_db ; cat *fsa >> plasmidfinder_db__2023-08-23.fasta

#3. Create a db folder for plasmidfinder in the abricate db:
cd /opt/anaconda/anaconda3/envs/abricate/db
mkdir plasmidfinder_db__2023-08-23 ; cd plasmidfinder_db__2023-08-23/ 
cp plasmidfinder_db__2023-08-23.fasta sequences 

#5. Change the FASTA headers to work with abricate:
#5.1. Add db name + tildes at the beginning of the header:
sed -i 's/>/>plasmidfinder~~~/g' sequences

#5.2. Replace underscores (_) with tildes (~) between gene and accession:
cat sequences | rev | sed 's/_/~~~/' | rev | sed 's/_~~~/~~~/g' >> sequences.temp ; rm sequences; mv sequences.temp sequences 

#5.3. Some accessions (those with underscores in the accession name) need more wrangling:
sed -i 's/__NC~~~/~~~NC_/g' sequences
sed -i 's/_NC~~~/~~~NC_/g' sequences
sed -i 's/_NZ~~~/~~~NZ_/g' sequences

#6. Setup the database so it can be used
abricate --setupdb #or run makeblastdb -in sequences -title plasmidfinder_db__2023-08-23 -dbtype nucl -hash_index
abricate --list

#Versions (as of 2023-08-23)
conda activate abricate ;
abricate --version #v1.0.1
plasmidfinder_db --version #v2023-01-18
conda deactivate ;

#Running
conda activate abricate ;
#Run abricate with plasmidfinder database
parallel --jobs 4 ' abricate --threads 4 --db=plasmidfinder_db__2023-08-23 {} >> {}_abricate.txt ' ::: $(ls *.fasta) ;

#Combine results from all assemblies to one file and remove individual files
cat *_abricate.txt | head -1 >> plFinder_abricate_fullgenes_KG_illumina__2023-08-23.txt ;
cat *_abricate.txt | grep -v "^#FILE" >> plFinder_abricate_fullgenes_KG_illumina__2023-08-23.txt ;
rm *_abricate.txt  ;

#Run abricate summary for pres/abs matrix
abricate --summary plFinder_abricate_fullgenes_KG_illumina__2023-08-23.txt >> plFinder_abricate_genes_KG_illumina__2023-08-23.txt ;
