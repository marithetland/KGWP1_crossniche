#!/usr/bin/env bash
#################################################
# KGWP1 cross-niche paper
# Commands to:
# - run panaroo pangenome estimation
#################################################

#Installation
conda create -n panaroo_v133_env 
conda activate panaroo_v133_env
conda install -c conda-forge -c bioconda -c defaults 'panaroo>=1.3'

#Versions
conda activate panaroo_v133_env
panaroo --version #panaroo 1.3.3

#Running
conda activate panaroo_v133_env ;

#Running without core gene alignment
panaroo -i gff3_paths.txt -o ./panaroo_no_alignment/ --clean-mode strict -t 40 --remove-invalid-genes | tee ./panaroo_no_alignment/panaroo_log.txt
#Was run 1) for all genomes and 2) only per niche (human, animal, marine) by KpSC.