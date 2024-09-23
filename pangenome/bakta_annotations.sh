#!/usr/bin/env bash
#################################################
# KGWP1 cross-niche paper
# Commands to:
# - run bakta annotations
#################################################

#Installation and database download
mamba create -n bakta_v18_env -c conda-forge -c bioconda bakta
conda activate bakta_v18_env

bakta_db download --output /opt/anaconda/anaconda3/envs/bakta_v18_env/lib/python3.10/site-packages/bakta/db

#Versions
conda activate bakta_v18_env ;
bakta --version #v1.8.1
bakta_db list #Required database schema version: 5

#Running
conda activate bakta_v18_env ;

#Running on Illumina-only assemblies:
parallel --jobs 5 'echo {} ; bakta --db /opt/anaconda/anaconda3/envs/bakta_v18_env/lib/python3.10/site-packages/bakta/db/db --threads 4 --prefix {} --output {}_bakta -v {}.fasta ' ::: $(ls *.fasta | sed 's/.fasta//g')   

#Running on draft hybrid genomes:
parallel --jobs 5 'echo {} ; bakta --db /opt/anaconda/anaconda3/envs/bakta_v18_env/lib/python3.10/site-packages/bakta/db/db --replicons unclosed_hybrid_data.tsv --threads 4 --prefix {} --output {}_bakta -v {}.fasta ' ::: $(ls *.fasta | sed 's/.fasta//g') 

#Running on closed hybrid genomes:
parallel --jobs 5 'echo {} ; bakta --db /opt/anaconda/anaconda3/envs/bakta_v18_env/lib/python3.10/site-packages/bakta/db/db --complete --replicons closed_replicon_data.tsv --threads 4 --prefix {} --output {}_bakta -v {}.fasta ' ::: $(ls *.fasta | sed 's/.fasta//g') 