#!/usr/bin/env bash
################################################
# KGWP1 cross-niche paper
# Commands to:
# - basecall fast5 files from ONT runs
################################################

#Set up environment
mamba create -n ontbasecalling_env -c bioconda filtlong ;
conda activate ontbasecalling_env ;
pip install h5py ;

#Set path for guppy_basecaller and fast_count 
conda activate ontbasecalling_env ;
export PATH=~/Programs/ont-guppy_v6.4.2/bin:~/Programs/MinION-desktop:$PATH ;

#check versions
guppy_basecaller #v6.4.2+97a7f06
filtlong #v0.2.1

################################################
#ONT basecalling 
################################################
#Loop over each run folder
cd ./ONT/run_folders/ ;
for run_folder in $(ls -d 2*) ; do
cd ${run_folder} ;
#Run basecalling script (from https://github.com/marithetland/ont-basecalling/) 

#For R9 w ligation
python ~/Programs/ont-basecalling/ont-basecalling.py --input_dir . --basecalling_model r9.4.1_sup --barcode_kit native_1-96 --key_file barcode_sample_key.csv   >> ${run_folder}.log 2>&1 ;
cd ./ONT/run_folders/ ;
touch ${run_folder}__basecalling_complete.txt ;
done ;

#For R10 w rapid
python ~/Programs/ont-basecalling/ont-basecalling.py --input_dir . --basecalling_model r10.4.1_sup --barcode_kit rapid_1-24 --key_file barcode_sample_key.csv   >> ${run_folder}.log 2>&1 ;
cd ./ONT/run_folders/ ;
touch ${run_folder}__basecalling_complete.txt ;
done ;