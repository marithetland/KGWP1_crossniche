#!/usr/bin/env bash
#################################################
# KGWP1 cross-niche paper
# Commands to:
# - install and run Kleborate w Kaptive
#################################################

#Install Kleborate
mamba create -n kleborate -c bioconda kleborate
cd ~/Programs/
git clone --recursive https://github.com/katholt/Kleborate.git
cd Kleborate/kaptive
git pull https://github.com/katholt/Kaptive master
cd ..
python3 setup.py install
kleborate -h

#Updating MLST scheme
cd Kleborate/scripts
python3 getmlst.py --species "Klebsiella pneumoniae"
mv Klebsiella_pneumoniae.fasta ../kleborate/data
mv profiles_csv ../kleborate/data/kpneumoniae.txt
rm -f ../kleborate/data/Klebsiella_pneumoniae.fasta.n*

#Versions
conda activate kleborate ; 
python ~/Programs/Kleborate/kleborate-runner.py --version #v2.4.0
python ~/Programs/Kleborate/kaptive/kaptive.py --version #v2.0.8
blastn -version #2.13.0+

#Run Kleborate
conda activate kleborate ; 
python ~/Programs/Kleborate/kleborate-runner.py --all -a *fasta -o Kleborate_results.txt
conda deactivate ;