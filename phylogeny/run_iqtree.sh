#!/bin/bash
#SBATCH --job-name=IQKp1
#SBATCH --partition=m3i,comp
#SBATCH --account=js66
#SBATCH --qos=normal
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=36
#SBATCH --mem=700G

# Job modules
module load python/3.8.5-gcc8
module load iqtree

# Job commands
cd ./KGWP1/run_iqtree/rdkg_KGWP1_Kp1__iqtree ;

iqtree -s SGH10_CP025080.1_alleles_var_cons0.95.mfasta --mem 400G --runs 1 -m TEST -B 1000 -nm 10000 -nt AUTO 