#!/usr/bin/env bash
#################################################
# KGWP1 cross-niche paper
# Commands to:
# - run pyseer w all genetic features and 
# phenotypes
#################################################

#Activate environment
conda activate pyseer_env
pyseer --version #1.3.11

cd ./GWAS/run_dir 

# Define the dataset as a multi-line string
dataset=$(cat <<EOF
BSIvGC,Kp1,BSI,COG,running_BSIvGC_Kp1_BSI_COG,BSIvGC_Kp1_BSI_COG,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,bsi_Kp1.pheno,BSIvGC_Kp1_gene_presence_absence.Rtab
BSIvGC,Kp1,GC,COG,running_BSIvGC_Kp1_GC_COG,BSIvGC_Kp1_GC_COG,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,gc_Kp1.pheno,BSIvGC_Kp1_gene_presence_absence.Rtab
KGWP1,Kp1,animal,COG,running_KGWP1_Kp1_animal_COG,KGWP1_Kp1_animal_COG,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_animal.pheno,KGWP1_Kp1_gene_presence_absence.Rtab
KGWP1,Kp1,human,COG,running_KGWP1_Kp1_human_COG,KGWP1_Kp1_human_COG,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_human.pheno,KGWP1_Kp1_gene_presence_absence.Rtab
BSIvGC,Kp1,BSI,sCOG,running_BSIvGC_Kp1_BSI_sCOG,BSIvGC_Kp1_BSI_sCOG,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,bsi_Kp1.pheno,BSIvGC_Kp1_struct_presence_absence.Rtab
BSIvGC,Kp1,GC,sCOG,running_BSIvGC_Kp1_GC_sCOG,BSIvGC_Kp1_GC_sCOG,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,gc_Kp1.pheno,BSIvGC_Kp1_struct_presence_absence.Rtab
KGWP1,Kp1,animal,sCOG,running_KGWP1_Kp1_animal_sCOG,KGWP1_Kp1_animal_sCOG,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_animal.pheno,KGWP1_Kp1_struct_presence_absence.Rtab
KGWP1,Kp1,human,sCOG,running_KGWP1_Kp1_human_sCOG,KGWP1_Kp1_human_sCOG,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_human.pheno,KGWP1_Kp1_struct_presence_absence.Rtab
BSIvGC,Kp1,BSI,unitig,running_BSIvGC_Kp1_BSI_unitig,BSIvGC_Kp1_BSI_unitig,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,bsi_Kp1.pheno,BSIvGC_Kp1_unitigcaller.pyseer.gz
BSIvGC,Kp1,GC,unitig,running_BSIvGC_Kp1_GC_unitig,BSIvGC_Kp1_GC_unitig,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,gc_Kp1.pheno,BSIvGC_Kp1_unitigcaller.pyseer.gz
KGWP1,Kp1,animal,unitig,running_KGWP1_Kp1_animal_unitig,KGWP1_Kp1_animal_unitig,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_animal.pheno,KGWP1_Kp1_unitigcaller.pyseer.gz
KGWP1,Kp1,human,unitig,running_KGWP1_Kp1_human_unitig,KGWP1_Kp1_human_unitig,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_human.pheno,KGWP1_Kp1_unitigcaller.pyseer.gz
BSIvGC,Kp1,BSI,SNP,running_BSIvGC_Kp1_BSI_SNP,BSIvGC_Kp1_BSI_SNP,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,bsi_Kp1.pheno,BSIvGC_Kp1_SNPs.vcf.gz
BSIvGC,Kp1,GC,SNP,running_BSIvGC_Kp1_GC_SNP,BSIvGC_Kp1_GC_SNP,BSIvGC_Kp1_iqtree_phylogeny_similarity.tsv,gc_Kp1.pheno,BSIvGC_Kp1_SNPs.vcf.gz
KGWP1,Kp1,animal,SNP,running_KGWP1_Kp1_animal_SNP,KGWP1_Kp1_animal_SNP,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_animal.pheno,KGWP1_Kp1_SNPs.vcf.gz
KGWP1,Kp1,human,SNP,running_KGWP1_Kp1_human_SNP,KGWP1_Kp1_human_SNP,KGWP1_Kp1_iqtree_phylogeny_similarity.tsv,KGWP1_Kp1_human.pheno,KGWP1_Kp1_SNPs.vcf.gz
EOF
)

# Iterate over each line in the dataset
echo "$dataset" | while IFS=, read -r project kp_type source_focus genetic_feature echo_text output_prefix similarity phenotype genetic_feature_file
do
    # Determine the --pres, --vcf, or --kmers parameter based on the genetic feature
    case $genetic_feature in
        COG|sCOG)
            genetic_flag="--pres"
            ;;
        SNP)
            genetic_flag="--vcf"
            ;;
        unitig)
            genetic_flag="--kmers"
            ;;
        *)
            echo "Unknown genetic feature: $genetic_feature"
            continue
            ;;
    esac

    #echo "$echo_text"
    echo $similarity
    echo $phenotype
    echo $genetic_flag
    echo $genetic_feature_file
    echo $output_prefix
    pyseer --lmm --similarity "$similarity" --phenotypes "$phenotype" $genetic_flag "$genetic_feature_file" --min-af 0.05 --max-af 0.95 --cpu 80 --output-patterns "${output_prefix}_patterns.tsv" --save-lmm "${output_prefix}_LMM" >> "${output_prefix}_output.tsv"
done

