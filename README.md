# KGWP1_crossniche
Code/commands used for analysis and to create figures for: 

Hetland MAK, Winkler MA, Kaspersen HP, et al. A genome-wide One Health study of <i>Klebsiella pneumoniae</i> in Norway reveals overlapping populations but few recent transmission events across reservoirs. Genome Medicine, 2025. https://doi.org/10.1186/s13073-025-01466-0

Contains code related to:
* ONT basecalling
  * [ONT basecalling](https://github.com/marithetland/KGWP1_crossniche/blob/main/basecalling/ont_basecalling.sh)
* Assembly
  * [Illumina assembly](https://github.com/marithetland/KGWP1_crossniche/blob/main/assembly/illumina_assembly.sh)
  * [Hybrid assembly](https://github.com/marithetland/KGWP1_crossniche/blob/main/assembly/hybrid_assembly.sh)
* Phylogeny
  * [IQtree](https://github.com/marithetland/KGWP1_crossniche/blob/main/phylogeny/run_iqtree.sh)
* Pangenome
  * [Annotations (Bakta)](https://github.com/marithetland/KGWP1_crossniche/blob/main/pangenome/bakta_annotations.sh)
  * [Pangenome (Panaroo)](https://github.com/marithetland/KGWP1_crossniche/blob/main/pangenome/panaroo.sh)
  * [Pangenome (Panstripe)](https://github.com/marithetland/KGWP1_crossniche/blob/main/pangenome/panstripe.sh)
* Genotyping
  * [Species, ST, capsule, AMR, virulence (Kleborate & Kaptive)](https://github.com/marithetland/KGWP1_crossniche/blob/main/genotyping/Kleborate_kaptive.sh)
  * [Plasmid replicons (PlasmidFinder)](https://github.com/marithetland/KGWP1_crossniche/blob/main/genotyping/plasmidfinder_abricate.sh)
  * [Thermo- and heavy metal resistance genes (AMRFinder plus)](https://github.com/marithetland/KGWP1_crossniche/blob/main/genotyping/hmrgs_ncbiamrfinder.sh)
* GWAS
  * [Population structure matrix](https://github.com/marithetland/KGWP1_crossniche/blob/main/gwas/get_popstructure_matrix.sh)
  * [Pyseer](https://github.com/marithetland/KGWP1_crossniche/blob/main/gwas/run_pyseer.sh)
* Strain-sharing
  * [SL specific SNP distances](https://github.com/marithetland/KGWP1_crossniche/blob/main/transmission_analyses/SL_alignments.sh)
  * [Filter recombinations with Verticall](https://github.com/marithetland/KGWP1_crossniche/blob/main/transmission_analyses/run_verticall.sh)
  * [Root2Tip regression](https://github.com/marithetland/KGWP1_crossniche/blob/main/transmission_analyses/roo2tip.Rmd)
  * [BactDating global](https://github.com/marithetland/KGWP1_crossniche/blob/main/transmission_analyses/BactDating_global.Rmd)
  * [BactDating local](https://github.com/marithetland/KGWP1_crossniche/blob/main/transmission_analyses/BactDating_local.Rmd)


* [Figures](https://github.com/marithetland/KGWP1_crossniche/blob/main/figures/)
