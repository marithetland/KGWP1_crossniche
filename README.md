# KGWP1_crossniche
Code/commands used for analysis and to create figures for the manuscript "A genome-wide One Health study of <i>Klebsiella pneumoniae</i> in Norway reveals overlapping populations but few recent transmission events across reservoirs".

Contains code related to:
* ONT basecalling
  * [ONT basecalling](https://github.com/marithetland/KGWP1_crossniche/blob/main/basecalling/basecalling.sh)
* Assembly (with read trimming and QC)
  * [illumina assembly](https://github.com/marithetland/KGWP1_crossniche/blob/main/assembly/illumina_assembly.sh)
  * [hybrid assembly](https://github.com/marithetland/KGWP1_crossniche/blob/main/assembly/hybrid_assembly.sh)
* Alignments
  * [Kp1 alignment](https://github.com/marithetland/KGWP1_crossniche/blob/main/alignment/kp1_alignment.sh)
  * [Kp3 alignment](https://github.com/marithetland/KGWP1_crossniche/blob/main/alignment/kp3_alignment.sh)
  * [iuc3 alignment](https://github.com/marithetland/KGWP1_crossniche/blob/main/alignment/iuc3_alignment.sh)
  * [colicin alignment](https://github.com/marithetland/KGWP1_crossniche/blob/main/alignment/colicin_alignment.sh)
* Phylogeny
  * [Kp1_phylogeny](https://github.com/marithetland/KGWP1_crossniche/blob/main/phylogeny/kp1_phylogeny.sh)
  * [Kp3_phylogeny](https://github.com/marithetland/KGWP1_crossniche/blob/main/phylogeny/kp3_phylogeny.sh)
* Pangenome
  * [Annotations (Bakta)](https://github.com/marithetland/KGWP1_crossniche/blob/main/pangenome/bakta_annotations.sh)
  * [Pangenome (Panaroo)](https://github.com/marithetland/KGWP1_crossniche/blob/main/pangenome/panaroo.sh)
  * [Pangenome (Panstripe)](https://github.com/marithetland/KGWP1_crossniche/blob/main/pangenome/panstripe.sh)
* Genotyping
  * [Species, ST, capsule, AMR, virulence (Kleborate & Kaptive)](https://github.com/marithetland/KGWP1_crossniche/blob/main/genotyping/kleborate_kaptive.sh)
  * [Plasmid replicons (PlasmidFinder)](https://github.com/marithetland/KGWP1_crossniche/blob/main/genotyping/plasmidfinder_abricate.sh)
  * [Thermo- and heavy metal resistance genes (AMRFinder plus)](https://github.com/marithetland/KGWP1_crossniche/blob/main/genotyping/hmrgs_ncbiamrfinder.sh)
* GWAS
  * [population structure matrix](https://github.com/marithetland/KGWP1_crossniche/blob/main/gwas/get_popstructure_matrix.sh)
  * [pyseer](https://github.com/marithetland/KGWP1_crossniche/blob/main/gwas/run_pyseer.sh)
* Strain-sharing
  * [pipeline to get SL specific SNP distances](https://github.com/marithetland/KGWP1_crossniche/blob/main/transmission_analyses/pipeline_aln_to_network.sh)
  * [script to plot and count transmission events](https://github.com/marithetland/KGWP1_crossniche/blob/main/transmission_analyses/transmission_events__plot_and_count.R)

* [Figures](https://github.com/marithetland/KGWP1_crossniche/blob/main/figures/)
