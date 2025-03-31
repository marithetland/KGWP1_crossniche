---
title: "runPanstripe.Rmd"
author: "MAK Hetland, marit.hetland[at]outlook.com"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  prettydoc::html_pretty:
    theme: cayman
    number_sections: FALSE
    fig_width: 10
    fig_height: 8
---

## Setup knitr
```{r setup, include=TRUE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = '/Users/marit/AMR-SUS Dropbox/Marit Hetland/PhD/Projects/KGWP1/actual_data_for_paper/pangenome/panstripe') 
#R -e "rmarkdown::render('script.Rmd',output_file='output.html')"
```


## Table of contents
- [Package versions](#loadpackages)  
- [Functions](#functions)  
- [Import and wrangle data](#importdata)  
- [Kp1: Run panstripe](#runpanstripeKp1)  
- [Kp1: Compare pangenomes](#comparepangenomesKp1)
- [Kp3: Run panstripe](#runpanstripeKp3)  
- [Kp3: Compare pangenomes](#comparepangenomesKp3)
- [Kp3: Accumulation curves ](#acummulationcurves)

In this script, phylogenies (from IQtree) and Rtabs (from panaroo) are used as input to panstripe. This analysis was done overall, and split by 
Kp1 and Kp3, as shown here. 

## Setup R environment
```{r setup_r_env}
# Set seed
set.seed(1234)
```


## Load packages and show versions {#loadpackages}
```{r load_packages}

# R version
R.Version()$version.string

# Reporting
library(panstripe)
packageVersion("panstripe")

library(ape)
packageVersion("ape")

library(patchwork)
packageVersion("patchwork")

library(readr)
packageVersion("readr")

library(data.table)
packageVersion("data.table")

library(VennDiagram)
packageVersion("VennDiagram")

```

## Import and wrangle data {#importdata}
```{r importdata}
#Load data ----
Kp1_tree <- read.tree("/Users/marit/AMR-SUS\ Dropbox/Marit\ 
Hetland/PhD/Projects/KGWP1/actual_data_for_paper/phylogeny/Kp1/SGH10_CP025080_CP025080.1_alleles_2552strains_var_cons0.95.mfasta.treefile_iqtree.tree") 
#IQtree Kp1 core SNP chr tree SGH10
Kp3_tree <- read.tree("/Users/marit/AMR-SUS\ Dropbox/Marit\ 
Hetland/PhD/Projects/KGWP1/actual_data_for_paper/phylogeny/Kp3/Kv342_CP000964_CP000964.1_alleles_532strains_var_cons0.95_fasttree.tree") #FastTree 
Kp3 core SNP chr tree Kv342

kp1_rtab <- fread("/Users/marit/AMR-SUS Dropbox/Marit 
Hetland/PhD/Projects/KGWP1/actual_data_for_paper/pangenome/panaroo_3255_best_assemblies_alignment/KGWP1__Kp1_gene_presence_absence.Rtab", 
sep="\t") #Read as csv to edit it

kp3_rtab <- fread("/Users/marit/AMR-SUS Dropbox/Marit 
Hetland/PhD/Projects/KGWP1/actual_data_for_paper/pangenome/panaroo_3255_best_assemblies_alignment/KGWP1__Kp3_gene_presence_absence.Rtab", 
sep="\t") #Read as csv to edit it

```

```{r subset_tree}
#Split into one tree per niche ----
#File path for each list to exclude (i.e. all non-Kp1 and those not related to each niche)
file_path_Kp1_nonhuman <- "Kp1_nonKpHUM_list.txt"
file_path_Kp1_nonanimal <- "Kp1_nonKpVET_list.txt"
file_path_Kp1_nonmarine <- "Kp1_nonKpMAR_list.txt"

file_path_Kp3_nonhuman <- "Kp3_nonKpHUM_list.txt"
file_path_Kp3_nonanimal <- "Kp3_nonKpVET_list.txt"
file_path_Kp3_nonmarine <- "Kp3_nonKpMAR_list.txt"

#Set the same tree for each
Kp1_tree_human <- Kp1_tree
Kp1_tree_animal <- Kp1_tree
Kp1_tree_marine <- Kp1_tree

Kp3_tree_human <- Kp3_tree
Kp3_tree_animal <- Kp3_tree
Kp3_tree_marine <- Kp3_tree

#Drop tips in trees for each niche Kp1
non_targeted_tip_names_human <- read_lines(file_path_Kp1_nonhuman)
for (tip in 1:length(non_targeted_tip_names_human)){
  Kp1_tree_human <- drop.tip(Kp1_tree_human, non_targeted_tip_names_human[tip])
}

non_targeted_tip_names_animal <- read_lines(file_path_Kp1_nonanimal)
for (tip in 1:length(non_targeted_tip_names_animal)){
  Kp1_tree_animal <- drop.tip(Kp1_tree_animal, non_targeted_tip_names_animal[tip])
}

non_targeted_tip_names_marine <- read_lines(file_path_Kp1_nonmarine)
for (tip in 1:length(non_targeted_tip_names_marine)){
  Kp1_tree_marine <- drop.tip(Kp1_tree_marine, non_targeted_tip_names_marine[tip])
}

#Drop tips in trees for each niche Kp3
non_targeted_tip_names_human <- read_lines(file_path_Kp3_nonhuman)
for (tip in 1:length(non_targeted_tip_names_human)){
  Kp3_tree_human <- drop.tip(Kp3_tree_human, non_targeted_tip_names_human[tip])
}

non_targeted_tip_names_animal <- read_lines(file_path_Kp3_nonanimal)
for (tip in 1:length(non_targeted_tip_names_animal)){
  Kp3_tree_animal <- drop.tip(Kp3_tree_animal, non_targeted_tip_names_animal[tip])
}

non_targeted_tip_names_marine <- read_lines(file_path_Kp3_nonmarine)
for (tip in 1:length(non_targeted_tip_names_marine)){
  Kp3_tree_marine <- drop.tip(Kp3_tree_marine, non_targeted_tip_names_marine[tip])
}

```

```{r subset_rtabs}
#Split into one Rtab file per niche ----
#Kp1
Kp1_file1_columns <- readLines("Kp1_KpHUM_list.txt")
Kp1_file2_columns <- readLines("Kp1_KGVET_list.txt")
Kp1_file3_columns <- readLines("Kp1_KGMAR_list.txt")

#Kp3
Kp3_file1_columns <- readLines("Kp3_KpHUM_list.txt")
Kp3_file2_columns <- readLines("Kp3_KGVET_list.txt")
Kp3_file3_columns <- readLines("Kp3_KGMAR_list.txt")

# Create data.tables for each file based on the specified column names
file1_dt <- kp1_rtab[, ..Kp1_file1_columns, with = FALSE]
file2_dt <- kp1_rtab[, ..Kp1_file2_columns, with = FALSE]
file3_dt <- kp1_rtab[, ..Kp1_file3_columns, with = FALSE]
file4_dt <- kp3_rtab[, ..Kp3_file1_columns, with = FALSE]
file5_dt <- kp3_rtab[, ..Kp3_file2_columns, with = FALSE]
file6_dt <- kp3_rtab[, ..Kp3_file3_columns, with = FALSE]

# Function to check if all values in a row (except the first column) are 0
all_values_zero <- function(dt) {
  apply(dt[, -1, with = FALSE], 1, function(row) all(row == 0))
}

# Remove rows where all values (except column 1) equal 0
file1_dt <- file1_dt[!all_values_zero(file1_dt)]
file2_dt <- file2_dt[!all_values_zero(file2_dt)]
file3_dt <- file3_dt[!all_values_zero(file3_dt)]
file4_dt <- file4_dt[!all_values_zero(file4_dt)]
file5_dt <- file5_dt[!all_values_zero(file5_dt)]
file6_dt <- file6_dt[!all_values_zero(file6_dt)]

# Save the data.tables to new TSV files
fwrite(file1_dt, "Kp1_human_gene_presence_absence.Rtab", sep="\t")
fwrite(file2_dt, "Kp1_animal_gene_presence_absence.Rtab", sep="\t")
fwrite(file3_dt, "Kp1_marine_gene_presence_absence.Rtab", sep="\t")
fwrite(file4_dt, "Kp3_human_gene_presence_absence.Rtab", sep="\t")
fwrite(file5_dt, "Kp3_animal_gene_presence_absence.Rtab", sep="\t")
fwrite(file6_dt, "Kp3_marine_gene_presence_absence.Rtab", sep="\t")

#Load updated Rtabs to use as input ----
pa_Kp1_human <- read_rtab("Kp1_human_gene_presence_absence.Rtab")
pa_Kp1_animal <- read_rtab("Kp1_animal_gene_presence_absence.Rtab")
pa_Kp1_marine <- read_rtab("Kp1_marine_gene_presence_absence.Rtab")
pa_Kp3_human <- read_rtab("Kp3_human_gene_presence_absence.Rtab")
pa_Kp3_animal <- read_rtab("Kp3_animal_gene_presence_absence.Rtab")
pa_Kp3_marine <- read_rtab("Kp3_marine_gene_presence_absence.Rtab")

```

## Run panstripe on Kp1 data {#runpanstripeKp1}
```{r runpanstripeKp1}
panstripe_Kp1_human <- panstripe(pa_Kp1_human, Kp1_tree_human, family='gaussian')
panstripe_Kp1_animal <- panstripe(pa_Kp1_animal, Kp1_tree_animal, family='gaussian')
panstripe_Kp1_marine <- panstripe(pa_Kp1_marine, Kp1_tree_marine, family='gaussian') 
```

```{r plotpanstripeKp1}
panstripe_Kp1_human$summary
panstripe_Kp1_animal$summary
panstripe_Kp1_marine$summary

plot_pangenome_params(panstripe_Kp1_human)
plot_pangenome_params(panstripe_Kp1_animal)
plot_pangenome_params(panstripe_Kp1_marine)

plot_pangenome_cumulative(panstripe_Kp1_human)
plot_pangenome_cumulative(panstripe_Kp1_animal)
plot_pangenome_cumulative(panstripe_Kp1_marine)

```

## Compare Kp1 pangenomes {#comparepangenomesKp1}
```{r comparepangenomesKp1}
#Compare the fits of the different datasets ----
compare_Kp1_hum_ani <- compare_pangenomes(panstripe_Kp1_human, panstripe_Kp1_animal, family='gaussian')
compare_Kp1_hum_mar <- compare_pangenomes(panstripe_Kp1_human, panstripe_Kp1_marine, family='gaussian')
compare_Kp1_ani_mar <- compare_pangenomes(panstripe_Kp1_animal, panstripe_Kp1_marine, family='gaussian')
```

```{r plotpanstripecomparisonsKp1}
#Look at comparisons
compare_Kp1_hum_ani$summary 
compare_Kp1_hum_mar$summary 
compare_Kp1_ani_mar$summary

#Plot comparisons
plot_Kp1_pangenome_comparison <- plot_pangenome_params(list(a_human = panstripe_Kp1_human, c_animal = panstripe_Kp1_animal, b_marine = 
panstripe_Kp1_marine), legend = FALSE) + 
  plot_pangenome_cumulative(list(a_human = panstripe_Kp1_human, c_animal = panstripe_Kp1_animal, b_marine = panstripe_Kp1_marine)) + 
plot_layout(nrow = 1)

plot_Kp1_pangenome_comparison <- plot_pangenome_cumulative(list(a_human = panstripe_Kp1_human, c_animal = panstripe_Kp1_animal, b_marine = 
panstripe_Kp1_marine)) + plot_layout(nrow = 1)


#Write to PDF
pdf_filename <-  "plot_Kp1_pangenome_comparison.pdf"
pdf(pdf_filename)
plot_Kp1_pangenome_comparison
dev.off()

```

## Run panstripe on Kp3 data {#runpanstripeKp3}
Note: Only human and animal compared as there was only one Kp3 marine genome.
```{r runpanstripeKp3}
panstripe_Kp3_human <- panstripe(pa_Kp3_human, Kp3_tree_human, family='gaussian')
panstripe_Kp3_animal <- panstripe(pa_Kp3_animal, Kp3_tree_animal, family='gaussian')
```

```{r plotpanstripeKp3}
panstripe_Kp3_human$summary
panstripe_Kp3_animal$summary

plot_pangenome_params(panstripe_Kp3_human)
plot_pangenome_params(panstripe_Kp3_animal)

plot_pangenome_cumulative(panstripe_Kp3_human)
plot_pangenome_cumulative(panstripe_Kp3_animal)

```

## Compare Kp3 pangenomes {#comparepangenomesKp3}
```{r comparepangenomesKp3}
#Compare the fits of the different datasets ----
compare_Kp3_hum_ani <- compare_pangenomes(panstripe_Kp3_human, panstripe_Kp3_animal, family='gaussian')

```

```{r plotpanstripecomparisonsKp3}
#Look at comparisons
compare_Kp3_hum_ani$summary 

#Plot comparisons
plot_Kp3_pangenome_comparison <- plot_pangenome_params(list(a_human = panstripe_Kp3_human, b_animal = panstripe_Kp3_animal), legend = FALSE) + 
  plot_pangenome_cumulative(list(a_human = panstripe_Kp3_human, b_animal = panstripe_Kp3_animal)) + plot_layout(nrow = 1)

plot_Kp3_pangenome_comparison <- plot_pangenome_cumulative(list(a_human = panstripe_Kp3_human, b_animal = panstripe_Kp3_animal)) + 
plot_layout(nrow = 1)

pdf_filename <-  "plot_Kp3_pangenome_comparison.pdf"
pdf(pdf_filename)
plot_Kp3_pangenome_comparison
dev.off()


```

## Plot accumulation curves for comparisons {#acummulationcurves}
```{r acummulationcurves}
#This is commented out as it takes a very long time to run
#Plot accumulation curve for comparison
#plot_Kp1_acc_curve <- plot_acc(list(a_human = panstripe_Kp1_human$pa, c_animal = panstripe_Kp1_animal$pa, b_marine = panstripe_Kp1_marine$pa))
#plot_Kp1_acc_curve

#pdf_filename <-  "plot_Kp1_acc_curve.pdf"
#pdf(pdf_filename)
#plot_Kp1_acc_curve
#dev.off()

#plot_Kp3_acc_curve <- plot_acc(list(a_human = panstripe_Kp3_human$pa, #c_animal = panstripe_Kp3_animal$pa))
#plot_Kp3_acc_curve

#pdf_filename <-  "plot_Kp3_acc_curve.pdf"
#pdf(pdf_filename)
#plot_Kp3_acc_curve
#dev.off()
```

