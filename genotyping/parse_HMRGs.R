##########################################
# Heavy metal and thermoresistance genes
# 2024-03-05
# Winkler MA
##########################################
#To adapt to other datasets, swap input on line 15 and if necessary change the defined gene prefixes on lines 26/28

# Load pkgs ----
library(tidyverse)
library(stringr)

# Load data ----
#Genomes were annotated using bakta_annotations.sh, and filtered for HRMGs with get_HMRGs_from_bakta.sh:
bakta_KGWP1 <- read.csv("./external_datasets/KGWP1_hmrgs.txt", sep = "\t", check.names = FALSE)

# Create wide pres/abs table of all genes ----
hmrgs_wide <- bakta_KGWP1 %>%
  select(strain, Gene) %>%
  pivot_wider(names_from = "Gene", values_from = "Gene", values_fill = 0, values_fn = function(x) 1) %>% #When multiple values exist for a single strain-Gene combo, this says: just return 1, effectively making a 1/0 pres/abs matrix.
  select(sort(names(.))) %>% #Sort column names alphabetically
  relocate(strain, .before = arsA) #Move strain back to be the first column

# Add new columns for each operon based on the present/absence of specified genes ---- 
# Define a list of gene prefixes
sort(unique(substr(bakta_KGWP1$Gene, 1, 3))) 
#Select which ones to keep
gene_prefixes <- c("ars", "cad", "chr", "clp", "cop", "cue", "cus", "czc", "fie", "hsp", "kil", "kla", "mer", "ncc", 
                   "ncr", "nir", "pbr", "pco", "rcn", "sil", "teh", "tel", "ter", "zit", "znt", "clpK", "hsp")
# Function to filter only genes to keep
for (prefix in gene_prefixes) {
  
  hmrgs_wide <- hmrgs_wide %>% 
    mutate(!!paste0(prefix, "_genes") :=  
             ifelse(
               rowSums(select(., starts_with({{prefix}})) == "1") > 0,
               apply(select(., starts_with({{prefix}})) == "1", 1,  
                     function(x) paste(colnames(select(., starts_with({{prefix}})))[x], 
                                       collapse = ", ")),
               "0"))
  
}

# Create pres/abs of operons ----
# HMRG gene and operons were defined as described in the Supplementary Methods in: https://www.biorxiv.org/content/10.1101/2024.09.11.612360v1

#Function to check if operons are present
has_all_genes <- function(gene_list, gene_column) {
  Reduce(`&`, lapply(gene_list, function(g) grepl(g, gene_column)))
}

#Run function per operon, store as 1/0 pres/abs
hmrg_operons <- hmrgs_wide %>%
  select(strain, ars_genes:clpK_genes) %>%
  mutate(
    arsABCDR__op = as.numeric(has_all_genes(c("arsA", "arsB", "arsC", "arsD", "arsR"), ars_genes)),
    cadABC__op = as.numeric(has_all_genes(c("cadA", "cadB", "cadC"), cad_genes)),
    chrA__op = as.numeric(has_all_genes(c("chrA"), chr_genes)),
    chrAB__op = as.numeric(has_all_genes(c("chrA", "chrB"), chr_genes)),
    copABCD__op = as.numeric(has_all_genes(c("copA", "copB", "copC", "copD"), cop_genes)),
    cueOR__op = as.numeric(has_all_genes(c("cueO", "cueR"), cue_genes)),
    cusABCF__op = as.numeric(has_all_genes(c("cusA", "cusB", "cusC", "cusF"), cus_genes)),
    czcABC__op = as.numeric(has_all_genes(c("czcA", "czcB", "czcC"), czc_genes)),
    merACP_AFP_APT__op = as.numeric(
      has_all_genes(c("merA", "merC", "merP"), mer_genes) |
        has_all_genes(c("merA", "merF", "merP"), mer_genes) |
        has_all_genes(c("merA", "merP", "merT"), mer_genes)),
    ncrABC__op = as.numeric(has_all_genes(c("ncrA", "ncrB", "ncrC"), ncr_genes)),
    pcoABCDRS__op = as.numeric(has_all_genes(c("pcoA", "pcoB", "pcoC", "pcoD", "pcoR", "pcoS"), pco_genes)),
    rcnAR__op = as.numeric(has_all_genes(c("rcnA", "rcnR"), rcn_genes)),
    silABCERS__op = as.numeric(has_all_genes(c("silA", "silB", "silC", "silE", "silR", "silS"), sil_genes)),
    terBCDE__op = as.numeric(has_all_genes(c("terB", "terC", "terD", "terE"), ter_genes)),
    kilA__op = as.numeric(has_all_genes(c("kilA"), kil_genes)),
    klaABC__op = as.numeric(has_all_genes(c("klaA", "klaB", "klaC"), kla_genes)),
    telAB__op = as.numeric(has_all_genes(c("telA", "telB"), tel_genes)),
    tehAB__op = as.numeric(has_all_genes(c("tehA", "tehB"), teh_genes)),
    zitB__op = as.numeric(has_all_genes(c("zitB"), zit_genes)),
    zntAR__op = as.numeric(has_all_genes(c("zntA", "zntR"), znt_genes)),
    clpK__op = as.numeric(has_all_genes(c("clpK"), clp_genes)),
    hsp20__op = as.numeric(has_all_genes(c("hsp20"), hsp_genes))  ) %>%
  relocate(arsABCDR__op, cadABC__op, chrA__op, chrAB__op, copABCD__op, cueOR__op,
           cusABCF__op, czcABC__op, merACP_AFP_APT__op, ncrABC__op, pcoABCDRS__op,
           rcnAR__op, silABCERS__op, terBCDE__op, kilA__op, klaABC__op, telAB__op,
           tehAB__op, zitB__op, zntAR__op, clpK__op, hsp20__op,
           .before = ars_genes) %>%
  mutate(num_HMRGs = rowSums(across(where(is.numeric)))) 

hmrg_operons

#write.csv(x = hmrg_operons, file = "./external_datasets/extData_HMRGs.tsv", row.names = FALSE)
