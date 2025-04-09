#KGWP1: R script for exploring combinations of clinical features
#2024-02-02

#Set env
set.seed(123) 

#Load pkgs
library(dplyr)
library(tidyverse)
library(pheatmap)

#Read data, select cols for heatmap
metadata_all <- read.csv("~/KGWP1/genotyping/metadata_all.csv", sep="\t", check.names = FALSE)

#Define colours
custom_colors <- c("white","blue" ,"blue", "red", "green", "#0DB14B","purple","black")
names(custom_colors) <- c(0, 0.3, 0.4, 0.7, 1)
niche_colours <- c("Animal"="#228B22","Human"="#CD2626","Marine"="#1874CD")
niche_colours_heatmap <- list(Niche = c("Marine" = "#1874CD","Animal" = "#228B22","Human" = "#CD2626"))

#Plot a heatmap of all clin features together ----

#First, select which plasmidfinder cols to keep based on num presence
pl_cols <- metadata_all %>% 
  select(c(Niche,`Col.BS512._1`:repB_KLEB_VIR)) %>%
  #select(-matches(c("_c_op$","_incom_op$","_incomp_op$"))) %>% #Remove unwanted plasmid columns
  pivot_longer(cols= `Col.BS512._1`:repB_KLEB_VIR, values_to = "n") %>% 
  filter(n>0) %>% #Select only presence rows for counting
  group_by(name) %>% #group_by(Niche,name)
  count() %>% arrange() %>%
  filter(n>20) %>% #Only those present in more than n genomes
  select(name)
col_names <- pl_cols$name

#Then set the rest of the columns and merge
k <- metadata_all %>% 
  #filter(Niche.x == "Marine") %>% #Uncomment to create a heatmap for only one niche
  select(c(Niche,Yersiniabactin,Colibactin,Aerobactin,Salmochelin,all_of(col_names),AGly_acquired:Bla_Carb_acquired,arsABCDR__op:hsp20__op)) %>%  #selects virulence, selected pl cols, amr and hmrg
  mutate(across(c(Yersiniabactin:Salmochelin,AGly_acquired:Bla_Carb_acquired), ~ if_else(. == "-", 0, 1))) %>% #Change to 0/1 matrix
  select(-matches(c("_c_op$","_incom_op$","_incomp_op$"))) # #Remove unwanted columns


# Define column categories
virulence_genes <- c("Yersiniabactin", "Colibactin", "Aerobactin", "Salmochelin")
plasmid_genes <- c("Col.MG828._1", "Col.Ye4449._1", "Col.pHAD28._1", "Col156_1", 
                   "Col3M_1", "Col440II_1", "Col440I_1", "Col8282_1", "ColRNAI_1", 
                   "ColpVC_1", "IncC_1", "IncFIA.HI1._1", "IncFIA_1", "IncFIB.AP001918._1", 
                   "IncFIB.K..pCAV1099.114._1", "IncFIB.K._1", "IncFIB.pB171._1", 
                   "IncFIB.pKPHS1._1", "IncFIB.pNDM.Mar._1", "IncFIB.pQil._1", "IncFIC.FII._1", 
                   "IncFII.Yp._1", "IncFII.pAMA1167.NDM.5._1", "IncFII.pCRY._1", "IncFII.pKP91._1", 
                   "IncFII.pRSB107._1", "IncHI1B.pNDM.MAR._1", "IncM1_1", "IncN_1", "IncR_1", 
                   "p0111_1", "pENTAS02_1", "repB.R1701._1", "repB_KLEB_VIR")
amr_genes <- c("AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired", "Gly_acquired", 
               "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired", "Tet_acquired", 
               "Tgc_acquired", "Tmt_acquired", "Bla_acquired", "Bla_inhR_acquired", "Bla_ESBL_acquired", 
               "Bla_ESBL_inhR_acquired", "Bla_Carb_acquired")
hmrg_genes <- c("arsABCDR__op", "cadABCR__op", "chrA__op", "copABCD__op", "cusABCF__op", "czcABC__op", 
                "merACP_AFP_APT__op", "ncrABC__op", "nirABCD__op", "pcoABCDRS__op", "rcnAR__op", 
                "silABCERS__op", "terBCDE__op", "klaABC__op", "telAB__op", "clpK__op", "hsp20__op")


# Combine all categories into a named list
categories <- list(
  Virulence = virulence_genes,
  Plasmid = plasmid_genes,
  AMR = amr_genes,
  HeavyMetal = hmrg_genes
)

# Define a function to map the gene categories to specific values
map_gene_to_value <- function(gene) {
  if (gene %in% virulence_genes) return(0.3)
  if (gene %in% plasmid_genes) return(1.0)
  if (gene %in% amr_genes) return(0.4)
  if (gene %in% hmrg_genes) return(0.7)
  return(0)
}

# Create a data matrix (without the first column (niche))
data_matrix <- as.matrix(k[,-1])

# Ensure all columns are numeric
data_matrix <- apply(data_matrix, 2, as.numeric)

# Check for any remaining NA/NaN/Inf values and handle them
data_matrix[is.na(data_matrix)] <- 0
data_matrix[is.nan(data_matrix)] <- 0
data_matrix[is.infinite(data_matrix)] <- 0

# Remove columns where all values are 0 to clean up fig
data_matrix <- data_matrix[, colSums(data_matrix != 0) > 0]

# Transform the data matrix to represent categories
for (gene in colnames(data_matrix)) {
  data_matrix[, gene] <- data_matrix[, gene] * map_gene_to_value(gene)
}

# Check the matrix again for any non-numeric or problematic values
if (any(is.na(data_matrix)) || any(is.nan(data_matrix)) || any(is.infinite(data_matrix))) {
  stop("Data matrix still contains NA, NaN, or Inf values after cleaning.")
}

# Create annotation data frame
annotation <- data.frame(Niche = k$Niche)
rownames(annotation) <- rownames(k)

# Ensure rownames of data_matrix match rownames of annotation
rownames(data_matrix) <- rownames(k)

# Create a heatmap with hierarchical clustering, custom colors, and annotations
FigS10 <- pheatmap(data_matrix, 
                  scale = "none", 
                  border_color = "white",
                  clustering_distance_rows = "binary",
                  clustering_distance_cols = "binary",
                  clustering_method = "complete",
                  show_rownames = TRUE, 
                  show_colnames = TRUE,
                  color = custom_colors,
                  annotation_row = annotation,
                  annotation_colors = niche_colours_heatmap)
FigS10

ggsave(filename= "~/KGWP1/figures/FigS10.pdf", plot = FigS10, dpi = 300, width = 12, height = 8, units = "in")
