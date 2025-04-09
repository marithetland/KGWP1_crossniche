#########################
# KGWP1 cross-niche paper
# 2023-03-30, Hetland MAK
# Figure 3: SNPs
#########################

#Set env ----
set.seed(123)
library(tidyverse)
library(eulerr)
library(data.table)

setwd("~/KGWP1/pangenome/panaroo_3255_best_assemblies_alignment/")

source_colors <- c("Human_infection"="#CC004C","Community_carriage"="#a38b46",
                   "Turkey"="#6460AA", "Broiler"="#0089D0", "Pig"="#FCD7DE",
                   "Dog" = "#0DB14B", "Bivalves" = "#00008B")

#Load data ----
KpSC_csv <- read.csv("gene_presence_absence.csv")
metadata_all <- read.csv("../../genotyping/metadata_all_2024-05-29.csv")

#Set functions ----
#Fig 3A: Function to calculate gene counts
get_gene_counts_for_pie <- function(file_path, niche_name = "Niche") {
  # Read the input gene presence/absence table
  my_data_table <- fread(file_path, sep="\t")
  
  # Calculate thresholds
  total_isolates <- ncol(my_data_table) - 1
  threshold_5 <- 0.05 * total_isolates
  threshold_95 <- 0.95 * total_isolates
  
  # Calculate gene categories
  counts <- rowSums(my_data_table[, -1])
  gene_5 <- sum(counts < threshold_5)
  gene_between <- sum(counts >= threshold_5 & counts <= threshold_95)
  gene_95 <- sum(counts > threshold_95)
  unique_genes <- length(unique(my_data_table$Gene))
  genes_in_all <- sum(counts == total_isolates)
  
  # Create a named vector of counts
  gene_counts <- c(
    "Genes <5%" = gene_5,
    "Genes 5–95%" = gene_between,
    "Genes >95%" = gene_95
  )
  
  # Create a dataframe for plotting
  pie_df <- data.frame(
    Category = names(gene_counts),
    Count = as.numeric(gene_counts)
  )
  
  # Create the pie chart using ggplot2
  pie_chart <- ggplot(pie_df, aes(x = "", y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(title = paste("Gene Category Distribution (", niche_name, ")", sep=""), x = NULL, y = NULL) +
    theme_void() +
    theme(legend.title = element_blank())
  
  # Return a list of results
  return(list(
    gene_counts = gene_counts,
    unique_genes = unique_genes,
    genes_in_all_isolates = genes_in_all,
    plot = pie_chart
  ))
}

#Jaccard: Function to identify and filter out core genes and to keep only genomes in Kp list
get_acc_genes <- function(df,kp_list) {
  #Identify core genes
  core_genes_list <- df %>%
    mutate(across(-c(Gene, `Non-unique Gene name`, Annotation), ~ if_else(. != "" & !is.na(.), 1, 0))) %>%
    column_to_rownames(var = "Gene") %>%
    filter(rowSums(select(., -c(`Non-unique Gene name`, Annotation))) > 0) %>%
    mutate(proportion = rowSums(select(., -c( `Non-unique Gene name`, Annotation))) / (ncol(.) - 2)) %>%
    filter(proportion > 0.95) %>%
    rownames() 
  
  #Remove core genes
  df_accessory <- df %>% 
    filter(!Gene %in% core_genes_list) %>%
    select(c(Gene, `Non-unique Gene name`, Annotation, all_of(unname(kp_list))))
  
  return(df_accessory)
}

#Jaccard: Function to combine data, calculate statistics and plot and save ggplot
get_acc_stats <- function(name,jaccard_dist,meta) {
  
  #Set output name
  outname=paste0("jaccard_dists/",name,"_accgenome_JacDist.pdf")
  
  #Merge distances with metadata for plotting
  jaccard_dist_meta <- jaccard_dist %>%
    left_join(meta, by = c("genome1" = "strain")) %>%
    rename_with(.fn = ~paste0(., "1"), .cols = -c(genome1, genome2, distance)) %>%
    left_join(meta, by = c("genome2" = "strain")) %>%
    rename_with(.fn = ~paste0(., "2"), .cols = -c(genome1, genome2, distance, ends_with("1")))
  
  #Plot results 
  plot <- jaccard_dist_meta %>%
    filter(source1 == source2) %>%
    mutate(source1 = factor(source1, levels = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler", "Bivalves", "Seawater"))) %>%  # Set the order of x-axis categories
    ggplot(aes(x = source1, y = distance, fill = source1)) + 
    geom_violin(trim = FALSE) +  # Draw violin plot
    stat_summary(fun = median, geom = "point", color = "black", size = 3, shape = 16) +  # Add median point
    labs(x = "", y = "Jaccard Distance") +
    theme_bw() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_fill_manual(values = source_colors) 
  
  
  #Save file
  ggsave(outname, plot = plot, dpi = 300, width = 20, height = 14, units = "cm")
  
  #Get mean and median values
  stats <- jaccard_dist_meta %>%
    group_by(source1) %>%
    summarize(median = median(distance, na.rm = TRUE),
              mean = mean(distance, na.rm = TRUE),
              min = min(distance, na.rm = TRUE),
              max = max(distance, na.rm = TRUE),
              num_pairs = n(),
              .groups = 'drop')
  
  return(list(stats = stats, plot = plot))
}

#Fig 3A Pie charts ----
KpSC_overall <- get_gene_counts_for_pie("gene_presence_absence.Rtab", "Overall")
KpSC_human <- get_gene_counts_for_pie("Kp1_human_gene_presence_absence.Rtab", "Human")
KpSC_animal <- get_gene_counts_for_pie("Kp1_animal_gene_presence_absence.Rtab", "Animal")
KpSC_marine <- get_gene_counts_for_pie("Kp1_marine_gene_presence_absence.Rtab", "Marine")

#Combine gene count results
results_df <- bind_rows(
  data.frame(Niche = "Overall", as.list(KpSC_overall$gene_counts), Unique_Genes = KpSC_overall$unique_genes, Genes_In_All_Isolates = KpSC_overall$genes_in_all_isolates),
  data.frame(Niche = "Human", as.list(KpSC_human$gene_counts), Unique_Genes = KpSC_human$unique_genes, Genes_In_All_Isolates = KpSC_human$genes_in_all_isolates),
  data.frame(Niche = "Animal", as.list(KpSC_animal$gene_counts), Unique_Genes = KpSC_animal$unique_genes, Genes_In_All_Isolates = KpSC_animal$genes_in_all_isolates),
  data.frame(Niche = "Marine", as.list(KpSC_marine$gene_counts), Unique_Genes = KpSC_marine$unique_genes, Genes_In_All_Isolates = KpSC_marine$genes_in_all_isolates)
)

#Plot and save pie charts 
ggsave("~/KGWP1/figures/KpSC_overall_pangenome_pie_2024-09-02.pdf", plot=KpSC_overall$plot, dpi = 300, width = 18, height = 18, units = "cm")
ggsave("~/KGWP1/figures/KpSC_human_pangenome_pie_2024-09-02.pdf", plot=KpSC_human$plot, dpi = 300, width = 18, height = 18, units = "cm")
ggsave("~/KGWP1/figures/KpSC_animal_pangenome_pie_2024-09-02.pdf", plot=KpSC_animal$plot, dpi = 300, width = 18, height = 18, units = "cm")
ggsave("~/KGWP1/figures/KpSC_marine_pangenome_pie_2024-09-02.pdf", plot=KpSC_marine$plot, dpi = 300, width = 18, height = 18, units = "cm")

#For Panstripe p-values: See script pangenom/panstripe.R

#Fig 3B Euler pangenome ----
rtab_long <- KpSC_overall %>% #Loaded in fig 3a chunk
  pivot_longer(-Gene, names_to = "strain", values_to = "presence") %>%
  left_join(metadata_all %>% select(strain, Niche), by = "strain") %>% #add metadata
  filter(presence == 1) %>% #Keep only those that are present
  distinct(Gene, Niche) %>% #Remove if duplicates
  #Summarise gene sets by niche
  group_by(Gene) %>%
  summarise(niche_combo = paste(sort(unique(Niche)), collapse = "&")) %>%
  ungroup() %>%
  count(niche_combo)

fit <- euler(setNames(rtab_long$n, rtab_long$niche_combo))
plot(fit, fills = list(fill = c("#CD2626","#228B22","#1874CD"), alpha = 0.7),
     edges = TRUE, quantities = TRUE, labels = TRUE)

#Fig 3C Jaccard dists ----
#Calculating Jaccard Distances to compare between niches (and sources). 
#Used jaccard_from_panaroo_p3.py to calculate the distances, and used only accessory genes (present in ≤95% of genomes). 

#Get accessory genes only and store
KpSC_acc__csv <- get_acc_genes(KpSC_csv,KpSC_list)

#Remove if any are 0 or 1 across for the whole gene (should be handled in function above)
KpSC_acc__csv <- KpSC_acc__csv %>%
  rowwise() %>%
  filter(!(all(c_across(4:ncol(KpSC_acc__csv)) == first(c_across(4:ncol(KpSC_acc__csv))))))
KpSC_acc__csv <- as.data.frame(KpSC_acc__csv)

write.csv(KpSC_acc__csv,"/jaccard_dists/panaroo_KpSC_acc.csv", quote = FALSE, row.names = FALSE)
# ---> Run jaccard_from_panaroo_p3.py on this output! Then read in output from that script:
KpSC_jaccard_dist <- read.csv("jaccard_dists/panaroo_KpSC_acc.csv__panaroo_jaccarddist.tsv", sep="\t", col.names = c("genome1","genome2","distance"))

#Plot by within-source
sub_metadata_source <- metadata_all %>% select(strain,Subtype) %>% rename(source=Subtype)
KpSC_stats <- get_acc_stats("KpSC",KpSC_jaccard_dist,sub_metadata_source) 
KpSC_stats$plot
KpSC_stats$stats

#Save
ggsave("~/KGWP1/figures/KpSC_Jaccard_within_2024-09-02.pdf", plot=KpSC_stats$plot, dpi = 300, width = 18, height = 18, units = "cm")

#Plot by between-source, comparing all sources to human infection
meta <- metadata_all %>% select(strain,Subtype)

KpSC_jaccard_between_source_plot <- jaccard_dist_meta <- KpSC_jaccard_dist %>%
  left_join(meta, by = c("genome1" = "strain")) %>%
  rename_with(.fn = ~paste0(., "1"), .cols = -c(genome1, genome2, distance)) %>%
  left_join(meta, by = c("genome2" = "strain")) %>%
  rename_with(.fn = ~paste0(., "2"), .cols = -c(genome1, genome2, distance, ends_with("1"))) %>%
  distinct() %>% #Remove duplicate rows
  #Make sure "Human_infection" is in "Subtype1" column
  mutate(temp_genome = ifelse(Subtype1 == "Human_infection", genome1, genome2),
         temp_subtype = ifelse(Subtype1 == "Human_infection", Subtype1, Subtype2),
         genome2 = ifelse(Subtype1 == "Human_infection", genome2, genome1),
         Subtype2 = ifelse(Subtype1 == "Human_infection", Subtype2, Subtype1),
         genome1 = temp_genome,
         Subtype1 = temp_subtype) %>%
  select(-temp_genome, -temp_subtype) %>%
  #Filter for pairs that include "Human_infection"
  filter(Subtype1 == "Human_infection" | Subtype2 == "Human_infection") %>%
  mutate(
    other_subtype = ifelse(Subtype1 == "Human_infection", Subtype2, Subtype1)
  ) %>%
  mutate(other_subtype = factor(other_subtype, levels = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler", "Bivalves", "Seawater"))) %>%  #Set the order of x-axis categories
  #Create the violin plot with color and median
  ggplot( aes(x = other_subtype, y = distance, fill = other_subtype)) + 
  geom_violin(trim = FALSE) +  
  stat_summary(fun = median, geom = "point", color = "black", size = 3, shape = 16) +  
  labs(x = "", y = "Jaccard Distance", title = "Genetic Distances Between Human Infection and Other Subtypes") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = source_colors)

# Save
ggsave("~/KGWP1/figures/KpSC_Jaccard_between_2024-09-02.pdf", plot=KpSC_jaccard_between_source_plot, dpi = 300, width = 18, height = 18, units = "cm")