#########################
# KGWP1 cross-niche paper
# 2023-03-30, Hetland MAK
# Figure S15A: SL clusters
#########################

#Load pkgs ----
library(tidyverse)
library(igraph)
library(ggnetwork)
library(ggplot2)
library(reshape2)

#set env ----
set.seed(123)
setwd("~/KGWP1/genotyping/")

subtype_colours = c("Human_infection"="#CC004C","Community_carriage"="#a38b46",
                    "Turkey"="#6460AA", "Broiler"="#0089D0", "Pig"="#FCD7DE",
                    "Dog" = "#0DB14B", "Bivalves" = "#00008B")

#Load and prep input files ----
metadata_all <- read.csv("metadata_all_2024-05-29.csv", sep=",", check.names = FALSE)

clstr_assignment <- read.csv("../22SNP_cluster_assignment.csv") %>%
  mutate(clusterSL = paste0(SL,"__",Cluster)) %>% select(clusterSL,Genome) 

data_long <- merge(clstr_assignment,metadata_all, by.x = "Genome", by.y= "strain") %>% group_by(clusterSL,source) %>% count()  

total_counts <- data_long %>%
  group_by(clusterSL) %>%
  summarise(Total = sum(n)) %>%
  arrange(desc(Total))

data_long$clusterSL <- factor(data_long$clusterSL, levels = rev(total_counts$clusterSL))
data_long$source <- factor(data_long$source, levels = c("Human_infection", "Community_carriage", "Dog","Pig","Turkey","Broiler","Bivalves","Seawater"))

#Plot FigS15A ----
SL_dist <- data_long %>% 
  pivot_wider(names_from=source,values_from=n,values_fill=0) %>%
  mutate(total=Human_infection+Community_carriage+Dog+Pig+Turkey+Broiler+Bivalves) %>%
  filter(total >1) %>% 
  mutate(tally = ifelse(total <3, "less","more")) %>% 
  filter(!(tally == "less" &(Pig ==2 | Human_infection ==2 | Broiler ==2 | Community_carriage ==2 | Turkey ==2 | Bivalves ==2 | Dog ==2) )) %>%
  #  filter(total >2) %>%
  pivot_longer(cols=Pig:Dog, names_to="source",values_to="n")  %>%
  ggplot(aes(x = clusterSL, y = n, fill = source)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_manual(values=subtype_colours) +
  labs(title = "Distribution of Genomes Across Niches by Cluster",
       x = "Cluster",
       y = "Count of Genomes",
       fill = "Niche") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()

SL_dist

ggsave("~/KGWP1/figures/FigS15A_SL_dist_barplot.pdf", plot=SL_dist, dpi = 300, width = 18, height = 18, units = "cm")


#Plot FigS15B ----
#Define SNP threshold (change if needed)
snp_threshold <- 22

#Define paths
path <- "~/KGWP1/transmission_analyses/run_transmission_comparison/"
snp_dist_files <- list.files(path, pattern = "_reddog_snpdists.csv$", full.names = TRUE)

#Load metadata 
metadata_all <- read.csv("metadata_all_2024-05-29.csv", sep=",", check.names = FALSE)
metadata <- metadata_all

# --- Create output folders ---
dir.create(file.path(path, "figures"), showWarnings = FALSE)
dir.create(file.path(path, "transmission_counts"), showWarnings = FALSE)

# --- Prepare output summary tables ---
cluster_summary <- tibble()
source_event_summary <- tibble()

# --- Loop through all SL files ---
for (file in snp_dist_files) {
  SL <- str_extract(basename(file), "SL[0-9]+")
  message("Processing ", SL)
  
  # --- Load SNP distances ---
  snp_dists <- read.csv(file, header = TRUE)
  colnames(snp_dists) <- c("to", "from", "snp_dist")
  
  # --- Filter for SNPs <= threshold ---
  edges <- snp_dists %>%
    filter(snp_dist <= snp_threshold, from != to) %>%
    rowwise() %>%
    mutate(pair_id = paste(sort(c(from, to)), collapse = "_")) %>%
    ungroup() %>%
    distinct(pair_id, .keep_all = TRUE) %>%
    mutate(snp_weight = 1)
  
  # --- Skip if no edges ---
  if (nrow(edges) == 0) {
    message("  No transmission links (≤", snp_threshold, " SNPs) for ", SL)
    next
  }
  
  # --- Prepare node list ---
  genomes_in_edges <- unique(c(edges$from, edges$to))
  nodes <- metadata %>%
    filter(strain %in% genomes_in_edges) %>%
    select(strain, Subtype, Niche)
  
  # --- Build network ---
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # --- Plot network ---
  g_plot <- ggnetwork(g)
  plot_file <- file.path(path, "figures", paste0(SL, "_network_", snp_threshold, "SNPs.pdf"))
  
  p <- ggplot(g_plot, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size = snp_weight), color = "black") +
    geom_nodes(aes(color = Subtype, shape = Niche), size = 6) +
    scale_color_manual(values = subtype_colours) +
    theme_void() +
    theme(legend.position = "right") +
    ggtitle(paste0(SL, " Transmission Network (≤", snp_threshold, " SNPs)")) +
    guides(size = "none")
  
  ggsave(plot_file, plot = p, width = 10, height = 7)
  
  # --- Identify clusters with >1 genome ---
  clusters <- components(g)$membership
  cluster_df <- tibble(Genome = names(clusters), cluster = clusters)
  cluster_sizes <- cluster_df %>% count(cluster) %>% rename(cluster_size = n)
  cluster_df <- cluster_df %>%
    left_join(cluster_sizes, by = "cluster") %>%
    filter(cluster_size > 1) %>%
    select(-cluster_size)
  
  num_genomes <- n_distinct(cluster_df$Genome)
  cluster_summary <- bind_rows(cluster_summary, tibble(SL = SL, Num_genomes_in_clusters = num_genomes))
  
  # --- Source-type transmission event summary (deduplicated per cluster) ---
  edges_in_clusters <- edges %>%
    filter(from %in% cluster_df$Genome & to %in% cluster_df$Genome)
  
  edges_with_clusters <- edges_in_clusters %>%
    left_join(cluster_df %>% rename(from = Genome), by = "from") %>%
    rename(cluster_from = cluster) %>%
    left_join(cluster_df %>% rename(to = Genome), by = "to") %>%
    rename(cluster_to = cluster) %>%
    filter(cluster_from == cluster_to) %>%
    rename(cluster = cluster_from)
  
  edges_with_sources <- edges_with_clusters %>%
    left_join(metadata, by = c("from" = "strain")) %>%
    rename(from_source = Subtype) %>%
    left_join(metadata, by = c("to" = "strain")) %>%
    rename(to_source = Subtype) %>%
    rowwise() %>%
    mutate(source_pair = paste(sort(c(from_source, to_source)), collapse = "-")) %>%
    ungroup() %>%
    distinct(cluster, source_pair, .keep_all = TRUE) %>%
    count(source_pair, name = "Count") %>%
    mutate(SL = SL)
  
  source_event_summary <- bind_rows(source_event_summary, edges_with_sources)
}

# --- Save outputs ---
write.csv(cluster_summary, file.path(path, "num_genomes_in_clusters_by_SL.csv"), row.names = FALSE)
write.csv(source_event_summary, file.path(path, "source_transmission_events_by_SL.csv"), row.names = FALSE)
