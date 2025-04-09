#####################################
#KGWP1 cross-niche paper
#Plot pairwise SNP distributions
# Fig 5B and S1A and S1B
#Hetland MAK, 2023-10-30
#####################################


#Load pkgs ----
library(scales)
library(dplyr)
library(ggplot2)
library(stringr)

#Load data ----
#NOTE: We decided to use the 107 SLs data instead of the SGH data as they produced nearly the same results, and as the niche-overlapping SLs were also used downstream for BactDating analyses.

#Load SNP dist file for full 3250 genomes (5 failed read mapping) + SL only
#snpdist_3250 <- read.csv("~/KGWP1/transmission_analyses/tranmission_reconstruction/snps/SGH10_CP025080_CP025080.1_alleles_3255strains_var_cons0.95.mfasta.snpdists.csv", header=FALSE)
snpdist_SLs <- read.csv("~/KGWP1/run_transmission_comparison/SLs_reddog_snpdists.tsv", header=TRUE)

#colnames(snpdist_3250) <- c("strain1","strain2","distance")
colnames(snpdist_SLs) <- c("strain1","strain2","distance")

#Load metadata 
metadata <- read.csv("~/KGWP1/metadata/sourceinfo__2023-10-12.tsv", sep="\t")
metadata <- metadata %>% select(strain,Niche,Subtype) %>% filter(Subtype != "Kres") 
colnames(metadata) <- c("strain","niche","source")

#Load kleborate 
kleborate<- read.csv("~/KGWP1/genotyping/kleborate_v236_kaptive_v207_KG_illumina_2023-07-10.txt", sep="\t")
kleborate <- kleborate %>% select(strain,species) %>% filter(!str_detect(strain, "^Kres"))  %>% filter(!str_detect(strain, "^SAM"))

#Prep data ----
#reduce the SNP-dist file, to remove self-interactions and double rows (where strain1 and strain2 are inverted)
filter_snpdist <- function(df) {
  df %>%
    # Step 1: Filter out rows where strain1 and strain2 are the same
    filter(strain1 != strain2) %>%
    # Step 2: Sort strain1 and strain2 alphabetically in each row
    rowwise() %>%
    mutate(strain1_sorted = min(strain1, strain2),
           strain2_sorted = max(strain1, strain2)) %>%
    # Step 3: Remove duplicated rows
    select(strain1_sorted, strain2_sorted, distance) %>%
    distinct(strain1_sorted, strain2_sorted, .keep_all = TRUE) %>%
    # (Optional) Renaming columns back to original names
    rename(strain1 = strain1_sorted,
           strain2 = strain2_sorted)
}

#snpdist_3250_filt <- filter_snpdist(snpdist_3250)
snpdist_SLs_filt <- filter_snpdist(snpdist_SLs)

#Match snpdists and metadata ----
join_and_rename <- function(df_filt, metadata, kleborate) {
  # Joining with metadata for strain1 and renaming columns
  df_filt <- df_filt %>%
    left_join(metadata, by = c("strain1" = "strain")) %>%
    rename(
      niche1 = niche,
      source1 = source
    )
  
  # Joining with metadata for strain2 and renaming columns
  df_filt <- df_filt %>%
    left_join(metadata, by = c("strain2" = "strain")) %>%
    rename(
      niche2 = niche,
      source2 = source
    )
  
  # Joining with kleborate for strain1 and renaming columns
  df_filt <- df_filt %>%
    left_join(kleborate, by = c("strain1" = "strain")) %>%
    rename(
      species1 = species
    )
  
  # Joining with kleborate for strain2 and renaming columns
  df_filt <- df_filt %>%
    left_join(kleborate, by = c("strain2" = "strain")) %>%
    rename(
      species2 = species
    )
  
  return(df_filt)
}

#snpdist_3250_filt <- join_and_rename(snpdist_3250_filt, metadata, kleborate)
snpdist_SLs_filt <- join_and_rename(snpdist_SLs_filt, metadata, kleborate)

#Set as tibble
#snpdist_3250_filt_pt <- as.tibble(snpdist_3250_filt)
snpdist_SLs_filt_pt <- as.tibble(snpdist_SLs_filt)

#Set columns for niche and source
# snpdist_3250_filt_pt <- snpdist_3250_filt_pt %>%
#   mutate(same_niche = if_else(niche1==niche2, 1, 0)) %>%
#   mutate(same_source = if_else(source1==source2, 1, 0)) %>%
#   mutate(same_species = if_else(species1==species2, 1, 0))

snpdist_SLs_filt_pt <- snpdist_SLs_filt_pt %>%
  mutate(same_niche = if_else(niche1==niche2, 1, 0)) %>%
  mutate(same_source = if_else(source1==source2, 1, 0)) %>%
  mutate(same_species = if_else(species1==species2, 1, 0))

#Plot distributions + p-value ----
#snpdist_3250_filt_pt %>%
snpdist_SLs_filt_pt %>%
  mutate(same_source = ifelse(same_source=="1","within-source","between-source")) %>%
  filter(distance <1000) %>%
  ggplot(aes(x = distance, color = as.factor(same_source))) + 
  geom_density(aes(fill = as.factor(same_source)), alpha = 0.5) + 
  labs(title = "Distribution of Distances",
    x = "Distance",
    y = "Density",
    color = "same_source",
    fill = "same_source") + 
  scale_colour_manual(values=c("between-source"="#624185", "within-source"="#ffa345" )) +
  scale_fill_manual(values=c("between-source"="#624185", "within-source"="#ffa345" )) +
  theme_minimal() +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)

ggsave("~/KGWP1/figures/pairwise_SNP_distribution.pdf", dpi = 300, width = 12, height = 8, units = "cm")


# Splitting the data based on the 'match' column, to run KS-test
within_niche <- snpdist_SLs_filt_pt %>% filter(same_niche == 1) %>% select(distance) %>% pull()
between_niche <- snpdist_SLs_filt_pt %>% filter(same_niche == 0) %>% select(distance) %>% pull()
# Performing the Kolmogorov-Smirnov test
ks_test <- ks.test(within_niche, between_niche) # p-value < 2.2e-16
ks_test
#D = 0.16181, p-value < 2.2e-16

#Plot by niche / source pairs ----
#Note not used in paper
snpdist_SLs_filt_pt %>% 
  mutate(pair = if_else(niche1 <= niche2, 
                        paste(niche1, niche2, sep = "-"), 
                        paste(niche2, niche1, sep = "-"))) %>%
  filter(distance <100) %>% #Show only up to 100 SNPs
  group_by(distance, pair) %>%  
  tally() %>%
  mutate(pair = ifelse(pair == "Animal-Human", "Human-Animal", pair)) %>%
  ggplot(aes(x = distance, y = n, fill = pair)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Count of pairs (niche1, niche2)",
       x = "Distance",
       fill = "Pair Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("Human-Human"="#E41A1C", 
                               "Human-Marine"="#377EB8", 
                               "Animal-Animal"="#4DAF4A", 
                               "Human-Animal"="#984EA3", 
                               "Marine-Marine"="#FF7F00", 
                               "Animal-Marine"="#FFFF33"))


#Plot Fig S1 ----
library(cutpointr) 

#SLs unfiltered - not used 
cp <- cutpointr(snpdist_SLs_filt_pt, distance, same_source, #or set same_niche
                method = maximize_metric, metric = sum_sens_spec, pos=1)

summary(cp)
plot(cp)

#SLs filtered 100 - used
snpdist_SLs_filt_pt_100 <- snpdist_SLs_filt_pt %>% filter(distance <= 100)
cp <- cutpointr(snpdist_SLs_filt_pt_100, distance, same_source,
                method = maximize_metric, metric = sum_sens_spec, pos=1)

summary(cp)
plot(cp) #This is Fig S1A - saved as PDF 

#Fig S1B: strain-sharing pairs
snpdist_SLs_filt_pt %>%
  filter(!is.na(distance)) %>%
  count(distance) %>%
  arrange(distance) %>%
  mutate(
    cumulative_pairs = cumsum(n),
    proportion_pairs = cumulative_pairs / sum(n)
  ) %>%
  filter(distance <= 100) %>%
  ggplot(aes(x = proportion_pairs, y = distance)) +
  geom_line(color = "black", size = 1.2) +
  geom_hline(yintercept = 22, linetype = "dotted", color = "red", linewidth = 1) +
  #scale_x_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(
    x = "Proportion of pairs",
    y = "Threshold",
    title = "Proportion of connected pairs by SNP threshold"
  )

