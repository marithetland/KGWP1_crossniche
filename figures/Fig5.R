#########################
# KGWP1 cross-niche paper
# 2023-03-30, Hetland MAK
# Figure 5: SNPs
#########################

library(tidyverse)
library(eulerr)
library(stats)

#set env ----
set.seed(123)
setwd("~/KGWP1/")

niche_colours <- c("Animal" = "#228B22", "Human" = "#CD2626", "Marine" = "#1874CD")

#Load and prep input files ----
metadata_all <- read.csv("./genotyping/metadata_all_2024-05-29.csv", sep="\t", check.names = FALSE)
snpdists <- read.csv("./transmission_analyses/SLs_reddog_snpdists.csv")

#Wrangle data for Figs 5B and 5C ----
snpdists_df <- snpdists %>%
  #Filter out rows where strain1 and strain2 are the same
  filter(strain1 != strain2) %>%
  #Sort strain1 and strain2 alphabetically in each row
  rowwise() %>%
  mutate(strain1 = min(strain1, strain2),
         strain2 = max(strain1, strain2)) %>%
  # Remove duplicated rows
  select(strain1, strain2, distance) %>%
  distinct(strain1, strain2, .keep_all = TRUE) %>%
  #Match snpdists and metadata
  left_join(metadata_all, by = c("strain1" = "strain")) %>%
  rename(niche1 = niche,source1 = source) %>%
  left_join(metadata_all, by = c("strain2" = "strain")) %>%
  rename(niche2 = niche, source2 = source) %>%
  #Sort naming
  mutate(source2 = if_else(source2 == "Blood" | source2 == "Urine" , "Human_infection", source2)) %>% 
  mutate(source1 = if_else(source1 == "Blood" | source1 == "Urine" , "Human_infection", source1)) %>%
  #Add 1/0 for if pairs are in same source/niche or not
  mutate(match = if_else(niche1==niche2, 1, 0),
         same_niche = if_else(niche1==niche2, 1, 0),
         same_source = if_else(source1==source2, 1, 0))

#Fig 5A ----
niche_counts <- metadata_all %>% 
  group_by(SL, niche) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = niche,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    present_Animal = Animal > 0,
    present_Human = Human > 0,
    present_Marine = Marine > 0,
    niche_combo = case_when(
      present_Animal & !present_Human & !present_Marine ~ "Animal only",
      !present_Animal & present_Human & !present_Marine ~ "Human only",
      !present_Animal & !present_Human & present_Marine ~ "Marine only",
      present_Animal & present_Human & !present_Marine ~ "Animal+Human",
      present_Animal & !present_Human & present_Marine ~ "Animal+Marine",
      !present_Animal & present_Human & present_Marine ~ "Human+Marine",
      present_Animal & present_Human & present_Marine ~ "All three",
      TRUE ~ "None"
    )
  ) %>% 
  count(niche_combo, name = "Num_SLs") %>%
  arrange(desc(Num_SLs))

fit <- euler(c(
  "Human"         = 673,   # Only Human
  "Animal&Human"  = 72,    # Animal + Human
  "Animal"        = 60,    # Only Animal
  "Human&Marine"  = 19,    # Human + Marine
  "Marine"        = 17,    # Only Marine
  "Animal&Human&Marine" = 12,  # All three
  "Animal&Marine" = 4     # Animal + Marine
))

plot(fit, fills = list(fill = c("#CD2626","#228B22","#1874CD"), alpha = 0.7),
     edges = TRUE, quantities = TRUE, labels = TRUE)

#Fig 5B ----
#See script Figs_5B_and_S1.R, but also:

#Split the data based on the 'match' column
within_niche <- snpdists_df %>% filter(match == 1) %>% select(distance) %>% pull()
between_niche <- snpdists_df %>% filter(match == 0) %>% select(distance) %>% pull()

#Kolmogorov-Smirnov test
ks_test <- ks.test(within_niche, between_niche) # p-value < 2.2e-16

#Plot
fig5B <- snpdists_df %>%
  mutate(match = ifelse(match=="1","within-niche","between-niche")) %>%
  filter(distance <1000) %>%
  ggplot(aes(x = distance, color = as.factor(match))) + 
  geom_density(aes(fill = as.factor(match)), alpha = 0.5) + 
  labs(
    title = "Distribution of Distances",
    x = "Distance",
    y = "Density",
    color = "Match",
    fill = "Match"
  ) + 
  scale_colour_manual(values=c("within-niche"="#ffa345", "between-niche"="#624185" )) +
  scale_fill_manual(values=c("within-niche"="#ffa345", "between-niche"="#624185" )) +
  theme_minimal() 

ggsave("~/KGWP1/figures/Fig5B.pdf", plot=fig5B, dpi = 300, width = 18, height = 18, units = "cm")


#Fig 5C ----
#Define mapping: source name in data -> column name you want in output
source_map <- c(
  Faeces = "human_carriage",
  Pig = "pig",
  Turkey = "turkey",
  Broiler = "broiler",
  Dog = "dog",
  Bivalves = "bivalve",
  Seawater = "seawater"
)

#Initialise output list
human_connect <- vector("list", length = 1001)

#Function to count human_infection connections to each source
count_connections <- function(data, source_type, max_dist) {
  s1 <- data %>%
    filter(as.numeric(distance) <= max_dist,
           source1 == "Human_infection", source2 == source_type) %>%
    pull(strain1)
  s2 <- data %>%
    filter(as.numeric(distance) <= max_dist,
           source2 == "Human_infection", source1 == source_type) %>%
    pull(strain2)
  length(unique(c(s1, s2)))
}

#Loop to count connections at SNP distance i
for (i in 0:1000) {
  # Connected human infection strains (to anything)
  h1 <- snpdists_df %>%
    filter(as.numeric(distance) <= i,
           strain1 != strain2, source1 == "Human_infection") %>%
    pull(strain1)
  h2 <- snpdists_df %>%
    filter(as.numeric(distance) <= i,
           strain1 != strain2, source2 == "Human_infection") %>%
    pull(strain2)
  hn <- length(unique(c(h1, h2)))
  
  #Get counts and rename
  counts <- map_dbl(names(source_map), ~count_connections(snpdists_df, .x, i))
  named_counts <- setNames(as.list(counts), source_map)
  
  #Store results for this distance
  human_connect[[i + 1]] <- tibble(x = i, n = hn, !!!named_counts)
}

# Combine results
human_connect <- bind_rows(human_connect)


#Plot - denominator = human samples
#x = threshold SNP distance
#y = of all human samples, what proportion have a connection to a non-human sample below this threshold distance

#Count number of human infection isolates, to use as denominator
n_human <- metadata_all %>%
  mutate(source = if_else(source == "Blood" | source == "Urine" , "Human_infection", source)) %>%
  filter(source == "Human_infection") %>%
  nrow()

sources_to_plot <- c("human_carriage", "pig", "turkey", "broiler", "dog", "bivalve", "seawater")

fig5C <- human_connect %>%
  mutate(across(all_of(sources_to_plot), ~ .x / n_human, .names = "{.col}_prop_dist")) %>%
  select(x, ends_with("_prop_dist")) %>%
  pivot_longer(cols = -x, names_to = "source", values_to = "probability") %>%
  #reshape2::melt(id.vars=c("x"), variable.name="source", value.name = "probability") %>%
  ggplot(aes(x = x, y = probability, color = source)) +
  geom_line(size = 1, linetype = "solid") +  
  #ggplot(aes(x = x, y = probability, col = source, linetype = source)) + 
  theme_bw() + 
  labs(x="genetic distance", y="Proportion of human infection samples\n that are connected to a non-human sample") +
  scale_color_manual(values = c(
    "human_carriage__prop_dist" = "#E63946",
    "pig_prop_dist" = "#F4A261",
    "turkey_prop_dist" = "#2A9D8F",
    "broiler_prop_dist" = "#264653",
    "dog_prop_dist" = "#2B9348",
    "bivalve_prop_dist" = "#FCCD04",
    "seawater_prop_dist" = "#9C89B8"))

fig5C 
ggsave("~/KGWP1/figures/Fig5C.pdf", plot=fig5C, dpi = 300, width = 12, height = 8, units = "cm")


#Fig 5D ----
#Load data generated in FigS15.R
count_source_events <- read.csv("~/KGWP1/transmission_analyses/run_transmission_comparison/source_transmission_events_by_SL.csv")

fig5D <- count_source_events %>%
  group_by(source_pair) %>%
  summarise(Total_Transmission_Count = sum(Count), .groups = "drop") %>%
  arrange(desc(Total_Transmission_Count)) %>%
  separate(source_pair, into = c("xaxis", "yaxis"), sep = "-") %>%
  ggplot(aes(x = xaxis, y = yaxis)) +
  geom_point(aes(size = Total_Transmission_Count), color = "black") +
  scale_x_discrete(
    name = "", 
    limits = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler", "Bivalves", "Seawater"),
  ) +
  scale_y_discrete(
    name = "", 
    limits = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler", "Bivalves", "Seawater"),
  ) +
  scale_size_continuous(range = c(2, 15), name = "Transmission count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Total_Transmission_Count), hjust = -2, size = 3) +
  coord_flip()

ggsave("~/KGWP1/figures/Fig5D.pdf", plot=fig5D, dpi = 300, width = 18, height = 18, units = "cm")
