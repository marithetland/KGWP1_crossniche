#########################
# KGWP1 cross-niche paper
# 2023-03-30, Hetland MAK
# Figure S14: SL, KL, OL
#########################

library(tidyverse)
library(dplyr)
library(ggplot2)

#set env ----
set.seed(123)
setwd("~/KGWP1/genotyping/")

#Load and prep input files ----
metadata_all <- read.csv("metadata_all_2024-05-29.csv", sep="\t", check.names = FALSE)

#Correct K and O locus calls (kaptive 2) ----
metadata_all <- metadata_all %>% mutate(K_locus = if_else(K_locus=='unknown (KL107)', "Unknown", K_locus),
                                        K_locus = if_else(str_detect(K_locus, 'unknown'), str_extract(K_locus, "(?<=\\().*(?=\\))"), K_locus),
                                        K_type = if_else(str_detect(K_type, 'unknown'), str_extract(K_type, "(?<=\\().*(?=\\))"), K_type),
                                        O_locus = if_else(str_detect(O_locus, 'unknown'), str_extract(O_locus, "(?<=\\().*(?=\\))"), O_locus),
                                        O_type = if_else(str_detect(O_type, 'unknown'), str_extract(O_type, "(?<=\\().*(?=\\))"), O_type)) 
#Set colours ----
niche_colours <- c("Animal" = "#228B22", "Human" = "#CD2626", "Marine" = "#1874CD")
OL_colours <- c("O1/O2v1"="#FB0045", "O1/O2v2"="#00B838", "O12"="#FEE103", "O3/O3a"="#9E17BB", "O3b"="#FF7B00", "O4"="#FF10ED", "O5"="#02F5F0", "OL102"="#ACFA00", "OL103"="#008280", "OL13"="#030F7B", "OL104"="#EEBBFE", "O1/O2v3"="#3A61E0")

#Set plot orders
order_niche <- c("Human","Animal","Marine")
order_source <- c("Human_infection","Community_carriage","Dog","Pig","Turkey","Broiler","Bivalves","Seawater")

#Plot K loci by source ----
#Select the 10 most prevalent KLs in human and non-human samples:
top_k_locus <- c("KL21","KL14","KL22","KL30","KL28","KL38",
                 "KL31","KL23","KL24","KL55","KL57",
                 "KL102","KL28","KL103","KL2","KL10",
                 "KL25","KL38","KL22","KL62","KL24")

#Set colours
num_colours <- length(top_k_locus)
bright_colours <- c("#FF0000", "#0000FF", "#00FF00", "#FFFF00", "#FFA500", "#800080", 
                   "#FFC0CB", "#00FFFF", "#FF00FF", "#FFD700", "#40E0D0", 
                   "#EE82EE", "#4B0082", "#FF7F50", "#DC143C", "#FF69B4", "#008080", 
                   "#FA8072", "#1E90FF", "#7FFF00","red")
random_colours <- sample(bright_colours, num_colours)
colour_mapping <- setNames(random_colours, top_k_locus)
colour_mapping <- c(colour_mapping, "Other" = "gray")
colour_mapping


KL_perc <- metadata_all %>%
  group_by(source, K_locus) %>% 
  count() %>%
  ungroup() %>%
  mutate(K_locus = ifelse(K_locus %in% top_k_locus, K_locus,"Other")) %>%
  select(source, K_locus, n) %>%
  group_by(source, K_locus) %>%
  summarise(n = sum(n), .groups = 'drop') %>%
  pivot_wider(names_from=K_locus,values_from=n,values_fill=0) %>%
  mutate(total = rowSums(select(., -source))) %>%
  pivot_longer(cols=-c(source,total), names_to="K_locus", values_to="n") %>%
  mutate(percent=(n/total)*100)  %>%
  mutate(source = factor(source, levels = order_source)) %>%
  ggplot( aes(x = source, y = percent, fill = K_locus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Source", y = "Percent", fill = "K_locus") +
  scale_fill_manual(values = colour_mapping) +
  theme_bw() +
  coord_flip()

KL_perc

ggsave("~/KGWP1/figures/FigS14A_KLs_by_source.pdf", plot=KL_perc, dpi = 300, width = 18, height = 18, units = "cm")

#Plot O loci source ----
OL_perc <- metadata_all %>%
  group_by(source, O_locus) %>% 
  count() %>% 
  ungroup() %>%
  pivot_wider(names_from=O_locus,values_from=n,values_fill=0) %>%
  mutate(total = rowSums(select(., -source))) %>%
  pivot_longer(cols=-c(source,total), names_to="O_locus", values_to="n") %>%
  mutate(percent=(n/total)*100) %>%
  mutate(source = factor(source, levels = order_source)) %>%
  ggplot( aes(x = source, y = percent, fill = O_locus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Source", y = "Percent", fill = "O_locus") +
  scale_fill_manual(values=OL_colours)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  coord_flip()

OL_perc

ggsave("~/KGWP1/figures/FigS14B_OLs_by_source.pdf", plot=OL_perc, dpi = 300, width = 18, height = 18, units = "cm")


#Plot SLs vs KL+OL combinations, by niche ----
KL_OL_SL <- metadata_all %>% select(strain,SL,niche,source,K_locus,O_locus) %>%
  mutate(KO_loci = paste(K_locus, O_locus, sep = " + ")) %>%
  mutate(SLKO = paste(SL, KO_loci, sep = " + ")) %>%
  group_by(SL, KO_loci, niche) %>%
  tally() %>% 
  pivot_wider(names_from=niche, values_from=n, values_fill = list(n = 0)) %>% 
  mutate(total = Human + Marine + Animal) %>%
  filter(total > 10) %>%
  mutate(KO_loci = factor(KO_loci, levels = ko_counts$KO_loci),
         SL = factor(SL, levels = sl_counts$SL)) %>%
  pivot_longer(cols = Human:Animal, names_to = "niche", values_to = "n") %>%
  mutate(niche = factor(niche, levels = c("Human", "Animal", "Marine"))) %>%
  ggplot(aes(x = KO_loci, y = SL, size = n, color = niche)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ niche) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_manual(values = custom_colors) +  
  labs(x = "KO_loci", y = "SL", size = "Niche Count", color = "Niche") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),  
    axis.text.y = element_text(size = 8))

KL_OL_SL

ggsave("~/KGWP1/figures/FigS14C_KL_OL_SL_by_niche.pdf", plot=KL_OL_SL, dpi = 300, width = 24, height = 12, units = "cm")
