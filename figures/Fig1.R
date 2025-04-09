#########################
# KGWP1 cross-niche paper
# 2023-03-30, Hetland MAK
# Figure: Dataset and 
# phylogroups
#########################

#Setwd ----
setwd("~/KGWP1/genotyping/")

#Load pkgs ----
library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(gplots)
library(dplyr)
library(data.table)

#Load  files ----
#Metadata, phylogroups and STs from Kleborate
metadata <- read.csv("../metadata/sourceinfo__2023-03-16.tsv", sep="\t", check.names = FALSE)
metadata_prevalence <- metadata
kleborate <- read.csv("Kleborate_v2.4.0_best3255__2024-01-09.txt", sep="\t", check.names = FALSE)
metadata_kleborate <- merge(metadata_prevalence,kleborate,by="strain") #Combine data
ont_sequenced <- read.csv("../metadata/list_ont_sequenced_genomes__2023-03-30.csv", sep=",", check.names = FALSE)
metadata_kleborate <- metadata_kleborate %>%
  left_join(ont_sequenced,by="strain") #Combine data with ONT data

#Panel A: Dataset with Illumina + ONT counts by source type
plot_panelA <- metadata_kleborate %>% 
  group_by(Subtype) %>%
  select(Subtype,ONT_sequenced, Niche) %>%
  mutate(ONT_sequenced = ifelse(is.na(ONT_sequenced), "No", "Yes")) %>%
  mutate(Subtype = ifelse(Subtype == "Urine" , "Blood", Subtype)) %>%
  mutate(ONT_sequenced=ifelse(Niche=="Human",gsub("No","No_human",ONT_sequenced),ONT_sequenced)) %>%
  mutate(ONT_sequenced=ifelse(Niche=="Human",gsub("Yes","AYes_human",ONT_sequenced),ONT_sequenced)) %>%
  mutate(ONT_sequenced=ifelse(Niche=="Animal",gsub("No","No_animal",ONT_sequenced),ONT_sequenced)) %>%
  mutate(ONT_sequenced=ifelse(Niche=="Animal",gsub("Yes","AYes_animal",ONT_sequenced),ONT_sequenced)) %>%
  mutate(ONT_sequenced=ifelse(Niche=="Marine",gsub("No","No_marine",ONT_sequenced),ONT_sequenced)) %>%
  mutate(ONT_sequenced=ifelse(Niche=="Marine",gsub("Yes","AYes_marine",ONT_sequenced),ONT_sequenced)) %>%
  group_by(Subtype,ONT_sequenced) %>%
  count() %>%
  mutate(freq = ifelse(Subtype == "Blood" & ONT_sequenced == "No_human", 700, n)) %>% #This reduces the size of the Illumina-only blood isolates. In AI, fix it so that there is a break.
  ggplot(aes(fill=ONT_sequenced, y=freq, x=Subtype)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("AYes_human" = "#8E0500", "No_human"= "#B30600", "AYes_animal" ="#29AF34", "No_animal"= "#00C914",  "AYes_marine"= "#2347BD", "No_marine"= "#3760E3")) +
  #scale_x_discrete(name ="", limits=c("Seawater", "Bivalve","Broiler","Turkey","Pig","Dog","Faeces","Urine","Blood") ) + #, labels=c("0_SNPs"="0 SNPs", "2_SNPs"="2 SNPs", "5_SNPs" = "5 SNPs", "10_SNPs" = "10 SNPs", "20_SNPs" = "20 SNPs"))  +
  scale_x_discrete(name ="", limits=c("Blood","Faeces","Dog","Pig","Turkey","Broiler","Bivalves","Seawater"), labels=c("Blood" = "Hospital infection","Faeces" = "Community carriage","Dog" = "Dog","Pig" = "Pig","Turkey" = "Turkey","Broiler" = "Broiler","Bivalves" = "Bivalves","Seawater" = "Seawater"))  +
  #xlab("SNP threshold") +
  ylab("Number of genomes") +
  coord_flip() +
  theme_bw()

plot_panelA
ggsave(filename="~/KGWP1/figures/Fig1A.pdf", width = 10, height = 8, dpi = 300)



#Plot B ----
#Panel b: Percentage dot plots of KpSC species by source type

plot_panelB <- metadata_kleborate %>% 
  group_by(Subtype,species,Niche) %>%
  select(Subtype,species, Niche) %>%
  mutate(Subtype = ifelse(Subtype == "Urine" , "Blood", Subtype)) %>%
  count() %>% 
  group_by(Subtype) %>%
  transmute(Subtype, species, percent = n/sum(n)*100) %>% #Percent of source
  ggplot(aes(x=Subtype, colour=species, y=species))+ 
  geom_point(aes(size=percent)) +
  scale_x_discrete(name ="", limits=c("Blood","Faeces","Dog","Pig","Turkey","Broiler","Bivalves","Seawater"), labels=c("Blood" = "Hospital infection","Faeces" = "Community carriage","Dog" = "Dog","Pig" = "Pig","Turkey" = "Turkey","Broiler" = "Broiler","Bivalves" = "Bivalves","Seawater" = "Seawater"))  +
  scale_y_discrete(name ="", limits=c("Klebsiella pneumoniae","Klebsiella variicola subsp. variicola","Klebsiella quasipneumoniae subsp. similipneumoniae","Klebsiella quasipneumoniae subsp. quasipneumoniae","Klebsiella quasivariicola","Klebsiella africana")) + #, labels=c("Blood" = "Hospital infection","Faeces" = "Community carriage","Dog" = "Dog","Pig" = "Pig","Turkey" = "Turkey","Broiler" = "Broiler","Bivalves" = "Bivalves","Seawater" = "Seawater"))  +
  scale_colour_manual(values=c("#FCB711","#6460AA","#F37021","#CC004C", "#0089D0","#0DB14B")) +
  coord_flip() +
  theme_bw() 
plot_panelB

ggsave(filename="~/KGWP1/figures/Fig1B.pdf", width = 16, height = 8, dpi = 300)


#Plot C ---- 
df <- read.csv("~/KGWP1/actual_data_for_paper/genotyping/Fig_clin_rel_features.tsv",sep="\t")
df <- df %>% filter(Feature != "Hypervirulence")%>% filter(Feature != "Virulence")

#Set orders
df$source <- factor(df$source, levels = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler", "Bivalves", "Seawater"))
df$Feature <- factor(df$Feature, levels = c("AMR", "MDR", "Carbapenemase", "Colistin", "ESBL", 
                                            "Yersiniabactin", "Aerobactin", "rmpADC", "iuc + rmpADC", "Heavy metal resistance (all)", 
                                            "Heavy metal resistance (plasmid)", "Thermoresistance", "Plasmid replicon markers"))
#Plot
df %>% 
  filter(Proportion_with_AMR > 0) %>%
  ggplot(aes(x = source, y = Feature, size = Proportion_with_AMR)) +
  geom_point(aes(color = Feature)) +
  geom_text(aes(label = paste0(round(Strains_with_AMR, 1))), 
            color = "black", 
            size = 3, 
            vjust = 0.5, 
            hjust = -0.5) +  
  scale_size_continuous(
    range = c(3, 20), 
    breaks = c(1, 2, 5, 10, 15, 25, 50, 75, 100),
    labels = c("1%", "2%", "5%", "10%", "15%", "25%", "50%", "75%", "100%")  ) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
  

ggsave(filename="~/KGWP1/figures/Fig1C.pdf", width = 16, height = 8, dpi = 300)
