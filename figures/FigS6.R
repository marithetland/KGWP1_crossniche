#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Plot virulence loci by source
#####################################

metadata_all <- read.csv("~/KGWP1/actual_data_for_paper/genotyping/metadata_all.csv", sep="\t", check.names = FALSE)

#Load pkgs ----
library(dplyr)
library(ggplot2)
library(tidyverse)

# Plot all VIR genes by source as bubble plot
overlapping_VIR_genes <- metadata_all %>% 
  select(strain,Subtype,Yersiniabactin,Colibactin,Aerobactin,Salmochelin,RmpADC) %>% 
  pivot_longer(cols = -c("strain","Subtype"), names_to = "VIR_acquired", values_to = "VIR_gene") %>%
  mutate(VIR_gene = sapply(str_split(VIR_gene, ";"), function(x) paste(unique(x), collapse = ";"))) %>% #Deduplicate if same gene found >1 in same isolates
  mutate(across('VIR_gene', str_replace, "\\^", ""))  %>%
  mutate(across('VIR_gene', str_replace, "\\*", ""))  %>%
  mutate(across('VIR_gene', str_replace, "\\?", "")) %>%
  mutate(across('VIR_gene', str_replace, '-', "0"))  %>% #change to only replace if whole string
  filter(VIR_gene >0) %>%
  select(-c("strain")) %>%
  #separate_rows(VIR_gene, sep=";") %>%#
  mutate(across('VIR_gene', str_replace, "\\^", ""))  %>%
  mutate(across('VIR_gene', str_replace, "\\*", ""))  %>%
  mutate(across('VIR_gene', str_replace, "\\?", "")) %>%
  mutate(across('VIR_gene', str_replace, '-', "0"))  %>%
  group_by(Subtype,VIR_acquired,VIR_gene) %>% #
  summarize(num_present = n(), .groups = "drop")  %>%
  #Uncomment below to get raw numbers
  #pivot_wider(names_from="Subtype",values_from=num_present,values_fill=0) # %>% filter(VIR_acquired == "Yersiniabactin") %>% arrange(VIR_acquired) %>% print(n=200)  %>% mutate(tot=Bivalves+Broiler +Community_carriage  + Dog +   Pig +Seawater +Turkey) %>% filter(tot >0) %>% print(n=200)
  mutate(class_AMR = paste(VIR_acquired,VIR_gene, sep="-")) %>% #
  arrange(VIR_acquired,VIR_gene) %>%#
  mutate(class_AMR = factor(class_AMR, levels= unique(class_AMR))) %>% #
  #select(Subtype,class_AMR,num_present) %>%
  pivot_wider(names_from=Subtype, values_from=num_present, values_fill = 0) %>%
  mutate(total_nonhuman=(Bivalves+Broiler+Dog+Pig+Turkey+Broiler+Seawater)) %>%
  filter(!(VIR_acquired == "Yersiniabactin" & total_nonhuman > 0)) %>%
  select(-c("total_nonhuman")) %>%
  pivot_longer(cols = -c("VIR_acquired","VIR_gene","class_AMR"), names_to = "Subtype", values_to = "num_present") %>%
  filter(num_present != 0) 

overlapping_VIR_genes

overlapping_VIR_genes$Subtype <- factor(overlapping_VIR_genes$Subtype, levels = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler",  "Bivalves", "Seawater"))

overlapping_VIR_genes %>%
  ggplot(aes(x=Subtype, y=class_AMR, fill=num_present)) + 
  geom_tile(color="white") +
  geom_text(aes(label= num_present), size=3) +
  scale_fill_gradient(low="white",  high="red" ) +
  theme_minimal()

vir_main_plot <- overlapping_VIR_genes %>%
  ggplot(aes(x=Subtype, y=class_AMR, fill=num_present)) + 
  geom_tile(color="white") +
  geom_text(aes(label= num_present), size=3) +
  scale_fill_gradientn(
    colors = c("white", "lightblue", "blue", "red"),
    values = scales::rescale(c(0, 1, 100, max(overlapping_VIR_genes$num_present))),
    name = "Number Present"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

vir_main_plot
ggsave(filename = "~/KGWP1/figures/VIR_pres.pdf", plot=main_plot ,width = 12, height = 8, dpi = 300)