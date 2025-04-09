#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Plot AMR determinants by source
#####################################

metadata_all <- read.csv("~/KGWP1/actual_data_for_paper/genotyping/metadata_all.csv", sep="\t", check.names = FALSE)

#Load pkgs ----
library(dplyr)
library(ggplot2)
library(tidyverse)

# Plot all AMR genes by source as bubble plot
overlapping_AMR_genes <- metadata_all %>% 
  select(strain,Subtype,AGly_acquired:Bla_Carb_acquired,Omp_mutations:Flq_mutations) %>% 
  pivot_longer(cols = -c("strain","Subtype"), names_to = "AMR_acquired", values_to = "AMR_gene") %>%
  mutate(AMR_gene = sapply(str_split(AMR_gene, ";"), function(x) paste(unique(x), collapse = ";"))) %>% #Deduplicate if same gene found >1 in same isolates
  mutate(across('AMR_gene', str_replace, "\\^", ""))  %>%
  mutate(across('AMR_gene', str_replace, "\\*", ""))  %>%
  mutate(across('AMR_gene', str_replace, "\\?", "")) %>%
  mutate(across('AMR_gene', str_replace, '-', "0"))  %>% #change to only replace if whole string
  filter(AMR_gene >0) %>%
  select(-c("strain")) %>%
  separate_rows(AMR_gene, sep=";") %>%#
  mutate(across('AMR_gene', str_replace, "\\^", ""))  %>%
  mutate(across('AMR_gene', str_replace, "\\*", ""))  %>%
  mutate(across('AMR_gene', str_replace, "\\?", "")) %>%
  mutate(across('AMR_gene', str_replace, '-', "0"))  %>%
  group_by(Subtype,AMR_acquired,AMR_gene) %>% #
  summarize(num_present = n(), .groups = "drop") %>%#
  # group_by(Subtype,AMR_gene) %>% 
  # count() %>%
  #mutate(across('AMR_gene', str_replace, "\\^", ""))  %>%
  #mutate(across('AMR_gene', str_replace, "\\*", ""))  %>%
  #mutate(across('AMR_gene', str_replace, "\\?", "")) %>%
  #TODO: Must combine the rows to get one row per num_present
  mutate(class_AMR = paste(AMR_acquired,AMR_gene, sep="-")) %>% #
  arrange(AMR_acquired,AMR_gene) %>%#
  mutate(class_AMR = factor(class_AMR, levels= unique(class_AMR))) %>% #
  #select(Subtype,class_AMR,num_present) %>%
  pivot_wider(names_from=Subtype, values_from=num_present, values_fill = 0) %>%
  mutate(total_nonhuman=(Bivalves+Broiler+Dog+Pig+Turkey+Broiler+Seawater)) %>%
  filter(total_nonhuman >0) %>%
  select(-c("total_nonhuman")) %>%
  pivot_longer(cols = -c("AMR_acquired","AMR_gene","class_AMR"), names_to = "Subtype", values_to = "num_present") %>%
  filter(num_present != 0)  %>%
  filter(AMR_gene != "MgrB00%")  #remove incorrect pig MGRB mutation

overlapping_AMR_genes$Subtype <- factor(overlapping_AMR_genes$Subtype, levels = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler",  "Bivalves", "Seawater"))


main_plot <- overlapping_AMR_genes %>%
  ggplot(aes(x=Subtype, y=class_AMR, fill=num_present)) + 
  geom_tile(color="white") +
  geom_text(aes(label= num_present), size=3) +
  scale_fill_gradientn(
    colors = c("white", "lightblue", "blue", "red"),
    values = scales::rescale(c(0, 1, 100, max(overlapping_AMR_genes$num_present))),
    name = "Number Present"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

main_plot
ggsave(filename = "~/KGWP1/figures/AMR_pres.pdf", plot=main_plot ,width = 12, height = 8, dpi = 300)