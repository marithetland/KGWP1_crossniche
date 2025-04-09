#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Fig S4
#####################################

#Plots for S4. See file FigS4_statistics.R for chi2/Kruskal Wallis/Mann Whitney

#Read data ----
metadata_all <- read.csv("~/KGWP1/genotyping/metadata_all.csv", sep="\t", check.names = FALSE)

#Load pkgs ----
library(tidyverse)
library(purrr)
library(tibble)

#Plot barplots ----

#AMR
metadata_all %>% 
  mutate(AMR = ifelse(num_resistance_classes >= 1, 1, 0)) %>%
  group_by(niche,AMR) %>% 
  count() %>% 
  pivot_wider(names_from=AMR, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="AMR") %>% 
  filter(AMR == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_AMR_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")

#MDR
metadata_all %>% 
  mutate(MDR = ifelse(num_resistance_classes >= 3, 1, 0)) %>%
  group_by(niche,MDR) %>% 
  count() %>% 
  pivot_wider(names_from=MDR, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="MDR") %>% 
  filter(MDR == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_MDR_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")


#Ybt
metadata_all %>% 
  mutate(Yersiniabactin = ifelse(Yersiniabactin != "-",1,0) )%>%
  group_by(niche,Yersiniabactin) %>% 
  count() %>% 
  pivot_wider(names_from=Yersiniabactin, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="Yersiniabactin") %>% 
  filter(Yersiniabactin == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_Ybt_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")


#Iuc
metadata_all %>% 
  mutate(Aerobactin = ifelse(Aerobactin != "-",1,0) )%>%
  group_by(niche,Aerobactin) %>% 
  count() %>% 
  pivot_wider(names_from=Aerobactin, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="Aerobactin") %>% 
  filter(Aerobactin == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_iuc_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")

#rmpadc
metadata_all %>% 
  mutate(RmpADC = ifelse(RmpADC != "-",1,0) )%>%
  group_by(niche,RmpADC) %>% 
  count() %>% 
  pivot_wider(names_from=RmpADC, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="RmpADC") %>% 
  filter(RmpADC == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_rmp_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")

#iuc_rmpADC
metadata_all %>% 
  mutate(RmpADC=ifelse(RmpADC == "-",0,1)) %>% 
  mutate(iuc=ifelse(Aerobactin == "-",0,1)) %>% 
  mutate(iuc_rmpADC=RmpADC+iuc) %>%
  mutate(iuc_rmpADC=ifelse(iuc_rmpADC == 2, 1, 0)) %>%
  group_by(niche,iuc_rmpADC) %>% 
  count() %>% 
  pivot_wider(names_from=iuc_rmpADC, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="iuc_rmpADC") %>% 
  filter(iuc_rmpADC == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_iucANDrmp_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")


#thermoresistance
metadata_all %>% 
  mutate(hsp20=ifelse(hsp20__op == 0,0,1)) %>% 
  mutate(clpK=ifelse(clpK__op == 0,0,1)) %>% 
  mutate(thermoresistance=hsp20+clpK) %>%
  mutate(thermoresistance=ifelse(thermoresistance == 2, 1, 0)) %>%
  group_by(niche,thermoresistance) %>% 
  count() %>% 
  pivot_wider(names_from=thermoresistance, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="thermoresistance") %>% 
  filter(thermoresistance == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_thermoresistance_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")


#plasmid replicon markers
metadata_all %>% 
  mutate(plmarkers=ifelse(NUM_FOUND == 0,0,1)) %>% 
  group_by(niche,plmarkers) %>% 
  count() %>% 
  pivot_wider(names_from=plmarkers, values_from=n, values_fill = 0) %>% 
  mutate(total=`0`+`1`) %>% pivot_longer(cols=c(`0`,`1`), values_to="n", names_to="plmarkers") %>% 
  filter(plmarkers == 1) %>% 
  mutate(Proportion=n/total) %>%
  ggplot(aes(x = niche, y = Proportion, fill = niche)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#CD2626","#228B22","#1874CD")) + 
  theme_bw() +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave("~/KGWP1/figures/prop_plasmidRepMarkers_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")


#Plot boxplots ----
#plasmid replicon markers -range
metadata_all %>%
  mutate(Niche = factor(Niche, levels = c("Human", "Animal", "Marine"))) %>%
  ggplot(aes(x = Niche, y = NUM_FOUND, fill = Niche)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
  theme_bw() +
  scale_fill_manual(values = c(Human = "#CD2626", Animal = "#228B22", Marine = "#1874CD")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("~/KGWP1/figures/range_plasmidRepMarkers_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")



#heavy metal operons (excluding thermoresistance) -range
metadata_all %>%
  mutate(sum_hmrgs = rowSums(across(arsABCDR__op:telAB__op))) %>%
  mutate(Niche = factor(Niche, levels = c("Human", "Animal", "Marine"))) %>%
  ggplot(aes(x = Niche, y = sum_hmrgs, fill = Niche)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.8) +
  theme_bw() +
  scale_fill_manual(values = c(Human = "#CD2626", Animal = "#228B22", Marine = "#1874CD")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("~/KGWP1/figures/range_HMRGoperons_by_niche.pdf", dpi = 300, width = 12, height = 8, units = "cm")


