#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Fig S5 - Global MDR/HV clones
#####################################

#Read data ----
metadata_all <- read.csv("~/KGWP1/genotyping/metadata_all.csv", sep="\t", check.names = FALSE)
#Set order
metadata_all$source <- factor(metadata_all$Subtype, levels = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey", "Broiler", "Bivalves", "Seawater"))

#Define MDR and HV SLs
MDR_SLs <- c("SL101","SL147","SL15","SL258","SL29","SL307","SL37","SL17")
vir_SLs <- c("SL23","SL25","SL380","SL66","SL86")

#Plot as numbers
mdr_SLs_plot <- metadata_all %>%
  filter(SL %in% SLs) %>%
  mutate(MDR = ifelse(num_resistance_classes < 3, 0, 1)) %>%
  group_by(SL, source, MDR) %>%
  count() %>%
  group_by(SL, source) %>%
  mutate(total_n = sum(n), proportion = n / total_n) %>%
  ggplot(aes(x = source, y = proportion, fill = as.factor(MDR))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values=c("0"="gray", "1"="red")) +
  facet_wrap(~SL) +
  labs(x = "Source", y = "Proportion", fill = "MDR") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5) 

ggsave(filename= "~/KGWP1/figures/mdr_SLs_plot.pdf", plot = mdr_SLs_plot, dpi = 300, width = 12, height = 8, units = "in")

#Virulent clones
#Plot as numbers
vir_SLs_plot <- metadata_all %>% 
  filter(SL %in% vir_SLs) %>% 
  mutate(Aerobactin = ifelse(Aerobactin == "-",0,1) )%>% 
  mutate(Salmochelin = ifelse(Salmochelin == "-",0,1) )%>% 
  mutate(RmpADC = ifelse(RmpADC == "-",0,1)) %>% 
  group_by(SL,source,Aerobactin,RmpADC) %>% count() %>% 
  mutate(hyp = ifelse((Aerobactin == 1 & RmpADC == 1),"Iuc+RmpADC",
                      (ifelse(Aerobactin == 1 & RmpADC == 0, "iuc_only",#
                        (ifelse(Aerobactin == 0 & RmpADC == 1, "rmpadc_only","none")))))) %>%
  ungroup() %>%
  select(SL,source,hyp,n) %>%
  group_by(SL, source, hyp) %>%
  summarise(total_n = sum(n), .groups = 'drop') %>% 
  ggplot( aes(x = source, y = total_n, fill = as.factor(hyp))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~SL, scales = "free_y") +
  labs(x = "SL", y = "Count", fill = "Hypervirulent") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = total_n), vjust = -0.5, color = "black", size = 3.5) 

vir_SLs_plot

ggsave(filename= "~/KGWP1/figures/vir_SLs_plot.pdf", plot = vir_SLs_plot, dpi = 300, width = 12, height = 8, units = "in")
