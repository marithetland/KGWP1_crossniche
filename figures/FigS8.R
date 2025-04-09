#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Fig S8
#####################################

#Read data ----
metadata_all <- read.csv("~/KGWP1/genotyping/metadata_all.csv", sep=",", check.names = FALSE)

#Load pkgs ----
library(tidyverse)
library(purrr)

#Set colours
source_colours = c("Human_infection" = "#CC004C", "Community_carriage" = "#a38b46", "Turkey" = "#6460AA", 
                   "Broiler" = "#0089D0", "Pig" = "#FCD7DE", "Dog" = "#0DB14B", "Bivalves" = "#00008B", "Seawater" = "#00004B") 

#Plot barplots ----
plot_operon_proportion_by_Subtype <- function(metadata_df, operon_col, output_dir = "~/KGWP1/figures/") {
  
  plot_hmrg <- metadata_df %>%
    rename(operon = !!sym(operon_col)) %>%
    group_by(Subtype, operon) %>%
    count() %>%
    pivot_wider(names_from = operon, values_from = n, values_fill = 0) %>%
    # Ensure both '0' and '1' columns exist
    mutate(
      `0` = if (!"0" %in% names(.)) 0 else `0`,
      `1` = if (!"1" %in% names(.)) 0 else `1`
    ) %>%
    mutate(total = `0` + `1`) %>%
    pivot_longer(cols = c(`0`, `1`), values_to = "n", names_to = "operon") %>%
    filter(operon == 1) %>%
    mutate(Proportion = n / total) %>%
    mutate(Subtype = factor(Subtype, levels = c("Human_infection", "Community_carriage", "Dog", "Pig", "Turkey","Broiler","Bivalves","Seawater"))) %>%
    ggplot(aes(x = Subtype, y = Proportion, fill = Subtype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = source_colours) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5)
  
  ggsave(
    filename = paste0(output_dir, operon_col, "_by_source.pdf"),
    plot = plot_hmrg,
    dpi = 300, width = 12, height = 8, units = "cm"
  )
  
}

operon_cols <- c(
  "arsABCDR__op", "cadABCR__op", "chrA__op", "copABCD__op", "cusABCF__op", "czcABC__op",
  "merACP_AFP_APT__op", "ncrABC__op", "nirABCD__op", "pcoABCDRS__op", "rcnAR__op",
  "silABCERS__op", "terBCDE__op", "klaABC__op", "telAB__op", "clpK__op", "hsp20__op"
)

for (op in operon_cols) {
  plot_operon_proportion_by_Subtype(metadata_all, op)
}
