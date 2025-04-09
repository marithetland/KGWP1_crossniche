#########################
# KGWP1 cross-niche paper
# 2023-03-30, Winkler MA
# Figure S8: HMRGs by replicon
#########################

hmrg_by_replicon_df <- read.csv("~/KGWP1/genotyping/hybrid_splitrep_hmrg_long__2024-07-17.csv")


location_labels = c("chr" = "Chromosome", "pl" = "Plasmid", "ctg" = "Contig", "both" = "Chromosome & Plasmid")
location_order <- c("pl", "both", "chr")

export_gene_labels = c("arsA", "arsB", "arsC", "arsD", "arsH", "arsR", "silA", "silB", "silC", "silE", "silF", "silP", "silR", "silS",
                       "chrA", "chrB1", "merA", "merB", "merC", "merD", "merE", "merF", "merP", "merR", "merT",
                       "copA", "copB", "copC", "copD", "copG", "copZ", "pcoA", "pcoB", "pcoC", "pcoD", "pcoE", "pcoR", "pcoS",
                       "terA", "terB", "terC", "terD", "terE", "terF", "terL", "terW", "terX", "terY", "terZ",
                       "rcnA", "rcnB", "rcnR", "ncrA", "ncrB", "ncrC", "ncrY", "clpK", "hsp20")

export_facet_labels<- c("antimony_arsenic" = "Antimony, Arsenic\n(arsABCDR)", "silver" = "Silver\n(silABCERS)", "chromium" = "Chromium\n(chrA, chrB1)", "mercury" = "Mercury\n(merACP, merAFP, merAPT)",
                        "copper_cop" = "Copper\n(copABCD)", "copper_pco" = "Copper\n(pcoABCDERS)", "tellurite" ="Tellurite\n(terBCDE)", "cobalt_nickel" = "Cobalt, Nickel\n(rcnAR)",
                        "nickel" = "Nickel\n(ncrABC)", "heat_shock" = "Heat Shock\n(clpK, hsp20)")

export_facet_order <- c("antimony_arsenic", "silver", "chromium", "mercury", "copper_cop", "copper_pco", "tellurite", "cobalt_nickel", "nickel", "heat_shock")


hmrg_by_replicon_df$location <- factor(hmrg_by_replicon_df$location, levels = location_order)
hmrg_by_replicon_df$metal_res <- factor(hmrg_by_replicon_df$metal_res, levels = export_facet_order)

ggplot(hmrg_by_replicon_df, aes(x = Gene, fill = location)) +
  geom_bar(position = "fill", aes(y = ..count.. / sum(..count..)), show.legend = FALSE) +
  scale_fill_manual(values = c("chr" = "blue", "both" = "green4", "pl" = "red"),
                    labels = location_labels) +
  labs(title = "Gene Presence on Chromosomes and Plasmids by Metal Resistance",
       x = "Gene",
       y = "Proportion of Count",
       fill = "Gene\nLocation") +
  facet_wrap(~ metal_res, scales = "free_x", ncol = 4, labeller = labeller(metal_res = as_labeller(export_facet_labels))) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.x = element_text(size = 10)) 

ggsave(filename="~/KGWP1/figures/FigS8.pdf", width = 10, height = 8, dpi = 300)
