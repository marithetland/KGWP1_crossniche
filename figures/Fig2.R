#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Fig 2: Plot and chi2 stats
#####################################

#Read data ----
metadata_all <- read.csv("~/KGWP1/genotyping/metadata_all.csv", sep=",", check.names = FALSE)

#Load pkgs ----
library(tidyverse)
library(purrr)

#Plot barplots ----
plot_operon_proportion_by_niche <- function(metadata_df, operon_col, output_dir = "~/KGWP1/figures/") {
  
  plot_hmrg <- metadata_df %>%
    rename(operon = !!sym(operon_col)) %>%
    group_by(niche, operon) %>%
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
    mutate(niche = factor(niche, levels = c("Human", "Animal", "Marine"))) %>%
    ggplot(aes(x = niche, y = Proportion, fill = niche)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(Human = "#CD2626", Animal = "#228B22", Marine = "#1874CD")) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    geom_text(aes(label = n), vjust = -0.5, color = "black", size = 3.5)

  ggsave(
    filename = paste0(output_dir, operon_col, "_by_niche.pdf"),
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
  plot_operon_proportion_by_niche(metadata_all, op)
}


#Stats ----
#Function to run chi2 overall and pairwise
run_chi2_tests <- function(data, binary_col) {
  
  col_sym <- rlang::ensym(binary_col)
  niche_pairs <- combn(unique(data$Niche), 2, simplify = FALSE)
  
  compute_chi2 <- function(df) {
    df %>%
      group_by(Niche, !!col_sym) %>%
      count() %>%
      pivot_wider(names_from = !!col_sym, values_from = n, values_fill = 0) %>%
      column_to_rownames("Niche") %>%
      as.matrix() %>%
      chisq.test()
  }
  
  interpret_p <- function(p) {
    case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  }
  
  overall_test <- data %>%
    group_by(Niche, !!col_sym) %>%
    count() %>%
    pivot_wider(names_from = !!col_sym, values_from = n, values_fill = 0) %>%
    column_to_rownames("Niche") %>%
    as.matrix() %>%
    chisq.test()
  
  pairwise_tests <- map(niche_pairs, function(pair) {
    subset_data <- data %>% filter(Niche %in% pair)
    test <- compute_chi2(subset_data)
    tibble(
      comparison = paste(pair, collapse = " vs "),
      p_value = test$p.value,
      statistic = unname(test$statistic),
      df = test$parameter,
      significance = interpret_p(test$p.value)
    )
  }) %>% bind_rows()
  
  overall_df <- tibble(
    comparison = "Overall",
    p_value = overall_test$p.value,
    statistic = unname(overall_test$statistic),
    df = overall_test$parameter,
    significance = interpret_p(overall_test$p.value)
  )
  
  bind_rows(overall_df, pairwise_tests)
}

#Run stats
hmrg_chi2_results <- map_df(operon_cols, function(op_col) {
  result <- run_chi2_tests(metadata_all, !!sym(op_col))
  result$operon <- op_col
  result
})

hmrg_chi2_results

