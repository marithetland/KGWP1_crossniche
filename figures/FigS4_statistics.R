#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Run chi2 tests (Fig S4)
#####################################

#Read data ----
metadata_all <- read.csv("~/KGWP1/actual_data_for_paper/genotyping/metadata_all.csv", sep="\t", check.names = FALSE)

#Load pkg s----
library(tidyverse)
library(purrr)
library(tibble)

#Functions ---
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

run_range_tests <- function(data, numeric_col) {
  
  col_sym <- rlang::ensym(numeric_col)
  niche_levels <- unique(data$Niche)
  niche_pairs <- combn(niche_levels, 2, simplify = FALSE)
  
  interpret_p <- function(p) {
    case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  }
  
  # Overall Kruskal-Wallis test
  formula <- as.formula(paste(rlang::as_string(col_sym), "~ Niche"))
  overall_test <- kruskal.test(formula, data = data)
  
  # Pairwise Wilcoxon rank-sum tests
  pairwise_tests <- map(niche_pairs, function(pair) {
    subset_data <- data %>% filter(Niche %in% pair)
    test <- wilcox.test(reformulate("Niche", response = rlang::as_string(col_sym)), data = subset_data)
    tibble(
      comparison = paste(pair, collapse = " vs "),
      p_value = test$p.value,
      statistic = unname(test$statistic),
      significance = interpret_p(test$p.value)
    )
  }) %>% bind_rows()
  
  overall_df <- tibble(
    comparison = "Overall",
    p_value = overall_test$p.value,
    statistic = NA_real_,  # Kruskal-Wallis uses chi-sq stat; not comparable to Wilcoxon
    significance = interpret_p(overall_test$p.value)
  )
  
  bind_rows(overall_df, pairwise_tests)
}

#Run chi2 tests ----

#AMR
metadata_all %>% 
  mutate(AMR = ifelse(num_resistance_classes >= 1, 1, 0)) %>%
  run_chi2_tests(binary_col = AMR)

#MDR
metadata_all %>% 
  mutate(MDR = ifelse(num_resistance_classes >= 1, 1, 0)) %>%
  run_chi2_tests(binary_col = MDR)

#ESBL
metadata_all %>% 
  mutate(ESBL=ifelse(Bla_ESBL_acquired == "-",0,1)) %>% 
  run_chi2_tests(binary_col = ESBL)

#Virulence
#YBT
metadata_all %>% 
  mutate(Ybt=ifelse(Yersiniabactin == "-",0,1)) %>% 
  run_chi2_tests(binary_col = Ybt)

#iuc
metadata_all %>% 
  mutate(iuc=ifelse(Aerobactin == "-",0,1)) %>% 
  run_chi2_tests(binary_col = iuc)

#iuc
metadata_all %>% 
  mutate(iro=ifelse(Salmochelin == "-",0,1)) %>% 
  run_chi2_tests(binary_col = iro)

#rmp
metadata_all %>% 
  mutate(RmpADC=ifelse(RmpADC == "-",0,1)) %>% 
  run_chi2_tests(binary_col = RmpADC)

#clb
metadata_all %>% 
  mutate(clb=ifelse(Colibactin == "-",0,1)) %>% 
  run_chi2_tests(binary_col = clb)

#Aero+rmp
metadata_all %>% 
  mutate(RmpADC=ifelse(RmpADC == "-",0,1)) %>% 
  mutate(iuc=ifelse(Aerobactin == "-",0,1)) %>% 
  mutate(iuc_rmpADC=RmpADC+iuc) %>%
  mutate(iuc_rmpADC=ifelse(iuc_rmpADC == 2, 1, 0)) %>%
  run_chi2_tests(binary_col = iuc_rmpADC)

#plasmid replicon markers
metadata_all %>% 
  mutate(NUM_FOUND=ifelse(NUM_FOUND == 0,0,1)) %>% 
  run_chi2_tests(binary_col = NUM_FOUND)

#HMRGs
metadata_all %>% 
  mutate(sum_hmrgs = rowSums(across(arsABCDR__op:telAB__op))) %>%
  mutate(sum_hmrgs=ifelse(sum_hmrgs == 0,0,1)) %>% 
  run_chi2_tests(binary_col = sum_hmrgs)

#Thermoresistance
metadata_all %>% 
  mutate(hsp20__op=ifelse(hsp20__op == 0,0,1)) %>% 
  run_chi2_tests(binary_col = hsp20__op)

metadata_all %>% 
  mutate(clpK__op=ifelse(clpK__op == 0,0,1)) %>% 
  run_chi2_tests(binary_col = clpK__op)


#Run Kruskal-Wallis and Wilcoxon tests ----
#Num AMR classes per genome 
metadata_all %>%
  run_range_tests(numeric_col = num_resistance_classes)

#Num AMR genes per genome 
metadata_all %>%
  run_range_tests(numeric_col = num_resistance_genes)

#Num Plasmid replicon markers per genome 
metadata_all %>%
  run_range_tests(numeric_col = NUM_FOUND)

#Num HMRGs per genome 
metadata_all %>%
  mutate(sum_hmrgs = rowSums(across(arsABCDR__op:telAB__op))) %>%
  run_range_tests(numeric_col = sum_hmrgs)