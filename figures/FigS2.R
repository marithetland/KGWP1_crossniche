#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Fig S2 - qqplots
#####################################

#Setwd
set.seed(123)
pwd="~/KGWP1/GWAS/pyseer_june_2024/"

#Load pyseer output files (Kp1 only) ----
#COGs - genes from panaroo
pyseer_Kp1_COG_human <- read.csv(paste0(pwd,"KGWP1_Kp1_human_COG_output.tsv"), sep="\t")
pyseer_Kp1_COG_animal <- read.csv(paste0(pwd,"KGWP1_Kp1_animal_COG_output.tsv"), sep="\t")

#sCOGs - structural genes from panaroo, 3 consecutive
pyseer_Kp1_sCOG_human <- read.csv(paste0(pwd,"KGWP1_Kp1_human_sCOG_output.tsv"), sep="\t")
pyseer_Kp1_sCOG_animal <- read.csv(paste0(pwd,"KGWP1_Kp1_animal_sCOG_output.tsv"), sep="\t")

#SNPs - from SGH10 RedDog alignment
pyseer_Kp1_SNP_human <- read.csv(paste0(pwd,"KGWP1_Kp1_human_SNP_output.tsv"), sep="\t")
pyseer_Kp1_SNP_animal <- read.csv(paste0(pwd,"KGWP1_Kp1_animal_SNP_output.tsv"), sep="\t")

#Unitigs - unitigcaller
pyseer_Kp1_unitig_human <- read.csv(paste0(pwd,"KGWP1_Kp1_human_unitig_output.tsv"), sep="\t")
pyseer_Kp1_unitig_animal <- read.csv(paste0(pwd,"KGWP1_Kp1_animal_unitig_output.tsv"), sep="\t")

#Functions ----
#QQplot with lines added for pyseer-suggested threshold (blue) and manual set threshold (red)
#Based on the qqplot script from https://github.com/mgalardini/pyseer/blob/master/scripts/qq_plot.py
plot_qqplot <- function(pvalues,pyseer_suggsted_threshold,manual_threshold) {
  #  'm' is your data vector in R
  m <- pvalues
  # Prepare the observed data
  y <- -log10(m)
  
  # Generate theoretical quantiles
  x <- -log10(runif(length(m), 0, 1))
  
  # Create a QQ plot
  qqplot(as.numeric(x), as.numeric(y), xlab = "Theoretical Quantiles", ylab = "Observed Quantiles")
  abline(0, 1) # adds a 45-degree reference line
  
  # Add a horizontal line
  # Replace 'intercept_value' with the y-axis value where you want the line
  abline(h = as.numeric(pyseer_suggsted_threshold), col = "blue", lwd = 2)
  abline(h = as.numeric(manual_threshold), col = "red", lwd = 2)
  
}

#Plot qq plots ----
plot_qqplot(pyseer_Kp1_COG_human$lrt.pvalue,5,24)
plot_qqplot(pyseer_Kp1_COG_animal$lrt.pvalue,5,24)
plot_qqplot(pyseer_Kp1_sCOG_human$lrt.pvalue,5,35) 
plot_qqplot(pyseer_Kp1_sCOG_animal$lrt.pvalue,5,35)
plot_qqplot(pyseer_Kp1_unitig_human$lrt.pvalue,8,48)
plot_qqplot(pyseer_Kp1_unitig_animal$lrt.pvalue,8,65)
