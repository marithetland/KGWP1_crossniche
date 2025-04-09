#####################################
# KGWP1 paper
# 2023-03-28, Hetland MAK
# Fig S13
#####################################

#Set env ----
set.seed(123)
pwd="~/KGWP1/BactDating_July_2024/"

#Load pkgs ----
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(lubridate)
library(tibble)


#Load tip dates ----
list_SLs <- c("SL107","SL3010","SL25","SL34","SL111","SL29","SL290","SL35","SL45","SL641")

metadata_all <- read.csv("~/KGWP1/genotyping/metadata_all_2024-05-29.csv")

SL_metadata <- metadata_all %>% filter(SL %in% list_SLs) %>% select(strain,SL,Collection_Year,Niche,source,virulence_score,resistance_score,num_resistance_classes)
names(SL_metadata) <- c("tip","SL","date","niche","source","virulence_score","resistance_score","num_resistance_classes")

#Load trees ----
tree_files <- c("SL107__BactDating_verticall.tree",
                #"SL200__BactDating_verticall.tree",
                "SL3010__BactDating_verticall.tree",
                #"SL3676__BactDating_verticall.tree",
                #"SL10__BactDating_verticall.tree",
                "SL25__BactDating_verticall.tree",
                "SL34__BactDating_verticall.tree",
                #"SL461__BactDating_verticall.tree",
                "SL111__BactDating_verticall.tree",
                "SL29__BactDating_verticall.tree",
                "SL290__BactDating_verticall.tree",
                #"SL35__BactDating_verticall.tree",
                "SL45__BactDating_verticall.tree") #,
                #"SL641__BactDating_verticall.tree")

#Set colour schemes ----

#Create tip colors with names
tip_colours <- c(
  "Human_infection"="black", "Community_carriage"="black",
  "Turkey"="black", "Broiler"="black", "Pig"="black",
  "Dog" = "black", "Bivalves" = "black")

#Create heatmap colors with names
heatmap_colours <- c(
  "Human_infection"="#CC004C", "Community_carriage"="#a38b46",
  "Turkey"="#6460AA", "Broiler"="#0089D0", "Pig"="#FCD7DE",
  "Dog" = "#0DB14B", "Bivalves" = "#00008B",
  #Set rainbow colours:
  "0"="white", "1"="#e60718", "2"="#e0863d", "3"="#F9C74F", "4"="#90BE6D",
  "5"="#56CFE1", "6"="#023E8A", "7"="#9C89B8", "8"="#000000")

#For each tree, plot
#Plot local trees
# Define the function
generate_tree_heatmap <- function(SL, tree, metadata, tip_colours, heatmap_colours, pwd, x_breaks, x_limits = NULL) {
  # Read metadata and edit tip colname
  tree_metadata <- metadata %>%
    slice(match(tree$tip.label, tip))  # Order same as tip labels, will also filter by those matching the tree tips
  
  # Join metadata to tree tips, to ensure correct labels
  tree_data <- as_tibble(tree)
  tree_data <- tree_data %>%
    left_join(tree_metadata, by = c("label" = "tip")) %>%
    filter(!is.na(label)) %>%  # Remove non-tips
    as.data.frame()
  
  # Plot
  p_tree <- ggtree(tree, right = TRUE, mrsd = date_decimal(max(tree_data$date))) %<+% 
    tree_data +
    geom_tippoint(aes(color = factor(source)), size = 2) +
    scale_color_manual(values = tip_colours, name = "source") + 
    theme_tree2() +  # Adds a timescale to the plot
    scale_x_continuous(breaks = x_breaks, limits = x_limits)
  
  # Prepare heatmap data
  heatmap_data <- tree_data %>%
    select(source, virulence_score, resistance_score, num_resistance_classes)
  rownames(heatmap_data) <- tree_data$label
  
  # Add heatmap to the tree plot
  p_heatmap <- gheatmap(p_tree, heatmap_data, offset = 0.005, width = 1, font.size = 3, color = "white", colnames_position = "bottom", colnames_angle = 90, hjust = 0) +
    scale_fill_manual(values = heatmap_colours, breaks = names(heatmap_colours))
  
  # Save the plot as PDF
  pdf_filename <- paste0(pwd, "local/", SL, "_tree_heatmap.pdf")
  ggsave(pdf_filename, plot=p_heatmap, dpi = 300, width = 18, height = 12, units = "cm")
  
  return(p_heatmap)
}

#Plot for each tree
## SL107
SL="SL107"
tree <- read.nexus(paste0(pwd,"/local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1970, 2020, by = 10))

## SL25
SL="SL25"
tree <- read.nexus(paste0(pwd,"local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1970, 2020, by = 10))

## SL3010
SL="SL3010"
tree <- read.nexus(paste0(pwd,"local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1630, 2020, by = 50))

## SL34
SL="SL34"
tree <- read.nexus(paste0(pwd,"local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1300, 2020, by = 50))

## SL29
SL="SL29"
tree <- read.nexus(paste0(pwd,"local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1810, 2020, by = 50))

## SL290
SL="SL290"
tree <- read.nexus(paste0(pwd,"local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1950, 2020, by = 10))

## SL111
SL="SL111"
tree <- read.nexus(paste0(pwd,"local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1700, 2020, by = 50))

## SL45
SL="SL45"
tree <- read.nexus(paste0(pwd,"local/",SL,"__BactDating_verticall.tree"))
generate_tree_heatmap(SL, tree, SL_metadata, tip_colours, heatmap_colours, pwd, seq(1300, 2020, by = 100))
