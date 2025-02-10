library(tidyverse)

# Define the directory where KmerFinder files are stored
kmerfinder_dir <- "kmerfinder_results"

# List all KmerFinder result files
kmerfinder_files <- list.files(kmerfinder_dir, pattern = "_kmerfinder_results.txt$", full.names = TRUE)

# Function to process a single KmerFinder file
process_kmerfinder_file <- function(file_path) {
  # Extract strain name from filename
  strain <- str_remove(basename(file_path), "_kmerfinder_results.txt")
  
  # Read in the file
  df <- read_delim(file_path, delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
    select(Query_Coverage, `TAXID Species`) %>%  # Keep only the relevant columns
    mutate(
      Query_Coverage = as.numeric(Query_Coverage),  # Ensure numeric type
      `TAXID Species` = str_squish(`TAXID Species`),  # Clean up species names
      genome = strain
    )
  
  # Check for any NA values in Query_Coverage due to conversion issues
  if (any(is.na(df$Query_Coverage))) {
    warning(paste("NA values found in Query_Coverage for", strain, "- check file formatting"))
  }
  
  # Summarize query coverage per species
  species_summary <- df %>%
    group_by(`TAXID Species`) %>%
    summarise(total_query_coverage = sum(Query_Coverage, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_query_coverage))
  
  # Identify the main species (highest query coverage)
  top_species_row <- species_summary %>% slice(1)
  top_species <- top_species_row$`TAXID Species`
  top_species_coverage <- top_species_row$total_query_coverage
  
  # Sum query coverage of all other species (contamination)
  num_contamination_coverage <- species_summary %>%
    filter(`TAXID Species` != top_species) %>%
    summarise(total = sum(total_query_coverage, na.rm = TRUE)) %>%
    pull(total)
  
  # Compute % contamination
  total_coverage <- top_species_coverage + num_contamination_coverage
  percent_contamination <- ifelse(total_coverage > 0, 
                                  (num_contamination_coverage / total_coverage) * 100, 
                                  0)
  
  # Return as a dataframe
  tibble(
    strain = strain,
    top_species = top_species,
    top_species_coverage = top_species_coverage,
    num_contamination_coverage = num_contamination_coverage,
    percent_contamination = percent_contamination
  )
}

# Process all files and combine results
kmerfinder_results_df <- map_dfr(kmerfinder_files, process_kmerfinder_file)
kmerfinder_results_df

write.csv(kmerfinder_results_df,"kmerfinder_results_df.csv")