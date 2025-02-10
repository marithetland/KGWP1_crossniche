library(tidyverse)

# Define the directory where Kraken2 files are stored
kraken2_dir <- "kraken2_results"

# List all Kraken2 result files
kraken2_files <- list.files(kraken2_dir, pattern = "__kraken2.txt$", full.names = TRUE)

# Function to process a single Kraken2 file
process_kraken2_file <- function(file_path) {
  # Extract strain name from filename
  strain <- str_remove(basename(file_path), "__kraken2.txt")
  
  # Read in the file with predefined column names
  df <- read_delim(file_path, delim = "\t", col_names = FALSE, show_col_types = FALSE) %>%
    rename(
      percent_fragments = X1,
      num_fragments_covered = X2,
      num_fragments_directly_assigned = X3,
      num_minimizers = X4,
      est_minimizers = X5,
      rank_code = X6,
      ncbi_taxID = X7,
      indented_scientific_name = X8
    ) %>%
    mutate(
      genome = strain,
      indented_scientific_name = str_squish(indented_scientific_name) # Remove extra spaces
    )
  
  # Filter for species-level classifications
  species_df <- df %>% filter(rank_code == "S")
  
  # Identify the top species (highest num_fragments_directly_assigned)
  top_species_row <- species_df %>%
    arrange(desc(num_fragments_directly_assigned)) %>%
    slice(1)
  
  # Get top species name and fragment count
  top_species <- top_species_row$indented_scientific_name
  top_species_num_fragments <- top_species_row$num_fragments_directly_assigned
  
  # Sum fragments assigned to other species (contamination)
  num_contamination_fragments <- species_df %>%
    filter(indented_scientific_name != top_species) %>%
    summarise(total = sum(num_fragments_directly_assigned, na.rm = TRUE)) %>%
    pull(total)
  
  # Compute % contamination
  total_species_fragments <- species_df %>%
    summarise(total = sum(num_fragments_directly_assigned, na.rm = TRUE)) %>%
    pull(total)
  
  percent_contamination <- ifelse(total_species_fragments > 0, 
                                  (num_contamination_fragments / total_species_fragments) * 100, 
                                  0)
  
  # Return as a dataframe
  tibble(
    strain = strain,
    top_species = top_species,
    top_species_num_fragments = top_species_num_fragments,
    num_contamination_fragments = num_contamination_fragments,
    percent_contamination = percent_contamination
  )
}

# Process all files and combine results
kraken2_results_df <- map_dfr(kraken2_files, process_kraken2_file)

write.csv(kraken2_results_df,"kraken2_results_df.csv")