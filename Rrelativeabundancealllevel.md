```shell
library(readr)
library(stringr)
library(httr)
library(XML)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

# Define the directory
dir_path <- 'Desktop/Coriolis_Greenhouse/greenhouse23h/readbased/'

# Get a list of all 'read_report.txt' files
files <- list.files(path = dir_path, pattern = '*read_report.txt$', full.names = TRUE)

# Barcode order
barcode_order <- c("barcode03", "barcode04", "barcode05", "barcode08", "barcode09", "barcode10", 
                   "barcode13", "barcode14", "barcode15", "barcode20", "barcode21", "barcode22", 
                   "barcode23", "barcode24", "barcode01")

# Loop through each file, read the data, and perform the analysis
all_reports_list <- lapply(files, function(file) {
  kraken_report <- read.table(file, header = F, sep = '\t')
  colnames(kraken_report) <- c('percentage', 'count', 'taxon_count', 'rank_code', 'NCBI_taxon_ID', 'name')
  kraken_report$relative_abundance <- kraken_report$count / sum(kraken_report$count)
  
  # Add a new column 'barcode' extracted from the filename
  kraken_report$barcode <- sub(".*\\/(barcode\\d+).*", "\\1", file)
  
  # Return the dataframe
  return(kraken_report)
})




# Define function for filtering and relative abundance calculation
filter_and_calculate_abundance <- function(df, rank_code) {
  filtered_df <- df[df$rank_code == rank_code,]
  filtered_df$relative_abundance <- filtered_df$count / sum(filtered_df$count)
  
  # Reorder barcodes according to given barcode order
  filtered_df$barcode <- factor(filtered_df$barcode, levels = barcode_order)
  
  return(filtered_df)
}

# Define rank codes for taxonomic levels
rank_codes <- c("D", "P", "C", "O", "F", "G", "S")

# Use lapply to apply the function to every data frame in our list
all_reports_list_filtered <- lapply(all_reports_list, function(df) {
  lapply(rank_codes, function(rank_code) {
    filter_and_calculate_abundance(df, rank_code)
  })
})

# Generate a color palette with 24 colors
my_palette <- c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1") )

# Add "Rare" color (gray) at the start of the palette
my_palette <- c("gray", my_palette)

# Function to perform operations
perform_operations <- function(reports_list, rank_code, level_name, top_n){
  
  # Filter dataframes for each taxonomic level
  level_list <- lapply(reports_list, function(df) df[df$rank_code == rank_code,])
  
  # Combine all 'level' dataframes from the list
  level_df <- do.call(rbind, level_list)
  
  # Recalculate relative abundance for each barcode
  level_df <- level_df %>%
    group_by(barcode) %>%
    mutate(relative_abundance = count / sum(count)) %>%
    ungroup()
  
  # Filter to only include the top_n most abundant levels for each barcode
  top_level <- level_df %>%
    group_by(barcode, name) %>%
    summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
    arrange(barcode, desc(relative_abundance)) %>%
    group_by(barcode) %>%
    slice_max(order_by = relative_abundance, n = top_n) %>%
    ungroup()
  
  # Collect remaining levels for each barcode as 'Rare'
  rare_level <- level_df %>%
    anti_join(top_level, by = c("barcode", "name")) %>%
    group_by(barcode) %>%
    summarise(name = "Rare", relative_abundance = sum(relative_abundance), .groups = "drop")
  
  # Bind top levels and rare
  final_level <- bind_rows(top_level, rare_level)
  
  # Trim whitespace from names
  final_level$name <- str_trim(final_level$name)
  
  # Ensure "Rare" comes first in the legend
  final_level$name <- factor(final_level$name, levels = c("Rare", unique(final_level$name[final_level$name != "Rare"])))
  
  # Create a stacked bar plot
  plot <- ggplot(final_level, aes(fill = name, y = relative_abundance, x = barcode)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    labs(x = "Barcode", y = "Relative Abundance (%)", fill = level_name) +
    scale_fill_manual(values = my_palette) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      legend.box = "horizontal"
    )
  
  print(plot)
}

# Apply the function to each taxonomic level
taxonomic_levels <- c("D", "P", "C", "O", "F", "G", "S")
level_names <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
top_n_values <- c(10, 10, 10, 10, 10, 15, 14)

for (i in 1:length(taxonomic_levels)) {
  perform_operations(all_reports_list, taxonomic_levels[i], level_names[i], top_n_values[i])
}
