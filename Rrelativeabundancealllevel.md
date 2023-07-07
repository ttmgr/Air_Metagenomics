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

# Define the strings to filter
strings_to_filter <- c("Homo", "Chordata", "Hominidae")

# Filter the dataframes and remove rows with the specified strings in the 'name' column
all_reports_list <- lapply(all_reports_list, function(df) df[!grepl(paste(strings_to_filter, collapse = "|"), df$name), ])



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

# Initialize an empty list for each rank_code
filtered_lists <- setNames(replicate(length(rank_codes), list(), simplify = FALSE), rank_codes)

# Loop through each rank_code
for (rank_code in rank_codes) {
  # Loop through each report in all_reports_list
  for (i in seq_along(all_reports_list)) {
    # Apply filter_and_calculate_abundance function and save the result to the corresponding list
    filtered_lists[[rank_code]][[i]] <- filter_and_calculate_abundance(all_reports_list[[i]], rank_code)
  }
}


# Generate a color palette with 24 colors
my_palette <- c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8,"Dark2") )

# Add "Rare" color (gray) at the start of the palette
my_palette <- c("gray", my_palette)

# Taxonomic level mapping
taxonomic_level_names <- c("D" = "Domain", "P" = "Phylum", "C" = "Class", "O" = "Order", 
                           "F" = "Family", "G" = "Genus", "S" = "Species")


# Define function to calculate and print median and standard deviation
calc_stats <- function(df){
  median_val <- median(df$relative_abundance)
  sd_val <- sd(df$relative_abundance)
  
  cat(paste0('Median of "Rare" group: ', median_val, '\n'))
  cat(paste0('Standard deviation of "Rare" group: ', sd_val, '\n'))
}

# Initialize an empty color mapping
# Initialize an empty color mapping with grey color for 'Rare'
color_mapping <- list("Rare" = "gray")
#color_mapping <- list()

# Update the perform_operations function
perform_operations <- function(reports_list, rank_code, top_n, include_rare = TRUE, x_tick_labels = NULL, plot_title = NULL){  
  
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
  if (include_rare) {
    rare_level <- level_df %>%
      anti_join(top_level, by = c("barcode", "name")) %>%
      group_by(barcode) %>%
      summarise(name = "Rare", relative_abundance = sum(relative_abundance), .groups = "drop")
    
    # Bind top levels and rare
    final_level <- bind_rows(top_level, rare_level)
  } else {
    final_level <- top_level
  }
  
  # Recalculate relative abundance for final_level
  final_level <- final_level %>%
    group_by(barcode) %>%
    mutate(relative_abundance = relative_abundance / sum(relative_abundance)) %>%
    ungroup()
  
  # Calculate and print statistics for "Rare" group
  if (include_rare) {
    calc_stats(final_level[final_level$name == "Rare",])
  }
  # Trim whitespace from names
  final_level$name <- str_trim(final_level$name)
  
  # Ensure "Rare" comes first in the legend
  final_level$name <- factor(final_level$name, levels = c("Rare", unique(final_level$name[final_level$name != "Rare"])))
  
  # Extract unique names in the current level excluding 'Rare'
  unique_names <- setdiff(unique(final_level$name), "Rare")
  
  # Assign colors to new names
  new_names <- setdiff(unique_names, names(color_mapping))
  if (length(new_names) > 0) {
    available_colors <- setdiff(my_palette, color_mapping)
    color_mapping[new_names] <- available_colors[seq_along(new_names)]
  }
  
  # Alphabetically sort names for legend (excluding 'Rare')
  sorted_names <- sort(setdiff(unique(final_level$name), "Rare"))
  
  # If 'Rare' should be included, it is added as first element
  if(include_rare) {
    sorted_names <- c("Rare", sorted_names)
  }
  
  # Relevel factor levels for ordered legend
  final_level$name <- factor(final_level$name, levels = sorted_names)
  
  # Create a stacked bar plot
  plot <- ggplot(final_level, aes(fill = name, y = relative_abundance, x = barcode)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    labs(x = "Barcode", y = "Relative Abundance (%)", fill = taxonomic_level_names[rank_code], title = plot_title) +
    scale_fill_manual(values = color_mapping) +
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
  
  # If x_tick_labels is not NULL, change the x-axis tick labels
  if (!is.null(x_tick_labels)) {
    plot <- plot + scale_x_discrete(labels = x_tick_labels)
  }
  
  print(plot)
}
# Call the function with the desired parameters
perform_operations(filtered_lists[["G"]], "G", 16, include_rare = T, x_tick_labels = c("BC01", "BC02", "BC03"), plot_title = "My Plot")

