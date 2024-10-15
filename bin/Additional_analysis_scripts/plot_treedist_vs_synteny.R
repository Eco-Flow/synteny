#!/usr/bin/env Rscript

# Load necessary libraries
if (!requireNamespace("lattice", quietly = TRUE)) {
  install.packages("lattice")
}
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(lattice)
library(gridExtra)
library(dplyr)


# Function to normalize comparison strings
normalize_comparison <- function(comparison) {
  species <- unlist(strsplit(comparison, "\\."))
  return(paste(sort(species), collapse = "."))
}


# Function to process each file
process_file <- function(file_path) {
  # Read the file, skipping the header
  data <- as.data.frame(read.table(file_path, header = FALSE, sep = "\t", skip = 1))
  
  # Extract the 8th (syntenic blocks), 17th (translocations minimum) columns and 19th (tree distance)
  extracted_data <- data[, c(1, 8, 17, 19)]

  # Add a combined name column
  extracted_data$Normalized_V1 <- sapply(extracted_data$V1, normalize_comparison)
  
  # Average pairs:
  averaged_df <- extracted_data %>%
  dplyr::group_by(Normalized_V1) %>%
  dplyr::summarise(
    Avg_V8 = mean(V8),
    Avg_V17 = mean(V17),
    Avg_V19 = mean(V19)
  ) %>%
  dplyr::select(Avg_V8, Avg_V17, Avg_V19)

  # Rename columns for clarity
  colnames(averaged_df) <- c("Syntenic_Breaks", "Translocation_Junctions", "Matrix_Value")
  
  return(averaged_df)
}

# Function to create the plots
create_plots <- function(files, colors = NULL, shapes = NULL) {
  # Initialize a list to hold data from each file
  all_data <- list()
  labels <- character(length(files))  # To store labels for the key
  
  # Process each file and store its data with a source identifier
  for (i in seq_along(files)) {
    file_data <- process_file(files[i])
    
    # Get the file name without the path and extension
    file_label <- gsub(".tsv$", "", basename(files[i]))
    labels[i] <- file_label
    
    file_data$Source <- file_label  # Add a column to indicate source file
    all_data[[i]] <- file_data
  }
  
  # Combine all the data into one data frame
  combined_data <- do.call(rbind, all_data)
  combined_data$Source <- factor(combined_data$Source, levels = unique(combined_data$Source))

  # Write the table with the Comparison column
  write.table(combined_data, "My_averaged_summary_with_comparison.tsv", sep = "\t", row.names = FALSE)
  
  # Create the xyplot with colors and shapes if provided
  plot_colors <- if (!is.null(colors)) strsplit(colors, ",")[[1]] else NULL
  plot_shapes <- if (!is.null(shapes)) as.numeric(strsplit(shapes, ",")[[1]]) else NULL
  
  # Define the settings for symbols
  symbol_settings <- list(col = plot_colors, pch = plot_shapes, fill = plot_colors)
  
  # Create the first plot
  plot1 <- xyplot(Syntenic_Breaks ~ Matrix_Value, data = combined_data,
                  groups = combined_data$Source,
                  auto.key = FALSE, # Disable auto key for custom key placement
                  xlab = "Tree Distance",
                  ylab = "Number of Syntenic Breaks",
                  main = "Tree Distance vs. Number of Syntenic Breaks",
                  par.settings = list(superpose.symbol = symbol_settings))
  
  # Create the second plot
  plot2 <- xyplot(Translocation_Junctions ~ Matrix_Value, data = combined_data,
                  groups = combined_data$Source,
                  auto.key = FALSE, # Disable auto key for custom key placement
                  xlab = "Tree Distance",
                  ylab = "Number of Translocation Junctions",
                  main = "Tree Distance vs. Number of Translocation Junctions",
                  par.settings = list(superpose.symbol = symbol_settings))
  
  # Open a PDF device for the plots
  pdf("plots_tree.pdf", width = 11.7, height = 8.3)  # A4 size in inches
  
  # Arrange the plots side by side
  grid.arrange(plot1, plot2, ncol = 2)
  
  # Close the PDF device for the plots
  dev.off()
  
  # Create a separate plot for the legend
  legend_data <- data.frame(
    Source = labels,
    Color = plot_colors,
    Shape = plot_shapes
  )
  
  legend_plot <- xyplot(1 ~ 1, type = "n", xlab = "", ylab = "",
                        xlim = c(0, 2), ylim = c(0, length(labels) + 1),
                        main = "Key",
                        panel = function(...) {
                          panel.grid()
                          for (i in seq_along(labels)) {
                            panel.points(0.5, length(labels) - i + 1, 
                                         col = legend_data$Color[i],
                                         pch = legend_data$Shape[i],
                                         cex = 2)  # Increase size of symbols
                            panel.text(1.5, length(labels) - i + 1,
                                       labels[i], pos = 4, cex = 1, font = 2)  # Make text bold
                          }
                        })
  
  # Open a PDF device for the legend
  pdf("legend_tree.pdf", width = 6, height = 6)  # Increase size if needed
  
  # Print the legend plot
  print(legend_plot)
  
  # Close the PDF device for the legend
  dev.off()
}

# Check if the script is run with at least one file argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide at least one input file as a command line argument.")
}

# Parse input arguments
input_files <- args[1:sum(file.exists(args))]  # Select only existing files
colors <- if (length(args) > length(input_files)) args[length(input_files) + 1] else NULL
shapes <- if (length(args) > length(input_files) + 1) args[length(input_files) + 2] else NULL

# Call the main function with the provided command line arguments
create_plots(input_files, colors, shapes)


