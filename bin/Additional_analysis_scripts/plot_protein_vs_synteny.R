#!/usr/bin/env Rscript

# Load necessary libraries
if (!requireNamespace("lattice", quietly = TRUE)) {
  install.packages("lattice")
}
library(lattice)

# Function to process each file
process_file <- function(file_path) {
  # Read the file, skipping the header
  data <- read.table(file_path, header = FALSE, sep = "\t", skip = 1)
  
  # Extract the 8th, 11th, and 17th columns
  extracted_data <- data[, c(8, 11, 17)]
  
  # Rename columns for clarity
  colnames(extracted_data) <- c("Syntenic_Breaks", "Protein_Identity", "Translocation_Junctions")
  
  return(extracted_data)
}

# Main function
plot_syntenic_vs_protein <- function(files, colors = NULL, shapes = NULL) {
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
  
  # Create the xyplot with colors and shapes if provided
  plot_colors <- if (!is.null(colors)) strsplit(colors, ",")[[1]] else NULL
  plot_shapes <- if (!is.null(shapes)) as.numeric(strsplit(shapes, ",")[[1]]) else NULL
  
  # Construct the plot calls conditionally adding colors and shapes
  plot1_call <- list(
    Syntenic_Breaks ~ Protein_Identity,
    data = combined_data,
    groups = combined_data$Source,
    auto.key = list(points = TRUE, text = labels, title = "Source File"),
    xlab = "Protein Identity",
    ylab = "Number of Syntenic Breaks",
    main = "Protein Identity vs. Number of Syntenic Breaks"
  )
  
  plot2_call <- list(
    Translocation_Junctions ~ Protein_Identity,
    data = combined_data,
    groups = combined_data$Source,
    auto.key = list(points = TRUE, text = labels, title = "Source File"),
    xlab = "Protein Identity",
    ylab = "Number of Translocation Junctions",
    main = "Protein Identity vs. Number of Translocation Junctions"
  )
  
  if (!is.null(plot_colors)) {
    plot1_call$col <- plot_colors
    plot2_call$col <- plot_colors
  }
  
  if (!is.null(plot_shapes)) {
    plot1_call$pch <- plot_shapes
    plot2_call$pch <- plot_shapes
  }
  
  # Open a PDF device for output
  pdf("protein_vs_synteny.pdf", width = 11.7, height = 8.3)  # A4 size in inches
  
  # Create the plots
  plot1 <- do.call(xyplot, plot1_call)
  plot2 <- do.call(xyplot, plot2_call)
  
  # Arrange the plots side by side
  print(plot1, position = c(0, 0, 0.5, 1), more = TRUE)
  print(plot2, position = c(0.5, 0, 1, 1))
  
  # Close the PDF device
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
plot_syntenic_vs_protein(input_files, colors, shapes)

