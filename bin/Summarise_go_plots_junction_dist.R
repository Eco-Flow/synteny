#!/usr/bin/Rscript

library(pheatmap)

# Check if the correct number of arguments are provided
if (length(commandArgs(TRUE)) != 1) {
  stop("Usage: script_name.R <file_name>")
}

# Read the file name argument from the command line
file_name <- commandArgs(TRUE)[1]

# Construct the full file paths
file_paths <- c(
  paste0("Go_summary_", file_name, "_closest.tsv"),
  paste0("Go_summary_", file_name, "_farthest.tsv"),
  paste0("Go_summary_", file_name, "_all_merged.tsv")
)

#print("Reached here\n")

# Iterate over each file path
for (file_path in file_paths) {
  # Read the data
  erefd <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Process data as per your original script
  newdata <- erefd[order(erefd$Count_significant, decreasing = TRUE), ]
  rownames(newdata) <- paste(newdata$GO_ID, newdata$GO_term)
  
  # Conditionally remove columns if they exist
  cols_to_remove <- c("GO_term", "GO_ID", "Count_significant", "Count_closest", "Count_farthest")
  cols_to_remove <- intersect(cols_to_remove, colnames(newdata))
  
  df <- newdata[, !(colnames(newdata) %in% cols_to_remove)]
  
  # Replace NAs with 1
  df[is.na(df)] <- 1  
  df2 <- as.matrix(df)

  # Check for NaN or Inf
  if (any(is.na(df2)) || any(is.infinite(df2))) {
    print("Matrix contains NA/NaN/Inf values:")
    print(df2)
  }
  
  # Check mean
  my_mean <- mean(df2)
  
  # Plot if mean is not 1
  if (my_mean != 1) {
    # Add a small constant to avoid log10(0) issues
    df2 <- df2 + 1e-10
    
    # Check again for invalid values after adjustment
    if (any(is.na(df2)) || any(is.infinite(df2))) {
      print("Matrix still contains NA/NaN/Inf values after adjustment:")
      print(df2)
      next  # Skip plotting for this file if invalid values exist
    }
    
    my_palette <- colorRampPalette(c("purple", "red", "orange", "yellow", "white"))(n = 20)
    pdf(paste0("Go_summary_", file_name, "_", gsub(".tsv", ".pdf", basename(file_path))),
        width = 6, height = 5)
    
    pheatmap::pheatmap(log10(head(df2, n = 30)), col = my_palette, cluster_rows = FALSE, treeheight_row = 0, treeheight_col = 0, legend = TRUE)
    dev.off()
  }
}

