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
  paste0("Go_summary_", file_name, "_topSynteny.tsv"),
  paste0("Go_summary_", file_name, "_botSynteny.tsv"),
  paste0("Go_summary_", file_name, "_averhigh.tsv"),
  paste0("Go_summary_", file_name, "_averlow.tsv"),
  paste0("Go_summary_", file_name, "_highScore.tsv"),
  paste0("Go_summary_", file_name, "_lowScore.tsv"),
  paste0("Go_summary_", file_name, "_top_orthologous.tsv"),
  paste0("Go_summary_", file_name, "_bot_orthologous.tsv")
)

# Iterate over each file path
for (file_path in file_paths) {
  # Read the data
  erefd <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Process data as per your original script
  newdata <- erefd[order(erefd$Count_significant, decreasing = TRUE), ]
  rownames(newdata) <- paste(newdata$GO_ID, newdata$GO_term)
  df <- subset(newdata, select = -c(GO_term, GO_ID, Count_significant))
  df[is.na(df)] <- 1
  df2 <- as.matrix(df)
  
  # Check mean
  my_mean <- mean(df2)
  
  # Plot if mean is not 1
  if (my_mean != 1) {
    my_palette <- colorRampPalette(c("purple", "red", "orange", "yellow", "white"))(n = 20)
    pdf(paste0("Go_summary_", file_name, "_", gsub(".tsv", ".pdf", basename(file_path))),
        width = 6, height = 5)
    pheatmap::pheatmap(log10(head(df2, n = 30)), col = my_palette, cluster_rows = FALSE, treeheight_row = 0, treeheight_col = 0, legend = TRUE)
    dev.off()
  }
}
