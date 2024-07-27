# Install Library

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}

# Load Library

library(SingleCellExperiment)

# process_panglaodb function definition

process_panglaodb <- function(data_path) {
  load(data_path)

  object <- get(ls()[1])

  # counts : sparse matrix
  # colData : metadata about cells
  # rowData : metadata about genes

  data_matrix <- object$counts
  col_data <- object$colData
  row_data <- object$rowData

  col_data_df <- as.data.frame(col_data)
  row_data_df <- as.data.frame(row_data)

  return(list(counts_matrix = data_matrix, col_data_df = col_data_df, row_data_df = row_data_df)) # nolint
}