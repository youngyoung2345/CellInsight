# Install Library

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}
if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
  BiocManager::install("GenomeInfoDbData")
}

# Load Library

library(GenomeInfoDb)
library(SingleCellExperiment)

# process_panglaodb function definition

process_panglaodb <- function(data_path) {
  loaded_objects <- load(data_path)

  object_name <- loaded_objects[1]
  object <- get(object_name)

  # counts : sparse matrix
  # colData : metadata about cells
  # rowData : metadata about genes

  if (inherits(object, "SingleCellExperiment")) {
    data_matrix <- object$counts
    col_data <- object$colData
    row_data <- object$rowData

    col_data_df <- as.data.frame(col_data)
    row_data_df <- as.data.frame(row_data)

    return(list(counts_matrix = data_matrix, col_data_df = col_data_df, row_data_df = row_data_df)) # nolint
  } else if (inherits(object, "dgCMatrix")) {
    data_matrix <- as.matrix(object)
    col_data_df <- as.data.frame(colnames(object))
    row_data_df <- as.data.frame(rownames(object))

    return(list(counts_matrix = data_matrix, col_data_df = col_data_df, row_data_df = row_data_df)) # nolint
  }
}