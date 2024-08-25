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
  if (!file.exists(data_path)) {
    stop("File does not exist: ", data_path)
  }
  loaded_objects <- load(data_path)
  
  if (length(loaded_objects) == 0) {
    stop("No objects loaded from RData file.")
  }
  
  # 로드된 객체의 이름과 구조를 확인
  # 로드된 객체의 이름과 구조를 확인
  object_name <- loaded_objects[1]  # 첫 번째 객체 이름 가져오기
  object <- get(object_name)  # 객체 가져오기

  print(paste("Loaded object name:", object_name))
  print(str(object))  # 객체의 구조 출력
  # counts : sparse matrix
  # colData : metadata about cells
  # rowData : metadata about genes
  # 객체가 SingleCellExperiment인지 확인
  if (inherits(object, "SingleCellExperiment")) {
    data_matrix <- object$counts
    col_data <- colData(object)
    row_data <- rowData(object)

    col_data_df <- as.data.frame(col_data)
    row_data_df <- as.data.frame(row_data)

    return(list(counts_matrix = data_matrix, col_data_df = col_data_df, row_data_df = row_data_df))
  } else if (inherits(object, "dgCMatrix")) {
    # dgCMatrix 객체를 처리하는 코드
    data_matrix <- as.matrix(object)  # dgCMatrix 객체를 numpy 배열로 변환
    col_data_df <- as.data.frame(colnames(object))
    row_data_df <- as.data.frame(rownames(object))
    
    return(list(counts_matrix = data_matrix, col_data_df = col_data_df, row_data_df = row_data_df))
  } else {
    stop("Unsupported object type: ", class(object))
  }
}