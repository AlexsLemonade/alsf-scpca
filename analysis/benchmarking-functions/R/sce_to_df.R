coldata_to_df <- function(sce) {
  df <- as.data.frame(colData(sce)) %>%
    tibble::rownames_to_column(var = "cell_id") %>%
    # add a column to get # of cells detected in that sample
    dplyr::mutate(cells_detected = ncol(sce))
  return(df)
}

rowdata_to_df <- function(sce) {
  df <- as.data.frame(rowData(sce)) %>%
    tibble::rownames_to_column(var = "gene_id")
  return(df)
}
