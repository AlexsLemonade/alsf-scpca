read_cellranger_counts <- function(quant_ids, data_dir){
  cellranger_sces <- quant_ids %>%
    purrr::map(
      ~ DropletUtils::read10xCounts(file.path(data_dir, .x,
                                    "outs", "filtered_feature_bc_matrix.h5"),
                                    sample.names = .x,
                                    col.names = TRUE)
    ) %>%
  purrr::map(
    # for consistency with other quantifiers:
    # change the column names just the barcode value, which is the first part of the barcode name
    # drop colData
    function(x) {
      colnames(x) <- stringr::str_extract(colnames(x), "^([ACGT]+)")
      colData(x) <- NULL
      return(x)
    }
  )
  names(cellranger_sces) <- quant_ids
  return(cellranger_sces)
}