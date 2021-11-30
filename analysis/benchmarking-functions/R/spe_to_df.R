#' Combine colData and in_tissue column of SpatialData from SpatialExperiment into a data.frame
#'
#' @param spe SpatialExperiment object with per spot QC metrics in colData and 
#'  in_tissue column present in spatialData
#'
#' @return data.frame containing both colData and in_tissue from spatialData of original spe
spatial_coldata_to_df <- function(spe){
  
  # grab coldata and add in tissue column
  col_df <- as.data.frame(colData(spe)) %>%
    tibble::rownames_to_column(var = "spot_id") %>%
    dplyr::mutate(in_tissue = spatialData(spe)$in_tissue)
}
