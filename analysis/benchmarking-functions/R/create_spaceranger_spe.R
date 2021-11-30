#' Import Spaceranger output for Spatial Transcriptomics library into a SpatialExperiment object.
#'
#' @param spaceranger_dir Path to folder output by Spaceranger.
#' @param sample_name Sample id.
#'
#' @return SpatialExperiment object.
create_spaceranger_spe <- function(spaceranger_dir, 
                                       sample_name){
  # path to counts dir
  counts_dir <- file.path(spaceranger_dir, "outs", "raw_feature_bc_matrix")
  
  # path to spatial data 
  spatial_dir <- file.path(spaceranger_dir, "outs", "spatial")
  
  # read in counts as an sce 
  sce <- DropletUtils::read10xCounts(counts_dir, col.names = TRUE)
  
  # read in image data
  image_data <- SpatialExperiment::readImgData(spatial_dir,
                                               sample_id = sample_name)
  
  # read in spatial coordinates 
  spatial_data_file <- file.path(spatial_dir, "tissue_positions_list.csv")
  spatial_data <- readr::read_csv(spatial_data_file, col_names = c("barcode", "in_tissue", "array_row",
                                                                   "array_col","pxl_row_in_fullres", "pxl_col_in_fullres"))
  
  # construct 'SpatialExperiment'
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = assay(sce)),
    colData = SingleCellExperiment::colData(sce), 
    rowData = DataFrame(symbol = rowData(sce)$Symbol),
    imgData = image_data,
    spatialData = DataFrame(spatial_data),
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    sample_id = sample_name)
  
  # return SpatialExperiment object
  return(spe)
}

