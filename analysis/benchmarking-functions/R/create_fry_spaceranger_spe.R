#' Import Alevin-fry and Spaceranger output for Spatial Transcriptomics library into a SpatialExperiment object.
#'
#' @param fry_dir Path to alevin-fry outputs directory. 
#' @param spaceranger_dir Path to folder output by Spaceranger. 
#' @param sample_name Sample id. 
#'
#' @return SpatialExperiment object with shared spots detected in Alevin-fry and Spaceranger.
#'
create_fry_spaceranger_spe <- function(fry_dir, 
                                       spaceranger_dir, 
                                       sample_name){
  # path to spatial output folder 
  spatial_dir <- file.path(spaceranger_dir, "outs", "spatial")
  
  # create SingleCellExperiment using scpcaTools
  fry_sce <- scpcaTools::read_alevin(fry_dir, 
                                     usa_mode = TRUE,
                                     which_counts = "spliced")
  
  # read in image data
  image_data <- SpatialExperiment::readImgData(spatial_dir,
                                               sample_id = sample_name)
  
  # read in spatial coordinates 
  spatial_data_file <- file.path(spatial_dir, "tissue_positions_list.csv")
  spatial_data <- readr::read_csv(spatial_data_file, col_names = c("barcode", "in_tissue", "array_row",
                                                                   "array_col","pxl_row_in_fullres",
                                                                   "pxl_col_in_fullres"))
  # find common cells between Alevin-fry and spaceranger 
  spatial_data <- spatial_data %>%
    dplyr::mutate(barcode = gsub("-1", "", barcode))
  common_cells <- intersect(colnames(fry_sce), spatial_data$barcode)
  
  # subset fry_sce by common cells
  fry_sce_subset <- fry_sce[,common_cells]
  spatial_data_subset <- data.frame(barcode = common_cells) %>%
    dplyr::left_join(spatial_data)
  
  # construct 'SpatialExperiment'
  fry_spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = assay(fry_sce_subset)),
    colData = SingleCellExperiment::colData(fry_sce_subset), 
    imgData = image_data,
    spatialData = DataFrame(spatial_data_subset),
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    sample_id = sample_name)
  
  # return SpatialExperiment object
  return(fry_spe)
}
