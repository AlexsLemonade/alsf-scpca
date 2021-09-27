#' Create a list of SingleCellExperiment objects using scpcaTools::import_quant_data. 
#'
#' @param info_df data.frame containing information for each sample to create a SingleCellExperiment object for. 
#'   Each row corresponds to an individual sample and columns contain information for that sample: 
#'   data_dir - Local directory where the output quants data for all samples is stored. 
#'   tool - which tool was used to process the samples. 
#'   which_counts - Which counts to keep in the counts matrix, cDNA only ("spliced") or cDNA and introns ("unspliced").
#'   intron_mode - Whether or not the data was aligned to the cDNA only (FALSE) or cDNA and introns index (TRUE).
#'   usa_mode - If the data was processed with alevin-fry, was it processed using USA mode. 
#'   filter - Logical indicating if the returned SingleCellExperiment should contain filtered (TRUE) or unfiltered (FALSE) data.
#' @param mito List of mitochondrial genes used to calculate per cell QC mitochondrial fraction. 
#' @param filename Optional filename used to save list of SingleCellExperiments as an .rds file.
#'
#' @return List of SingleCellExperiments.
#'
make_sce_list <- function(info_df, mito, filename = NULL){
  # for each row in the info data frame, create a SingleCellExperiment object and save as a list
  sce_list <- mapply(scpcaTools::import_quant_data,
                     info_df$data_dir,
                     info_df$tool,
                     info_df$which_counts,
                     FALSE,
                     info_df$intron_mode,
                     info_df$usa_mode,
                     info_df$filter) 
  # calcluate per cell and per feature qc 
  sce_list <- sce_list %>%
    purrr::map(scpcaTools::add_cell_mito_qc, mito) %>%
    purrr::map(
      scater::addPerFeatureQC
    )
  # use the unique quant directory (contains sample + run parameters as the name) to name each sce in the list
  names(sce_list) <- info_df$quant_dir
  
  # if the filename is present, save as rds file
  if(!is.null(filename)){
    readr::write_rds(sce_list, file = filename)
  }
  return(sce_list)
}
