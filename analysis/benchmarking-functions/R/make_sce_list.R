addCellQC <- function(sce, mito, threshold = 0){
  scater::addPerCellQC(sce, 
                       subsets = list(mito = mito[mito %in% rownames(sce)]),
                       threshold = threshold)
}


make_sce_list <- function(info_df, filename, mito){
  sce_list <- mapply(scpcaTools::import_quant_data,
                     info_df$data_dir,
                     info_df$tool,
                     info_df$which_counts,
                     info_df$intron_mode,
                     info_df$usa_mode,
                     info_df$filter) 
  sce_list <- sce_list %>%
    purrr::map(addCellQC, mito) %>%
    purrr::map(
      scater::addPerFeatureQC
    )
  names(sce_list) <- info_df$quant_dir
  readr::write_rds(sce_list, filename)
  return(sce_list)
}
