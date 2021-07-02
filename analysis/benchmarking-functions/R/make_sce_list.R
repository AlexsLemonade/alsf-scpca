make_sce_list <- function(info_df, filename, mito){
  sce_list <- mapply(scpcaTools::import_quant_data,
                     info_df$data_dir,
                     info_df$tool,
                     info_df$which_counts,
                     info_df$intron_mode,
                     info_df$usa_mode,
                     info_df$filter) 
  sce_list <- sce_list %>%
    purrr::map(
      ~ scater::addPerCellQC(.x,
                             subsets = list(mito = mito[mito %in% rownames(.x)]))
    ) %>%
    purrr::map(
      ~ scater::addPerFeatureQC(.x)
    )
  readr::write_rds(sce_list, filename)
  return(sce_list)
}