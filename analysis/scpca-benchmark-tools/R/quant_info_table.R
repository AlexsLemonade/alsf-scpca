quant_info_table <- function(data_dir, tools, samples){
  
  for (tool in tools){
    
    tool_data_dir <- file.path(data_dir, tool)
    
    if (tool == "cellranger") {
      cellranger_quant_info <- data.frame(tool = "cellranger", 
                                          quant_dir = list.dirs(tool_data_dir, 
                                                                recursive = FALSE, 
                                                                full.names = FALSE)) %>%
        tidyr::separate(quant_dir, sep = "[-]",
                        into = c("sample", "index_type"),
                        remove = FALSE) %>%
        dplyr::filter(sample %in% samples) %>%
        dplyr::mutate(index_type = ifelse(index_type == "cdnapre_mRNA", "splici", "cDNA"))
    } if (tool %in% c("alevin-fry-unfiltered", "alevin-fry-knee")) {
      alevin_fry_quant_info <- data.frame(tool = "alevin-fry", 
                                          quant_dir = list.dirs(tool_data_dir, 
                                                                recursive = FALSE, 
                                                                full.names = FALSE)) %>%
        tidyr::separate(quant_dir, sep = "[-]",
                        into = c("sample", "index_type", "alevin_alignment", "alevin_resolution"),
                        extra = "drop", 
                        fill = "right") %>%
        dplyr::filter(sample %in% samples) %>%
        dplyr::mutate(index_type = ifelse(index_type == "spliced_intron_txome_k31", "splici", "cDNA"))
    } %>% 
  }
}