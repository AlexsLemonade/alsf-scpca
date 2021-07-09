quant_info_table <- function(data_dir, tools, samples){

  quant_info_list <- list()

  for (tool in tools){

    tool_data_dir <- file.path(data_dir, tool)

    if (tool %in% c("cellranger", "alevin", "kallisto")) {

      tool_quant_info <- data.frame(tool = tool,
                                          quant_dir = list.dirs(tool_data_dir,
                                                                recursive = FALSE,
                                                                full.names = FALSE)) %>%
        tidyr::separate(quant_dir, sep = "[-]",
                        into = c("sample", "index_type"),
                        remove = FALSE) %>%
        dplyr::filter(sample %in% samples) %>%
        dplyr::mutate(index_type = dplyr::case_when(index_type %in% c("cdnapre_mRNA", "spliced_intron_txome_k31_full_sa",
                                                               "spliced_intron_txome_k31") ~ 'splici',
                                             index_type %in% c("cdna", "spliced_txome_k31_full_sa",
                                                               "txome_k31_full_sa", "cdna_k31_full_sa",
                                                               "cdna_k31", "txome_k31") ~ 'cDNA'),
                      data_dir = file.path(tool_data_dir, quant_dir),
                      filter = ifelse(tool == "kallisto", TRUE, FALSE))


      quant_info_list[[tool]] <- tool_quant_info

    } else if (tool %in% c("alevin-fry-unfiltered", "alevin-fry-knee")) {

      alevin_fry_quant_info <- data.frame(tool = tool,
                                          quant_dir = list.dirs(tool_data_dir,
                                                                recursive = FALSE,
                                                                full.names = FALSE)) %>%
        tidyr::separate(quant_dir, sep = "[-]",
                        into = c("sample", "index_type", "alevin_alignment", "alevin_resolution"),
                        extra = "drop",
                        fill = "right",
                        remove = FALSE) %>%
        dplyr::mutate(filter_strategy = dplyr::case_when(tool == "alevin-fry-unfiltered" ~ "unfiltered",
                                                         tool == "alevin-fry-knee" ~ "knee"),
                      tool = stringr::str_extract(tool, "alevin-fry")) %>%
        dplyr::filter(sample %in% samples) %>%
        dplyr::mutate(index_type = ifelse(index_type == "spliced_intron_txome_k31", "splici", "cDNA"),
                      data_dir = file.path(tool_data_dir, quant_dir),
                      filter = ifelse(alevin_alignment == "knee", FALSE, TRUE))

      quant_info_list[[tool]] <- alevin_fry_quant_info
    }
  }

  all_quant_info <- dplyr::bind_rows(quant_info_list) %>%
    dplyr::mutate(alevin_alignment = tidyr::replace_na(alevin_alignment, "not_alevin"),
                  alevin_resolution = tidyr::replace_na(alevin_resolution, "not_alevin"),
                  usa_mode = index_type == "splici" & alevin_resolution == "cr" &  
                  intron_mode = index_type == "splici" & tool != "cellranger"  & !usa_mode
  return(all_quant_info)
}
