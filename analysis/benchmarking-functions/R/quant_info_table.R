#' Make metadata table with information about the run for each sample and tool combination. 
#'
#' @param data_dir Path to the local directory where quants output files are stored.
#' @param tools Tools used to produce the quants output for each sample. 
#'   Each tool should have a subdirectory within the data_dir and 
#'   each sample should have a subdirectory within each tool.
#' @param samples List of samples to keep in the final metadata table.
#'
#' @return Metadata table with each row corresponding to a sample and information 
#'   about the quants file and run parameters used for that sample and tool combination.
#'   The columns in the table include: 
#'   data_dir - Local directory where the output quants data for all samples is stored. 
#'   tool - which tool was used to process the samples. 
#'   which_counts - Which counts to keep in the counts matrix, cDNA only ("spliced") or cDNA and introns ("unspliced").
#'   intron_mode - Whether or not the data was aligned to the cDNA only (FALSE) or cDNA and introns index (TRUE).
#'   usa_mode - If the data was processed with alevin-fry, was it processed using USA mode. 
#'   filter - Logical indicating if the returned SingleCellExperiment should contain filtered (TRUE) or unfiltered (FALSE) data.
#'
#' @examples
#' \dontrun{
#' # make a table for alevin-fry-unfiltered
#' quant_info_table(data_dir = data_dir,
#'                  tools = "alevin-fry-unfiltered",
#'                  samples = c("SCPCR000006", "SCPCR000220")
#')
#' 
#' }
quant_info_table <- function(data_dir, tools, samples){

  # intialize a list to be populated with a smaller metadata table for each tool 
  quant_info_list <- list()
  
  # check that the data directory is not empty 
  if(length(list.files(data_dir)) == 0){
    stop("Data directory is empty, no files present.")
  }
  
  # check that tools are input
  if(length(tools) == 0){
    stop("Missing input for tools used.") 
  }
  
  # check for samples input 
  if(length(samples) == 0){
    stop("No sample ids input.")
  }
  
  # make a list of tool names to check against 
  tool_names <- c("alevin", 
                  "alevin-fry", 
                  "alevin-fry-unfiltered", 
                  "alevin-fry-knee", 
                  "cellranger", 
                  "kallisto")
  
  # group possible index types into splici or cDNA  
  splici_types <- c("cdnapre_mRNA", 
                    "spliced_intron_txome_k31_full_sa",
                    "spliced_intron_txome_k31",
                    "cdna-pre_mRNA")
  
  cdna_types <- c("cdna", 
                  "spliced_txome_k31_full_sa",
                  "txome_k31_full_sa", 
                  "cdna_k31_full_sa",
                  "cdna_k31", 
                  "txome_k31", 
                  "spliced_txome_k31")
  
  # make a list of tools that use knee filtering
  knee_tools <- c("cellranger",
                  "kallisto")

  for (tool in tools){
    
    # check to make sure tool is a valid tool option before proceeding
    if(!(tool %in% tool_names)){
      stop("Invalid tool name.")
    }
    
    tool_data_dir <- file.path(data_dir, tool)
    
    # initiliaze the metadata table for each tool using the list of directories within the tool directory
    # the list of directories within each tool's directory should correspond to the sample names + run configurations
    tool_quant_info <- data.frame(tool = tool,
                                  quant_dir = list.dirs(tool_data_dir,
                                                        recursive = FALSE,
                                                        full.names = FALSE)) 
    
    # grab the run parameters from the directory names (quant_dir) and group parameters into their respective columns
    # this will be dependent on the tool that was used 
    
    if (tool %in% c("cellranger", "alevin", "kallisto")) {
      
      # split directory name by [-] to get sample name and index type 
      tool_quant_info <- tool_quant_info %>%
        tidyr::separate(quant_dir, sep = "[-]",
                        extra = "merge",
                        into = c("sample", "index_type"),
                        remove = FALSE) %>%
        dplyr::mutate(data_dir = file.path(tool_data_dir, quant_dir))


    } else if (tool %in% c("alevin-fry-unfiltered", "alevin-fry-knee")) {
      
      # split directory name by [-] to get sample name, index type, alignment strategy, and resolution used
      tool_quant_info <- tool_quant_info %>%
        tidyr::separate(quant_dir, sep = "[-]",
                        into = c("sample", "index_type", "alevin_alignment", "alevin_resolution"),
                        extra = "drop",
                        fill = "right",
                        remove = FALSE) %>%
        # create filter strategy column based on tool used for alevin fry
        dplyr::mutate(filter_strategy = dplyr::case_when(tool == "alevin-fry-unfiltered" ~ "unfiltered",
                                                         tool == "alevin-fry-knee" ~ "knee"),
                      tool = stringr::str_extract(tool, "alevin-fry"), 
                      data_dir = file.path(tool_data_dir, quant_dir))

    }
    # add each smaller metadata to named list of metadata tables 
    quant_info_list[[tool]] <- tool_quant_info
  }
# merge into one large metadata table for all samples and all tool configurations used 
  all_quant_info <- dplyr::bind_rows(quant_info_list) %>%
    # filter for only specific samples that are of interest
    dplyr::filter(sample %in% samples) %>%
    # convert index type to either splici or cDNA
    dplyr::mutate(index_type = dplyr::case_when(index_type %in% splici_types ~ 'splici',
                                                index_type %in% cdna_types ~ 'cDNA'),
                  # create columns with information about each parameter used for each run 
                  filter = ifelse(alevin_alignment == "knee" | tool %in% knee_tools, FALSE, TRUE), 
                  alevin_alignment = tidyr::replace_na(alevin_alignment, "not_alevin"),
                  alevin_resolution = tidyr::replace_na(alevin_resolution, "not_alevin"),
                  usa_mode = index_type == "splici" & alevin_resolution == "cr",
                  intron_mode = index_type == "splici" & tool != "cellranger"  & !usa_mode)
  return(all_quant_info)
}
