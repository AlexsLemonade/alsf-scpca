
import_quant_data <- function(quant_ids, tool = tool, data_dir, intron_metadata) {
  # if tool is kallisto, then read in using the read_kallisto_counts
  if(tool == "kallisto") {
    base_file <- file.path(data_dir, quant_ids, "counts", "gene_count")
    counts_list <- base_file %>%
      purrr::map(read_kallisto_counts)
    # otherwise read in directly using tximport
    # we are only using this for alevin or kallisto based tools, so I did hardcode it for alevin in tximport
  } else {
    counts_list <- quant_ids %>%
      purrr::map(
        ~ tximport::tximport(file.path(data_dir, .x, "alevin", "quants_mat.gz"),
                   type = "alevin")) %>%
      purrr::map(~ SingleCellExperiment::SingleCellExperiment(list(counts = .x$counts))) %>%
      purrr::map(counts)
  }
  # once we have a list of counts matrices we want to take any that were aligned to the spliced_intron index and collapse the counts
  names(counts_list) <- quant_ids
  counts_list_nuclei <- counts_list[grep("spliced_intron", names(counts_list))] %>%
    purrr::map(
      ~ collapse_intron_counts(.x, intron_metadata)
    )
  
  # replace the sces in the original list with the collapsed sces
  counts_list[names(counts_list_nuclei)] <- counts_list_nuclei
  
  # convert counts list to sces and return 
  sces <- counts_list %>%
    purrr::map(~ SingleCellExperiment::SingleCellExperiment(list(counts = .x)))
  return(sces)
  
}