#' Import Spaceranger output for Spatial Transcriptomics library into a SpatialExperiment object.
#'
#' @param spaceranger_outs_dir Path to "outs" folder output by Spaceranger.
#' @param sample_name Sample id.
#'
#' @return SpatialExperiment object.
create_spaceranger_spe <- function(spaceranger_outs_dir, sample_name){
  # read in spaceranger outputs directly using read10XVisium
  spaceranger_spe <- SpatialExperiment::read10xVisium(spaceranger_outs_dir,
                                                      sample_id = sample_name,
                                                      type = "sparse",
                                                      data = "raw",
                                                      images = "lowres", 
                                                      load = FALSE)
}
