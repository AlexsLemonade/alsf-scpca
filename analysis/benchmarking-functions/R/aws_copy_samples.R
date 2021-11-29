#' Copy samples from AWS S3 to a local directory
#'
#' @param local_dir Path to local directory where quants files will be copied to.
#' @param s3_dir S3 bucket where quants files are stored. 
#' @param samples List of samples to include when copying over quants files.
#' @param tools List of tools used for producing quants files, 
#'   dictating which files will be included in copying.
#'
#'
#' @examples
#' \dontrun{
#' # copy over files from alevin-fry-unfiltered directory on S3 to local output folder
#' aws_copy_samples( local_dir = outputs_directory,
#'                   s3_dir = "s3://nextflow-ccdl-results/scpca",
#'                   samples = "SCPCR000006",
#'                   tool = "alevin-fry-unfiltered"
#'
#' )
aws_copy_samples <- function(local_dir, s3_dir, samples, tools) {
  
  tool_names <- c("alevin", "alevin-fry", "alevin-fry-unfiltered", 
                  "alevin-fry-knee", "cellranger", "kallisto")
  
  if(length(samples) == 0){
    stop("No sample ids input.")
  }
  
  if(length(tools) == 0){
   stop("Missing input for tools used.") 
  }
  
  # make one character string with --include {sample} for all samples to be copied 
  includes <- paste("--include '", samples, "*'", sep = '', collapse = ' ')

  for (tool in tools) {
    
    if(!(tool %in% tool_names)){
      stop("Invalid tool name.")
    }
    
    # create a directory for the tool that is being used if not already existent
    local_tool_dir = file.path(local_dir, tool)
    if(!dir.exists(local_tool_dir)){
      dir.create(local_tool_dir, recursive = TRUE, showWarnings = FALSE) 
    }
    s3_tool_dir = glue::glue("{s3_dir}/{tool}-quant/")
    
    # based on each tool create the aws s3 cp call used to copy files from S3 to the local directory 
    # for each tool exclude certain large files 
    if (tool == 'alevin') {
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*"', includes, '--recursive', sep = " ")
    }
    if (tool %in% c('alevin-fry', 'alevin-fry-unfiltered', 'alevin-fry-knee')) {
      # exclude large rad files 
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*"', includes, '--exclude "*.rad"',
                         '--recursive')
    }
    if (tool == 'kallisto') {
      # excludes large outputs in bus directory for each sample
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*"', includes, '--exclude "*/bus/*"',
                         '--recursive')
    }
    if (tool == 'cellranger') {
      # exclude large bam files 
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*"', includes, 
                         '--exclude "*/SC_RNA_COUNTER_CS/*"',
                         '--exclude "*/SPATIAL_RNA_COUNTER_CS/*"',
                         '--exclude "*.bam"', '--exclude "*.bam.bai"',
                         '--recursive')
    }
    system(sync_call, ignore.stdout = TRUE)
  }
}
