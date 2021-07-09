aws_copy_samples <- function(local_dir, s3_dir, samples, tools) {
  
  ## need to add in check for tools 
  
  includes <- paste('--include "', samples, '*"', sep = '', collapse = ' ')
  
  for (tool in tools) {
    local_tool_dir = file.path(local_dir, tool)
    s3_tool_dir = glue::glue("{s3_dir}/{tool}-quant/")
    dir.create(local_tool_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (tool == 'alevin') {
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*/"', includes, '--recursive')
      system(sync_call, ignore.stdout = TRUE)
    }
    if (tool %in% c('alevin-fry', 'alevin-fry-unfiltered', 'alevin-fry-knee')) {
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*/"', includes, '--exclude "*.rad"',
                         '--recursive')
      system(sync_call, ignore.stdout = TRUE)
    }
    if (tool == 'kallisto') {
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*/"', includes, '--exclude "*/bus/*"',
                         '--recursive')
      system(sync_call, ignore.stdout = TRUE)
    }
    if (tool == 'cellranger') {
      sync_call <- paste('aws s3 cp', s3_tool_dir, local_tool_dir, 
                         '--exclude "*/"', includes, 
                         '--exclude "*/SC_RNA_COUNTER_CS/*"',
                         '--exclude "*.bam"', '--exclude "*.bam.bai"',
                         '--recursive')
      system(sync_call, ignore.stdout = TRUE)
    }
  }
}
