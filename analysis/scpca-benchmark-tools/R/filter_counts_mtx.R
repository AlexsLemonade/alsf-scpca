filter_counts_mtx <- function(sce){
  # grab counts from single cell experiment
  counts <- counts(sce)
  # calculate probability of being an empty droplet
  empty_df <- DropletUtils::emptyDrops(counts)
  cells <- rownames(empty_df[which(empty_df$Limited == "TRUE" & empty_df$FDR <= 0.01),])
  # subset original counts matrix by cells that pass filter
  counts <- counts[, cells]
  return(counts)
}