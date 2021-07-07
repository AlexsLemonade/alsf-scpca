perCellQCMetrics_threshold <- function(sce, mito = mito_genes, threshold){
  # add per cell QC with mitochondrial genes separately for later comparisons
  scater::perCellQCMetrics(
    sce,
    threshold = threshold
  )
}