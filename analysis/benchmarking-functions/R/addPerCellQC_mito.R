addPerCellQC_mito <- function(sce, mito = mito_genes){
  # add per cell QC with mitochondrial genes separately for later comparisons
  scater::addPerCellQC(
    sce, 
    subsets = list(mito = mito[mito %in% rownames(sce)])
  )
}