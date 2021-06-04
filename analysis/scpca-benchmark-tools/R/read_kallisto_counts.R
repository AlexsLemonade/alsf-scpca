read_kallisto_counts <- function(base) {
  counts <- Matrix::readMM(paste0(base,".mtx"))%>%
    t() %>% # transpose to gene x cell orientation
    as("dgCMatrix") # compress sparse matrix
  dimnames(counts) <- list(readLines(paste0(base,".", "genes.txt")),
                           readLines(paste0(base,".barcodes.txt")))
  return(counts)
}