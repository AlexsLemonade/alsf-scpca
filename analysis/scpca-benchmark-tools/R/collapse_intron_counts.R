# for every combination of mrna + intronic region create a collapsed counts matrix that contains the sum of the total counts for that gene (mrna + intron)
# we could use tximeta::splitSE but it just creates two separate matrices rather than one combined one 
collapse_intron_counts <- function(counts, intron_metadata) {
  # find the shared genes in that counts matrix 
  shared_genes <- intersect(row.names(counts), rownames(intron_metadata))
  # replace row names with -I appended with corresponding spliced gene 
  row.names(counts)[which(row.names(counts) %in% shared_genes)] <- intron_metadata[shared_genes, "spliced"]
  # aggregate Matrix counts by gene name 
  Matrix.utils::aggregate.Matrix(counts, row.names(counts))
}