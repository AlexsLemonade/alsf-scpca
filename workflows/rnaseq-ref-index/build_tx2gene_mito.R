#!/usr/bin/env Rscript

# This script uses Bioconductor tools to download ensembl annotation data,
# create a transcript2gene mapping file for use in salmon alevin,
# and a output a list of mitochondrial genes for use in QC.

# Currently hard coded to Ensembl v100 Homo sapiens, mostly because 
# AnnotationHub is not as simple to search efficiently as I might like.

library(AnnotationHub)
library(magrittr)
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-t", "--t2g_out"),
    type = "character",
    default = "annotation/Homo_sapiens.ensembl.100.tx2gene.tsv",
    help = "File path for tx2gene output tsv file.",
  ),
  make_option(
    opt_str = c("-m", "--mito_out"),
    type = "character",
    default = "annotation/Homo_sapiens.ensembl.100.mitogenes.txt",
    help = "File path mitochondrial genes output.",
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

hub = AnnotationHub()
# Ensembl v100 Homo Sapiens is AH79689
ensdb = hub[["AH79689"]]

# get tx to gene table
tx <- transcriptsBy(ensdb, "gene")
# use the versioned transcript ids
as.data.frame(tx) %>%
  dplyr::select(tx_id_version, gene_id) %>%
  readr::write_tsv(opt$t2g_out,
                   col_names = FALSE)

# get genes and select mitochondrial only
ensg <- genes(ensdb)
mitogenes <- ensg[seqnames(ensg) == 'MT']

# write out mitochondrial gene list
writeLines(mitogenes$gene_id, opt$mito_out)

