#!/usr/bin/env Rscript

# This script takes the input reference fasta of the genome and the 
# corresponding gtf file, identifies the regions of interest corresponding
# to spliced and unspliced transcripts, then subsets the genome and gtf for 
# only those particular regions, and finally outputs an expanded.fa, 
# expanded.gtf, and expanded.tx2gene.tsv all corresponding to the genomic 
# regions with spliced and unspliced transcripts.

# Additionally this script creates a transcript2gene mapping file 
# for use in salmon alevin, and a output a list of mitochondrial 
# genes for use in QC. 

# load needed packages
library(Biostrings)
library(BSgenome)
library(eisaR)
library(GenomicFeatures)
library(magrittr)
library(optparse)
library(tidyverse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-f", "--gtf"),
    type = "character",
    default = "annotation/Homo_sapiens.GRCh38.103.gtf.gz",
    help = "File path for input gtf file",
  ),
  make_option(
    opt_str = c("-g", "--genome"),
    type = "character",
    default = "fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    help = "File path for reference fasta file",
  ), 
  make_option(
    opt_str = c("-O", "--output_dir"),
    type = "character",
    default = "annotation",
    help = "Directory to write output files",
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Make output paths 
expanded_fasta <- file.path(opt$output_dir, 
                            "Homo_sapiens.ensembl.103.expanded.fa")
expanded_gtf <- file.path(opt$output_dir, 
                          "Homo_sapiens.ensembl.103.expanded.gtf")
expanded_tx2gene <- file.path(opt$output_dir, 
                              "Homo_sapiens.ensembl.103.tx2gene_expanded.tsv")
tx2_gene_out <- file.path(opt$output_dir, 
                      "Homo_sapiens.ensembl.103.tx2gene.tsv")
mito_out <- file.path(opt$output_dir,
                      "Homo_sapiens.ensembl.103.mitogenes.txt")

# Check for output directory 
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir)
}

# extract GRanges object containing genomic coordinates of each 
# annotated transcript - both spliced and unspliced transcripts
grl <- eisaR::getFeatureRanges(
  gtf = opt$gtf,
  featureType = c("spliced", "unspliced"), 
  verbose = TRUE
)

# load in primary genome
genome <- Biostrings::readDNAStringSet(opt$genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)

# extract unspliced and spliced sequences
# genomic regions defined above
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)

# write unspliced and spliced sequences to fasta file
Biostrings::writeXStringSet(
  seqs, filepath = expanded_fasta
)

# write the associated annotations to gtf 
eisaR::exportToGtf(
  grl, 
  filepath = expanded_gtf
)

# create text file mapping transcript and intron identifiers to 
# corresponding gene identifiers
# first make Tx2Gene for all spliced and unspliced sequences
full_tx2gene <- eisaR::getTx2Gene(
  grl, filepath = expanded_tx2gene
)

# next make Tx2gene for only spliced transcripts 

# get list of all spliced transcript ID's
spliced_genes = metadata(grl)$corrtx %>%
  pull(spliced)

# subset sequences for only 
splice_grl = grl[spliced_genes]

# export Tx2Gene for spliced transcripts
spliced_tx2gene <- eisaR::getTx2Gene(
  splice_grl, filepath = tx2_gene_out
)

# reimport gtf to get list of mito genes 
gtf <- rtracklayer::import(expanded_gtf)
mitogenes <- gtf[seqnames(gtf) == 'MT']

# write out mitochondrial gene list
writeLines(mitogenes$gene_id, mito_out)