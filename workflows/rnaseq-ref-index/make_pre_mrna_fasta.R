#!/usr/bin/env Rscript

# This script takes the input reference fasta of the genome and the 
# corresponding gtf file, identifies the regions of interest corresponding
# to spliced transcripts and introns, then subsets the genome and gtf for 
# only those particular regions, and finally outputs an expanded.fa, 
# expanded.gtf, and expanded.tx2gene.tsv all corresponding to the genomic 
# regions with spliced transcripts and introns

# This script then subsets for only spliced transcripts to 
# create a corresponding fasta, gtf, and transcript2gene mapping file 
# for use with index generation for scRNA-seq 

# finally, a list of mitochondrial genes is output for QC

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
    opt_str = c("-a", "--annotation_output"),
    type = "character",
    default = "annotation",
    help = "Directory to write output files",
  ), 
  make_option(
    opt_str = c("-o", "--fasta_output"),
    type = "character",
    default = "fasta",
    help = "Directory to write output files",
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Get base file name from input gtf file
file_base <- stringr::str_replace(basename(opt$gtf), "\\.gtf(\\.gz)?$", "")

# main files with only spliced regions
txome_fasta_file <- paste0(file_base, ".spliced.txome.fa.gz")
txome_gtf_file <- paste0(file_base, ".spliced.txome.gtf")
txome_tx2gene_file <- paste0(file_base, ".spliced.tx2gene.tsv")

## files with spliced + intron regions
spliced_intron_fasta_file <- paste0(file_base,".", "spliced_intron.txome.fa.gz")
spliced_intron_gtf_file <- paste0(file_base, ".", "spliced_intron.txome.gtf")
spliced_intron_tx2gene_file <- paste0(file_base, ".", 
                                      "spliced_intron.tx2gene.tsv")
spliced_intron_metadata_file <- paste0(file_base, ".", 
                                       "spliced_intron.metadata.tsv")

mito_file <- paste0(file_base, ".", "mitogenes.txt")

# make final output file names needed
txome_fasta <- file.path(opt$fasta_output, txome_fasta_file)
txome_gtf <- file.path(opt$annotation_output, txome_gtf_file)
txome_tx2gene_file <- file.path(opt$annotation_output, txome_tx2gene_file)

spliced_intron_fasta <- file.path(opt$fasta_output, spliced_intron_fasta_file)
spliced_intron_gtf <- file.path(opt$annotation_output, spliced_intron_gtf_file)
spliced_intron_tx2gene <- file.path(opt$annotation_output, 
                                    spliced_intron_tx2gene_file)
spliced_intron_metadata <- file.path(opt$annotation_output, 
                                     spliced_intron_metadata_file)

mito_out <- file.path(opt$annotation_output, mito_file)

# Check for output directory 
if (!dir.exists(opt$fasta_output)) {
  dir.create(opt$fasta_output)
}

if (!dir.exists(opt$annotation_output)) {
  dir.create(opt$annotation_output)
}


# extract GRanges object containing genomic coordinates of each 
# annotated transcript - both spliced transcripts and intronic regions
grl <- eisaR::getFeatureRanges(
  gtf = opt$gtf,
  featureType = c("spliced", "intron"), 
  verbose = TRUE
)

# load in primary genome
genome <- Biostrings::readDNAStringSet(opt$genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)

# extract spliced and intron sequences
# genomic regions defined above
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)

# write spliced and intron sequences to fasta file
Biostrings::writeXStringSet(
  seqs, filepath = spliced_intron_fasta, compress = TRUE
)

# write the associated annotations to gtf 
eisaR::exportToGtf(
  grl, 
  filepath = spliced_intron_gtf
)

# create text file mapping transcript and intron identifiers to 
# corresponding gene identifiers
# first make Tx2Gene for all spliced and intron sequences
full_tx2gene <- eisaR::getTx2Gene(
  grl, filepath = spliced_intron_tx2gene
)

# next make Tx2gene for only spliced transcripts 
## need to write out to metadata 
readr::write_tsv(metadata(grl)$corrgene, spliced_intron_metadata)

# get list of all spliced transcript ID's
spliced_genes = metadata(grl)$featurelist$spliced

# subset sequences for only spliced cDNA
splice_grl = grl[spliced_genes]

splice_seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = splice_grl
)

# write spliced sequences only to fasta file
Biostrings::writeXStringSet(
  splice_seqs, filepath = txome_fasta, compress = TRUE
)

# write the associated annotations to gtf 
eisaR::exportToGtf(
  splice_grl, 
  filepath = txome_gtf
)

# export Tx2Gene for spliced transcripts
spliced_tx2gene <- eisaR::getTx2Gene(
  splice_grl, filepath = txome_tx2gene_file
)

# reimport gtf to get list of mito genes 
gtf <- rtracklayer::import(txome_gtf)
mitogenes <- gtf[seqnames(gtf) == 'MT']

# write out mitochondrial gene list
writeLines(mitogenes$gene_id, mito_file)
