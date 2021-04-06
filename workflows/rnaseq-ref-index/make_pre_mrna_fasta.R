#!/usr/bin/env Rscript

# This script takes the input reference fasta of the genome and the 
# corresponding gtf file, identifies the regions of interest corresponding
# to spliced transcripts and unspliced transcripts, then subsets the genome 
# and gtf for only those particular regions, and finally outputs an expanded.fa, 
# expanded.gtf, and expanded.tx2gene.tsv all corresponding to the genomic 
# regions with spliced and unspliced transcripts. 

# load needed packages
library(magrittr)
library(optparse)
library(eisaR)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome)

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
    default = "expanded_annotation",
    help = "Directory to write output files",
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Make output paths 
expanded_fasta <- file.path(opt$output_dir, 
                            "ensembl.GRCh38.annotation.expanded.fa")
expanded_gtf <- file.path(opt$output_dir, 
                          "ensembl.GRCh38.annotation.expanded.gtf")
expanded_tx2gene <- file.path(opt$output_dir, 
                              "ensembl.GRCh38.annotation.expanded.tx2gene.tsv")

# Check for output directory 
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir)
}

# extract GRanges object containing genomic coordinates of each 
# annotated transcript and intron 
gtf <- opt$gtf
grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 90L, 
  joinOverlappingIntrons = FALSE,
  verbose = TRUE
)
#grl[4:6]

# load in primary genome
genome <- Biostrings::readDNAStringSet(opt$genome)
head(genome)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)

# extract sequences for features of interest based on spliced and intronic
# genomic regions defined above
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)

# write regions of interest to fasta file
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
df <- eisaR::getTx2Gene(
  grl, filepath = expanded_tx2gene
)