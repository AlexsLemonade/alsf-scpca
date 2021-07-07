# This script takes data from S3, imports it using the scpcaTools package, calculates QC metrics, 
# and then creates two data frames with per cell and per gene metrics for all samples of interest. 

# import libraries 
library(magrittr)
library(scater)
library(optparse)

# import benchmarking functions
function_path <- file.path(".." ,"benchmarking-functions", "R")
miceadds::source.all(function_path)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_ids"),
    type = "character",
    help = "list of scpca_run_ids used for benchmarking"
  ),
  make_option(
    opt_str = c("-t", "--tools"),
    type = "character",
    help = "tools used for benchmarking"
  ), 
  make_option(
    opt_str = c("-s", "--quant_s3"),
    type = "character",
    default = "s3://nextflow-ccdl-results/scpca",
    help = "storage location for quant output files on S3"
  ), 
  make_option(
    opt_str = c("-a", "--annotation_files_s3"),
    type = "character",
    default = "s3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103/annotation/",
    help = "storage location for annotation info data on S3, must have *.mitogenes.txt"
  ), 
  make_option(
    opt_str = c("-o", "--output_dir"),
    type = "character",
    default = "results",
    help = "path to folder where output files should be stored"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# file paths for data
base_dir <- here::here()
data_dir <- file.path(base_dir, 'data', 'quants')

# output directory used to store sce objects and qc dataframes 
if(!dir.exists(opt$output_dir)){
  dir.create(opt$output_dir, recursive = TRUE)
}

# copy data from AWS
aws_copy_samples(local_dir = data_dir,
                 s3_dir = opt$quant_s3,
                 samples = opt$sample_ids,
                 tools = opt$tools)

# get library data to add cell or nucleus information for each sample
library_data_dir <- file.path(base_dir, 'sample-info')
dir.create(library_data_dir, recursive = TRUE, showWarnings = FALSE)
library_data_s3 <- 's3://ccdl-scpca-data/sample_info'

# grab library metadata from location in s3
sync_call <- paste('aws s3 sync', library_data_s3, library_data_dir,
                   '--exclude "*"', 
                   '--include "*scpca-library-metadata.tsv"')
system(sync_call, ignore.stdout = TRUE)

# get mitochondrial gene list from s3
sync_call <- paste('aws s3 sync', opt$annotation_files_s3, library_data_dir,
                   '--exclude "*"', 
                   '--include "*.mitogenes.txt"')
system(sync_call, ignore.stdout = TRUE)

# Generate quant_info table for all benchmarking samples
quant_info <- quant_info_table(data_dir = data_dir,
                               tools = opt$tools,
                               samples = opt$sample_ids)

# read in sample metadata
library_df <- readr::read_tsv(file.path(library_data_dir, "scpca-library-metadata.tsv"))
select_metadata_df <- library_df %>%
  dplyr::select(scpca_run_id, seq_unit)

# add sample metadata to quant_info
quant_info <- quant_info %>% 
  dplyr::left_join(select_metadata_df, by = c("sample" = "scpca_run_id")) %>%
  dplyr::mutate(which_counts = dplyr::case_when(seq_unit == "cell" ~ "spliced",
                                                seq_unit == "nucleus" ~ "unspliced"))

# read in mito gene list
mito_file <- file.path(library_data_dir, "Homo_sapiens.GRCh38.103.mitogenes.txt")
mito_genes <- readr::read_tsv(mito_file, col_names = "gene_id")
mito_genes <- mito_genes %>%
  dplyr::pull(gene_id) %>%
  unique()

# get sces for each tool

# load in data for non-fry tools first
# process alevin-fry tools separately since there are multiple different configurations possible 

sces_list <- list()
non_fry_tools = opt$tools[which(opt$tools %in% c("cellranger", "kallisto", "alevin"))]
if(length(non_fry_tools) > 0) {
  for(i in length(non_fry_tools)){
    tool_info <- quant_info %>%
      dplyr::filter(tool %in% non_fry_tools[i])
    
    sces_list[[non_fry_tools[[i]]]] <- make_sce_list(tool_info,
                                                     filename = file.path(opt$output_dir,
                                                                          paste(non_fry_tools[[i]], "sces.rds", sep = "_")),
                                                     mito = mito_genes)
  }
}

# make coldata and rowdata dataframes for non-fry tools
i = 1
coldata_non_fry_df_list <- list()
for(list in sces_list){
  coldata_non_fry_df_list[[non_fry_tools[[i]]]] <- list %>%
    purrr::map_df(coldata_to_df, .id = "quant_id")
  i = i + 1
}

coldata_non_fry_df <- coldata_non_fry_df_list %>%
  dplyr::bind_rows(.id = "tool")


i = 1
rowdata_non_fry_df_list <- list()
for(list in sces_list){
  rowdata_non_fry_df_list[[non_fry_tools[[i]]]] <- list %>%
    purrr::map_df(rowdata_to_df, .id = "quant_id")
  i = i + 1
}

rowdata_non_fry_df <- rowdata_non_fry_df_list %>%
  dplyr::bind_rows(.id = "tool")

# process alevin fry data

if("alevin-fry-unfiltered" %in% opt$tools) {

  # make separate tables to save each condition of alevin-fry as its own sce 
  # making smaller sces makes it easier to work with if we ever want to go back to these later 
  
  alevin_info <- quant_info %>%
    dplyr::filter(tool == "alevin-fry") %>%
    dplyr::mutate(group_name = paste(index_type, alevin_alignment, alevin_resolution, sep = "_")) %>%
    dplyr::arrange(group_name) 
  
  # make a list of split data frames based on different configurations
  alevin_split_info <- alevin_info %>%
    dplyr::group_split(group_name)
  
  names(alevin_split_info) <- alevin_info %>%
    dplyr::pull(group_name) %>%
    unique()
  
  # get sces
  alevin_fry_sces <- list()
  i = 1
  for(df in alevin_split_info){
    alevin_fry_sces[[names(alevin_split_info)[i]]] <- make_sce_list(df,
                                                                    filename = file.path(opt$output_dir,
                                                                                         paste(names(alevin_split_info)[i], "sces.rds", sep = "_")),
                                                                    mito = mito_genes)
    i = i+1
  }
  
  # make coldata qc df
  i = 1
  coldata_fry_unfiltered_qc_list <- list()
  for(sce_list in alevin_fry_sces){
    coldata_fry_unfiltered_qc_list[[i]] <- sce_list %>%
      purrr::map_df(coldata_to_df, .id = "quant_id")
    i = i + 1
  }
  coldata_fry_unfiltered_df <- coldata_fry_unfiltered_qc_list %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(tool = "alevin-fry-unfiltered")
  
  # make rowdata qc df
  
  i = 1
  rowdata_fry_unfiltered_qc_list <- list()
  for(sce_list in alevin_fry_sces){
    rowdata_fry_unfiltered_qc_list[[i]] <- sce_list %>%
      purrr::map_df(rowdata_to_df, .id = "quant_id")
    i = i + 1
  }
  rowdata_fry_unfiltered_df <- rowdata_fry_unfiltered_qc_list %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(tool = "alevin-fry-unfiltered")
  

}


# if alevin-fry-unfiltered and non-alevin tools, join together 

if("alevin-fry-unfiltered" %in% opt$tools & length(non_fry_tools) > 0){
  coldata_qc <- dplyr::bind_rows(
    list(coldata_fry_unfiltered_df,
         coldata_non_fry_df)
  )
  rowdata_qc <- dplyr::bind_rows(
    list(rowdata_fry_unfiltered_df,
         rowdata_non_fry_df)
  )
} else if ("alevin-fry-unfiltered" %in% opt$tools & length(non_fry_tools) == 0) {
  coldata_qc <- coldata_fry_unfiltered_df
  rowdata_qc <- rowdata_non_fry_df
} else {
  coldata_qc <- coldata_non_fry_df
  rowdata_qc <- rowdata_non_fry_df
}

# save quant_info, coldata and rowdata qc
readr::write_tsv(quant_info, file.path(opt$output_dir, "quant_info.tsv"))
readr::write_tsv(coldata_qc, file.path(opt$output_dir, "coldata_qc.tsv"))
readr::write_tsv(rowdata_qc, file.path(opt$output_dir, "rowdata_qc.tsv"))


