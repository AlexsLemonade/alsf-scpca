# This script takes data from S3, imports it using the scpcaTools package, calculates QC metrics,
# and then creates two data frames with per cell and per gene metrics for all samples of interest.

# import libraries
library(magrittr)
library(scater)
library(optparse)

# import aws_copy_samples, make_sce_list, and quant_info_table functions
function_path <- file.path(".." ,"benchmarking-functions", "R")
file.path(function_path, list.files(function_path, pattern = "*.R$")) %>%
  purrr::walk(source)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_file"),
    type = "character",
    default = "sample-file.tsv",
    help = "path to tab separated file with one column with list of sample ID's used"
  ),
  make_option(
    opt_str = c("-t", "--tool_list"),
    type = "character",
    default = "tools.tsv",
    help = "path to tab separated file with one column with list of tools used"
  ),
  make_option(
    opt_str = c("-q", "--quant_s3"),
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
    help = "path to folder where all output files should be stored"
  ),
  make_option(
    opt_str = c("--save"),
    action = "store_true",
    type = "logical",
    help = "option to save individual rds files containing SingleCellExperiment objects."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check that sample file and tool list exists
if(!file.exists(opt$sample_file)){
  stop("Sample file does not exist.")
}

if(!file.exists(opt$tool_list)){
  stop("File containing tool list does not exist.")
}

# read in list of sample ID's and tools
sample_ids <- readr::read_tsv(opt$sample_file, col_names = c("sample_id")) %>%
  dplyr::pull("sample_id")
tools <- readr::read_tsv(opt$tool_list, col_names = c("tools")) %>%
  dplyr::pull("tools")

# output directory used to store aws data, sce objects and qc dataframes
if(!dir.exists(opt$output_dir)){
  dir.create(opt$output_dir, recursive = TRUE)
}

# directory within output directory to store data copied over from AWS S3
data_dir <- file.path(opt$output_dir, "data", "quants")

# create results directory to write output files
results_dir <- file.path(opt$output_dir, "results")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
}

# copy data from AWS
# grab data from S3 for each tool and sample
aws_copy_samples(local_dir = data_dir,
                 s3_dir = opt$quant_s3,
                 samples = sample_ids,
                 tools = tools)

# get library data to add cell or nucleus information for each sample
library_data_dir <- file.path(opt$output_dir, 'sample-info')
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

# Generate metadata table for all benchmarking samples containing run parameters for each sample
quant_info <- quant_info_table(data_dir = data_dir,
                               tools = tools,
                               samples = sample_ids)

# read in sample metadata and obtain sample id and sequencing unit
library_df <- readr::read_tsv(file.path(library_data_dir, "scpca-library-metadata.tsv"))
select_metadata_df <- library_df %>%
  dplyr::select(scpca_run_id, seq_unit)

# add sample metadata to quant_info
quant_info <- quant_info %>%
  dplyr::left_join(select_metadata_df, by = c("sample" = "scpca_run_id")) %>%
  # get which_counts column based on seq_unit to be used to import quants data into R
  dplyr::mutate(which_counts = dplyr::case_when(seq_unit == "cell" ~ "spliced",
                                                seq_unit == "nucleus" ~ "unspliced"))

# get mitochondrial file name
mito_file <- list.files(library_data_dir, pattern = ".mitogenes.txt$", full.names = TRUE)

# error if there is more than one mitochondrial file detected
if(length(mito_file) != 1){
  stop("Check that only 1 mitochondrial file with the pattern *.mitogenes.txt is
       in the sample-info folder in the designated output directory.")
}

# read in mito gene list
mito_genes <- readr::read_tsv(mito_file, col_names = "gene_id")
mito_genes <- mito_genes %>%
  dplyr::pull(gene_id) %>%
  unique()


# make a list of list of sce tools for each tool + run configuration
# this list will then be collapsed into two dataframes by extracting the rowData and colData from each sce
alevin_fry_tools <- c("alevin-fry", "alevin-fry-unfiltered", "alevin-fry-knee")

quant_info <- quant_info %>%
  # add group name so that different configurations of alevin fry are in separate sub lists of sces
  # doing this means that sces.rds objects are saved in smaller groups based on both tool + run parameters
  dplyr::mutate(group_name = ifelse(tool %in% alevin_fry_tools,
                       paste(index_type, alevin_alignment, alevin_resolution, sep = "_"),
                       tool)) %>%
  # put them in order for easy naming of the list
  dplyr::arrange(group_name)

# if --save is used, pass filename through with make_sce_list() to save each list of sces for each
# unique combo of tool+run parameters as its own rds file
 all_sces_list <- quant_info %>%
    # split into list of metadata tables
    dplyr::group_split(group_name, .keep = TRUE) %>%
    # make individual sce lists for each group based on tool + run configurations
    purrr::map(
      ~ make_sce_list(
        info_df = .x,
        mito = mito_genes,
        filename = ifelse(opt$save,
                          file.path(
                            results_dir,
                            paste0(.x$group_name[1], "_sces.rds")),
                          NULL)
      ))

# name list of list of sces based on tool + run configuration stored in group_name column
names(all_sces_list) <- quant_info %>%
  dplyr::pull(group_name) %>%
  unique()

# extract colData from all sces and combine into data frame with each cell as a row
coldata_df <- purrr::map_df(
  all_sces_list,
  ~ purrr::map_df(.x, scpcaTools::coldata_to_df, .id = "quant_id"),
  .id = "tool"
)

# extract rowData from all sces and combine into data frame with each gene as a row
rowdata_df <- purrr::map_df(
  all_sces_list,
  ~ purrr::map_df(.x, scpcaTools::rowdata_to_df, .id = "quant_id"),
  .id = "tool"
)


# save quant_info, coldata and rowdata qc
readr::write_tsv(quant_info, file.path(results_dir, "quant_info.tsv"))
readr::write_tsv(coldata_df, file.path(results_dir, "coldata_qc.tsv"))
readr::write_tsv(rowdata_df, file.path(results_dir, "rowdata_qc.tsv"))
