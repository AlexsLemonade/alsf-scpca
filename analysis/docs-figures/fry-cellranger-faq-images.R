## This is the script used to make the figures that are found in scpca-docs comparing Alevin-fry to Cell Ranger 
## The sces used as input were generated during benchmarking using 
## quantifier-comparisons/benchmarking_generate_qc_df.R and saved in data/results

library(magrittr)
library(ggplot2)
library(SingleCellExperiment)
library(patchwork)
library(ggpubr)

# path to results files with sces and qc dataframes from running benchmarking_generate_qc_df.R
base_dir <- here::here()
sces_dir <- file.path(base_dir, "data", "results")

# sce objects
cr_like_file <- file.path(sces_dir, "splici_salign_cr_sces.rds")
cellranger_file <- file.path(sces_dir, "cellranger_sces.rds")

# plots to save
umi_density <- file.path("plots", "total_umi_per_cell.png")
gene_density <- file.path("plots", "total_genes_per_cell.png")
gene_exp_correlation <- file.path("plots", "gene_exp_correlation.png")

# read in mito gene list
library_data_dir <- file.path(base_dir, 'sample-info')
mito_file <- file.path(library_data_dir, "Homo_sapiens.GRCh38.103.mitogenes.txt")
mito_genes <- readr::read_tsv(mito_file, col_names = "gene_id")
mito_genes <- mito_genes %>%
  dplyr::pull(gene_id) %>%
  unique()

# grab library metadata from AWS 
library_data_s3 <- 's3://ccdl-scpca-data/sample_info'

# grab library metadata from location in s3
sync_call <- paste('aws s3 sync', library_data_s3, library_data_dir,
                   '--exclude "*"',
                   '--include "*scpca-library-metadata.tsv"')
system(sync_call, ignore.stdout = TRUE)

# read in library metadata to grab SCPCL id later 
library_df <- readr::read_tsv(file.path(library_data_dir, "scpca-library-metadata.tsv"))
select_metadata_df <- library_df %>%
  dplyr::select(scpca_run_id, scpca_library_id)

# read in cr_like_sce
cr_like_sce <- readr::read_rds(cr_like_file) 

# remove 006, 118, and 119
cr_like_sce <- cr_like_sce[!(names(cr_like_sce) %in% c("SCPCR000006-spliced_intron_txome_k31-salign-cr-like", 
                                                                "SCPCR000118-spliced_intron_txome_k31-salign-cr-like",
                                                                "SCPCR000119-spliced_intron_txome_k31-salign-cr-like"))]

# read in cellranger 
cellranger_sce <- readr::read_rds(cellranger_file)

# make combined list of all sces 
all_sces_list <- list(cr_like_sce,cellranger_sce)
names(all_sces_list) <- c("Alevin-fry", "Cell Ranger")

# make list of samples to use for assigning samples to cell or nuclei 
single_cell <- c("SCPCR000126", "SCPCR000127")
single_nuclei <- c("SCPCR000220", "SCPCR000221")

# grab coldata from all samples and create a dataframe used for plotting
coldata_df <- purrr::map_df(
  all_sces_list,
  ~ purrr::map_df(.x, scpcaTools::coldata_to_df, .id = "quant_id"),
  .id = "tool"
) %>% 
  # create new columns with just sample ID and then seq unit 
  dplyr::mutate(sample = stringr::word(quant_id, 1, sep = "-"),
                seq_unit = dplyr::case_when(sample %in% single_cell ~ "Cell",
                                            sample %in% single_nuclei ~ "Nuclei")) %>%
  # join with metadata to get corresponding library ID 
  dplyr::left_join(select_metadata_df, by = c("sample" = "scpca_run_id")) %>%
  dplyr::mutate(plot_id = paste(scpca_library_id, seq_unit, sep = "-"))

# filter for cells that are found in both alevin + cellranger
cell_counts <- coldata_df %>%  
  dplyr::count(cell_id, sample)

common_cells <- cell_counts %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(cell_id)

coldata_common <- coldata_df %>%
  dplyr::filter(cell_id %in% common_cells) 

# create two separate dataframes used for plotting
cell_coldata_common <- coldata_common %>%
  dplyr::filter(seq_unit == "Cell")

nuclei_coldata_common <- coldata_common %>%
  dplyr::filter(seq_unit == "Nuclei")

# create combined UMI per cell plot 
umi_cell <- ggplot(cell_coldata_common, aes(x = sum, color = tool)) + 
  geom_density() + 
  facet_wrap(~ plot_id) +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "#56B4E9"),
        axis.text.x = element_text(size = 6)) +
  labs(x = "",
       y = "") +
  scale_color_manual(values = c("#D55E00", "#0072B2"))

umi_nuclei <- ggplot(nuclei_coldata_common, aes(x = sum, color = tool)) + 
  geom_density() + 
  facet_wrap(~ plot_id) +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "#E69F00"),
        axis.text.x = element_text(size = 6)) +
  labs(x = expression(paste(Log[10], " Total UMI per cell")),
       y = "") +
  scale_color_manual(values = c("#D55E00", "#0072B2"))

# combine plots and add combined legend 
combined_umi_plot <- umi_cell + umi_nuclei + plot_layout(guides = "collect", nrow = 2) & 
  theme(legend.position = "top", legend.title = element_blank())

# write combined plot to pdf
ggsave(umi_density, plot = combined_umi_plot, width = 8, height = 5)

# do the same thing for genes detected per cell 
genes_cell <- ggplot(cell_coldata_common, aes(x = detected, color = tool)) + 
  geom_density() + 
  facet_wrap(~ plot_id) +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "#56B4E9"),
        axis.text.x = element_text(size = 6)) +
  labs(x = "",
       y = "") +
  scale_color_manual(values = c("#D55E00", "#0072B2"))

genes_nuclei <- ggplot(nuclei_coldata_common, aes(x = detected, color = tool)) + 
  geom_density() + 
  facet_wrap(~ plot_id) +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "#E69F00"),
        axis.text.x = element_text(size = 6)) +
  labs(x = expression(paste(Log[10], " total genes detected per cell")),
       y = "") +
  scale_color_manual(values = c("#D55E00", "#0072B2"))

combined_gene_plot <- genes_cell + genes_nuclei + plot_layout(guides = "collect", nrow = 2) & 
  theme(legend.position = "top", legend.title = element_blank())

ggsave(gene_density, plot = combined_gene_plot, width = 8, height = 5)

## Mean gene expression correlation  
# grab rowdata from filtered sces
rowdata_df <- all_sces_list %>% 
  purrr::map_df(
    ~ purrr::map_df(.x, scpcaTools::rowdata_to_df, .id = "quant_id"),
    .id = "tool"
  ) %>%
  dplyr::mutate(sample = stringr::word(quant_id, 1, sep = "-"),
                seq_unit = dplyr::case_when(sample %in% single_cell ~ "Cell",
                                            sample %in% single_nuclei ~ "Nuclei")) %>%
  # join with metadata to get corresponding library ID 
  dplyr::left_join(select_metadata_df, by = c("sample" = "scpca_run_id")) %>%
  dplyr::mutate(plot_id = paste(scpca_library_id, seq_unit, sep = "-"))
  
# filter genes with low detection 
gene_counts <- rowdata_df %>% 
  # remove genes that have a low frequency of being detected
  dplyr::filter(detected >= 5.0) %>%
  dplyr::count(gene_id, sample)

# filter to only genes that are in both tools
common_genes <- gene_counts %>%
  dplyr::filter(n == 2) %>%
  dplyr::pull(gene_id)

rowdata_df_common <- rowdata_df %>%
  dplyr::filter(
    (gene_id %in% common_genes) 
  )

# spread table to put mean expression for Alevin-fry and Cell ranger in its own columns for plotting
rowdata_cor <- rowdata_df_common %>%
  dplyr::select(tool, gene_id, plot_id, mean) %>%
  # spread the mean expression stats to one column per caller
  tidyr::pivot_wider(id_cols = c(gene_id, plot_id),
                     names_from = c("tool"),
                     values_from = mean) %>%
  # drop rows with NA values to ease correlation calculations
  tidyr::drop_na()

# plot correlation of alevin-fry to Cellranger for each sample 
ggplot(rowdata_cor, aes(x = `Alevin-fry`, y = `Cell Ranger`)) +
  geom_point(size = 0.5, alpha = 0.1) + 
  facet_wrap(~ plot_id) + 
  scale_x_log10(labels = scales::label_number()) + 
  scale_y_log10(labels = scales::label_number()) + 
  labs(x = "Cell Ranger mean gene expression", y = "Alevin-fry mean gene expression") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 6),
        strip.background = element_rect())+
  stat_cor(aes(label = ..r.label..), method = "pearson", size = 2)
  
ggsave(gene_exp_correlation)
