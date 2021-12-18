#!/usr/bin/env Rscript

# Run cellhash demux notebook for each sample

rmd_template <- "15-cellhash-demux.Rmd"
outdir <- file.path("results", "notebooks")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# samples have odd lib ids between 533 and 567
lib_ids <- paste0("SCPCL000", seq(533, 567, by = 2))

lib_ids |>
  purrr::walk(
    ~ try(rmarkdown::render(rmd_template,
                        output_dir = outdir,
                        output_file = paste0(., "_demux.html"),
                        params = list(lib_id = .),
                        envir = new.env()))
  )
