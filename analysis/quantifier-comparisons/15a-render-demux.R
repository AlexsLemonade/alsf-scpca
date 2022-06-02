#!/usr/bin/env Rscript

# This little script runs the demux comparisons notebook for each library 
# present in the data folder. Note that the main text will not change, 
# even though results can vary among libraries, sometimes widely.

rmd_template <- "15-demux-comparisons.Rmd"
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
