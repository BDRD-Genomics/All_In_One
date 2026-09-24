#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
  library(fs)
})

option_list <- list(
  make_option(c("-e","--run_id"),
              type="character", help="Project identifier (used in title and default output name)"),
  make_option(c("-p","--qc_plot"),
              type="character", help="Path to the QC JPEG produced by qc_stats_plot.R"),
  make_option(c("-r","--rmd"),
              type="character", default="ngs_summary_tabs.rmd",
              help="Path to the Rmd template (default: ngs_summary_tabs.rmd in CWD)"),
  make_option(c("-o","--out_html"),
              type="character", default=NULL,
              help="Output HTML path (default: <out_dir>/<run_id>_ngs_summary.html)"),
  make_option(c("-d","--out_dir"),
              type="character", default=".",
              help="Output directory (default: current directory)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate inputs 
if (is.null(opt$run_id)) {
  stop("Missing --run_id")
}
if (is.null(opt$qc_plot)) {
  stop("Missing --qc_plot (path to the QC JPEG)")
}
if (!file_exists(opt$rmd)) {
  stop(sprintf("Rmd not found: %s", opt$rmd))
}
if (!file_exists(opt$qc_plot)) {
  stop(sprintf("QC plot not found: %s", opt$qc_plot))
}

dir_create(opt$out_dir)

# Default output file if not provided
if (is.null(opt$out_html) || nchar(opt$out_html) == 0) {
  opt$out_html <- path(opt$out_dir, paste0(opt$run_id, "_ngs_summary.html"))
}

# ---- Render ----
rmarkdown::render(
  input  = opt$rmd,
  params = list(
    run_id = opt$run_id,
    qc_plot_path  = opt$qc_plot
  ),
  output_file = opt$out_html,
  quiet = TRUE
)

cat("Wrote HTML:", opt$out_html, "\n")
