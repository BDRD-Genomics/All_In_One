#!/usr/bin/env Rscript

library(optparse)
library(stringr)

# ─── Option Parsing ─────────────────────────────────────────────
option_list <- list(
  make_option(c("-i", "--in_dir"), type = "character", default = NULL,
              help = "Input directory of FASTQ files", metavar = "character"),
  make_option(c("-o", "--out_dir"), type = "character", default = NULL,
              help = "Output directory for samplesheet", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

fastq_dir <- args$in_dir
outdir    <- args$out_dir

# ─── List FASTQ Files ───────────────────────────────────────────
files <- list.files(fastq_dir, pattern = "\\.fastq\\.gz$", full.names = FALSE)

# ─── Extract Unique Sample IDs ─────────────────────────────────
samples <- unique(str_replace(
  files,
  "_S\\d+_L\\d{3}_R[12]_001\\.fastq\\.gz$|_R[12]_001\\.fastq\\.gz$|_R[12]\\.fastq\\.gz$|\\.fastq\\.gz$",
  ""
))

# ─── Collect Reads Per Sample ──────────────────────────────────
sample    <- list()
fastq_1   <- list()
fastq_2   <- list()
long_read <- list()

for (sam in samples) {
  sample[[sam]]    <- sam
  #fastq_1[[sam]]   <- files[grepl(paste0("^", sam, "(_S([1-9][0-9]?|100)__R1_001|_R1_001|_R1)\\.fastq\\.gz$"), files)]
  #fastq_1[[sam]]   <- files[grepl(paste0("^", sam, "(_S([1-9][0-9]?|100)_R1_001|_R1_001|_R1)\\.fastq\\.gz$"), files)]
  #fastq_2[[sam]]   <- files[grepl(paste0("^", sam, "(_S[123456789]_L001_R2_001|_R2_001|_R2)\\.fastq\\.gz$"), files)]
  fastq_1[[sam]]   <- files[grepl(paste0("^", sam, ".*\\_R1.*\\.fastq\\.gz$"), files)]
  fastq_2[[sam]]   <- files[grepl(paste0("^", sam, ".*\\_R2.*\\.fastq\\.gz$"), files)]
  long_read[[sam]] <- files[grepl(paste0("^", sam, "\\.fastq\\.gz$"), files)]
}

# ─── Replace Missing Reads with Placeholder ────────────────────
fastq_1   <- lapply(fastq_1, function(x) if (length(x) == 0) "No_Read_R1" else x)
fastq_2   <- lapply(fastq_2, function(x) if (length(x) == 0) "No_Read_R2" else x)
long_read <- lapply(long_read, function(x) if (length(x) == 0) "No_Read" else x)

# ─── Build Final DataFrame ─────────────────────────────────────
df <- data.frame(
  sample_id       = unlist(sample),
  fastq_1         = file.path(fastq_dir, unlist(fastq_1)),
  fastq_2         = file.path(fastq_dir, unlist(fastq_2)),
  long_read       = file.path(fastq_dir, unlist(long_read))
  stringsAsFactors = FALSE
)

# ─── Write Output ──────────────────────────────────────────────
out_file <- file.path(outdir, "samplesheet.csv")
write.table(df, file = out_file, row.names = FALSE, sep = ",", quote = FALSE)

cat("Samplesheet written to:", out_file, "\n")

