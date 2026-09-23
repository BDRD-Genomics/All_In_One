#!/usr/bin/env Rscript

library(tidyr)
library(optparse)
library(stringr)

# ─── Option Parsing ─────────────────────────────────────────────
option_list <- list(
  make_option(c("-i", "--in_dir"), type = "character", default = NULL,
              help = "Input directory of assembly files", metavar = "character"),
  make_option(c("-o", "--out_dir"), type = "character", default = NULL,
              help = "Output directory for samplesheet", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

fasta_dir <- args$in_dir
outdir    <- args$out_dir

# ─── List FASTQ Files ───────────────────────────────────────────
files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = FALSE)

# ─── Parse File Names ───────────────────────────────────────────
parsed_data <- data.frame(full_name = files ) %>%
    separate(full_name, 
        into=c("sample_id", "assembler"),
        sep="_",
        extra="drop"
)
parsed_data[["fasta"]] = files

write.csv(parsed_data, "contig_samplesheet.csv", row.names=FALSE)
