#!/bin/bash R

library(optparse)
option_list <- list(
  make_option(opt_str = c("-i","--in_dir"),
              default = NULL,
              help = "Input directory fastq files",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_dir"),
              type = "character",
              default = NULL,
              help = "output for samplesheet",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)


pod5_dir <- args$in_dir
outdir   <- args$out_dir

files   = list.files(pod5_dir)
files_df <- as.data.frame(files)
samples = unique(stringr::str_remove(files,".pod5$"))
sample  = list()
pod5    = list()

for (sam in samples){
    sample[[sam]] = sam
    pod5[[sam]] = files_df[match(paste0(sam,".pod5"),files_df$files),]
}

df <- data.frame(sample_id   = unlist(sample),
                 pod5_file   = paste0(pod5_dir,unlist(pod5)))

write.table(df,paste0(outdir,'/',"pod5_list.csv"),row.names=FALSE,sep=",",quote=FALSE)
