#!/bin/bash R

library(optparse)
option_list <- list(
  make_option(opt_str = c("-i","--in_dir"),
              default = NULL,
              help = "Input directory ffor trimmed interleaved files",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_dir"),
              type = "character",
              default = NULL,
              help = "output for samplesheet",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)


fastq_dir <- args$in_dir
outdir    <- args$out_dir

#files = list.files(fastq_dir)
files = list.files(fastq_dir,pattern="\\.fastq.gz$")


#samples = unique(stringr::str_remove(files,"_R1.fastq.gz|_R2.fastq.gz|.fastq.gz|_[0-9]+_S.*_L.*_R.*\\.fastq\\.gz$|_L001_R.*_001.*\\.fastq\\.gz$"))
#samples = unique(stringr::str_remove(files,"_R1_001.fastq.gz|_R2_001.fastq.gz|_R1.fastq.gz|_R2.fastq.gz|.fastq.gz|_L001_R.*_001.*\\.fastq\\.gz$"))
samples = unique(stringr::str_remove(files,"_trimmed_R1.fastq.gz|_trimmed_R2.fastq.gz|_trimmed.fastq.gz|.fastq.gz|_L001_R.*_001.*\\.fastq\\.gz$"))


sample    = list()
interleaved   = list()
for (sam in samples){
    sample[[sam]]  = sam
    interleaved[[sam]] = files[grepl(paste0("^",sam,".*_trimmed.fastq.gz$"),files)] 
    # If there is no read change character(0) to "No_Read"

}

data_dir <- fastq_dir
df <- data.frame(sample_id   = unlist(sample),
                 interleaved   = paste0(fastq_dir,unlist(interleaved)),
#                 all_fastq_files = paste0(fastq_dir,'[!Undetermined_]*.fastq.gz'),
                 data_dir = paste0(data_dir))

write.table(df,paste0(outdir,'/',"interleaved_samplesheet.csv"),row.names=FALSE,sep=",",quote=FALSE)
