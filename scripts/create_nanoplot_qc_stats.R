#!/usr/bin/env Rscript

### NanoPlot Parser Script
### Creates qc_stats_final.xlsx from NanoPlot NanoStats.txt files

library(optparse)
library(scales)
library(openxlsx)

option_list <- list(
  make_option(opt_str = c("-i", "--in_pre_trim"),
              default = NULL,
              help = "Input directory containing pretrim NanoPlot NanoStats.txt files",
              metavar = "character"),
  make_option(opt_str = c("-p", "--in_post_trim"),
              default = NULL,
              help = "Input directory containing posttrim NanoPlot NanoStats.txt files",
              metavar = "character"),
  make_option(opt_str = c("-o", "--out_dir"),
              type = "character",
              default = NULL,
              help = "Output directory for qc_stats_final.xlsx",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

pretrim_dir  <- args$in_pre_trim
posttrim_dir <- args$in_post_trim
outdir       <- args$out_dir

if (is.null(pretrim_dir) || is.null(posttrim_dir) || is.null(outdir)) {
  stop("ERROR: --in_pre_trim, --in_post_trim, and --out_dir are required")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

find_nanostats <- function(x) {
  files <- list.files(
    path = x,
    pattern = "NanoStats\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop(paste0("ERROR: No NanoStats.txt files found in: ", x))
  }

  return(files)
}

clean_number <- function(x) {
  x <- gsub(",", "", x)
  x <- gsub("%", "", x)
  x <- trimws(x)
  suppressWarnings(as.numeric(x))
}

get_value <- function(lines, patterns) {
  for (p in patterns) {
    hit <- grep(p, lines, ignore.case = TRUE, value = TRUE)

    if (length(hit) > 0) {
      value <- sub(".*:\\s*", "", hit[1])
      return(clean_number(value))
    }
  }

  return(NA_real_)
}

get_sample_name <- function(file) {
  parent <- basename(dirname(file))
  sample <- parent

  sample <- gsub("_pretrim_nanoplot$", "", sample)
  sample <- gsub("_posttrim_nanoplot$", "", sample)
  sample <- gsub("_nanoplot$", "", sample)
  sample <- gsub("_pretrim$", "", sample)
  sample <- gsub("_posttrim$", "", sample)
  sample <- gsub("_post_trim$", "", sample)

  return(sample)
}

parse_nanostats_file <- function(file) {
  lines <- readLines(file, warn = FALSE)

  sample <- get_sample_name(file)

  num_reads <- get_value(lines, c(
    "^Number of reads",
    "^Number of Reads",
    "^Reads"
  ))

  total_bases <- get_value(lines, c(
    "^Total bases",
    "^Total Bases",
    "^Total basepairs",
    "^Total Basepairs"
  ))

  mean_len <- get_value(lines, c(
    "^Mean read length",
    "^Mean Read Length",
    "^Mean length"
  ))

  median_len <- get_value(lines, c(
    "^Median read length",
    "^Median Read Length",
    "^Median length"
  ))

  n50 <- get_value(lines, c(
    "^Read length N50",
    "^N50"
  ))

  mean_qual <- get_value(lines, c(
    "^Mean read quality",
    "^Mean Read Quality",
    "^Mean quality"
  ))

  median_qual <- get_value(lines, c(
    "^Median read quality",
    "^Median Read Quality",
    "^Median quality"
  ))

  df <- data.frame(
    Sample = sample,
    Num_Reads = num_reads,
    Sum_Length = total_bases,
    Avg_Sequence_Length = mean_len,
    Median_Sequence_Length = median_len,
    Read_Length_N50 = n50,
    Mean_Read_Quality = mean_qual,
    Median_Read_Quality = median_qual,
    stringsAsFactors = FALSE
  )

  return(df)
}

pre_files  <- find_nanostats(pretrim_dir)
post_files <- find_nanostats(posttrim_dir)

pre_df  <- do.call(rbind, lapply(pre_files, parse_nanostats_file))
post_df <- do.call(rbind, lapply(post_files, parse_nanostats_file))

merged <- merge(
  pre_df,
  post_df,
  by = "Sample",
  suffixes = c("_Raw", "_Trimmed"),
  all = TRUE
)

merged$Reads_Removed <- merged$Num_Reads_Raw - merged$Num_Reads_Trimmed
merged$Percent_Reads_Removed <- ifelse(
  is.na(merged$Num_Reads_Raw) | merged$Num_Reads_Raw == 0,
  NA,
  round((merged$Reads_Removed / merged$Num_Reads_Raw) * 100, 2)
)

merged$Bases_Removed <- merged$Sum_Length_Raw - merged$Sum_Length_Trimmed
merged$Percent_Bases_Removed <- ifelse(
  is.na(merged$Sum_Length_Raw) | merged$Sum_Length_Raw == 0,
  NA,
  round((merged$Bases_Removed / merged$Sum_Length_Raw) * 100, 2)
)

merged$Mean_Length_Delta <- merged$Avg_Sequence_Length_Trimmed - merged$Avg_Sequence_Length_Raw
merged$N50_Delta <- merged$Read_Length_N50_Trimmed - merged$Read_Length_N50_Raw
merged$Mean_Quality_Delta <- merged$Mean_Read_Quality_Trimmed - merged$Mean_Read_Quality_Raw

total_df <- data.frame(
  Sample = merged$Sample,

  Raw_Num_Reads = comma(merged$Num_Reads_Raw),
  Raw_Sum_Length = comma(merged$Sum_Length_Raw),
  Raw_Avg_Length = comma(round(merged$Avg_Sequence_Length_Raw, digits = 2)),
  Raw_Median_Length = comma(round(merged$Median_Sequence_Length_Raw, digits = 2)),
  Raw_N50 = comma(merged$Read_Length_N50_Raw),
  Raw_Mean_Quality = round(merged$Mean_Read_Quality_Raw, digits = 2),

  Trimmed_Num_Reads = comma(merged$Num_Reads_Trimmed),
  Trimmed_Sum_Length = comma(merged$Sum_Length_Trimmed),
  Trimmed_Avg_Length = comma(round(merged$Avg_Sequence_Length_Trimmed, digits = 2)),
  Trimmed_Median_Length = comma(round(merged$Median_Sequence_Length_Trimmed, digits = 2)),
  Trimmed_N50 = comma(merged$Read_Length_N50_Trimmed),
  Trimmed_Mean_Quality = round(merged$Mean_Read_Quality_Trimmed, digits = 2),

  Reads_Removed = comma(merged$Reads_Removed),
  Percent_Reads_Removed = merged$Percent_Reads_Removed,

  Bases_Removed = comma(merged$Bases_Removed),
  Percent_Bases_Removed = merged$Percent_Bases_Removed,

  Mean_Length_Delta = round(merged$Mean_Length_Delta, digits = 2),
  N50_Delta = comma(merged$N50_Delta),
  Mean_Quality_Delta = round(merged$Mean_Quality_Delta, digits = 2),

  stringsAsFactors = FALSE
)

write.xlsx(
  total_df,
  file = file.path(outdir, "qc_stats_final.xlsx"),
  overwrite = TRUE
)

write.table(total_df,paste0(outdir,"qc_stats_final.tsv"), sep='\t')

