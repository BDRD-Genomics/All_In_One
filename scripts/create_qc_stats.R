#!/bin/bash R

### Multiqc Parser Script

library(optparse)
library(scales)
library(openxlsx)
option_list <- list(
  make_option(opt_str = c("-i","--in_pre_trim"),
              default = NULL,
              help = "Input csv file pretrim multiqc csv",
              metavar = "character"),
  make_option(opt_str = c("-p","--in_post_trim"),
              default = NULL,
              help = "Input csv file posttrim multiqc csv",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_dir"),
              type = "character",
              default = NULL,
              help = "output for compiled_qc_stat_sheet",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

pre_trim_multiqc_csv     <- args$in_pre_trim
post_trim_multiqc_csv    <- args$in_post_trim
outdir                   <- args$out_dir

pre_trim_multiqc_csv <- read.csv(pre_trim_multiqc_csv)
post_trim_multiqc_csv <- read.csv(post_trim_multiqc_csv)
post_trim_multiqc_csv <- subset(post_trim_multiqc_csv,!grepl("_fastp_R1",Sample))
post_trim_multiqc_csv <- subset(post_trim_multiqc_csv,!grepl("_fastp_R2",Sample))
colnames(pre_trim_multiqc_csv) <- c("Sample", "Percent_Duplicates", "Percent_GC", "Avg_Sequence_Length","Median_Sequence_Length","Percent_Fails","Total_Sequences")
colnames(post_trim_multiqc_csv) <- c("Sample", "Percent_Duplicates", "Percent_GC", "Avg_Sequence_Length","Median_Sequence_Length","Percent_Fails","Total_Sequences")

# Num Reads
pre_trim_multiqc_csv$Num_Reads <- pre_trim_multiqc_csv$Total_Sequences * 10^6
post_trim_multiqc_csv$Num_Reads <- post_trim_multiqc_csv$Total_Sequences * 10^6
# Sum length
pre_trim_multiqc_csv$Sum_Length  <- (pre_trim_multiqc_csv$Num_Reads * pre_trim_multiqc_csv$Avg_Sequence_Length)
post_trim_multiqc_csv$Sum_Length <- (post_trim_multiqc_csv$Num_Reads * post_trim_multiqc_csv$Avg_Sequence_Length)
# Removed during QC
Full_Reads   <- (post_trim_multiqc_csv$Num_Reads - pre_trim_multiqc_csv$Num_Reads)
perc_Reads   <- Full_Reads / (post_trim_multiqc_csv$Num_Reads + pre_trim_multiqc_csv$Num_Reads )
Content      <- (post_trim_multiqc_csv$Sum_Length - pre_trim_multiqc_csv$Sum_Length)
perc_Content <- post_trim_multiqc_csv$Sum_Length / (post_trim_multiqc_csv$Sum_Length + pre_trim_multiqc_csv$Sum_Length)
GC_Delta     <- post_trim_multiqc_csv$Percent_GC - pre_trim_multiqc_csv$Percent_GC
#total_reads <- (post_trim_multiqc_csv$Num_Reads + pre_trim_multiqc_csv$Num_Reads)
#perc_Reads <- (Full_Reads/total_reads) * 100
total_df <- as.data.frame(cbind(pre_trim_multiqc_csv$Sample,comma(pre_trim_multiqc_csv$Num_Reads),comma(pre_trim_multiqc_csv$Sum_Length),
                                comma(round(pre_trim_multiqc_csv$Avg_Sequence_Length,digits=2)),pre_trim_multiqc_csv$Percent_GC,
                                comma(post_trim_multiqc_csv$Num_Reads),comma(post_trim_multiqc_csv$Sum_Length),
                                comma(round(post_trim_multiqc_csv$Avg_Sequence_Length,digits=2)),post_trim_multiqc_csv$Percent_GC,
                                comma(as.numeric(Full_Reads)),perc_Reads,Content,perc_Content,GC_Delta))
# rename columns
colnames(total_df) <- c("Sample","Raw_Num_Reads","Sum_Length","Avg_Length","GC",
                        "Trimmed_Num_Reads","Sum_Length","Avg_Length","GC",
                        "Full_Reads","% Reads","Content","% Content","GC Delta")
  
#row.names(total_df) <- pre_trim_multiqc_csv$Sample

write.table(total_df,paste0(outdir,"qc_stats_final.tsv"),quote=FALSE,row.names=FALSE,sep='\t')
write.xlsx(total_df,paste0(outdir,"qc_stats_final.xlsx"))
#write.table(total_df,paste0(outdir,"qc_stats_final.tsv"), sep='\t')
