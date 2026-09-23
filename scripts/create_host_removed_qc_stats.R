#!/bin/bash R

### Multiqc Parser Script

library(optparse)
library(scales)
library(openxlsx)
option_list <- list(
 # make_option(opt_str = c("-i","--qc_stats_pretrim_post_trim"),
 #             default = NULL,
 #             help = "Input csv file pretrim multiqc csv",
 #             metavar = "character"),
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

#pre_trim_post_trim_qc_stats_excel <- args$in_pre_trim
host_removed_multiqc_csv           <- args$in_post_trim
outdir                             <- args$out_dir

#pre_trim_post_trim_qc_stats_excel <- openxlsx::read.xlsx(pre_trim_post_trim_qc_stats_excel)
host_removed_multiqc_csv <- read.csv(host_removed_multiqc_csv)
host_removed_samples_idx <- grepl("_host_removed", host_removed_multiqc_csv$Sample, ignore.case=TRUE)
host_removed_samples <- host_removed_multiqc_csv[host_removed_samples_idx,]
host_mapped_samples_idx  <- !grepl("_host_removed", host_removed_multiqc_csv$Sample, ignore.case=TRUE)
host_mapped_samples <- host_removed_multiqc_csv[host_mapped_samples_idx,]
# strip host from sample name
host_mapped_samples$Sample  <- sub("_host.*", "", host_mapped_samples$Sample )
host_removed_samples$Sample <- sub("_host.*", "", host_removed_samples$Sample )

colnames(host_mapped_samples) <- c("Sample", "Percent_Duplicates", "Percent_GC", "Avg_Sequence_Length","Median_Sequence_Length","Percent_Fails","Total_Sequences")
colnames(host_removed_samples) <- c("Sample", "Percent_Duplicates", "Percent_GC", "Avg_Sequence_Length","Median_Sequence_Length","Percent_Fails","Total_Sequences")

# Num Reads
host_mapped_samples$Num_Reads <- host_mapped_samples$Total_Sequences * 10^6
host_removed_samples$Num_Reads <- host_removed_samples$Total_Sequences * 10^6
# Sum length
host_mapped_samples$Sum_Length  <- (host_mapped_samples$Num_Reads * host_mapped_samples$Avg_Sequence_Length)
host_removed_samples$Sum_Length <- (host_removed_samples$Num_Reads * host_removed_samples$Avg_Sequence_Length)
# Removed during QC
Full_Reads   <- (host_removed_samples$Num_Reads - host_mapped_samples$Num_Reads)
perc_Reads   <- Full_Reads / (host_removed_samples$Num_Reads + host_mapped_samples$Num_Reads )
Content      <- (host_removed_samples$Sum_Length - host_mapped_samples$Sum_Length)
perc_Content <- host_removed_samples$Sum_Length / (host_removed_samples$Sum_Length + host_mapped_samples$Sum_Length)
GC_Delta     <- host_removed_samples$Percent_GC - host_mapped_samples$Percent_GC
#total_reads <- (post_trim_multiqc_csv$Num_Reads + pre_trim_multiqc_csv$Num_Reads)
#perc_Reads <- (Full_Reads/total_reads) * 100
total_df <- as.data.frame(cbind(host_mapped_samples$Sample,comma(host_mapped_samples$Num_Reads),comma(host_mapped_samples$Sum_Length),
                                comma(round(host_mapped_samples$Avg_Sequence_Length,digits=2)),host_mapped_samples$Percent_GC,
                                comma(host_removed_samples$Num_Reads),comma(host_removed_samples$Sum_Length),
                                comma(round(host_removed_samples$Avg_Sequence_Length,digits=2)),host_removed_samples$Percent_GC,
                                comma(as.numeric(Full_Reads)),perc_Reads,Content,perc_Content,GC_Delta))
# rename columns
colnames(total_df) <- c("Sample","Host_Num_Reads","Sum_Length","Avg_Length","GC",
                        "Non_Host_Num_Reads","Sum_Length","Avg_Length","GC",
                        "Full_Reads","% Reads","Content","% Content","GC Delta")
  
#row.names(total_df) <- pre_trim_multiqc_csv$Sample

#write.table(total_df,paste0(outdir,"qc_stats_final.tsv"),quote=FALSE,row.names=FALSE,sep='\t')
write.xlsx(total_df,paste0(outdir,"qc_stats_host_removed.xlsx"))
