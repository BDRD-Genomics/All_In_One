#!/bin/bash R

### Multiqc Parser Script

library(optparse)
library(scales)
library(openxlsx)
option_list <- list(
  make_option(opt_str = c("-r","--in_rRNA_removed_csv"),
              default = NULL,
              help = "Input csv rRNA Removed multiqc csv",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_dir"),
              type = "character",
              default = NULL,
              help = "output for compiled_qc_stat_sheet",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

rRNA_removed_multiqc_csv           <- args$in_rRNA_removed_csv
outdir                             <- args$out_dir
# rRNA Removed
rRNA_removed_multiqc_csv <- read.csv(rRNA_removed_multiqc_csv)
rRNA_removed_idx <- grepl("_host_contaminant_", rRNA_removed_multiqc_csv$Sample, ignore.case=TRUE)
rRNA_removed_samples <- rRNA_removed_multiqc_csv[rRNA_removed_idx,]
# rRNA 
rRNA_reads_idx  <- !grepl("_host_contaminant_", rRNA_removed_multiqc_csv$Sample, ignore.case=TRUE)
rRNA_reads <- rRNA_removed_multiqc_csv[rRNA_reads_idx,]
# strip host from sample name
rRNA_removed_samples$Sample  <- sub("_host.*", "", rRNA_removed_samples$Sample )
rRNA_reads$Sample <- sub("_rRNA.*", "", rRNA_reads$Sample )

colnames(rRNA_removed_samples) <- c("Sample", "Percent_Duplicates", "Percent_GC", "Avg_Sequence_Length","Median_Sequence_Length","Percent_Fails","Total_Sequences")
colnames(rRNA_reads) <- c("Sample", "Percent_Duplicates", "Percent_GC", "Avg_Sequence_Length","Median_Sequence_Length","Percent_Fails","Total_Sequences")

# Num Reads
rRNA_reads$Num_Reads <- rRNA_reads$Total_Sequences * 10^6
rRNA_removed_samples$Num_Reads <- rRNA_removed_samples$Total_Sequences * 10^6
# Sum length
rRNA_reads$Sum_Length  <- (rRNA_reads$Num_Reads * rRNA_reads$Avg_Sequence_Length)
rRNA_removed_samples$Sum_Length <- (rRNA_removed_samples$Num_Reads * rRNA_removed_samples$Avg_Sequence_Length)
# Removed during QC
Full_Reads   <- (rRNA_removed_samples$Num_Reads - rRNA_reads$Num_Reads)
perc_Reads   <- Full_Reads / (rRNA_removed_samples$Num_Reads + rRNA_reads$Num_Reads )
Content      <- (rRNA_removed_samples$Sum_Length - rRNA_reads$Sum_Length)
perc_Content <- rRNA_removed_samples$Sum_Length / (rRNA_removed_samples$Sum_Length + rRNA_reads$Sum_Length)
GC_Delta     <- rRNA_removed_samples$Percent_GC - rRNA_reads$Percent_GC
#total_reads <- (post_trim_multiqc_csv$Num_Reads + pre_trim_multiqc_csv$Num_Reads)
#perc_Reads <- (Full_Reads/total_reads) * 100
total_df <- as.data.frame(cbind(rRNA_reads$Sample,comma(rRNA_reads$Num_Reads),comma(rRNA_reads$Sum_Length),
                                comma(round(rRNA_reads$Avg_Sequence_Length,digits=2)),rRNA_reads$Percent_GC,
                                comma(rRNA_removed_samples$Num_Reads),comma(rRNA_removed_samples$Sum_Length),
                                comma(round(rRNA_removed_samples$Avg_Sequence_Length,digits=2)),rRNA_removed_samples$Percent_GC,
                                comma(as.numeric(Full_Reads)),perc_Reads,Content,perc_Content,GC_Delta))
# rename columns
colnames(total_df) <- c("Sample","rRNA_Num_Reads","Sum_Length","Avg_Length","GC",
                        "Non_rRNA_Num_Reads","Sum_Length","Avg_Length","GC",
                        "Full_Reads","% Reads","Content","% Content","GC Delta")
  
#row.names(total_df) <- pre_trim_multiqc_csv$Sample

#write.table(total_df,paste0(outdir,"qc_stats_final.tsv"),quote=FALSE,row.names=FALSE,sep='\t')
write.xlsx(total_df,paste0(outdir,"qc_stats_rRNA_removed.xlsx"))
