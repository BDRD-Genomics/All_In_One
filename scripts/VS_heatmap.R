#MetaDetector summary heatmaps

### 
# Example Run
# Rscript VS_heatmaps.R -r TEST_total_hits_merged.csv -g Family -o output/dir
#
#
#
#
library(ggplot2)
library(reshape2)
library(plyr)
library(pheatmap)
library(scales)
library(optparse)

#process shotgun data for pneumonia samples

option_list <- list(
  make_option(opt_str = c("-r","--read_counts"),
              default = NULL,
              help = "csv file containing VS output for read counts",
              metavar = "character"),
  make_option(opt_str = c("-g","--group"),
              default = NULL,
              help = "which of the three read counts file is this?",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_dir"),
              type = "character",
              default = NULL,
              help = "output for heatmaps",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)

group <- args$group
outdir <- args$out_dir
read_counts_input <- args$read_counts
read_counts <- read.csv(read_counts_input, sep = ",", header = TRUE)
read_counts$X <- NULL

#DNA_reads_agg <- aggregate(cbind(name_list) ~ Family, data = DNA_reads, FUN = sum, na.rm = T)
DNA_reads_agg <- ddply(read_counts, .(Family), numcolwise(sum, na.rm=T))
#DNA_reads_agg <- aggregate(cbind(AFI_GEO_Sandflies_001_hostRem_normalizedFamily, AFI_GEO_Sandflies_002_hostRem_normalizedFamily) ~ Family, data = DNA_reads, FUN = sum, na.rm = TRUE)
#name_list <- colnames(DNA_reads_agg)
tmpcol1 <- gsub("_normalizedFamily","",as.character(colnames(DNA_reads_agg)))
tmpcol2 <- gsub("_normalizedRPM","",tmpcol1)
tmpcol3 <- gsub("_total","",tmpcol2)
colnames(DNA_reads_agg) <- tmpcol3
rownames(DNA_reads_agg) <- DNA_reads_agg$Family
#DNA_reads_agg
DNA_reads_agg$Family <- NULL
#DNA_reads_agg
#DNA_reads_agg <- DNA_reads_agg[rownames(DNA_reads_agg) != "Unclassified",]
#DNA_reads_agg

giveNAs <- which(is.na(as.matrix(dist(DNA_reads_agg))), arr.ind = T)
#head(giveNAs)


palette <- colour_ramp(c("white", "lightyellow",  "darkblue"))
#palette(seq(0,1,length.out=50))
#col <- heat.colors(30)

color_var <- palette(seq(0,1,length.out=50))
breaks_var <- seq(min(DNA_reads_agg),max(DNA_reads_agg),length.out=51)


#DNA_reads_agg
#min(DNA_reads_agg)
#max(DNA_reads_agg)

scaled_data <- scale(DNA_reads_agg)

png(paste0(outdir,group,"_vs_heatmap.png"),width=1200,height=800)
pheatmap(DNA_reads_agg, 
         #scale="row", # Filtered gene x sample matrix
         cluster_rows=F, # Cluster based on genes
         cluster_cols=F, # Cluster based on samples
         #annotation_col = my_sample_col,
         #color=inferno(10),
         color=color_var,
         breaks=breaks_var,
         #scale = "column",
         #display_numbers = T,
         show_colnames=T, # Show the sample names
         #annotation_colors=annot_colors,
         angle_col = 90,
         heatmap_legend_param = list(title = "Counts", # Legend for genes
             #direction = "horizontal",
             title_position = "topcenter"))
dev.off()

png(paste0(outdir,group,"_vs_scaleByFamily_heatmap.png"),width=1200,height=800)
pheatmap(DNA_reads_agg, 
         scale="row", # Filtered gene x sample matrix
         cluster_rows=F, # Cluster based on genes
         cluster_cols=F, # Cluster based on samples
         #annotation_col = my_sample_col,
         #color=inferno(10),
         color=color_var,
         breaks=breaks_var,
         #scale = "column",
         #display_numbers = T,
         show_colnames=T, # Show the sample names
         #annotation_colors=annot_colors,
         angle_col = 90,
         heatmap_legend_param = list(title = "Counts", # Legend for genes
                                     #direction = "horizontal",
                                     title_position = "topcenter"))
dev.off()

png(paste0(outdir,group,"_vs_scaleBySample_heatmap.png"),width=1200,height=800)
pheatmap(DNA_reads_agg, 
         #scale="row", # Filtered gene x sample matrix
         cluster_rows=F, # Cluster based on genes
         cluster_cols=F, # Cluster based on samples
         #annotation_col = my_sample_col,
         #color=inferno(10),
         color=color_var,
         breaks=breaks_var,
         scale = "column",
         #display_numbers = T,
         show_colnames=T, # Show the sample names
         #annotation_colors=annot_colors,
         angle_col = 90,
         heatmap_legend_param = list(title = "Counts", # Legend for genes
                                     #direction = "horizontal",
                                     title_position = "topcenter"))
dev.off()
