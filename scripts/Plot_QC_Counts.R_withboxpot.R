#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggpubr)
})

# CLI options 
opt_list <- list(
  make_option(c("-i","--input"),   type="character", help="Input Excel file (.xlsx) [required]"),
  make_option(c("-s","--sheet"),   type="character", default=NULL, help="Worksheet name or index (default: first)"), # If you have a multisheeted excel file
  make_option(c("-o","--outdir"),  type="character", default=".",   help="Output directory (default: .)"),
  make_option(c("-f","--outfile"), type="character", default="qc_stats_raw_vs_trimmed.jpeg", help="Output JPEG filename"),
  make_option(c("--width"),        type="double",    default=12,    help="Plot width (inches)"),
  make_option(c("--height"),       type="double",    default=7,     help="Plot height (inches)"),
  make_option(c("--dpi"),          type="integer",   default=300,   help="DPI")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input)) stop("Missing --input")

# Read Excel 
sheet_arg <- if (is.null(opt$sheet)) 1 else opt$sheet
df <- readxl::read_excel(opt$input, sheet = sheet_arg)
Samples <- df$Sample
# Ensure required columns exist
needed <- c("Sample","Raw_Num_Reads","Trimmed_Num_Reads")
missing <- setdiff(needed, names(df))
if (length(missing) > 0) stop(sprintf("Missing required columns: %s", paste(missing, collapse=", ")))

# Select & clean 
plot_df <- df %>%
  select(Sample, Raw = Raw_Num_Reads, Trimmed = Trimmed_Num_Reads) %>%
  mutate(
    Raw     = as.numeric(gsub(",", "", as.character(Raw))),
    Trimmed = as.numeric(gsub(",", "", as.character(Trimmed)))
  )

# Long format for ggplot
plot_long <- plot_df %>%
  pivot_longer(cols = c(Raw, Trimmed), names_to = "Type", values_to = "Count") %>%
  mutate(
    Type   = factor(Type, levels = c("Raw","Trimmed")),
    Sample = factor(Sample, levels = unique(Sample))
  )
if (length(Samples) > 100 ) {
  #  Plot 
  p <- ggplot(plot_long, aes(x = Sample, y = Count, fill = Type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    # put labels on top of each bar
    geom_text(
      aes(label = scales::label_number(big.mark = ",")(Count)),
      position = position_dodge(width = 0.8),
      vjust = -0.25,
      size = 3
    ) +
    scale_y_continuous(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0.02, 0.10))  # add headroom for labels
    ) +
    scale_fill_manual(values = c(Raw = "#A6CEE3", Trimmed = "#1F78B4")) +
    labs(
      title = "Total Read Count: Raw Reads vs. Trimmed Reads",
      x = "Sample",
      y = "Read Count",
      fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title  = element_text(hjust = 0.5, face = "bold")
    )
} else {
  p <- ggboxplot(plot_long, "QC_state", "Read Counts", label = "Sample", label.select = list(top.up = 10), repel = TRUE,
                           fill = "QC_state", palette = c("#A6CEE3", "#1F78B4"), 
          title = "Total Read Count: Raw Reads vs. Trimmed Reads") + 
  theme(legend.title = element_blank(), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::label_number(scale_cut = cut_short_scale()))
}

# Save 
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(opt$outdir, opt$outfile)
ggsave(outfile, p, width = opt$width, height = opt$height, dpi = opt$dpi, device = "jpeg")
message("Wrote plot: ", outfile)
