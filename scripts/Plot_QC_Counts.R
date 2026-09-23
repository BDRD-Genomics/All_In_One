#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# CLI options 
opt_list <- list(
  make_option(c("-i","--input"),   type="character", help="Input Excel file (.xlsx) [required]"),
  make_option(c("-s","--sheet"),   type="character", default=NULL, help="Worksheet name or index (default: first)"), # If you have a multisheeted excel file
  make_option(c("-o","--outdir"),  type="character", default=".",   help="Output directory (default: .)"),
  make_option(c("-f","--outfile"), type="character", default="qc_stats_raw_vs_trimmed_w_numbers.jpeg", help="Output JPEG filename"),
  make_option(c("-r","--outfile2"), type="character", default="qc_stats_raw_vs_trimmed_wo_numbers.jpeg", help="Output JPEG filename"),
  make_option(c("--width"),        type="double",    default=12,    help="Plot width (inches)"),
  make_option(c("--height"),       type="double",    default=7,     help="Plot height (inches)"),
  make_option(c("--dpi"),          type="integer",   default=300,   help="DPI")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input)) stop("Missing --input")

# Read Excel 
sheet_arg <- if (is.null(opt$sheet)) 1 else opt$sheet
df <- readxl::read_excel(opt$input, sheet = sheet_arg)

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

#  Plot 
p <- ggplot(plot_long, aes(x = Sample, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  # put labels on top of each bar
  geom_text(
    aes(label = scales::label_number(big.mark = ",")(Count)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    hjust = 0.01,
    angle = 45,
    size = 3
  ) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0, 0.10))  # add headroom for labels
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


#  Plot 
p2 <- ggplot(plot_long, aes(x = Sample, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  # (removed geom_text() labels over bars)
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale(), big.mark = ","),
    breaks = scales::pretty_breaks(n = 10),    
    expand = expansion(mult = c(0, 0.05))    
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



# Save 
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(opt$outdir, opt$outfile)
outfile2 <- file.path(opt$outdir, opt$outfile2)

ggsave(outfile, p, width = opt$width, height = opt$height, dpi = opt$dpi, device = "jpeg")
ggsave(outfile2, p2, width = opt$width, height = opt$height, dpi = opt$dpi, device = "jpeg")
message("Wrote plot: ", outfile)
message("Wrote plot: ", outfile2)

