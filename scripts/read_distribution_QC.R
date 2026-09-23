suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(optparse)
})

# Define Options
option_list <- list(
    make_option(c("-a", "--pre"), type = "character", help = "pre-trimmed fastq(.gz) file", metavar = "file"),
    make_option(c("-b", "--post"), type = "character", help = "post-trimmed fastq(.gz) file", metavar = "file"),
    make_option(c("-p", "--prefix"), type = "character", default = "output", help = "Output file prefix [default = %default]"),
    make_option(c("-o", "--outdir"), type="character", default = ".", help = "Output directory [default = current directory]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$pre) || is.null(opt$post)) {
    stop("Both --pre and --post FASTQ files are required.", call. = FALSE)
}

if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
}
# function to extract read length from fastq file
get_read_lengths_vroom <- function(file,label) {
    message("Reading file: ", file)

    lines <- vroom::vroom_lines(file)

    seq_lines <- lines[seq(2, length(lines), by = 4)]

    tibble(ReadLength = nchar(seq_lines), Stage = label)
}

# Summarize length

summarize_lengths <- function(df, label) {
    df %>%
        summarize(
            Label = label,
            Total_Reads = n(),
            Mean_Length = mean(ReadLength),
            Median_Length = median(ReadLength),
            Q1 = quantile(ReadLength, 0.25),
            Q2 = quantile(ReadLength, 0.5),
            Q3 = quantile(ReadLength, 0.75),
            Q4 = quantile(ReadLength, 1),
            Min_Length = min(ReadLength),
            Max_Length = max(ReadLength),
            Std_Dev = sd(ReadLength)
        )
}

# Run
cat("Processing files...\n")
pre_df <- get_read_lengths_vroom(opt$pre, "Pre-Trimmed")
post_df <- get_read_lengths_vroom(opt$post, "Post-Trimmed")
combined_df <- bind_rows(pre_df,post_df)

# calculate median lines for each stage
medians <- combined_df %>%
		group_by(Stage) %>%
		summarize(median = median(ReadLength), .groups = "drop")
# reference break points for min frag length
breakpoints <- c(500, 1500, 3000)

# output path
plot_path <- file.path(opt$outdir, paste0(opt$prefix, "_read_length_distribution.png"))
summary_path <- file.path(opt$outdir, paste0(opt$prefix, "_read_length_summary.csv"))

#max_density <- combined_df %>% 
#		group_by(Stage) %>%
#		summarise(max_dens = max(density(ReadLength)$y)) %>%
#		pull(max_dens) %>%
#		max()
max_y <- combined_df %>% 
	group_by(Stage) %>%
	summarise(d = max(density(ReadLength)$y)) %>% 
	pull(d) %>%
	max()
max_density <- combined_df %>% 
	group_by(Stage) %>%
	summarise(max_y = max(density(ReadLength)$y)) %>%
	pull(max_y) %>%
	max()
label_height <- max_density * 1.1
# plots
p <- ggplot(combined_df, aes(x = ReadLength, fill = Stage)) +
    geom_histogram(aes(y = ..density..), position = "identity", bins = 100,
        alpha = 0.15, color = "black") + 
    geom_density(alpha = 0.1, adjust = 1.5) +  
    geom_vline(xintercept = breakpoints, linetype = "dotted", color = "gray30", linewidth = 1.1) +
    annotate("text", x = breakpoints, y = label_height, vjust = 2, angle = 90,
             label = paste(breakpoints, "bp"), size = 4, color = "gray30", fontface = "bold", hjust = 0) +
    scale_fill_manual(values = c("Pre-Trimmed" = "#1f77b4", "Post-Trimmed" = "#ff7f0e")) +
    scale_color_manual(values = c("Pre-Trimmed" = "#1f77b4", "Post-Trimmed" = "#ff7f0e")) +
    labs(
        title = "Read Length Distribution with Filtering Reference Lines",
        x = "Read Length (bp)",
        y = "Density"
    ) +
    theme_bw() + 
    theme(
        legend.title = element_blank(),
        legend.position = "top"
    ) 
ggsave(plot_path, plot = p, width = 14, height = 16, dpi = 300)

# Summary
summary_df <- bind_rows(
    summarize_lengths(pre_df, "Pre-Trimmed"),
    summarize_lengths(post_df, "Post-Trimmed")
)
write_csv(summary_df, summary_path)

cat("Done.\n")
cat("Saved Plot to:", plot_path, "\n")
cat("Saved summary to:", summary_path, "\n")
