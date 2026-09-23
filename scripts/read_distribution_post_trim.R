suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(optparse)
})


# --------------------
# Command-line options
# --------------------

option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "post-trimmed fastq(.gz) file", metavar = "file"),
    make_option(c("-p", "--prefix"), type = "character", default = "output", help = "Output file prefix [default = %default]"),
    make_option(c("-o", "--outdir"), type="character", default = ".", help = "Output directory [default = current directory]")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("Both --pre and --post FASTQ files are required.", call. = FALSE)
if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
}

# --------------------
# Function to extract read lengths 
# --------------------
get_read_lengths <- function(file) {
    lines <- vroom:vroom_lines(file)
    seq_lines <- lines[seq(2, length(lines), by = 4)]
    tibble(ReadLength = nchar(seq_lines))
}

# --------------------
# Load and summarize
# --------------------
read_lengths <- get_read_lengths(opt$input)

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

# save summary
write_csv(summary_df, file.path(opt$outdir, paste0(opt$prefix, "_read_length_summary.csv")))

# --------------------
# Plot histogram with lines
# --------------------

med <- summary_df$Median_Length
breakpoints <- c(500, 1500, 3000)

max_density <- max(density(read_lengths$ReadLength)$y)

p <- ggplot(read_lengths, aes(x = ReadLength)) +
    geom_histogram(aes(y = ..density..), position = "identity", bins = 100, fill = "#2c7bb6",
        alpha = 0.15, color = "black") + 
    geom_density(alpha = 0.1, adjust = 1.5) +  
    geom_vline(xintercept=med,color="black", linewidth=1.2) + 
    geom_text(aes(x = med, y = max_density * 1.05), label = paste0("Median: ", round(med)), angle = 90,
              fontface = "bold", color = "black", size = 4, hjust = 0) + 
    geom_vline(xintercept =  breakpoints, linetype = "dotted", color = "read", linewidth = 0.8) + 
    annotate("text", x = breakpoints, y = max_density * 1.05, 
             label = paste(breakpoints, "bp"),
             angle = 90, color = "red", fontface = "bold", size = 3.5, hjust = 0) + 
    labs(title = "Post-Trimmed Read Length Distribution",
         x = "Read Length (bp)", y = "Density") + 
    theme_bw()

# save plot
ggsave(file.path(opt$outdir, paste0(opt$prefix, "_read_length_distribution.png")),
       plot = p, width = 12, height = 5, dpi = 300)
cat("Done. Plot and summary saved in: ", opt$outdir, "\n ")