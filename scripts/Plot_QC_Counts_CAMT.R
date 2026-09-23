##########################################################
# Script for plotting the read counts - pre or post QC ###
# Chanel Mosby-Tourtellot                              ###
##########################################################


# Load library
library(tidyr)
library(stringr)
library(readr)
library(ggpubr)
library(scales)

# Set working directory to your desired location
setwd("")

# Input file looks like this:
# Samples       `Raw reads` `Trimmed reads`
# <chr>              <dbl>           <dbl>
# 1 ISB_779_2022    41649056        27863072
# 2 ISB_787_2022    39127844        27193496
# 3 ISB_788_2022    45643686        34600884


QC_DF <- read_delim("QC_DCPHD_Hyb.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

QC_DF_long <- QC_DF %>% pivot_longer(
  cols = !Sample, names_to = "QC_state", values_to = "Read Counts")

#Adjust names as needed
QC_DF_long$Sample <- str_remove(QC_DF_long$Sample, "AFI_KEN_USAMRD_A_")

# Bar graph version
QC_barplot <- ggbarplot(QC_DF_long, "Sample", "Read Counts",
                        fill = "QC_state", color = "QC_state", palette = "Paired", 
                        position = position_dodge(0.8), title = "Read Counts for 75 USAFSAM Swabs - Semi-Agnostic") + 
  rotate_x_text() + font("x.text", size = 9) + theme(legend.title = element_blank(), legend.position = "right", 
                                                     plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::label_number(scale_cut = cut_short_scale()))


# Boxplot version
QC_barplot <- ggboxplot(QC_DF_long, "QC_state", "Read Counts", label = "Sample", label.select = list(top.up = 10), repel = TRUE,
                           fill = "QC_state", palette = c("#00AFBB", "#E7B800"), 
          title = "MIDRP EID Read counts: 397 samples + 17 controls") + 
  theme(legend.title = element_blank(), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::label_number(scale_cut = cut_short_scale()))
