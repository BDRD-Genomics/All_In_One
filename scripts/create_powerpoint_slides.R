#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(officer)
  library(magick)    # to read image dimensions
  library(fs)
})

# CLI 
opt_list <- list(
  make_option(c("-p","--run_id"), type="character", help="Run ID to show on title slide (required)"),
  make_option(c("-i","--image"),      type="character", help="Path to JPEG(s)"),
  make_option(c("-o","--out"),        type="character", default=NULL,
              help="Output .pptx path (default: ./<run_id>_QC_Summary.pptx)"),
  make_option(c("-t","--template"),   type="character", default=NULL,
              help="Optional PowerPoint template .pptx to use (default: blank Office Theme)") # Need to find NMRC template slides
)
opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$run_id)) stop("Missing --run_id")
if (is.null(opt$image))      stop("Missing --image (path to JPEG)")

if (!file_exists(opt$image)) stop(sprintf("Image not found: %s", opt$image))

if (is.null(opt$out) || nchar(opt$out) == 0) {
  opt$out <- path(".", sprintf("%s_QC_Summary.pptx", opt$run_id))
}
dir_create(path_dir(opt$out))

# Helper functions
add_title_slide <- function(doc, title_text) {
  layouts <- layout_summary(doc)$layout
  layout <- if ("Title Slide" %in% layouts) "Title Slide"
        else if ("Title and Content" %in% layouts) "Title and Content"
        else if ("Title Only" %in% layouts) "Title Only"
        else "Blank"

  doc <- add_slide(doc, layout = layout, master = "Office Theme")

  ph_types <- c("ctrTitle", "title")
  placed <- FALSE
  for (tp in ph_types) {
    try({
      doc <- ph_with(doc, value = title_text, location = ph_location_type(type = tp))
      placed <- TRUE
    }, silent = TRUE)
    if (placed) break
  }

  if (!placed) {
    sz <- slide_size(doc)
    w <- sz$width * 0.9
    h <- 1.5
    left <- (sz$width - w) / 2
    top  <- (sz$height - h) / 2
    doc <- ph_with(
      doc,
      value = fpar(ftext(title_text, fp_text(font.size = 36, bold = TRUE))),
      location = ph_location(left = left, top = top, width = w, height = h)
    )
  }

  doc
}

add_image_slide <- function(doc, img_path, slide_title = "Read Counts: Raw Vs. Trimmed",
                            layout = "Title and Content") {
  layouts <- layout_summary(doc)$layout
  layout <- if (layout %in% layouts) layout else if ("Title Only" %in% layouts) "Title Only" else "Blank"
  doc <- add_slide(doc, layout = layout, master = "Office Theme")

  
  placed_title <- FALSE
  for (tp in c("title", "ctrTitle")) {
    try({
      doc <- ph_with(doc, value = slide_title, location = ph_location_type(type = tp))
      placed_title <- TRUE
    }, silent = TRUE)
    if (placed_title) break
  }

  sz <- slide_size(doc)   # inches
  info <- magick::image_info(magick::image_read(img_path))
  img_w_in <- info$width / 96   # assume 96 dpi
  img_h_in <- info$height / 96

  top_margin <- if (placed_title) 1.0 else 0.25
  margin <- 0.25
  max_w <- sz$width  - 2*margin
  max_h <- sz$height - (top_margin + margin)

  scale <- min(max_w / img_w_in, max_h / img_h_in)
  w <- img_w_in * scale
  h <- img_h_in * scale

  left <- (sz$width - w) / 2
  top  <- top_margin + (max_h - h)/2

  doc <- ph_with(
    doc,
    value = external_img(img_path, width = w, height = h),
    location = ph_location(left = left, top = top, width = w, height = h)
  )

  doc
}

# Build PPTX 
doc <- if (!is.null(opt$template) && file_exists(opt$template)) read_pptx(opt$template) else read_pptx()

doc <- add_title_slide(doc, opt$run_id)
doc <- add_image_slide(doc, opt$image, slide_title = "Read Counts: Raw Vs. Trimmed")

print(doc, target = opt$out)
cat("Wrote PowerPoint:", opt$out, "\n")
