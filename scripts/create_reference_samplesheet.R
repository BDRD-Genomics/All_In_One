# create samplesheet of references

library(optparse)
option_list <- list(
  make_option(opt_str = c("-i","--in_dir"),
              default = NULL,
              help = "Input directory fastq files",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_dir"),
              type = "character",
              default = NULL,
              help = "output for samplesheet",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)


ref_dir <- args$in_dir
outdir    <- args$out_dir

#files = list.files(ref_dir)
files = list.files(ref_dir,pattern="\\.fasta$|\\.fna$|\\.fa$|\\.fsa$")


#samples = unique(stringr::str_remove(files,"_R1.fastq.gz|_R2.fastq.gz|.fastq.gz|_[0-9]+_S.*_L.*_R.*\\.fastq\\.gz$|_L001_R.*_001.*\\.fastq\\.gz$"))
accessions = unique(stringr::str_remove(files,"\\.fasta$|\\.fna$|\\.fa$|\\.fsa$"))


acc_id    = list()
for (acc in accessions){
  acc_id[[acc]]  = acc
  }
reference_dir   = list()

files
for (f in files){
  reference_dir[[f]] = f
}
reference_dir
#test <- lapply(long_read, function(x) if(identical(x,character(0))) "No_Read" else x)
#data_dir <- ref_dir
df <- data.frame(acc_id   = unlist(acc_id),
                 reference_dir   = paste0(ref_dir,unlist(reference_dir)))


write.table(df,paste0(outdir,'/',"reference_samplesheet.csv"),row.names=FALSE,sep=",",quote=FALSE)
