#!/usr/bin/env Rscript
# Author: Adam Sorbie
# Date: 21/11/22
# Version 1.0.0

if(!require("pacman")){
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(dada2, optparse, parallel, ggplot2, tictoc)


option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="dada2_qc", 
              help="output folder", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Functions

RIGHT <- function(x,n){
  substring(x,nchar(x)-n+1)
}

# Start of script

# check all paths have trailing forward slash 
if (RIGHT(opt$path, 1) != "/") {
  opt$path <- paste0(opt$path, "/")
}
if (RIGHT(opt$out, 1) != "/") {
  opt$out <- paste0(opt$out, "/")
}


if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least three arguments must be supplied (path, trimF, trimR)", 
       call.=FALSE)
}

# time analysis
tic()
print(paste("QC STARTING", Sys.time(), sep=" "))

path <- opt$path 
setwd(opt$path)


fnFs <- sort(list.files(pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pattern="_R2_001.fastq.gz", full.names = TRUE))

# get sample names from forward reads 
sample.names <- sapply(strsplit(basename(fnFs), "_S1_L001_R1"), `[`, 1)

dir.create(opt$out)

print("Plotting quality profiles")
qc_F <- plotQualityProfile(fnFs, aggregate=TRUE)

ggsave(paste(opt$out,"quality_profile.pdf", sep="/"), qc_F, device = "pdf", dpi = 300)

qc_R <- plotQualityProfile(fnRs, aggregate=TRUE)
ggsave(paste(opt$out,"quality_profile_r.pdf", sep="/"), qc_R, device = "pdf", dpi = 300)

print(paste("QC COMPLETED", Sys.time(), sep=" "))

sink(paste(opt$out, "parameter_log.txt", sep="/"))
print(paste0("Filepath: ", opt$path))
print(paste0("Output: ", opt$out))
sessionInfo()
toc()
closeAllConnections() 
