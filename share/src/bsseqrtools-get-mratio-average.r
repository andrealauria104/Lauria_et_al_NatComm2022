#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

if(length(args)<3) {
  message("\n[!] Missing required arguments")
  message("Usage: bsseqrtools-get-mratio-average.r <path_mcall_ratio> <path_metadata> <metadata_features_summarize> <outfile> (optional)")
  quit(save = "no", status = 0, runLast = TRUE)
}

suppressWarnings(suppressMessages(library(BSseqRtools)))

path_mcall_ratio <- args[1]
path_metadata    <- args[2]
metadata_features_summarize <- args[3]

if(length(args)==3) {
  outfile <- gsub(".ratio.rds$",".ratio.average.rds",path_mcall_ratio)
} else {
  outfile <- args[4]
}

mcall_ratio <- readRDS(path_mcall_ratio)
metadata    <- read.delim(path_metadata)
mcall_ratio_average <- BSseqRtools::get_average_methylation_2(mdata_ratio = mcall_ratio
                                                              , metadata = metadata
                                                              , metadata_features_summarize = metadata_features_summarize)

message("-- saving to: ", outfile)
saveRDS(mcall_ratio_average, file = outfile)