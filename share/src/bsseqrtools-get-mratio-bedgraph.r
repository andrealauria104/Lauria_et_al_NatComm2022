#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = T)

if(length(args)<1) {
  message("\n[!] Missing required arguments")
  message("Usage: bsseqrtools-get-mratio-bedgraph.r <path_mcall_ratio> <outdir> (optional)")
  quit(save = "no", status = 0, runLast = TRUE)
}

suppressWarnings(suppressMessages(library(BSseqRtools)))

path_mcall_ratio <- args[1]

if(length(args)==1) {
  outdir <- "."
} else {
  outdir <- args[2]
}

mcall_ratio <- readRDS(path_mcall_ratio)
for(i in colnames(mcall_ratio)) {
  outfile <- paste0(outdir, "/",i,".",gsub(".rds","",basename(path_mcall_ratio)),".bedGraph")  
  BSseqRtools::create_begraph(mtr.ratio = mcall_ratio, outfile = outfile, samp = i, make.zero.based=F)
  
}
