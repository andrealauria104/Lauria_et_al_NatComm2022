#!/usr/bin/env Rscript
# # # # # # # # # # # # # # # # # # # # #
#                                       #
#  Bisulfite sequencing data analysis   #
#                                       #
# # # # # # # # # # # # # # # # # # # # #
#
# bsseqrtools-get-mtile.r
#
# Get methylation scores in tiling windows
#
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(docopt)))

'Bisulfite sequencing data analysis
 
 bsseqrtools-get-mtile.r
 
 Get methylation scores in tiling windows.

Usage:
   bsseqrtools-get-mtile.r [-m <mcall> -o <outdir> -w <window> -s <step> -c <covbases> -r <ratio> -f <file> -t <threads>]

Options:
   -m, --mcall Path to mcall data (required). 
   -o, --outdir Path to output directory [default: .].
   -w, --window Window size [default: 1000].
   -s, --step Step size [default: 1000].
   -c, --covbases Minimum number of bases to be covered in window [default: 0].
   -r, --ratio Call methylation ratio [default: TRUE].
   -f, --file Output file name [default: mtile.rds].
   -t, --threads Number of threads [default: 1].

Author:
   Andrea Lauria' -> doc

opts <- docopt(doc)
# required arguments ---
required_args <- opts[1]
if(any(sapply(required_args, is.null))) {
  missing_idx <- sapply(required_args, is.null)
  missing_args <- gsub("--"," ",names(required_args)[missing_idx])
  message("\n[!] Missing required arguments: ",missing_args,"\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

# resources ----
suppressWarnings(suppressMessages(library(BSseqRtools)))

# 1. Read and process methylation data ----
mcall <- readRDS(opts$mcall)

# Genomic Tiles ---
mtile <- methylKit::tileMethylCounts(mcall
                                     , win.size = as.integer(opts$window)
                                     , step.size = as.integer(opts$step)
                                     , cov.bases = as.integer(opts$covbases)
                                     , mc.cores = as.integer(opts$threads))

outfile <- paste0(opts$outdir,"/",opts$file)
message(" -- saving to: ", outfile)
saveRDS(mtile, file = outfile)

# Methylation percentage 
if(as.logical(opts$ratio)) {
  mtile_ratio <- methylKit::percMethylation(mtile, rowids = T)
  outfile <- gsub(".rds$",".ratio.rds",outfile)
  message(" -- saving to: ", outfile)
  saveRDS(mtile_ratio, file = outfile)
}
