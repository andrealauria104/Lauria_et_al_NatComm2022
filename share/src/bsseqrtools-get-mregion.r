#!/usr/bin/env Rscript
# # # # # # # # # # # # # # # # # # # # #
#                                       #
#  Bisulfite sequencing data analysis   #
#                                       #
# # # # # # # # # # # # # # # # # # # # #
#
# bsseqrtools-get-mregion.r
#
# Get methylation scores in predefined region
#
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(docopt)))

'Bisulfite sequencing data analysis
 
 bsseqrtools-get-region-methylation.r
 
 Get methylation scores in predefined region

Usage:
   bsseqrtools-get-mregion.r [-m <mcall> -r <region> -b <covbases> -c <coverage> -p <ratio> -w <writebed> -f <file> -t <threads>]

Options:
   -m, --mcall Path to mcall data (required).
   -r, --region Region bed file (required).
   -b, --covbases Minimum number of bases to be covered in region [default: 0].
   -c, --coverage Minimum coverage threshold to retain region [default: 8].
   -p, --ratio Call methylation ratio [default: TRUE].
   -w, --writebed Write filtered bed [default: TRUE].
   -f, --file Output file name [default: mregion.rds].
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
mcall   <- readRDS(opts$mcall)

regions <- BSseqRtools::build_other_annotation(opts$region)[[1]]
mregion <- methylKit::regionCounts(object = mcall
                                   , regions = regions
                                   , cov.bases = as.integer(opts$covbases)
                                   , mc.cores = as.integer(opts$threads))
if(as.integer(opts$coverage)!=0) {
  cov_df <- getData(mregion)[,mregion@coverage.index]  
  min_cov_idx <- which(rowSums(cov_df>=as.integer(opts$coverage))==ncol(cov_df))
  mregion <- mregion[min_cov_idx,]
}

if(as.logical(opts$writebed) & (as.integer(opts$coverage)!=0 | as.integer(opts$covbases)!=0)) {
  bed_df <- getData(mregion)[,c("chr","start","end")]
  utilsRtools::write_bedfile(bed_df, bedfile = gsub(".bed$",".coverage_filt.bed",opts$region))
}

outfile <- opts$file
message(" -- saving to: ", outfile)
saveRDS(mregion, file = outfile)

# 2. Methylation percentage ----
if(as.logical(opts$ratio)) {
  mregion_ratio <- methylKit::percMethylation(mregion, rowids = T)
  outfile <- gsub(".rds$",".ratio.rds",outfile)
  message(" -- saving to: ", outfile)
  saveRDS(mregion_ratio, file = outfile)
}
