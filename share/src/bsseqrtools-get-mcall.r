#!/usr/bin/env Rscript
# # # # # # # # # # # # # # # # # # # # #
#                                       #
#  Bisulfite sequencing data analysis   #
#                                       #
# # # # # # # # # # # # # # # # # # # # #
#
# bsseqrtools-get-mcall.r
#
# Create data structure object for downstream analysis 
# with methylKit at single base resolution (mcall) 
#
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(docopt)))

'Bisulfite sequencing data analysis
 
 bsseqrtools-get-mcall.r
 
 Create data structure (mcall) object for downstream analysis with methylKit

Usage:
   bsseqrtools-get-mcall.r [-d <datadir> -m <metadata> -o <outdir> -c <coverage> -n <normalize> -p <pipeline> -g <genome> -h <hicov> -s <samplecov> -r <ratio> -f <file>]

Options:
   -d, --datadir Path to datadir (required).
   -m, --metadata Path to experiment metadata (required). 
   -o, --outdir Path to output directory [default: .].
   -c, --coverage Minimum coverage threshold to retain base [default: 5].
   -n, --normalize Normalize coverage [default: TRUE].
   -p, --pipeline Pipeline used for data processing (bismarkCoverage | bsmap) [default: bismarkCoverage].
   -g, --genome Genome assembly [default: mm10].
   -h, --hicov High coverage threshold (reccomended 99.9 for RRBS) [default: NULL].
   -s, --samplecov Minimum coverage in sample per condition [default: 1].
   -r, --ratio Call methylation ratio [default: TRUE].
   -f, --file Output file name [default: mcall.rds].

Author:
   Andrea Lauria' -> doc

opts <- docopt(doc)
# required arguments ---
required_args <- opts[1:2]
if(any(sapply(required_args, is.null))) {
  missing_idx <- sapply(required_args, is.null)
  missing_args <- gsub("--"," ",names(required_args)[missing_idx])
  message("\n[!] Missing required arguments: ",missing_args,"\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

# resources ----
suppressWarnings(suppressMessages(library(BSseqRtools)))

# 1. List files ----
if(as.character(opts$pipeline) == "bismarkCoverage") {
  fl <- list.files(opts$datadir, pattern = ".merged_CpG_evidence.cov.gz", full.names = T)
  names(fl) <- gsub(".trimmed_bismark.*","", basename(fl))
} else if(as.character(opts$pipeline) == "bsmap") {
  fl <- list.files(opts$datadir, pattern = ".txt(.gz)?", full.names = T)
  names(fl) <- gsub(".txt(.gz)?","", basename(fl))
}

if(length(fl)==0) {
  message("\n[!] No methylation data in the provided folder. Check pipeline settings.\n")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}
# 2. Read and process methylation data ----
metadata <- read.delim(as.character(opts$metadata))
ord_idx  <- metadata$sample.id[order(metadata$treatment)]
fl       <- fl[ord_idx]

mcall <- BSseqRtools::read_mcall(mfiles = fl
                                 , pipeline = opts$pipeline
                                 , experimental_info = metadata
                                 , cv = as.integer(opts$coverage)
                                 , normalize = as.logical(opts$normalize)
                                 , hiTh = opts$hicov
                                 , min = as.integer(opts$samplecov)
                                 , assembly = opts$genome)

outfile <- paste0(opts$outdir,"/",opts$file)
message(" -- saving to: ", outfile)
saveRDS(mcall, file = outfile)

# Methylation percentage 
if(as.logical(opts$ratio)) {
  mcall_ratio <- methylKit::percMethylation(mcall, rowids = T)
  outfile <- gsub(".rds$",".ratio.rds",outfile)
  message(" -- saving to: ", outfile)
  saveRDS(mcall_ratio, file = outfile)
}