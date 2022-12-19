#!/usr/bin/env Rscript
#
#  Call Differentially Methylated Regions (DMR)
#
# 0. Resources ----
suppressWarnings(suppressMessages(library(docopt)))

'Call Differentially Methylated Regions (DMR) with DSS

Usage:
   call-DMR.R [-d <dss>] [-o <outdir>] [--dmr_delta <dmr_delta>] [--dmr_p_threshold <dmr_p_threshold>] [--dmr_minlen <dmr_minlen>] [--dmr_minCG <dmr_minCG>] [--dmr_dis_merge <dmr_dis_merge>] [--dmr_pct_sig <dmr_pct_sig>]

Options:
   -d, --dss Path to DSS analysis object (as .rds file) (required)
   -o, --outdir Path to output directory [default: .]
   --dmr_delta DSS parameter dmr_delta [default: 0.2]
   --dmr_p_threshold DSS parameter dmr_p_threshold [default: 0.05]
   --dmr_minlen DSS parameter dmr_minlen [default: 50]
   --dmr_minCG DSS parameter dmr_minCG [default: 3]
   --dmr_dis_merge DSS parameter dmr_dis_merge [default: 100]
   --dmr_pct_sig DSS parameter dmr_pct_sig [default: 0.5]
' -> doc

opts <- docopt(doc)
# required arguments ---
if(is.null(opts$dss)) {
  message("\n [!] Missing required arguments \n ")
  message(doc)
  quit(save = "no", status = 0, runLast = TRUE)
}

# read data ----
message(" -- reading DSS analysis results ...")

dss_analysis <- readRDS(opts$dss)
dmr_analysis <- list()

# call DMR ----
for(i in names(dss_analysis)) {

  message("\n[+] callDMR, parameters:")
  message(" -- dmr.delta = ", as.numeric(opts$dmr_delta))
  message(" -- dmr.p.threshold  = ", as.numeric(opts$dmr_p_threshold))
  message(" -- dmr.minlen = ", as.integer(opts$dmr_minlen))
  message(" -- dmr.minCG = ", as.integer(opts$dmr_minCG))
  message(" -- dmr.dis.merge = ", as.integer(opts$dmr_dis_merge))
  message(" -- dmr.pct.sig = ", as.numeric(opts$dmr_pct_sig))
  
  dmrs <- DSS::callDMR(dss_analysis[[i]]$dmlTest
	, delta = as.numeric(opts$dmr_delta)
	, p.threshold = as.numeric(opts$dmr_p_threshold)
	, minlen = as.integer(opts$dmr_minlen)
	, minCG = as.integer(opts$dmr_minCG)
	, dis.merge = as.integer(opts$dmr_dis_merge)
	, pct.sig = as.numeric(opts$dmr_pct_sig))
  
  options(scipen = 999)
  output_suffix <- paste0("delta_", as.numeric(opts$dmr_delta), "_p_", as.numeric(opts$dmr_p_threshold)
                    ,".",as.integer(opts$dmr_minlen),"_",as.integer(opts$dmr_minCG)
                    ,"_",as.integer(opts$dmr_dis_merge),"_",as.numeric(opts$dmr_pct_sig))

  output_suffix <- paste0(i, "_", output_suffix)

  outfile <- paste0(opts$outdir, "/dssDMR_", output_suffix, ".txt")
  message("[+] writing DMR result to file: ", outfile)
  write.table(dmrs, file = outfile, row.names = F, col.names = T, sep = "\t", quote = F)
  options(scipen = 0)
}

