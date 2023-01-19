#!/usr/bin/env Rscript
# # # # # # # # # # # # # # # # # #
#                                 #
#  RNA-sequencing data analysis   #
#                                 #
# # # # # # # # # # # # # # # # # #
#
# rnaseqrtools-dge-edger.r
#
# Perfom Differential Expression Analysis with edgeR
# 
# 0. Resources, variables and parameters ----
#
# Mandatory variables to be assigned before running the script:
#   - path_dgelist_obj <- ASSIGN_PATH_DGELIST_OBJ
#     path to dgelist object with raw counts and metadata (e.g. obtained from rnaseq-get_dgelist_obj.r script)
#   - path_results    <- ASSIGN_PATH_RESULTS
#     path to results folder
# 
suppressWarnings(suppressMessages(library(RNAseqRtools)))

# paths ---
path_dgelist_obj <- "data/dgelist_obj.rds" # ASSIGN_PATH_DGELIST_OBJ
path_results <- "results/dge-edger" # ASSIGN_PATH_RESULTS

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# Differential gene expression analysis ---
de_method = "qlf" # qlf, lrt, exact
de_robust.dispersion = TRUE
de_return.y = FALSE
de_formula = "~0+StageOverexpression"
de_cf = NULL
de_contrast  = list("EpiLC_3B_vs_none"=c("StageOverexpression","EpiLC_3B","EpiLC_none")
                    ,"ME24h_3B_vs_none"=c("StageOverexpression","ME24h_3B","ME24h_none")
                    ,"ME48h_3B_vs_none"=c("StageOverexpression","ME48h_3B","ME48h_none"))

de_anovalike = FALSE
de_analysis  = "timecourse_contrasts" # dge

fdrTh  = 0.05
fcTh   = .5
save_table = TRUE
save_excel = FALSE

# read object ---
dgelist_obj <- readRDS(path_dgelist_obj)

s33d = 1991
# 1. Differential expression analysis ----
if(!is.null(de_contrast) && grepl("~0+",de_formula)) {
  # test defined contrasts (parametrization: ~ 0 + group ...) ---
  if(!is.list(de_contrast)) de_contrast <- list("de_contrast"=de_contrast)
  de <- lapply(de_contrast, function(x) {
    RNAseqRtools::calculateDiffExprEdgeR(y = dgelist_obj
                                         , design = de_formula
                                         , cf = de_cf
                                         , method = de_method
                                         , contrast = x
                                         , robust.dispersion = de_robust.dispersion
                                         , return.y = de_return.y
                                         , anovalike = de_anovalike)
  })
} else {
  # test model coefficients (parametrization: ~ group ... )---
  if(!is.list(de_cf)) de_cf <- list("de_cf"=de_cf)
  de <- lapply(de_cf, function(x) {
    RNAseqRtools::calculateDiffExprEdgeR(y = dgelist_obj
                                         , design = de_formula
                                         , cf = x
                                         , method = de_method
                                         , contrast = de_contrast
                                         , robust.dispersion = de_robust.dispersion
                                         , return.y = de_return.y
                                         , anovalike = de_anovalike)
  })
}

# Get differentially expressed genes ---
de <- lapply(de, function(x) {
  x$sig <- tryCatch(RNAseqRtools::getDEgenes(x
                                             , fdrTh = fdrTh
                                             , fcTh = fcTh)
                    , error = function(e) {message(e);return(NA)})
  return(x)
})

# Save results ---
if(save_table) {
  for(i in names(de)) {
    outfile <- paste0(path_results,"/dge-edger.",de_method,".",i,".table.txt.gz")
    message(" -- writing to: ", outfile)
    write.table(de[[i]]$table
                , file = gzfile(outfile)
                , row.names = F
                , col.names = T
                , sep = "\t"
                , quote = F)
  }
}

na_idx <- lapply(de, function(x) !is.data.frame(x$sig) && is.na(x$sig))
na_idx <- which(unlist(na_idx))
if(length(na_idx)!=0) de <- de[-na_idx] # filter non significant results

for(i in names(de)) {
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",i,".fcTh_",fcTh,".fdrTh_",fdrTh,".txt.gz")
  message(" -- writing to: ", outfile)
  write.table(de[[i]]$sig
              , file = gzfile(outfile)
              , row.names = F
              , col.names = T
              , sep = "\t"
              , quote = F)
}

if(save_excel) {
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",de_analysis,".fcTh_",fcTh,".fdrTh_",fdrTh,".xlsx")
  RNAseqRtools::saveXLSresEdgeR(res = lapply(de,"[[","sig")
                                , outfile = outfile
                                , force = T)
  
}
