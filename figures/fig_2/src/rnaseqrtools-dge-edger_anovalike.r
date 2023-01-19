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
path_results <- "results/dge-edger_anovalike" # ASSIGN_PATH_RESULTS

if(!dir.exists(path_results)) dir.create(path_results, recursive = T)

# Differential gene expression analysis ---
de_method = "qlf" # qlf, lrt, exact
de_robust.dispersion = TRUE
de_return.y = FALSE
de_formula = "~Genotype*Stage" # de_formula = "~time"
de_cf = list("timecourse_anovalike"=c(2:8))
de_contrast  = NULL
de_anovalike = TRUE
de_analysis  = "timecourse_anovalike" # dge

fdrTh  = 0.001
fcTh   = 1.5
save_table = TRUE
save_excel = TRUE

s33d = 1991
# read object ---
dgelist_obj <- readRDS(path_dgelist_obj)
# reorder factors
dgelist_obj$samples$Stage    <- factor(dgelist_obj$samples$Stage, levels = c("ESC","EpiLC","ME24h","ME48h"))
dgelist_obj$samples$Genotype <- factor(dgelist_obj$samples$Genotype, levels = c("WT","3BKO"))

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
  x$sig <- RNAseqRtools::getDEgenes(x
                                    , fdrTh = fdrTh
                                    , fcTh = fcTh)
  return(x)
})

# Save results ---
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

if(save_excel) {
  outfile <- paste0(path_results,"/dge-edger.",de_method,".",de_analysis,".fcTh_",fcTh,".fdrTh_",fdrTh,".xlsx")
  RNAseqRtools::saveXLSresEdgeR(res = lapply(de,"[[","sig")
                                , outfile = outfile
                                , force = T)
  
}
