#!/usr/bin/env Rscript
#
# collect-dmr-info.r
#
# 
args <- commandArgs(trailingOnly = T)

dmr_pattern <- as.character(args[1])

if(length(args)!=1) {
  message("[!] Invalid argument")
  message("Usage: collect-dmr-info.r <dmr_pattern>")
  quit(save = "no", status = 0, runLast = TRUE)
}
path_dss <- paste0(dmr_pattern,".txt")
path_region_gene_map <- paste0(dmr_pattern,"_rGREAT/region_gene_map_default.gz")
path_region_ccre_map <- paste0(dmr_pattern,".ccREs_map.gz")
path_hyper_names     <- paste0(dmr_pattern,".hyper.names")
path_hypo_names      <- paste0(dmr_pattern,".hypo.names")

hyper_names <- read.delim(path_hyper_names, header = F, col.names = "idx")
hyper_names$status <- "hyper"
hypo_names <- read.delim(path_hypo_names, header = F, col.names = "idx")
hypo_names$status <- "hypo"

dmr <- rbind.data.frame(hyper_names,hypo_names)
# dss
dss <- read.delim(path_dss, header = T)
dss$idx <- with(dss, paste0(chr,":",start,"-",end))
dmr <- merge(dmr, dss[,-c(1:3)], by = "idx", all.x = T)
# gene map
region_gene_map <- read.delim(path_region_gene_map, header = T)
dmr <- merge(dmr, region_gene_map, by = "idx", all = T)
# ccre map
region_ccre_map <- read.delim(path_region_ccre_map, header = T, stringsAsFactors = F)
region_ccre_map$ccRE_grouped <- region_ccre_map$ccRE
region_ccre_map$ccRE_grouped[grep("PLS",region_ccre_map$ccRE_grouped)] <- "Promoters (PLS)"
region_ccre_map$ccRE_grouped[grep("pELS",region_ccre_map$ccRE_grouped)] <- "Enhancers (pELS)"
region_ccre_map$ccRE_grouped[grep("dELS",region_ccre_map$ccRE_grouped)] <- "Enhancers (dELS)"
region_ccre_map$ccRE_grouped[grep("CTCF",region_ccre_map$ccRE_grouped)] <- "CTCF-bound"
dmr <- merge(dmr, region_ccre_map, by = "idx", all = T)

outfile <- paste0(dmr_pattern,".dmr_info.txt.gz")
message(" -- saving to: ", outfile)
write.table(dmr
            , file = gzfile(outfile)
            , row.names = F
            , col.names = T
            , sep = "\t"
            , quote = F)
